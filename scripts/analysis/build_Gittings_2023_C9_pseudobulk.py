#!/usr/bin/env python3

"""Build donor x region x major-cell-type pseudobulk from raw.X counts."""

from __future__ import annotations

import csv
import datetime as dt
import os
import pathlib

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse


ROOT = pathlib.Path(
    os.environ.get("ALS_FTD_PROJECT_ROOT", pathlib.Path(__file__).resolve().parents[2])
).expanduser().resolve()
RAW = ROOT / "data/raw/Gittings_2023_C9_snRNA"
OUT = ROOT / "results/Gittings_2023_C9_snRNA"
MANIFEST = RAW / "cellxgene_h5ad_manifest.tsv"
STATUS = OUT / "status.tsv"


def timestamp():
    return dt.datetime.now(dt.timezone.utc).astimezone().strftime("%Y-%m-%dT%H:%M:%S%z")


def status(stage, state, detail, progress):
    with STATUS.open("a", encoding="utf-8", newline="") as handle:
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerow(
            [timestamp(), stage, state, detail, f"{progress:.1f}"]
        )
    print(f"[{timestamp()}] {stage} | {state} | {detail} | {progress:.1f}%", flush=True)


def diagnosis(value):
    text = str(value).strip().lower()
    if text == "normal":
        return "Control"
    if "frontotemporal dementia" in text and "amyotrophic" not in text:
        return "C9_FTD"
    if "with or without frontotemporal dementia" in text or "als/ftd" in text:
        return "C9_ALS_FTD"
    if "amyotrophic lateral sclerosis" in text:
        return "C9_ALS"
    raise ValueError(f"Unmapped disease: {value}")


def major_celltype(value):
    text = str(value).strip().lower()
    if text.startswith("excitatory"):
        return "Exc"
    if text.startswith("inhibitory"):
        return "Inh"
    if any(token in text for token in (
        "intratelencephalic", "corticothalamic", "near projecting",
        "extratelencephalic",
    )):
        return "Exc"
    if "interneuron" in text or "lamp5" in text:
        return "Inh"
    if text == "astrocytes":
        return "Astro"
    if text == "opcs":
        return "OPC"
    if text == "oligodendrocytes":
        return "Oligo"
    if text == "microglia":
        return "Micro"
    if "endo" in text or "pericyte" in text:
        return "Vascular"
    raise ValueError(f"Unmapped cluster: {value}")


def age_number(value):
    text = str(value)
    for token in text.replace("-year-old stage", "").split():
        try:
            return float(token)
        except ValueError:
            pass
    return np.nan


def main():
    manifest = pd.read_csv(MANIFEST, sep="\t", dtype=str, keep_default_na=False)
    complete = manifest.apply(
        lambda row: pathlib.Path(row["local_path"]).exists()
        and pathlib.Path(row["local_path"]).stat().st_size == int(row["expected_bytes"]), axis=1
    )
    if complete.sum() != 4:
        raise RuntimeError(f"Pseudobulk requires 4 complete assets; found {complete.sum()}")
    status("pseudobulk", "started", "Aggregate integer raw.X by donor-region-major cell type", 71)

    all_counts = []
    all_samples = []
    all_gene_ids = []
    symbol_by_gene_id = {}
    for file_index, (_, record) in enumerate(manifest.iterrows(), 1):
        path = pathlib.Path(record["local_path"])
        data = ad.read_h5ad(path, backed="r")
        if data.raw is None:
            raise RuntimeError(f"raw.X missing: {path}")
        obs = data.obs.copy()
        required = {"donor_id", "disease", "sex", "development_stage"}
        missing = required.difference(obs.columns)
        if missing:
            raise RuntimeError(f"obs missing {sorted(missing)} in {path}")
        cluster_column = "cluster" if "cluster" in obs.columns else "cluster_layers"
        if cluster_column not in obs.columns:
            raise RuntimeError(f"author cluster/cluster_layers missing in {path}")
        obs["diagnosis"] = obs["disease"].map(diagnosis)
        obs["cell_type"] = obs[cluster_column].map(major_celltype)
        obs["region"] = record["region"]
        obs["partition"] = record["cohort_partition"]
        obs["donor_id"] = obs["donor_id"].astype(str)
        obs["sex"] = obs["sex"].astype(str)
        obs["age"] = obs["development_stage"].map(age_number)
        group_columns = ["donor_id", "diagnosis", "region", "cell_type", "partition", "sex", "age"]
        group_index = pd.MultiIndex.from_frame(obs[group_columns])
        codes, levels = pd.factorize(group_index, sort=True)
        n_groups = len(levels)

        gene_ids = np.asarray(data.raw.var_names.astype(str))
        # In these CELLxGENE exports raw.var_names are Ensembl IDs while the
        # unfortunately named `gene_ids` column contains HGNC-style symbols.
        symbol_column = "gene_ids" if "gene_ids" in data.raw.var.columns else "feature_name"
        symbols = np.asarray(data.raw.var[symbol_column].astype(str))
        if len(set(gene_ids)) != len(gene_ids):
            raise RuntimeError(f"Duplicate raw gene IDs in {path}")
        all_gene_ids.append(gene_ids)
        for gene_id, symbol in zip(gene_ids, symbols):
            previous = symbol_by_gene_id.setdefault(gene_id, symbol)
            if previous != symbol:
                raise RuntimeError(
                    f"Conflicting symbol for {gene_id}: {previous} vs {symbol} in {path}"
                )

        group_counts = np.zeros((n_groups, len(gene_ids)), dtype=np.int64)
        group_nuclei = np.bincount(codes, minlength=n_groups).astype(np.int64)
        matrix = data.raw.X
        for start in range(0, data.n_obs, 4096):
            end = min(start + 4096, data.n_obs)
            block = matrix[start:end, :]
            if not sparse.issparse(block):
                block = sparse.csr_matrix(block)
            if block.data.size and not np.allclose(block.data, np.round(block.data)):
                raise RuntimeError(f"Noninteger raw.X values in {path} rows {start}:{end}")
            local_codes = codes[start:end]
            selector = sparse.csr_matrix(
                (np.ones(end - start), (local_codes, np.arange(end - start))),
                shape=(n_groups, end - start),
            )
            group_counts += np.asarray((selector @ block).todense(), dtype=np.int64)

        level_frame = levels.to_frame(index=False)
        level_frame.columns = group_columns
        level_frame["n_nuclei"] = group_nuclei
        level_frame["total_counts"] = group_counts.sum(axis=1)
        level_frame["dataset_version_id"] = record["dataset_version_id"]
        level_frame["sample_id"] = level_frame.apply(
            lambda row: "|".join([
                row["donor_id"], row["diagnosis"], row["region"], row["cell_type"]
            ]), axis=1
        )
        if level_frame["sample_id"].duplicated().any():
            raise RuntimeError(f"Duplicate pseudobulk IDs within {path}")
        all_counts.append(group_counts)
        all_samples.append(level_frame)
        data.file.close()
        status(
            "pseudobulk", "running",
            f"aggregated {file_index}/4 {record['dataset_version_id']}; groups={n_groups}",
            71 + 4 * file_index / 4,
        )

    samples = pd.concat(all_samples, ignore_index=True)
    # Each CELLxGENE asset can omit genes that are all-zero within that asset.
    # Align raw-count pseudobulks on the union of raw.var identifiers and use
    # zero for an identifier absent from one asset.
    union_gene_ids = np.asarray(sorted(set().union(*map(set, all_gene_ids))))
    union_index = {gene_id: index for index, gene_id in enumerate(union_gene_ids)}
    aligned_counts = []
    for group_counts, gene_ids in zip(all_counts, all_gene_ids):
        aligned = np.zeros((group_counts.shape[0], len(union_gene_ids)), dtype=np.int64)
        columns = np.fromiter(
            (union_index[gene_id] for gene_id in gene_ids), dtype=np.int64,
            count=len(gene_ids),
        )
        aligned[:, columns] = group_counts
        aligned_counts.append(aligned)
    counts = np.vstack(aligned_counts)
    if samples["sample_id"].duplicated().any():
        duplicates = samples.loc[samples["sample_id"].duplicated(False), "sample_id"].tolist()
        raise RuntimeError(f"Cross-asset duplicate pseudobulk IDs: {duplicates[:10]}")
    if counts.shape != (len(samples), len(union_gene_ids)):
        raise RuntimeError("Pseudobulk dimension mismatch")

    samples.to_csv(OUT / "pseudobulk_sample_metadata.tsv", sep="\t", index=False)
    union_symbols = np.asarray([symbol_by_gene_id[gene_id] for gene_id in union_gene_ids])
    genes = pd.DataFrame({"gene_id": union_gene_ids, "gene": union_symbols})
    genes.to_csv(OUT / "pseudobulk_gene_metadata.tsv", sep="\t", index=False)
    with h5py.File(OUT / "pseudobulk_raw_counts.h5", "w") as handle:
        handle.create_dataset(
            "counts", data=counts.T, compression="gzip", compression_opts=4,
            shuffle=True, chunks=(min(1024, counts.shape[1]), min(16, counts.shape[0])),
        )
        handle["counts"].attrs["orientation"] = "genes_by_samples"
        handle.create_dataset("gene_id", data=np.asarray(union_gene_ids, dtype="S"))
        handle.create_dataset("gene_symbol", data=np.asarray(union_symbols, dtype="S"))
        handle.create_dataset("sample_id", data=np.asarray(samples["sample_id"], dtype="S"))
    qc = pd.DataFrame([{
        "assets": 4, "nuclei": int(samples["n_nuclei"].sum()),
        "pseudobulk_samples": len(samples), "donors": samples["donor_id"].nunique(),
        "genes": len(union_gene_ids), "duplicate_sample_ids": int(samples["sample_id"].duplicated().sum()),
        "asset_gene_counts": ";".join(map(str, map(len, all_gene_ids))),
        "gene_alignment": "union_by_raw_var_gene_id_absent_asset_gene_zero_filled",
        "all_counts_integer": bool(np.issubdtype(counts.dtype, np.integer)),
        "minimum_nuclei": int(samples["n_nuclei"].min()),
        "maximum_nuclei": int(samples["n_nuclei"].max()),
    }])
    qc.to_csv(OUT / "pseudobulk_build_QC.tsv", sep="\t", index=False)
    status(
        "pseudobulk", "complete",
        f"nuclei={qc.iloc[0]['nuclei']};samples={len(samples)};donors={samples['donor_id'].nunique()};genes={len(union_gene_ids)}",
        76,
    )


if __name__ == "__main__":
    main()
