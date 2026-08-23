# 冻结注册清单（Frozen Registration Manifest）

**冻结时间戳：2026-08-14 22:43:59 CST**（Asia/Shanghai，UTC+8）
**用途：上传 OSF/Zenodo 作为时间戳预注册（建议 OSF registries → 'Frozen Analysis Plan'）**

## 1. 冻结的检验全集（在外部队列检验前固定）

### 7 个预冻结机制家族（成员定义）
- 文件：results/GSE219280_same_celltype_fixed_programs/predefined_gene_set_membership.tsv
### 17 个预冻结候选基因
BAG3, HSP90AA1, HSP90AB1, CHORDC1, ARHGAP35, PARK7, CCT5, ARHGAP39, COX5B, COX7C, MAP1B, MAPT, BCR, CCT8, COA3, DOCK7, STARD13
### 9 个 MAGMA 基因集（7 家族 + CAND17 + LE_strict13 并集 + pool3）
- 文件：external/magma/program_gene_sets.gsa

## 2. 关键内容 SHA256 哈希（冻结时刻）
- `results/GSE219280_same_celltype_fixed_programs/predefined_gene_set_membership.tsv`：`2fcdf4a7e9cab7e295c08c4de8ac43965e52c9290513e2b075ba911164717542`
- `external/magma/program_gene_sets.gsa`：`04123d2462c2df32c6d7497dc27df54f7f29dbc970f9df91b1ff8a03d7b9c7fa`
- `scripts/binomial_concordance_tests.R`：`73a917718e007cadccb1cb90e0ac0af66ea0bcc37b9a1f393bf2b07494b9d245`
- `scripts/ruf_gate_waterfall.R`：`0a389c54415ef9629ad1a19acf60db207c6cb2e1b32226f276d1187950d7ea1d`
- `scripts/gate_power_and_tCI.R`：`1977169c57e9d4b84641a3257fe36f14c50b59db07b9e592733d7512bf6b7723`
- `scripts/tdp43_perturbation_validation.R`：`f8ab0f3784ede099fbbaae7b5f793fe93f8979f491b8cfee4e3a7fc29cc9ea54`
- `scripts/run_magma_als_gwas.R`：`c6e5b7428863f936e4314cdceaaf5a823bfc1a356a281579a9cd5ff9c03989e5`
- `scripts/build_magma_inputs.R`：`22b082453cd4ef687215b8675c2e770146097c752d4922f03722df8fe0b77e57`

## 3. 冻结的统计门槛
- Ruf 五重门槛：全局家族FDR<0.05 → 供体LOO≥90% → 技术+年龄同向 → SVA12同向 → SVA12全局FDR<0.05
- 三队列完整门槛：三队列同向 + 各自全局FDR + 稳健性（仅 21 条胶质 ALS 轴可估计）
- 严格复制定义：双队列同向 + 双队列全局FDR + LOO≥90%
- MAGMA 基因集：竞争性检验，Bonferroni 0.05/17（单基因）

## 4. 声明
以上哈希对应文件在冻结时刻后未经修改；后续任何修改须在新版本文件中进行并重新注册。
