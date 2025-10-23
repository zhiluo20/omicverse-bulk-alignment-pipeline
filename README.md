# omicverse-bulk-alignment-pipeline
#  RNA-seq Bulk Alignment Pipeline — Class 集成与批处理支持

## 摘要（Summary）
本 PR 将原有分散脚本封装为 `omicverse.bulk.Alignment` 统一入口，打通 **SRA → FASTQ → QC(fastp) → STAR → featureCounts** 全链条，支持**并发**、**幂等跳过**与**一致的产出结构**，并新增：
- `Alignment.fetch_metadata()`：GEO/ENA 元数据获取与 RunInfo 生成  
- `Alignment.prefetch()`：多路 SRA 并发下载（带进度、断点与官方校验）  
- `Alignment.fasterq()`：`fasterq-dump` 并行（按 SRR 隔离输出与 tmp）  
- `Alignment.fastp()`：批处理 QC，缓存检测与报告收集  
- `Alignment.star_align()`：批量 STAR（索引自动选择/缓存、BAM 幂等、SRR 命名软链接）  
- `Alignment.featurecounts()`：批量计数与**自动合并矩阵**（自动推断 GTF；返回矩阵路径）

> 产物命名标准化：  
> - STAR 目录中保留官方名 `Aligned.sortedByCoord.out.bam`，同时**额外暴露** `SRR.sorted.bam`（软链接/拷贝）以避免合并矩阵时列名冲突。  
> - featureCounts 单样本表统一为 `<counts_root>/<SRR>/<SRR>.counts.txt`；合并矩阵为 `<counts_root>/matrix.<by>.csv`。

---

## 变更范围（Scope）
- 新增/调整的模块（示例，实际以本 PR diff 为准）：
  - `bulk/alignment.py`（核心类与适配）
  - `bulk/sra_prefetch.py`（并发 prefetch 与进度条最小侵入改造）
  - `bulk/sra_fasterq.py`（批处理 fasterq；幂等与重试）
  - `bulk/qc_fastp.py`（fastp 批处理与产物检测）
  - `bulk/star_step.py`（基于现有 `star_tools` 的批处理适配；SRR 软链接）
  - `bulk/count_step.py` / `bulk/count_tools.py`（批量 featureCounts 与合并矩阵；列名统一为 SRR）
  - `bulk/tools_check.py`（`which_or_find`、`merged_env` 等工具）
- `bulk/__init__.py` 导出 `Alignment`, `AlignmentConfig`

---

## 依赖（Dependencies）

### 系统 & 生信工具
- **sra-tools**（需要 `prefetch`, `vdb-validate`, `fasterq-dump`）
- **samtools**（BAM index）
- **STAR**
- **subread**（提供 `featureCounts`）
- **fastp**
- 建议：`pigz`（gzip 加速）、`aria2`（下载加速，可选）

### Python 包
- 必需：`pandas`, `tqdm`, `numpy`, `requests`, `lxml`
- 可能用到（视 meta 抓取实现）：`biopython`

**见本 PR 附带的 `environment.yml`（bioconda 优先）**。

---

## 目录结构（Outputs Layout）
运行后默认输出如下（以 `work/` 为根）：
```
work/
 ├── prefetch/         # prefetch 的 .sra
 │   └── SRRxxxxxx/SRRxxxxxx.sra
 ├── fasterq/
 │   └── SRRxxxxxx/
 │       ├── SRRxxxxxx_1.fastq.gz
 │       └── SRRxxxxxx_2.fastq.gz
 ├── fastp/
 │   └── SRRxxxxxx/
 │       ├── SRRxxxxxx_clean_1.fastq.gz
 │       ├── SRRxxxxxx_clean_2.fastq.gz
 │       ├── SRRxxxxxx.fastp.json
 │       └── SRRxxxxxx.fastp.html
 ├── star/
 │   └── SRRxxxxxx/
 │       ├── Aligned.sortedByCoord.out.bam      # 官方文件名（保留）
 │       ├── Aligned.sortedByCoord.out.bam.bai
 │       ├── SRRxxxxxx.sorted.bam               # SRR 命名（软链接/拷贝）
 │       └── SRRxxxxxx.sorted.bam.bai
 └── counts/
     ├── SRRxxxxxx/SRRxxxxxx.counts.txt
     └── matrix.auto.csv                        # 合并矩阵（行=gene_id, 列=SRR）
```

---

## 环境变量 / 配置（Configuration）
- `NCBI_SETTINGS`（可选）：SRA 工具配置路径  
- `TMPDIR`（可选）：大文件临时目录  
- `FC_GTF_HINT`（可选）：当无法从 STAR index 推断 GTF 时，提供 GTF 路径提示

`AlignmentConfig` 关键字段（有默认）：
- `work_root`, `prefetch_root`, `fasterq_root`, `fastp_root`,  
  `star_index_root`, `star_align_root`, `counts_root`
- `threads`（并发控制，与内部每任务线程需平衡）
- `memory`（传递给 fasterq 的 `--mem`，如 `"8G"`）
- `gzip_fastq`（fasterq 输出是否压缩）

---

## 端到端用法（Usage）
```python
from omicverse.bulk import Alignment, AlignmentConfig

cfg = AlignmentConfig(
    work_root="work",
    prefetch_root="work/prefetch",
    fasterq_root="work/fasterq",
    fastp_root="work/fastp",
    star_index_root="index",
    star_align_root="work/star",
    counts_root="work/counts",
    threads=16,
    memory="8G",
    gzip_fastq=True,
)

aln = Alignment(cfg)

meta = aln.fetch_metadata("GSE157103")
sra_paths = aln.prefetch(meta["srr_list"], max_concurrent=4)
fq_pairs = aln.fasterq(meta["srr_list"])
qc_results = aln.fastp(fq_pairs)
pairs_for_star = [(srr, c1, c2) for (srr, c1, c2, _, _) in qc_results]

bam_triples = aln.star_align(
    pairs_for_star,
    gencode_release="v44",
    sjdb_overhang=149,
    accession_for_species=None,
    max_workers=2,
)

fc_out = aln.featurecounts(
    bam_triples,
    simple=True,
    by="auto",
    threads=8,
)
print("merged matrix:", fc_out.get("matrix"))
```

---

## 复现实验（Repro Checklist）
- [ ] 同一批次 SRR **重复运行**：所有阶段出现 `[SKIP]`，不重复生成产物  
- [ ] fasterq 失败重试有效（网络不稳时自动切换本地 `.sra` 输入）  
- [ ] STAR 产物包含**官方名**与**SRR 命名软链接**，二者指向同一数据  
- [ ] featureCounts 合并矩阵列名为 **SRR**，无重复列名冲突  
- [ ] `counts/matrix.auto.csv` 行数 ≥ 单样本 gene 行数上限，且 `gene_id` 非空  

---

## 并发与性能（Tuning）
- **机器核数 N**：`max_workers × per-sample threads ≤ N`（含超线程时适当打折）  
- `prefetch`：外层并发（例 `max_concurrent=4`），单个下载内部 0.25s 轮询进度  
- `fasterq-dump`：建议 `--mem 8–16G`、`threads_per_job 12–24`；并发样本数谨慎  
- `STAR`：**内存敏感**；大型基因组建议单样本 8–16 线程，并发 1–2  
- `featureCounts`：`-T` 适中（8–16），I/O 是主要瓶颈  

---

## 向后兼容性（Compatibility）
- 原有脚本可继续**独立调用**；类方法只是薄适配，不改变原业务逻辑  
- STAR 输出增加了 SRR 软链接，不影响原有消费方；有助于矩阵列名唯一化  

---

## 测试（Test Plan）
- [ ] 小样本（2–3 SRR）本地端到端测试  
- [ ] 中等批量（8–12 SRR）并发参数  
- [ ] 断点续跑（kill 后重启）产物与日志一致性  
- [ ] 手动删除某一步产物，仅重做该步（其余 `[SKIP]`）  
- [ ] 仅保留 STAR 官方 BAM 名时，`_normalize_bam` 能补齐 `SRR.sorted.bam` 软链接  
- [ ] `FC_GTF_HINT` 指向 GTF 时，可绕过自动推断  

---

## 未来工作（Next）
- [ ] 增加 **CLI**（`ov-bulk align ...`）封装类方法  
- [ ] 支持 **STAR 索引自动下载/构建**（gencode/ensembl）与缓存记录  
- [ ] 整合 **salmon**/**kallisto** 作为可选轻量计数  
- [ ] 增加 **md5/sha256** 与产物 manifest，便于审计与复现  
- [ ] GitHub Actions 做最小 CI：lint + 单元测试 + 小数据集集成测试  

---

## 故障排查（Troubleshooting）
- **MergeError: duplicate columns** → 已在 `count_tools.py` 中合并前将计数列名改为 **SRR**  
- **fasterq 退出码 3 / 输出缺失** → 网络不稳或 S3 超时；已加入重试与切换本地 `.sra`  
- **STAR 内存不足** → 降低 `threads` 或减少 `max_workers`；必要时调整 `--limitGenomeGenerateRAM`  
- **找不到 GTF** → 显式传 `gtf=` 或设 `FC_GTF_HINT`；或确保 STAR index 上级 `_cache/` 下有 gtf  
- **权限/软链接问题** → 若文件系统不支持 symlink，代码自动回退为拷贝策略  

---

## Checklist（提交前）
- [ ] 代码通过 `flake8`/`black`（或项目既有规范）  
- [ ] 大文件未纳入仓库（.sra/.bam/.fastq.gz 等）  
- [ ] 文档与示例路径与默认配置一致  
- [ ] 在 `CHANGELOG.md` 或本 PR 中清楚记录变更与迁移说明  
