# WORF Pipeline（独立自包含版本）

整合 `worf_benchmark`（fastq → QC → 比对 → counts）与 `worfscore2`（counts → diff → 图）
两大模块，一条命令全流程，可续跑、可预览，输出统一到固定目录。

> **自包含**：本文件夹即可独立运行，不依赖任何目录外代码。

## 结构

```
Worf_pipeline/
├── run_pipeline.py    # 主编排入口（唯一需要手动调用的脚本）
├── config.json        # 全局配置（路径 / 实验列表 / window / 环境）
├── targets/           # 靶点定义 on.csv（85322 on-target） / off.csv（135638 off-target）
├── modules/           # step1(BAM→parquet) / step2(parquet→target_counts)
└── plots/             # plot_by_n_s(散点) / plot_population(群体统计)
```

## 全流程

```
raw fastq ─fastp→ clean fastq ─minimap2(hg38)→ SAM ─samtools sort→ BAM
   → step1(BAM→parquet) → step2(on/off.csv + window→ target_counts)
   → diff(sample−ctrl, 排序, center, 合并 on/off, filtered)
   → plot_by_n_s(散点) + plot_population(群体统计)
```

## 统一输出目录（workflow_root，默认 `/data/lulab_commonspace/guozehua/worf_pipeline`）

```
{workflow_root}/
├── 1.QC/          fastp: {exp}/{sample}_clean_{1,2}.fq.gz, fastp.json/html
├── 2.Alignment/   minimap2+samtools: {exp}/{sample}_aligned.sorted.bam (+.bai)
├── 3.Counts/      {exp}/{sample}_reads.parquet,
│                  {exp}/{kind}_target_counts_{frac}.csv  (kind=on/off, frac=7 种饱和度)
├── 4.Diff/        {exp}/{exp}_{kind}_diff[{_frac}].csv  (按实验归档)
│                  {exp}/{exp}_diff[{_frac}].csv, {exp}_diff[{_frac}]_filtered.csv
└── 5.Plots/       scatter_slim/{exp}/{frac}/ 与 population/{exp}/{frac}/
                   （无时间戳，直接覆盖同名文件）
```

> 每个分部互不污染；`4.Diff` 按实验归档；`5.Plots` 按 大类/样品/depth 组织，**直接覆盖**（方便反复微调）。
> counts/diff/plot 均按 **7 种饱和度深度（frac）** 分别计算：`0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 1`（frac=1 为全量，文件名无后缀）。
> 每个步骤幂等：输出已存在则跳过，可用 `--force` 强制重跑。

## 用法

```bash
cd ~/Projects/Lab/Worf_pipeline
conda activate worf_env

# 全部实验（config.experiments）全流程
python run_pipeline.py

# 扫描 raw_root 下全部实验
python run_pipeline.py --all

# 指定实验 + 指定步骤（续跑）
python run_pipeline.py --exp 20260508WORFT5-SNAP --steps qc,align,count

# 强制重算某实验的 diff 与图
python run_pipeline.py --exp 20260508WORFT5-SNAP --steps diff,plot --force

# 预览将执行的命令（不实际执行）
python run_pipeline.py --exp 20260508WORFT5-SNAP --dry-run

# 换 window / 绘图参数
python run_pipeline.py --window 30 --exclude-chrm --cap-auto

# 为某实验单独指定绘图 N 值（覆盖 config.n_rows）
python run_pipeline.py --exp 20260508WORFT5-SNAP --steps plot --n-rows-per-exp 20260508WORFT5-SNAP:10000
```

### 参数

| 参数 | 说明 |
|------|------|
| `--steps qc,align,count,diff,plot` | 要执行的步骤（默认全部） |
| `--exp NAME [NAME...]` | 指定实验；对照组始终自动包含 |
| `--all` | 扫描 `raw_root` 下全部实验目录 |
| `--window N` | 覆盖 step2 扩窗 bp（默认 150） |
| `--force` | 覆盖已存在输出 |
| `--n-rows N` | 全局绘图 N 值（默认取 `config.n_rows`） |
| `--n-rows-per-exp EXP:N` | 按实验指定绘图 N 值，可多次，优先级最高 |
| `--exclude-chrm` / `--cap-auto` | 绘图选项（排除 chrM / 自动裁剪） |
| `--y-cap iqr\|p99\|max` | 散点图稳健 y 上限策略（默认 `iqr` = Q3+1.5×IQR） |
| `--no-population` | 跳过群体统计图 |
| `--keep-sam` | 保留中间 SAM（默认排序后删除） |
| `--dry-run` | 只打印命令，不执行 |

> **N 值说明**：N 只影响第 5 步绘图（每 on/off 组取前 N 行），前 4 步（qc/align/count/diff）完全不用 N。每个实验的 N 在 `config.n_rows` 中定义，可用 `--n-rows-per-exp EXP:N` 覆盖（同实验可多次生成不同 N 的图）。plot 输出目录带 `_N{n}` 后缀，如 `5.Plots/scatter_slim/{exp}/{frac}_N10000/`，不同 N 的图互不覆盖。

> **y 轴 scale 说明**：scatter 一律线性轴，ctrl/sample 共用同一 y 上限（默认 `--y-cap iqr` = Q3+1.5×IQR，可换 `p99` / `max`）。IQR 法能避免零星大点把 y 轴撑开——**数据不删除，仅裁剪显示**（超上限的点被 ylim 裁掉不画）。每个子文件夹的 `max50.log` 记录完整分布统计（mean/median/Q1/Q3/IQR/P99/abs_max + ctrl/sample 各自 top50 + 实际 ylim），用于诊断 scale 是否合理。

## 配置（config.json）

| 键 | 说明 |
|----|------|
| `conda_env` / `conda_prefix` | 运行环境（worf_env：含 fastp/minimap2/samtools/matplotlib/scipy） |
| `ref_mmi` | 参考基因组 minimap2 索引（hg38.mmi） |
| `raw_root` | 原始数据根（含各实验目录，fastq 在 `{exp}/00.mergeRawFq/{sample}/`） |
| `control_exp` | 固定对照组（20260127humangenecontrol） |
| `experiments` | 要分析的实验组列表 |
| `targets.on/off` | 靶点定义 CSV（相对本文件夹，可自定） |
| `window` | step2 匹配扩窗 bp（默认 150） |
| `mapq_threshold` | step2 MAPQ 阈值（默认 30） |
| `threads` | minimap2 / samtools 线程数（默认 8） |

## 运行时环境

- conda env **`worf_env`**（具备全部依赖：fastp, minimap2, samtools,
  pandas, pyarrow, pysam, matplotlib, scipy）
- 参考基因组索引：`/data/lulab_commonspace/guozehua/Worf/references/hg38.mmi`
- 原始 fastq 约定：`{raw_root}/{exp}/00.mergeRawFq/{sample}/{sample}_raw_{1,2}.fq.gz`

## 数据语义（重要）

- **on / off 是“靶点类型”**：on-target 基因位点（85322 个）vs off-target 全基因组背景位点（135638 个），二者位点几乎不重叠。
- **对照 / 实验是“样本（日期）”维度**：对照组固定 `20260127humangenecontrol`，实验组为其后各批次。
- diff 公式：`diff = read_count_sample − read_count_ctrl`，按 target_id 对齐；ON 按 diff 降序、OFF 按 diff 升序；合并后加 `on_target` 列；filtered 版删除 OFF 组中 `ctrl > 10` 的背景行。
