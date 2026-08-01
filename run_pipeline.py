#!/usr/bin/env python3
"""
WORF Targeted-Sequencing Unified Pipeline
=========================================
把 worf_benchmark (fastq→QC→align→counts) 与 worfscore2 (counts→diff→plots) 两大
模块整合为一条可编排、可续跑、输出统一的全流程。

【自包含】本文件夹即可独立运行，不依赖任何目录外代码：
  Worf_pipeline/
  ├── run_pipeline.py    # 主编排入口
  ├── config.json        # 全局配置
  ├── targets/           # on.csv / off.csv 靶点定义
  ├── modules/           # step1(BAM→parquet) / step2(parquet→target_counts)
  └── plots/             # plot_by_n_s(散点) / plot_population(群体统计)

工作目录（固定）:
  {workflow_root}/
  ├── 1.QC/           fastp 质控: {exp}/{sample}_clean_{1,2}.fq.gz + fastp.json
  ├── 2.Alignment/    minimap2+samtools: {exp}/{sample}_aligned.sorted.bam (+.bai)
  ├── 3.Counts/       step1(parquet)+step2(on/off target_counts): {exp}/
  ├── 4.Diff/         {exp}_on_diff.csv / {exp}_off_diff.csv / {exp}_diff.csv
  │                          / {exp}_diff_filtered.csv
  └── 5.Plots/        scatter_slim/{ts}/ 与 population/{ts}/

用法示例:
  # 全部实验全流程
  python run_pipeline.py --all --steps qc,align,count,diff,plot
  # 单实验续跑（只跑缺失步骤）
  python run_pipeline.py --exp 20260508WORFT5-SNAP
  # 只重算 diff + 图
  python run_pipeline.py --exp 20260508WORFT5-SNAP --steps diff,plot --force
  # 预览将执行的命令（不实际执行）
  python run_pipeline.py --exp 20260508WORFT5-SNAP --dry-run
"""
import argparse
import json
import os
import subprocess
import sys
from datetime import datetime
from pathlib import Path

PIPELINE_ROOT = Path(__file__).resolve().parent
CONFIG_PATH = PIPELINE_ROOT / "config.json"

STEPS_ORDER = ["qc", "align", "count", "diff", "plot"]

# saturation 深度比例（counts 按此 7 种深度分别计算 diff / plot）
# 对应 3.Counts/{exp}/{kind}_target_counts_{frac}.csv
FRACTIONS = ["0.001", "0.01", "0.05", "0.1", "0.2", "0.5", "1"]

args_singleton = None  # 全局 CLI 参数（供步骤函数读取），在 main 中赋值


# ---------------------------------------------------------------- helpers
def log(msg: str):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def load_config(path: str = None) -> dict:
    cfg_path = Path(path) if path else CONFIG_PATH
    with open(cfg_path, "r", encoding="utf-8") as f:
        cfg = json.load(f)
    # 解析相对路径（相对 pipeline 根目录）
    for key in ("targets",):
        for k, v in cfg[key].items():
            if not os.path.isabs(v):
                cfg[key][k] = str((PIPELINE_ROOT / v).resolve())
    return cfg


def env_bin(cfg: dict, tool: str) -> str:
    """返回 conda 环境内工具/解释器绝对路径。"""
    base = Path(cfg["conda_prefix"]) / "envs" / cfg["conda_env"] / "bin"
    return str(base / tool)


def run(cmd: list, cwd: str = None, desc: str = "") -> bool:
    log(f"$ {' '.join(cmd)}" if not desc else f"{desc}\n$ {' '.join(cmd)}")
    if args_singleton and getattr(args_singleton, "dry_run", False):
        log("[dry-run] 跳过执行")
        return True
    try:
        subprocess.run(cmd, cwd=cwd, check=True)
        return True
    except subprocess.CalledProcessError as e:
        log(f"[ERROR] 命令失败 (exit={e.returncode}): {' '.join(cmd)}")
        return False


# ---------------------------------------------------------------- discovery
# 历史命名既有 "00.mergeRawFq"（0508/0628/control），也有 "00.mergerRawFq"（0615，多一个 r）
RAW_FQ_DIR_NAMES = ("00.mergeRawFq", "00.mergerRawFq")


def find_samples(cfg: dict, exp: str) -> list:
    """在 {raw_root}/{exp}/00.merge( r)?RawFq/ 下寻找 sample（每个 sample 目录含 *_raw_*.fq.gz）。"""
    raw_root = Path(cfg["raw_root"])
    exp_dir = None
    for name in RAW_FQ_DIR_NAMES:
        cand = raw_root / exp / name
        if cand.is_dir():
            exp_dir = cand
            break
    if exp_dir is None:
        log(f"[WARN] 未找到实验目录: {raw_root / exp}（无 00.mergeRawFq / 00.mergerRawFq）")
        return []
    samples = []
    for d in sorted(exp_dir.iterdir()):
        if not d.is_dir():
            continue
        fq1 = sorted(d.glob("*_raw_1.fq*")) + sorted(d.glob("*_1.fq*"))
        fq2 = sorted(d.glob("*_raw_2.fq*")) + sorted(d.glob("*_2.fq*"))
        if fq1 and fq2:
            samples.append({"name": d.name, "dir": d, "fq1": fq1[0], "fq2": fq2[0]})
    return samples


def resolve_experiments(cfg: dict, args) -> list:
    """确定要处理的实验列表（含对照）。返回 [(exp, is_control), ...]"""
    if args.all:
        raw_root = Path(cfg["raw_root"])
        exps = sorted(p.name for p in raw_root.iterdir() if p.is_dir())
        exps = [e for e in exps if not e.startswith((".", "_"))]
    else:
        exps = list(args.exp) if args.exp else list(cfg["experiments"])

    control = cfg["control_exp"]
    result = []
    if control not in exps:
        exps = [control] + exps  # 对照组始终最先处理
    for e in exps:
        result.append((e, e == control))
    return result


# ---------------------------------------------------------------- steps
def step_qc(cfg: dict, exp: str, sample: dict, force: bool) -> bool:
    out_dir = Path(cfg["workflow_root"]) / "1.QC" / exp
    out_dir.mkdir(parents=True, exist_ok=True)
    clean1 = out_dir / f"{sample['name']}_clean_1.fq.gz"
    clean2 = out_dir / f"{sample['name']}_clean_2.fq.gz"
    if not force and clean1.exists() and clean2.exists():
        log(f"[skip] QC 已存在: {clean1}")
        return True
    fastp = env_bin(cfg, "fastp")
    extra = cfg.get("fastp_extra_args", "").split()
    cmd = [fastp,
           "-i", str(sample["fq1"]), "-I", str(sample["fq2"]),
           "-o", str(clean1), "-O", str(clean2),
           "-j", str(out_dir / "fastp.json"),
           "-h", str(out_dir / "fastp.html")] + extra
    return run(cmd, desc=f"[qc] {exp}/{sample['name']} fastp")


def step_align(cfg: dict, exp: str, sample: dict, force: bool) -> bool:
    out_dir = Path(cfg["workflow_root"]) / "2.Alignment" / exp
    out_dir.mkdir(parents=True, exist_ok=True)
    bam = out_dir / f"{sample['name']}_aligned.sorted.bam"
    bai = out_dir / f"{sample['name']}_aligned.sorted.bam.bai"
    if not force and bam.exists() and bai.exists():
        log(f"[skip] BAM 已存在: {bam}")
        return True

    qc_dir = Path(cfg["workflow_root"]) / "1.QC" / exp
    clean1 = qc_dir / f"{sample['name']}_clean_1.fq.gz"
    clean2 = qc_dir / f"{sample['name']}_clean_2.fq.gz"
    dry = bool(args_singleton and getattr(args_singleton, "dry_run", False))
    if not clean1.exists() and not dry:
        log(f"[ERROR] 缺少 QC 输出: {clean1}，请先跑 qc 步骤")
        return False

    sam = out_dir / f"{sample['name']}_aligned.sam"
    minimap2, samtools = env_bin(cfg, "minimap2"), env_bin(cfg, "samtools")
    threads = cfg.get("threads", 8)

    if not force and sam.exists():
        log(f"[skip] SAM 已存在: {sam}")
    else:
        cmd = [minimap2, "-ax", "sr", "-t", str(threads),
               cfg["ref_mmi"], str(clean1), str(clean2)]
        log(f"[align] {exp} minimap2\n$ {' '.join(cmd)}")
        if dry:
            log("[dry-run] 跳过执行")
        else:
            with open(sam, "w") as fout:
                try:
                    subprocess.run(cmd, stdout=fout, check=True)
                except subprocess.CalledProcessError:
                    log("[ERROR] minimap2 失败")
                    return False

    if not run([samtools, "sort", "-@", str(threads), "-o", str(bam), str(sam)],
               desc=f"[align] {exp} samtools sort"):
        return False
    if not run([samtools, "index", str(bam)], desc=f"[align] {exp} samtools index"):
        return False
    # 排序完成后删除中间 SAM（省空间，用 --keep-sam 保留）
    if not (args_singleton and getattr(args_singleton, "keep_sam", False)):
        try:
            sam.unlink()
            log(f"[align] 已删除中间 SAM: {sam}")
        except OSError:
            pass
    return True


def step_count(cfg: dict, exp: str, sample: dict, force: bool) -> bool:
    out_dir = Path(cfg["workflow_root"]) / "3.Counts" / exp
    out_dir.mkdir(parents=True, exist_ok=True)
    parquet = out_dir / f"{sample['name']}_reads.parquet"
    python = env_bin(cfg, "python")
    step1 = PIPELINE_ROOT / "modules" / "step1_extract_reads.py"
    step2 = PIPELINE_ROOT / "modules" / "step2_1_match_targets.py"
    bam = Path(cfg["workflow_root"]) / "2.Alignment" / exp / f"{sample['name']}_aligned.sorted.bam"
    dry = bool(args_singleton and getattr(args_singleton, "dry_run", False))
    if not bam.exists() and not dry:
        log(f"[ERROR] 缺少 BAM: {bam}，请先跑 align 步骤")
        return False
    window = cfg.get("window", 150)

    # step1: BAM -> parquet
    if not force and parquet.exists():
        log(f"[skip] parquet 已存在: {parquet}")
    else:
        cmd = [python, str(step1),
               "--bam", str(bam),
               "--output", str(parquet),
               "--min-mapq", str(cfg.get("min_mapq_step1", 0)),
               "--compression", "snappy",
               "--target-csv", cfg["targets"]["on"],
               "--window", str(window)]
        if not run(cmd, desc=f"[count] {exp} step1 (BAM→parquet)"):
            return False
        # step1 实际会写到 {out_dir}/{task_id}/ 子目录，归位到顶层规范路径
        if not dry and not parquet.exists():
            candidates = sorted(out_dir.rglob("*.parquet"))
            if not candidates:
                log(f"[ERROR] step1 未生成 parquet: {out_dir}")
                return False
            src = candidates[0]
            log(f"[count] {exp} step1 parquet 归位: {src.name} → {parquet.name}")
            src.rename(parquet)

    # step2: parquet -> on/off target_counts
    for kind in ("on", "off"):
        tgt_csv = cfg["targets"][kind]
        out_csv = out_dir / f"{kind}_target_counts.csv"
        if not force and out_csv.exists():
            log(f"[skip] {kind} counts 已存在: {out_csv}")
            continue
        cmd = [python, str(step2),
               "--reads-parquet", str(parquet),
               "--target-csv", tgt_csv,
               "--output", str(out_csv),
               "--window", str(window),
               "--mapq-threshold", str(cfg.get("mapq_threshold", 30))]
        if not run(cmd, desc=f"[count] {exp} step2 {kind}"):
            return False
    return True


def _load_counts(path: Path) -> dict:
    """读取 target_counts.csv（保留基因注释列）。"""
    import pandas as pd
    return pd.read_csv(path).set_index("target_id")


def _normalize_n_values(value):
    """把绘图 N 值（int / str / list）规范化为 int 列表。

    同一实验可配置多个 N（如 [10000, 80000]），绘图阶段会对每个 N 各生成一批图
    （输出目录 {frac}_N{n}/，互不覆盖）。
    """
    if value is None:
        return [80000]
    if isinstance(value, (list, tuple)):
        out = []
        for v in value:
            try:
                out.append(int(v))
            except (TypeError, ValueError):
                continue
        return out or [80000]
    try:
        return [int(value)]
    except (TypeError, ValueError):
        return [80000]


def _apply_diff_blacklist(df, blacklist: dict):
    """按 gene_name / target_id 黑名单过滤计数表。"""
    if not blacklist:
        return df

    target_ids = {str(x) for x in blacklist.get("target_ids", [])}
    gene_names = {str(x).strip() for x in blacklist.get("gene_names", [])}
    if not target_ids and not gene_names:
        return df

    mask = False
    if "gene_name" in df.columns:
        mask = df["gene_name"].astype(str).isin(gene_names)
    if target_ids:
        mask = mask | df.index.astype(str).isin(target_ids)

    return df.loc[~mask]


def _make_diff(exp_df, ctrl_df, out_path: Path, ascend: bool, blacklist: dict = None):
    """按 target_id 对齐实验组与对照组，计算 diff，排序并输出。

    在计算前先滤掉来自 chrM 的条目（线粒体基因组计数会对全图造成干扰）。
    并可根据黑名单移除指定 gene_name / target_id。
    """
    import pandas as pd
    # 滤掉 chrM（on/off counts 中都可能有）
    exp_df = exp_df[~exp_df["chromosome"].astype(str).str.startswith("chrM")]
    ctrl_df = ctrl_df[~ctrl_df["chromosome"].astype(str).str.startswith("chrM")]

    exp_df = _apply_diff_blacklist(exp_df, blacklist)
    ctrl_df = _apply_diff_blacklist(ctrl_df, blacklist)

    df = exp_df.join(ctrl_df["read_count"].rename("ctrl"), how="left")
    df["read_count_ctrl"] = df["ctrl"].fillna(0).astype(int)
    df["read_count_sample"] = exp_df["read_count"].astype(int)
    df["diff"] = df["read_count_sample"] - df["read_count_ctrl"]
    df = df.drop(columns=["ctrl"])

    out = df[["read_count_ctrl", "read_count_sample", "diff",
              "gene_name", "chromosome", "start_bp", "end_bp"]].copy()
    out["center"] = ((out["start_bp"] + out["end_bp"]) // 2)
    out = out.sort_values("diff", ascending=ascend)
    out.insert(0, "target_id", out.index)
    out.to_csv(out_path, index=False)
    return out


def step_diff(cfg: dict, exps: list, force: bool, apply_blacklist: bool = False) -> bool:
    diff_root = Path(cfg["workflow_root"]) / "4.Diff"
    diff_root.mkdir(parents=True, exist_ok=True)

    control = cfg["control_exp"]
    counts_base = Path(cfg["workflow_root"]) / "3.Counts"
    ctrl_paths = {
        "on": counts_base / control / "on_target_counts.csv",
        "off": counts_base / control / "off_target_counts.csv",
    }
    dry = bool(args_singleton and getattr(args_singleton, "dry_run", False))
    if not all(p.exists() for p in ctrl_paths.values()):
        if dry:
            log(f"[dry-run] 缺对照组 counts（{control}），将需要先跑 count 步骤")
        else:
            log(f"[ERROR] 缺少对照组 counts: {ctrl_paths}，请先对对照组跑 count 步骤")
            return False

    # 每种饱和度（frac）分别计算 diff；每个实验归档到 4.Diff/{exp}/ 子目录
    for exp, is_ctrl in exps:
        if is_ctrl:
            continue  # 对照组自身不算 diff
        diff_dir = diff_root / exp
        diff_dir.mkdir(parents=True, exist_ok=True)
        for frac in FRACTIONS:
            for kind, ascend in (("on", False), ("off", True)):
                exp_csv = counts_base / exp / f"{kind}_target_counts_{frac}.csv"
                ctrl_csv = counts_base / control / f"{kind}_target_counts_{frac}.csv"
                suffix = "" if frac == "1" else f"_{frac}"
                out_csv = diff_dir / f"{exp}_{kind}_diff{suffix}.csv"
                if not exp_csv.exists() or not ctrl_csv.exists():
                    if dry:
                        log(f"[dry-run] 将生成 {out_csv.name}（需 {exp}/{control} {kind} counts {frac}）")
                    else:
                        log(f"[WARN] 缺少 {kind} counts（frac={frac}）: {exp_csv} 或 {ctrl_csv}，跳过")
                    continue
                if not force and out_csv.exists():
                    log(f"[skip] diff 已存在: {out_csv}")
                    continue
                if dry:
                    log(f"[dry-run] 将生成 {out_csv.name}（diff: {exp} vs {control} @{frac}）")
                    continue
                exp_df = _load_counts(exp_csv)
                ctrl_df = _load_counts(ctrl_csv)
                blacklist = cfg.get("diff_blacklist") if apply_blacklist else None
                _make_diff(exp_df, ctrl_df, out_csv, ascend, blacklist)
                log(f"[diff] {exp} {kind} @{frac}: → {out_csv}")

            # 合并 on/off -> {exp}_diff{suffix}.csv
            merged_csv = diff_dir / f"{exp}_diff{suffix}.csv"
            on_csv = diff_dir / f"{exp}_on_diff{suffix}.csv"
            off_csv = diff_dir / f"{exp}_off_diff{suffix}.csv"
            if on_csv.exists() and off_csv.exists() and (force or not merged_csv.exists()):
                import pandas as pd
                df_on = pd.read_csv(on_csv)
                df_on["on_target"] = True
                df_off = pd.read_csv(off_csv)
                df_off["on_target"] = False
                df_merged = pd.concat([df_on, df_off], ignore_index=True)
                df_merged.to_csv(merged_csv, index=False)
                log(f"[diff] {exp} 合并 on/off @{frac} → {merged_csv}")

                # filtered 版本：删除 OFF 组中对照 read_count > 10 的“背景信号”行
                filtered_csv = diff_dir / f"{exp}_diff{suffix}_filtered.csv"
                df_filtered = df_merged[~((df_merged["on_target"] == False)
                                          & (df_merged["read_count_ctrl"] > 10))]
                df_filtered.to_csv(filtered_csv, index=False)
                log(f"[diff] {exp} filtered @{frac}: {len(df_merged)} → {len(df_filtered)} 行")
    return True


def step_plot(cfg: dict, exps: list, args) -> bool:
    """调用绘图脚本。输出结构：5.Plots/{scatter_slim|population}/{exp}/{frac}_N{n}/。

    直接覆盖同名文件（无时间戳），方便反复微调。N 值按实验分别取 args.n_rows_map。
    """
    plots_root = Path(cfg["workflow_root"]) / "5.Plots"
    plots_root.mkdir(parents=True, exist_ok=True)
    python = env_bin(cfg, "python")
    diff_root = Path(cfg["workflow_root"]) / "4.Diff"
    dates = [e for e, is_ctrl in exps if not is_ctrl]

    scatter_script = PIPELINE_ROOT / "plots" / "plot_by_n_s.py"
    pop_script = PIPELINE_ROOT / "plots" / "plot_population.py"

    for frac in FRACTIONS:
        # scatter_slim：每个实验单独调用，data-dir 指向 4.Diff/{exp}/
        # 同一实验可配置多个 N，plot_by_n_s 一次调用生成 {frac}_N{n} 多批图
        if scatter_script.exists():
            scatter_root = plots_root / "scatter_slim"
            for exp in dates:
                n_list = _normalize_n_values(args.n_rows_map.get(exp, 80000))
                cmd = [python, str(scatter_script),
                       "--n-rows"] + [str(n) for n in n_list] + [
                       "--data-dir", str(diff_root / exp),
                       "--output-root", str(scatter_root),
                       "--frac", frac,
                       "--y-cap", args.y_cap,
                       "--dates", exp]
                if args.exclude_chrm:
                    cmd.append("--exclude-chrm")
                if args.cap_auto:
                    cmd.append("--cap-auto")
                n_tag = "/".join(str(n) for n in n_list)
                if not run(cmd, desc=f"[plot] scatter_slim {exp} @{frac} (N={n_tag})"):
                    return False
        else:
            log("[warn] 跳过 scatter_slim：脚本缺失")

        # population：每个实验单独调用（脚本只接受单个 N，多 N 时逐个调用）
        if pop_script.exists() and not args.no_population:
            pop_root = plots_root / "population"
            for exp in dates:
                n_list = _normalize_n_values(args.n_rows_map.get(exp, 80000))
                for n_rows in n_list:
                    cmd = [python, str(pop_script),
                           "--data-dir", str(diff_root / exp),
                           "--output-root", str(pop_root),
                           "--n-rows", str(n_rows),
                           "--frac", frac,
                           "--dates", exp]
                    if args.exclude_chrm:
                        cmd.append("--exclude-chrm")
                    if not run(cmd, desc=f"[plot] population {exp} @{frac} (N={n_rows})"):
                        return False
    return True


# ---------------------------------------------------------------- main
def parse_args():
    p = argparse.ArgumentParser(description="WORF 统一 pipeline 编排（自包含）")
    p.add_argument("--config", default=None, help="config.json 路径")
    p.add_argument("--steps", default=",".join(STEPS_ORDER),
                   help=f"要执行的步骤，逗号分隔: {','.join(STEPS_ORDER)}（默认全部）")
    p.add_argument("--exp", nargs="+", help="指定实验（可多个）；默认用 config.experiments")
    p.add_argument("--all", action="store_true", help="扫描 raw_root 下全部实验目录")
    p.add_argument("--window", type=int, default=None, help="覆盖 step2 扩窗 bp")
    p.add_argument("--force", action="store_true", help="强制重跑（覆盖已存在输出）")
    p.add_argument("--n-rows", type=int, default=None,
                   help="绘图每 on/off 组取行数（默认用 config.n_rows，可被 --n-rows-per-exp 覆盖）")
    p.add_argument("--n-rows-per-exp", action="append", default=[], metavar="EXP:N[,N...]",
                   help="按实验指定绘图 N 值，可多次，支持多值，如 --n-rows-per-exp 20260508WORFT5-SNAP:10000,80000")
    p.add_argument("--exclude-chrm", action="store_true", help="绘图排除 chrM")
    p.add_argument("--cap-auto", action="store_true", help="绘图自动 outlier 裁剪")
    p.add_argument("--apply-diff-blacklist", action="store_true",
                   help="启用 config.json 中的 diff_blacklist，真正从 diff 输出中移除这些条目")
    p.add_argument("--y-cap", type=str, default="topk", choices=["topk", "iqr", "p99", "max"],
                   help="散点图稳健 y 上限策略（默认 topk = 只保留顶部少数点到微缩轴；iqr / p99 / max 也可选）")
    p.add_argument("--no-population", action="store_true", help="跳过群体统计图")
    p.add_argument("--keep-sam", action="store_true", help="保留中间 SAM 文件")
    p.add_argument("--dry-run", action="store_true", help="只打印将执行的命令，不实际执行")
    return p.parse_args()


def main():
    global args_singleton
    args = parse_args()
    args_singleton = args
    cfg = load_config(args.config)
    if args.window is not None:
        cfg["window"] = args.window

    exps = resolve_experiments(cfg, args)
    log(f"实验列表: {[f'{e}(ctrl)' if c else e for e, c in exps]}")

    # ---- 解析按实验的 N 值（仅 plot 用；前 4 步不涉及 N）----
    # 优先级: --n-rows-per-exp EXP:N[,N...] > --n-rows > config.n_rows[exp] > 80000
    # N 支持多个（如 10000,80000），同一实验会生成多批图
    per_exp_n = {}
    for item in args.n_rows_per_exp:
        if ":" not in item:
            log(f"[WARN] --n-rows-per-exp 格式应为 EXP:N 或 EXP:N1,N2，忽略: {item}")
            continue
        exp_key, n_str = item.split(":", 1)
        vals = []
        for part in n_str.split(","):
            part = part.strip()
            if not part:
                continue
            try:
                vals.append(int(part))
            except ValueError:
                log(f"[WARN] 无效 N 值: {part}")
        if vals:
            per_exp_n[exp_key.strip()] = vals
        else:
            log(f"[WARN] --n-rows-per-exp 无有效 N，忽略: {item}")

    cfg_n_rows = cfg.get("n_rows", {})
    global_n = args.n_rows if args.n_rows is not None else None
    args.n_rows_map = {}
    for exp, _ in exps:
        if exp in per_exp_n:
            n = per_exp_n[exp]
        elif global_n is not None:
            n = _normalize_n_values(global_n)
        else:
            n = _normalize_n_values(cfg_n_rows.get(exp, 80000))
        args.n_rows_map[exp] = n
    log("按实验绘图 N 值: " + ", ".join(f"{e}={n}" for e, n in args.n_rows_map.items()))

    steps = [s.strip() for s in args.steps.split(",") if s.strip()]
    unknown = set(steps) - set(STEPS_ORDER)
    if unknown:
        log(f"[ERROR] 未知步骤: {unknown}")
        sys.exit(1)

    # 逐实验执行 qc/align/count
    for exp, is_ctrl in exps:
        samples = find_samples(cfg, exp)
        if not samples:
            log(f"[WARN] {exp}: 未发现 fastq，跳过")
            continue
        for s in samples:
            if "qc" in steps and not step_qc(cfg, exp, s, args.force):
                sys.exit(1)
            if "align" in steps and not step_align(cfg, exp, s, args.force):
                sys.exit(1)
            if "count" in steps and not step_count(cfg, exp, s, args.force):
                sys.exit(1)

    if "diff" in steps:
        if not step_diff(cfg, exps, args.force, apply_blacklist=args.apply_diff_blacklist):
            sys.exit(1)

    if "plot" in steps:
        if not step_plot(cfg, exps, args):
            sys.exit(1)

    log("✅ 全部完成")


if __name__ == "__main__":
    main()
