#!/usr/bin/env python3
"""
生成 WORF pipeline 文件索引 JSON（file_index.json）

作用：
- 记录每个样品在每个步骤中的文件是否完整（存在 / 大小 / 生成时间）
- 记录每个步骤的原文件来源目录（source_dir）
- 记录时间戳，可作为后续脚本断点续传的参考

用法：
    python build_index.py            # 生成到 ./file_index.json
    python build_index.py -o x.json  # 指定输出路径
"""
import argparse
import json
import os
from datetime import datetime
from pathlib import Path

PIPELINE_ROOT = Path(__file__).resolve().parent
CONFIG_PATH = PIPELINE_ROOT / "config.json"
DEFAULT_OUT = PIPELINE_ROOT / "file_index.json"

# 旧流程中 counts 的源目录（worf_benchmark saturation 全量 frac_1）
SAT_RUNS_ROOT = Path("/data/lulab_commonspace/guozehua/worf_benchmark/saturation_runs")

# 7 种饱和度（深度比例），与 run_pipeline.py 的 FRACTIONS 保持一致
FRACTIONS = ["0.001", "0.01", "0.05", "0.1", "0.2", "0.5", "1"]


def load_config():
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        cfg = json.load(f)
    for key in ("targets",):
        for k, v in cfg[key].items():
            if not os.path.isabs(v):
                cfg[key][k] = str((PIPELINE_ROOT / v).resolve())
    return cfg


def file_status(path: Path):
    """返回文件状态：exists / size_bytes / mtime（ISO 时间）"""
    if not path.exists():
        return {"exists": False, "size_bytes": None, "mtime": None}
    st = path.stat()
    return {
        "exists": True,
        "size_bytes": st.st_size,
        "mtime": datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%dT%H:%M:%S"),
    }


def step_entry(source_dir, file_paths, check_dir=False):
    """构造一个步骤条目。

    source_dir: 原文件地址目录
    file_paths: 该步骤在 workflow 下的目标文件路径列表
    """
    files = {str(p): file_status(p) for p in file_paths}
    # 判断 complete：必须有声明的文件，且全部存在
    complete = bool(files) and all(v["exists"] for v in files.values())
    # 步骤时间：取存在文件的最早 mtime；没有则用当前时间
    mtimes = [v["mtime"] for v in files.values() if v["mtime"]]
    timestamp = min(mtimes) if mtimes else datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    return {
        "source_dir": str(source_dir),
        "timestamp": timestamp,
        "complete": complete,
        "files": files,
    }


def resolve_raw_dir(cfg, exp):
    """返回该实验的 00.mergeRawFq / 00.mergerRawFq 目录（任一存在）。"""
    raw_root = Path(cfg["raw_root"])
    for name in ("00.mergeRawFq", "00.mergerRawFq"):
        cand = raw_root / exp / name
        if cand.is_dir():
            return cand
    return raw_root / exp


def build_index(cfg):
    wf = Path(cfg["workflow_root"])
    control = cfg["control_exp"]
    exps = list(cfg["experiments"]) + [control]

    experiments = {}
    for exp in exps:
        raw_dir = resolve_raw_dir(cfg, exp)
        is_ctrl = (exp == control)

        # ---- 扫描样品 ----
        samples = {}
        if raw_dir.is_dir():
            for sample_dir in sorted(raw_dir.iterdir()):
                if not sample_dir.is_dir():
                    continue
                sname = sample_dir.name
                src_res = sample_dir / "results"

                # qc: 源 = results/01_qc，目标 = 1.QC/{exp}/
                qc_src = src_res / "01_qc"
                qc_files = [
                    wf / "1.QC" / exp / f"{sname}_clean_1.fq.gz",
                    wf / "1.QC" / exp / f"{sname}_clean_2.fq.gz",
                ]
                # align: 源 = results/03_bam，目标 = 2.Alignment/{exp}/
                align_src = src_res / "03_bam"
                align_files = [
                    wf / "2.Alignment" / exp / f"{sname}_aligned.sorted.bam",
                    wf / "2.Alignment" / exp / f"{sname}_aligned.sorted.bam.bai",
                ]
                # count: 源 = saturation_runs/{exp}_{kind}/frac_1/
                count_src = SAT_RUNS_ROOT / f"{exp}_on" / "frac_1"
                count_files = [
                    wf / "3.Counts" / exp / f"{sname}_reads.parquet",
                    wf / "3.Counts" / exp / "on_target_counts.csv",
                    wf / "3.Counts" / exp / "off_target_counts.csv",
                ]

                samples[sname] = {
                    "raw_fq_dir": str(sample_dir),
                    "steps": {
                        "qc": step_entry(qc_src, qc_files),
                        "align": step_entry(align_src, align_files),
                        "count": step_entry(count_src, count_files),
                    },
                }

        # ---- 实验级步骤（diff / plot，不按样品；按 7 种饱和度 frac 分别记录）----
        # diff 归档到 4.Diff/{exp}/ 子目录
        diff_src = wf / "4.Diff" / exp
        diff_files = []
        for frac in FRACTIONS:
            suffix = "" if frac == "1" else f"_{frac}"
            diff_files += [
                wf / "4.Diff" / exp / f"{exp}_on_diff{suffix}.csv",
                wf / "4.Diff" / exp / f"{exp}_off_diff{suffix}.csv",
                wf / "4.Diff" / exp / f"{exp}_diff{suffix}.csv",
                wf / "4.Diff" / exp / f"{exp}_diff{suffix}_filtered.csv",
            ]
        # plot: 输出在 5.Plots/{scatter_slim,population}/{exp}/{frac}/（无时间戳，直接覆盖）
        plot_src = wf / "5.Plots"
        plot_files = []
        plots_root = wf / "5.Plots"
        if plots_root.is_dir():
            for category in ("scatter_slim", "population"):
                cat_dir = plots_root / category / exp
                if cat_dir.is_dir():
                    for frac_dir in sorted(cat_dir.iterdir()):
                        for f in sorted(frac_dir.rglob("*.png")):
                            plot_files.append(f)

        exp_steps = {"diff": step_entry(diff_src, diff_files)}
        # plot 若无输出文件，也记录（complete=false）
        exp_steps["plot"] = step_entry(plot_src, plot_files)

        experiments[exp] = {
            "is_control": is_ctrl,
            "raw_dir": str(raw_dir),
            "samples": samples,
            "steps": exp_steps,
        }

    return {
        "generated_at": datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),
        "workflow_root": str(wf),
        "pipeline_root": str(PIPELINE_ROOT),
        "window": cfg.get("window", 150),
        "mapq_threshold": cfg.get("mapq_threshold", 30),
        "experiments": experiments,
    }


def main():
    p = argparse.ArgumentParser(description="生成 WORF pipeline 文件索引 JSON")
    p.add_argument("-o", "--output", default=str(DEFAULT_OUT), help="输出 JSON 路径")
    p.add_argument("--config", default=None, help="config.json 路径")
    args = p.parse_args()

    cfg = load_config() if not args.config else json.load(open(args.config, encoding="utf-8"))
    index = build_index(cfg)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:
        json.dump(index, f, ensure_ascii=False, indent=2)
    print(f"✅ 索引已写入: {out}")
    # 简单摘要
    for exp, e in index["experiments"].items():
        n_ok = 0
        n_steps = 0
        for s in e["samples"].values():
            for st in s["steps"].values():
                n_steps += 1
                if st["complete"]:
                    n_ok += 1
        for st in e["steps"].values():
            n_steps += 1
            if st["complete"]:
                n_ok += 1
        print(f"  {exp}: 步骤完成 {n_ok}/{n_steps} (ctrl={e['is_control']})")


if __name__ == "__main__":
    main()
