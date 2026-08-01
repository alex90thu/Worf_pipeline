#!/usr/bin/env python3
"""
Step 2.1: Match extracted reads (Parquet) against targets with containment rules.

Rules (on already chromosome-matched candidates):
1) If read length > target length: match iff read covers target.
   a <= target_start and b >= target_end
2) If read length < target length: match iff read is inside target.
   target_start <= a and b <= target_end
3) If read length == target length: not matched (strict > / < only).
"""
import argparse
import csv
import sys
import os
import json
from typing import Sequence
from collections import defaultdict
from intervaltree import IntervalTree
from datetime import datetime

import pandas as pd


def parse_int_field(val: str | None) -> int | None:
    if not val:
        return None
    try:
        return int(float(str(val).strip().replace(",", "")))
    except ValueError:
        return None


def find_column(fieldnames: Sequence[str], candidates: Sequence[str]) -> str | None:
    lowered = {f.lower().strip(): f for f in fieldnames}
    for c in candidates:
        key = c.lower().strip()
        if key in lowered:
            return lowered[key]
    return None


def load_targets(target_csv: str, window: int):
    """Load targets and build interval trees for candidate lookup."""
    targets = []
    interval_trees = defaultdict(IntervalTree)

    with open(target_csv, "r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            return targets, interval_trees

        fieldnames = reader.fieldnames
        chrom_col = find_column(fieldnames, ["chromosome", "chrom", "chr", "染色体"])
        start_col = find_column(fieldnames, ["start", "start_bp", "pos_start", "起始", "起点", "开始"])
        end_col = find_column(fieldnames, ["end", "end_bp", "pos_end", "终止", "终点", "结束"])
        pos_col = find_column(fieldnames, ["position", "pos", "site", " 位点", "坐标", "position (bp)"])

        if not chrom_col:
            raise ValueError("target CSV missing chromosome column")
        if not ((start_col and end_col) or pos_col):
            raise ValueError("target CSV needs start+end columns, or a position column")

        gene_col = find_column(fieldnames, ["gene_name", "gene", "name", "symbol", "target", "靶点", "基因"])

        for idx, row in enumerate(reader, start=1):
            chrom_raw = (row.get(chrom_col) or "").strip()
            if not chrom_raw:
                continue

            if start_col and end_col:
                start_bp = parse_int_field(row.get(start_col))
                end_bp = parse_int_field(row.get(end_col))
                if start_bp == end_bp and start_bp is not None:
                    start_bp = max(0, start_bp - window)
                    end_bp = end_bp + window
            elif pos_col:
                pos_bp = parse_int_field(row.get(pos_col))
                if pos_bp is not None:
                    start_bp = max(0, pos_bp - window)
                    end_bp = pos_bp + window
                else:
                    start_bp = None
                    end_bp = None
            else:
                start_bp = None
                end_bp = None

            if start_bp is None or end_bp is None:
                continue

            if start_bp > end_bp:
                start_bp, end_bp = end_bp, start_bp

            gene_name = (row.get(gene_col) or "NA").strip() if gene_col else "NA"

            t = {
                "target_id": idx,
                "gene_name": gene_name,
                "chromosome": chrom_raw,
                "start_bp": start_bp,
                "end_bp": end_bp,
                "target_len": end_bp - start_bp,
                "read_count": 0,
            }

            targets.append(t)
            interval_trees[chrom_raw].addi(start_bp, end_bp, idx - 1)

    return targets, interval_trees


def is_step21_match(read_start: int, read_end: int, target_start: int, target_end: int) -> bool:
    read_len = read_end - read_start
    target_len = target_end - target_start

    if read_len > target_len:
        return read_start <= target_start and read_end >= target_end

    if read_len < target_len:
        return target_start <= read_start and read_end <= target_end

    # Strict rules: length equality is not matched.
    return False


def main():
    parser = argparse.ArgumentParser(description="Step2.1 containment matcher for extracted reads (Parquet)")
    parser.add_argument("--reads-parquet", required=True, help="Input reads Parquet file (from step1)")
    parser.add_argument("--target-csv", required=True, help="Target CSV file")
    parser.add_argument("--output", required=True, help="Output hit counts CSV")
    parser.add_argument("--window", type=int, default=500, help="Flanking window around single-point targets (default=500)")
    parser.add_argument("--mapq-threshold", type=int, default=30, help="Minimum MAPQ threshold (default=30)")
    args = parser.parse_args()

    print(f"[INFO] Loading targets from {args.target_csv}...")
    try:
        targets, interval_trees = load_targets(args.target_csv, args.window)
    except Exception as e:
        print(f"[ERROR] Failed to load target CSV: {e}", file=sys.stderr)
        return 1

    if not targets:
        print("[ERROR] No valid targets found in target CSV", file=sys.stderr)
        return 1

    print(f"[INFO] Loaded {len(targets)} targets with +/- {args.window}bp window.")
    chroms_with_targets = set(interval_trees.keys())
    print(f"[INFO] Target chromosomes: {len(chroms_with_targets)}")

    print(f"[INFO] Loading reads from {args.reads_parquet}...")

    target_counts = [0] * len(targets)

    total_reads = 0
    matched_reads = 0

    df = pd.read_parquet(
        args.reads_parquet,
        columns=["read_chromosome", "read_start_bp", "read_end_bp", "read_mapq", "duplicate_count"],
    )
    print(f"[INFO] Total reads in parquet: {len(df):,}")

    for chrom_name, chrom_df in df.groupby("read_chromosome", observed=True):
        if chrom_name not in interval_trees:
            continue

        chrom_df = chrom_df[chrom_df["read_mapq"] > args.mapq_threshold]
        if len(chrom_df) == 0:
            continue

        tree = interval_trees[chrom_name]

        for _, row in chrom_df.iterrows():
            total_reads += 1
            start = int(row["read_start_bp"])
            end = int(row["read_end_bp"])
            dup_count = int(row.get("duplicate_count", 1))
            if pd.isna(dup_count) or dup_count < 1:
                dup_count = 1

            candidates = tree.overlap(start, end)
            any_hit = False
            for interval in candidates:
                target = targets[interval.data]
                if is_step21_match(start, end, target["start_bp"], target["end_bp"]):
                    target_counts[interval.data] += dup_count
                    any_hit = True

            if any_hit:
                matched_reads += 1

            if total_reads % 5000000 == 0:
                print(f"[INFO] Processed {total_reads:,} reads, {matched_reads:,} matched...")

    for i, t in enumerate(targets):
        t["read_count"] = target_counts[i]

    print(f"[INFO] Writing {len(targets)} output rows to {args.output}")
    with open(args.output, "w", newline="", encoding="utf-8") as f_out:
        writer = csv.writer(f_out)
        writer.writerow(["target_id", "gene_name", "chromosome", "start_bp", "end_bp", "read_count"])
        for t in targets:
            writer.writerow([
                t["target_id"],
                t["gene_name"],
                t["chromosome"],
                t["start_bp"],
                t["end_bp"],
                t["read_count"],
            ])

    print(f"[SUCCESS] Wrote {len(targets)} records to {args.output}")
    print(f"[SUMMARY] Total reads processed: {total_reads:,}, Matched: {matched_reads:,}")

    try:
        json_path = args.reads_parquet.replace('.parquet', '_step1_reads_summary.json')
        if os.path.exists(json_path):
            with open(json_path, 'r', encoding='utf-8') as jf:
                payload = json.load(jf)
        else:
            payload = {}

        if payload:
            task_key = next(iter(payload.keys()))
            record = payload.get(task_key, {})
        else:
            task_key = os.path.splitext(os.path.basename(args.reads_parquet))[0]
            record = {}

        params = record.setdefault('params', {})
        params['step2'] = {
            'matcher': 'step2.1',
            'reads_parquet': args.reads_parquet,
            'target_csv': args.target_csv,
            'output': args.output,
            'window': args.window,
            'mapq_threshold': args.mapq_threshold,
            'argv': sys.argv,
        }

        artifacts = record.setdefault('artifacts', {})
        artifacts['step2_output'] = args.output
        record['updated_at'] = datetime.now().strftime('%Y-%m-%dT%H:%M:%S')
        payload[task_key] = record

        with open(json_path, 'w', encoding='utf-8') as jf:
            json.dump(payload, jf, ensure_ascii=False, indent=2)
        print(f"[INFO] Updated metadata JSON: {json_path}")
    except Exception as e:
        print(f"[WARN] Failed to update metadata JSON in step2.1: {e}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
