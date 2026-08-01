#!/usr/bin/env python3
"""
Step 1: Export read genomic intervals from BAM to Parquet
Uses pysam for direct BAM reading (no samtools required)
"""
import argparse
import os
import re
import sys
from collections import defaultdict, OrderedDict
import json

import pandas as pd
import pysam
from datetime import datetime

def cigar_ref_length(cigar: str) -> int:
    """Calculate reference-consuming length from CIGAR string."""
    if not cigar or cigar == "*":
        return 0
    ref_ops = {"M", "D", "N", "=", "X"}
    total = 0
    for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        if op in ref_ops:
            total += int(length)
    return total


def sanitize_token(s: str) -> str:
    """Convert arbitrary text into a short filesystem-safe token."""
    s = (s or "").strip()
    s = re.sub(r"\.[^.]+$", "", s)  # remove extension
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s or "NA"


def extract_sample_name_from_bam(bam_path: str) -> str:
    """Extract sample name from path layout, fallback to BAM basename.

    Preferred: parent before "00.mergeRawFq" (project/sample-level identifier).
    Fallback: BAM basename leading token, e.g. UDI001_aligned.sorted.bam -> UDI001.
    """
    marker = "00.mergeRawFq"
    norm = bam_path.replace("\\", "/")
    parts = [p for p in norm.split("/") if p]
    if marker in parts:
        idx = parts.index(marker)
        if idx >= 1:
            return sanitize_token(parts[idx - 1])

    base = os.path.basename(bam_path)
    stem = os.path.splitext(base)[0]
    if "_" in stem:
        stem = stem.split("_", 1)[0]
    return sanitize_token(stem)


def extract_task_id_from_path(path: str, timestamp_source: str = "now") -> str:
    """Extract task_id using the rule: find '00.mergeRawFq' in the path,
        take its parent directory name and append a timestamp as YYYYMMDDHHMMSS.
    If the marker isn't found, fall back to filename (without extension).

        timestamp_source:
            - now (default): use current run time
            - parent_mtime: use the marker parent directory mtime (legacy behavior)
            - bam_mtime: use BAM file mtime
    """
    marker = "00.mergeRawFq"
    # normalize separators
    norm = path.replace("\\", "/")
    if marker not in norm:
        return os.path.splitext(os.path.basename(path))[0]

    parts = [p for p in norm.split("/") if p != ""]
    try:
        idx = parts.index(marker)
    except ValueError:
        return os.path.splitext(os.path.basename(path))[0]

    # parent directory name (the TASK_ID according to your rule)
    if idx >= 1:
        parent_name = parts[idx - 1]
        # reconstruct parent full path to get timestamp
        parent_path = "/".join(norm.split("/")[:idx])
        if parent_path == "":
            parent_path = os.path.dirname(path)
    else:
        parent_name = os.path.splitext(os.path.basename(path))[0]
        parent_path = os.path.dirname(path)

    # Choose timestamp source. Default is current run time to make each run identifiable.
    if timestamp_source == "parent_mtime":
        try:
            mtime = os.path.getmtime(parent_path)
            ts = datetime.fromtimestamp(mtime).strftime("%Y%m%d%H%M%S")
        except Exception:
            ts = datetime.now().strftime("%Y%m%d%H%M%S")
    elif timestamp_source == "bam_mtime":
        try:
            mtime = os.path.getmtime(path)
            ts = datetime.fromtimestamp(mtime).strftime("%Y%m%d%H%M%S")
        except Exception:
            ts = datetime.now().strftime("%Y%m%d%H%M%S")
    else:
        ts = datetime.now().strftime("%Y%m%d%H%M%S")

    return f"{parent_name}_{ts}"

def main():
    parser = argparse.ArgumentParser(
        description="Export read genomic intervals from BAM to Parquet (using pysam)"
    )
    parser.add_argument("--bam", required=True, help="Input BAM file path")
    parser.add_argument(
        "--output",
        default="reads_positions.parquet",
        help="Output Parquet file path (default: reads_positions.parquet)",
    )
    parser.add_argument(
        "--min-mapq",
        type=int,
        default=0,
        help="Minimum MAPQ to keep (default: 0)",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=1000000,
        help="Print progress every N processed records",
    )
    parser.add_argument(
        "--compression",
        default="snappy",
        choices=["snappy", "gzip", "brotli", "none"],
        help="Parquet compression codec (default: snappy)",
    )
    parser.add_argument(
        "--target-csv",
        default="",
        help="Target CSV path used to build task folder name",
    )
    parser.add_argument(
        "--window",
        default="NA",
        help="Window value used to build task folder name",
    )
    parser.add_argument(
        "--task-suffix",
        default="",
        help="Optional suffix appended to the inferred task_id to separate runs from the same BAM",
    )
    parser.add_argument(
        "--task-timestamp-source",
        default="now",
        choices=["now", "parent_mtime", "bam_mtime"],
        help="Timestamp source for inferred task_id (default: now)",
    )
    parser.add_argument(
        "--step2-mode",
        default="step2_1",
        choices=["step2", "step2_1", "step2.1"],
        help="Default step2 matcher route recorded in metadata (default: step2_1)",
    )
    parser.add_argument(
        "--bam-name",
        default="",
        help="BAM identifier used for task naming when --bam is stdin (-)",
    )
    args = parser.parse_args()

    is_stdin = (args.bam == "-")
    if not is_stdin and not os.path.exists(args.bam):
        print(f"[ERROR] BAM file not found: {args.bam}", file=sys.stderr)
        return 1

    step2_mode = "step2_1" if args.step2_mode == "step2.1" else args.step2_mode

    # New naming rule: sample_[TARGET_CSV]_[window]
    # Example: UDI001_unmarked_sites_annotated_with_gtf_plus_150
    # When reading from stdin, use --bam-name for sample name extraction
    bam_ref_for_naming = args.bam_name if (is_stdin and args.bam_name) else args.bam
    sample_name = extract_sample_name_from_bam(bam_ref_for_naming)
    target_csv_token = sanitize_token(os.path.basename(args.target_csv)) if args.target_csv else "TARGET_CSV"
    window_token = sanitize_token(str(args.window))
    task_id = f"{sample_name}_{target_csv_token}_{window_token}"
    print(f"[INFO] Detected task_id: {task_id}", file=sys.stderr)

    # Determine base output directory and final output paths.
    # If user passed a .parquet filename, keep that name inside the task dir;
    # otherwise treat args.output as a directory and use the default parquet name.
    arg_output = args.output
    if arg_output.lower().endswith('.parquet'):
        abs_output_path = os.path.abspath(arg_output)
        base_output_dir = os.path.dirname(abs_output_path) or os.getcwd()
        output_basename = os.path.basename(abs_output_path)
    else:
        base_output_dir = os.path.abspath(arg_output)
        output_basename = os.path.basename(args.output) if os.path.splitext(args.output)[1] else 'reads_positions.parquet'

    # create task-specific directory under base_output_dir
    task_dir = os.path.join(base_output_dir, task_id)
    os.makedirs(task_dir, exist_ok=True)

    final_output = os.path.join(task_dir, output_basename)
    # Ensure directory for the final output exists (redundant but safe)
    os.makedirs(os.path.dirname(final_output), exist_ok=True)

    processed = 0
    mapped_alignment_count = 0
    high_mapq_alignment_count = 0

    # Use pysam to read BAM directly; "-" reads BAM from stdin
    try:
        bam_file = pysam.AlignmentFile(args.bam, "rb")
    except Exception as e:
        print(f"[ERROR] Failed to open BAM file: {e}", file=sys.stderr)
        return 1

    # Use lists for faster appending than DataFrame
    rows = []
    buffer_rows = {}
    buffer_counts = defaultdict(int)
    last_flush_chrom = None
    last_flush_pos = -1

    for read in bam_file:
        processed += 1

        # Skip unmapped reads
        if read.is_unmapped:
            continue

        mapped_alignment_count += 1
        mapq = read.mapq
        if mapq > 30:
            high_mapq_alignment_count += 1
            
        if mapq < args.min_mapq:
            continue

        # Get read info
        read_name = read.query_name
        chromosome = bam_file.get_reference_name(read.reference_id)
        pos = read.reference_start  # 0-based
        cigar = read.cigarstring if read.cigarstring else ""
        
        is_reverse = read.is_reverse
        ref_len = cigar_ref_length(cigar)
        end_bp = pos + ref_len if ref_len > 0 else pos + 1
        strand = "-" if is_reverse else "+"

        row = [
            read_name,
            chromosome,
            pos,
            end_bp,
            mapq,
            strand,
            cigar,
        ]

        pos_key = (chromosome, strand, pos, end_bp)
        
        # Check if we should flush buffer
        if chromosome != last_flush_chrom or pos > last_flush_pos + 10000:
            keys_to_flush = []
            for bk in list(buffer_rows.keys()):
                if bk[0] != chromosome or pos > bk[2] + 10000:
                    b_row = buffer_rows[bk]
                    b_row.append(buffer_counts[bk])
                    rows.append(b_row)
                    keys_to_flush.append(bk)
                    
            for bk in keys_to_flush:
                del buffer_rows[bk]
                del buffer_counts[bk]
            
            last_flush_chrom = chromosome
            last_flush_pos = pos
        
        buffer_counts[pos_key] += 1
        if pos_key not in buffer_rows:
            buffer_rows[pos_key] = row

        if args.progress_every > 0 and processed % args.progress_every == 0:
            print(f"[INFO] Processed {processed:,} records", file=sys.stderr)

    # Flush any remaining rows
    for bk in buffer_rows.keys():
        b_row = buffer_rows[bk]
        b_row.append(buffer_counts[bk])
        rows.append(b_row)
    buffer_rows.clear()
    buffer_counts.clear()

    bam_file.close()

    # Create DataFrame and save as Parquet
    columns = [
        "read_name",
        "read_chromosome",
        "read_start_bp",
        "read_end_bp",
        "read_mapq",
        "read_strand",
        "read_cigar",
        "duplicate_count",
    ]
    df = pd.DataFrame(rows, columns=columns)

    # attach task_id column
    try:
        df['task_id'] = task_id
    except Exception:
        df['task_id'] = None

    # Determine compression
    compression = args.compression if args.compression != "none" else None
    
    df.to_parquet(final_output, compression=compression, index=False)
    print(f"[SUCCESS] Export completed into {final_output}")
    print(f"[INFO] Total rows: {len(df):,}")

    # Save summary statistics as separate metadata file
    summary_path = final_output.replace(".parquet", "_summary.txt")
    map_ratio = (mapped_alignment_count / processed) * 100 if processed > 0 else 0.0
    high_mapq_ratio = (high_mapq_alignment_count / mapped_alignment_count) * 100 if mapped_alignment_count > 0 else 0.0
    
    with open(summary_path, "w") as f:
        f.write(f"Task_ID\t{task_id}\n")
        f.write(f"Total_Processed_Reads\t{processed}\n")
        f.write(f"Total_Mapped_Reads\t{mapped_alignment_count}\n")
        f.write(f"High_MAPQ_Reads_(>30)\t{high_mapq_alignment_count}\n")
        f.write(f"Mapping_Ratio\t{map_ratio:.2f}%\n")
        f.write(f"High_MAPQ_Ratio_(vs_Mapped)\t{high_mapq_ratio:.2f}%\n")

    # Create JSON metadata with task_id as the first key
    json_path = final_output.replace('.parquet', '_step1_reads_summary.json')
    metadata = OrderedDict()
    metadata[task_id] = {
        "bam": args.bam,
        "output_parquet": final_output,
        "params": {
            "step1": {
                "bam": args.bam,
                "output": args.output,
                "min_mapq": args.min_mapq,
                "progress_every": args.progress_every,
                "compression": args.compression,
                "target_csv": args.target_csv,
                "window": args.window,
                "step2_mode": step2_mode,
                "resolved_compression": compression,
                "base_output_dir": base_output_dir,
                "task_dir": task_dir,
                "final_output": final_output,
                "argv": sys.argv,
            }
        },
        "stats": {
            "total_processed_reads": processed,
            "total_mapped_reads": mapped_alignment_count,
            "high_mapq_reads_gt30": high_mapq_alignment_count,
            "mapping_ratio_percent": round(map_ratio, 2),
            "high_mapq_ratio_percent": round(high_mapq_ratio, 2),
        },
        "created_at": datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),
    }

    try:
        with open(json_path, 'w', encoding='utf-8') as jf:
            json.dump(metadata, jf, ensure_ascii=False, indent=2)
        print(f"[INFO] JSON metadata saved to {json_path}")
    except Exception as e:
        print(f"[WARN] Failed to write JSON metadata: {e}", file=sys.stderr)
    print(f"[INFO] Summary statistics saved to {summary_path}")

if __name__ == "__main__":
    sys.exit(main())