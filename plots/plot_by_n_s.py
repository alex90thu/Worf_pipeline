#!/usr/bin/env python3
"""
plot_by_n 的窄版变体（s = slim），适合 A4 纸排版：
- 图宽 18 cm（~7.09 英寸）
- 顶刊标准无背景交替配色：深蓝/浅蓝 (on-target)、深红/浅红 (off-target)，纯白底
- X 轴标签去除 chr 前缀，直接展示 1, 2, ... X，彻底解决拥挤问题
"""
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
import argparse

# ---- matplotlib 全局设置 ----
mpl.rcParams['font.family'] = 'Arial'

# 图宽 18 cm（A4 宽 21 cm，留边距后正文约 18 cm）
CM_TO_INCH = 1 / 2.54
FIG_WIDTH_CM = 18
FIG_WIDTH_INCH = FIG_WIDTH_CM * CM_TO_INCH  # ≈ 7.09"

chrom_order = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY', 'chrM']
chrom_order_no_m = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']


def chrom_center_to_x(df, include_chrM=True):
    chrom_col = 'chromosome_x' if 'chromosome_x' in df.columns else 'chromosome'
    center_col = 'center'

    chroms = chrom_order if include_chrM else chrom_order_no_m
    df = df[df[chrom_col].isin(chroms)].copy()

    df['chrom_num'] = df[chrom_col].str.replace('chr', '').replace({'X': '23', 'Y': '24', 'M': '25'})
    df['chrom_num'] = pd.to_numeric(df['chrom_num'], errors='coerce')
    df = df.dropna(subset=['chrom_num'])
    df = df.sort_values(by=['chrom_num', center_col]).reset_index(drop=True)

    chrom_lengths = df.groupby('chrom_num', observed=False)[center_col].max().fillna(0)
    chrom_offsets = chrom_lengths.shift(1).fillna(0).cumsum()
    df['offset'] = df['chrom_num'].map(chrom_offsets)
    df['absolute_bp'] = df[center_col] + df['offset']

    return df, chrom_offsets, chrom_lengths


def auto_select_scale(df, y_col):
    sample_data = df[y_col].dropna()
    if len(sample_data) == 0:
        return True

    min_val = sample_data.min()
    max_val = sample_data.max()

    if max_val > 0 and min_val > 0:
        log_range = np.log10(max_val) - np.log10(min_val)
        if log_range > 2:
            return True

    return False


def auto_select_scale_values(values):
    """基于一维 values（已 dropna）判断是否用 log scale。"""
    values = np.asarray(values)
    if len(values) == 0:
        return True
    min_val = values.min()
    max_val = values.max()
    if max_val > 0 and min_val > 0:
        log_range = np.log10(max_val) - np.log10(min_val)
        if log_range > 2:
            return True
    return False


def write_max50_log(output_dir, label, frac, n, df, ylim_shared, scale_str):
    """在绘图子文件夹写 max50.log，记录该组别 top 50 的最大 read counts。

    作用：当散点太多看不清零星大点时，用此文件诊断 y 轴 scale 是否正确。
    - 分别列出 ctrl 和 sample 各自的 top 50 read_count（原始值，非 +1）
    - 附带本组绘图实际使用的 ylim（y_max = max*1.1），方便对照
    """
    lines = [
        f"# max50.log — top 50 read counts diagnostic",
        f"# experiment: {label}",
        f"# depth(frac): {frac}",
        f"# N: {n}",
        f"# scale: {scale_str}",
        f"# ylim used for BOTH ctrl & sample: y_min={ylim_shared[0]:.1f}, y_max={ylim_shared[1]:.1f} (max*1.1)",
        "",
    ]
    for col in ("read_count_ctrl", "read_count_sample"):
        vals = df[col].dropna().astype(int).sort_values(ascending=False).head(50)
        lines.append(f"[{col}] top 50 (max={vals.max() if len(vals) else 0}):")
        for i, v in enumerate(vals, start=1):
            lines.append(f"  {i:2d}. {v}")
        lines.append("")

    log_path = os.path.join(output_dir, "max50.log")
    with open(log_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"  max50.log 已写入: {log_path}")


def auto_cap(values):
    """基于 min / median / mean 自动估算 ctrl 背景 outlier 裁剪阈值。

    公式：cap = median + factor * (median - min)
    其中 factor = max(2, 6 / skew)，skew = mean / median。
    分布越对称 → factor 越大 → cap 越宽松；右尾越重 → cap 越紧。
    median==0 时回退为 3*mean。
    """
    vmin = np.min(values)
    vmed = np.median(values)
    vavg = np.mean(values)

    if vmed == 0:
        # 中位数为零（如 sample 大半为零）：以均值为中心，3 倍均值截顶
        return max(3.0 * vavg, 1.0)

    spread = vmed - vmin          # 主体数据半宽
    skew = vavg / vmed            # >1 右偏，≈1 对称
    factor = max(2.0, 6.0 / skew) # 对称 → factor≈6，严重右偏 → factor→2
    cap = vmed + factor * spread
    return max(cap, 1.0)


def compute_robust_ylim(df, cap_mode="iqr"):
    """基于 ctrl+sample 合并的 read_count 计算稳健 y 上限，并返回统计量。

    参数:
        df: 含 read_count_ctrl / read_count_sample 列的 DataFrame（已按需过滤 chrM）
        cap_mode: y 上限策略
            - "iqr": y_max = Q3 + 1.5*IQR（经典稳健法，推荐，容忍零星大点）
            - "p99": y_max = P99
            - "max": y_max = max（旧行为，backward compat）

    返回:
        (ylim, stats): ylim=(y_min, y_max)；stats 为 dict，含 mean/median/Q1/Q3/
                       P99/abs_max/top50 等，供写 max50.log 诊断。
    注意: 返回的 y_max 不 +1（+1 在绘图内部做）；此处统计基于原始 read_count。
    """
    combined = pd.concat([
        df['read_count_ctrl'].dropna(),
        df['read_count_sample'].dropna(),
    ]).astype(int)

    if len(combined) == 0:
        return (0.0, 1.0), {}

    q1, q3 = combined.quantile([0.25, 0.75])
    iqr = q3 - q1
    p99 = combined.quantile(0.99)
    abs_max = combined.max()
    n = len(combined)

    if cap_mode == "iqr":
        y_max_raw = q3 + 1.5 * iqr
        # 兜底：若 IQR 法上限过小（如分布集中在低值区），至少覆盖到 P99
        y_max_raw = max(y_max_raw, p99)
        y_max_raw = max(y_max_raw, 1.0)
    elif cap_mode == "p99":
        y_max_raw = max(p99, 1.0)
    else:  # "max"
        y_max_raw = max(abs_max, 1.0)

    y_max = y_max_raw + 1  # 与绘图内部 +1 保持一致
    data_min_p1 = combined.min() + 1
    y_min = max(0.0, data_min_p1 * 0.9) if data_min_p1 > 0 else 0.0

    stats = {
        "n": n,
        "mean": round(float(combined.mean()), 2),
        "median": int(combined.median()),
        "q1": int(q1),
        "q3": int(q3),
        "iqr": int(iqr),
        "p99": int(p99),
        "abs_max": int(abs_max),
        "y_max_raw": y_max_raw,
        "top50_ctrl": df['read_count_ctrl'].dropna().astype(int)
                      .sort_values(ascending=False).head(50).tolist(),
        "top50_sample": df['read_count_sample'].dropna().astype(int)
                        .sort_values(ascending=False).head(50).tolist(),
    }
    return (y_min, y_max), stats


def write_max50_log(output_dir, label, frac, n, stats, ylim_shared, scale_str):
    """在绘图子文件夹写 max50.log，记录该组别完整统计量（供诊断 y 轴 scale）。

    当散点太多看不清零星大点时，用此文件判断 y 轴 scale 是否合理。
    - 列出 ctrl/sample 各自的 top 50 read_count（原始值，非 +1）
    - 列出 mean/median/Q1/Q3/IQR/P99/abs_max 与绘图实际使用的 ylim
    """
    lines = [
        f"# max50.log — read count distribution diagnostic",
        f"# experiment: {label}",
        f"# depth(frac): {frac}",
        f"# N: {n}",
        f"# scale: {scale_str}",
        f"# ylim used for BOTH ctrl & sample: y_min={ylim_shared[0]:.1f}, y_max={ylim_shared[1]:.1f}",
        "",
    ]
    if stats:
        lines += [
            f"[overall] n={stats['n']}, mean={stats['mean']}, median={stats['median']}, "
            f"Q1={stats['q1']}, Q3={stats['q3']}, IQR={stats['iqr']}",
            f"[overall] P99={stats['p99']}, abs_max={stats['abs_max']}",
            f"[overall] y_max_raw (before +1) = {stats['y_max_raw']}",
            "",
        ]
        for col_key, col_name in (("top50_ctrl", "read_count_ctrl"),
                                  ("top50_sample", "read_count_sample")):
            vals = stats[col_key]
            lines.append(f"[{col_name}] top 50 (max={vals[0] if vals else 0}):")
            for i, v in enumerate(vals, start=1):
                lines.append(f"  {i:2d}. {v}")
            lines.append("")

    log_path = os.path.join(output_dir, "max50.log")
    with open(log_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print(f"  max50.log 已写入: {log_path}")


def plot_scatter_slim(df, y_col, output_path, title, use_log=False, include_chrM=True, ylim=None):
    """绘制窄版散点图：18 cm 宽，高级全染色体标注与交替背景。

    一律使用线性 y 轴（不接受 log scale）。
    """
    df, chrom_offsets, chrom_lengths = chrom_center_to_x(df.copy(), include_chrM)

    df_on = df[df['on_target'] == True]
    df_off = df[df['on_target'] == False]

    # ---- 18 cm 宽的窄图 ----
    fig_height = FIG_WIDTH_INCH * 0.42  # 高宽比约 0.42
    fig, ax = plt.subplots(figsize=(FIG_WIDTH_INCH, fig_height), dpi=300)

    y_on = df_on[y_col] + 1
    y_off = df_off[y_col] + 1

    if ylim:
        ax.set_ylim(ylim[0], ylim[1])

    # ---- 颜色定义 (顶刊标准：纯净白底 + 染色体交替配色) ----
    ON_COLORS  = ['#2874A6', '#5DADE2']   # 深蓝 / 浅蓝 交替
    OFF_COLORS = ['#B03A2E', '#EC7063']   # 深红 / 浅红 交替

    # ---- 生成刻度标签与位置 (不再绘制 axvspan 背景带) ----
    tick_labels = []
    tick_positions = []
    chrom_centers = chrom_offsets + (chrom_lengths / 2)
    chrom_nums = sorted([c for c in chrom_lengths.index if chrom_lengths[c] > 0])

    for cnum in chrom_nums:
        # 收集刻度中心位置
        tick_positions.append(chrom_centers[cnum])
        # 生成精简版标签 (去除 'chr')
        if cnum == 23:
            tick_labels.append('X')
        elif cnum == 24:
            tick_labels.append('Y')
        elif cnum == 25:
            tick_labels.append('M')
        else:
            tick_labels.append(str(int(cnum)))

    # ---- 按染色体交替绘制散点 (替代原 axvspan 背景带) ----
    for i, cnum in enumerate(chrom_nums):
        on_color = ON_COLORS[i % 2]
        off_color = OFF_COLORS[i % 2]

        # 当前染色体的 on_target 点
        mask_on = (df_on['chrom_num'] == cnum)
        if mask_on.any():
            ax.scatter(df_on.loc[mask_on, 'absolute_bp'], y_on[mask_on],
                       c=on_color, alpha=0.35, s=2, zorder=3, edgecolors='none')

        # 当前染色体的 off_target 点
        mask_off = (df_off['chrom_num'] == cnum)
        if mask_off.any():
            ax.scatter(df_off.loc[mask_off, 'absolute_bp'], y_off[mask_off],
                       c=off_color, alpha=0.35, s=2, zorder=2, edgecolors='none')

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, fontsize=7)
    
    # ---- 排版收尾 ----
    ax.set_xlim(0, chrom_offsets.iloc[-1] + chrom_lengths.iloc[-1])
    ax.set_xlabel('Chromosomes', fontsize=9, fontweight='bold')
    ax.set_ylabel('Read counts', fontsize=9, fontweight='bold')
    ax.set_title(title, fontsize=10, pad=10)
    
    # 极简边框处理
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # 自定义图例（代理艺术家，避免多染色体重复标签）
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=ON_COLORS[0],
               markersize=7, label=f'on (n={len(df_on)})'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=OFF_COLORS[0],
               markersize=7, label=f'off (n={len(df_off)})'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=6, frameon=False)
    
    plt.tight_layout(pad=0.5)
    plt.savefig(output_path, dpi=300, bbox_inches='tight', transparent=False)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='plot_by_n_s: slim scatter plots (18 cm wide, aesthetic axis)')
    parser.add_argument('--n-rows', nargs='+', type=int, default=[500, 1000, 5000, 10000, 20000, 30000, 40000])
    parser.add_argument('--exclude-chrm', action='store_true', help='Exclude chrM')
    parser.add_argument('--linear', action='store_true', help='Use linear scale (default is log)')
    parser.add_argument('--filtered', action='store_true', help='Use filtered diff files (_diff_filtered.csv)')
    parser.add_argument('--data-dir', type=str, default=None, help='Data directory')
    parser.add_argument('--csv', type=str, default=None, help='Direct CSV file path (bypasses --dates and --data-dir)')
    parser.add_argument('--cap-pct', type=float, default=None, metavar='PCT',
                        help='Cap read_count_ctrl/sample at given percentile (e.g. 99) to suppress outliers')
    parser.add_argument('--cap-auto', action='store_true',
                        help='Auto-estimate cap from min/median/mean and filter outliers above it')
    parser.add_argument('--dates', nargs='+', type=str, default=['20260404', '20260405', '20260407'],
                        help='Date prefixes (default: 20260404 20260405 20260407)')
    parser.add_argument('--frac', type=str, default='1',
                        help='Saturation depth fraction (default: 1 = full depth; '
                             'diff file suffix _<frac>, frac=1 uses no suffix)')
    parser.add_argument('--y-cap', type=str, default='iqr', choices=['iqr', 'p99', 'max'],
                        help='Robust y-max strategy (default: iqr = Q3+1.5*IQR, tolerate sparse outliers; '
                             'p99 = P99; max = old max*1.1 behavior). Data is never removed, only y-axis clipped.')
    parser.add_argument('--output-root', type=str, default=None,
                        help='Override output root dir (default: {WORKSPACE}/output_scatter_slim)')
    args = parser.parse_args()

    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

    if args.data_dir:
        WORKSPACE = args.data_dir
    else:
        WORKSPACE = os.path.dirname(SCRIPT_DIR)

    # 若 data-dir 本身直接包含 *_diff.csv，则直接作为 diff 目录
    if os.path.isdir(WORKSPACE) and any(f.endswith('_diff.csv') for f in os.listdir(WORKSPACE)):
        DATA_DIR = WORKSPACE
    elif os.path.isdir(os.path.join(WORKSPACE, 'diff')):
        DATA_DIR = os.path.join(WORKSPACE, 'diff')
    else:
        DATA_DIR = os.path.join(WORKSPACE, 'data/diff')

    if args.output_root:
        BASE_OUTPUT = args.output_root
    else:
        BASE_OUTPUT = os.path.join(WORKSPACE, 'output_scatter_slim')

    # 不再使用时间戳子目录，直接写入 output_root（反复执行直接覆盖）

    # frac 标签：frac=1 用目录名 "1"（depth 文件夹统一带深度名）
    print(f"Using DATA_DIR: {DATA_DIR}")
    print(f"Output directory: {BASE_OUTPUT}")
    print(f"Frac: {args.frac}")

    # ---- 确定输入文件列表 ----
    if args.csv:
        # --csv 模式：直接使用指定文件
        if not os.path.exists(args.csv):
            print(f"ERROR: CSV file not found: {args.csv}")
            return
        file_list = [(args.csv, os.path.splitext(os.path.basename(args.csv))[0])]
        print(f"Using CSV: {args.csv}")
    else:
        # 日期模式：按 --dates + --frac 拼接路径（frac=1 无后缀，兼容旧命名）
        print(f"Using DATA_DIR: {DATA_DIR}")
        frac_suffix = "" if args.frac == "1" else f"_{args.frac}"
        file_list = []
        for date in args.dates:
            if args.filtered:
                fpath = f'{DATA_DIR}/{date}_diff{frac_suffix}_filtered.csv'
            else:
                fpath = f'{DATA_DIR}/{date}_diff{frac_suffix}.csv'
            if not os.path.exists(fpath):
                print(f"  {date}: SKIP (file not found: {fpath})")
                continue
            file_list.append((fpath, date))

    if not file_list:
        print("No input files found. Exiting.")
        return

    configs = [
        {'use_log': args.linear, 'include_chrM': not args.exclude_chrm, 'filtered': args.filtered},
    ]

    for cfg in configs:
        for n in args.n_rows:
            print(f"\n=== N={n}, frac={args.frac}, chrM={'include' if cfg['include_chrM'] else 'exclude'}, "
                  f"filtered={cfg['filtered']}, scale=auto ===")

            for diff_file, label in file_list:
                df = pd.read_csv(diff_file)
                df_on = df[df['on_target'] == True].head(n)
                df_off = df[df['on_target'] == False].head(n)
                df_subset = pd.concat([df_on, df_off], ignore_index=True)
                print(f"  {label}: on={len(df_on)}, off={len(df_off)}")

                # 输出目录：{output_root}/{样品}/{frac}_N{n}/（depth + N 后缀）
                OUTPUT_DIR = f"{BASE_OUTPUT}/{label}/{args.frac}_N{n}"
                os.makedirs(OUTPUT_DIR, exist_ok=True)

                # ---- outlier 裁剪 ----
                if args.cap_auto:
                    # 仅用 ctrl 数据估算 cap（背景异常值集中在 ctrl），ctrl/sample 统一阈值
                    cap_val = auto_cap(df_subset['read_count_ctrl'].dropna().values)
                    n_before = len(df_subset)
                    df_subset = df_subset[df_subset['read_count_ctrl'] <= cap_val]
                    df_subset['read_count_sample'] = df_subset['read_count_sample'].clip(upper=cap_val)
                    n_filtered = n_before - len(df_subset)
                    n_on = (df_subset['on_target'] == True).sum()
                    n_off = (df_subset['on_target'] == False).sum()
                    print(f"  cap-auto: ctrl_cap={cap_val:.0f} → "
                          f"filtered {n_filtered} rows, kept on={n_on} off={n_off}")
                elif args.cap_pct is not None:
                    cap_val = np.percentile(df_subset['read_count_ctrl'].dropna(), args.cap_pct)
                    n_capped = (df_subset['read_count_ctrl'] > cap_val).sum()
                    df_subset = df_subset[df_subset['read_count_ctrl'] <= cap_val]
                    df_subset['read_count_sample'] = df_subset['read_count_sample'].clip(upper=cap_val)
                    n_on = (df_subset['on_target'] == True).sum()
                    n_off = (df_subset['on_target'] == False).sum()
                    print(f"  cap-pct={args.cap_pct}: ctrl capped at {cap_val:.0f} → "
                          f"filtered {n_capped} rows, kept on={n_on} off={n_off}")

                filter_str = ' (filtered)' if cfg['filtered'] else ''

                # y 轴范围：基于过滤 chrM 后的数据计算，与绘图一致
                if not cfg['include_chrM']:
                    df_ylim = df_subset[~df_subset['chromosome'].str.startswith('chrM')]
                else:
                    df_ylim = df_subset

                # ---- 统一 y 轴 scale（ctrl 与 sample 共用）----
                # 一律使用线性轴；基于 ctrl+sample 合并数据计算稳健 y 上限（默认 IQR 法），
                # 避免零星大点把 y 轴撑开。数据不删除，仅裁剪显示（set_ylim 天然行为）。
                use_log = False  # 一律线性轴
                ylim_shared, stats = compute_robust_ylim(df_ylim, cap_mode=args.y_cap)
                scale_str = 'linear'

                # ctrl 图
                title = f'{label} - read_count_ctrl (N={n}, {scale_str}){filter_str}'
                plot_scatter_slim(df_subset, 'read_count_ctrl',
                                  f'{OUTPUT_DIR}/read_count_ctrl.png',
                                  title, use_log=use_log,
                                  include_chrM=cfg['include_chrM'], ylim=ylim_shared)

                # sample 图（与 ctrl 共用相同 scale）
                title = f'{label} - read_count_sample (N={n}, {scale_str}){filter_str}'
                plot_scatter_slim(df_subset, 'read_count_sample',
                                  f'{OUTPUT_DIR}/read_count_sample.png',
                                  title, use_log=use_log,
                                  include_chrM=cfg['include_chrM'], ylim=ylim_shared)

                # ---- 写 max50.log（诊断 y 轴 scale 用，含完整统计量）----
                write_max50_log(OUTPUT_DIR, label, args.frac, n,
                                stats, ylim_shared, scale_str)

    print("\n✅ Done!")


if __name__ == '__main__':
    main()
