#!/usr/bin/env python3
"""
plot_by_n 的窄版变体（s = slim），适合 A4 纸排版：
- 图宽 18 cm（~7.09 英寸）
- 引入不对称断轴 (1:5) 结合特殊高亮标注 (Ceiling Annotation) 解决 Outlier 压缩问题
- 顶刊标准无背景交替配色：深蓝/浅蓝 (on-target)、深红/浅红 (off-target)，纯白底
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


def compute_robust_ylim(df, cap_mode="topk", topk=5):
    """基于 ctrl+sample 合并的 read_count 计算稳健 y 上限，并返回统计量。

    cap_mode:
        - "topk"（默认）: y_max = 合并数据的第 K 大值，保证顶部微缩轴内
          （read_count > y_max 的点）数量 <= K（期望 ≤5，至多 10）。这是
          为断轴图设计的：顶部只高亮少数几个真正的 outlier。
        - "iqr": y_max = Q3 + 1.5*IQR（兜底 >= P99）
        - "p99": y_max = P99
        - "max": y_max = 绝对最大值（backward compat）

    topk: 顶部允许的最大点数（默认 5）。
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

    if cap_mode == "topk":
        # 第 K 大值：保证 read_count > y_max_raw 的点数 <= K
        sorted_desc = combined.sort_values(ascending=False).values
        if n <= topk:
            y_max_raw = max(abs_max, 1.0)
        else:
            y_max_raw = int(sorted_desc[topk - 1])
        # 兜底：不能为 0 / 负
        y_max_raw = max(y_max_raw, 1.0)
        n_above = int((combined > y_max_raw).sum())
    elif cap_mode == "iqr":
        y_max_raw = q3 + 1.5 * iqr
        y_max_raw = max(y_max_raw, p99)
        y_max_raw = max(y_max_raw, 1.0)
        n_above = int((combined > y_max_raw).sum())
    elif cap_mode == "p99":
        y_max_raw = max(p99, 1.0)
        n_above = int((combined > y_max_raw).sum())
    else:  # "max"
        y_max_raw = max(abs_max, 1.0)
        n_above = 0

    y_max = y_max_raw + 1
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
        "n_above_topk": n_above,
        "cap_mode": cap_mode,
        "topk": topk,
        "top50_ctrl": df['read_count_ctrl'].dropna().astype(int).sort_values(ascending=False).head(50).tolist(),
        "top50_sample": df['read_count_sample'].dropna().astype(int).sort_values(ascending=False).head(50).tolist(),
    }
    return (y_min, y_max), stats


def write_max50_log(output_dir, label, frac, n, stats, ylim_shared, scale_str):
    """写入 max50.log，记录统计量与断轴位置判断信息"""
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
            f"[break-axis] cap_mode={stats.get('cap_mode')}, topk={stats.get('topk')}, "
            f"n_above (points above y_max) = {stats.get('n_above_topk')}",
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


def plot_scatter_slim(df, y_col, output_path, title, use_log=False, include_chrM=True, ylim=None, global_abs_max=None):
    """绘制窄版散点图，采用 1:5 不对称断轴高亮极限 Outlier"""
    df, chrom_offsets, chrom_lengths = chrom_center_to_x(df.copy(), include_chrM)

    df_on = df[df['on_target'] == True]
    df_off = df[df['on_target'] == False]

    y_on = df_on[y_col] + 1
    y_off = df_off[y_col] + 1
    
    # 获取全局绝对最大值
    abs_max_p1 = (global_abs_max + 1) if global_abs_max is not None else max(y_on.max() if not y_on.empty else 0, y_off.max() if not y_off.empty else 0)
    
    # 判定是否需要启动断轴：如果最大值突破了 IQR 上限的 1.5 倍
    needs_break = ylim is not None and abs_max_p1 > ylim[1] * 1.5

    fig_height = FIG_WIDTH_INCH * 0.42
    if needs_break:
        fig = plt.figure(figsize=(FIG_WIDTH_INCH, fig_height), dpi=300)
        # 核心逻辑：1:5 断轴，顶部只占不到 20% 空间，且由于刻度跨度大，视觉上起到了强压缩作用
        gs = fig.add_gridspec(2, 1, height_ratios=[1, 5], hspace=0.1)
        ax_top = fig.add_subplot(gs[0])
        ax_bot = fig.add_subplot(gs[1])
    else:
        fig, ax_bot = plt.subplots(figsize=(FIG_WIDTH_INCH, fig_height), dpi=300)
        ax_top = None

    # ---- 颜色定义 (顶刊标准) ----
    ON_COLORS  = ['#2874A6', '#5DADE2']   # 深蓝 / 浅蓝
    OFF_COLORS = ['#B03A2E', '#EC7063']   # 深红 / 浅红

    tick_labels = []
    tick_positions = []
    chrom_centers = chrom_offsets + (chrom_lengths / 2)
    chrom_nums = sorted([c for c in chrom_lengths.index if chrom_lengths[c] > 0])

    for cnum in chrom_nums:
        tick_positions.append(chrom_centers[cnum])
        if cnum == 23:
            tick_labels.append('X')
        elif cnum == 24:
            tick_labels.append('Y')
        elif cnum == 25:
            tick_labels.append('M')
        else:
            tick_labels.append(str(int(cnum)))

    # 动态阈值：只对真正脱颖而出的大点打字标号（如超出天花板的两倍以上），避免文字重叠
    text_threshold = max(ylim[1] * 2, abs_max_p1 * 0.5) if needs_break else float('inf')

    for i, cnum in enumerate(chrom_nums):
        on_color = ON_COLORS[i % 2]
        off_color = OFF_COLORS[i % 2]

        mask_on = (df_on['chrom_num'] == cnum)
        mask_off = (df_off['chrom_num'] == cnum)

        # 1. 正常绘制底部主群体
        if mask_on.any():
            ax_bot.scatter(df_on.loc[mask_on, 'absolute_bp'], y_on[mask_on],
                           c=on_color, alpha=0.35, s=2, zorder=3, edgecolors='none')
        if mask_off.any():
            ax_bot.scatter(df_off.loc[mask_off, 'absolute_bp'], y_off[mask_off],
                           c=off_color, alpha=0.35, s=2, zorder=2, edgecolors='none')

        # 2. 如果存在断轴，针对顶部区间绘制特殊形状与数值高亮
        if needs_break:
            outlier_on = mask_on & (y_on > ylim[1])
            outlier_off = mask_off & (y_off > ylim[1])

            if outlier_on.any():
                # 改用醒目的向上的三角形 `^`，加大字号并提高不透明度
                ax_top.scatter(df_on.loc[outlier_on, 'absolute_bp'], y_on[outlier_on],
                               c=on_color, alpha=0.9, s=25, marker='^', zorder=4, edgecolors='none')
                for bp, val in zip(df_on.loc[outlier_on, 'absolute_bp'], y_on[outlier_on]):
                    if val > text_threshold:
                        ax_top.annotate(f"{int(val-1)}", xy=(bp, val), xytext=(3, 3),
                                        textcoords='offset points', fontsize=6.0,
                                        color=on_color, ha='left', va='bottom', fontweight='bold')

            if outlier_off.any():
                ax_top.scatter(df_off.loc[outlier_off, 'absolute_bp'], y_off[outlier_off],
                               c=off_color, alpha=0.9, s=25, marker='^', zorder=4, edgecolors='none')
                for bp, val in zip(df_off.loc[outlier_off, 'absolute_bp'], y_off[outlier_off]):
                    if val > text_threshold:
                        ax_top.annotate(f"{int(val-1)}", xy=(bp, val), xytext=(3, 3),
                                        textcoords='offset points', fontsize=6.0,
                                        color=off_color, ha='left', va='bottom', fontweight='bold')

    # ---- 坐标轴限定与排版修饰 ----
    xlim_max = chrom_offsets.iloc[-1] + chrom_lengths.iloc[-1]
    ax_bot.set_xticks(tick_positions)
    ax_bot.set_xticklabels(tick_labels, fontsize=6)
    ax_bot.tick_params(axis='y', labelsize=6)
    ax_bot.set_xlim(0, xlim_max)

    if needs_break:
        ax_top.set_xlim(0, xlim_max)
        ax_top.tick_params(axis='y', labelsize=6)
        ax_bot.set_ylim(ylim[0], ylim[1])
        # 给顶部的文字预留一点空隙
        ax_top.set_ylim(ylim[1], abs_max_p1 + (abs_max_p1 - ylim[1]) * 0.45) 

        # 隐藏断口处的实线边框
        ax_top.spines['bottom'].set_visible(False)
        ax_top.spines['top'].set_visible(False)
        ax_top.spines['right'].set_visible(False)
        ax_top.tick_params(bottom=False, labelbottom=False)

        ax_bot.spines['top'].set_visible(False)
        ax_bot.spines['right'].set_visible(False)

        # 绘制断轴的对角斜线 `//`
        d = 0.012
        kwargs = dict(transform=ax_top.transAxes, color='#333333', clip_on=False, lw=0.8)
        ax_top.plot((-d, +d), (-d*5, +d*5), **kwargs) # 左上方断口
        kwargs.update(transform=ax_bot.transAxes)
        ax_bot.plot((-d, +d), (1-d, 1+d), **kwargs)   # 左下方断口

        ax_top.set_title(title, fontsize=10, pad=10)
        ax_bot.set_xlabel('Chromosomes', fontsize=8, fontweight='bold')
        ax_bot.set_ylabel('Read counts', fontsize=8, fontweight='bold')
        # 把 Y 轴标题文本向右移开，避免与坐标轴线重叠
        ax_bot.yaxis.set_label_coords(-0.05, 0.6)
    else:
        if ylim:
            ax_bot.set_ylim(ylim[0], ylim[1])
        ax_bot.spines['top'].set_visible(False)
        ax_bot.spines['right'].set_visible(False)
        ax_bot.set_title(title, fontsize=10, pad=10)
        ax_bot.set_xlabel('Chromosomes', fontsize=8, fontweight='bold')
        ax_bot.set_ylabel('Read counts', fontsize=8, fontweight='bold')

    # ---- 图例生成 ----
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=ON_COLORS[0], markersize=7, label=f'WORF-SEQ (n={len(df_on)})'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=OFF_COLORS[0], markersize=7, label=f'Control (n={len(df_off)})'),
    ]
    if needs_break:
        ax_top.legend(handles=legend_elements, loc='upper right', fontsize=6, frameon=False)
    else:
        ax_bot.legend(handles=legend_elements, loc='upper right', fontsize=6, frameon=False)
    
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
    parser.add_argument('--csv', type=str, default=None, help='Direct CSV file path')
    parser.add_argument('--dates', nargs='+', type=str, default=['20260404', '20260405', '20260407'])
    parser.add_argument('--frac', type=str, default='1')
    parser.add_argument('--y-cap', type=str, default='topk', choices=['topk', 'iqr', 'p99', 'max'],
                        help='Robust y-max strategy (default: topk = 第 K 大值, 顶部微缩轴 <=K 个点; '
                             'iqr = Q3+1.5*IQR; p99 = P99; max = absolute max)')
    parser.add_argument('--topk', type=int, default=5,
                        help='Break-axis top-region max points (default 5, used with --y-cap topk)')
    parser.add_argument('--output-root', type=str, default=None)
    args = parser.parse_args()

    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

    if args.data_dir:
        WORKSPACE = args.data_dir
    else:
        WORKSPACE = os.path.dirname(SCRIPT_DIR)

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

    print(f"Using DATA_DIR: {DATA_DIR}")
    print(f"Output directory: {BASE_OUTPUT}")

    if args.csv:
        if not os.path.exists(args.csv):
            print(f"ERROR: CSV file not found: {args.csv}")
            return
        file_list = [(args.csv, os.path.splitext(os.path.basename(args.csv))[0])]
    else:
        frac_suffix = "" if args.frac == "1" else f"_{args.frac}"
        file_list = []
        for date in args.dates:
            fpath = f'{DATA_DIR}/{date}_diff{frac_suffix}_filtered.csv' if args.filtered else f'{DATA_DIR}/{date}_diff{frac_suffix}.csv'
            if not os.path.exists(fpath):
                print(f"  {date}: SKIP (file not found: {fpath})")
                continue
            file_list.append((fpath, date))

    if not file_list:
        print("No input files found. Exiting.")
        return

    configs = [{'include_chrM': not args.exclude_chrm, 'filtered': args.filtered}]

    for cfg in configs:
        for n in args.n_rows:
            for diff_file, label in file_list:
                df = pd.read_csv(diff_file)
                df_on = df[df['on_target'] == True].head(n)
                df_off = df[df['on_target'] == False].head(n)
                df_subset = pd.concat([df_on, df_off], ignore_index=True)
                
                OUTPUT_DIR = f"{BASE_OUTPUT}/{label}/{args.frac}_N{n}"
                os.makedirs(OUTPUT_DIR, exist_ok=True)

                filter_str = ' (filtered)' if cfg['filtered'] else ''
                df_ylim = df_subset[~df_subset['chromosome'].str.startswith('chrM')] if not cfg['include_chrM'] else df_subset

                # 计算稳健 y 上限（默认 topk：顶部微缩轴 <=K 个点）
                ylim_shared, stats = compute_robust_ylim(df_ylim, cap_mode=args.y_cap, topk=args.topk)
                scale_str = 'linear'

                # 传入 global_abs_max 触发极值画图策略
                title_ctrl = 'Control'
                plot_scatter_slim(df_subset, 'read_count_ctrl', f'{OUTPUT_DIR}/read_count_ctrl.png',
                                  title_ctrl, include_chrM=cfg['include_chrM'], ylim=ylim_shared, global_abs_max=stats['abs_max'])

                title_sample = 'WORF-SEQ'
                plot_scatter_slim(df_subset, 'read_count_sample', f'{OUTPUT_DIR}/read_count_sample.png',
                                  title_sample, include_chrM=cfg['include_chrM'], ylim=ylim_shared, global_abs_max=stats['abs_max'])

                write_max50_log(OUTPUT_DIR, label, args.frac, n, stats, ylim_shared, scale_str)

    print("\n✅ Done!")

if __name__ == '__main__':
    main()