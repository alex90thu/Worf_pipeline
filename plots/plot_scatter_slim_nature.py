#!/usr/bin/env python3
"""
plot_by_n_s 的 Nature 排版草稿版（v1，不覆盖原图）。

排版规格（用户指定）：
- 单图尺寸 120 mm × 40 mm
- 全部字体 Arial regular 6 pt（标题/轴标签/刻度/图例，无粗体）
- 图内不写 title（Nature 的标题放 figure caption），靠文件名区分 ctrl/sample
- 保留 on(蓝)/off(红) 分色与断轴高亮逻辑（数据不丢），断口低调处理
- 输出 PDF(矢量) + PNG(600 dpi)

用法：
    python plots/plot_scatter_slim_nature.py \
        --exp 20260508WORFT5-SNAP --frac 1 --n 10000
    # 默认即上述 example，输出到 output/5.Plots/nature_draft/
"""
import argparse
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from plot_by_n_s import chrom_center_to_x, compute_robust_ylim, CHIP_COL_BY_N

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DIFF_DIR = os.path.join(ROOT, 'output', '4.Diff')
OUT_BASE = os.path.join(ROOT, 'output', '5.Plots', 'nature_draft')

# ---- Nature 排版规格 ----
WIDTH_MM, HEIGHT_MM = 120.0, 40.0
MM2IN = 1 / 25.4
DPI = 600

mpl.rcParams.update({
    'font.family': 'Arial',
    'font.size': 6,
    'font.weight': 'normal',
    'axes.linewidth': 0.4,
    'axes.labelsize': 6,
    'axes.titlesize': 6,
    'axes.labelweight': 'normal',
    'xtick.major.width': 0.4,
    'ytick.major.width': 0.4,
    'xtick.major.size': 2.0,
    'ytick.major.size': 2.0,
    'xtick.labelsize': 6,
    'ytick.labelsize': 6,
    'legend.fontsize': 6,
    'legend.frameon': False,
    'lines.linewidth': 0.4,
    # Illustrator 友好：TrueType 内嵌（Type 42）保持文字为连续可编辑文本，
    # 而非默认 Type 3 字形（会把单词拆成单个字母）
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'svg.fonttype': 'none',
})

# slim scatter 经典配色：on 蓝 / off 红（交替深浅，便于连续染色体区分）
ON_COLORS  = ['#2874A6', '#5DADE2']
OFF_COLORS = ['#B03A2E', '#EC7063']


def plot_nature(df, y_col, out_stem, ylim, abs_max, alt_x=True,
                width_mm=120.0, height_mm=40.0, rasterize=False):
    """一张 Nature 版散点图（含可选低调断轴）。

    alt_x: 横轴交替标注（只标奇数 1,3,5,…21 与 X/Y/M，偶数留空），减少拥挤。
    width_mm/height_mm: 画布物理尺寸（默认 120×40mm）。
    rasterize: 数据点栅格化（嵌为位图），文本/轴/刻度保持矢量可编辑。
              期刊投稿常用：散点数量大，栅格化可大幅减小 PDF 体积。
    """
    df, chrom_offsets, chrom_lengths = chrom_center_to_x(df.copy(), include_chrM=True)

    df_on = df[df['on_target'] == True]
    df_off = df[df['on_target'] == False]

    y_on = df_on[y_col] + 1
    y_off = df_off[y_col] + 1

    abs_max_p1 = abs_max + 1
    has_upper_outliers = ((y_on > ylim[1]).any() or (y_off > ylim[1]).any())
    needs_break = bool(has_upper_outliers or abs_max_p1 > ylim[1])

    fig = plt.figure(figsize=(width_mm * MM2IN, height_mm * MM2IN))
    if needs_break:
        gs = fig.add_gridspec(2, 1, height_ratios=[1, 6], hspace=0.14,
                              left=0.14, right=0.985, top=0.94, bottom=0.30)
        ax_top = fig.add_subplot(gs[0])
        ax_bot = fig.add_subplot(gs[1])
    else:
        ax_bot = fig.add_subplot(1, 1, 1)
        fig.subplots_adjust(left=0.14, right=0.985, top=0.94, bottom=0.30)
        ax_top = None

    # ---- 逐染色体交替配色散点 ----
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
        elif alt_x:
            # 交替标注：只标奇数，偶数留空
            tick_labels.append(str(int(cnum)) if int(cnum) % 2 == 1 else '')
        else:
            tick_labels.append(str(int(cnum)))

    for i, cnum in enumerate(chrom_nums):
        on_color = ON_COLORS[i % 2]
        off_color = OFF_COLORS[i % 2]
        mask_on = (df_on['chrom_num'] == cnum)
        mask_off = (df_off['chrom_num'] == cnum)

        if mask_on.any():
            coll = ax_bot.scatter(df_on.loc[mask_on, 'absolute_bp'], y_on[mask_on],
                                  c=on_color, alpha=0.35, s=2, zorder=3, edgecolors='none')
            coll.set_rasterized(rasterize)
        if mask_off.any():
            coll = ax_bot.scatter(df_off.loc[mask_off, 'absolute_bp'], y_off[mask_off],
                                  c=off_color, alpha=0.35, s=2, zorder=2, edgecolors='none')
            coll.set_rasterized(rasterize)

        if needs_break:
            outlier_on = mask_on & (y_on > ylim[1])
            outlier_off = mask_off & (y_off > ylim[1])
            if outlier_on.any():
                coll = ax_top.scatter(df_on.loc[outlier_on, 'absolute_bp'], y_on[outlier_on],
                                      c=on_color, alpha=0.9, s=12, marker='^', zorder=4,
                                      edgecolors='none')
                coll.set_rasterized(rasterize)
            if outlier_off.any():
                coll = ax_top.scatter(df_off.loc[outlier_off, 'absolute_bp'], y_off[outlier_off],
                                      c=off_color, alpha=0.9, s=12, marker='^', zorder=4,
                                      edgecolors='none')
                coll.set_rasterized(rasterize)

    xlim_max = chrom_offsets.iloc[-1] + chrom_lengths.iloc[-1]
    ax_bot.set_xlim(0, xlim_max)
    ax_bot.set_xticks(tick_positions)
    ax_bot.set_xticklabels(tick_labels)
    ax_bot.set_ylim(ylim[0], ylim[1])

    # ---- 轴与断口 ----
    ax_bot.set_xlabel('Chromosome')
    ax_bot.set_ylabel('Read counts')
    ax_bot.tick_params(axis='x', which='major', length=1.5)
    ax_bot.tick_params(axis='y', which='major', length=1.5)

    if needs_break:
        ax_top.set_xlim(0, xlim_max)
        ax_top.set_ylim(ylim[1], abs_max_p1 + (abs_max_p1 - ylim[1]) * 0.45)
        ax_top.tick_params(axis='y', which='major', length=1.5)
        ax_top.tick_params(axis='x', bottom=False, labelbottom=False)

        # 去掉上下面板多余边框（只留左 + 各自外侧）
        ax_top.spines['bottom'].set_visible(False)
        ax_top.spines['top'].set_visible(False)
        ax_top.spines['right'].set_visible(False)
        ax_bot.spines['top'].set_visible(False)
        ax_bot.spines['right'].set_visible(False)

        # 低调断口斜线
        d = 0.006
        kw = dict(transform=ax_top.transAxes, color='#333333', clip_on=False, lw=0.4)
        ax_top.plot((-d, +d), (-d * 6, +d * 6), **kw)
        kw.update(transform=ax_bot.transAxes)
        ax_bot.plot((-d, +d), (1 - d, 1 + d), **kw)
    else:
        ax_bot.spines['top'].set_visible(False)
        ax_bot.spines['right'].set_visible(False)

    # ---- 图例（<2 序列必须有 legend）----
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=ON_COLORS[0],
               markersize=4, label=f'WORF-SEQ (n={len(df_on)})'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=OFF_COLORS[0],
               markersize=4, label=f'Control (n={len(df_off)})'),
    ]
    ax_bot.legend(handles=legend_elements, loc='upper right', fontsize=6,
                  frameon=False, handlelength=1.2, handletextpad=0.4,
                  borderaxespad=0.2, labelspacing=0.2)

    # ---- 保存：PDF 矢量 + PNG 600 dpi ----
    # 不用 bbox_inches='tight'，保证画布物理尺寸精确 = 120mm × 40mm
    # PDF 用 transparent=True：去掉白色背景色块（方便 Illustrator 叠底）
    fig.savefig(out_stem + '.pdf', dpi=DPI, transparent=True)
    fig.savefig(out_stem + '.png', dpi=DPI, transparent=False)
    print(f'  saved: {out_stem}.pdf / .png')
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description='Nature-style slim scatter (draft v1)')
    ap.add_argument('--exp', default='20260508WORFT5-SNAP')
    ap.add_argument('--frac', default='1', help='depth fraction; 1 = full depth')
    ap.add_argument('--n', type=int, default=10000)
    ap.add_argument('--out-root', default=None, help='override output root')
    ap.add_argument('--alt-x-labels', action=argparse.BooleanOptionalAction, default=True,
                    help='横轴交替标注（默认：只标奇数 1,3,5,… 与 X/Y/M，偶数留空）；'
                         '加 --no-alt-x-labels 全标')
    ap.add_argument('--width-mm', type=float, default=120.0, help='画布宽 mm（默认 120）')
    ap.add_argument('--height-mm', type=float, default=40.0, help='画布高 mm（默认 40）')
    ap.add_argument('--rasterize', action=argparse.BooleanOptionalAction, default=False,
                    help='数据点栅格化（默认关闭=全矢量；开启后散点嵌为位图、文本保持矢量可编辑）')
    args = ap.parse_args()

    frac_suffix = '' if str(args.frac) == '1' else f'_{args.frac}'
    diff_path = os.path.join(DIFF_DIR, args.exp, f'{args.exp}_diff{frac_suffix}.csv')
    print(f'load: {diff_path}')
    df = pd.read_csv(diff_path)

    # 芯片过滤：小芯片只画本芯片内位点（与 pipeline plot 口径一致）
    chip_col = CHIP_COL_BY_N.get(args.n)
    if chip_col and chip_col in df.columns:
        n_before = len(df)
        df = df[df[chip_col]].copy()
        print(f'  芯片过滤({chip_col}): {n_before} → {len(df)} 行')
    df_on = df[df['on_target'] == True].head(args.n)
    df_off = df[df['on_target'] == False].head(args.n)
    df_subset = pd.concat([df_on, df_off], ignore_index=True)
    print(f'  on={len(df_on)}, off={len(df_off)}')

    ylim, stats = compute_robust_ylim(df_subset, cap_mode='topk', topk=5)
    print(f'  ylim={ylim[0]:.1f}~{ylim[1]:.1f}, abs_max={stats["abs_max"]}, '
          f'break_axis={bool(stats["break_axis_enabled"])}')

    # 输出目录：120×40 默认 → nature_draft/；其他尺寸 → nature_draft_{宽}mm/；
    # 栅格化版追加 _raster 后缀，互不覆盖。
    if args.out_root:
        out_base = args.out_root
    else:
        base_name = ('nature_draft' if (args.width_mm == 120.0 and args.height_mm == 40.0)
                     else f'nature_draft_{int(args.width_mm)}mm')
        if args.rasterize:
            base_name += '_raster'
        out_base = os.path.join(ROOT, 'output', '5.Plots', base_name)
    out_dir = os.path.join(out_base, args.exp, f'{args.frac}_N{args.n}')
    os.makedirs(out_dir, exist_ok=True)

    for y_col, tag in (('read_count_ctrl', 'read_count_ctrl'),
                       ('read_count_sample', 'read_count_sample')):
        plot_nature(df_subset, y_col, os.path.join(out_dir, tag), ylim, stats['abs_max'],
                    alt_x=args.alt_x_labels,
                    width_mm=args.width_mm, height_mm=args.height_mm,
                    rasterize=args.rasterize)

    print(f'\n✅ done → {out_dir}')


if __name__ == '__main__':
    main()
