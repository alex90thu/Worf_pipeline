#!/usr/bin/env python3
"""
整合图：比较同一 depth（frac）下不同 N 值（Multiplexing Density / Panel Size）的 average read count。

- 统计口径沿用 plot_population：ON/OFF 分组，均值 + 标准误（SEM）
- 配色沿用 slim scatter：On-target 深海蓝 #2874A6，Off-target 深红 #B03A2E
- 上子图：Control（read_count_ctrl）
- 下子图：WORF-SEQ（read_count_sample）
- 横轴按 N 降序排列（N 大 → N 小），默认 N 集合 {100, 1,000, 10,000, 80,000}
- 同时输出 未过滤 / 过滤（OFF 组 ctrl>10 背景行）两个版本
- 输出：output/5.Plots/bar_contrast/{depth}_{not_filtered|filtered}_contrastbar.{png,pdf}
- WORF-SEQ 子图默认自动判轴：ON 均值跨 >2 个数量级时切 log（可用 --scale 覆盖）

用法：
    python plots/plot_n_escalation.py --frac 1
    python plots/plot_n_escalation.py --frac 0.5 \
        --exp-n 20260622chip-100:100 20260615sc-WORF:1000 20260508WORFT5-SNAP:10000 20260508WORFT5-SNAP:80000
"""
import argparse
import json
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# ==========================================
# 1. 全局样式设置 (NBT 极简美学)
# ==========================================
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['font.size'] = 9
mpl.rcParams['axes.linewidth'] = 0.8
mpl.rcParams['xtick.major.width'] = 0.8
mpl.rcParams['ytick.major.width'] = 0.8

# slim scatter 经典配色（on 蓝 / off 红）
COLOR_ON = '#2874A6'   # 深海蓝 (On-target)
COLOR_OFF = '#B03A2E'  # 深红 (Off-target)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DIFF_DIR = os.path.join(ROOT, 'output', '4.Diff')
OUT_DIR = os.path.join(ROOT, 'output', '5.Plots')

# 默认 (实验, N) 对（来自 config.n_rows；WORFT5-SNAP 同时保留 N=10,000 与 N=80,000）
DEFAULT_EXP_N = [
    ('20260622chip-100', 100),
    ('20260615sc-WORF', 1000),
    ('20260508WORFT5-SNAP', 10000),
    ('20260508WORFT5-SNAP', 80000),
]

# 绘图 N 值 → 需过滤的芯片标记列（diff 文件由 pipeline 附带 in_*_chip 列）。
# 80k 芯片是全量位点，不过滤；其余小芯片只统计本芯片内的位点。
CHIP_COL_BY_N = {100: 'in_100_chip', 1000: 'in_1k_chip', 10000: 'in_10k_chip'}


def frac_suffix(frac):
    """frac=1 无后缀（兼容旧命名），其他 frac 加 _<frac>"""
    return '' if str(frac) == '1' else f'_{frac}'


def load_group_stats(exp, n, frac, filter_suffix=''):
    """读取 {exp}_diff[{frac_suf}][{filter_suffix}].csv，on/off 各取前 N 行，返回 mean + SEM。

    filter_suffix: ''（未过滤）或 '_filtered'（过滤 OFF 组 ctrl>10 的背景行）。
    返回 {(col, 'on'/'off'): (mean, sem)}，col ∈ {read_count_ctrl, read_count_sample}
    """
    path = os.path.join(DIFF_DIR, exp, f'{exp}_diff{frac_suffix(frac)}{filter_suffix}.csv')
    if not os.path.exists(path):
        raise FileNotFoundError(f'missing diff file: {path}')
    df = pd.read_csv(path)
    # 芯片过滤：小芯片只统计本芯片内的位点（diff 已附带 in_*_chip 列）。
    # 过滤发生在取 top N 之前，避免混入非本芯片位点；80k 全量不过滤。
    chip_col = CHIP_COL_BY_N.get(n)
    if chip_col and chip_col in df.columns:
        n_before = len(df)
        df = df[df[chip_col]].copy()
        print(f'  芯片过滤({exp}, {chip_col}): {n_before} → {len(df)} 行')
    on = df[df['on_target'] == True].head(n)
    off = df[df['on_target'] == False].head(n)

    out = {}
    for col in ('read_count_ctrl', 'read_count_sample'):
        for key, sub in (('on', on), ('off', off)):
            v = sub[col].dropna()
            out[(col, key)] = (v.mean(), stats.sem(v) if len(v) else 0.0)
    return out


def plot_grouped_bars(ax, n_labels, on_m, on_sem, off_m, off_sem, title,
                      is_bottom=False, use_log=False):
    """极简分组柱状图核心函数（on 蓝 / off 红，mean ± SEM）"""
    x = np.arange(len(n_labels))
    width = 0.35
    err_kw = dict(lw=1, capsize=3, capthick=1, ecolor='#333333')

    if use_log:
        # log 轴：下误差棒封底（保证 value - low > 0），上误差棒取 SEM
        off_err = [np.minimum(off_sem, off_m * 0.9), off_sem]
        on_err = [np.minimum(on_sem, on_m * 0.9), on_sem]
        y_annot_off = off_m + off_sem
        y_annot_on = on_m + on_sem
    else:
        off_err, on_err = off_sem, on_sem
        y_annot_off, y_annot_on = off_m + off_sem, on_m + on_sem

    rects_off = ax.bar(x - width / 2, off_m, width, yerr=off_err,
                       color=COLOR_OFF, alpha=0.9, edgecolor='none',
                       error_kw=err_kw, label='Off-target')
    rects_on = ax.bar(x + width / 2, on_m, width, yerr=on_err,
                      color=COLOR_ON, alpha=0.9, edgecolor='none',
                      error_kw=err_kw, label='On-target')

    # 数值标签（极简风）
    for rect, val, err in zip(rects_off, off_m, off_sem):
        ax.annotate(f'{val:.1f}',
                    xy=(rect.get_x() + rect.get_width() / 2, val + err),
                    xytext=(0, 3), textcoords="offset points",
                    ha='center', va='bottom', fontsize=7.5, color='#333333')
    for rect, val, err in zip(rects_on, on_m, on_sem):
        ax.annotate(f'{val:.1f}',
                    xy=(rect.get_x() + rect.get_width() / 2, val + err),
                    xytext=(0, 3), textcoords="offset points",
                    ha='center', va='bottom', fontsize=7.5, color='#333333')

    # 轴线和网格修饰
    ax.set_ylabel('Read counts', fontweight='bold')
    ax.set_title(title, pad=12, fontweight='bold')
    ax.grid(axis='y', color='#E0E0E0', linestyle='--', linewidth=0.5, alpha=0.7)
    ax.set_axisbelow(True)          # 网格沉在柱子后方
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    if use_log:
        ax.set_yscale('log')

    # X 轴刻度处理
    if is_bottom:
        ax.set_xticks(x)
        ax.set_xticklabels(n_labels, fontweight='bold')
        ax.set_xlabel('Multiplexing Density (Panel Size)', fontweight='bold', labelpad=10)
    else:
        ax.tick_params(bottom=False)  # 隐藏顶部图 x 轴刻度小短线


def main():
    parser = argparse.ArgumentParser(description='N escalation integrated bar plot')
    parser.add_argument('--frac', type=str, default='1',
                        help='Saturation depth fraction (default: 1 = full depth)')
    parser.add_argument('--exp-n', nargs='+', type=str, default=None,
                        help='(Experiment:N) list, overrides defaults, e.g. 20260622chip-100:100 '
                             '20260615sc-WORF:1000 20260508WORFT5-SNAP:10000 20260508WORFT5-SNAP:80000')
    parser.add_argument('--scale', type=str, choices=['auto', 'linear', 'log'], default='auto',
                        help='WORF-SEQ subplot y-scale (default: auto → log if ON range >2 orders)')
    parser.add_argument('--output', type=str, default=None,
                        help='Output file stem base (default: output/5.Plots/bar_contrast/<frac>_<state>_contrastbar)')
    args = parser.parse_args()

    # ---- (实验, N) 映射（--exp-n 可整体覆盖默认）----
    pairs = list(DEFAULT_EXP_N)
    if args.exp_n:
        pairs = []
        for pair in args.exp_n:
            exp, n = pair.rsplit(':', 1)
            pairs.append((exp, int(n)))
    # 横轴降序：N 大 -> N 小
    pairs.sort(key=lambda kv: kv[1], reverse=True)
    n_labels = [f'N={n:,}' for _, n in pairs]
    print('N order (desc):', list(zip(n_labels, [e for e, _ in pairs])))

    # ---- 过滤状态：未过滤 / 过滤（OFF 组 ctrl>10 的背景行）----
    states = [('not_filtered', ''), ('filtered', '_filtered')]

    for state_tag, filter_suffix in states:
        print(f'\n[filter state] {state_tag}')

        # ---- 聚合统计 ----
        on_m, on_sem = {}, {}
        off_m, off_sem = {}, {}
        for col in ('read_count_ctrl', 'read_count_sample'):
            on_m[col], on_sem[col] = [], []
            off_m[col], off_sem[col] = [], []

        for exp, n in pairs:
            st = load_group_stats(exp, n, args.frac, filter_suffix)
            for col in ('read_count_ctrl', 'read_count_sample'):
                on_m[col].append(st[(col, 'on')][0])
                on_sem[col].append(st[(col, 'on')][1])
                off_m[col].append(st[(col, 'off')][0])
                off_sem[col].append(st[(col, 'off')][1])

        for col in ('read_count_ctrl', 'read_count_sample'):
            on_m[col] = np.array(on_m[col]);  on_sem[col] = np.array(on_sem[col])
            off_m[col] = np.array(off_m[col]); off_sem[col] = np.array(off_sem[col])

        # ---- WORF-SEQ 子图判轴 ----
        samp_vals = np.concatenate([on_m['read_count_sample'], off_m['read_count_sample']])
        if args.scale == 'log':
            use_log = True
        elif args.scale == 'linear':
            use_log = False
        else:  # auto
            pos = samp_vals[samp_vals > 0]
            use_log = bool(len(pos) >= 2 and np.log10(pos.max() / pos.min()) > 2)
        print(f'  WORF-SEQ subplot scale: {"log" if use_log else "linear"}')

        # ---- 绘图 ----
        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(6, 7), dpi=300, sharex=True)
        fig.subplots_adjust(hspace=0.25)

        # 上层：Control
        plot_grouped_bars(ax1, n_labels,
                          on_m['read_count_ctrl'], on_sem['read_count_ctrl'],
                          off_m['read_count_ctrl'], off_sem['read_count_ctrl'],
                          'Control', is_bottom=False)
        ax1.legend(frameon=False, loc='upper left', fontsize=8)

        # 下层：WORF-SEQ
        plot_grouped_bars(ax2, n_labels,
                          on_m['read_count_sample'], on_sem['read_count_sample'],
                          off_m['read_count_sample'], off_sem['read_count_sample'],
                          'WORF-SEQ', is_bottom=True, use_log=use_log)

        # 输出：output/5.Plots/bar_contrast/{depth}_{state}_contrastbar.{png,pdf}
        if args.output:
            stem = f'{args.output}_{state_tag}'
        else:
            stem = os.path.join(OUT_DIR, 'bar_contrast', f'{args.frac}_{state_tag}_contrastbar')
        os.makedirs(os.path.dirname(stem), exist_ok=True)
        plt.savefig(stem + '.pdf', transparent=True, bbox_inches='tight')
        plt.savefig(stem + '.png', dpi=300, bbox_inches='tight')
        print(f'  saved: {stem}.png / {stem}.pdf')
        plt.close()


if __name__ == '__main__':
    main()
