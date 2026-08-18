#!/usr/bin/env python3
"""
生成群体统计图：
- 柱状图：均值 + standard error
- Violin plot
- Box plot：max, 25%, median, 25%, min
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
from scipy import stats

# 绘图 N 值 → 需过滤的芯片标记列（diff 文件由 pipeline 附带 in_*_chip 列）。
# 80k 芯片是全量位点，不过滤；其余小芯片只画本芯片内的位点。
CHIP_COL_BY_N = {100: 'in_100_chip', 1000: 'in_1k_chip', 10000: 'in_10k_chip'}


def plot_bar_with_error(df, output_path, title):
    """柱状图：ON vs OFF 对比（ctrl 和 sample 分开展示），均值 + standard error
    
    修正：ON 和 OFF 构成对比，用颜色区分；ctrl 和 sample 分组
    """
    df_off = df[df['on_target'] == False]
    df_on = df[df['on_target'] == True]
    
    categories = ['OFF', 'ON']
    positions = [0, 1]  # OFF=0, ON=1
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=150)
    
    # ===== 左图：ctrl 的 ON vs OFF =====
    ctrl_means_off = df_off['read_count_ctrl'].mean()
    ctrl_means_on = df_on['read_count_ctrl'].mean()
    ctrl_sems_off = stats.sem(df_off['read_count_ctrl'].dropna())
    ctrl_sems_on = stats.sem(df_on['read_count_ctrl'].dropna())
    
    # 红色=OFF, 绿色=ON
    axes[0].bar(0, ctrl_means_off, yerr=ctrl_sems_off, capsize=5, color='red', alpha=0.7, edgecolor='black', label='OFF')
    axes[0].bar(1, ctrl_means_on, yerr=ctrl_sems_on, capsize=5, color='green', alpha=0.7, edgecolor='black', label='ON')
    
    # 添加数值标签
    max_ctrl_err = max(ctrl_sems_off, ctrl_sems_on) if ctrl_sems_off and ctrl_sems_on else 1
    axes[0].text(0, ctrl_means_off + max_ctrl_err * 1.5, f'{ctrl_means_off:.1f}', ha='center', va='bottom', fontsize=9)
    axes[0].text(1, ctrl_means_on + max_ctrl_err * 1.5, f'{ctrl_means_on:.1f}', ha='center', va='bottom', fontsize=9)
    
    axes[0].set_ylabel('Read Count', fontsize=12, fontweight='bold')
    axes[0].set_title('read_count_ctrl', fontsize=12, pad=10)
    axes[0].set_xticks(positions)
    axes[0].set_xticklabels(categories)
    axes[0].legend()
    axes[0].grid(axis='y', alpha=0.3)
    
    # ===== 右图：sample 的 ON vs OFF =====
    sample_means_off = df_off['read_count_sample'].mean()
    sample_means_on = df_on['read_count_sample'].mean()
    sample_sems_off = stats.sem(df_off['read_count_sample'].dropna())
    sample_sems_on = stats.sem(df_on['read_count_sample'].dropna())
    
    axes[1].bar(0, sample_means_off, yerr=sample_sems_off, capsize=5, color='red', alpha=0.7, edgecolor='black', label='OFF')
    axes[1].bar(1, sample_means_on, yerr=sample_sems_on, capsize=5, color='green', alpha=0.7, edgecolor='black', label='ON')
    
    max_sample_err = max(sample_sems_off, sample_sems_on) if sample_sems_off and sample_sems_on else 1
    axes[1].text(0, sample_means_off + max_sample_err * 1.5, f'{sample_means_off:.1f}', ha='center', va='bottom', fontsize=9)
    axes[1].text(1, sample_means_on + max_sample_err * 1.5, f'{sample_means_on:.1f}', ha='center', va='bottom', fontsize=9)
    
    axes[1].set_ylabel('Read Count', fontsize=12, fontweight='bold')
    axes[1].set_title('read_count_sample', fontsize=12, pad=10)
    axes[1].set_xticks(positions)
    axes[1].set_xticklabels(categories)
    axes[1].legend()
    axes[1].grid(axis='y', alpha=0.3)
    
    fig.suptitle(title, fontsize=14)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"  ✓ {title}")


def plot_violin(df, output_path, title):
    """Violin plot: ON vs OFF 对比（ctrl 和 sample 分开展示）
    
    修正：ON 和 OFF 构成对比，用颜色区分；ctrl 和 sample 分组
    """
    df_off = df[df['on_target'] == False]
    df_on = df[df['on_target'] == True]
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=150)
    
    # ===== 左图：ctrl 的 ON vs OFF =====
    data_ctrl_off = df_off['read_count_ctrl'].dropna()
    data_ctrl_on = df_on['read_count_ctrl'].dropna()
    parts1 = axes[0].violinplot([data_ctrl_off, data_ctrl_on], positions=[1, 2], showmeans=True, showmedians=True)
    # 红色=OFF, 绿色=ON
    parts1['bodies'][0].set_facecolor('red')
    parts1['bodies'][0].set_alpha(0.7)
    parts1['bodies'][1].set_facecolor('green')
    parts1['bodies'][1].set_alpha(0.7)
    axes[0].set_xticks([1, 2])
    axes[0].set_xticklabels(['OFF', 'ON'])
    axes[0].set_ylabel('Read Count', fontsize=12, fontweight='bold')
    axes[0].set_title('read_count_ctrl', fontsize=12)
    axes[0].grid(axis='y', alpha=0.3)
    
    # ===== 右图：sample 的 ON vs OFF =====
    data_sample_off = df_off['read_count_sample'].dropna()
    data_sample_on = df_on['read_count_sample'].dropna()
    parts2 = axes[1].violinplot([data_sample_off, data_sample_on], positions=[1, 2], showmeans=True, showmedians=True)
    # 红色=OFF, 绿色=ON
    parts2['bodies'][0].set_facecolor('red')
    parts2['bodies'][0].set_alpha(0.7)
    parts2['bodies'][1].set_facecolor('green')
    parts2['bodies'][1].set_alpha(0.7)
    axes[1].set_xticks([1, 2])
    axes[1].set_xticklabels(['OFF', 'ON'])
    axes[1].set_ylabel('Read Count', fontsize=12, fontweight='bold')
    axes[1].set_title('read_count_sample', fontsize=12)
    axes[1].grid(axis='y', alpha=0.3)
    
    # 添加图例
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='red', alpha=0.7, label='OFF'),
                       Patch(facecolor='green', alpha=0.7, label='ON')]
    fig.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(0.99, 0.99))
    
    fig.suptitle(title, fontsize=14)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"  ✓ {title}")


def plot_box(df, output_path, title):
    """Box plot: ON vs OFF 对比（ctrl 和 sample 分开展示）
    
    修正：ON 和 OFF 构成对比，用颜色区分；ctrl 和 sample 分组
    """
    df_off = df[df['on_target'] == False]
    df_on = df[df['on_target'] == True]
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), dpi=150)
    
    # ===== 左图：ctrl 的 ON vs OFF =====
    data_ctrl_off = df_off['read_count_ctrl'].dropna()
    data_ctrl_on = df_on['read_count_ctrl'].dropna()
    bp1 = axes[0].boxplot([data_ctrl_off, data_ctrl_on], positions=[1, 2], patch_artist=True,
                          showfliers=True, showmeans=True,
                          meanprops=dict(marker='D', markerfacecolor='yellow', markeredgecolor='black', markersize=8))
    # 红色=OFF, 绿色=ON
    for patch, color in zip(bp1['boxes'], ['red', 'green']):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    axes[0].set_xticks([1, 2])
    axes[0].set_xticklabels(['OFF', 'ON'])
    axes[0].set_ylabel('Read Count', fontsize=12, fontweight='bold')
    axes[0].set_title('read_count_ctrl', fontsize=12)
    axes[0].grid(axis='y', alpha=0.3)
    
    # 添加统计信息
    stats1 = f"OFF: min={data_ctrl_off.min():.0f}, max={data_ctrl_off.max():.0f}\nON: min={data_ctrl_on.min():.0f}, max={data_ctrl_on.max():.0f}"
    axes[0].text(0.02, 0.98, stats1, transform=axes[0].transAxes, fontsize=9,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # ===== 右图：sample 的 ON vs OFF =====
    data_sample_off = df_off['read_count_sample'].dropna()
    data_sample_on = df_on['read_count_sample'].dropna()
    bp2 = axes[1].boxplot([data_sample_off, data_sample_on], positions=[1, 2], patch_artist=True,
                          showfliers=True, showmeans=True,
                          meanprops=dict(marker='D', markerfacecolor='yellow', markeredgecolor='black', markersize=8))
    # 红色=OFF, 绿色=ON
    for patch, color in zip(bp2['boxes'], ['red', 'green']):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    axes[1].set_xticks([1, 2])
    axes[1].set_xticklabels(['OFF', 'ON'])
    axes[1].set_ylabel('Read Count', fontsize=12, fontweight='bold')
    axes[1].set_title('read_count_sample', fontsize=12)
    axes[1].grid(axis='y', alpha=0.3)
    
    stats2 = f"OFF: min={data_sample_off.min():.0f}, max={data_sample_off.max():.0f}\nON: min={data_sample_on.min():.0f}, max={data_sample_on.max():.0f}"
    axes[1].text(0.02, 0.98, stats2, transform=axes[1].transAxes, fontsize=9,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # 添加图例
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='red', alpha=0.7, label='OFF'),
                       Patch(facecolor='green', alpha=0.7, label='ON')]
    fig.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(0.99, 0.99))
    
    fig.suptitle(title, fontsize=14)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"  ✓ {title}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data-dir', type=str, required=True, help='Data directory (e.g., /root/.openclaw/workspace/worfscore_w15_new)')
    parser.add_argument('--n-rows', type=int, default=80000, help='Number of rows to load per group')
    parser.add_argument('--output-root', type=str, default=None,
                        help='Override output root dir (default: {data_dir}/output_scatter)')
    parser.add_argument('--dates', nargs='+', type=str, default=None,
                        help='Experiment/date prefixes; default: auto-discover from diff dir')
    parser.add_argument('--frac', type=str, default='1',
                        help='Saturation depth fraction (default: 1 = full depth; '
                             'diff file suffix _<frac>, frac=1 uses no suffix)')
    args = parser.parse_args()

    # 若 data-dir 本身直接包含 *_diff.csv，则直接作为 diff 目录
    if any(f.endswith('_diff.csv') for f in os.listdir(args.data_dir)):
        DATA_DIR = args.data_dir
    elif os.path.isdir(os.path.join(args.data_dir, 'diff')):
        DATA_DIR = os.path.join(args.data_dir, 'diff')
    else:
        DATA_DIR = os.path.join(args.data_dir, 'data/diff')

    if args.output_root:
        BASE_OUTPUT = args.output_root
    else:
        BASE_OUTPUT = os.path.join(args.data_dir, 'output_scatter')

    # 不再使用时间戳子目录，直接写入 output_root（反复执行直接覆盖）

    # frac 标签：frac=1 无后缀（兼容旧命名），其他 frac 加 _<frac>
    frac_suffix = "" if args.frac == "1" else f"_{args.frac}"
    
    print(f"Using DATA_DIR: {DATA_DIR}")
    print(f"Output directory: {BASE_OUTPUT}")
    print(f"Frac: {args.frac} (diff suffix: '{frac_suffix}')")
    
    # 指定配置：原始 + 过滤 两种
    configs = [
        {'include_chrM': False, 'filtered': False},
        {'include_chrM': False, 'filtered': True},
    ]
    
    # 确定实验/日期列表：优先 --dates，否则从 diff 目录自动发现合并后的 *_diff{suffix}.csv
    if args.dates:
        dates = list(args.dates)
    else:
        # 文件名形如 {exp}_diff[_<frac>].csv / {exp}_diff[_<frac>]_filtered.csv
        # 提取 {exp}：匹配 '_diff' 前缀前的部分
        dates = sorted({
            f.split('_diff')[0]
            for f in os.listdir(DATA_DIR)
            if (f.endswith('_diff.csv') or f.endswith('_diff_filtered.csv'))
            and '_on_diff' not in f and '_off_diff' not in f
            and (frac_suffix == "" or frac_suffix in f)
        })
        print(f"Auto-discovered experiments: {dates}")

    n = args.n_rows
    for cfg in configs:
        # filtered 用文件名后缀区分（同一 depth 目录内）
        filter_tag = "_filtered" if cfg['filtered'] else ""

        for date in dates:
            # 输出目录：{output_root}/{样品}/{frac}_N{n}/
            POP_OUTPUT_DIR = f"{BASE_OUTPUT}/{date}/{args.frac}_N{n}"
            os.makedirs(POP_OUTPUT_DIR, exist_ok=True)
            
            # 选择原始文件或过滤文件
            if cfg['filtered']:
                diff_file = f'{DATA_DIR}/{date}_diff{frac_suffix}_filtered.csv'
            else:
                diff_file = f'{DATA_DIR}/{date}_diff{frac_suffix}.csv'
            if not os.path.exists(diff_file):
                print(f"  {date}: SKIP (file not found: {diff_file})")
                continue
            
            df = pd.read_csv(diff_file)

            # 芯片过滤：N 对应小芯片时，只保留本芯片内的位点（diff 已附带 in_*_chip 列）。
            # 过滤发生在取 top N 之前，避免混入非本芯片位点；无标记列（旧 diff）时跳过。
            chip_col = CHIP_COL_BY_N.get(n)
            if chip_col and chip_col in df.columns:
                n_before = len(df)
                df = df[df[chip_col]].copy()
                print(f"  芯片过滤({chip_col}): {n_before} → {len(df)} 行")

            # 过滤染色体
            chrom_order = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY', 'chrM']
            chrom_order_no_m = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
            chroms = chrom_order if cfg['include_chrM'] else chrom_order_no_m
            df = df[df['chromosome'].isin(chroms)]
            
            # 取前 N 行
            df_on = df[df['on_target'] == True].head(n)
            df_off = df[df['on_target'] == False].head(n)
            df_subset = pd.concat([df_on, df_off], ignore_index=True)
            
            print(f"  {date} @{args.frac}: on={len(df_on)}, off={len(df_off)}, filtered={cfg['filtered']}")
            
            # 生成三种图（每张图包含 ctrl 和 sample 对比）
            # 1. 柱状图
            plot_bar_with_error(df_subset,
                               f"{POP_OUTPUT_DIR}/bar{filter_tag}.png",
                               f"{date} - Bar Plot (N={n})")
            
            # 2. Violin plot
            plot_violin(df_subset,
                       f"{POP_OUTPUT_DIR}/violin{filter_tag}.png",
                       f"{date} - Violin Plot (N={n})")
            
            # 3. Box plot
            plot_box(df_subset,
                    f"{POP_OUTPUT_DIR}/box{filter_tag}.png",
                    f"{date} - Box Plot (N={n})")
    
    print("\n✅ Done!")


if __name__ == '__main__':
    main()