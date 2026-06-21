# utils/data_quality.py
# 数据质量分析绘图模块
# 修改：
#   1. 强制统一全局字体为 Arial，字号一致（通过 sns.set_theme + rcParams）
#   2. 样本相关性热图已移除颜色条标签和标题
#   3. 所有图表字体大小完全一致（标题12、轴标签10、刻度9、图例9）
#   4. 详细调试信息

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from matplotlib_venn import venn2, venn3
from matplotlib.patches import Ellipse, Patch
import re
import time
import sys

# ==============================
# 全局字体与样式设置（最高优先级）
# ==============================
print("[DEBUG] 开始配置全局字体为 Arial，字号统一...")
# 使用 seaborn 的主题，同时指定字体
sns.set_theme(style="white", font="Arial")
# 补充 matplotlib 的 rcParams，确保所有元素一致
mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans', 'Helvetica', 'sans-serif'],
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.titlesize': 12
})
print("[DEBUG] 全局字体设置完成：Arial, 标题12, 轴标签10, 刻度9, 图例9")

print(f"[LOAD] data_quality.py loaded from {__file__}", flush=True)

# ===================== 工具函数 =====================
def parse_sample_name(col_name):
    short = col_name.replace('LFQ intensity ', '')
    parts = short.split('-')
    try:
        if parts[0].upper() == 'WT':
            return ('WT', 0, int(parts[1]))
        else:
            return (parts[0], int(parts[1]), int(parts[2]))
    except:
        return ('zzz', 0, 0)

def sort_lfq_cols(lfq_cols):
    group_order = {'WT': 0, '100': 1, '200': 2}
    parsed = []
    for c in lfq_cols:
        g, t, r = parse_sample_name(c)
        parsed.append((group_order.get(g, 99), t, r, c))
    parsed.sort()
    return [x[-1] for x in parsed]

def extract_pure_group_prefix(sample_name):
    if sample_name.startswith('LFQ intensity '):
        short = sample_name[len('LFQ intensity '):]
    else:
        short = sample_name
    base = re.sub(r'[_-]\d+$', '', short)
    if base == short:
        parts = re.split(r'[_-]', short)
        if len(parts) > 1 and parts[-1].isdigit():
            base = '_'.join(parts[:-1])
    return base

def get_group_color(prefix):
    if '100' in prefix: return '#FF6347'
    elif '200' in prefix: return '#4682B4'
    elif 'WT' in prefix.upper(): return '#32CD32'
    else: return '#D3D3D3'

def format_group_display_name(prefix):
    if prefix.upper() == 'WT': return 'WT'
    parts = prefix.split('-')
    if len(parts) == 2: return f"{parts[0]}mM {parts[1]}h"
    return prefix

# ===================== 1. 缺失值热图 =====================
def plot_missing_heatmap(expr_df, lfq_cols):
    print("[DEBUG] 1. 缺失值热图 开始")
    if expr_df is None or lfq_cols is None: return None
    sorted_lfq = sort_lfq_cols(lfq_cols)
    mat = expr_df[sorted_lfq].apply(pd.to_numeric, errors='coerce')
    miss = (mat.isna() | (mat == 0)).astype(int)
    fig, ax = plt.subplots(figsize=(12, 8))
    sns.heatmap(miss.T, cmap=['#3498db', 'white'], cbar=False, xticklabels=False, yticklabels=True, ax=ax)
    ax.set_title("Missing Value Heatmap (Blue=Detected, White=Missing)")
    ax.set_xlabel("Proteins")
    ax.set_ylabel("Samples")
    plt.tight_layout()
    print("[DEBUG] 1. 缺失值热图 完成")
    return fig

# ===================== 2. 有效值数量 =====================
def plot_valid_values_per_sample(expr_df, lfq_cols):
    print("[DEBUG] 2. 有效值数量 开始")
    sorted_lfq = sort_lfq_cols(lfq_cols)
    mat = expr_df[sorted_lfq].apply(pd.to_numeric, errors='coerce')
    valid = (mat.notna() & (mat > 0)).sum()
    colors = [get_group_color(extract_pure_group_prefix(c)) for c in sorted_lfq]
    fig, ax = plt.subplots(figsize=(14, 6))
    ax.bar(sorted_lfq, valid.values, color=colors)
    ax.set_xlabel('Sample')
    ax.set_ylabel('Number of Valid Values')
    ax.set_title('Valid Values per Sample')
    ax.tick_params(axis='x', rotation=45)
    import matplotlib.patches as mpatches
    ax.legend(handles=[
        mpatches.Patch(color='#FF6347', label='100 group'),
        mpatches.Patch(color='#4682B4', label='200 group'),
        mpatches.Patch(color='#32CD32', label='WT'),
    ])
    plt.tight_layout()
    print("[DEBUG] 2. 有效值数量 完成")
    return fig

# ===================== 3. 韦恩图 =====================
def plot_venn_by_group(expr_df, lfq_cols):
    print("[DEBUG] 3. 韦恩图 开始")
    groups = {}
    for c in lfq_cols:
        prefix = extract_pure_group_prefix(c)
        groups.setdefault(prefix, []).append(c)
    order = ['WT', '100-6', '100-12', '100-24', '200-6', '200-12', '200-24']
    sorted_group_names = [gn for gn in order if gn in groups]
    figs = []
    for gname in sorted_group_names:
        cols = groups[gname]
        n = len(cols)
        if n < 2 or n > 3:
            continue
        print(f"[DEBUG] 韦恩图 - 组 {gname}，样本数 {n}")
        sets = {c: set(expr_df.index[expr_df[c].notna() & (expr_df[c] > 0)]) for c in cols}
        total = len(set.union(*sets.values()))
        title = f"{format_group_display_name(gname)} ({total})"
        labels = [f"R{i+1}" for i in range(n)]
        fig, ax = plt.subplots(figsize=(6, 6))
        try:
            if n == 2:
                A, B = sets[cols[0]], sets[cols[1]]
                venn2(subsets=(len(A), len(B), len(A & B)), set_labels=labels, ax=ax)
            elif n == 3:
                A, B, C = sets[cols[0]], sets[cols[1]], sets[cols[2]]
                subsets = (len(A), len(B), len(C),
                           len(A & B), len(B & C), len(A & C),
                           len(A & B & C))
                venn3(subsets=subsets, set_labels=labels, ax=ax)
            ax.set_title(title, fontweight='bold')
        except Exception as e:
            print(f"[ERROR] 韦恩图失败: {e}")
            ax.text(0.5, 0.5, f"Error: {e}", ha='center')
        figs.append((gname, fig))
    print("[DEBUG] 3. 韦恩图 完成")
    return figs

# ===================== 4. 肽段序列 =====================
def get_peptide_sequences_table(expr_df, protein_ids=None):
    print("[DEBUG] 4. 肽段序列 开始")
    peptide_col = next((c for c in expr_df.columns if 'peptide sequences' in c.lower()), None)
    if peptide_col is None:
        print("[DEBUG] 未找到 Peptide sequences 列")
        return None, None
    seqs = expr_df[peptide_col].dropna()
    if protein_ids is not None:
        seqs = seqs.loc[seqs.index.intersection(protein_ids)]
    merged = seqs.reset_index()
    merged.columns = ['Protein', 'Peptide Sequences']
    expanded_rows = []
    for _, row in merged.iterrows():
        if isinstance(row['Peptide Sequences'], str):
            for pep in row['Peptide Sequences'].split(';'):
                expanded_rows.append({'Protein': row['Protein'], 'Peptide Sequence': pep.strip()})
    expanded = pd.DataFrame(expanded_rows) if expanded_rows else pd.DataFrame()
    print(f"[DEBUG] 肽段序列: 合并 {len(merged)} 行, 展开 {len(expanded)} 行")
    return merged, expanded

# ===================== 5. 肽段长度分布 =====================
def plot_peptide_length_histogram(peptide_sequences_merged):
    print("[DEBUG] 5. 肽段长度分布 开始")
    lengths = []
    if peptide_sequences_merged is None or peptide_sequences_merged.empty:
        st.warning("无肽段序列数据")
        return None
    for seqs in peptide_sequences_merged['Peptide Sequences'].dropna():
        if isinstance(seqs, str):
            for pep in seqs.split(';'):
                lengths.append(len(pep.strip()))
    if not lengths:
        st.warning("未找到有效肽段序列")
        return None
    fig, ax = plt.subplots()
    n, bins, patches = ax.hist(lengths, bins='fd', color='#3498db', edgecolor='white')
    ax.set_xlabel('Peptide Length')
    ax.set_ylabel('Frequency')
    ax.set_title('Distribution of peptide lengths for all detected proteins in the selected samples.')
    print(f"[DEBUG] 肽段长度分布: {len(lengths)} 条, bins='fd', 实际 bins 数={len(n)}")
    return fig

# ===================== 6. 样本相关性热图（无标题，无色条标签，统一字体） =====================
def plot_sample_correlation_heatmap(expr_df, lfq_cols):
    print("[DEBUG] 6. 样本相关性热图（pandas corr, RdYlBu, 无标题/无颜色条标签）开始")
    if expr_df is None or lfq_cols is None or len(lfq_cols) < 2:
        print("[DEBUG] 数据不足，退出")
        return None

    # 数据预处理...
    mat = expr_df[lfq_cols].apply(pd.to_numeric, errors='coerce')
    valid_counts = mat.notna().sum(axis=1)
    mat = mat[valid_counts > 1]
    print(f"[DEBUG] 过滤后蛋白数: {mat.shape[0]} (至少2个样本有值)")

    if mat.isna().any().any():
        fill_value = np.nanquantile(mat.values, 0.01)
        mat = mat.fillna(fill_value)
        print(f"[DEBUG] 1%分位数填充值: {fill_value:.4f}")
    else:
        print("[DEBUG] 无缺失值，跳过填充")

    log_mat = np.log2(mat.values + 1)
    with np.errstate(divide='ignore', invalid='ignore'):
        z_mat = (log_mat - log_mat.mean(axis=1, keepdims=True)) / log_mat.std(axis=1, keepdims=True)
    z_mat = np.nan_to_num(z_mat, nan=0.0, posinf=0.0, neginf=0.0)
    row_vars = np.var(z_mat, axis=1)
    if np.sum(row_vars > 0) < 10:
        print("[DEBUG] 变异蛋白不足")
        return None
    top_idx = np.argsort(row_vars)[-500:]
    z_mat = z_mat[top_idx, :]
    n_prots, n_samps = z_mat.shape
    print(f"[DEBUG] 选择了 {n_prots} 个高变异蛋白，样本数={n_samps}")

    short_names = [c.replace('LFQ intensity ', '') for c in lfq_cols]
    z_df = pd.DataFrame(z_mat.T, columns=[f'Protein_{i}' for i in range(z_mat.T.shape[1])])
    z_df.index = short_names

    start_time = time.time()
    corr_df = z_df.T.corr()
    elapsed = time.time() - start_time
    print(f"[DEBUG] Pearson 计算完成，耗时 {elapsed:.2f}s")
    corr_values = corr_df.values[np.triu_indices(n_samps, k=1)]
    if len(corr_values) > 0:
        print(f"[DEBUG] Pearson 统计: min={np.min(corr_values):.3f}, max={np.max(corr_values):.3f}, "
              f"median={np.median(corr_values):.3f}, mean={np.mean(corr_values):.3f}")

    # 分组颜色
    color_list = []
    for c in lfq_cols:
        short = c.replace('LFQ intensity ', '')
        if short.startswith('WT'):
            color_list.append('#32CD32')
        elif short.startswith('100'):
            color_list.append('#FF6347')
        elif short.startswith('200'):
            color_list.append('#4682B4')
        else:
            color_list.append('#D3D3D3')
    row_colors = pd.Series(color_list, index=short_names, name='Group')
    print(f"[DEBUG] 分组颜色数: {len(row_colors)}")

    # 绘制聚类热图（无标题，无色条标签）
    print("[DEBUG] 绘制聚类热图 (seaborn clustermap) ...")
    g = sns.clustermap(corr_df,
                       cmap='RdYlBu',
                       center=0,
                       vmin=-1, vmax=1,
                       row_cluster=True, col_cluster=True,
                       row_colors=row_colors,
                       col_colors=row_colors,
                       figsize=(14, 12),
                       dendrogram_ratio=0.08,
                       linewidths=0.5,
                       linecolor='gray',
                       annot=True,
                       fmt='.2f',
                       annot_kws={'size': 7},
                       cbar_kws={'label': '', 'shrink': 0.8},
                       xticklabels=1, yticklabels=1)
    # 移除标题（不设置任何标题）
    g.ax_heatmap.set_title("")
    print("[DEBUG] 已移除热图标题和颜色条标签")

    # 分组图例
    legend_elements = [Patch(facecolor='#32CD32', label='WT'),
                       Patch(facecolor='#FF6347', label='100'),
                       Patch(facecolor='#4682B4', label='200')]
    g.fig.legend(handles=legend_elements, loc='upper right',
                 bbox_to_anchor=(1.15, 0.9), title='Group', frameon=True)
    g.fig.subplots_adjust(right=0.85)
    print("[DEBUG] 6. 样本相关性热图 完成")
    return g

# ===================== 7. PCA（统一字体） =====================
def plot_pca_raw_by_group(expr_df, lfq_cols, sample_info):
    print("[DEBUG] 7. PCA 开始")
    data = expr_df[lfq_cols].apply(pd.to_numeric, errors='coerce')
    q_val = data.stack().quantile(0.01)
    print(f"[DEBUG] 1%分位数填充值: {q_val:.4f}")
    data_filled = data.fillna(q_val)

    log_data = np.log2(data_filled.values + 1)
    X = StandardScaler().fit_transform(log_data.T)

    pca = PCA(n_components=2)
    scores = pca.fit_transform(X)
    var = pca.explained_variance_ratio_ * 100

    if sample_info is not None and 'SubGroup' in sample_info.columns:
        info_short = [str(idx).replace('LFQ intensity ', '') for idx in sample_info.index]
        subgroups = []
        for c in lfq_cols:
            short = c.replace('LFQ intensity ', '')
            match = [info_short[i] for i, x in enumerate(info_short) if x == short]
            if match:
                idx = info_short.index(match[0])
                subgroups.append(sample_info.iloc[idx]['SubGroup'])
            else:
                subgroups.append(extract_pure_group_prefix(c))
    else:
        subgroups = [extract_pure_group_prefix(c) for c in lfq_cols]

    plot_df = pd.DataFrame({
        'PC1': scores[:, 0],
        'PC2': scores[:, 1],
        'Group': subgroups
    })

    fig, ax = plt.subplots()
    group_colors = plt.cm.Set2(np.linspace(0, 1, len(plot_df['Group'].unique())))
    for g, col in zip(plot_df['Group'].unique(), group_colors):
        subset = plot_df[plot_df['Group'] == g]
        ax.scatter(subset['PC1'], subset['PC2'], label=g, color=col, s=60)
        if len(subset) >= 3:
            cov = np.cov(subset['PC1'], subset['PC2'])
            mean_vals = subset[['PC1', 'PC2']].mean(axis=0).to_numpy()
            vals, vecs = np.linalg.eigh(cov)
            angle = np.degrees(np.arctan2(*vecs[:, 1][::-1]))
            w, h = 2 * np.sqrt(vals * 5.991)
            ell = Ellipse(xy=mean_vals, width=w, height=h, angle=angle,
                          edgecolor=col, facecolor='none', linestyle='--', linewidth=1.5)
            ax.add_patch(ell)
    ax.set_xlabel(f'PC1 ({var[0]:.1f}%)')
    ax.set_ylabel(f'PC2 ({var[1]:.1f}%)')
    ax.set_title('PCA (Raw Data, 1% quantile imputation)')
    ax.legend()
    print("[DEBUG] 7. PCA 完成")
    return fig
