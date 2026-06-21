# utils/visualization.py
# 火山图及组合图绘制模块
# 修改：
#   1. 修正列名：识别 '-log10P' 或 'log10P'
#   2. 组合图布局、样式与 R 版本完全一致
#   3. 详细调试信息

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import matplotlib.gridspec as gridspec
import streamlit as st
import warnings
import time

# ---------------------------- 颜色映射（与R一致） ----------------------------
default_colors = {
    'Up': '#FF0000',
    'Down': '#0000FF',
    'Increase': '#C00000',
    'Decrease': '#0945A5',
    'NS': '#7f7e83'
}

# ---------------------------- 火山图（交互式，Plotly） ----------------------------
def plot_volcano(de_results, fc_up=1.2, fc_down=0.84, p_cut=0.05,
                 title='Volcano Plot', width=800, height=600):
    """使用 Plotly 绘制交互式火山图（用于页面展示）"""
    if de_results is None or de_results.empty:
        st.warning("No DE results to plot.")
        return go.Figure()

    fig = px.scatter(
        de_results, x='log2FC', y='-log10P', color='regulation',
        color_discrete_map=default_colors,
        hover_data=['Master protein IDs', 'FC', 'Pvalue', 'regulation_note'],
        labels={'log2FC': 'log2(Fold Change)', '-log10P': '-log10(P-value)'},
        title=title,
        category_orders={'regulation': ['Up', 'Down', 'Increase', 'Decrease', 'NS']}
    )
    fig.add_vline(x=np.log2(fc_up), line_dash="dash", line_color="gray")
    fig.add_vline(x=np.log2(fc_down), line_dash="dash", line_color="gray")
    fig.add_hline(y=-np.log10(p_cut), line_dash="dash", line_color="gray")
    fig.update_layout(width=width, height=height)
    return fig

# ---------------------------- 单图下载用（matplotlib，带标注） ----------------------------
def plot_volcano_single_annotated(de_results, fc_up=1.2, fc_down=0.84, p_cut=0.05,
                                  cols=None, point_size=4, title=''):
    """返回 matplotlib figure，带上/下调数目注释，无刻度"""
    if de_results is None or de_results.empty:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, 'No data', ha='center', va='center')
        return fig

    if cols is None:
        cols = default_colors

    fig, ax = plt.subplots()
    # 统一列名：可能为 -log10P 或 log10P
    log10_col = '-log10P' if '-log10P' in de_results.columns else 'log10P'
    for reg, color in cols.items():
        subset = de_results[de_results['regulation'] == reg]
        if not subset.empty:
            ax.scatter(subset['log2FC'], subset[log10_col], c=color, label=reg, alpha=0.6, s=point_size)

    ax.axvline(np.log2(fc_up), color='gray', linestyle='dashed')
    ax.axvline(np.log2(fc_down), color='gray', linestyle='dashed')
    ax.axhline(-np.log10(p_cut), color='gray', linestyle='dashed')

    counts = de_results['regulation'].value_counts()
    up = counts.get('Up', 0)
    down = counts.get('Down', 0)
    inc = counts.get('Increase', 0)
    dec = counts.get('Decrease', 0)
    y_annot = de_results[log10_col].max() * 1.05 if not de_results[log10_col].isnull().all() else 5
    if up > 0 or inc > 0:
        ax.text(4.5, y_annot, f'{up} ({inc})', color=cols['Up'], fontweight='bold', ha='right')
    if down > 0 or dec > 0:
        ax.text(-4.5, y_annot, f'{down} ({dec})', color=cols['Down'], fontweight='bold', ha='left')

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_title(title)
    return fig

# ---------------------------- 核心子图（用于组合图） ----------------------------
def plot_volcano_core_combined(df, fc_up, fc_down, p_cut, cols, point_size=1.8):
    """绘制无边框、无刻度的火山图子图，返回 matplotlib Figure"""
    print(f"[DEBUG] plot_volcano_core_combined 开始，数据行数={len(df)}")
    if df is None or df.empty:
        print("[DEBUG] df 为空")
        return None

    # 自动识别 -log10P 列
    if '-log10P' in df.columns:
        log10_col = '-log10P'
    elif 'log10P' in df.columns:
        log10_col = 'log10P'
    else:
        print("[DEBUG] 缺少 -log10P 或 log10P 列")
        return None

    required = {'log2FC', log10_col, 'regulation'}
    if not required.issubset(df.columns):
        print(f"[DEBUG] 缺少必要列: {required - set(df.columns)}")
        return None

    fig, ax = plt.subplots()
    for reg, color in cols.items():
        subset = df[df['regulation'] == reg]
        if not subset.empty:
            ax.scatter(subset['log2FC'], subset[log10_col], c=color, s=point_size, alpha=0.6)

    # 阈值线
    ax.axvline(np.log2(fc_up), color='gray', linestyle='dashed', linewidth=0.4)
    ax.axvline(np.log2(fc_down), color='gray', linestyle='dashed', linewidth=0.4)
    ax.axhline(-np.log10(p_cut), color='gray', linestyle='dashed', linewidth=0.4)

    ax.set_xlim(-8.5, 8.5)
    ax.set_ylim(0, 8)

    # 计数标注
    up = (df['regulation'] == 'Up').sum()
    down = (df['regulation'] == 'Down').sum()
    inc = (df['regulation'] == 'Increase').sum()
    dec = (df['regulation'] == 'Decrease').sum()
    y_text = 7.8
    if up > 0 or inc > 0:
        ax.text(4.5, y_text, f'{up}', color=cols['Up'], fontsize=10, fontweight='bold', ha='right')
        ax.text(4.6, y_text, f'({inc})', color=cols['Increase'], fontsize=10, fontweight='bold', ha='left')
    if down > 0 or dec > 0:
        ax.text(-4.5, y_text, f'{down}', color=cols['Down'], fontsize=10, fontweight='bold', ha='left')
        ax.text(-4.6, y_text, f'({dec})', color=cols['Decrease'], fontsize=10, fontweight='bold', ha='right')

    # 移除所有可视元素
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('')
    ax.set_ylabel('')
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_facecolor('none')
    fig.patch.set_facecolor('none')
    ax.set_position([0, 0, 1, 1])
    print("[DEBUG] plot_volcano_core_combined 完成")
    return fig

# ---------------------------- 构建组合图（完全模仿R） ----------------------------
def get_optimal_layout(n_plots):
    """根据子图数量计算最佳行列布局（与R一致）"""
    if n_plots <= 0: return (1, 1)
    if n_plots == 1: return (1, 1)
    if n_plots == 2: return (2, 1)
    if n_plots == 3: return (3, 1)
    if n_plots == 4: return (2, 2)
    if n_plots == 5: return (3, 2)
    if n_plots == 6: return (3, 2)
    if n_plots <= 8: return (4, 2)
    if n_plots <= 9: return (3, 3)
    if n_plots <= 12: return (4, 3)
    ncol = int(np.ceil(np.sqrt(n_plots)))
    nrow = int(np.ceil(n_plots / ncol))
    return (ncol, nrow)

def build_combined_plot(results_list, comp_names, main_title, fc_up=1.2, fc_down=0.84, p_cut=0.05,
                        sub_titles=None, point_size=1.8, cols=None):
    """
    生成与 R 版本 export_server.R 中 build_combined_plot 完全一致的组合火山图。
    results_list: DE结果 DataFrames 列表
    comp_names:   比较名称列表
    sub_titles:   自定义子标题列表（可选）
    """
    print(f"[DEBUG] build_combined_plot 开始，比较数量={len(results_list)}")
    if not results_list:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, 'No comparisons available', ha='center', va='center')
        return fig

    if cols is None:
        cols = default_colors

    n = len(results_list)
    ncol, nrow = get_optimal_layout(n)
    print(f"[DEBUG] 布局: {ncol} 列 x {nrow} 行")

    # 尺寸与 R 一致
    base_w = 3.5
    base_h = 3.5
    fig_width = ncol * base_w + 1.5
    fig_height = nrow * base_h + 2.0

    fig = plt.figure(figsize=(fig_width, fig_height), facecolor='white')
    gs = gridspec.GridSpec(nrow, ncol, figure=fig,
                           left=0.12, right=0.88, top=0.9, bottom=0.12,
                           wspace=0.4, hspace=0.4)

    # 确保所有输入数据都有正确的列名（统一处理）
    for i, res in enumerate(results_list):
        if '-log10P' in res.columns:
            res = res.rename(columns={'-log10P': 'log10P'})
            results_list[i] = res

    for i, (res, name) in enumerate(zip(results_list, comp_names)):
        print(f"[DEBUG] 绘制子图 {i+1}/{n}: {name}")
        ax = fig.add_subplot(gs[i // ncol, i % ncol])
        # 直接在当前ax上绘制核心内容，避免创建新figure
        if res is not None and not res.empty and 'log2FC' in res.columns and 'log10P' in res.columns and 'regulation' in res.columns:
            for reg, color in cols.items():
                subset = res[res['regulation'] == reg]
                if not subset.empty:
                    ax.scatter(subset['log2FC'], subset['log10P'], c=color, s=point_size, alpha=0.6)
            ax.axvline(np.log2(fc_up), color='gray', linestyle='dashed', linewidth=0.4)
            ax.axvline(np.log2(fc_down), color='gray', linestyle='dashed', linewidth=0.4)
            ax.axhline(-np.log10(p_cut), color='gray', linestyle='dashed', linewidth=0.4)
            ax.set_xlim(-8.5, 8.5)
            ax.set_ylim(0, 8)
            up = (res['regulation'] == 'Up').sum()
            down = (res['regulation'] == 'Down').sum()
            inc = (res['regulation'] == 'Increase').sum()
            dec = (res['regulation'] == 'Decrease').sum()
            y_text = 7.8
            if up > 0 or inc > 0:
                ax.text(4.5, y_text, f'{up}', color=cols['Up'], fontsize=10, fontweight='bold', ha='right')
                ax.text(4.6, y_text, f'({inc})', color=cols['Increase'], fontsize=10, fontweight='bold', ha='left')
            if down > 0 or dec > 0:
                ax.text(-4.5, y_text, f'{down}', color=cols['Down'], fontsize=10, fontweight='bold', ha='left')
                ax.text(-4.6, y_text, f'({dec})', color=cols['Decrease'], fontsize=10, fontweight='bold', ha='right')
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
            ax.set_facecolor('none')
        else:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center')

        # 子标题
        if sub_titles and i < len(sub_titles):
            stitle = sub_titles[i] if sub_titles[i] else name
        else:
            stitle = name
        ax.set_title(stitle, fontsize=10, pad=3)

    # 隐藏多余子图
    for j in range(n, nrow * ncol):
        ax = fig.add_subplot(gs[j // ncol, j % ncol])
        ax.axis('off')

    # 公共轴标签（位置与R一致）
    fig.text(0.5, 0.02, 'log2(Fold Change)', ha='center', fontsize=12)
    fig.text(0.02, 0.5, '-log10(P-value)', va='center', rotation='vertical', fontsize=12)

    # 主标题
    fig.suptitle(main_title, fontsize=14, fontweight='bold', y=0.97)

    print("[DEBUG] build_combined_plot 完成")
    return fig

# ---------------------------- 热图（保持不变） ----------------------------
def plot_heatmap(data_matrix, row_labels=None, col_labels=None, title='Heatmap'):
    print("[DEBUG] plot_heatmap 被调用")
    if data_matrix is None or data_matrix.empty:
        st.warning("No data for heatmap.")
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, 'No data', ha='center', va='center')
        return fig

    numeric_cols = data_matrix.select_dtypes(include=[np.number]).columns
    if len(numeric_cols) == 0:
        st.warning("热图数据中没有数值列。")
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, 'No numeric columns', ha='center', va='center')
        return fig

    mat = data_matrix[numeric_cols].apply(pd.to_numeric, errors='coerce').values
    if np.all(np.isnan(mat)):
        st.warning("处理后无可用的热图数据。")
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, 'No valid data', ha='center', va='center')
        return fig

    mat = mat[~np.all(np.isnan(mat), axis=1), :]
    mat = mat[:, ~np.all(np.isnan(mat), axis=0)]
    if mat.shape[0] < 2 or mat.shape[1] < 2:
        st.warning("数据量不足以绘制热图。")
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, 'Insufficient data', ha='center', va='center')
        return fig

    mat = np.nan_to_num(mat, nan=0.0, posinf=0.0, neginf=0.0)
    row_vars = np.var(mat, axis=1)
    mat = mat[row_vars > 0, :]
    if mat.shape[0] < 2:
        st.warning("所有行的方差为 0，无法聚类。")
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, 'No variable rows', ha='center', va='center')
        return fig

    with np.errstate(divide='ignore', invalid='ignore'):
        z_scored = (mat - mat.mean(axis=1, keepdims=True)) / mat.std(axis=1, keepdims=True)
    z_scored = np.nan_to_num(z_scored, nan=0.0, posinf=0.0, neginf=0.0)

    z_df = pd.DataFrame(z_scored, index=data_matrix.index[:len(mat)], columns=numeric_cols)
    g = sns.clustermap(
        z_df, cmap='RdBu_r', standard_scale=None,
        row_cluster=True, col_cluster=True,
        xticklabels=col_labels if col_labels else True,
        yticklabels=row_labels if row_labels else True,
        figsize=(10, 8)
    )
    g.ax_heatmap.set_title(title)
    print("[DEBUG] heatmap 完成")
    return g.fig

# ---------------------------- PCA ----------------------------
def plot_pca(data, lfq_cols, groups=None, batch=None):
    print("[DEBUG] plot_pca 被调用")
    if data is None or not lfq_cols:
        st.warning("No data for PCA.")
        return go.Figure()

    mat = data[lfq_cols].apply(pd.to_numeric, errors='coerce').fillna(0).values
    if mat.shape[0] < 2 or mat.shape[1] < 2:
        st.warning("Not enough data for PCA.")
        return go.Figure()

    X = StandardScaler().fit_transform(mat.T)
    pca = PCA(n_components=2)
    scores = pca.fit_transform(X)
    var = pca.explained_variance_ratio_ * 100

    plot_df = pd.DataFrame({
        'PC1': scores[:, 0],
        'PC2': scores[:, 1],
        'Sample': lfq_cols
    })
    if groups is not None:
        plot_df['Group'] = groups
    if batch is not None:
        plot_df['Batch'] = batch

    color_col = 'Group' if groups is not None else None
    fig = px.scatter(
        plot_df, x='PC1', y='PC2', color=color_col,
        hover_data=['Sample'],
        labels={'PC1': f'PC1 ({var[0]:.1f}%)', 'PC2': f'PC2 ({var[1]:.1f}%)'},
        title='PCA Plot'
    )
    print("[DEBUG] plot_pca 完成")
    return fig

# ---------------------------- Venn ----------------------------
def plot_venn(sets_dict, title='Venn Diagram'):
    print("[DEBUG] plot_venn 被调用")
    if len(sets_dict) < 2 or len(sets_dict) > 3:
        st.warning("Venn 图仅支持 2-3 个集合")
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, 'Need 2-3 sets', ha='center', va='center')
        return fig

    from matplotlib_venn import venn2, venn3
    fig, ax = plt.subplots()
    names = list(sets_dict.keys())
    sets = [set(v) for v in sets_dict.values()]

    if len(sets) == 2:
        venn2(subsets=(len(sets[0]), len(sets[1]), len(sets[0] & sets[1])), set_labels=names, ax=ax)
    elif len(sets) == 3:
        venn3(subsets=[len(s) for s in sets], set_labels=names, ax=ax)
    ax.set_title(title)
    print("[DEBUG] plot_venn 完成")
    return fig

# ---------------------------- UpSet ----------------------------
def plot_upset(sets_dict):
    print("[DEBUG] plot_upset 被调用")
    if not sets_dict:
        st.warning("No sets for UpSet plot.")
        return None

    from upsetplot import from_contents, UpSet
    contents = {}
    for name, ids in sets_dict.items():
        for i in ids:
            if i not in contents:
                contents[i] = []
            contents[i].append(name)

    series = pd.Series(contents)
    upset = UpSet(from_contents(series), subset_size='count', sort_by='cardinality')
    print("[DEBUG] plot_upset 完成")
    return upset.plot()
