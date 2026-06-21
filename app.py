# app.py
# Streamlit 主程序 - 蛋白质组学分析平台
# 修改：
#   1. Data Upload 页面同时显示原始数据维度和清洗后维度
#   2. Preprocessing 页面默认基线样本为 WT-1，插补方法为 quantile 时可自定义分位数
#   3. 大量调试信息
#   4. STRING PPI Network 页面已替换为调试好的版本

import streamlit as st
import pandas as pd
import numpy as np
import warnings
import plotly.express as px
import tempfile
import os
import requests
import gzip
import shutil
import re
import json
import time

# 工具模块
from utils.data_io import load_expression_data, load_sample_info
from utils.preprocessing import apply_missing_filter, impute_missing, run_combat_py
from utils.normalization import total_intensity_normalize
from utils.de_analysis import run_de_analysis
from utils.visualization import (
    plot_volcano, plot_heatmap, plot_venn, plot_upset,
    plot_volcano_single_annotated, build_combined_plot
)
from utils.cd_search import batch_cd_search, parse_cd_tsv
from utils.export_utils import create_excel_report
from utils.data_quality import (
    plot_missing_heatmap, plot_valid_values_per_sample,
    plot_venn_by_group, get_peptide_sequences_table,
    plot_peptide_length_histogram, plot_sample_correlation_heatmap,
    plot_pca_raw_by_group
)
from utils.eggnog_annotation import (
    run_eggnog_annotation_local,
    parse_eggnog_manual_file,
    EMAPPER_PATH
)
from utils.esm_search import (
    load_esm_model, build_reference_library, get_esm_embedding,
    find_similar_proteins, get_sequence_for_id, is_pretrained_available,
    ESM_AVAILABLE, batch_search_top1, fetch_uniprot_annotations_batch,
    parse_fasta, compute_fasta_hash, is_uniprot_accession, extract_accession
)
from utils.enrichment import get_go_background, enrich_go
from utils.llm_summary import generate_function_summary

# STRING 模块
from utils.string_ppi import (
    extract_uniprot_from_eggnog,
    run_blast_mapping,
    call_string_api,
    build_ppi_network_html
)

warnings.filterwarnings("ignore", category=RuntimeWarning)

st.set_page_config(page_title="Proteomics Platform", layout="wide")

# ==================== 辅助函数 ====================

def parse_manual_uniprot_ids(text: str) -> list:
    """解析手动输入的 UniProt ID 列表"""
    ids = [
        x.strip()
        for x in text.replace(",", "\n").replace(";", "\n").splitlines()
        if x.strip()
    ]

    if any(x.startswith(">") for x in ids):
        st.error("这里需要输入 UniProt Accession，不是 FASTA 序列。FASTA 请粘贴到上方 FASTA 输入框。")
        st.stop()

    uniprot_pattern = re.compile(
        r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9](?:[A-Z][A-Z0-9]{2}[0-9])?)$"
    )

    invalid_ids = [x for x in ids if not uniprot_pattern.match(x)]

    if invalid_ids:
        st.error("以下内容不是有效的 UniProt Accession，请检查输入：")
        st.code("\n".join(invalid_ids[:20]))
        st.stop()

    return ids


# ==================== 初始化 session state ====================
if 'raw_expr_df' not in st.session_state:
    st.session_state.raw_expr_df = None       # 新增：原始未清洗数据
if 'expr_df' not in st.session_state:
    st.session_state.expr_df = None
if 'lfq_cols' not in st.session_state:
    st.session_state.lfq_cols = None
if 'sample_names' not in st.session_state:
    st.session_state.sample_names = None
if 'sample_info' not in st.session_state:
    st.session_state.sample_info = None
if 'processed' not in st.session_state:
    st.session_state.processed = None
if 'norm_data' not in st.session_state:
    st.session_state.norm_data = None
if 'groups' not in st.session_state:
    st.session_state.groups = {}
if 'comparisons' not in st.session_state:
    st.session_state.comparisons = []
if 'de_results' not in st.session_state:
    st.session_state.de_results = {}
if 'page' not in st.session_state:
    st.session_state.page = "Home"
if 'cleaning_stats' not in st.session_state:
    st.session_state.cleaning_stats = None
if 'cd_result' not in st.session_state:
    st.session_state.cd_result = None
if 'eggnog_result' not in st.session_state:
    st.session_state.eggnog_result = None

# ==================== 侧边栏导航 ====================
st.sidebar.title("Navigation")
pages = [
    "Home", "Data Upload", "Data Quality", "Preprocessing",
    "Define Groups & Comparisons", "Differential Analysis",
    "Visualization", "CD-Search", "Export", "EggNOG Annotation",
    "ESM2 Similarity Search", "STRING PPI Network"
]
current_page = pages.index(st.session_state.page) if st.session_state.page in pages else 0
page = st.sidebar.radio("Go to", pages, index=current_page)
st.session_state.page = page

# ==================== 页面逻辑 ====================

# ---------- Home ----------
if page == "Home":
    st.title("Proteomics Differential Analysis Platform")
    st.write("从 MaxQuant 数据到差异分析、注释和报告。")
    if st.button("开始分析 → 上传数据"):
        st.session_state.page = "Data Upload"
        st.rerun()

# ---------- Data Upload ----------
elif page == "Data Upload":
    st.header("数据上传")
    expr_file = st.file_uploader("上传 proteinGroups.txt", type=['txt'])
    sample_file = st.file_uploader("上传样本信息 (CSV 或 Excel)", type=['csv', 'xlsx'])
    
    if expr_file is not None and sample_file is not None:
        df_clean, lfq_cols, sample_names, col_type, clean_stats, df_raw = load_expression_data(expr_file)
        si = load_sample_info(sample_file)
        if df_clean is not None and si is not None:
            st.session_state.raw_expr_df = df_raw          # 原始未清洗数据
            st.session_state.expr_df = df_clean
            st.session_state.lfq_cols = lfq_cols
            st.session_state.sample_names = sample_names
            st.session_state.sample_info = si
            st.session_state.cleaning_stats = clean_stats
            print(f"[DEBUG] 数据上传成功: 原始行数={df_raw.shape[0]}, 清洗后行数={df_clean.shape[0]}, LFQ列数={len(lfq_cols)}")
            st.success("数据加载成功！")
            # 显示两种维度
            st.write(f"**原始数据维度**（上传文件直接解析）: {df_raw.shape[0]} 行 × {df_raw.shape[1]} 列")
            st.write(f"**清洗后数据维度**（去除Reverse/污染物等）: {df_clean.shape[0]} 行 × {df_clean.shape[1]} 列")
            st.write("前5行（清洗后）:", df_clean.head())
            st.write("样本信息:", si.head())

# ---------- Data Quality ----------
elif page == "Data Quality":
    st.header("数据质量分析")
    if st.session_state.expr_df is None:
        st.warning("请先上传数据")
    else:
        expr_df = st.session_state.expr_df
        lfq_cols = st.session_state.lfq_cols
        sample_info = st.session_state.sample_info

        st.subheader("缺失值热图")
        print("[DEBUG] Data Quality: 绘制缺失值热图")
        fig_miss = plot_missing_heatmap(expr_df, lfq_cols)
        if fig_miss:
            st.pyplot(fig_miss)

        st.subheader("每个样本的有效值数量")
        print("[DEBUG] Data Quality: 绘制 Valid Values per Sample")
        fig_valid = plot_valid_values_per_sample(expr_df, lfq_cols)
        if fig_valid:
            st.pyplot(fig_valid)

        st.subheader("按处理组的韦恩图")
        print("[DEBUG] Data Quality: 绘制分组韦恩图")
        venn_figs = plot_venn_by_group(expr_df, lfq_cols)
        if venn_figs:
            cols = st.columns(3)
            for i, (gname, fig) in enumerate(venn_figs):
                with cols[i % 3]:
                    st.pyplot(fig)
        else:
            st.info("没有可绘制的韦恩图（每组需要2-3个样本）")

        st.subheader("肽段长度分布")
        merged, _ = get_peptide_sequences_table(expr_df)
        if merged is not None:
            print("[DEBUG] Data Quality: 绘制肽段长度分布")
            fig_len = plot_peptide_length_histogram(merged)
            if fig_len:
                st.pyplot(fig_len)
        else:
            st.info("无肽段序列数据")

        st.subheader("样本相关性热图")
        print("[DEBUG] Data Quality: 绘制样本相关性热图")
        fig_cor = plot_sample_correlation_heatmap(expr_df, lfq_cols)
        if fig_cor:
            st.pyplot(fig_cor)

        st.subheader("PCA 分析 (原始数据)")
        print("[DEBUG] Data Quality: 绘制 PCA")
        fig_pca = plot_pca_raw_by_group(expr_df, lfq_cols, sample_info)
        if fig_pca:
            st.pyplot(fig_pca)

# ---------- Preprocessing ----------
elif page == "Preprocessing":
    st.header("预处理")
    if st.session_state.expr_df is None:
        st.warning("请先上传数据")
    else:
        if st.session_state.cleaning_stats is not None:
            st.subheader("Data Cleaning Summary")
            cs = st.session_state.cleaning_stats
            st.write(f"原始蛋白数: {cs['original']}")
            st.write(f"移除 Reverse hits: {cs['reverse_removed']}")
            st.write(f"移除 Potential contaminants: {cs['contaminant_removed']}")
            st.write(f"移除 CON_ contaminants: {cs['con_removed']}")
            st.write(f"清洗后保留蛋白数: {cs['retained']}")
            st.write("---")

        max_missing = st.slider("最大缺失比例", 0.0, 1.0, 0.5)
        impute_method = st.selectbox("插补方法", ["none", "knn", "ppca", "quantile", "minvalue"], index=3)
        
        # ---------- 分位数自定义（当选择 quantile 时） ----------
        quantile_prob = 0.01
        if impute_method == "quantile":
            st.info("默认使用 1% 分位数（0.01）填充缺失值，可修改：")
            quantile_prob = st.number_input("分位数", min_value=0.001, max_value=0.5, value=0.01, step=0.01,
                                            help="用于替换缺失值的分位数，越小则填充值越低")
            print(f"[DEBUG] 用户设置 quantile 分位数: {quantile_prob}")

        # ---------- 基线样本默认选择 WT-1 ----------
        lfq_cols = st.session_state.lfq_cols
        # 默认索引：寻找 LFQ intensity WT-1，找不到则用第一个
        default_baseline = "LFQ intensity WT-1"
        if default_baseline in lfq_cols:
            default_idx = lfq_cols.index(default_baseline)
        else:
            default_idx = 0
            print(f"[DEBUG] 未找到 {default_baseline}，默认使用第一个列作为基线")
        baseline_sample = st.selectbox("归一化基线样本", lfq_cols, index=default_idx)
        print(f"[DEBUG] 选择基线样本: {baseline_sample}")

        if st.button("运行预处理"):
            before_filter = st.session_state.expr_df.shape[0]
            print(f"[DEBUG] 用户点击运行预处理: 缺失阈值={max_missing}, 插补={impute_method}, 基线={baseline_sample}")
            df = st.session_state.expr_df.copy()
            df = apply_missing_filter(df, st.session_state.lfq_cols, max_missing, st.session_state.sample_info)
            after_filter = df.shape[0]
            print(f"[DEBUG] 缺失过滤后保留蛋白数: {after_filter} (移除 {before_filter - after_filter})")
            # 将分位数参数传递给插补函数
            df = impute_missing(df, st.session_state.lfq_cols, impute_method, quantile_prob=quantile_prob, k=10)
            df = total_intensity_normalize(df, st.session_state.lfq_cols, baseline_sample, st.session_state.sample_names)
            st.session_state.processed = df
            final_na = df[st.session_state.lfq_cols].isna().sum().sum()
            print(f"[DEBUG] 预处理完成后 LFQ 列缺失值总数: {final_na}")
            st.success(f"预处理完成！缺失过滤前 {before_filter} 个蛋白，过滤后保留 {after_filter} 个。插补后缺失值数量: {final_na}")
            st.write("处理后数据 (仅显示前5行, 可横向滚动):")
            st.dataframe(df.head())

# ---------- Define Groups & Comparisons ----------
elif page == "Define Groups & Comparisons":
    st.header("定义分组与比较")
    
    if st.session_state.sample_info is not None and st.session_state.lfq_cols is not None:
        si = st.session_state.sample_info
        lfq_cols = st.session_state.lfq_cols
        
        short_from_cols = []
        for c in lfq_cols:
            if c.startswith('LFQ intensity '):
                short_from_cols.append(c[len('LFQ intensity '):])
            else:
                short_from_cols.append(c)
        short_std = [s.replace('.', '_').strip() for s in short_from_cols]
        col_std = [c.replace('.', '_').strip() for c in lfq_cols]
        info_index_std = [str(idx).replace('.', '_').strip() for idx in si.index]
        
        print(f"[DEBUG] 短名示例 (前3): {short_std[:3]}")
        
        match_map = {}
        for i, (full_col, short, full_std) in enumerate(zip(lfq_cols, short_std, col_std)):
            if short in info_index_std:
                match_map[full_col] = info_index_std[info_index_std.index(short)]
            elif full_std in info_index_std:
                match_map[full_col] = full_std
            else:
                for idx_name in info_index_std:
                    if short in idx_name or idx_name in short:
                        match_map[full_col] = idx_name
                        break
        
        st.info(f"匹配到 {len(match_map)} 个样本")
        
        groups = {}
        used_column = None
        if 'SubGroup' in si.columns:
            for full_col, matched_idx in match_map.items():
                subgroup_val = si.loc[si.index[info_index_std.index(matched_idx)], 'SubGroup']
                groups.setdefault(subgroup_val, []).append(full_col)
            used_column = 'SubGroup'
        elif 'Group' in si.columns:
            for full_col, matched_idx in match_map.items():
                group_val = si.loc[si.index[info_index_std.index(matched_idx)], 'Group']
                groups.setdefault(group_val, []).append(full_col)
            used_column = 'Group'
        else:
            prefixes = {}
            for s in short_std:
                parts = s.split('_')
                if len(parts) >= 2:
                    prefix = '_'.join(parts[:-1])
                else:
                    prefix = s
                prefixes.setdefault(prefix, []).append(s)
            for prefix, short_list in prefixes.items():
                full_list = [lfq_cols[short_std.index(s)] for s in short_list if s in short_std]
                if full_list:
                    groups[prefix] = full_list
            used_column = 'prefix'
        
        st.session_state.groups = groups
        
        if groups:
            st.write(f"### 分组结果 (基于 {used_column})")
            for gname, cols in groups.items():
                st.write(f"**{gname}** ({len(cols)} 样本): {', '.join(cols)}")
        else:
            st.warning("未能生成分组，请检查样本信息表。")
            st.write("### 备选方案：根据样本名前缀自动分组")
            prefixes = {}
            for s in short_std:
                parts = s.split('_')
                if len(parts) >= 2:
                    prefix = '_'.join(parts[:-1])
                else:
                    prefix = s
                prefixes.setdefault(prefix, []).append(s)
            st.json(prefixes)
            if st.button("使用前缀分组"):
                new_groups = {}
                for prefix, short_list in prefixes.items():
                    full_list = [lfq_cols[short_std.index(s)] for s in short_list if s in short_std]
                    if full_list:
                        new_groups[prefix] = full_list
                st.session_state.groups = new_groups
                st.rerun()
    
    if st.session_state.groups:
        group_list = list(st.session_state.groups.keys())
        st.write("---")
        st.subheader("手动添加比较")
        col1, col2 = st.columns(2)
        with col1:
            treat = st.selectbox("处理组", group_list, key="treat")
        with col2:
            ctrl = st.selectbox("对照组", group_list, key="ctrl")
        if st.button("添加比较"):
            if treat != ctrl:
                comp_name = f'{treat} vs {ctrl}'
                if comp_name not in [c['name'] for c in st.session_state.comparisons]:
                    st.session_state.comparisons.append({
                        'treat': treat,
                        'ctrl': ctrl,
                        'name': comp_name
                    })
                    st.success(f"比较 '{comp_name}' 已添加")
                else:
                    st.warning("该比较已存在")
            else:
                st.error("处理组和对照组不能相同")
        
        st.subheader("一键添加所有 vs WT 比较")
        wt_candidates = [g for g in group_list if 'wt' in g.lower() or g.lower() == 'control']
        if wt_candidates:
            wt_group = wt_candidates[0]
            other_groups = [g for g in group_list if g != wt_group]
            if other_groups:
                if st.button(f"添加所有其他组 vs {wt_group} 的比较"):
                    added = 0
                    for treat_grp in other_groups:
                        comp_name = f'{treat_grp} vs {wt_group}'
                        if comp_name not in [c['name'] for c in st.session_state.comparisons]:
                            st.session_state.comparisons.append({
                                'treat': treat_grp,
                                'ctrl': wt_group,
                                'name': comp_name
                            })
                            added += 1
                    if added > 0:
                        st.success(f"已添加 {added} 个比较")
                    else:
                        st.info("所有比较已存在")
            else:
                st.info("没有其他组可比较")
        else:
            st.warning("未找到 WT/Control 组，无法批量添加。请手动添加比较。")
    
    if st.session_state.comparisons:
        st.write("---")
        st.subheader("当前比较")
        for comp in st.session_state.comparisons:
            treat_samples = st.session_state.groups.get(comp['treat'], [])
            ctrl_samples = st.session_state.groups.get(comp['ctrl'], [])
            st.write(f"- **{comp['name']}**：处理组 **{comp['treat']}** ({len(treat_samples)} 样本) vs 对照组 **{comp['ctrl']}** ({len(ctrl_samples)} 样本)")
    else:
        st.info("尚未定义任何比较")

# ---------- Differential Analysis ----------
elif page == "Differential Analysis":
    st.header("差异分析")
    if not st.session_state.comparisons:
        st.warning("请先在 'Define Groups & Comparisons' 页面添加比较")
    elif st.session_state.processed is None:
        st.warning("请先运行预处理")
    else:
        method = st.selectbox("统计方法", ["t-test", "wilcoxon", "limma"])
        fc_up = st.number_input("FC up >", value=1.2)
        fc_down = st.number_input("FC down <", value=0.84)
        p_cut = st.number_input("P-value 阈值", value=0.05)
        
        st.write("---")
        st.subheader("有效重复数阈值")
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            min_treat_valid = st.number_input("Treatment group min valid replicates", min_value=1, value=2)
        with col2:
            min_ctrl_valid = st.number_input("Control group min valid replicates", min_value=1, value=2)
        with col3:
            min_rep_ttest = st.number_input("t-test 最少重复数", min_value=1, value=2)
        with col4:
            pass
        
        col5, col6 = st.columns(2)
        with col5:
            min_rep_inc = st.number_input("标记为 Increase 所需处理组最少重复数", min_value=1, value=2)
        with col6:
            min_rep_dec = st.number_input("标记为 Decrease 所需对照组最少重复数", min_value=1, value=2)
        
        st.write("---")
        st.subheader("蛋白质过滤")
        min_unique_pep = st.number_input("最小 Unique Peptides (≥ 此值)", min_value=1, value=2, step=1)
        use_pep_filter = st.checkbox("启用 Unique Peptides 过滤", value=True)
        
        if st.button("运行所有比较"):
            norm_df = st.session_state.processed.copy()
            if use_pep_filter and 'Unique peptides' in norm_df.columns:
                print(f"[DEBUG] 差异分析前 Unique peptide 过滤: 阈值 >= {min_unique_pep}")
                before = norm_df.shape[0]
                norm_df['Unique peptides'] = pd.to_numeric(norm_df['Unique peptides'], errors='coerce')
                norm_df = norm_df.dropna(subset=['Unique peptides'])
                norm_df = norm_df[norm_df['Unique peptides'] >= min_unique_pep]
                print(f"[DEBUG] 过滤后剩余蛋白质: {norm_df.shape[0]} (从 {before})")
            
            results = {}
            for comp in st.session_state.comparisons:
                treat_cols = st.session_state.groups[comp['treat']]
                ctrl_cols = st.session_state.groups[comp['ctrl']]
                print(f"[DEBUG] 差异分析: {comp['name']}, 方法={method}")
                res = run_de_analysis(norm_df, treat_cols, ctrl_cols,
                                      fc_up=fc_up, fc_down=fc_down, p_cut=p_cut,
                                      method=method,
                                      min_treat_valid=min_treat_valid,
                                      min_ctrl_valid=min_ctrl_valid,
                                      min_rep_ttest=min_rep_ttest,
                                      min_rep_inc=min_rep_inc,
                                      min_rep_dec=min_rep_dec)
                if res is not None:
                    results[comp['name']] = res
            st.session_state.de_results = results
            st.success("差异分析完成！")
            for name, res in results.items():
                total = len(res)
                up = (res['regulation']=='Up').sum()
                down = (res['regulation']=='Down').sum()
                inc = (res['regulation']=='Increase').sum()
                dec = (res['regulation']=='Decrease').sum()
                st.write(f"**{name}**：总蛋白 {total}，Up {up}，Down {down}，Increase {inc}，Decrease {dec}")

# ---------- Visualization ----------
elif page == "Visualization":
    st.header("可视化")
    if not st.session_state.de_results:
        st.warning("请先运行差异分析")
    else:
        comp_name = st.selectbox("选择比较", list(st.session_state.de_results.keys()))
        res = st.session_state.de_results[comp_name]
        fig_vol = plot_volcano(res)
        st.plotly_chart(fig_vol)
        
        if st.button("生成热图"):
            sig_prots = res[res['regulation'].isin(['Up','Down','Increase','Decrease'])]['Master protein IDs']
            if not sig_prots.empty:
                full_mat = st.session_state.processed.set_index('Master protein IDs')
                sel_mat = full_mat.loc[full_mat.index.isin(sig_prots)]
                sel_mat = sel_mat.select_dtypes(include=[np.number]).apply(pd.to_numeric, errors='coerce')
                if not sel_mat.empty:
                    heatmap_fig = plot_heatmap(sel_mat)
                    st.pyplot(heatmap_fig)
                else:
                    st.warning("未找到差异蛋白对应的数值表达数据")
            else:
                st.warning("没有差异蛋白可绘制热图")

# ---------- CD-Search (保存结果到 session_state) ----------
elif page == "CD-Search":
    st.header("NCBI CD-Search")
    st.subheader("方式一：自动搜索（粘贴 FASTA）")
    fasta_text = st.text_area("粘贴 FASTA 序列", height=200)
    if st.button("运行自动 CD-Search"):
        if fasta_text.strip():
            with st.spinner("正在搜索..."):
                print("[DEBUG] 开始自动 CD-Search")
                result = batch_cd_search(fasta_text)
                if result is not None:
                    st.session_state.cd_result = result
                    st.dataframe(result)
                else:
                    st.error("自动搜索失败")
        else:
            st.warning("请输入 FASTA 序列")

    st.subheader("方式二：手动上传 CD-Search 结果文件")
    uploaded_cd = st.file_uploader("上传 CD-Search 结果 (TSV/TXT)", type=['tsv', 'txt', 'csv'])
    if uploaded_cd is not None:
        try:
            print("[DEBUG] 开始解析手动上传的 CD-Search 文件")
            content = uploaded_cd.read().decode('utf-8')
            df = parse_cd_tsv(content)
            if df is not None and not df.empty:
                st.session_state.cd_result = df
                st.success(f"成功解析 {len(df)} 条 CD-Search 记录")
                st.dataframe(df, height=600)
            else:
                st.error("解析失败，请检查文件格式")
        except Exception as e:
            st.error(f"读取文件出错: {e}")

    if st.session_state.cd_result is not None:
        st.info(f"当前已加载 CD-Search 注释：{st.session_state.cd_result.shape[0]} 行")
    else:
        st.info("尚未加载任何 CD-Search 注释")

# ---------- Export ----------
elif page == "Export":
    st.header("导出")
    if st.session_state.processed is None or not st.session_state.de_results:
        st.warning("请先完成预处理和差异分析")
    else:
        st.subheader("单图下载")
        selected_comp = st.selectbox("选择比较", list(st.session_state.de_results.keys()))
        plot_format = st.selectbox("图片格式", ["svg", "tiff"])
        col1, col2 = st.columns(2)
        with col1:
            width = st.number_input("宽度 (inch)", min_value=5, max_value=30, value=10)
        with col2:
            height = st.number_input("高度 (inch)", min_value=5, max_value=30, value=8)
        custom_title = st.text_input("自定义标题 (留空使用默认)", value="")
        
        if st.button("下载单图"):
            print(f"[DEBUG] 下载单图: {selected_comp}, 格式={plot_format}, 尺寸={width}x{height}")
            res = st.session_state.de_results[selected_comp]
            title = custom_title if custom_title else selected_comp
            fig = plot_volcano_single_annotated(res, fc_up=1.2, fc_down=0.84, p_cut=0.05, title=title)
            if fig:
                with tempfile.NamedTemporaryFile(suffix=f".{plot_format}", delete=False) as tmp:
                    fig.savefig(tmp.name, dpi=300, bbox_inches='tight')
                    with open(tmp.name, 'rb') as f:
                        st.download_button(
                            label="点击下载单图",
                            data=f,
                            file_name=f"volcano_{selected_comp}.{plot_format}",
                            mime=f"image/{plot_format}"
                        )
        
        st.subheader("组合图下载")
        combined_title = st.text_input("组合图主标题", value="Combined Volcano Plots")
        sub_titles_input = st.text_area("子图标题 (每行一个，按比较顺序)", value="")
        if st.button("下载组合图"):
            print(f"[DEBUG] 下载组合图，标题={combined_title}")
            results_list = [st.session_state.de_results[comp['name']] for comp in st.session_state.comparisons]
            comp_names = [comp['name'] for comp in st.session_state.comparisons]
            if sub_titles_input.strip():
                sub_titles = [s.strip() for s in sub_titles_input.split('\n') if s.strip()]
                while len(sub_titles) < len(comp_names):
                    sub_titles.append('')
            else:
                sub_titles = comp_names
            fig = build_combined_plot(results_list, comp_names, main_title=combined_title,
                                      fc_up=1.2, fc_down=0.84, p_cut=0.05,
                                      sub_titles=sub_titles)
            if fig:
                with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
                    fig.savefig(tmp.name, dpi=300, bbox_inches='tight')
                    with open(tmp.name, 'rb') as f:
                        st.download_button(
                            label="点击下载组合图",
                            data=f,
                            file_name="combined_volcano.png",
                            mime="image/png"
                        )
        
        st.subheader("Excel 报告导出")
        if st.button("生成 Excel 报告"):
            print("[DEBUG] 开始生成 Excel 报告")
            wb = create_excel_report(
                raw_data=st.session_state.raw_expr_df,   # 使用真正的原始数据
                clean_data=st.session_state.expr_df,
                norm_data=st.session_state.processed,
                de_results=st.session_state.de_results,
                comparisons=st.session_state.comparisons,
                groups=st.session_state.groups,
                sample_info=st.session_state.sample_info
            )
            with tempfile.NamedTemporaryFile(delete=False, suffix='.xlsx') as tmp:
                wb.save(tmp.name)
                with open(tmp.name, 'rb') as f:
                    st.download_button('下载 Excel 报告', f, file_name='proteomics_report.xlsx')

# ---------- EggNOG Annotation (含数据库下载助手) ----------
elif page == "EggNOG Annotation":
    st.header("EggNOG Annotation")
    st.markdown("使用本地 eggnog-mapper 或手动上传结果文件进行功能注释。")
    
    # ========== 数据库下载工具 ==========
    st.subheader("📦 数据库下载")
    st.info("数据库只需下载一次。如果已下载，可跳过此步骤。")
    
    data_dir = os.path.join(os.path.dirname(EMAPPER_PATH), "data")
    db_ready = os.path.exists(os.path.join(data_dir, "eggnog.db")) and \
               os.path.exists(os.path.join(data_dir, "eggnog_proteins.dmnd"))
    
    if db_ready:
        st.success("✅ 数据库已就绪，可以直接使用本地注释。")
    else:
        st.warning("⚠️ 数据库尚未下载或未完整解压。请使用下方按钮下载。")
        
        col1, col2 = st.columns(2)
        with col1:
            if st.button("下载 eggnog.db (约 8 GB)", use_container_width=True):
                url = "https://downloads.eggnogdb.org/emapper/emapperdb-5.0.2/eggnog.db.gz"
                dest = os.path.join(data_dir, "eggnog.db.gz")
                os.makedirs(data_dir, exist_ok=True)
                
                with st.spinner("正在下载 eggnog.db，请耐心等待（约 8 GB）..."):
                    try:
                        with requests.get(url, stream=True, timeout=30) as r:
                            r.raise_for_status()
                            total_size = int(r.headers.get('content-length', 0))
                            progress_bar = st.progress(0)
                            downloaded = 0
                            with open(dest, 'wb') as f:
                                for chunk in r.iter_content(chunk_size=8192):
                                    if chunk:
                                        f.write(chunk)
                                        downloaded += len(chunk)
                                        if total_size > 0:
                                            progress_bar.progress(min(downloaded / total_size, 1.0))
                        st.success("下载完成！请点击旁边的解压按钮。")
                    except Exception as e:
                        st.error(f"下载失败: {e}")
        
        with col2:
            if st.button("解压 eggnog.db.gz", use_container_width=True):
                gz_path = os.path.join(data_dir, "eggnog.db.gz")
                db_path = os.path.join(data_dir, "eggnog.db")
                if not os.path.exists(gz_path):
                    st.error("找不到 eggnog.db.gz，请先下载。")
                else:
                    with st.spinner("正在解压..."):
                        try:
                            with gzip.open(gz_path, 'rb') as src, open(db_path, 'wb') as dst:
                                shutil.copyfileobj(src, dst)
                            os.remove(gz_path)
                            st.success("解压成功！")
                        except Exception as e:
                            st.error(f"解压失败: {e}")
        
        col3, col4 = st.columns(2)
        with col3:
            if st.button("下载 eggnog_proteins.dmnd (约 4 GB)", use_container_width=True):
                url = "https://downloads.eggnogdb.org/emapper/emapperdb-5.0.2/eggnog_proteins.dmnd.gz"
                dest = os.path.join(data_dir, "eggnog_proteins.dmnd.gz")
                os.makedirs(data_dir, exist_ok=True)
                
                with st.spinner("正在下载 eggnog_proteins.dmnd (约 4 GB)..."):
                    try:
                        with requests.get(url, stream=True, timeout=30) as r:
                            r.raise_for_status()
                            total_size = int(r.headers.get('content-length', 0))
                            progress_bar = st.progress(0)
                            downloaded = 0
                            with open(dest, 'wb') as f:
                                for chunk in r.iter_content(chunk_size=8192):
                                    if chunk:
                                        f.write(chunk)
                                        downloaded += len(chunk)
                                        if total_size > 0:
                                            progress_bar.progress(min(downloaded / total_size, 1.0))
                        st.success("下载完成！请点击旁边的解压按钮。")
                    except Exception as e:
                        st.error(f"下载失败: {e}")
        
        with col4:
            if st.button("解压 eggnog_proteins.dmnd.gz", use_container_width=True):
                gz_path = os.path.join(data_dir, "eggnog_proteins.dmnd.gz")
                dmnd_path = os.path.join(data_dir, "eggnog_proteins.dmnd")
                if not os.path.exists(gz_path):
                    st.error("找不到 eggnog_proteins.dmnd.gz，请先下载。")
                else:
                    with st.spinner("正在解压..."):
                        try:
                            with gzip.open(gz_path, 'rb') as src, open(dmnd_path, 'wb') as dst:
                                shutil.copyfileobj(src, dst)
                            os.remove(gz_path)
                            st.success("解压成功！")
                        except Exception as e:
                            st.error(f"解压失败: {e}")

    st.markdown("---")
    
    # ========== 本地运行 ==========
    tab_local, tab_upload = st.tabs(["本地运行 (推荐)", "手动上传"])
    
    with tab_local:
        st.subheader("粘贴 FASTA 序列，使用本地 emapper.py 进行注释")
        fasta_text = st.text_area("粘贴 FASTA 序列", height=200)
        if st.button("开始本地注释"):
            if fasta_text.strip():
                df = run_eggnog_annotation_local(fasta_text)
                if df is not None:
                    st.success("注释完成！")
                    st.session_state.eggnog_result = df
                    st.dataframe(df, height=600)
                    csv = df.to_csv(index=False).encode('utf-8')
                    st.download_button(
                        label="下载注释结果 (CSV)",
                        data=csv,
                        file_name="eggnog_annotations.csv",
                        mime="text/csv"
                    )
                else:
                    st.error("注释失败，请检查序列或重试。")
            else:
                st.warning("请输入 FASTA 序列")
        st.info("首次使用可能需要下载数据库（约20GB），请耐心等待。")
    
    # ========== 手动上传 ==========
    with tab_upload:
        st.subheader("手动上传 eggNOG 结果文件")
        st.markdown("""
        1. 访问 [eggNOG-mapper 网站](http://eggnog-mapper.embl.de/) 提交序列。
        2. 下载 **annotations** 文件（TSV 格式）。
        3. 在此上传。
        """)
        uploaded_eggnog = st.file_uploader("上传 eggNOG 注释结果 (TSV)", type=['tsv', 'txt', 'csv'])
        if uploaded_eggnog is not None:
            try:
                content = uploaded_eggnog.getvalue().decode('utf-8')
                df = parse_eggnog_manual_file(content)
                if df is not None and not df.empty:
                    st.session_state.eggnog_result = df
                    st.success(f"成功解析 {len(df)} 条注释记录")
                    st.dataframe(df, height=600)
                    csv = df.to_csv(index=False).encode('utf-8')
                    st.download_button(
                        label="下载注释结果 (CSV)",
                        data=csv,
                        file_name="eggnog_annotations.csv",
                        mime="text/csv"
                    )
                else:
                    st.error("解析失败，请检查文件格式")
            except Exception as e:
                st.error(f"读取文件出错: {e}")

    if st.session_state.eggnog_result is not None:
        st.info(f"当前已加载 eggNOG 注释：{st.session_state.eggnog_result.shape[0]} 行")
    else:
        st.info("尚未加载任何 eggNOG 注释")

# ---------- ESM2 Similarity Search ----------
elif page == "ESM2 Similarity Search":
    st.header("🔬 ESM2 蛋白质功能相似性搜索")
    st.caption("功能1: 单条蛋白相似检索 → 功能3: AI摘要生成 | 功能2: 批量注释与富集分析")
    
    if not ESM_AVAILABLE:
        st.error("**ESM2 功能不可用**，请安装 PyTorch 和 fair-esm：")
        st.code("pip install torch fair-esm", language="bash")
        st.info("安装完成后请重启 Streamlit 应用。")
        st.stop()

    model_state_key = "esm_model"
    if model_state_key not in st.session_state:
        with st.spinner("正在加载 ESM2 模型（首次加载可能需要下载，约 140MB）..."):
            model, alphabet, batch_converter = load_esm_model()
            if model is not None:
                st.session_state.esm_model = model
                st.session_state.esm_alphabet = alphabet
                st.session_state.esm_batch_converter = batch_converter
                st.success("模型加载成功！")
            else:
                st.error("模型加载失败，请检查控制台输出。")
                st.stop()

    model = st.session_state.esm_model
    batch_converter = st.session_state.esm_batch_converter

    st.subheader("📚 参考库管理")
    pretrained_available = is_pretrained_available()
    use_pretrained = st.checkbox("使用预构建 Swiss-Prot 参考库", value=pretrained_available, disabled=not pretrained_available)
    
    custom_fasta = None
    if not use_pretrained:
        custom_fasta = st.file_uploader("上传自定义参考蛋白 FASTA 文件", type=["fasta", "fa", "txt"], key="global_ref_fasta")
        if custom_fasta is not None:
            fasta_text = custom_fasta.getvalue().decode("utf-8")
            hash_val = compute_fasta_hash(fasta_text)
            cache_exists = os.path.exists(f"data/custom_cache_{hash_val}.npz")
            if cache_exists:
                st.info("💾 检测到缓存，构建时将直接加载。")
            if st.button("构建/更新自定义嵌入库"):
                print("[DEBUG] 开始构建自定义嵌入库...")
                with st.spinner("正在处理..."):
                    emb_dict, ids = build_reference_library(fasta_text, model, batch_converter, use_cache=True)
                    st.session_state.custom_library = emb_dict
                    st.session_state.custom_ids = ids
                    if len(emb_dict) > 0:
                        st.success(f"自定义库已就绪，共 {len(emb_dict)} 条蛋白。")
                    else:
                        st.error("构建失败。")
    else:
        if "custom_library" in st.session_state:
            del st.session_state.custom_library

    if use_pretrained:
        st.info("当前使用预构建 Swiss-Prot 库。")
    elif "custom_library" in st.session_state and st.session_state.custom_library:
        st.info(f"当前使用自定义库，包含 {len(st.session_state.custom_library)} 个蛋白。")
    else:
        st.info("尚未构建或选择任何参考库。")

    st.markdown("---")

    with st.expander("🔎 功能1：单条蛋白相似性检索 + AI摘要", expanded=True):
        st.markdown("输入一条蛋白序列或蛋白ID，在参考库中查找功能最相似的蛋白。")
        col1, col2 = st.columns(2)
        with col1:
            query_id = st.text_input("Master Protein ID", placeholder="Pt_Chr0100002")
            upload_query_fasta = st.file_uploader("上传查询蛋白 FASTA 文件（用于提取ID序列）", type=["fasta", "fa", "txt"], key="query_fasta_upload")
        with col2:
            query_seq = st.text_area("或直接粘贴蛋白序列", height=120, placeholder="MKGAKSK...")
        top_n = st.slider("返回结果数", 5, 50, 10, key="top_n_single")
        
        if st.button("🔍 搜索相似蛋白", key="single_search"):
            final_seq = None
            if query_seq.strip():
                lines = query_seq.strip().splitlines()
                seq_parts = [l for l in lines if not l.startswith(">")]
                final_seq = "".join(seq_parts)
                print(f"[DEBUG] 功能1: 使用粘贴序列，长度={len(final_seq)}")
            elif query_id.strip():
                upload_text = None
                if upload_query_fasta is not None:
                    upload_text = upload_query_fasta.getvalue().decode("utf-8")
                final_seq = get_sequence_for_id(query_id.strip(), fasta_text=upload_text)
                if final_seq is None:
                    st.error(f"无法获取蛋白 {query_id} 的序列。")
                    st.stop()
                print(f"[DEBUG] 功能1: 获取到序列，长度={len(final_seq)}")
            else:
                st.error("请至少输入蛋白ID或粘贴序列。")
                st.stop()
            
            custom_lib = st.session_state.get("custom_library", None)
            result_df = find_similar_proteins(
                final_seq, top_n=top_n,
                model=model, batch_converter=batch_converter,
                use_pretrained=use_pretrained,
                custom_library=custom_lib
            )
            if "Error" in result_df.columns:
                st.error(result_df["Error"][0])
            else:
                st.dataframe(result_df, width='stretch')
                st.session_state.last_search_result = result_df
                
                st.markdown("---")
                st.subheader("🤖 功能3：AI 功能摘要")
                st.caption("选择一个蛋白，生成其功能描述。")
                selected_idx = st.selectbox(
                    "选择蛋白：",
                    range(len(result_df)),
                    format_func=lambda i: f"#{i+1}: {result_df.iloc[i]['Protein_ID']} (相似度 {result_df.iloc[i]['Similarity']:.3f})",
                    key="summary_select"
                )
                use_api = st.checkbox("使用 DeepSeek AI（需 API Key）", value=False, key="use_api_summary")
                if st.button("生成摘要", key="gen_summary"):
                    row = result_df.iloc[selected_idx]
                    pid = row['Protein_ID']
                    sim = row['Similarity']
                    pname = row.get('Protein_Name', '')
                    gos = row.get('GO', '')
                    ecs = row.get('EC', '')
                    print(f"[DEBUG] 生成摘要: {pid}, use_api={use_api}")
                    with st.spinner("正在生成..."):
                        summary = generate_function_summary(
                            protein_id=pid,
                            protein_name=pname,
                            similarity_score=sim,
                            similar_protein=result_df.iloc[0]['Protein_ID'],
                            go_terms=gos,
                            ec_numbers=ecs,
                            use_api=use_api
                        )
                        st.session_state.current_summary = summary
                if 'current_summary' in st.session_state:
                    st.success(st.session_state.current_summary)

    st.markdown("---")

    with st.expander("🧬 功能2：批量AI注释 + 富集分析", expanded=False):
        st.markdown("输入多条蛋白序列，对每条序列找到参考库中最相似的蛋白，并对这些蛋白进行GO富集分析。")
        st.warning("富集分析要求参考库具有 UniProt 注释（标准 Accession），否则将跳过富集分析。")
        
        input_mode = st.radio("输入方式", ["粘贴序列", "上传FASTA文件"], horizontal=True, key="batch_mode")
        batch_seqs = []
        if input_mode == "粘贴序列":
            txt = st.text_area("每条序列一行（可包含FASTA头）", height=150, key="batch_text")
            if txt:
                lines = txt.strip().split('\n')
                seq = ""
                for line in lines:
                    if line.startswith('>'):
                        if seq:
                            batch_seqs.append(seq); seq = ""
                    else:
                        seq += line.strip()
                if seq:
                    batch_seqs.append(seq)
        else:
            batch_file = st.file_uploader("上传FASTA文件", type=["fasta", "fa", "txt"], key="batch_file")
            if batch_file:
                content = batch_file.read().decode('utf-8')
                batch_seqs = list(parse_fasta(content).values())
                st.success(f"解析到 {len(batch_seqs)} 条序列")
        
        if st.button("📥 为当前自定义库加载 UniProt 注释", disabled=(use_pretrained or "custom_library" not in st.session_state)):
            lib_ids = list(st.session_state.custom_library.keys())
            if lib_ids:
                sample_ids = lib_ids[:20]
                valid = sum(is_uniprot_accession(extract_accession(x)) for x in sample_ids)
                if valid < len(sample_ids) * 0.1:
                    st.warning("当前库的蛋白ID不是标准UniProt格式，无法获取注释。请使用标准Swiss-Prot蛋白测试富集。")
                else:
                    with st.spinner("获取注释中..."):
                        annots = fetch_uniprot_annotations_batch(lib_ids)
                        st.session_state.custom_annotations = annots
                        st.success(f"已获取 {len(annots)} 个蛋白的注释。")
            else:
                st.warning("参考库为空。")
        
        if st.button("🚀 运行批量注释与富集", disabled=(len(batch_seqs) == 0)):
            if "custom_library" not in st.session_state or not st.session_state.custom_library:
                st.error("请先构建参考库。")
            else:
                lib = st.session_state.custom_library
                with st.spinner(f"正在处理 {len(batch_seqs)} 条序列..."):
                    results = batch_search_top1(batch_seqs, model, batch_converter, lib)
                    target_ids = [r['best_id'] for r in results if r['best_id']]
                    st.subheader("📋 最佳匹配蛋白")
                    display_df = pd.DataFrame({
                        "查询序号": [r['query_index']+1 for r in results],
                        "最佳匹配蛋白": [r['best_id'] if r['best_id'] else '失败' for r in results],
                        "相似度": [round(r['similarity'], 4) for r in results]
                    })
                    st.dataframe(display_df, width='stretch')
                    
                    annotations = st.session_state.get('custom_annotations', {})
                    if annotations and len(target_ids) > 0:
                        st.subheader("📊 富集分析")
                        go_bg = get_go_background(annotations)
                        enrich_df = enrich_go(target_ids, annotations, go_bg)
                        if len(enrich_df) > 0:
                            fig = px.scatter(
                                enrich_df,
                                x='enrichment_ratio',
                                y='GO',
                                size='count',
                                color='adjusted_p',
                                color_continuous_scale='Reds_r',
                                title="富集的GO术语（超几何检验）",
                                labels={'enrichment_ratio': '富集因子', 'GO': 'GO术语', 'count': '目标基因数'}
                            )
                            st.plotly_chart(fig, width='stretch')
                            st.dataframe(enrich_df)
                        else:
                            st.info("未检测到显著富集的GO术语 (p<0.05)。")
                    else:
                        st.warning("缺少注释信息，无法进行富集分析。请先点击按钮加载注释。")

# ---------- STRING PPI Network (调试好的版本) ----------
elif page == "STRING PPI Network":
    st.header("STRING PPI Network")

    st.info(
        "STRING 网络通过在线 STRING API 获取。"
        "UniProt/EBI BLAST 也是在线服务，可能排队较久。"
        "推荐优先使用 eggNOG 注释或手动输入 UniProt ID。"
    )

    fasta_text = st.text_area(
        "FASTA 输入框：粘贴 FASTA 格式的蛋白序列",
        height=200,
        placeholder=">protein_1\nMKT...",
        help="这里填写 FASTA。只有选择 eggNOG 注释或 BLAST 映射时才需要。",
    )

    if fasta_text.strip():
        fasta_seqs = parse_fasta(fasta_text)
        original_ids = list(fasta_seqs.keys())
        st.write(f"解析到 {len(original_ids)} 条序列")
    else:
        fasta_seqs = {}
        original_ids = []

    st.subheader("映射到 UniProt / STRING")

    mapping_method = st.radio(
        "映射方法",
        [
            "优先使用 eggNOG 注释",
            "手动输入 UniProt ID",
            "使用序列 BLAST (UniProt，较慢)",
        ],
        index=1,
        help="推荐优先使用 eggNOG 注释或手动输入 UniProt ID；在线 BLAST 可能等待很久。",
    )

    manual_uniprot_text = ""
    if mapping_method == "手动输入 UniProt ID":
        manual_uniprot_text = st.text_area(
            "UniProt ID 输入框：每行一个 UniProt Accession",
            height=160,
            placeholder="P93004\nQ9S7W5\nA0A178WKB2",
            help="这里只能填写 UniProt Accession，不能填写 FASTA 序列。",
        )
        st.caption("示例格式：每行一个 accession，例如 P93004、Q9S7W5。不要粘贴 FASTA。")

    st.subheader("STRING 网络参数")

    species_options = {
        "拟南芥 Arabidopsis thaliana (3702)": 3702,
        "水稻 Oryza sativa japonica (39947)": 39947,
        "水稻 Oryza sativa (4530)": 4530,
        "人类 Homo sapiens (9606)": 9606,
    }

    species_label = st.selectbox(
        "选择物种",
        list(species_options.keys()),
        index=0,
    )
    species = species_options[species_label]

    required_score = st.slider(
        "相互作用阈值 required_score",
        0,
        1000,
        400,
    )

    if st.button("获取 STRING 网络"):
        id_mapping = {}

        with st.spinner("正在准备 UniProt ID..."):
            if mapping_method == "手动输入 UniProt ID":
                manual_ids = parse_manual_uniprot_ids(manual_uniprot_text)

                if not manual_ids:
                    st.error("请先输入 UniProt Accession。")
                    st.stop()

                id_mapping = {uid: uid for uid in manual_ids}

            elif mapping_method == "优先使用 eggNOG 注释":
                eggnog_df = st.session_state.get("eggnog_result", None)

                if eggnog_df is None:
                    st.error("尚未加载 eggNOG 注释。请先在 EggNOG 页面运行注释，或改用手动输入 UniProt ID。")
                    st.stop()

                if not original_ids:
                    st.warning("未粘贴 FASTA，将尝试从 eggNOG 表格中直接提取 ID。")
                    original_ids = list(eggnog_df.iloc[:, 0].astype(str))

                id_mapping = extract_uniprot_from_eggnog(eggnog_df, original_ids)

                if not id_mapping:
                    st.error("没有从 eggNOG 注释中提取到 UniProt ID。建议改用手动输入 UniProt ID。")
                    st.stop()

                st.success(f"从 eggNOG 提取到 {len(id_mapping)} 个映射")

            elif mapping_method == "使用序列 BLAST (UniProt，较慢)":
                if not fasta_seqs:
                    st.error("请先在 FASTA 输入框粘贴 FASTA 序列。")
                    st.stop()

                st.warning("在线 BLAST 可能排队较久，建议只测试少量序列。")
                id_mapping = run_blast_mapping(fasta_seqs)

                if not id_mapping:
                    st.error("BLAST 没有得到 UniProt 映射结果。")
                    st.stop()

        uniprot_ids = list(dict.fromkeys(id_mapping.values()))

        st.info(f"最终获得 {len(uniprot_ids)} 个 UniProt ID")

        map_df = pd.DataFrame(
            list(id_mapping.items()),
            columns=["Original_ID", "UniProt_ID"],
        )
        st.dataframe(map_df, use_container_width=True)

        with st.spinner("正在查询 STRING 在线接口..."):
            network_json = call_string_api(
                uniprot_ids,
                species=species,
                required_score=required_score,
            )

        if network_json is None:
            st.error("STRING API 调用失败，请检查网络、物种 ID 或 UniProt ID 是否匹配。")
            st.stop()

        if len(network_json) == 0:
            st.warning("STRING 返回为空。可能是物种选择不匹配，或这些蛋白之间没有满足阈值的互作。")
            st.stop()

        st.success(f"收到 STRING 网络数据，共 {len(network_json)} 条边")

        html = build_ppi_network_html(network_json, id_mapping)
        st.components.v1.html(html, height=800, scrolling=True)

        st.download_button(
            "下载原始网络 JSON",
            data=json.dumps(network_json, indent=2),
            file_name="string_network.json",
            mime="application/json",
        )
