# utils/data_io.py
import pandas as pd
import re
import streamlit as st
import numpy as np

def load_expression_data(file, clean=False):
    """
    读取 MaxQuant proteinGroups.txt。
    如果 clean=True，则进行完整清洗（移除 Reverse/Contaminant/CON_）。
    否则仅进行基本转换（数值化、零值替换、列重命名）。
    返回 (df, lfq_cols, sample_names, col_type, cleaning_stats, raw_df)
    """
    try:
        df = pd.read_csv(file, sep='\t', dtype=str)
        print(f"[DEBUG] load_expression_data: 原始数据形状={df.shape}")
    except Exception as e:
        st.error(f"无法读取文件: {e}")
        return None, None, None, None, None, None

    # 保存原始副本
    raw_df = df.copy()

    # 寻找强度列
    lfq_pattern = r'^LFQ intensity '
    intensity_pattern = r'^Intensity '
    lfq_cols = [c for c in df.columns if re.match(lfq_pattern, c)]
    intensity_cols = [c for c in df.columns if re.match(intensity_pattern, c)]

    if lfq_cols:
        sel_cols = lfq_cols
        col_type = 'LFQ'
    elif intensity_cols:
        sel_cols = intensity_cols
        col_type = 'Intensity'
    else:
        st.error("未找到强度列（LFQ intensity 或 Intensity）")
        return None, None, None, None, None, None

    print(f"[DEBUG] 检测到 {len(sel_cols)} 个 {col_type} 列")

    # 数值转换，零值替换为 NaN
    for c in sel_cols:
        df[c] = pd.to_numeric(df[c], errors='coerce')
    zero_count = (df[sel_cols] == 0).sum().sum()
    df[sel_cols] = df[sel_cols].replace(0, np.nan)
    print(f"[DEBUG] 将 {zero_count} 个零值替换为 NaN")

    # 清洗统计
    original_rows = df.shape[0]
    reverse_removed = 0
    contam_removed = 0
    con_removed = 0

    if clean:
        # 只在 clean=True 时才过滤
        if 'Reverse' in df.columns:
            reverse_mask = df['Reverse'].str.contains(r'\+', na=False)
            reverse_removed = reverse_mask.sum()
            df = df[~reverse_mask]
            print(f"[DEBUG] 移除 {reverse_removed} 行 Reverse hits")
        if 'Potential contaminant' in df.columns:
            contam_mask = df['Potential contaminant'].str.contains(r'\+', na=False)
            contam_removed = contam_mask.sum()
            df = df[~contam_mask]
            print(f"[DEBUG] 移除 {contam_removed} 行 Potential contaminants")
        if 'Protein IDs' in df.columns:
            con_mask = df['Protein IDs'].str.startswith('CON_', na=False)
            con_removed = con_mask.sum()
            df = df[~con_mask]
            print(f"[DEBUG] 移除 {con_removed} 行 CON_ contaminants")
        # 生成 Master protein IDs（仅在过滤后需要）
        if 'Protein IDs' in df.columns:
            df['Master protein IDs'] = df['Protein IDs'].str.split(';').str[0]
    else:
        # 即使不过滤，也生成 Master protein IDs 以便后续使用
        if 'Protein IDs' in df.columns:
            df['Master protein IDs'] = df['Protein IDs'].str.split(';').str[0]

    # 标准化列名：点转下划线
    rename_dict = {}
    for c in sel_cols:
        new_c = c.replace('.', '_')
        if new_c != c:
            rename_dict[c] = new_c
    if rename_dict:
        df.rename(columns=rename_dict, inplace=True)
        sel_cols = [rename_dict.get(c, c) for c in sel_cols]

    # 提取短样本名
    prefix = 'LFQ intensity ' if col_type == 'LFQ' else 'Intensity '
    sample_names = []
    for c in sel_cols:
        name = c[len(prefix):].replace('.', '_')
        sample_names.append(name)

    final_rows = df.shape[0]
    cleaning_stats = {
        'original': original_rows,
        'reverse_removed': reverse_removed,
        'contaminant_removed': contam_removed,
        'con_removed': con_removed,
        'retained': final_rows
    }
    print(f"[DEBUG] 清洗后保留 {final_rows} 行, {len(sel_cols)} 个强度列, 缺失值总数 {df[sel_cols].isna().sum().sum()}")
    return df, sel_cols, sample_names, col_type, cleaning_stats, raw_df


def load_sample_info(file):
    """读取样本信息表，支持 CSV 和 Excel"""
    try:
        filename = file.name.lower()
        if filename.endswith('.csv'):
            si = pd.read_csv(file, index_col=0)
        elif filename.endswith(('.xlsx', '.xls')):
            si = pd.read_excel(file, index_col=0)
        else:
            st.error("不支持的文件格式，请上传 CSV 或 Excel 文件")
            return None
        # 标准化行名（点转下划线）
        si.index = si.index.astype(str).str.replace('.', '_', regex=False)
        print(f"[DEBUG] 样本信息加载成功: {si.shape}")
        return si
    except Exception as e:
        st.error(f"读取样本信息失败: {e}")
        return None
