import pandas as pd
import numpy as np
import streamlit as st

def total_intensity_normalize(data, lfq_cols, baseline_sample=None, sample_names=None):
    """总强度归一化，异常数据自动跳过"""
    print(f"[DEBUG] total_intensity_normalize: 基线样本={baseline_sample}")
    mat = data[lfq_cols].copy()
    missing_before = mat.isna().sum().sum()
    print(f"[DEBUG] 归一化前缺失值总数: {missing_before}")
    
    # 检查数据完整性
    if mat.sum().sum() == 0:
        st.error("所有样本总强度为零，可能数据异常，跳过归一化")
        print("[ERROR] 所有样本总强度为零，跳过归一化")
        return data
    
    # 确定基线样本
    if baseline_sample is None:
        baseline_sample = lfq_cols[0]
    
    actual_baseline_col = None
    if baseline_sample in lfq_cols:
        actual_baseline_col = baseline_sample
    else:
        short_to_full = {n: c for n, c in zip(sample_names, lfq_cols) if n and c}
        actual_baseline_col = short_to_full.get(baseline_sample)
    
    if actual_baseline_col is None:
        raise ValueError(f"基线样本 {baseline_sample} 未找到")
    
    baseline_total = mat[actual_baseline_col].sum()
    print(f"[DEBUG] 基线样本 {actual_baseline_col} 总强度: {baseline_total}")
    
    if baseline_total == 0:
        st.warning(f"基线样本 {baseline_sample} 总强度为零，正在自动选择其他样本...")
        for col in lfq_cols:
            total = mat[col].sum()
            if total > 0:
                actual_baseline_col = col
                baseline_total = total
                st.info(f"已自动切换基线样本为: {col}")
                print(f"[DEBUG] 自动切换基线样本为: {col}, 总强度: {baseline_total}")
                break
        if baseline_total == 0:
            st.error("所有样本总强度均为零，无法归一化，保留原始数据")
            print("[ERROR] 所有样本总强度为零，归一化跳过")
            return data
    
    # 执行归一化
    norm_mat = mat.div(mat.sum(axis=0), axis=1) * baseline_total
    norm_cols = ['Norm_' + c for c in lfq_cols]
    for orig_col, norm_col in zip(lfq_cols, norm_cols):
        data[norm_col] = norm_mat[orig_col]
    
    missing_after = data[norm_cols].isna().sum().sum()
    print(f"[DEBUG] 归一化后 Norm_ 列缺失值总数: {missing_after}")
    return data
