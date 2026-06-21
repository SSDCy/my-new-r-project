import numpy as np
import pandas as pd
from sklearn.impute import KNNImputer
import subprocess, tempfile, os, time
import streamlit as st

# ==================== 缺失过滤 ====================
def apply_missing_filter(data, lfq_cols, threshold, sample_info=None):
    """分组缺失过滤"""
    print(f"[DEBUG] apply_missing_filter: 输入数据形状={data.shape}, 阈值={threshold}")
    mat = data[lfq_cols]
    if sample_info is None or 'Group' not in sample_info.columns:
        missing_frac = mat.isna().mean(axis=1)
        keep = missing_frac <= threshold
        print(f"[DEBUG] 全局过滤: 保留 {keep.sum()} 行")
        return data.loc[keep]
    
    short_to_group = {}
    for idx, row in sample_info.iterrows():
        short_name = idx.replace('LFQ intensity ', '').replace('.', '_')
        short_to_group[short_name] = row.get('Group', 'Unknown')
    
    sample_short = [c[len('LFQ intensity '):].replace('.', '_') for c in lfq_cols]
    group_vec = [short_to_group.get(s, 'Unknown') for s in sample_short]
    groups = list(set(group_vec))
    print(f"[DEBUG] 分组过滤: 检测到组 {groups}")
    
    keep = np.zeros(len(data), dtype=bool)
    for g in groups:
        cols_in_group = [c for c, grp in zip(lfq_cols, group_vec) if grp == g]
        if not cols_in_group:
            continue
        miss = mat[cols_in_group].isna().mean(axis=1)
        keep |= (miss <= threshold).values
        print(f"[DEBUG] 组 '{g}': {keep.sum()} 行通过")
    print(f"[DEBUG] 分组过滤: 最终保留 {keep.sum()} 行")
    return data.loc[keep]

# ==================== 缺失值插补 ====================
def impute_missing(data, lfq_cols, method='quantile', k=10, min_value=1e-4, quantile_prob=0.01):
    """缺失值插补"""
    print(f"[DEBUG] impute_missing: 方法={method}, 数据形状={data.shape}")
    missing_before = data[lfq_cols].isna().sum().sum()
    print(f"[DEBUG] 插补前缺失值总数: {missing_before}")
    
    if method == 'none':
        print("[DEBUG] 跳过插补")
        return data
    elif method == 'knn':
        print(f"[DEBUG] KNN插补: k={k}")
        imputer = KNNImputer(n_neighbors=k)
        imputed = imputer.fit_transform(data[lfq_cols])
        data[lfq_cols] = imputed
    elif method == 'ppca':
        # 先使用 KNN 预填充，确保无缺失，然后 PPCA 优化
        print("[DEBUG] PPCA 插补：先执行 KNN 预填充")
        imputer = KNNImputer(n_neighbors=k)
        knn_imputed = imputer.fit_transform(data[lfq_cols])
        data[lfq_cols] = knn_imputed
        try:
            data = _impute_ppca_safe(data, lfq_cols)
        except Exception as e:
            print(f"[WARNING] PPCA 优化失败: {e}，将保留 KNN 结果")
            st.warning("PPCA 优化失败，保留 KNN 插补结果")
    elif method == 'minvalue':
        data = _impute_minvalue(data, lfq_cols, min_value)
    elif method == 'quantile':
        print(f"[DEBUG] 分位数插补: q={quantile_prob}")
        mat = data[lfq_cols].copy()
        for col in lfq_cols:
            q = mat[col].quantile(quantile_prob)
            data[col] = mat[col].fillna(q)
    else:
        raise ValueError(f"未知的插补方法: {method}")
    
    missing_after = data[lfq_cols].isna().sum().sum()
    print(f"[DEBUG] 插补后缺失值总数: {missing_after}")
    if missing_after > 0:
        print("[WARNING] 最终仍存在缺失值，请检查数据")
    else:
        print("[DEBUG] 所有缺失值已填充")
    return data

def _impute_minvalue(data, lfq_cols, factor):
    """最低值填充"""
    print(f"[DEBUG] minvalue: factor={factor}")
    mat = data[lfq_cols].copy()
    all_vals = mat.values.flatten()
    valid_vals = all_vals[(~np.isnan(all_vals)) & (all_vals > 0)]
    if len(valid_vals) == 0:
        global_min = 1e-4
        print("[WARNING] minvalue: 无有效非零值，使用默认最小值 1e-4")
    else:
        global_min = np.min(valid_vals)
        print(f"[DEBUG] minvalue: 全局最小非零值={global_min}")
    fill_val = global_min * factor
    print(f"[DEBUG] minvalue: 填充值={fill_val}")
    
    for col in lfq_cols:
        before_na = data[col].isna().sum()
        data[col] = data[col].fillna(fill_val)
        after_na = data[col].isna().sum()
        if before_na > 0:
            print(f"[DEBUG] 列 '{col}': 填充前缺失 {before_na}, 填充后缺失 {after_na}")
    return data

def _impute_ppca_safe(data, lfq_cols):
    """
    PPCA 优化：输入数据已经过 KNN 填充，无缺失值。
    调用 R 的 pcaMethods::ppca，仅用于改善分布，不引入新缺失。
    """
    print("[DEBUG] 执行 PPCA 优化（输入已无缺失）")
    expr = data[lfq_cols].T   # 样本×蛋白
    expr.index.name = 'Sample'
    
    # 保存临时 CSV
    expr_file = tempfile.NamedTemporaryFile(suffix='.csv', delete=False)
    expr.to_csv(expr_file.name)
    expr_file.close()
    expr_path = expr_file.name.replace('\\', '/')
    
    out_file = tempfile.NamedTemporaryFile(suffix='.csv', delete=False)
    out_file.close()
    out_path = out_file.name.replace('\\', '/')
    
    # R 代码：直接对完整矩阵执行 PPCA
    r_code = f"""
    library(pcaMethods)
    mat <- as.matrix(read.csv("{expr_path}", row.names=1))
    # 矩阵已经无缺失，但以防万一
    if(any(is.na(mat))) mat[is.na(mat)] <- 0
    # PPCA
    pc <- pca(mat, method="ppca", nPcs=min(2, ncol(mat)), scale="uv", center=TRUE)
    imputed <- completeObs(pc)
    write.csv(imputed, "{out_path}", row.names=TRUE)
    """
    
    r_script = tempfile.NamedTemporaryFile(suffix='.R', delete=False, mode='w')
    r_script.write(r_code)
    r_script.close()
    r_script_path = r_script.name
    
    try:
        start_time = time.time()
        result = subprocess.run(
            ["Rscript", r_script_path], check=True, capture_output=True, text=True,
            encoding='utf-8', errors='replace'
        )
        elapsed = time.time() - start_time
        print(f"[DEBUG] PPCA R 执行成功，耗时 {elapsed:.2f} 秒")
        if result.stderr:
            print(f"[DEBUG] R 输出:\n{result.stderr[:500]}")
        
        corrected = pd.read_csv(out_file.name, index_col=0).T
        # 清理临时文件
        for f in [expr_file.name, out_file.name, r_script_path]:
            if os.path.exists(f):
                os.unlink(f)
        
        # 列对齐：只更新存在的列
        common_cols = [c for c in lfq_cols if c in corrected.columns]
        if len(common_cols) != len(lfq_cols):
            print(f"[WARNING] PPCA 结果列数不一致，保留 KNN 结果")
            # 不更新 data，返回即可
            return data
        
        # 检查是否引入了异常值（全零行等）
        if corrected[common_cols].sum().sum() == 0:
            print("[WARNING] PPCA 结果全为零，丢弃，保留 KNN 结果")
            return data
        
        data[common_cols] = corrected[common_cols].values
        print("[DEBUG] PPCA 优化完成")
        return data
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] PPCA R 调用失败: {e.stderr}")
        for f in [expr_file.name, out_file.name, r_script_path]:
            if os.path.exists(f):
                os.unlink(f)
        raise RuntimeError(f"Rscript 失败: {e.stderr}")

# ==================== 批次校正 ====================
def run_combat_py(data, lfq_cols, batch_vector):
    """ComBat 包装"""
    return _run_combat_r(data, lfq_cols, batch_vector)

def _run_combat_r(data, lfq_cols, batch_vector):
    """通过 Rscript 调用 sva::ComBat"""
    print("[DEBUG] 开始 ComBat 批次校正")
    expr = data[lfq_cols].copy().T
    expr.index.name = 'Sample'
    
    expr_file = tempfile.NamedTemporaryFile(suffix='.csv', delete=False)
    expr.to_csv(expr_file.name)
    expr_file.close()
    expr_path = expr_file.name.replace('\\', '/')

    batch_file = tempfile.NamedTemporaryFile(suffix='.txt', delete=False, mode='w')
    batch_file.write('\n'.join(batch_vector))
    batch_file.close()
    batch_path = batch_file.name.replace('\\', '/')

    out_file = tempfile.NamedTemporaryFile(suffix='.csv', delete=False)
    out_file.close()
    out_path = out_file.name.replace('\\', '/')

    r_code = f"""
    library(sva)
    expr <- as.matrix(read.csv("{expr_path}", row.names=1))
    batch <- readLines("{batch_path}")
    corrected <- ComBat(dat=expr, batch=batch)
    write.csv(corrected, "{out_path}", row.names=TRUE)
    """

    r_script = tempfile.NamedTemporaryFile(suffix='.R', delete=False, mode='w')
    r_script.write(r_code)
    r_script.close()
    r_script_path = r_script.name

    try:
        subprocess.run(["Rscript", r_script_path], check=True, capture_output=True, text=True,
                       encoding='utf-8', errors='replace')
        corrected = pd.read_csv(out_file.name, index_col=0).T
        for f in [expr_file.name, batch_file.name, out_file.name, r_script_path]:
            if os.path.exists(f):
                os.unlink(f)
        for col in lfq_cols:
            if col in corrected.columns:
                data[col] = corrected[col]
        print("[DEBUG] ComBat 完成")
        return data
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] ComBat 失败: {e.stderr}")
        st.error(f"ComBat 失败: {e.stderr}")
        return data
