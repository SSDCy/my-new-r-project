import numpy as np
import pandas as pd
from sklearn.impute import KNNImputer
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects.conversion import localconverter

def apply_missing_filter(data, lfq_cols, threshold, sample_info=None):
    """分组缺失过滤（如无样本信息则全局过滤）"""
    mat = data[lfq_cols]
    if sample_info is None or 'Group' not in sample_info.columns:
        # 全局过滤
        missing_frac = mat.isna().mean(axis=1)
        keep = missing_frac <= threshold
        return data.loc[keep]
    
    # 分组过滤
    # 构建样本到组的映射（通过短名匹配）
    short_to_group = {}
    for idx, row in sample_info.iterrows():
        short_name = idx.replace('LFQ intensity ', '').replace('.', '_')
        short_to_group[short_name] = row.get('Group', 'Unknown')
    
    sample_short = [c[len('LFQ intensity '):].replace('.', '_') for c in lfq_cols]
    group_vec = [short_to_group.get(s, 'Unknown') for s in sample_short]
    groups = list(set(group_vec))
    
    keep = np.zeros(len(data), dtype=bool)
    for g in groups:
        cols_in_group = [c for c, grp in zip(lfq_cols, group_vec) if grp == g]
        if not cols_in_group:
            continue
        miss = mat[cols_in_group].isna().mean(axis=1)
        keep |= (miss <= threshold).values
    return data.loc[keep]

def impute_missing(data, lfq_cols, method='quantile', k=10, min_value=1e-4, quantile_prob=0.01):
    """缺失值插补"""
    mat = data[lfq_cols].copy()
    if method == 'none':
        return data
    elif method == 'knn':
        imputer = KNNImputer(n_neighbors=k)
        imputed = imputer.fit_transform(mat)
        data[lfq_cols] = imputed
        return data
    elif method == 'minvalue':
        global_min = mat.min().min() if (mat > 0).any().any() else min_value
        fill_val = global_min * min_value
        data[lfq_cols] = mat.fillna(fill_val)
        return data
    elif method == 'quantile':
        for col in lfq_cols:
            q = mat[col].quantile(quantile_prob)
            data[col] = mat[col].fillna(q)
        return data
    else:
        raise ValueError(f"未知的插补方法: {method}")

def run_combat_py(data, lfq_cols, batch_vector):
    """通过 rpy2 调用 sva::ComBat"""
    with localconverter(ro.default_converter + pandas2ri.converter):
        # 准备 R 数据
        r_mat = ro.conversion.py2rpy(data[lfq_cols].T)  # R 需要基因×样本，我们转置
        ro.globalenv['expr_mat'] = r_mat
        ro.globalenv['batch'] = ro.StrVector(batch_vector)
    
    ro.r('''
    library(sva)
    combat_res <- ComBat(dat = as.matrix(expr_mat), batch = batch)
    ''')
    
    with localconverter(ro.default_converter + pandas2ri.converter):
        corrected = ro.globalenv['combat_res']
        corrected = pd.DataFrame(corrected, index=lfq_cols).T
    data[lfq_cols] = corrected
    return data
