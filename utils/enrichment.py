# utils/enrichment.py
"""
GO 富集分析模块
- 超几何检验
- Benjamini-Hochberg FDR 校正
- 生成气泡图数据
"""

import numpy as np
from scipy.stats import hypergeom
import pandas as pd
from typing import Dict, List          # 修复：导入类型提示

DEBUG = True
def debug(msg: str):
    if DEBUG:
        print(f"[ENRICHMENT DEBUG] {msg}")

def get_go_background(annotation_dict: Dict) -> Dict[str, int]:
    """
    从注释字典中统计每个 GO 术语在背景（全体参考蛋白）中的出现次数。
    annotation_dict: {protein_id: {'go': [...], 'ec': [...]}, ...}
    返回：{GO_term: frequency}
    """
    debug("计算背景 GO 频率...")
    go_counts = {}
    for acc, ann in annotation_dict.items():
        for go in ann.get('go', []):
            go_counts[go] = go_counts.get(go, 0) + 1
    debug(f"背景 GO 术语总数：{len(go_counts)}")
    return go_counts

def enrich_go(target_ids: List[str], annotation_dict: Dict, go_background: Dict[str, int],
              p_cutoff: float = 0.05) -> pd.DataFrame:
    """
    对目标 ID 列表进行 GO 富集分析（超几何检验）。
    参数：
        target_ids: 目标蛋白 ID 列表（如 Top1 匹配结果）
        annotation_dict: 所有蛋白的注释字典
        go_background: 背景 GO 频率
        p_cutoff: 显著性阈值
    返回：
        DataFrame，包含 GO, p_value, adjusted_p, count, background_count, enrichment_ratio
    """
    debug(f"开始富集分析，目标蛋白数={len(target_ids)}，背景蛋白数={len(annotation_dict)}")
    M = len(annotation_dict)          # 背景总蛋白数
    N = len(target_ids)               # 目标蛋白数

    # 统计目标中每个 GO 出现次数
    target_go_counts = {}
    for tid in target_ids:
        if tid in annotation_dict:
            for go in annotation_dict[tid].get('go', []):
                target_go_counts[go] = target_go_counts.get(go, 0) + 1

    results = []
    for go, k in target_go_counts.items():
        n = go_background.get(go, 0)   # 背景中该 GO 的总数
        if n == 0:
            continue
        # 超几何检验：在 M 个蛋白中随机抽取 N 个，至少观测到 k 个的概率
        p_val = hypergeom.sf(k - 1, M, n, N)
        if p_val < p_cutoff:
            results.append({
                'GO': go,
                'p_value': p_val,
                'count': k,
                'background_count': n,
                'enrichment_ratio': (k / N) / (n / M)
            })

    df = pd.DataFrame(results)
    if len(df) > 0:
        # Benjamini-Hochberg 校正
        df['adjusted_p'] = p_adjust_bh(df['p_value'].values)
        df = df.sort_values('p_value')
    debug(f"富集分析完成，显著 GO 数量：{len(df)}")
    return df

def p_adjust_bh(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR 校正"""
    pvals = np.array(pvals, dtype=float)
    n = len(pvals)
    if n == 0:
        return np.array([])
    sorted_idx = np.argsort(pvals)
    sorted_p = pvals[sorted_idx]
    adj_p = np.empty(n)
    for i, p in enumerate(sorted_p):
        adj_p[sorted_idx[i]] = min(p * n / (i + 1), 1.0)
    # 确保单调性
    for i in range(n - 2, -1, -1):
        adj_p[sorted_idx[i]] = min(adj_p[sorted_idx[i]], adj_p[sorted_idx[i + 1]])
    return adj_p
