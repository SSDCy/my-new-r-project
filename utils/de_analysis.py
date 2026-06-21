import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, mannwhitneyu
import subprocess, tempfile, os
import streamlit as st

def run_de_analysis(data, treat_cols, ctrl_cols, fc_up=1.2, fc_down=0.84, p_cut=0.05,
                    method='t-test',
                    min_treat_valid=2, min_ctrl_valid=2,
                    min_rep_ttest=2, min_rep_inc=2, min_rep_dec=2):
    """
    差异分析，完全对齐原 R 代码逻辑。
    返回结果包含分类 regulation 及注释 regulation_note。
    """
    treat_mat = data[treat_cols].values.astype(float)
    ctrl_mat = data[ctrl_cols].values.astype(float)

    # 计算有效值（非 NaN）数量
    n_treat = np.sum(~np.isnan(treat_mat), axis=1)
    n_ctrl = np.sum(~np.isnan(ctrl_mat), axis=1)

    # 均值
    mean_treat = np.nanmean(treat_mat, axis=1)
    mean_ctrl = np.nanmean(ctrl_mat, axis=1)

    # Fold Change
    with np.errstate(divide='ignore', invalid='ignore'):
        fc = mean_treat / mean_ctrl
        log2fc = np.log2(fc)

    # P 值初始化
    p_values = np.full(len(data), np.nan)

    if method == 'limma':
        try:
            r_res = _run_limma_r(data, treat_cols, ctrl_cols)
            if r_res is not None:
                p_values = r_res['P.Value'].reindex(data.index).values
                log2fc = r_res['logFC'].reindex(data.index).values
                fc = 2 ** log2fc
        except Exception as e:
            st.error(f"limma 失败: {e}，回退到 t-test")
            method = 't-test'

    if method in ['t-test', 'wilcoxon']:
        for i in range(len(data)):
            t_vals = treat_mat[i, ~np.isnan(treat_mat[i,:])]
            c_vals = ctrl_mat[i, ~np.isnan(ctrl_mat[i,:])]
            # 必须两组都满足最小有效值要求，且至少各有2个观测才做检验（R 代码中逻辑相同）
            if len(t_vals) < min_treat_valid or len(c_vals) < min_ctrl_valid:
                continue
            if len(t_vals) < 2 or len(c_vals) < 2:
                continue
            try:
                if method == 't-test':
                    _, p = ttest_ind(t_vals, c_vals, equal_var=True)
                else:
                    _, p = mannwhitneyu(t_vals, c_vals, alternative='two-sided')
                p_values[i] = p
            except:
                pass

    # ---------- 分类规则（严格对齐 R） ----------
    regulation = []
    regulation_note = []
    for i in range(len(data)):
        nt = n_treat[i]
        nc = n_ctrl[i]
        fc_i = fc[i]
        p_i = p_values[i]

        # 情况1：处理组有值、对照组完全无值 → Increase（前提：处理组有效值 >= min_rep_inc）
        if nt > 0 and nc == 0 and nt >= min_rep_inc:
            regulation.append('Increase')
            regulation_note.append(f'仅处理组检测到 (n≥{min_rep_inc})')
        # 情况2：对照组有值、处理组完全无值 → Decrease（前提：对照组有效值 >= min_rep_dec）
        elif nc > 0 and nt == 0 and nc >= min_rep_dec:
            regulation.append('Decrease')
            regulation_note.append(f'仅对照组检测到 (n≥{min_rep_dec})')
        # 情况3：两组都有值，且满足检验所需最小重复数
        elif nt >= min_rep_ttest and nc >= min_rep_ttest:
            if np.isnan(p_i):
                regulation.append('NS')
                regulation_note.append('统计检验不可用')
            elif fc_i > fc_up and p_i < p_cut:
                regulation.append('Up')
                regulation_note.append(f'显著上调 (FC>{fc_up:.2f}, p<{p_cut})')
            elif fc_i < fc_down and p_i < p_cut:
                regulation.append('Down')
                regulation_note.append(f'显著下调 (FC<{fc_down:.2f}, p<{p_cut})')
            else:
                regulation.append('NS')
                regulation_note.append('未显著')
        else:
            # 重复数不足
            regulation.append('NS')
            reasons = []
            if nt < min_rep_ttest:
                reasons.append(f'处理组重复数 {nt}<{min_rep_ttest}')
            if nc < min_rep_ttest:
                reasons.append(f'对照组重复数 {nc}<{min_rep_ttest}')
            regulation_note.append('重复数不足: ' + '; '.join(reasons))

    # 为 Increase/Decrease 生成模拟坐标（仅用于火山图可视化，不影响计数）
    inc_idx = [i for i, reg in enumerate(regulation) if reg == 'Increase']
    dec_idx = [i for i, reg in enumerate(regulation) if reg == 'Decrease']
    if inc_idx:
        log2fc[inc_idx] = np.random.uniform(5.5, 7.5, size=len(inc_idx))
        p_values[inc_idx] = np.random.uniform(1e-6, 1e-2, size=len(inc_idx))
        print(f"[DEBUG] 为 {len(inc_idx)} 个 Increase 蛋白生成了模拟坐标")
    if dec_idx:
        log2fc[dec_idx] = np.random.uniform(-7.5, -5.5, size=len(dec_idx))
        p_values[dec_idx] = np.random.uniform(1e-6, 1e-2, size=len(dec_idx))
        print(f"[DEBUG] 为 {len(dec_idx)} 个 Decrease 蛋白生成了模拟坐标")

    log10p = -np.log10(np.clip(p_values, 1e-10, None))

    result = data[['Master protein IDs']].copy()
    result['n_treat'] = n_treat
    result['n_ctrl'] = n_ctrl
    result['mean_treat'] = mean_treat
    result['mean_ctrl'] = mean_ctrl
    result['FC'] = fc
    result['log2FC'] = log2fc
    result['Pvalue'] = p_values
    result['-log10P'] = log10p
    result['regulation'] = regulation
    result['regulation_note'] = regulation_note

    # 输出汇总调试信息
    total = len(result)
    up = (result['regulation']=='Up').sum()
    down = (result['regulation']=='Down').sum()
    inc = (result['regulation']=='Increase').sum()
    dec = (result['regulation']=='Decrease').sum()
    ns = (result['regulation']=='NS').sum()
    print(f"[DEBUG] 差异分析结果: 总蛋白 {total}, Up {up}, Down {down}, Increase {inc}, Decrease {dec}, NS {ns}")
    return result


def _run_limma_r(data, treat_cols, ctrl_cols):
    """通过 Rscript 调用 limma"""
    expr = pd.concat([data[treat_cols], data[ctrl_cols]], axis=1)
    with tempfile.NamedTemporaryFile(suffix='.csv', delete=False) as f:
        expr.to_csv(f.name)
        input_path = f.name
    input_path_clean = input_path.replace('\\', '/')
    r_code = f"""
    library(limma)
    df <- read.csv("{input_path_clean}", row.names=1)
    mat <- as.matrix(df)
    group <- factor(c(rep("Treat", {len(treat_cols)}), rep("Ctrl", {len(ctrl_cols)})), levels=c("Ctrl","Treat"))
    design <- model.matrix(~group)
    fit <- lmFit(mat, design)
    fit <- eBayes(fit)
    res <- topTable(fit, coef="groupTreat", number=Inf, sort.by="none")
    write.csv(res, "{input_path_clean}.limma.csv", row.names=TRUE)
    """
    with tempfile.NamedTemporaryFile(suffix='.R', delete=False, mode='w') as f:
        f.write(r_code)
        r_script_path = f.name
    try:
        subprocess.run(["Rscript", r_script_path], check=True, capture_output=True, text=True,
                       encoding='utf-8', errors='replace')
        result = pd.read_csv(input_path + '.limma.csv', index_col=0)
        os.unlink(input_path)
        os.unlink(input_path + '.limma.csv')
        os.unlink(r_script_path)
        return result
    except Exception as e:
        print(f"limma 错误: {e}")
        return None
