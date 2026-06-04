# server/data_upload_quality.R
message("[DEBUG] data_upload_quality.R loading...")

# ================== 数据质量分析辅助函数 ==================
calculate_data_quality_score <- function(expr_matrix) {
  message("[DEBUG] calculate_data_quality_score: starting, dim = ", 
          if (!is.null(expr_matrix)) paste(nrow(expr_matrix), "x", ncol(expr_matrix)) else "NULL")
  if (is.null(expr_matrix) || nrow(expr_matrix) == 0 || ncol(expr_matrix) == 0) {
    message("[DEBUG] calculate_data_quality_score: empty or NULL matrix, returning score=0")
    return(list(score = 0, grade = "Poor", details = list()))
  }
  
  total_values <- nrow(expr_matrix) * ncol(expr_matrix)
  missing_values <- sum(is.na(expr_matrix))
  missing_ratio <- missing_values / total_values
  missing_score <- max(0, 30 * (1 - missing_ratio * 2))
  message(sprintf("[DEBUG] calculate_data_quality_score: missing_ratio=%.3f, missing_score=%.1f", missing_ratio, missing_score))
  
  original_na <- sum(is.na(expr_matrix))
  message("[DEBUG] calculate_data_quality_score: original NA count = ", original_na)
  
  filled <- tryCatch({
    message("[DEBUG] calculate_data_quality_score: attempting KNN imputation for correlation")
    impute_missing_values(as.data.frame(expr_matrix), method = "knn")
  }, error = function(e) {
    message("[DEBUG] calculate_data_quality_score: KNN failed (", e$message, "), using simple min fill")
    expr_matrix[is.na(expr_matrix)] <- 1e-4
    return(expr_matrix)
  })
  filled <- as.matrix(filled)
  after_na <- sum(is.na(filled))
  message("[DEBUG] calculate_data_quality_score: after imputation, NA count = ", after_na)
  
  sample_cor <- cor(filled, use = "complete.obs")
  diag(sample_cor) <- NA
  avg_cor <- mean(sample_cor, na.rm = TRUE)
  message(sprintf("[DEBUG] calculate_data_quality_score: avg_cor (post-imputation) = %.3f", avg_cor))
  
  consistency_score <- max(0, 40 * pmin(1, avg_cor / 0.8))
  message(sprintf("[DEBUG] calculate_data_quality_score: consistency_score=%.1f", consistency_score))
  
  protein_valid <- rowSums(!is.na(expr_matrix)) >= 2
  protein_valid_ratio <- mean(protein_valid)
  protein_score <- max(0, 30 * protein_valid_ratio)
  message(sprintf("[DEBUG] calculate_data_quality_score: valid_ratio=%.3f, protein_score=%.1f", protein_valid_ratio, protein_score))
  
  total_score <- round(missing_score + consistency_score + protein_score, 1)
  
  grade <- if (total_score >= 90) "Excellent"
  else if (total_score >= 80) "Good"
  else if (total_score >= 60) "Fair"
  else "Poor"
  
  details <- list(
    missing_ratio = round(missing_ratio * 100, 2),
    missing_score = round(missing_score, 1),
    avg_correlation = round(avg_cor, 3),
    consistency_score = round(consistency_score, 1),
    protein_valid_ratio = round(protein_valid_ratio * 100, 2),
    protein_score = round(protein_score, 1)
  )
  message("[DEBUG] calculate_data_quality_score: total=", total_score, ", grade=", grade)
  list(score = total_score, grade = grade, details = details)
}

generate_quality_report <- function(quality_score, expr_matrix, sample_info = NULL) {
  message("[DEBUG] generate_quality_report: starting report generation")
  details <- quality_score$details
  key_findings <- list()
  recommendations <- list()
  special_note <- ""
  
  # 缺失率评估
  if (details$missing_ratio < 10) {
    key_findings <- c(key_findings, list(list(
      type = "success",
      title = "✅ 低缺失率",
      content = paste0("数据缺失率仅为 ", details$missing_ratio, "%，完整性非常好，可直接进行下游分析。\n→ 对应预处理操作：无需特殊处理，或设置缺失值过滤阈值为 0.3 进行轻度过滤")
    )))
  } else if (details$missing_ratio < 20) {
    key_findings <- c(key_findings, list(list(
      type = "success",
      title = "✅ 较低缺失率",
      content = paste0("数据缺失率为 ", details$missing_ratio, "%，整体完整性较好，不影响主要分析。\n→ 对应预处理操作：建议设置缺失值过滤阈值为 0.5，并执行 KNN 填充")
    )))
  } else if (details$missing_ratio < 40) {
    key_findings <- c(key_findings, list(list(
      type = "warning",
      title = "⚠️ 中等缺失率",
      content = paste0("数据缺失率为 ", details$missing_ratio, "%，会降低部分统计方法的效力。\n→ 对应预处理操作：设置缺失值过滤阈值为 0.5，并执行 KNN 或 PPCA 填充")
    )))
  } else {
    key_findings <- c(key_findings, list(list(
      type = "danger",
      title = "高缺失率",
      content = paste0("数据缺失率高达 ", details$missing_ratio, "%，会直接影响后续差异分析、聚类的可靠性。\n→ 对应预处理操作：缺失值过滤（阈值 0.5）+ 缺失值填充（KNN/PPCA）")
    )))
  }
  
  # 样本一致性
  if (details$avg_correlation > 0.9) {
    key_findings <- c(key_findings, list(list(
      type = "success", title = "✅ 优秀的样本一致性",
      content = paste0("样本间平均相关性为 ", round(details$avg_correlation, 3), "，实验重复性极佳。\n→ 对应预处理操作：无需额外处理")
    )))
  } else if (details$avg_correlation > 0.8) {
    key_findings <- c(key_findings, list(list(
      type = "success", title = "✅ 良好的样本一致性",
      content = paste0("样本间平均相关性为 ", round(details$avg_correlation, 3), "，实验重复性较好。\n→ 对应预处理操作：可考虑执行批次校正以进一步提升一致性")
    )))
  } else if (details$avg_correlation > 0.7) {
    key_findings <- c(key_findings, list(list(
      type = "warning", title = "⚠️ 一般的样本一致性",
      content = paste0("样本间平均相关性为 ", round(details$avg_correlation, 3), "，可能存在技术波动。\n→ 对应预处理操作：检查并移除异常样本，执行缺失值填充")
    )))
  } else {
    key_findings <- c(key_findings, list(list(
      type = "danger", title = "较差的样本一致性",
      content = paste0("样本间平均相关性仅为 ", round(details$avg_correlation, 3), "，样本重复性较差，会降低后续分析的统计效力。\n→ 对应预处理操作：异常样本移除 + 缺失值填充")
    )))
  }
  
  # 蛋白有效检出率
  if (details$protein_valid_ratio > 80) {
    key_findings <- c(key_findings, list(list(
      type = "success", title = "蛋白有效检出率高",
      content = paste0("超过 ", details$protein_valid_ratio, "% 的蛋白在至少 2 个样本中被检出，蛋白整体质量优秀。\n→ 对应预处理操作：最小强度过滤阈值可设为 100000，无需大幅调整")
    )))
  } else if (details$protein_valid_ratio > 60) {
    key_findings <- c(key_findings, list(list(
      type = "warning", title = "⚠️ 蛋白有效检出率一般",
      content = paste0("约 ", details$protein_valid_ratio, "% 的蛋白在至少 2 个样本中被检出，部分蛋白检出率偏低。\n→ 对应预处理操作：适当降低强度过滤阈值，保留更多蛋白")
    )))
  } else {
    key_findings <- c(key_findings, list(list(
      type = "danger", title = "ⓧ 蛋白有效检出率低",
      content = paste0("仅有 ", details$protein_valid_ratio, "% 的蛋白在至少 2 个样本中被检出，数据可靠性较低。\n→ 对应预处理操作：建议放宽强度过滤阈值，并核对实验流程")
    )))
  }
  
  # PCA 相关的深入分析
  message("[DEBUG] generate_quality_report: running PCA analysis...")
  pca_result <- calculate_pca(expr_matrix, sample_info)
  if (!is.null(pca_result)) {
    pca_df <- pca_result$pca_df
    pc1 <- pca_df$PC1
    pc2 <- pca_df$PC2
    samples_pca <- pca_df$Sample
    
    # 批次效应检测
    if (!is.null(sample_info) && "Batch" %in% colnames(sample_info)) {
      message("[DEBUG] generate_quality_report: checking batch effect...")
      si_names <- standardize_sample_name(rownames(sample_info))
      batch_all <- sample_info$Batch
      batch_vec <- rep(NA_character_, nrow(pca_df))
      for (k in seq_len(nrow(pca_df))) {
        idx <- which(si_names == samples_pca[k])
        if (length(idx) == 1) batch_vec[k] <- batch_all[idx]
      }
      valid_batch <- !is.na(batch_vec)
      if (sum(valid_batch) >= 2 && length(unique(batch_vec[valid_batch])) >= 2) {
        pc1_valid <- pc1[valid_batch]
        batch_valid <- batch_vec[valid_batch]
        batch_means <- tapply(pc1_valid, batch_valid, mean)
        within_vars <- tapply(pc1_valid, batch_valid, var)
        within_var <- mean(within_vars, na.rm = TRUE)
        batch_var <- var(batch_means, na.rm = TRUE)
        if (!is.na(batch_var) && !is.na(within_var) && within_var > 0) {
          f_ratio <- batch_var / within_var
          if (f_ratio > 3) {
            key_findings <- c(key_findings, list(list(
              type = "warning", title = "⚠️ 检测到潜在批次效应",
              content = "PCA 显示样本按批次聚集，可能引入系统性偏差。\n→ 对应预处理操作：在预处理中启用批次校正（ComBat）"
            )))
          }
        }
      }
    }
    
    # 异常样本检测
    z1 <- abs((pc1 - mean(pc1)) / sd(pc1))
    z2 <- abs((pc2 - mean(pc2)) / sd(pc2))
    outlier_mask <- (z1 > 3 | z2 > 3)
    if (any(outlier_mask)) {
      outlier_samples <- samples_pca[outlier_mask]
      key_findings <- c(key_findings, list(list(
        type = "danger", title = "检测到异常样本",
        content = paste0("PCA 发现 ", length(outlier_samples), " 个可能的异常样本：", paste(outlier_samples, collapse = ", "), "，严重拉低样本一致性。\n→ 对应预处理操作：核对实验记录后，从样本信息表中移除该样本")
      )))
    }
    
    # 分组效果
    if (!is.null(pca_df) && "Group" %in% colnames(pca_df) && length(unique(pca_df$Group)) >= 2) {
      groups_pca <- pca_df$Group
      group_means <- tapply(pc1, groups_pca, mean)
      group_sd <- sd(pc1)
      if (length(group_means) >= 2 && !is.na(group_sd) && group_sd > 0) {
        group_diff <- max(group_means) - min(group_means)
        if (group_diff > group_sd) {
          key_findings <- c(key_findings, list(list(
            type = "success", title = "✅ 良好的分组效果",
            content = "不同处理组在 PCA 中明显分离，实验效应显著。\n→ 对应预处理操作：可直接进行差异分析"
          )))
        } else {
          key_findings <- c(key_findings, list(list(
            type = "warning", title = "⚠️ 分组效果不显著",
            content = "PCA 图中各组间未明显分开，可能处理效应较弱或噪声较大。\n→ 对应预处理操作：提高强度过滤阈值或移除异常样本，以增强组间差异"
          )))
        }
      }
    }
  } else {
    message("[DEBUG] generate_quality_report: PCA returned NULL, skipping deeper checks.")
  }
  
  # 推荐操作
  if (details$missing_ratio > 20) {
    recommendations <- c(recommendations, list(list(
      title = "缺失值处理", tag = "推荐", tag_type = "danger",
      items = c("设置 Max missing fraction = 0.5", "选择 KNN 或 PPCA 填充")
    )))
  } else {
    recommendations <- c(recommendations, list(list(
      title = "缺失值处理", tag = "可选", tag_type = "secondary",
      items = c("当前缺失率较低，可选用 KNN 填充（推荐）或保留缺失值。")
    )))
  }
  
  if (!is.null(sample_info) && "Batch" %in% colnames(sample_info)) {
    recommendations <- c(recommendations, list(list(
      title = "批次校正", tag = "推荐", tag_type = "warning",
      items = c("如果存在批次效应，请勾选“启用 ComBat 批次校正”。")
    )))
  }
  
  special_note <- "请根据上方的关键发现和操作建议逐步优化数据，以获得最佳分析结果。"
  if (details$missing_ratio > 40) {
    special_note <- paste0("您的数据缺失率较高（", details$missing_ratio, "%），但通过合理的过滤和填充，仍可获得可靠结果。请务必关注异常样本。")
  }
  
  message("[DEBUG] generate_quality_report: report generated with ", length(key_findings), " findings.")
  list(key_findings = key_findings, recommendations = recommendations, special_note = special_note)
}

calculate_missing_stats <- function(expr_matrix) {
  message("[DEBUG] calculate_missing_stats: dim = ", nrow(expr_matrix), "x", ncol(expr_matrix))
  if (is.null(expr_matrix) || nrow(expr_matrix) == 0) return(list())
  protein_missing <- rowMeans(is.na(expr_matrix))
  protein_missing_stats <- quantile(protein_missing, c(0, 0.25, 0.5, 0.75, 1))
  sample_missing <- colMeans(is.na(expr_matrix))
  sample_missing_stats <- quantile(sample_missing, c(0, 0.25, 0.5, 0.75, 1))
  list(protein_missing = protein_missing,
       protein_missing_stats = round(protein_missing_stats * 100, 2),
       sample_missing = sample_missing,
       sample_missing_stats = round(sample_missing_stats * 100, 2),
       total_missing_ratio = round(mean(is.na(expr_matrix)) * 100, 2))
}

calculate_sample_correlation <- function(expr_matrix) {
  message("[DEBUG] calculate_sample_correlation: dim = ", nrow(expr_matrix), "x", ncol(expr_matrix))
  if (is.null(expr_matrix) || nrow(expr_matrix) < 2 || ncol(expr_matrix) < 2) {
    message("[DEBUG] calculate_sample_correlation: insufficient data, return NULL")
    return(NULL)
  }
  
  original_na_count <- sum(is.na(expr_matrix))
  message("[DEBUG] calculate_sample_correlation: original NA count = ", original_na_count)
  
  filled <- tryCatch({
    message("[DEBUG] calculate_sample_correlation: attempting KNN imputation")
    impute_missing_values(as.data.frame(expr_matrix), method = "knn")
  }, error = function(e) {
    message("[DEBUG] calculate_sample_correlation: KNN failed (", e$message, "), using simple min fill")
    expr_matrix[is.na(expr_matrix)] <- 1e-4
    return(expr_matrix)
  })
  filled <- as.matrix(filled)
  new_na_count <- sum(is.na(filled))
  message("[DEBUG] calculate_sample_correlation: after imputation, NA count = ", new_na_count)
  
  log_expr <- log2(filled + 1)
  row_vars <- apply(log_expr, 1, var, na.rm = TRUE)
  n_keep <- min(500, nrow(log_expr))
  top_var <- order(row_vars, decreasing = TRUE)[1:n_keep]
  log_expr <- log_expr[top_var, , drop = FALSE]
  if (nrow(log_expr) < 2) {
    message("[DEBUG] calculate_sample_correlation: not enough variable rows, return NULL")
    return(NULL)
  }
  cor_matrix <- cor(log_expr, use = "complete.obs")
  na_count <- sum(is.na(cor_matrix))
  message(sprintf("[DEBUG] calculate_sample_correlation: cor matrix computed, NA count=%d", na_count))
  if (na_count > 0) {
    message("[DEBUG] calculate_sample_correlation: replacing NAs with 0")
    cor_matrix[is.na(cor_matrix)] <- 0
  }
  return(cor_matrix)
}

calculate_pca <- function(expr_matrix, sample_info = NULL) {
  message("[DEBUG] calculate_pca (data_upload) starting, dim = ", nrow(expr_matrix), "x", ncol(expr_matrix))
  if (is.null(expr_matrix) || nrow(expr_matrix) < 2 || ncol(expr_matrix) < 2) {
    message("[DEBUG] calculate_pca: insufficient data, return NULL")
    return(NULL)
  }
  
  filled <- tryCatch({
    message("[DEBUG] calculate_pca: attempting impute_missing_values (global) with KNN")
    impute_missing_values(as.data.frame(expr_matrix), method = "knn")
  }, error = function(e) {
    message("[DEBUG] calculate_pca: impute_missing_values failed (", e$message, "), using simple min fill")
    expr_matrix[is.na(expr_matrix)] <- 1e-4
    return(expr_matrix)
  })
  filled <- as.matrix(filled)
  message("[DEBUG] calculate_pca: missing values processed, remaining NAs: ", sum(is.na(filled)))
  
  log_expr <- log2(filled + 1)
  row_vars <- apply(log_expr, 1, var)
  log_expr <- log_expr[row_vars > 1e-6, ]
  row_unique <- apply(log_expr, 1, function(x) length(unique(x)))
  log_expr <- log_expr[row_unique > 1, ]
  if (nrow(log_expr) < 2) {
    message("[DEBUG] calculate_pca: not enough variable rows after filtering, return NULL")
    return(NULL)
  }
  
  tryCatch({
    pca_result <- prcomp(t(log_expr), scale. = TRUE)
    var_explained <- round(pca_result$sdev^2 / sum(pca_result$sdev^2) * 100, 1)
    pca_df <- as.data.frame(pca_result$x[, 1:2])
    pca_df$Sample <- rownames(pca_df)
    
    if (!is.null(sample_info) && "Group" %in% colnames(sample_info)) {
      sample_info_short <- sample_info
      rownames(sample_info_short) <- standardize_sample_name(rownames(sample_info_short))
      pca_sample_std <- standardize_sample_name(pca_df$Sample)
      common_idx <- match(pca_sample_std, rownames(sample_info_short))
      if (any(!is.na(common_idx))) {
        pca_df$Group <- sample_info_short$Group[common_idx]
      } else {
        pca_df$Group <- "All"
      }
    } else {
      pca_df$Group <- "All"
    }
    
    message("[DEBUG] calculate_pca: successful, PC1=", var_explained[1], "%, PC2=", var_explained[2], "%")
    list(pca_df = pca_df, var_explained = var_explained,
         pc1_var = var_explained[1], pc2_var = var_explained[2])
  }, error = function(e) {
    message("[DEBUG] calculate_pca: prcomp error - ", e$message)
    NULL
  })
}

# ================== 数据质量反应式（供 data_quality_plots.R 使用） ==================
dq_expr_matrix <- reactive({
  message("[DEBUG] dq_expr_matrix triggered")
  req(rv$clean_data, rv$lfq_cols)
  df <- rv$clean_data[, rv$lfq_cols, drop = FALSE]
  df <- suppressWarnings(as.data.frame(lapply(df, as.numeric)))
  df[df == 0] <- NA
  rownames(df) <- rv$clean_data$`Master protein IDs`
  colnames(df) <- rv$sample_names
  df
}) %>% bindCache(rv$clean_data, rv$lfq_cols)

dq_quality_score <- reactive({
  message("[DEBUG] dq_quality_score triggered")
  req(dq_expr_matrix())
  calculate_data_quality_score(dq_expr_matrix())
}) %>% bindCache(dq_expr_matrix())

dq_missing_stats <- reactive({
  message("[DEBUG] dq_missing_stats triggered")
  req(dq_expr_matrix())
  calculate_missing_stats(dq_expr_matrix())
}) %>% bindCache(dq_expr_matrix())

message("[DEBUG] data_upload_quality.R loaded successfully.")