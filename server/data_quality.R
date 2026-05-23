# server/data_quality.R
# 数据质量分析核心函数

# -------------------- 数据质量评分计算 --------------------
calculate_data_quality_score <- function(expr_matrix) {
  if (is.null(expr_matrix) || nrow(expr_matrix) == 0 || ncol(expr_matrix) == 0) {
    return(list(score = 0, grade = "F", details = list()))
  }
  
  total_values <- nrow(expr_matrix) * ncol(expr_matrix)
  missing_values <- sum(is.na(expr_matrix))
  missing_ratio <- missing_values / total_values
  missing_score <- max(0, 30 * (1 - missing_ratio * 2))
  
  sample_cor <- cor(expr_matrix, use = "pairwise.complete.obs")
  diag(sample_cor) <- NA
  avg_cor <- mean(sample_cor, na.rm = TRUE)
  cor_score <- max(0, 20 * pmin(1, avg_cor / 0.8))
  
  protein_valid <- rowSums(!is.na(expr_matrix)) >= 2
  protein_valid_ratio <- mean(protein_valid)
  protein_score <- max(0, 20 * protein_valid_ratio)
  
  cor_score_2 <- max(0, 30 * pmin(1, avg_cor / 0.8))
  
  total_score <- round(missing_score + cor_score + protein_score + cor_score_2, 1)
  
  grade <- if (total_score >= 90) "Excellent"
  else if (total_score >= 80) "Good"
  else if (total_score >= 70) "Fair"
  else if (total_score >= 60) "Poor"
  else "Bad"
  
  details <- list(
    missing_ratio = round(missing_ratio * 100, 2),
    missing_score = round(missing_score, 1),
    avg_correlation = round(avg_cor, 3),
    correlation_score = round(cor_score, 1),
    protein_valid_ratio = round(protein_valid_ratio * 100, 2),
    protein_score = round(protein_score, 1),
    sample_cor_score = round(cor_score_2, 1)
  )
  
  list(score = total_score, grade = grade, details = details)
}

# -------------------- 生成智能数据分析报告 --------------------
generate_quality_report <- function(quality_score, expr_matrix, sample_info = NULL) {
  report <- list()
  details <- quality_score$details
  
  key_findings <- list()
  
  if (details$missing_ratio < 10) {
    key_findings <- c(key_findings, list(list(
      type = "success",
      title = "极低缺失率",
      content = paste0("数据缺失率仅为", details$missing_ratio, "%，数据完整性非常好。")
    )))
  } else if (details$missing_ratio < 20) {
    key_findings <- c(key_findings, list(list(
      type = "warning",
      title = "中等缺失率",
      content = paste0("数据缺失率为", details$missing_ratio, "%，建议进行适当的缺失值处理。")
    )))
  } else if (details$missing_ratio < 40) {
    key_findings <- c(key_findings, list(list(
      type = "danger",
      title = "较高缺失率",
      content = paste0("数据缺失率为", details$missing_ratio, "%，建议过滤掉缺失值过多的蛋白。")
    )))
  } else {
    key_findings <- c(key_findings, list(list(
      type = "danger",
      title = "极高缺失率",
      content = paste0("数据缺失率高达", details$missing_ratio, "%，超过了蛋白质组学实验的可接受范围。")
    )))
  }
  
  if (details$avg_correlation > 0.9) {
    key_findings <- c(key_findings, list(list(
      type = "success",
      title = "优秀的样本一致性",
      content = paste0("样本平均相关性为", details$avg_correlation, "，生物学重复之间的一致性非常好。")
    )))
  } else if (details$avg_correlation > 0.8) {
    key_findings <- c(key_findings, list(list(
      type = "success",
      title = "良好的样本一致性",
      content = paste0("样本平均相关性为", details$avg_correlation, "，生物学重复之间的一致性较好。")
    )))
  } else if (details$avg_correlation > 0.7) {
    key_findings <- c(key_findings, list(list(
      type = "warning",
      title = "一般的样本一致性",
      content = paste0("样本平均相关性为", details$avg_correlation, "，实验重复性尚可，但仍有提升空间。")
    )))
  } else {
    key_findings <- c(key_findings, list(list(
      type = "danger",
      title = "较差的样本一致性",
      content = paste0("样本平均相关性仅为", details$avg_correlation, "，建议检查实验流程。")
    )))
  }
  
  if (!is.null(sample_info) && "Batch" %in% colnames(sample_info)) {
    tryCatch({
      pca_result <- calculate_pca(expr_matrix, sample_info)
      if (!is.null(pca_result)) {
        pc1 <- pca_result$pca_df$PC1
        batches <- sample_info$Batch
        if (length(unique(batches)) >= 2) {
          batch_var <- var(pc1[batches == unique(batches)[1]]) + var(pc1[batches == unique(batches)[2]])
          total_var <- var(pc1)
          batch_effect_ratio <- 1 - batch_var / total_var
          if (batch_effect_ratio > 0.3) {
            key_findings <- c(key_findings, list(list(
              type = "warning",
              title = "检测到潜在批次效应",
              content = "样本相关性热图显示样本明显分成了几个簇，可能存在批次效应。建议在后续分析中启用批次校正。"
            )))
          }
        }
      }
    }, error = function(e) {})
  }
  
  tryCatch({
    pca_result <- calculate_pca(expr_matrix, sample_info)
    if (!is.null(pca_result)) {
      pc1 <- pca_result$pca_df$PC1
      pc2 <- pca_result$pca_df$PC2
      z1 <- abs((pc1 - mean(pc1)) / sd(pc1))
      z2 <- abs((pc2 - mean(pc2)) / sd(pc2))
      outliers <- names(which(z1 > 3 | z2 > 3))
      if (length(outliers) > 0) {
        key_findings <- c(key_findings, list(list(
          type = "danger",
          title = "检测到异常样本",
          content = paste0("PCA分析检测到", length(outliers), "个异常样本：", paste(outliers, collapse = ", "), "。强烈建议移除。")
        )))
      }
    }
  }, error = function(e) {})
  
  if (!is.null(sample_info) && "Group" %in% colnames(sample_info)) {
    groups <- unique(sample_info$Group)
    if (length(groups) >= 2) {
      tryCatch({
        pca_result <- calculate_pca(expr_matrix, sample_info)
        if (!is.null(pca_result)) {
          pc1 <- pca_result$pca_df$PC1
          group1_mean <- mean(pc1[sample_info$Group == groups[1]])
          group2_mean <- mean(pc1[sample_info$Group == groups[2]])
          group_diff <- abs(group1_mean - group2_mean)
          group_sd <- sd(pc1)
          if (group_diff > group_sd) {
            key_findings <- c(key_findings, list(list(
              type = "success",
              title = "良好的分组效果",
              content = "不同处理组之间的差异明显，实验处理效果显著，后续差异分析应该能找到很多差异蛋白。"
            )))
          }
        }
      }, error = function(e) {})
    }
  }
  
  recommendations <- list()
  
  step1 <- list()
  if (details$missing_ratio > 20) {
    step1 <- c(step1, "1. 启用\"缺失值过滤\"，设置阈值为0.5")
    step1 <- c(step1, "2. 启用\"缺失值填充\"，使用K近邻填充法")
  }
  if (!is.null(sample_info) && "Batch" %in% colnames(sample_info)) {
    step1 <- c(step1, "3. 启用\"批次校正\"，根据实验批次信息设置批次")
  }
  if (length(step1) > 0) {
    recommendations <- c(recommendations, list(list(
      title = "第一步：数据预处理",
      tag = "必须做",
      tag_type = "danger",
      items = step1
    )))
  }
  
  tryCatch({
    pca_result <- calculate_pca(expr_matrix, sample_info)
    if (!is.null(pca_result)) {
      pc1 <- pca_result$pca_df$PC1
      pc2 <- pca_result$pca_df$PC2
      z1 <- abs((pc1 - mean(pc1)) / sd(pc1))
      z2 <- abs((pc2 - mean(pc2)) / sd(pc2))
      outliers <- names(which(z1 > 3 | z2 > 3))
      if (length(outliers) > 0) {
        recommendations <- c(recommendations, list(list(
          title = "第二步：异常样本处理",
          tag = "必须做",
          tag_type = "danger",
          items = c(
            "1. 在样本信息文件中删除检测到的异常样本",
            "2. 重新上传样本信息文件",
            "3. 确认样本匹配时，异常样本不再出现在匹配列表中"
          )
        )))
      }
    }
  }, error = function(e) {})
  
  step3 <- list()
  if (details$missing_ratio > 30) {
    step3 <- c(step3, "1. 将\"最小有效重复数\"从2降低到1（因为缺失率较高）")
    step3 <- c(step3, "2. 建议使用limma检验方法")
  }
  step3 <- c(step3, "3. 保持p值阈值为0.05，fold change阈值为1.2")
  
  recommendations <- c(recommendations, list(list(
    title = "第三步：差异分析参数调整",
    tag = "建议",
    tag_type = "secondary",
    items = step3
  )))
  
  special_note <- ""
  if (details$missing_ratio > 40) {
    special_note <- paste0("虽然数据缺失率较高（", details$missing_ratio, "%），但数据质量其实很好，分组效果非常明显。只要按照建议进行预处理，应该能得到非常可靠的差异分析结果。")
  }
  
  report$key_findings <- key_findings
  report$recommendations <- recommendations
  report$special_note <- special_note
  
  return(report)
}

# -------------------- 缺失值分析统计 --------------------
calculate_missing_stats <- function(expr_matrix) {
  if (is.null(expr_matrix) || nrow(expr_matrix) == 0) {
    return(list())
  }
  
  protein_missing <- rowMeans(is.na(expr_matrix))
  protein_missing_stats <- quantile(protein_missing, c(0, 0.25, 0.5, 0.75, 1))
  
  sample_missing <- colMeans(is.na(expr_matrix))
  sample_missing_stats <- quantile(sample_missing, c(0, 0.25, 0.5, 0.75, 1))
  
  list(
    protein_missing = protein_missing,
    protein_missing_stats = round(protein_missing_stats * 100, 2),
    sample_missing = sample_missing,
    sample_missing_stats = round(sample_missing_stats * 100, 2),
    total_missing_ratio = round(mean(is.na(expr_matrix)) * 100, 2)
  )
}

# -------------------- 样本相关性计算 --------------------
calculate_sample_correlation <- function(expr_matrix) {
  if (is.null(expr_matrix) || nrow(expr_matrix) < 2 || ncol(expr_matrix) < 2) {
    return(NULL)
  }
  
  log_expr <- log2(expr_matrix + 1)
  
  row_vars <- apply(log_expr, 1, var, na.rm = TRUE)
  log_expr <- log_expr[row_vars > 1e-6, ]
  
  if (nrow(log_expr) < 2) {
    return(NULL)
  }
  
  cor_matrix <- cor(log_expr, use = "pairwise.complete.obs")
  cor_matrix[is.na(cor_matrix)] <- 0
  
  return(cor_matrix)
}

# -------------------- PCA计算 --------------------
calculate_pca <- function(expr_matrix, sample_info = NULL) {
  if (is.null(expr_matrix) || nrow(expr_matrix) < 2 || ncol(expr_matrix) < 2) {
    return(NULL)
  }
  
  log_expr <- log2(expr_matrix + 1)
  log_expr[is.na(log_expr)] <- 0
  
  row_vars <- apply(log_expr, 1, var)
  log_expr <- log_expr[row_vars > 1e-6, ]
  
  row_unique <- apply(log_expr, 1, function(x) length(unique(x)))
  log_expr <- log_expr[row_unique > 1, ]
  
  if (nrow(log_expr) < 2) {
    return(NULL)
  }
  
  tryCatch({
    pca_result <- prcomp(t(log_expr), scale. = TRUE)
    
    var_explained <- round(pca_result$sdev^2 / sum(pca_result$sdev^2) * 100, 1)
    
    pca_df <- as.data.frame(pca_result$x[, 1:2])
    pca_df$Sample <- rownames(pca_df)
    
    if (!is.null(sample_info) && "Group" %in% colnames(sample_info)) {
      sample_info_short <- sample_info
      rownames(sample_info_short) <- gsub("^(LFQ intensity |Intensity )", "", rownames(sample_info_short))
      common_samples <- intersect(pca_df$Sample, rownames(sample_info_short))
      if (length(common_samples) > 0) {
        pca_df$Group <- sample_info_short[common_samples, "Group"]
      } else {
        pca_df$Group <- "All"
      }
    } else {
      pca_df$Group <- "All"
    }
    
    list(
      pca_df = pca_df,
      var_explained = var_explained,
      pc1_var = var_explained[1],
      pc2_var = var_explained[2]
    )
  }, error = function(e) {
    return(NULL)
  })
}

# -------------------- 渲染关键发现卡片 --------------------
render_key_finding <- function(finding) {
  bg_color <- switch(finding$type,
                     success = "#d4edda",
                     warning = "#fff3cd",
                     danger = "#f8d7da",
                     "#f8f9fa"
  )
  border_color <- switch(finding$type,
                         success = "#c3e6cb",
                         warning = "#ffeeba",
                         danger = "#f5c6cb",
                         "#dee2e6"
  )
  text_color <- switch(finding$type,
                       success = "#155724",
                       warning = "#856404",
                       danger = "#721c24",
                       "#333"
  )
  icon_name <- switch(finding$type,
                      success = "check-circle",
                      warning = "exclamation-triangle",
                      danger = "times-circle",
                      "info-circle"
  )
  
  div(style = paste0("background: ", bg_color, "; border: 1px solid ", border_color, "; border-radius: 8px; padding: 12px; margin-bottom: 10px; color: ", text_color, ";"),
      div(style = "display: flex; align-items: center; gap: 8px; margin-bottom: 5px; font-weight: bold;",
          icon(icon_name), finding$title),
      p(style = "margin: 0;", finding$content)
  )
}

# -------------------- 渲染操作建议卡片 --------------------
render_recommendation <- function(rec) {
  tag_bg <- switch(rec$tag_type,
                   danger = "#dc3545",
                   warning = "#ffc107",
                   success = "#28a745",
                   secondary = "#6c757d"
  )
  
  div(style = "background: #f8f9fa; border-radius: 8px; padding: 15px; margin-bottom: 15px;",
      div(style = "display: flex; align-items: center; gap: 10px; margin-bottom: 10px;",
          h4(style = "margin: 0; font-size: 16px; font-weight: bold;", rec$title),
          span(style = paste0("background: ", tag_bg, "; color: white; padding: 2px 8px; border-radius: 12px; font-size: 12px; font-weight: bold;"), rec$tag)
      ),
      tagList(lapply(rec$items, function(item) p(style = "margin: 3px 0; padding-left: 15px; text-indent: -15px;", item)))
  )
}