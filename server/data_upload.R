# server/data_upload.R

# ================== 全局辅助函数 ==================
extract_sample_names <- function(cols) {
  short <- sub("^(LFQ intensity |Intensity )", "", cols, ignore.case = TRUE)
  still_full <- which(short == cols)
  if (length(still_full) > 0) {
    short[still_full] <- sub("^(LFQ[._]?intensity[._]?|Intensity[._]?)", "", cols[still_full], ignore.case = TRUE)
  }
  short <- sub("^[^[:alnum:]]+", "", short)
  short <- sub("[^[:alnum:]]+$", "", short)
  short <- gsub("-", ".", short)
  short
}

get_raw_prefix <- function(type = input$intensity_type) {
  if (type == "LFQ") "LFQ intensity " else "Intensity "
}

get_norm_prefix <- function(type = input$intensity_type) {
  if (type == "LFQ") "Norm_LFQ intensity " else "Norm_Intensity "
}

extract_group_name <- function(x) {
  x <- sub("^Norm_([A-Za-z]+ intensity) ", "", x)
  sapply(strsplit(x, "-"), `[`, 1)
}

get_group_colors <- function(groups) {
  if (length(groups) == 0) return(character(0))
  pal <- c("#FFB6C1","#90EE90","#87CEEB","#DDA0DD","#FFD700",
           "#FFA07A","#98FB98","#B0C4DE","#FFB347","#C9A0DC")
  if (length(groups) > length(pal)) {
    pal <- colorRampPalette(pal)(length(groups))
  }
  setNames(pal[1:length(groups)], groups)
}

# ================== 数据质量分析辅助函数 ==================
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

# ================== 优化后的关键发现报告生成函数 ==================
generate_quality_report <- function(quality_score, expr_matrix, sample_info = NULL) {
  details <- quality_score$details
  key_findings <- list()
  recommendations <- list()
  special_note <- ""
  
  # ---- 缺失率发现 ----
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
  
  # ---- 样本一致性发现 ----
  if (details$avg_correlation > 0.9) {
    key_findings <- c(key_findings, list(list(
      type = "success",
      title = "✅ 优秀的样本一致性",
      content = paste0("样本间平均相关性为 ", round(details$avg_correlation, 3), "，实验重复性极佳。\n→ 对应预处理操作：无需额外处理")
    )))
  } else if (details$avg_correlation > 0.8) {
    key_findings <- c(key_findings, list(list(
      type = "success",
      title = "✅ 良好的样本一致性",
      content = paste0("样本间平均相关性为 ", round(details$avg_correlation, 3), "，实验重复性较好。\n→ 对应预处理操作：可考虑执行批次校正以进一步提升一致性")
    )))
  } else if (details$avg_correlation > 0.7) {
    key_findings <- c(key_findings, list(list(
      type = "warning",
      title = "⚠️ 一般的样本一致性",
      content = paste0("样本间平均相关性为 ", round(details$avg_correlation, 3), "，可能存在技术波动。\n→ 对应预处理操作：检查并移除异常样本，执行缺失值填充")
    )))
  } else {
    key_findings <- c(key_findings, list(list(
      type = "danger",
      title = "较差的样本一致性",
      content = paste0("样本间平均相关性仅为 ", round(details$avg_correlation, 3), "，样本重复性较差，会降低后续分析的统计效力。\n→ 对应预处理操作：异常样本移除 + 缺失值填充")
    )))
  }
  
  # ---- 蛋白有效检出率 ----
  if (details$protein_valid_ratio > 80) {
    key_findings <- c(key_findings, list(list(
      type = "success",
      title = "蛋白有效检出率高",
      content = paste0("超过 ", details$protein_valid_ratio, "% 的蛋白在至少 2 个样本中被检出，蛋白整体质量优秀。\n→ 对应预处理操作：最小强度过滤阈值可设为 100000，无需大幅调整")
    )))
  } else if (details$protein_valid_ratio > 60) {
    key_findings <- c(key_findings, list(list(
      type = "warning",
      title = "⚠️ 蛋白有效检出率一般",
      content = paste0("约 ", details$protein_valid_ratio, "% 的蛋白在至少 2 个样本中被检出，部分蛋白检出率偏低。\n→ 对应预处理操作：适当降低强度过滤阈值，保留更多蛋白")
    )))
  } else {
    key_findings <- c(key_findings, list(list(
      type = "danger",
      title = "ⓧ 蛋白有效检出率低",
      content = paste0("仅有 ", details$protein_valid_ratio, "% 的蛋白在至少 2 个样本中被检出，数据可靠性较低。\n→ 对应预处理操作：建议放宽强度过滤阈值，并核对实验流程")
    )))
  }
  
  # ---- PCA 分析（批次效应、异常样本、分组效果）----
  pca_result <- calculate_pca(expr_matrix, sample_info)
  pca_df <- if (!is.null(pca_result)) pca_result$pca_df else NULL
  
  batch_effect_detected <- FALSE
  outlier_samples <- character(0)
  group_separation_good <- FALSE
  has_groups <- !is.null(pca_df) && "Group" %in% colnames(pca_df) && length(unique(pca_df$Group)) >= 2
  
  if (!is.null(pca_df) && nrow(pca_df) >= 2) {
    pc1 <- pca_df$PC1
    pc2 <- pca_df$PC2
    samples_pca <- pca_df$Sample
    
    if (!is.null(sample_info) && "Batch" %in% colnames(sample_info)) {
      si_names <- gsub("-", ".", rownames(sample_info))
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
          if (f_ratio > 3) batch_effect_detected <- TRUE
        }
      }
    }
    
    z1 <- abs((pc1 - mean(pc1)) / sd(pc1))
    z2 <- abs((pc2 - mean(pc2)) / sd(pc2))
    outlier_mask <- (z1 > 3 | z2 > 3)
    if (any(outlier_mask)) {
      outlier_samples <- samples_pca[outlier_mask]
    }
    
    if (has_groups) {
      groups_pca <- pca_df$Group
      group_means <- tapply(pc1, groups_pca, mean)
      group_sd <- sd(pc1)
      if (length(group_means) >= 2 && !is.na(group_sd) && group_sd > 0) {
        group_diff <- max(group_means) - min(group_means)
        if (group_diff > group_sd) group_separation_good <- TRUE
      }
    }
  }
  
  if (batch_effect_detected) {
    key_findings <- c(key_findings, list(list(
      type = "warning",
      title = "⚠️ 检测到潜在批次效应",
      content = "PCA 显示样本按批次聚集，可能引入系统性偏差。\n→ 对应预处理操作：在预处理中启用批次校正（ComBat）"
    )))
  }
  if (length(outlier_samples) > 0) {
    key_findings <- c(key_findings, list(list(
      type = "danger",
      title = "检测到异常样本",
      content = paste0("PCA 发现 ", length(outlier_samples), " 个可能的异常样本：", paste(outlier_samples, collapse = ", "), "，严重拉低样本一致性。\n→ 对应预处理操作：核对实验记录后，从样本信息表中移除该样本")
    )))
  }
  if (group_separation_good && !batch_effect_detected) {
    key_findings <- c(key_findings, list(list(
      type = "success",
      title = "✅ 良好的分组效果",
      content = "不同处理组在 PCA 中明显分离，实验效应显著。\n→ 对应预处理操作：可直接进行差异分析"
    )))
  } else if (!group_separation_good && has_groups && !batch_effect_detected) {
    key_findings <- c(key_findings, list(list(
      type = "warning",
      title = "⚠️ 分组效果不显著",
      content = "PCA 图中各组间未明显分开，可能处理效应较弱或噪声较大。\n→ 对应预处理操作：提高强度过滤阈值或移除异常样本，以增强组间差异"
    )))
  }
  
  # ---- 操作建议（保持不变） ----
  if (details$missing_ratio > 20) {
    recommendations <- c(recommendations, list(list(
      title = "缺失值处理",
      tag = "推荐",
      tag_type = "danger",
      items = c(
        paste0("当前缺失率 ", details$missing_ratio, "% > 20%，建议："),
        "1. 设置 Max missing fraction = 0.5（过滤缺失过半的蛋白）",
        "2. 选择 KNN 或 PPCA 填充",
        "3. 如缺失率过高，考虑剔除相应样本"
      )
    )))
  } else {
    recommendations <- c(recommendations, list(list(
      title = "缺失值处理",
      tag = "可选",
      tag_type = "secondary",
      items = c("当前缺失率较低，可选用 KNN 填充（推荐）或保留缺失值。")
    )))
  }
  
  if (batch_effect_detected) {
    recommendations <- c(recommendations, list(list(
      title = "批次校正",
      tag = "强烈推荐",
      tag_type = "danger",
      items = c("检测到批次效应，请务必在预处理中勾选 'Perform batch correction' 并选择批次列。")
    )))
  } else {
    recommendations <- c(recommendations, list(list(
      title = "批次校正",
      tag = "可选",
      tag_type = "secondary",
      items = c("未检测到显著批次效应，但如有实验批次设计，仍可尝试校正。")
    )))
  }
  
  if (length(outlier_samples) > 0) {
    recommendations <- c(recommendations, list(list(
      title = "异常样本处理",
      tag = "必须执行",
      tag_type = "danger",
      items = c(
        paste0("检测到以下异常样本：", paste(outlier_samples, collapse = ", ")),
        "1. 核对实验记录",
        "2. 从样本信息表中删除这些样本，重新上传后再分析"
      )
    )))
  }
  
  if (details$missing_ratio > 30) {
    recommendations <- c(recommendations, list(list(
      title = "统计方法",
      tag = "建议",
      tag_type = "warning",
      items = c("缺失率较高，建议使用 limma 检验，并适当降低最低重复数。")
    )))
  } else {
    recommendations <- c(recommendations, list(list(
      title = "统计方法",
      tag = "默认",
      tag_type = "secondary",
      items = c("可使用 t-test 或 Wilcoxon，保留默认参数即可。")
    )))
  }
  
  if (details$missing_ratio > 40) {
    special_note <- paste0("您的数据缺失率较高（", details$missing_ratio, "%），但通过合理的过滤和填充，仍可获得可靠结果。请务必关注异常样本。")
  } else if (details$missing_ratio < 10 && details$avg_correlation > 0.9) {
    special_note <- "您的数据质量非常高，可以直接进行分析，预期结果可靠。"
  } else if (length(outlier_samples) > 0 && batch_effect_detected) {
    special_note <- "数据中存在潜在异常样本和批次效应，请优先处理这两项，否则后续分析可能受到较大干扰。"
  } else if (group_separation_good && !batch_effect_detected) {
    special_note <- "实验处理效果显著，PCA 组间分离清晰，有望筛选出大量差异蛋白。"
  } else {
    special_note <- "请根据上方的关键发现和操作建议逐步优化数据，以获得最佳分析结果。"
  }
  
  list(key_findings = key_findings, recommendations = recommendations, special_note = special_note)
}

calculate_missing_stats <- function(expr_matrix) {
  if (is.null(expr_matrix) || nrow(expr_matrix) == 0) return(list())
  protein_missing <- rowMeans(is.na(expr_matrix))
  protein_missing_stats <- quantile(protein_missing, c(0, 0.25, 0.5, 0.75, 1))
  sample_missing <- colMeans(is.na(expr_matrix))
  sample_missing_stats <- quantile(sample_missing, c(0, 0.25, 0.5, 0.75, 1))
  list(protein_missing = protein_missing, protein_missing_stats = round(protein_missing_stats * 100, 2),
       sample_missing = sample_missing, sample_missing_stats = round(sample_missing_stats * 100, 2),
       total_missing_ratio = round(mean(is.na(expr_matrix)) * 100, 2))
}

calculate_sample_correlation <- function(expr_matrix) {
  if (is.null(expr_matrix) || nrow(expr_matrix) < 2 || ncol(expr_matrix) < 2) return(NULL)
  log_expr <- log2(expr_matrix + 1)
  row_vars <- apply(log_expr, 1, var, na.rm = TRUE)
  n_keep <- min(500, nrow(log_expr))
  top_var <- order(row_vars, decreasing = TRUE)[1:n_keep]
  log_expr <- log_expr[top_var, , drop = FALSE]
  if (nrow(log_expr) < 2) return(NULL)
  cor_matrix <- cor(log_expr, use = "pairwise.complete.obs")
  cor_matrix[is.na(cor_matrix)] <- 0
  return(cor_matrix)
}

calculate_pca <- function(expr_matrix, sample_info = NULL) {
  if (is.null(expr_matrix) || nrow(expr_matrix) < 2 || ncol(expr_matrix) < 2) return(NULL)
  log_expr <- log2(expr_matrix + 1)
  log_expr[is.na(log_expr)] <- 0
  row_vars <- apply(log_expr, 1, var)
  n_keep <- min(500, nrow(log_expr))
  top_var <- order(row_vars, decreasing = TRUE)[1:n_keep]
  log_expr <- log_expr[top_var, , drop = FALSE]
  row_unique <- apply(log_expr, 1, function(x) length(unique(x)))
  log_expr <- log_expr[row_unique > 1, , drop = FALSE]
  if (nrow(log_expr) < 2) return(NULL)
  tryCatch({
    pca_result <- prcomp(t(log_expr), scale. = TRUE)
    var_explained <- round(pca_result$sdev^2 / sum(pca_result$sdev^2) * 100, 1)
    pca_df <- as.data.frame(pca_result$x[, 1:2])
    pca_df$Sample <- rownames(pca_df)
    if (!is.null(sample_info) && "Group" %in% colnames(sample_info)) {
      sample_info_short <- sample_info
      rownames(sample_info_short) <- gsub("^(LFQ intensity |Intensity )", "", rownames(sample_info_short))
      common_samples <- intersect(pca_df$Sample, rownames(sample_info_short))
      if (length(common_samples) > 0) pca_df$Group <- sample_info_short[common_samples, "Group"]
      else pca_df$Group <- "All"
    } else {
      pca_df$Group <- "All"
    }
    list(pca_df = pca_df, var_explained = var_explained, pc1_var = var_explained[1], pc2_var = var_explained[2])
  }, error = function(e) NULL)
}

render_key_finding <- function(finding) {
  bg_color <- switch(finding$type, success = "#d4edda", warning = "#fff3cd", danger = "#f8d7da", "#f8f9fa")
  border_color <- switch(finding$type, success = "#c3e6cb", warning = "#ffeeba", danger = "#f5c6cb", "#dee2e6")
  text_color <- switch(finding$type, success = "#155724", warning = "#856404", danger = "#721c24", "#333")
  icon_name <- switch(finding$type, success = "check-circle", warning = "exclamation-triangle", danger = "times-circle", "info-circle")
  div(style = paste0("background: ", bg_color, "; border: 1px solid ", border_color, "; border-radius: 8px; padding: 12px; margin-bottom: 10px; color: ", text_color, ";"),
      div(style = "display: flex; align-items: center; gap: 8px; margin-bottom: 5px; font-weight: bold;", icon(icon_name), finding$title),
      p(style = "margin: 0; white-space: pre-line;", finding$content))
}

render_recommendation <- function(rec) {
  tag_bg <- switch(rec$tag_type, danger = "#dc3545", warning = "#ffc107", success = "#28a745", secondary = "#6c757d")
  div(style = "background: #f8f9fa; border-radius: 8px; padding: 15px; margin-bottom: 15px;",
      div(style = "display: flex; align-items: center; gap: 10px; margin-bottom: 10px;",
          h4(style = "margin: 0; font-size: 16px; font-weight: bold;", rec$title),
          span(style = paste0("background: ", tag_bg, "; color: white; padding: 2px 8px; border-radius: 12px; font-size: 12px; font-weight: bold;"), rec$tag)),
      tagList(lapply(rec$items, function(item) p(style = "margin: 3px 0; padding-left: 15px; text-indent: -15px;", item))))
}

# ================== 数据上传与预处理 ==================
output$download_sample_template <- downloadHandler(
  filename = function() { "sample_info_template.xlsx" },
  content = function(file) {
    prefix <- get_raw_prefix()
    template <- data.frame(
      SampleName = paste0(prefix, c("L2.1.1", "L2.1.2", "L2.1.3", "L2.2.1", "L2.2.2", "L2.2.3")),
      Group = c("Control", "Control", "Control", "Treatment", "Treatment", "Treatment"),
      Batch = c("Batch1", "Batch1", "Batch2", "Batch1", "Batch2", "Batch2"),
      Note = c("", "", "", "", "", "")
    )
    writexl::write_xlsx(template, file)
  }
)

read_sample_info <- function(file_path) {
  ext <- tools::file_ext(file_path)
  if (ext %in% c("csv", "txt")) {
    df <- read.csv(file_path, header = TRUE, row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
  } else if (ext %in% c("xlsx", "xls")) {
    df <- as.data.frame(readxl::read_excel(file_path, col_names = TRUE))
  } else {
    stop("Unsupported file format for sample info.")
  }
  rownames(df) <- as.character(df[[1]])
  df <- df[, -1, drop = FALSE]
  return(df)
}

cached_sample_info <- reactiveValues(LFQ = NULL, Intensity = NULL)

observeEvent(input$intensity_type, {
  if (!is.null(rv$sample_info)) {
    if (input$intensity_type == "LFQ") cached_sample_info$Intensity <- rv$sample_info
    else cached_sample_info$LFQ <- rv$sample_info
  }
  target_cache <- if (input$intensity_type == "LFQ") cached_sample_info$LFQ else cached_sample_info$Intensity
  if (!is.null(target_cache)) rv$sample_info <- target_cache
  else { rv$sample_info <- NULL; reset("sample_info_file") }
  if (!is.null(rv$raw_data)) {
    new_cols <- grep(paste0("^", get_raw_prefix()), colnames(rv$raw_data), value = TRUE)
    if (length(new_cols) > 0) {
      rv$lfq_cols <- new_cols
      rv$sample_names <- extract_sample_names(new_cols)
      updateSelectInput(session, "baseline_sample", choices = c("Auto", rv$sample_names), selected = "Auto")
    } else {
      showNotification("No matching intensity columns found for the selected type.", type = "error", duration = 5)
    }
  }
}, ignoreInit = TRUE)

observeEvent(input$expression_file, {
  req(input$expression_file)
  rv$raw_data <- NULL; rv$clean_data <- NULL; rv$lfq_cols <- NULL
  rv$sample_names <- NULL; rv$groups <- list(); rv$comparisons <- list()
  rv$analysis_results <- NULL
  rv$comp_id_counter <- 0; rv$group_id_counter <- 0
  rv$group_id_map <- list(); rv$current_profile_protein <- NULL
  rv$batch_vector <- NULL; rv$sample_info <- NULL
  manual_sort_active(FALSE)
  updateSelectInput(session, "comp_treat", choices = character(0))
  updateSelectInput(session, "comp_ctrl", choices = character(0))
  updateSelectInput(session, "selected_comparison", choices = character(0))
  updateSelectInput(session, "baseline_sample", choices = c("Auto"), selected = "Auto")
  updateSelectInput(session, "batch_ref_group", choices = character(0))
  updateSelectizeInput(session, "venn_comparisons_select", choices = character(0), selected = character(0))
  updateCheckboxGroupInput(session, "venn_comparisons_checkbox", choices = character(0), selected = character(0))
  tryCatch({
    file_path <- input$expression_file$datapath
    data <- fread(file_path, sep = "\t", stringsAsFactors = FALSE, data.table = FALSE, check.names = FALSE, colClasses = "character")
    for (cn in names(data)) {
      if (grepl("^(LFQ intensity |Intensity )", cn)) data[[cn]] <- as.numeric(data[[cn]])
    }
    lfq_cols <- grep(paste0("^", get_raw_prefix()), colnames(data), value = TRUE)
    if (length(lfq_cols) == 0) {
      showNotification("No matching intensity columns found. Please check intensity type.", type = "error", duration = 5)
      return()
    }
    sample_names <- extract_sample_names(lfq_cols)
    clean_data <- data
    if ("Reverse" %in% colnames(clean_data)) clean_data <- filter(clean_data, is.na(Reverse) | Reverse != "+")
    if ("Potential contaminant" %in% colnames(clean_data)) clean_data <- filter(clean_data, is.na(`Potential contaminant`) | `Potential contaminant` != "+")
    if ("Protein IDs" %in% colnames(clean_data)) clean_data <- filter(clean_data, !grepl("^CON_", `Protein IDs`))
    if ("Protein IDs" %in% colnames(clean_data)) clean_data <- mutate(clean_data, `Master protein IDs` = sub(";.*", "", `Protein IDs`), .after = `Majority protein IDs`)
    rv$raw_data <- data
    if ("Protein IDs" %in% colnames(rv$raw_data)) rv$raw_data$`Master protein IDs` <- sub(";.*", "", rv$raw_data$`Protein IDs`)
    rv$clean_data <- clean_data
    rv$lfq_cols <- lfq_cols
    rv$sample_names <- sample_names
    updateSelectInput(session, "baseline_sample", choices = c("Auto", sample_names), selected = "Auto")
    showNotification("Expression matrix uploaded successfully!", type = "message", duration = 3)
  }, error = function(e) {
    showNotification(paste("Error reading expression file:", e$message), type = "error", duration = 5)
  })
})

observeEvent(input$sample_info_file, {
  req(input$sample_info_file)
  tryCatch({
    df <- read_sample_info(input$sample_info_file$datapath)
    rv$sample_info <- df
    if (input$intensity_type == "LFQ") cached_sample_info$LFQ <- df
    else cached_sample_info$Intensity <- df
    showNotification("Sample info uploaded successfully!", type = "message", duration = 3)
  }, error = function(e) {
    showNotification(paste("Error reading sample info:", e$message), type = "error", duration = 5)
  })
})

sample_match_validation <- reactive({
  if (is.null(rv$lfq_cols) || is.null(rv$sample_info)) {
    return(list(status = "waiting", message = "Please upload both expression matrix and sample information.",
                matched = character(0), unmatched_info = character(0), unmatched_expr = character(0)))
  }
  expr_col_full <- rv$lfq_cols
  info_names_full <- rownames(rv$sample_info)
  matched_full <- intersect(info_names_full, expr_col_full)
  unmatched_info_full <- setdiff(info_names_full, expr_col_full)
  unmatched_expr_full <- setdiff(expr_col_full, info_names_full)
  matched <- extract_sample_names(matched_full)
  unmatched_info <- extract_sample_names(unmatched_info_full)
  unmatched_expr <- extract_sample_names(unmatched_expr_full)
  if (length(unmatched_info) == 0 && length(unmatched_expr) == 0) {
    return(list(status = "success", message = paste0("All ", length(matched), " samples are successfully matched!"),
                matched = matched, unmatched_info = character(0), unmatched_expr = character(0)))
  } else {
    return(list(status = "warning", message = paste0(length(matched), " samples matched. ",
                                                     length(unmatched_info), " sample(s) in info but not in expression; ",
                                                     length(unmatched_expr), " sample(s) in expression but not in info."),
                matched = matched, unmatched_info = unmatched_info, unmatched_expr = unmatched_expr))
  }
}) %>% bindCache(rv$lfq_cols, rv$sample_info)

output$sample_match_hint <- renderUI({
  req(rv$sample_info)
  div(style = "margin-top: 8px; padding: 8px 12px; background: #fff3cd; border-radius: 6px; color: #856404; font-weight: bold;",
      icon("info-circle"), " Green highlighted samples are matched with the uploaded sample info. Samples without fill color are not matched.")
})

# ==================== 修复：正确的 renderDataTable 参数位置 ====================
# 注意：server = FALSE 必须作为 renderDataTable() 的顶层参数，
# 而不是放在 DT::datatable() 内部。
output$upload_preview <- DT::renderDataTable({
  message("[DEBUG] upload_preview: rv$lfq_cols length = ", length(rv$lfq_cols))
  req(rv$clean_data, rv$lfq_cols)
  df <- rv$clean_data[, rv$lfq_cols, drop = FALSE]
  DT::datatable(df,
                options = list(pageLength = 10, scrollX = TRUE),
                rownames = FALSE)
}, server = FALSE)   # 关键修改：将 server 参数移到这里

output$sample_info_preview <- DT::renderDataTable({
  message("[DEBUG] sample_info_preview triggered")
  req(rv$sample_info)
  df <- rv$sample_info
  df_display <- data.frame(SampleName = rownames(df), df, check.names = FALSE, stringsAsFactors = FALSE)
  DT::datatable(df_display,
                options = list(pageLength = 10, scrollX = TRUE),
                rownames = FALSE)
}, server = FALSE)   # 关键修改：同样移到这里

output$data_summary_ui <- renderUI({
  req(rv$raw_data)
  type_label <- if (input$intensity_type == "LFQ") "LFQ intensity" else "Intensity"
  div(
    p(strong("Dimensions:"), sprintf("%d rows × %d columns", nrow(rv$raw_data), ncol(rv$raw_data))),
    p(strong(paste(type_label, "columns:", sep = " ")), length(rv$lfq_cols)),
    p(strong("Samples:"), length(rv$sample_names))
  )
})

output$detected_samples_ui <- renderUI({
  req(rv$sample_names)
  type_label <- if (input$intensity_type == "LFQ") "LFQ intensity" else "Intensity"
  samples <- rv$sample_names
  validation <- sample_match_validation()
  matched <- if (validation$status %in% c("success", "warning")) validation$matched else character(0)
  tagList(
    h4(icon("vial"), paste(" Detected Samples (", type_label, ")", sep = "")),
    div(style = "max-height: 200px; overflow-y: auto;",
        lapply(samples, function(s) {
          if (s %in% matched) {
            div(class = "sample-item", style = "background: #d4edda; border-color: #c3e6cb;", icon("vial"), " ", s)
          } else {
            div(class = "sample-item", icon("vial"), " ", s)
          }
        })
    )
  )
})

observeEvent(input$reset_all, {
  raw <- rv$raw_data; clean <- rv$clean_data; lfq <- rv$lfq_cols; sn <- rv$sample_names; si <- rv$sample_info
  rv$groups <- list(); rv$comparisons <- list(); rv$analysis_results <- NULL
  rv$pending_duplicate <- NULL; rv$reset_counter <- rv$reset_counter + 1
  rv$comp_id_counter <- 0; rv$group_id_counter <- 0; rv$group_id_map <- list()
  rv$current_profile_protein <- NULL; rv$batch_vector <- NULL; rv$sample_info <- NULL
  manual_sort_active(FALSE)
  cached_sample_info$LFQ <- NULL; cached_sample_info$Intensity <- NULL
  for (name in names(subplot_old_values)) subplot_old_values[[name]] <- NULL
  updateSelectInput(session, "comp_treat", choices = character(0))
  updateSelectInput(session, "comp_ctrl", choices = character(0))
  updateSelectInput(session, "selected_comparison", choices = character(0))
  updateSelectInput(session, "batch_ref_group", choices = character(0))
  updateSelectizeInput(session, "venn_comparisons_select", choices = character(0), selected = character(0))
  updateCheckboxGroupInput(session, "venn_comparisons_checkbox", choices = character(0), selected = character(0))
  updateNumericInput(session, "fc_up", value = 1.2)
  updateNumericInput(session, "fc_down", value = 0.84)
  updateSelectInput(session, "p_cut", selected = "0.05")
  updateNumericInput(session, "min_treat_valid", value = 2)
  updateNumericInput(session, "min_ctrl_valid", value = 2)
  updateNumericInput(session, "min_rep_ttest", value = 2)
  updateNumericInput(session, "min_rep_inc", value = 2)
  updateNumericInput(session, "min_rep_dec", value = 2)
  updateNumericInput(session, "min_unique_pep", value = 2)
  updateNumericInput(session, "point_size", value = 1.8)
  cols <- default_colors()
  colourpicker::updateColourInput(session, "color_up", value = cols$Up)
  colourpicker::updateColourInput(session, "color_down", value = cols$Down)
  colourpicker::updateColourInput(session, "color_increase", value = cols$Increase)
  colourpicker::updateColourInput(session, "color_decrease", value = cols$Decrease)
  colourpicker::updateColourInput(session, "color_ns", value = cols$NS)
  updateRadioButtons(session, "stat_method", selected = "t-test")
  updateNumericInput(session, "replicate_fill_all", value = 2)
  updateTextInput(session, "download_single_title", value = "")
  updateTextInput(session, "combined_plot_title", value = "Combined Volcano Plots")
  updateTextInput(session, "single_plot_title", value = "")
  shinyjs::reset("plot_format")
  updateTextInput(session, "plot_width", value = "10")
  updateTextInput(session, "plot_height", value = "8")
  sub_ids <- grep("^subplot_title_", names(input), value = TRUE)
  for (sid in sub_ids) updateTextInput(session, sid, value = "")
  updateRadioButtons(session, "heatmap_data_source", selected = "LFQ")
  heatmap_raw_groups(NULL)
  rv$raw_data <- raw; rv$clean_data <- clean; rv$lfq_cols <- lfq; rv$sample_names <- sn; rv$sample_info <- si
  if (!is.null(sn) && length(sn) > 0) {
    updateSelectInput(session, "baseline_sample", choices = c("Auto", sn), selected = "Auto")
  } else {
    updateSelectInput(session, "baseline_sample", choices = c("Auto"), selected = "Auto")
  }
  showNotification("All settings reset.", type = "message", duration = 2)
})

get_base_sample <- function() {
  user_sel <- input$baseline_sample
  if (!is.null(user_sel) && user_sel != "Auto" && user_sel %in% rv$sample_names) return(user_sel)
  groups <- rv$groups
  if (length(groups) > 0) {
    wt_names <- grep("WT|Control|CK", names(groups), value = TRUE, ignore.case = TRUE)
    if (length(wt_names) > 0 && length(groups[[wt_names[1]]]) > 0) return(groups[[wt_names[1]]][1])
  }
  comps <- rv$comparisons
  if (length(comps) > 0) {
    ctrl_group <- comps[[1]]$ctrl
    if (!is.null(ctrl_group) && ctrl_group %in% names(groups) && length(groups[[ctrl_group]]) > 0)
      return(groups[[ctrl_group]][1])
  }
  if (length(rv$sample_names) > 0) return(rv$sample_names[1])
  NULL
}

current_baseline <- reactive({ get_base_sample() })

raw_totals <- reactive({
  req(rv$raw_data, rv$lfq_cols)
  totals <- sapply(rv$lfq_cols, function(col) sum(as.numeric(rv$raw_data[[col]]), na.rm = TRUE))
  names(totals) <- extract_sample_names(rv$lfq_cols)
  totals
})

get_analysis_matrix <- reactive({
  proc <- tryCatch(processed_data(), error = function(e) NULL)
  if (!is.null(proc)) return(proc) else return(expression_data())
})

norm_data_before_batch <- reactive({
  mat <- get_analysis_matrix()
  if (is.null(mat)) { showNotification("No expression data available.", type = "error"); return(NULL) }
  base_sample <- current_baseline()
  if (is.null(base_sample)) { showNotification("Unable to determine baseline sample.", type = "error"); return(NULL) }
  sample_short <- extract_sample_names(colnames(mat))
  base_idx <- which(sample_short == base_sample)
  if (length(base_idx) == 0) { showNotification(paste0("Baseline sample '", base_sample, "' not found."), type = "error"); return(NULL) }
  base_sum <- sum(mat[, base_idx], na.rm = TRUE)
  if (base_sum <= 0) { showNotification("Baseline sample total intensity is zero.", type = "error"); return(NULL) }
  norm_mat <- mat
  for (i in seq_len(ncol(mat))) {
    s <- sum(mat[, i], na.rm = TRUE)
    if (s > 0) norm_mat[, i] <- mat[, i] * base_sum / s
    else norm_mat[, i] <- mat[, i]
  }
  norm_prefix <- get_norm_prefix()
  colnames(norm_mat) <- paste0(norm_prefix, sample_short)
  norm_df <- as.data.frame(norm_mat)
  norm_df$`Master protein IDs` <- rownames(norm_mat)
  clean <- rv$clean_data
  if (!is.null(clean)) {
    extra_cols <- intersect(c("Protein IDs", "Majority protein IDs", "Unique peptides", "Fasta headers"), colnames(clean))
    if (length(extra_cols) > 0) {
      idx <- match(norm_df$`Master protein IDs`, clean$`Master protein IDs`)
      for (col in extra_cols) norm_df[[col]] <- clean[[col]][idx]
    }
  }
  norm_df
})

norm_totals <- reactive({
  nd <- norm_data_before_batch()
  if (is.null(nd)) return(NULL)
  norm_prefix <- get_norm_prefix()
  norm_cols <- grep(paste0("^", norm_prefix), colnames(nd), value = TRUE)
  if (length(norm_cols) == 0) return(NULL)
  totals <- sapply(norm_cols, function(col) sum(as.numeric(nd[[col]]), na.rm = TRUE))
  names(totals) <- gsub(paste0("^", norm_prefix), "", norm_cols)
  totals
})

norm_data_full <- reactive({
  norm_data_before_batch()
})

output$norm_comparison_plot <- renderPlotly({
  req(raw_totals())
  raw <- raw_totals()
  nrt <- tryCatch(norm_totals(), error = function(e) NULL)
  samples <- names(raw)
  baseline <- current_baseline()
  p <- plot_ly()
  p <- add_bars(p, x = samples, y = raw[samples], name = "Raw",
                marker = list(color = "steelblue"),
                hovertemplate = paste0("Sample: %{x}<br>Raw: %{y:.0f}<extra></extra>"))
  if (!is.null(nrt)) {
    p <- add_bars(p, x = samples, y = nrt[samples], name = "Normalized",
                  marker = list(color = "darkorange"),
                  hovertemplate = paste0("Sample: %{x}<br>Normalized: %{y:.0f}<extra></extra>"))
    title_text <- paste0("Total Intensity: Raw vs Normalized (Baseline: ", baseline, ")")
  } else {
    title_text <- paste0("Total Intensity (Raw only) (Baseline: ", baseline %||% "N/A", ")")
  }
  p %>% layout(title = title_text, yaxis = list(title = "Total Intensity"),
               xaxis = list(tickangle = -45), legend = list(title = list(text = "Type")), barmode = "group")
})

output$upload_status_ui <- renderUI({
  if (is.null(rv$raw_data)) {
    div(class = "status-badge status-warning", icon("exclamation-triangle"), " No expression file uploaded")
  } else {
    div(class = "status-badge status-success", icon("check-circle"), " Expression matrix uploaded!")
  }
})

observeEvent(input$reset_color, {
  cols <- default_colors()
  colourpicker::updateColourInput(session, "color_up", value = cols$Up)
  colourpicker::updateColourInput(session, "color_down", value = cols$Down)
  colourpicker::updateColourInput(session, "color_increase", value = cols$Increase)
  colourpicker::updateColourInput(session, "color_decrease", value = cols$Decrease)
  colourpicker::updateColourInput(session, "color_ns", value = cols$NS)
  showNotification("Colors reset to defaults.", type = "message", duration = 2)
})

# ==================== 核心：expression_data 反应式（唯一入口） ====================
expression_data <- reactive({
  message("[DEBUG] expression_data triggered from server/data_upload.R")
  req(rv$clean_data)
  if (is.null(rv$lfq_cols) || length(rv$lfq_cols) == 0) {
    validate(need(FALSE, "No intensity columns found. Please upload data first."))
  }
  df <- rv$clean_data
  if (!"Master protein IDs" %in% colnames(df)) {
    validate(need(FALSE, "Master protein IDs column not found in cleaned data."))
  }
  rownames(df) <- as.character(df$`Master protein IDs`)
  df <- df[, rv$lfq_cols, drop = FALSE]
  df <- suppressWarnings(as.data.frame(lapply(df, as.numeric)))
  df[df == 0] <- NA
  if (ncol(df) == 0) {
    validate(need(FALSE, "No intensity columns found."))
  }
  if (nrow(df) == 0) {
    validate(need(FALSE, "No protein rows found."))
  }
  message(sprintf("[DEBUG] expression_data: %d proteins, %d samples", nrow(df), ncol(df)))
  df
})

# 数据质量专用矩阵（使用简短样本名）
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