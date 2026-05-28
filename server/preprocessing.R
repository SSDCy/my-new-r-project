# server/preprocessing.R
# ============================================================
# 注意：expression_data 反应式已在 server/data_upload.R 中定义，
# 本文件直接使用 expression_data()，不再重复定义。
# 样本名标准化统一使用 standardize_sample_name()（定义在 data_upload.R）
# 修复：impute_missing_values 已移至 global.R，此处不再定义。
# 修复：imputation_comparison_data 和 filter_comparison_data 增加一致性检查，
#       但 filter_comparison_data 的检查移到过滤之后，避免误判。
# 修复：batch_pca_plot 增加显式颜色映射，防止点全灰。
# 修复：上传新数据时重置批次校正参数，避免显示旧结果。
# 修复：Comparison 选项卡增加预处理状态判断，显示更明确的提示。
# 修复：缺失值过滤分组模式下，若部分样本未匹配到组信息，自动回退全局模式并提醒用户，
#       同时更新UI下拉框。
# 修复：回退时记录实际过滤模式（本次修复）。
# 修复：切换强度类型或上传新文件时，标记预处理失效，避免热图使用旧数据。
# ============================================================

# ---------- 显示样本信息各列的概要 ----------
output$batch_info_preview <- renderPrint({
  req(rv$sample_info)
  col_data <- rv$sample_info[["Batch"]]
  unique_vals <- unique(col_data)
  cat("样本信息表“Batch”列概览：\n")
  cat("样本总数：", length(col_data), "\n")
  cat("唯一值数量：", length(unique_vals), "\n")
  cat("唯一值列表（前10个）：\n")
  if (length(unique_vals) <= 10) {
    cat(paste(unique_vals, collapse = "\n"), "\n")
  } else {
    cat(paste(head(unique_vals, 10), collapse = "\n"), "\n... 还有", length(unique_vals) - 10, "个\n")
  }
})

# ---------- 快速预设按钮 ----------
observeEvent(input$preset_missing_0.3, {
  updateSliderInput(session, "max_missing_fraction", value = 0.3)
})
observeEvent(input$preset_missing_0.5, {
  updateSliderInput(session, "max_missing_fraction", value = 0.5)
})
observeEvent(input$preset_missing_0.7, {
  updateSliderInput(session, "max_missing_fraction", value = 0.7)
})

# ---------- 动态缺失值过滤效果 ----------
output$missing_filter_effect <- renderPrint({
  req(expression_data())
  missing_per_protein <- rowMeans(is.na(expression_data()))
  filtered <- sum(missing_per_protein > input$max_missing_fraction)
  retained <- nrow(expression_data()) - filtered
  filtered_percent <- round(filtered / nrow(expression_data()) * 100, 1)
  cat("Single-step predicted removal (only missing rate):", filtered, "proteins (", filtered_percent, "%)\n", sep = "")
  cat("Single-step predicted retained (only missing rate):", retained, "proteins\n", sep = "")
  cat("\nNote: This prediction only considers the missing rate filter, without prior Inf/Intensity filters. It is intended as a reference for threshold estimation.\n")
})

# ---------- 强度过滤效果 ----------
output$intensity_filter_effect <- renderPrint({
  message("[DEBUG] intensity_filter_effect called")
  req(expression_data(), input$max_missing_fraction, input$min_intensity)
  data <- expression_data()
  missing_frac <- rowMeans(is.na(data))
  keep <- missing_frac <= input$max_missing_fraction
  data <- data[keep, , drop = FALSE]
  max_int <- apply(data, 1, max, na.rm = TRUE)
  keep_finite <- is.finite(max_int)
  total_before_intensity <- sum(keep_finite)
  data <- data[keep_finite, , drop = FALSE]
  max_int <- max_int[keep_finite]
  min_samples <- input$min_samples_above_intensity %||% 1
  above_thresh_counts <- apply(data, 1, function(x) sum(x > input$min_intensity, na.rm = TRUE))
  below_thresh <- sum(above_thresh_counts < min_samples)
  retained <- total_before_intensity - below_thresh
  threshold_log <- log10(input$min_intensity + 1)
  cat("After missing and Inf filter, ", below_thresh, " proteins are below the intensity threshold in < ", min_samples, " samples (threshold = ",
      input$min_intensity, ", log10 = ", round(threshold_log, 2), ") and will be removed. ",
      "All ", retained, " proteins are retained.\n", sep = "")
  if (below_thresh > 0) {
    cat(sprintf("本次过滤移除了 %d 个蛋白，保留了 %d 个蛋白。\n", below_thresh, retained))
  } else {
    cat("本次过滤未造成有效蛋白损失，所有符合信号质量要求的蛋白均被保留。\n")
  }
})

# ---------- 强度分布图 ----------
output$intensity_dist_plot <- renderPlot({
  req(expression_data(), input$max_missing_fraction)
  threshold <- if (is.null(input$min_intensity) || is.na(input$min_intensity)) 0 else input$min_intensity
  
  data <- expression_data()
  missing_frac <- rowMeans(is.na(data))
  keep <- missing_frac <= input$max_missing_fraction
  data <- data[keep, , drop = FALSE]
  max_int <- apply(data, 1, max, na.rm = TRUE)
  max_int <- max_int[is.finite(max_int)]
  if (length(max_int) == 0) return(NULL)
  
  log_int <- log10(max_int + 1)
  df <- data.frame(log_int = log_int)
  
  total_count <- length(max_int)
  retained_count <- sum(max_int > threshold)
  min_samples <- input$min_samples_above_intensity %||% 1
  
  subtitle_text <- paste0(
    "After missing & Inf filter: ", total_count, " proteins\n",
    "Retained if at least ", min_samples, " sample(s) > ", threshold, " (max-based ref): ", retained_count
  )
  
  p <- ggplot(df, aes(x = log_int)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "steelblue", alpha = 0.6) +
    geom_density(color = "darkorange", linewidth = 1.2) +
    labs(title = "Protein Max Intensity Distribution",
         subtitle = subtitle_text,
         x = "log10(Max Intensity + 1)", y = "Density") +
    theme_bw() +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9),
      plot.margin = margin(t = 10, r = 10, b = 5, l = 10)
    )
  
  if (threshold > 0) {
    threshold_log <- log10(threshold + 1)
    p <- p + geom_vline(xintercept = threshold_log, color = "red", linetype = "dashed", linewidth = 1) +
      annotate("text", x = threshold_log, y = Inf,
               label = paste0("Threshold = ", threshold, " (log10 = ", round(threshold_log, 2), ")"),
               vjust = 2, hjust = -0.1, color = "red", size = 3.5)
  }
  p
})

# ---------- 强度统计 ----------
intensity_stats <- reactive({
  message("[DEBUG] intensity_stats called, using expression_data from data_upload.R")
  req(expression_data(), input$max_missing_fraction)
  data <- expression_data()
  missing_per_protein <- rowMeans(is.na(data))
  keep <- missing_per_protein <= input$max_missing_fraction
  data <- data[keep, , drop = FALSE]
  if (nrow(data) == 0) {
    validate(need(FALSE, "No proteins remain after missing value filter. Lower the threshold."))
  }
  max_intensities <- apply(data, 1, max, na.rm = TRUE)
  keep_finite <- is.finite(max_intensities)
  if (sum(keep_finite) == 0) {
    validate(need(FALSE, "All max intensities are non-finite after filtering."))
  }
  max_intensities <- max_intensities[keep_finite]
  quantiles <- c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1)
  thresholds <- quantile(max_intensities, quantiles, na.rm = TRUE)
  filtered_counts <- sapply(thresholds, function(t) sum(max_intensities < t))
  recommended_threshold <- quantile(max_intensities, 0.15, na.rm = TRUE)
  list(max_intensities = max_intensities,
       thresholds = thresholds,
       quantiles = quantiles,
       filtered_counts = filtered_counts,
       recommended_threshold = recommended_threshold,
       total_proteins = length(max_intensities),
       original_total = nrow(expression_data()),
       missing_filtered = nrow(expression_data()) - nrow(data),
       inf_filtered = sum(!keep_finite))
})

observe({
  req(intensity_stats())
  updateNumericInput(session, "calc_threshold",
                     value = round(intensity_stats()$recommended_threshold))
})

calc_result <- reactive({
  req(intensity_stats(), input$calc_threshold)
  stats <- intensity_stats()
  threshold <- input$calc_threshold
  filtered <- sum(stats$max_intensities < threshold)
  retained <- stats$total_proteins - filtered
  filtered_percent <- round(filtered / stats$total_proteins * 100, 1)
  list(filtered = filtered, retained = retained, filtered_percent = filtered_percent,
       missing_filtered = stats$missing_filtered, inf_filtered = stats$inf_filtered,
       original_total = stats$original_total)
})

output$calc_result <- renderPrint({
  req(calc_result())
  res <- calc_result()
  cat("Original proteins:", res$original_total, "\n")
  cat("Missing value filter removed:", res$missing_filtered, "\n")
  if (res$inf_filtered > 0) cat("Inf value filter removed:", res$inf_filtered, "\n")
  cat("Intensity filter removed:", res$filtered, "(", res$filtered_percent, "%)\n")
  cat("Final predicted retained:", res$retained, "\n")
  cat("\nNote: This calculation already includes all current filter effects.\n")
})

output$pre_raw_summary <- renderPrint({
  tryCatch({
    req(expression_data())
    stats <- intensity_stats()
    cat("Raw data dimensions:", nrow(expression_data()), "proteins,", ncol(expression_data()), "samples\n")
    cat("Total missing values:", sum(is.na(expression_data())), "\n")
    cat("Overall missing ratio:", round(sum(is.na(expression_data()))/(nrow(expression_data())*ncol(expression_data()))*100, 2), "%\n")
    cat("Proteins with any missing:", sum(rowSums(is.na(expression_data())) > 0), "\n\n")
    cat("====================\n")
    cat("Current missing filter setting: max allowed missing =", input$max_missing_fraction, "\n")
    cat("Predicted removal by missing filter:", stats$missing_filtered, "proteins\n")
    if (stats$inf_filtered > 0) cat("Predicted Inf removal:", stats$inf_filtered, "proteins\n")
    cat("After missing/Inf filter:", stats$total_proteins, "proteins\n\n")
    cat("====================\n")
    cat("Intensity distribution (per protein max intensity):\n")
    cat("Min:", format(min(stats$max_intensities, na.rm = TRUE), scientific = FALSE, big.mark = ","), "\n")
    cat("25%:", format(quantile(stats$max_intensities, 0.25, na.rm = TRUE), scientific = FALSE, big.mark = ","), "\n")
    cat("Median:", format(median(stats$max_intensities, na.rm = TRUE), scientific = FALSE, big.mark = ","), "\n")
    cat("75%:", format(quantile(stats$max_intensities, 0.75, na.rm = TRUE), scientific = FALSE, big.mark = ","), "\n")
    cat("Max:", format(max(stats$max_intensities, na.rm = TRUE), scientific = FALSE, big.mark = ","), "\n\n")
    cat("Different quantile thresholds and predicted removal:\n")
    for (i in seq_along(stats$thresholds)) {
      q_percent <- stats$quantiles[i] * 100
      threshold <- format(stats$thresholds[i], scientific = FALSE, big.mark = ",")
      filtered <- stats$filtered_counts[i]
      filtered_percent <- round(filtered / stats$total_proteins * 100, 1)
      retained <- stats$total_proteins - filtered
      cat(q_percent, "% quantile =", threshold,
          ": remove", filtered, "proteins (", filtered_percent, "%), keep", retained, "\n")
    }
    cat("\nRecommended minimum intensity threshold:", format(round(stats$recommended_threshold), scientific = FALSE, big.mark = ","), "\n")
    cat("(This threshold will filter out approximately 15% of low-quality proteins)\n")
  }, error = function(e) {
    cat("Error generating summary:\n", e$message, "\n")
  })
})

output$pre_raw_missing_plot <- renderPlot({
  tryCatch({
    req(expression_data())
    missing_per_protein <- rowMeans(is.na(expression_data()))
    hist(missing_per_protein, breaks = 20,
         main = "Raw Data: Missing Value Proportion Distribution",
         xlab = "Missing Proportion", ylab = "Number of Proteins",
         col = "lightblue", border = "white")
    abline(v = input$max_missing_fraction, col = "red", lwd = 2, lty = 2)
    text(input$max_missing_fraction + 0.02, par("usr")[4]*0.9,
         paste("Threshold =", input$max_missing_fraction), col = "red")
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("Error drawing histogram:\n", e$message), cex = 1.2)
  })
})

# ---------- ComBat 批次校正 ----------
combat_correction <- function(data, batch) {
  if (!requireNamespace("sva", quietly = TRUE))
    stop("sva package required. Run BiocManager::install('sva')")
  batch_factor <- as.factor(batch)
  if (length(levels(batch_factor)) < 2) {
    stop("批次校正需要至少两个不同的批次值。")
  }
  data_matrix <- as.matrix(data)
  message("[DEBUG] combat_correction: input dim = ", nrow(data_matrix), " x ", ncol(data_matrix))
  message("[DEBUG] combat_correction: batches = ", paste(levels(batch_factor), collapse = ", "))
  message("[DEBUG] combat_correction: table(batch) = ", paste(capture.output(print(table(batch_factor))), collapse = "\n"))
  
  if (any(is.na(data_matrix))) {
    message("[DEBUG] ComBat: Detected missing values, applying KNN imputation temporarily for batch correction.")
    data_matrix <- impute_missing_values(data, method = "knn")
    data_matrix <- as.matrix(data_matrix)
    message("[DEBUG] ComBat: missing values after imputation = ", sum(is.na(data_matrix)))
  }
  row_vars <- apply(data_matrix, 1, var)
  zero_var_rows <- sum(row_vars == 0, na.rm = TRUE)
  if (zero_var_rows > 0) {
    data_matrix <- data_matrix[row_vars > 0, , drop = FALSE]
    message("[DEBUG] ComBat: removed ", zero_var_rows, " zero-variance rows")
  }
  
  message("[DEBUG] ComBat: running ComBat...")
  corrected <- sva::ComBat(dat = data_matrix, batch = batch_factor)
  nan_count <- sum(is.na(corrected))
  inf_count <- sum(!is.finite(as.matrix(corrected)))
  if (nan_count + inf_count > 0) {
    message("[DEBUG] ComBat: produced ", nan_count, " NAs and ", inf_count, " Infs; replacing with 1e-4")
    corrected[is.na(corrected)] <- 1e-4
    corrected[!is.finite(as.matrix(corrected))] <- 1e-4
  }
  message("[DEBUG] combat_correction: completed successfully")
  result <- as.data.frame(corrected)
  rownames(result) <- rownames(data_matrix)
  return(result)
}

# ---------- 预处理参数记录 ----------
preprocessing_params <- reactiveValues(
  imputation_method = NULL,
  last_run_time = NULL,
  inf_filtered_count = 0,
  inf_filtered_proteins = character(0),
  batch_performed = FALSE,
  batch_corrected_cols = NULL,
  batch_uncorrected_cols = NULL,
  batch_match_summary = NULL,
  pre_batch_data = NULL,
  post_batch_data = NULL,
  missing_filter_mode = "global",
  missing_filtered_by_group = 0,
  knn_k = 10,
  min_value = 1e-4,
  quantile_prob = 0.01,
  intensity_min_samples = 1,
  missing_filter_fallback = FALSE,
  missing_filter_fallback_unmatched = 0,
  intensity_type_used = NULL
)

# ============ 智能批次诊断 ============
batch_diagnostic <- reactive({
  if (is.null(rv$sample_info) || !"Batch" %in% colnames(rv$sample_info)) {
    return(list(status = "warning", message = "⚠️ 批次诊断失败：样本信息表中缺少 Batch 列。"))
  }
  if (is.null(expression_data())) {
    return(list(status = "info", message = "Expression data not loaded."))
  }
  
  tryCatch({
    mat <- expression_data()
    if (input$max_missing_fraction < 1) {
      missing_frac <- rowMeans(is.na(mat))
      mat <- mat[missing_frac <= input$max_missing_fraction, , drop = FALSE]
    }
    max_int <- apply(mat, 1, max, na.rm = TRUE)
    keep_finite <- is.finite(max_int)
    mat <- mat[keep_finite, , drop = FALSE]
    if (!is.null(input$min_intensity) && !is.na(input$min_intensity) && input$min_intensity > 0) {
      max_int_finite <- apply(mat, 1, max, na.rm = TRUE)
      mat <- mat[max_int_finite > input$min_intensity, , drop = FALSE]
    }
    if (nrow(mat) < 2 || ncol(mat) < 2) {
      return(list(status = "warning", message = "⚠️ 批次诊断失败：过滤后数据量不足。"))
    }
    suppressMessages({
      filled <- impute_missing_values(mat, method = "knn")
    })
    filled <- as.matrix(filled)
    mat_t <- t(log2(filled + 1))
    mat_t <- mat_t[, apply(mat_t, 2, var) > 1e-12, drop = FALSE]
    if (ncol(mat_t) < 2) {
      return(list(status = "warning", message = "⚠️ 批次诊断失败：可变蛋白数量不足。"))
    }
    pca <- prcomp(mat_t, scale. = TRUE)
    pc1 <- pca$x[,1]
    sample_names <- rownames(mat_t)
    info_rows <- standardize_sample_name(rownames(rv$sample_info))
    sample_norm <- standardize_sample_name(sample_names)
    idx <- match(sample_norm, info_rows)
    if (all(is.na(idx))) {
      return(list(status = "warning", message = "⚠️ 批次诊断失败：无法匹配样本到样本信息表。"))
    }
    batch <- rv$sample_info[idx, "Batch"]
    batch <- factor(batch)
    if (length(levels(batch)) < 2) {
      return(list(status = "info", message = "✅ 批次效应检测：仅检测到一个批次，无需校正。"))
    }
    anova_p <- tryCatch(summary(aov(pc1 ~ batch))[[1]]$`Pr(>F)`[1], error = function(e) NA)
    if (!is.na(anova_p) && anova_p < 0.05) {
      return(list(status = "significant", message = "⚠️ 检测到显著批次效应（ANOVA p < 0.05），建议启用 ComBat 校正。"))
    } else {
      return(list(status = "success", message = "✅ 批次效应检测完成：您的数据无显著批次分离，校正前后样本分布无明显变化，无需校正。"))
    }
  }, error = function(e) {
    return(list(status = "warning", message = "⚠️ 批次诊断失败：请检查样本信息表中 Batch 列是否存在，或每个批次样本数量是否≥3 个。"))
  })
})

observe({
  diag <- batch_diagnostic()
  req(diag)
  if (identical(input$imputation_method, "none")) {
    message("[DEBUG] batch auto-enable blocked: imputation method is 'none'")
    return()
  }
  if (diag$status == "significant") {
    if (!isTRUE(input$perform_batch_correction)) {
      message("[DEBUG] auto-enabling batch correction due to significant effect")
      updateCheckboxInput(session, "perform_batch_correction", value = TRUE)
      showNotification("⚠️ 检测到显著批次效应，已自动启用 ComBat 校正。", type = "warning", duration = 6, id = "auto_batch_notif")
    }
  }
})

output$batch_diagnostic_message <- renderPrint({
  cat(batch_diagnostic()$message)
})

output$batch_help_text <- renderUI({
  diag <- batch_diagnostic()
  base_text <- "采用 ComBat 算法校正实验批次偏差。批次列将自动识别为样本信息表中的 'Batch' 列。"
  if (diag$status == "success") {
    HTML(paste(base_text, "<br><b>👉 结合您的数据：校正前后 PCA 显示，批次无显著分离，无明显批次效应，建议取消勾选，避免过度校正。</b>"))
  } else if (diag$status == "significant") {
    HTML(paste(base_text, "<br><b>👉 注意：已自动启用校正，建议运行预处理后查看 Batch Correction 选项卡中的 PCA 对比图。</b>"))
  } else if (diag$status == "warning") {
    HTML(paste(base_text, "<br><b>", diag$message, "</b>"))
  } else {
    HTML(base_text)
  }
})

output$batch_correction_performed <- reactive({
  !is.null(processed_data()) && preprocessing_params$batch_performed
})
outputOptions(output, "batch_correction_performed", suspendWhenHidden = FALSE)

observe({
  message("[DEBUG] conflict check: imputation = ", input$imputation_method, ", batch = ", input$perform_batch_correction)
  if (identical(input$imputation_method, "none")) {
    if (isTRUE(input$perform_batch_correction)) {
      updateCheckboxInput(session, "perform_batch_correction", value = FALSE)
      showNotification("批次校正需要完整数据，已自动取消（不填充时无法校正）。", type = "warning", duration = 6, id = "conflict_notif")
    }
    shinyjs::disable("perform_batch_correction")
  } else {
    shinyjs::enable("perform_batch_correction")
  }
})

apply_missing_filter <- function(data, threshold, mode, sample_info = NULL, sample_names_short = NULL) {
  if (threshold >= 1) return(data)
  fallback_triggered <- FALSE
  unmatched_count <- 0
  if (mode == "group" && !is.null(sample_info) && "Group" %in% colnames(sample_info) && !is.null(sample_names_short)) {
    si <- sample_info
    si$short <- extract_sample_names(rownames(si))
    group_vec <- si$Group[match(sample_names_short, si$short)]
    na_mask <- is.na(group_vec)
    if (any(na_mask)) {
      unmatched_count <- sum(na_mask)
      message(sprintf("[DEBUG] apply_missing_filter: %d samples not matched to any group, falling back to global mode", unmatched_count))
      missing_frac <- rowMeans(is.na(data))
      keep <- missing_frac <= threshold
      fallback_triggered <- TRUE
    } else {
      groups <- unique(group_vec)
      keep <- rep(FALSE, nrow(data))
      for (g in groups) {
        cols_in_group <- which(group_vec == g)
        if (length(cols_in_group) > 0) {
          missing_frac_group <- rowMeans(is.na(data[, cols_in_group, drop = FALSE]))
          keep <- keep | (missing_frac_group <= threshold)
        }
      }
      message("[DEBUG] apply_missing_filter (group): kept ", sum(keep), " out of ", nrow(data))
    }
  } else {
    missing_frac <- rowMeans(is.na(data))
    keep <- missing_frac <= threshold
    message("[DEBUG] apply_missing_filter (global): kept ", sum(keep), " out of ", nrow(data))
  }
  result <- data[keep, , drop = FALSE]
  attr(result, "fallback_triggered") <- fallback_triggered
  attr(result, "unmatched_count") <- unmatched_count
  return(result)
}

apply_intensity_filter <- function(data, threshold, min_samples) {
  if (threshold <= 0) return(data)
  above_counts <- apply(data, 1, function(x) sum(x > threshold, na.rm = TRUE))
  keep <- above_counts >= min_samples
  message("[DEBUG] apply_intensity_filter: threshold = ", threshold, ", min_samples = ", min_samples, ", kept = ", sum(keep), " out of ", nrow(data))
  data[keep, , drop = FALSE]
}

# ============ 核心预处理反应式（修复回退时参数记录） ============
processed_data <- eventReactive(input$run_preprocessing, {
  showNotification("Running preprocessing...", type = "message", duration = NULL, id = "preprocess_notif")
  tryCatch({
    data <- expression_data()
    message("[DEBUG] processed_data: starting preprocessing")
    
    preprocessing_params$intensity_type_used <- input$intensity_type
    message("[DEBUG] processed_data: intensity type recorded: ", input$intensity_type)
    
    filter_mode <- input$missing_filter_mode %||% "global"
    preprocessing_params$missing_filter_mode <- filter_mode
    preprocessing_params$missing_filter_fallback <- FALSE
    preprocessing_params$missing_filter_fallback_unmatched <- 0
    message("[DEBUG] missing filter mode: ", filter_mode)
    
    # 1. 缺失值过滤（可能触发回退）
    data <- apply_missing_filter(
      data, 
      threshold = input$max_missing_fraction,
      mode = filter_mode,
      sample_info = rv$sample_info,
      sample_names_short = rv$sample_names
    )
    fallback_triggered <- attr(data, "fallback_triggered")
    unmatched_count <- attr(data, "unmatched_count")
    if (isTRUE(fallback_triggered)) {
      preprocessing_params$missing_filter_fallback <- TRUE
      preprocessing_params$missing_filter_fallback_unmatched <- unmatched_count
      # ---- 修复：同步更新实际执行的过滤模式 ----
      preprocessing_params$missing_filter_mode <- "global"
      message("[DEBUG] processed_data: missing filter fell back to global, preprocessing_params$missing_filter_mode set to 'global'")
      showNotification(
        paste0("⚠️ 分组缺失值过滤：", unmatched_count, " 个样本未匹配到任何组，已自动切换为全局模式过滤。"),
        type = "warning", duration = 10, id = "missing_fallback_notif"
      )
      # 更新UI下拉框为 global
      updateSelectInput(session, "missing_filter_mode", selected = "global")
      message("[DEBUG] processed_data: UI updated to 'global'")
    }
    
    if (nrow(data) == 0) stop("No proteins left after missing value filter. Relax the threshold.")
    
    # 2. 强度过滤
    max_int <- apply(data, 1, max, na.rm = TRUE)
    keep_finite <- is.finite(max_int)
    preprocessing_params$inf_filtered_count <- sum(!keep_finite)
    preprocessing_params$inf_filtered_proteins <- rownames(data)[!keep_finite]
    data <- data[keep_finite, , drop = FALSE]
    if (nrow(data) == 0) stop("No proteins left after removing Inf values.")
    
    if (!is.null(input$min_intensity) && !is.na(input$min_intensity) && input$min_intensity > 0) {
      min_samples <- input$min_samples_above_intensity %||% 1
      if (is.na(min_samples) || min_samples < 1) min_samples <- 1
      preprocessing_params$intensity_min_samples <- min_samples
      data <- apply_intensity_filter(data, input$min_intensity, min_samples)
      if (nrow(data) == 0) stop("No proteins left after intensity filter. Lower the threshold or reduce the required number of samples.")
    }
    
    # 3. 缺失值填补
    knn_k <- input$knn_k
    if (is.null(knn_k) || is.na(knn_k) || knn_k < 1) knn_k <- 10
    preprocessing_params$knn_k <- knn_k
    
    if (input$imputation_method == "none") {
      preprocessing_params$imputation_method <- "none"
      message("[DEBUG] imputation: skipped (none)")
    } else {
      data <- impute_missing_values(data, method = input$imputation_method, k = knn_k)
      actual_method <- attr(data, "actual_method")
      if (is.null(actual_method)) actual_method <- input$imputation_method
      preprocessing_params$imputation_method <- actual_method
      message("[DEBUG] imputation: performed with method recorded as '", actual_method, "' (user selected: ", input$imputation_method, ")")
    }
    
    # 4. 批次校正
    preprocessing_params$batch_performed <- FALSE
    preprocessing_params$batch_corrected_cols <- NULL
    preprocessing_params$batch_uncorrected_cols <- NULL
    preprocessing_params$batch_match_summary <- NULL
    
    if (isTRUE(input$perform_batch_correction)) {
      message("[DEBUG] batch correction: data is imputed, checking matching...")
      if (!is.null(rv$sample_info) && "Batch" %in% colnames(rv$sample_info)) {
        expr_cols <- colnames(data)
        info_rows <- rownames(rv$sample_info)
        expr_norm <- standardize_sample_name(expr_cols)
        info_norm <- standardize_sample_name(info_rows)
        
        message("[DEBUG] batch correction: matching samples...")
        message("[DEBUG] expr_norm (first 5): ", paste(head(expr_norm, 5), collapse = ", "))
        message("[DEBUG] info_norm (first 5): ", paste(head(info_norm, 5), collapse = ", "))
        
        match_idx <- match(expr_norm, info_norm)
        n_total <- length(expr_cols)
        n_matched <- sum(!is.na(match_idx))
        n_unmatched <- n_total - n_matched
        matched_samples <- expr_cols[!is.na(match_idx)]
        unmatched_samples <- if (n_unmatched > 0) expr_cols[is.na(match_idx)] else character(0)
        
        match_summary <- list(
          total = n_total,
          matched = n_matched,
          unmatched = n_unmatched,
          matched_samples = matched_samples,
          unmatched_samples = unmatched_samples
        )
        preprocessing_params$batch_match_summary <- match_summary
        
        message("[DEBUG] batch correction: total samples = ", n_total, ", matched = ", n_matched, ", unmatched = ", n_unmatched)
        if (n_unmatched > 0) {
          if (length(unmatched_samples) > 5) {
            message("[DEBUG] unmatched samples (first 5): ", paste(head(unmatched_samples, 5), collapse = ", "), " ... and ", length(unmatched_samples)-5, " more")
          } else {
            message("[DEBUG] unmatched samples: ", paste(unmatched_samples, collapse = ", "))
          }
        }
        
        if (n_matched == 0) {
          showNotification(
            paste0("⚠️ 批次校正跳过：没有样本能匹配到样本信息表。请检查样本名是否一致。"),
            type = "error", duration = 10, id = "batch_match_notif"
          )
        } else if (n_unmatched > 0) {
          sample_list <- if (length(unmatched_samples) > 5) {
            paste0(paste(head(unmatched_samples, 5), collapse = ", "), " ...等", length(unmatched_samples), "个")
          } else {
            paste(unmatched_samples, collapse = ", ")
          }
          showNotification(
            paste0("⚠️ 批次校正：", n_unmatched, " 个样本无法匹配（共 ", n_total, " 个）。未匹配样本：", sample_list),
            type = "warning", duration = 10, id = "batch_match_notif"
          )
        }
        
        if (n_matched > 0) {
          keep_samples <- !is.na(match_idx)
          batch_info <- rv$sample_info[match_idx[keep_samples], "Batch"]
          unique_vals <- unique(batch_info)
          if (length(unique_vals) < 2) {
            showNotification("批次校正跳过：匹配后的样本仅有一个批次值。", type = "warning", duration = 8)
          } else {
            preprocessing_params$pre_batch_data <- data[, keep_samples, drop = FALSE]
            data_to_correct <- data[, keep_samples, drop = FALSE]
            corrected_part <- combat_correction(data_to_correct, batch = batch_info)
            corrected_part[is.na(corrected_part)] <- 1e-4
            corrected_part[!is.finite(as.matrix(corrected_part))] <- 1e-4
            data[, keep_samples] <- corrected_part
            preprocessing_params$batch_performed <- TRUE
            preprocessing_params$batch_corrected_cols <- expr_cols[keep_samples]
            preprocessing_params$batch_uncorrected_cols <- if (n_unmatched > 0) unmatched_samples else NULL
            preprocessing_params$post_batch_data <- data
            showNotification(
              paste0("批次校正完成：校正了 ", n_matched, " 个样本", if (n_unmatched > 0) paste0("，", n_unmatched, " 个样本未校正")),
              type = "message", duration = 5, id = "batch_done_notif"
            )
          }
        }
      } else {
        showNotification("批次校正跳过：样本信息表中缺少 'Batch' 列。", type = "warning")
      }
    } else {
      message("[DEBUG] batch correction: not requested")
    }
    
    preprocessing_params$last_run_time <- Sys.time()
    removeNotification("preprocess_notif")
    showNotification("Preprocessing completed! Redirecting to Analysis & Export...", type = "message", duration = 3)
    updateNavbarPage(session, "main_navbar", selected = "plots")
    message("[DEBUG] processed_data: preprocessing finished successfully")
    return(data)
  }, error = function(e) {
    removeNotification("preprocess_notif")
    showNotification(paste("Preprocessing failed:", e$message), type = "error", duration = 10)
    return(NULL)
  })
})

observeEvent(input$expression_file, {
  message("[DEBUG] expression file uploaded, resetting batch params and invalidating processed data")
  preprocessing_params$batch_performed <- FALSE
  preprocessing_params$batch_corrected_cols <- NULL
  preprocessing_params$pre_batch_data <- NULL
  preprocessing_params$post_batch_data <- NULL
  preprocessing_params$batch_match_summary <- NULL
  preprocessing_params$intensity_type_used <- NULL
})

observeEvent(input$intensity_type, {
  message("[DEBUG] intensity type changed, resetting batch params and invalidating processed data")
  preprocessing_params$batch_performed <- FALSE
  preprocessing_params$batch_corrected_cols <- NULL
  preprocessing_params$pre_batch_data <- NULL
  preprocessing_params$post_batch_data <- NULL
  preprocessing_params$batch_match_summary <- NULL
  preprocessing_params$intensity_type_used <- NULL
})

output$pre_processed_summary <- renderPrint({
  req(processed_data())
  cat("Processed data dimensions:", nrow(processed_data()), "proteins,", ncol(processed_data()), "samples\n")
  cat("Remaining missing values:", sum(is.na(processed_data())), "\n")
  cat("Proteins with any missing:", sum(rowSums(is.na(processed_data())) > 0), "\n")
  
  cat("\n==================== 异常值(Inf)过滤详情 ====================\n")
  inf_count <- preprocessing_params$inf_filtered_count
  inf_ids <- preprocessing_params$inf_filtered_proteins
  if (inf_count == 0) {
    cat("未检测到 Inf 异常蛋白\n")
  } else {
    cat("共过滤掉 Inf 异常蛋白：", inf_count, " 个\n", sep = "")
    cat("被过滤的蛋白 ID 列表（前30个）：\n")
    if (length(inf_ids) <= 30) {
      cat(paste(inf_ids, collapse = "\n"), "\n")
    } else {
      cat(paste(head(inf_ids, 30), collapse = "\n"))
      cat("\n...（共计 ", length(inf_ids), " 个）\n", sep = "")
    }
    cat("说明：这些蛋白在所有样本中的强度值均为无穷大（Inf），已自动移除。\n")
  }
  
  if (!is.null(preprocessing_params$imputation_method)) {
    cat("\n====================\n")
    method_display <- switch(preprocessing_params$imputation_method,
                             knn = "K近邻填补法",
                             ppca = "概率主成分分析填补法",
                             minvalue = paste0("最小值填充（固定值 ", preprocessing_params$min_value, "）"),
                             quantile = paste0("分位数填充（概率 ", preprocessing_params$quantile_prob, "）"),
                             "knn (fallback from ppca)" = "K近邻填补法（因PPCA失败自动回退）",
                             none = "无（跳过填充）")
    cat("缺失值处理方式：", method_display, "\n")
    if (grepl("knn", preprocessing_params$imputation_method)) {
      cat("KNN 参数 k = ", preprocessing_params$knn_k, "\n")
    }
    cat("最后运行时间：", format(preprocessing_params$last_run_time, "%Y-%m-%d %H:%M:%S"), "\n")
    
    cat("\n预处理步骤顺序：\n")
    cat("1. 缺失值过滤 (max fraction = ", input$max_missing_fraction, ")\n", sep = "")
    if (preprocessing_params$missing_filter_mode == "group") {
      cat("   过滤模式：分组内缺失率")
      if (isTRUE(preprocessing_params$missing_filter_fallback)) {
        cat("（警告：自动回退为全局模式，", preprocessing_params$missing_filter_fallback_unmatched, " 个样本未匹配到组）")
      }
      cat("\n")
    } else {
      cat("   过滤模式：全局缺失率\n")
    }
    cat("2. 异常值(Inf)过滤（始终执行）\n")
    cat("3. 强度过滤 (min intensity = ", input$min_intensity, ", at least ", preprocessing_params$intensity_min_samples, " samples)\n", sep = "")
    cat("4. 缺失值填补 (", method_display, ")\n", sep = "")
    if (preprocessing_params$batch_performed) {
      cat("5. 批次校正\n")
      if (!is.null(preprocessing_params$batch_corrected_cols)) {
        cat("   已校正的样本列数量：", length(preprocessing_params$batch_corrected_cols), "\n")
      }
      if (!is.null(preprocessing_params$batch_match_summary)) {
        ms <- preprocessing_params$batch_match_summary
        cat("   样本匹配详情：总数 = ", ms$total, ", 已匹配 = ", ms$matched, ", 未匹配 = ", ms$unmatched, "\n")
        if (ms$unmatched > 0) {
          cat("   未匹配样本：", paste(ms$unmatched_samples, collapse = ", "), "\n")
        }
      }
    } else {
      cat("5. 批次校正（未执行）\n")
    }
  }
})

output$pre_processed_missing_plot <- renderPlot({
  req(processed_data())
  missing_per_protein <- rowMeans(is.na(processed_data()))
  hist(missing_per_protein, breaks = 20,
       main = "Processed Data: Missing Value Proportion Distribution",
       xlab = "Missing Proportion", ylab = "Number of Proteins",
       col = "lightgreen", border = "white")
})

output$pre_processed_table <- DT::renderDT({
  req(processed_data())
  df <- processed_data()
  df <- cbind(`Master Protein ID` = rownames(df), df)
  rownames(df) <- NULL
  dt <- DT::datatable(
    df,
    options = list(pageLength = 10, scrollX = TRUE, searchHighlight = TRUE, server = TRUE),
    rownames = FALSE, filter = "top"
  )
  corrected <- preprocessing_params$batch_corrected_cols
  if (!is.null(corrected) && length(corrected) > 0) {
    for (col in corrected) {
      dt <- DT::formatStyle(dt, columns = col, backgroundColor = "#e6f2ff")
    }
  }
  dt
})

# ============ 缺失值过滤前后对比 ============

output$preprocessing_done <- reactive({ !is.null(processed_data()) })
outputOptions(output, "preprocessing_done", suspendWhenHidden = FALSE)

output$imputation_skipped <- reactive({
  !is.null(processed_data()) && preprocessing_params$imputation_method == "none"
})
outputOptions(output, "imputation_skipped", suspendWhenHidden = FALSE)

filter_comparison_data <- reactive({
  req(processed_data())
  before <- expression_data()
  after <- processed_data()
  
  message("[DEBUG] filter_comparison_data: raw before dim = ", nrow(before), "x", ncol(before), 
          "; after dim = ", nrow(after), "x", ncol(after))
  
  filter_mode <- preprocessing_params$missing_filter_mode
  threshold <- input$max_missing_fraction
  message("[DEBUG] filter_comparison_data: mode = ", filter_mode, ", threshold = ", threshold)
  
  # 对原始数据应用相同的过滤
  before <- apply_missing_filter(
    before, threshold, filter_mode,
    sample_info = rv$sample_info,
    sample_names_short = rv$sample_names
  )
  max_int <- apply(before, 1, max, na.rm = TRUE)
  keep_finite <- is.finite(max_int)
  before <- before[keep_finite, , drop = FALSE]
  if (!is.null(input$min_intensity) && !is.na(input$min_intensity) && input$min_intensity > 0) {
    min_samples <- input$min_samples_above_intensity %||% 1
    before <- apply_intensity_filter(before, input$min_intensity, min_samples)
  }
  
  message("[DEBUG] filter_comparison_data: after filters, before dim = ", nrow(before), "x", ncol(before))
  
  # 提取蛋白ID
  protein_ids <- rv$clean_data$`Master protein IDs`
  if (is.null(protein_ids) || length(protein_ids) != nrow(before)) {
    protein_ids <- rownames(before)
  } else {
    protein_ids <- as.character(protein_ids)
  }
  
  before_missing_rate <- rowMeans(is.na(before))
  before_max <- apply(before, 1, max, na.rm = TRUE)
  before_mean <- rowMeans(before, na.rm = TRUE)
  before_median <- apply(before, 1, median, na.rm = TRUE)
  
  pass_inf <- is.finite(before_max)
  pass_missing <- before_missing_rate <= threshold
  pass_intensity <- ifelse(pass_inf, before_max > input$min_intensity, FALSE)
  retained <- pass_inf & pass_missing & pass_intensity
  
  detailed <- data.frame(
    Protein_ID = protein_ids,
    Missing_Rate_Before = round(before_missing_rate, 4),
    Max_Intensity_Before = before_max,
    Mean_Intensity_Before = before_mean,
    Median_Intensity_Before = before_median,
    Pass_Inf_Filter = pass_inf,
    Pass_Missing_Filter = pass_missing,
    Pass_Intensity_Filter = pass_intensity,
    Retained_After_Filter = retained,
    Missing_Fraction_Threshold = threshold,
    Intensity_Threshold = input$min_intensity,
    stringsAsFactors = FALSE
  )
  
  raw_for_export <- cbind(Protein_ID = protein_ids, before)
  
  # 填充缺失值以便可视化比较
  before_impute <- impute_missing_values(before, method = "knn")
  after_impute <- after
  
  message("[DEBUG] filter_comparison_data: before_impute dim = ", nrow(before_impute), "x", ncol(before_impute),
          "; after_impute dim = ", nrow(after_impute), "x", ncol(after_impute))
  
  # 检查过滤并填充后的数据是否维度一致
  if (nrow(before_impute) != nrow(after_impute) || ncol(before_impute) != ncol(after_impute)) {
    message("[DEBUG] filter_comparison_data: imputed dimension mismatch! Check preprocessing parameters.")
    return(NULL)
  }
  
  common_cols <- intersect(colnames(before_impute), colnames(after_impute))
  before_impute <- before_impute[, common_cols, drop = FALSE]
  after_impute <- after_impute[, common_cols, drop = FALSE]
  
  list(
    detailed = detailed,
    raw_before = raw_for_export,
    before_impute = before_impute,
    after_impute = after_impute
  )
})

output$filter_boxplot <- renderPlot({
  if (is.null(processed_data())) {
    plot.new(); text(0.5, 0.5, "Please run preprocessing first.", cex = 1.2)
    return()
  }
  dat <- filter_comparison_data()
  if (is.null(dat)) {
    plot.new(); text(0.5, 0.5, "Data mismatch. Please re-run preprocessing with current parameters.", cex = 1.2)
    return()
  }
  before_log <- log2(as.matrix(dat$before_impute) + 1)
  after_log <- log2(as.matrix(dat$after_impute) + 1)
  
  before_df <- data.frame(Sample = rep(colnames(before_log), each = nrow(before_log)),
                          Intensity = as.vector(before_log),
                          Stage = "Before Filtering")
  after_df <- data.frame(Sample = rep(colnames(after_log), each = nrow(after_log)),
                         Intensity = as.vector(after_log),
                         Stage = "After Filtering")
  plot_df <- rbind(before_df, after_df)
  plot_df$Stage <- factor(plot_df$Stage, levels = c("Before Filtering", "After Filtering"))
  
  ggplot(plot_df, aes(x = Sample, y = Intensity, fill = Stage)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
    scale_fill_manual(values = c("Before Filtering" = "#3498db", "After Filtering" = "#2ecc71")) +
    labs(title = "Intensity Distribution (log2) Before and After Filtering",
         y = "log2(Intensity)", x = "") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
})

safe_pca <- function(mat, scale = TRUE) {
  mat <- as.matrix(mat)
  mat[!is.finite(mat)] <- NA
  mat <- mat[complete.cases(mat), , drop = FALSE]
  if (nrow(mat) < 2) return(NULL)
  row_vars <- apply(mat, 1, var, na.rm = TRUE)
  keep <- row_vars > 1e-12
  mat <- mat[keep, , drop = FALSE]
  if (nrow(mat) < 2) return(NULL)
  tryCatch(prcomp(mat, scale. = scale), error = function(e) NULL)
}

output$filter_pca_plot <- renderPlot({
  if (is.null(processed_data())) {
    plot.new(); text(0.5, 0.5, "Please run preprocessing first.", cex = 1.2)
    return()
  }
  dat <- filter_comparison_data()
  if (is.null(dat)) {
    plot.new(); text(0.5, 0.5, "Data mismatch. Please re-run preprocessing.", cex = 1.2)
    return()
  }
  after_ids <- rownames(dat$after_impute)
  before_sub <- dat$before_impute[after_ids, , drop = FALSE]
  after_sub <- dat$after_impute[after_ids, , drop = FALSE]
  
  if (nrow(before_sub) < 3) {
    plot.new()
    text(0.5, 0.5, "Not enough common proteins for PCA.")
    return()
  }
  
  before_log <- log2(as.matrix(before_sub) + 1)
  after_log <- log2(as.matrix(after_sub) + 1)
  
  pca_before <- safe_pca(before_log)
  pca_after <- safe_pca(after_log)
  
  if (is.null(pca_before) || is.null(pca_after)) {
    plot.new()
    text(0.5, 0.5, "PCA failed due to insufficient variability or missing values.")
    return()
  }
  
  var_before <- round(pca_before$sdev^2 / sum(pca_before$sdev^2) * 100, 1)
  var_after <- round(pca_after$sdev^2 / sum(pca_after$sdev^2) * 100, 1)
  
  df_before <- data.frame(PC1 = pca_before$x[,1], PC2 = pca_before$x[,2], Stage = "Before")
  df_after <- data.frame(PC1 = pca_after$x[,1], PC2 = pca_after$x[,2], Stage = "After")
  pca_df <- rbind(df_before, df_after)
  
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Stage)) +
    geom_point(alpha = 0.6, size = 2) +
    stat_ellipse(type = "norm", level = 0.95) +
    scale_color_manual(values = c("Before" = "#3498db", "After" = "#2ecc71")) +
    labs(title = "PCA: Before vs After Filtering (common proteins)",
         x = paste0("PC1 (Before: ", var_before[1], "%, After: ", var_after[1], "%)"),
         y = paste0("PC2 (Before: ", var_before[2], "%, After: ", var_after[2], "%)")) +
    theme_bw() +
    theme(legend.position = "bottom")
})

output$filter_summary_table <- DT::renderDT({
  if (is.null(processed_data())) {
    return(DT::datatable(data.frame(Message = "Please run preprocessing first.")))
  }
  dat <- filter_comparison_data()
  if (is.null(dat)) {
    return(DT::datatable(data.frame(Message = "Data mismatch. Please re-run preprocessing.")))
  }
  detailed <- dat$detailed
  DT::datatable(detailed, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
})

output$download_filter_table <- downloadHandler(
  filename = function() paste0("Filter_Comparison_", Sys.Date(), ".xlsx"),
  content = function(file) {
    if (is.null(processed_data())) {
      showNotification("Please run preprocessing first.", type = "error")
      return()
    }
    dat <- filter_comparison_data()
    if (is.null(dat)) {
      showNotification("Data mismatch. Please re-run preprocessing.", type = "error")
      return()
    }
    detailed <- dat$detailed
    raw_before <- dat$raw_before
    
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "Raw Data Before Filtering")
    openxlsx::writeData(wb, "Raw Data Before Filtering", raw_before)
    
    total_proteins <- nrow(detailed)
    inf_removed <- sum(!detailed$Pass_Inf_Filter)
    missing_removed <- sum(detailed$Pass_Inf_Filter & !detailed$Pass_Missing_Filter)
    intensity_removed <- sum(detailed$Pass_Inf_Filter & detailed$Pass_Missing_Filter & !detailed$Pass_Intensity_Filter)
    retained <- sum(detailed$Retained_After_Filter)
    
    summary_rows <- data.frame(
      Metric = c("Total Proteins (Before Filtering)", "Inf Value Filter Removed",
                 "Missing Rate Filter Removed", "Intensity Filter Removed", "Final Retained Proteins"),
      Count = c(total_proteins, inf_removed, missing_removed, intensity_removed, retained),
      stringsAsFactors = FALSE
    )
    
    openxlsx::addWorksheet(wb, "Protein Details")
    openxlsx::writeData(wb, "Protein Details", summary_rows, startRow = 1, startCol = 1, colNames = TRUE)
    detail_start_row <- nrow(summary_rows) + 3
    openxlsx::writeData(wb, "Protein Details", detailed, startRow = detail_start_row, startCol = 1, colNames = TRUE)
    
    now_time <- Sys.time()
    uploaded_name <- if (!is.null(input$expression_file)) input$expression_file$name else "Unknown"
    n_samples <- length(rv$sample_names)
    if (is.null(n_samples)) n_samples <- ncol(dat$before_impute)
    
    miss_threshold <- input$max_missing_fraction
    inten_threshold <- input$min_intensity
    
    log_items <- c(
      "Experiment Name", "Number of Samples", "Analysis Time", "Analyst", "",
      "--- Filtering Procedure (applied in order) ---",
      "Step 1: Inf/Abnormal Value Filter", "  Rule", "  Filtered Protein Count",
      "Step 2: Missing Rate Filter", "  Formula (Missing Rate = missing samples / total samples)",
      "  Threshold (max allowed fraction)", "  Filtered Protein Count",
      "Step 3: Intensity Filter", "  Formula (Intensity = max expression value across all samples)",
      "  Threshold (min intensity)", "  Filtered Protein Count", "",
      "--- Column Definitions ---",
      "Missing_Rate_Before: Missing rate = missing sample count / total sample count",
      "Max_Intensity_Before: Maximum intensity among all samples for each protein",
      "Mean_Intensity_Before: Average intensity across samples",
      "Median_Intensity_Before: Median intensity across samples",
      "Pass_Inf_Filter: TRUE if max intensity is finite (not Inf/NaN)",
      "Pass_Missing_Filter: TRUE if missing rate <= threshold",
      "Pass_Intensity_Filter: TRUE if max intensity > intensity threshold",
      "Retained_After_Filter: TRUE if protein meets ALL three conditions",
      "Missing_Fraction_Threshold: user-set missing fraction cutoff",
      "Intensity_Threshold: user-set minimum intensity cutoff"
    )
    
    log_values <- c(
      uploaded_name, as.character(n_samples), format(now_time, "%Y-%m-%d %H:%M:%S"), "Not provided", "", "", "",
      "Proteins with non-finite (Inf/NaN) maximum intensity are removed",
      as.character(inf_removed), "", "", as.character(miss_threshold), as.character(missing_removed),
      "", "", as.character(inten_threshold), as.character(intensity_removed), "", "",
      "", "", "", "", "", "", "", "", "", ""
    )
    
    stopifnot(length(log_items) == length(log_values))
    log_df <- data.frame(Item = log_items, Value = log_values, stringsAsFactors = FALSE)
    
    openxlsx::addWorksheet(wb, "Filtering Log")
    openxlsx::writeData(wb, "Filtering Log", log_df)
    openxlsx::setColWidths(wb, "Filtering Log", cols = 1, widths = 60)
    openxlsx::setColWidths(wb, "Filtering Log", cols = 2, widths = 40)
    
    openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
  }
)

# ============ 缺失值填补前后对比 ============

pre_imputation_matrix <- reactive({
  req(processed_data())
  data <- expression_data()
  filter_mode <- preprocessing_params$missing_filter_mode
  threshold <- input$max_missing_fraction
  message("[DEBUG] pre_imputation_matrix: mode = ", filter_mode, ", threshold = ", threshold)
  data <- apply_missing_filter(
    data, threshold, filter_mode,
    sample_info = rv$sample_info,
    sample_names_short = rv$sample_names
  )
  max_int <- apply(data, 1, max, na.rm = TRUE)
  keep_finite <- is.finite(max_int)
  data <- data[keep_finite, , drop = FALSE]
  if (nrow(data) > 0 && !is.null(input$min_intensity) && !is.na(input$min_intensity) && input$min_intensity > 0) {
    min_samples <- input$min_samples_above_intensity %||% 1
    data <- apply_intensity_filter(data, input$min_intensity, min_samples)
  }
  data
})

imputation_comparison_data <- reactive({
  req(processed_data(), pre_imputation_matrix())
  before_imp <- pre_imputation_matrix()
  after_imp <- processed_data()
  
  if (nrow(before_imp) != nrow(after_imp) || ncol(before_imp) != ncol(after_imp)) {
    message("[DEBUG] imputation_comparison_data: dimension mismatch. Before: ", nrow(before_imp), "x", ncol(before_imp), "; After: ", nrow(after_imp), "x", ncol(after_imp))
    return(NULL)
  }
  
  message("[DEBUG] imputation_comparison_data: before rows = ", nrow(before_imp), ", after rows = ", nrow(after_imp))
  
  all_ids <- rv$clean_data$`Master protein IDs`
  expr_ids <- rownames(expression_data())
  before_ids <- rownames(before_imp)
  idx <- match(before_ids, expr_ids)
  protein_ids <- all_ids[idx]
  if (any(is.na(protein_ids))) {
    protein_ids <- before_ids
  }
  
  missing_before <- rowSums(is.na(before_imp))
  missing_rate_before <- rowMeans(is.na(before_imp))
  mean_before <- rowMeans(before_imp, na.rm = TRUE)
  median_before <- apply(before_imp, 1, median, na.rm = TRUE)
  
  mean_after <- rowMeans(as.matrix(after_imp), na.rm = TRUE)
  median_after <- apply(after_imp, 1, median, na.rm = TRUE)
  
  detailed <- data.frame(
    Protein_ID = protein_ids,
    Missing_Count_Before = missing_before,
    Missing_Rate_Before = round(missing_rate_before, 4),
    Mean_Intensity_Before = mean_before,
    Median_Intensity_Before = median_before,
    Mean_Intensity_After = mean_after,
    Median_Intensity_After = median_after,
    stringsAsFactors = FALSE
  )
  
  before_vis <- impute_missing_values(before_imp, method = "knn")
  after_vis <- after_imp
  common_cols <- intersect(colnames(before_vis), colnames(after_vis))
  before_vis <- before_vis[, common_cols, drop = FALSE]
  after_vis <- after_vis[, common_cols, drop = FALSE]
  
  total_missing_before <- sum(missing_before)
  missing_after <- sum(is.na(after_imp))
  method <- preprocessing_params$imputation_method
  params <- if (grepl("knn", method)) {
    paste0("k = ", preprocessing_params$knn_k)
  } else if (method == "ppca") {
    "nPcs = 2, scale = 'uv', center = TRUE"
  } else {
    method
  }
  
  list(
    detailed = detailed,
    before_vis = before_vis,
    after_vis = after_vis,
    total_missing_before = total_missing_before,
    missing_after = missing_after,
    method = method,
    params = params,
    n_proteins = nrow(before_imp),
    n_samples = ncol(before_imp)
  )
})

output$imputation_stats_text <- renderPrint({
  if (is.null(processed_data())) {
    cat("Please run preprocessing first.\n")
    return()
  }
  dat <- imputation_comparison_data()
  if (is.null(dat)) {
    cat("Data mismatch. Please re-run preprocessing.\n")
    return()
  }
  cat("Total missing values before imputation:", dat$total_missing_before, "\n")
  cat("Missing values after imputation:", dat$missing_after, "\n")
  cat("Imputation method:", dat$method, "\n")
  cat("Parameters:", dat$params, "\n")
})

output$imputation_boxplot <- renderPlot({
  if (is.null(processed_data())) {
    plot.new(); text(0.5, 0.5, "Please run preprocessing first.", cex = 1.2)
    return()
  }
  dat <- imputation_comparison_data()
  if (is.null(dat)) {
    plot.new(); text(0.5, 0.5, "Data mismatch. Please re-run preprocessing.", cex = 1.2)
    return()
  }
  before_log <- log2(as.matrix(dat$before_vis) + 1)
  after_log <- log2(as.matrix(dat$after_vis) + 1)
  
  before_df <- data.frame(Sample = rep(colnames(before_log), each = nrow(before_log)),
                          Intensity = as.vector(before_log),
                          Stage = "Before Imputation")
  after_df <- data.frame(Sample = rep(colnames(after_log), each = nrow(after_log)),
                         Intensity = as.vector(after_log),
                         Stage = "After Imputation")
  plot_df <- rbind(before_df, after_df)
  plot_df$Stage <- factor(plot_df$Stage, levels = c("Before Imputation", "After Imputation"))
  
  sample_names <- colnames(before_log)
  group_levels <- gsub("\\..*", "", sample_names)
  plot_df$Group <- rep(group_levels, each = nrow(before_log))
  
  n_proteins <- nrow(before_log)
  n_samples <- length(sample_names)
  method_display <- switch(dat$method, knn = "KNN", ppca = "PPCA", "knn (fallback from ppca)" = "KNN (fallback)", none = "None", dat$method)
  subtitle_text <- paste0(method_display, " Imputation | ", n_proteins, " proteins | ", n_samples, " samples")
  
  ggplot(plot_df, aes(x = Sample, y = Intensity, fill = Stage)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
    labs(title = "Intensity Distribution (log2) Before and After Imputation",
         subtitle = subtitle_text,
         y = "log2(Intensity)", x = "") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          legend.position = "bottom")
})

output$imputation_pca_plot <- renderPlot({
  if (is.null(processed_data())) {
    plot.new(); text(0.5, 0.5, "Please run preprocessing first.", cex = 1.2)
    return()
  }
  dat <- imputation_comparison_data()
  if (is.null(dat)) {
    plot.new(); text(0.5, 0.5, "Data mismatch. Please re-run preprocessing.", cex = 1.2)
    return()
  }
  before_sub <- dat$before_vis
  after_sub <- dat$after_vis
  
  if (nrow(before_sub) < 3) {
    plot.new()
    text(0.5, 0.5, "Not enough proteins for PCA.")
    return()
  }
  
  before_log <- log2(as.matrix(before_sub) + 1)
  after_log <- log2(as.matrix(after_sub) + 1)
  
  pca_before <- safe_pca(before_log)
  pca_after <- safe_pca(after_log)
  
  if (is.null(pca_before) || is.null(pca_after)) {
    plot.new()
    text(0.5, 0.5, "PCA failed due to insufficient variability or missing values.")
    return()
  }
  
  var_before <- round(pca_before$sdev^2 / sum(pca_before$sdev^2) * 100, 1)
  var_after <- round(pca_after$sdev^2 / sum(pca_after$sdev^2) * 100, 1)
  
  df_before <- data.frame(PC1 = pca_before$x[,1], PC2 = pca_before$x[,2], Stage = "Before Imputation")
  df_after <- data.frame(PC1 = pca_after$x[,1], PC2 = pca_after$x[,2], Stage = "After Imputation")
  pca_df <- rbind(df_before, df_after)
  
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Stage)) +
    geom_point(alpha = 0.6, size = 2) +
    stat_ellipse(type = "norm", level = 0.95) +
    scale_color_manual(values = c("Before Imputation" = "#e67e22", "After Imputation" = "#9b59b6")) +
    labs(title = "PCA: Before vs After Imputation",
         x = paste0("PC1 (Before: ", var_before[1], "%, After: ", var_after[1], "%)"),
         y = paste0("PC2 (Before: ", var_before[2], "%, After: ", var_after[2], "%)")) +
    theme_bw() +
    theme(legend.position = "bottom")
})

output$imputation_qq_plot <- renderPlot({
  if (is.null(processed_data())) {
    plot.new(); text(0.5, 0.5, "Please run preprocessing first.", cex = 1.2)
    return()
  }
  dat <- imputation_comparison_data()
  if (is.null(dat)) {
    plot.new(); text(0.5, 0.5, "Data mismatch. Please re-run preprocessing.", cex = 1.2)
    return()
  }
  before_log <- log2(as.matrix(dat$before_vis) + 1)
  after_log <- log2(as.matrix(dat$after_vis) + 1)
  
  qq_before <- data.frame(Intensity = as.vector(before_log), Stage = "Before Imputation")
  qq_after <- data.frame(Intensity = as.vector(after_log), Stage = "After Imputation")
  qq_df <- rbind(qq_before, qq_after)
  
  ggplot(qq_df, aes(sample = Intensity, color = Stage)) +
    stat_qq(size = 0.5, alpha = 0.5) +
    stat_qq_line() +
    scale_color_manual(values = c("Before Imputation" = "#e67e22", "After Imputation" = "#9b59b6")) +
    labs(title = "Q-Q Plot: Normality Check Before and After Imputation",
         subtitle = paste("Log2 Intensity Distribution"),
         x = "Theoretical Quantiles",
         y = "Sample Quantiles") +
    theme_bw() +
    theme(legend.position = "bottom")
})

output$imputation_summary_table <- DT::renderDT({
  if (is.null(processed_data())) {
    return(DT::datatable(data.frame(Message = "Please run preprocessing first.")))
  }
  dat <- imputation_comparison_data()
  if (is.null(dat)) {
    return(DT::datatable(data.frame(Message = "Data mismatch. Please re-run preprocessing.")))
  }
  DT::datatable(dat$detailed,
                options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
})

output$download_imputation_table <- downloadHandler(
  filename = function() paste0("Imputation_Comparison_", Sys.Date(), ".xlsx"),
  content = function(file) {
    if (is.null(processed_data())) {
      showNotification("Please run preprocessing first.", type = "error")
      return()
    }
    dat <- imputation_comparison_data()
    if (is.null(dat)) {
      showNotification("Data mismatch. Please re-run preprocessing.", type = "error")
      return()
    }
    detailed <- dat$detailed
    
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "Imputation Details")
    
    total <- nrow(detailed)
    na_count <- sum(detailed$Missing_Count_Before > 0)
    na_pct <- round(na_count / total * 100, 1)
    total_missing <- dat$total_missing_before
    total_cells <- dat$n_proteins * dat$n_samples
    missing_pct <- round(total_missing / total_cells * 100, 1)
    
    summary_labels <- c(
      paste0("Total Proteins: ", total),
      paste0("Proteins with NA before imputation: ", na_count, " (", na_pct, "%)"),
      paste0("Total missing values before imputation: ", total_missing, " (", total_missing, "/(", dat$n_proteins, "×", dat$n_samples, ") = ", missing_pct, "%)"),
      paste0("Missing values after imputation: ", dat$missing_after),
      paste0("Imputation method: ", dat$method),
      paste0("Parameters: ", dat$params)
    )
    summary_rows <- data.frame(Metric = summary_labels, stringsAsFactors = FALSE)
    
    openxlsx::writeData(wb, "Imputation Details", summary_rows, startRow = 1, startCol = 1, colNames = TRUE)
    openxlsx::writeData(wb, "Imputation Details", detailed, startRow = nrow(summary_rows) + 3, startCol = 1, colNames = TRUE)
    
    now_time <- Sys.time()
    log_items <- c(
      "Experiment Name", "Analysis Time", "Imputation Method", "Parameters",
      "", "--- Column Definitions ---",
      "Missing_Count_Before: Number of missing values per protein before imputation",
      "Missing_Rate_Before: Missing rate before imputation",
      "Mean_Intensity_Before: Average intensity before imputation (NA ignored)",
      "Median_Intensity_Before: Median intensity before imputation (NA ignored)",
      "Mean_Intensity_After: Average intensity after imputation",
      "Median_Intensity_After: Median intensity after imputation"
    )
    log_values <- c(
      if (!is.null(input$expression_file)) input$expression_file$name else "Unknown",
      format(now_time, "%Y-%m-%d %H:%M:%S"),
      dat$method,
      dat$params,
      "", "", "", "", "", "", "", ""
    )
    log_df <- data.frame(Item = log_items, Value = log_values, stringsAsFactors = FALSE)
    
    openxlsx::addWorksheet(wb, "Imputation Log")
    openxlsx::writeData(wb, "Imputation Log", log_df)
    openxlsx::setColWidths(wb, "Imputation Log", cols = 1, widths = 50)
    openxlsx::setColWidths(wb, "Imputation Log", cols = 2, widths = 40)
    
    openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
  }
)

# ============ 跳过填充时的专属图表 ============

output$missing_heatmap_skipped <- renderPlot({
  req(pre_imputation_matrix())
  mat <- pre_imputation_matrix()
  miss_mat <- is.na(as.matrix(mat)) * 1
  rownames(miss_mat) <- rownames(mat)
  colnames(miss_mat) <- colnames(mat)
  
  sample_groups <- gsub("\\..*", "", colnames(mat))
  sample_groups <- factor(sample_groups, levels = unique(sample_groups))
  ann_col <- data.frame(Group = sample_groups, row.names = colnames(mat))
  
  group_levels <- levels(sample_groups)
  if (length(group_levels) <= 8) {
    group_color_vec <- RColorBrewer::brewer.pal(length(group_levels), "Set2")
  } else {
    group_color_vec <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(group_levels))
  }
  names(group_color_vec) <- group_levels
  ann_colors <- list(Group = group_color_vec)
  
  if (nrow(miss_mat) > 1000) {
    set.seed(123)
    miss_mat <- miss_mat[sample(1:nrow(miss_mat), 1000), ]
  }
  
  pheatmap::pheatmap(miss_mat,
                     color = c("#3498db", "#e74c3c"),
                     legend_breaks = c(0, 1),
                     legend_labels = c("Present", "Missing"),
                     cluster_rows = TRUE,
                     cluster_cols = TRUE,
                     show_rownames = FALSE,
                     show_colnames = TRUE,
                     annotation_col = ann_col,
                     annotation_colors = ann_colors,
                     main = "Missing Value Heatmap (Blue = Present, Red = Missing)",
                     fontsize_col = 8)
})

output$valid_barplot_skipped <- renderPlot({
  req(pre_imputation_matrix())
  mat <- pre_imputation_matrix()
  valid_pct <- (1 - colMeans(is.na(mat))) * 100
  df <- data.frame(Sample = names(valid_pct), ValidPct = valid_pct)
  df$Sample <- factor(df$Sample, levels = df$Sample)
  
  ggplot(df, aes(x = Sample, y = ValidPct)) +
    geom_col(fill = "#3498db", alpha = 0.8) +
    geom_hline(yintercept = 70, color = "red", linetype = "dashed", linewidth = 1) +
    annotate("text", x = 1, y = 72, label = "70% threshold", color = "red", hjust = 0) +
    labs(title = "Valid Values per Sample", y = "Valid Percentage (%)", x = "") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
})

output$missing_summary_table_skipped <- renderTable({
  req(pre_imputation_matrix())
  mat <- pre_imputation_matrix()
  n_prot <- nrow(mat)
  n_samp <- ncol(mat)
  overall_na <- round(mean(is.na(mat)) * 100, 1)
  na_proteins <- sum(rowSums(is.na(mat)) > 0)
  avg_valid <- round(mean(1 - colMeans(is.na(mat))) * 100, 1)
  
  data.frame(
    Metric = c("Total Proteins", "Total Samples", "Overall Missing Rate",
               "Proteins with Missing Values", "Average Valid Value % per Sample"),
    Value = c(n_prot, n_samp, paste0(overall_na, "%"), na_proteins, paste0(avg_valid, "%")),
    Description = c("Number of proteins after filtering",
                    "Number of samples",
                    "Percentage of all protein-sample pairs with NA",
                    "Number of proteins with at least one missing value",
                    "Mean of valid value percentages across samples"),
    stringsAsFactors = FALSE
  )
}, striped = TRUE, bordered = TRUE, width = "100%")

# ============ 批次校正前后 PCA 对比（修复颜色） ============
batch_comparison_pca <- reactive({
  req(processed_data(), preprocessing_params$batch_performed)
  before <- preprocessing_params$pre_batch_data
  after <- preprocessing_params$post_batch_data
  if (is.null(before) || is.null(after)) return(NULL)
  
  common_cols <- intersect(colnames(before), colnames(after))
  before <- before[, common_cols, drop = FALSE]
  after <- after[, common_cols, drop = FALSE]
  
  before_t <- t(log2(as.matrix(before) + 1))
  after_t <- t(log2(as.matrix(after) + 1))
  
  before_ok <- apply(before_t, 2, function(x) all(is.finite(x)))
  after_ok <- apply(after_t, 2, function(x) all(is.finite(x)))
  keep_cols <- which(before_ok & after_ok)
  if (length(keep_cols) < 2) return(NULL)
  before_t <- before_t[, keep_cols, drop = FALSE]
  after_t <- after_t[, keep_cols, drop = FALSE]
  
  before_var <- apply(before_t, 2, var)
  after_var <- apply(after_t, 2, var)
  keep_var <- which(before_var > 1e-12 & after_var > 1e-12)
  if (length(keep_var) < 2) return(NULL)
  before_t <- before_t[, keep_var, drop = FALSE]
  after_t <- after_t[, keep_var, drop = FALSE]
  
  pca_before <- tryCatch(prcomp(before_t, scale. = TRUE), error = function(e) NULL)
  pca_after <- tryCatch(prcomp(after_t, scale. = TRUE), error = function(e) NULL)
  if (is.null(pca_before) || is.null(pca_after)) return(NULL)
  
  sample_info_short <- rv$sample_info
  rownames(sample_info_short) <- standardize_sample_name(rownames(sample_info_short))
  before_norm <- standardize_sample_name(colnames(before))
  after_norm <- standardize_sample_name(colnames(after))
  batch_before <- sample_info_short[before_norm, "Batch"]
  batch_after <- sample_info_short[after_norm, "Batch"]
  
  list(
    pca_before = pca_before,
    pca_after = pca_after,
    batch_before = batch_before,
    batch_after = batch_after,
    var_before = round(pca_before$sdev^2 / sum(pca_before$sdev^2) * 100, 1),
    var_after = round(pca_after$sdev^2 / sum(pca_after$sdev^2) * 100, 1)
  )
})

output$batch_pca_plot <- renderPlot({
  dat <- batch_comparison_pca()
  if (is.null(dat)) {
    plot.new(); text(0.5, 0.5, "PCA not available")
    return()
  }
  df_before <- data.frame(PC1 = dat$pca_before$x[,1], PC2 = dat$pca_before$x[,2],
                          Batch = dat$batch_before, Stage = "Before Correction")
  df_after <- data.frame(PC1 = dat$pca_after$x[,1], PC2 = dat$pca_after$x[,2],
                         Batch = dat$batch_after, Stage = "After Correction")
  pca_df <- rbind(df_before, df_after)
  pca_df$Stage <- factor(pca_df$Stage, levels = c("Before Correction", "After Correction"))
  
  batches <- unique(na.omit(c(as.character(dat$batch_before), as.character(dat$batch_after))))
  if (length(batches) > 8) {
    batch_colors <- setNames(rainbow(length(batches)), batches)
  } else {
    batch_colors <- setNames(RColorBrewer::brewer.pal(length(batches), "Set1"), batches)
  }
  
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(type = "norm", level = 0.95) +
    scale_color_manual(values = batch_colors) +
    facet_wrap(~ Stage) +
    labs(title = "PCA Before and After ComBat Batch Correction",
         x = paste0("PC1 (Before: ", dat$var_before[1], "%, After: ", dat$var_after[1], "%)"),
         y = paste0("PC2 (Before: ", dat$var_before[2], "%, After: ", dat$var_after[2], "%)")) +
    theme_bw() + theme(legend.position = "bottom")
})

output$batch_pca_interpretation <- renderUI({
  diag <- batch_diagnostic()
  if (diag$status == "significant") {
    p("图注：红色 = Batch1，青色 = Batch2。校正前后批次样本分布出现分离趋势，提示存在批次效应，建议启用校正。")
  } else {
    p("图注：红色 = Batch1，青色 = Batch2。校正前后批次样本分布无明显聚类，说明无显著批次效应，无需校正。")
  }
})

output$download_batch_pca <- downloadHandler(
  filename = function() "batch_correction_pca.png",
  content = function(file) {
    dat <- batch_comparison_pca()
    if (is.null(dat)) return()
    df_before <- data.frame(PC1 = dat$pca_before$x[,1], PC2 = dat$pca_before$x[,2],
                            Batch = dat$batch_before, Stage = "Before Correction")
    df_after <- data.frame(PC1 = dat$pca_after$x[,1], PC2 = dat$pca_after$x[,2],
                           Batch = dat$batch_after, Stage = "After Correction")
    pca_df <- rbind(df_before, df_after)
    pca_df$Stage <- factor(pca_df$Stage, levels = c("Before Correction", "After Correction"))
    
    batches <- unique(na.omit(c(as.character(dat$batch_before), as.character(dat$batch_after))))
    if (length(batches) > 8) {
      batch_colors <- setNames(rainbow(length(batches)), batches)
    } else {
      batch_colors <- setNames(RColorBrewer::brewer.pal(length(batches), "Set1"), batches)
    }
    
    plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
      geom_point(size = 3, alpha = 0.8) +
      stat_ellipse(type = "norm", level = 0.95) +
      scale_color_manual(values = batch_colors) +
      facet_wrap(~ Stage) +
      labs(title = "PCA Before and After ComBat Batch Correction",
           x = paste0("PC1 (Before: ", dat$var_before[1], "%, After: ", dat$var_after[1], "%)"),
           y = paste0("PC2 (Before: ", dat$var_before[2], "%, After: ", dat$var_after[2], "%)")) +
      theme_bw() + theme(legend.position = "bottom")
    ggsave(file, plot = plot, width = 12, height = 6, dpi = 150)
  }
)