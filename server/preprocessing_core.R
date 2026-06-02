# server/preprocessing_core.R
message("[DEBUG] preprocessing_core.R loaded")

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
    
    # 1. 缺失值过滤
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
      preprocessing_params$missing_filter_mode <- "global"
      message("[DEBUG] processed_data: missing filter fell back to global")
      showNotification(
        paste0("⚠️ 分组缺失值过滤：", unmatched_count, " 个样本未匹配到任何组，已自动切换为全局模式过滤。"),
        type = "warning", duration = 10, id = "missing_fallback_notif"
      )
      updateSelectInput(session, "missing_filter_mode", selected = "global")
    }
    
    if (nrow(data) == 0) stop("No proteins left after missing value filter. Relax the threshold.")
    
    # 2. 强度过滤（包含 Inf 过滤）
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
      if (nrow(data) == 0) stop("No proteins left after intensity filter.")
    }
    
    # 3. 缺失值填补
    knn_k <- input$knn_k
    if (is.null(knn_k) || is.na(knn_k) || knn_k < 1) knn_k <- 10
    preprocessing_params$knn_k <- knn_k
    
    if (input$imputation_method == "none") {
      preprocessing_params$imputation_method <- "none"
    } else {
      data <- impute_missing_values(data, method = input$imputation_method, k = knn_k)
      actual_method <- attr(data, "actual_method")
      if (is.null(actual_method)) actual_method <- input$imputation_method
      preprocessing_params$imputation_method <- actual_method
    }
    
    # 4. 批次校正
    preprocessing_params$batch_performed <- FALSE
    preprocessing_params$batch_corrected_cols <- NULL
    preprocessing_params$batch_uncorrected_cols <- NULL
    preprocessing_params$batch_match_summary <- NULL
    
    if (isTRUE(input$perform_batch_correction)) {
      if (!is.null(rv$sample_info) && "Batch" %in% colnames(rv$sample_info)) {
        expr_cols <- colnames(data)
        info_rows <- rownames(rv$sample_info)
        expr_norm <- standardize_sample_name(expr_cols)
        info_norm <- standardize_sample_name(info_rows)
        match_idx <- match(expr_norm, info_norm)
        n_total <- length(expr_cols)
        n_matched <- sum(!is.na(match_idx))
        n_unmatched <- n_total - n_matched
        matched_samples <- expr_cols[!is.na(match_idx)]
        unmatched_samples <- if (n_unmatched > 0) expr_cols[is.na(match_idx)] else character(0)
        
        preprocessing_params$batch_match_summary <- list(total = n_total, matched = n_matched, unmatched = n_unmatched,
                                                         matched_samples = matched_samples, unmatched_samples = unmatched_samples)
        if (n_matched > 0) {
          keep_samples <- !is.na(match_idx)
          batch_info <- rv$sample_info[match_idx[keep_samples], "Batch"]
          if (length(unique(batch_info)) >= 2) {
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
          }
        }
      }
    }
    
    preprocessing_params$last_run_time <- Sys.time()
    removeNotification("preprocess_notif")
    showNotification("Preprocessing completed!", type = "message", duration = 3)
    updateNavbarPage(session, "main_navbar", selected = "plots")
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

# ============ 原始 pre_raw_summary（保留，但不在 UI 中显示） ============
output$pre_raw_summary <- renderPrint({
  tryCatch({
    req(expression_data())
    stats <- intensity_stats()
    detail <- missing_filter_prediction_detail()
    
    message("[DEBUG] pre_raw_summary: mode = ", detail$mode)
    
    mode <- input$missing_filter_mode
    threshold <- input$max_missing_fraction
    data <- expression_data()
    filtered_data <- apply_missing_filter(data, threshold, mode, rv$sample_info, rv$sample_names)
    max_int <- apply(filtered_data, 1, max, na.rm = TRUE)
    inf_removed <- sum(!is.finite(max_int))
    after_inf <- sum(is.finite(max_int))
    inf_ids <- rownames(filtered_data)[!is.finite(max_int)]
    message("[DEBUG] pre_raw_summary: inf removed = ", inf_removed, ", after missing+inf = ", after_inf)
    message("[DEBUG] pre_raw_summary: first few inf IDs = ", paste(head(inf_ids, 5), collapse = ", "))
    
    cat("Raw data dimensions:", nrow(expression_data()), "proteins,", ncol(expression_data()), "samples\n")
    cat("Total missing values:", sum(is.na(expression_data())), "\n")
    cat("Overall missing ratio:", round(sum(is.na(expression_data()))/(nrow(expression_data())*ncol(expression_data()))*100, 2), "%\n")
    cat("Proteins with any missing:", sum(rowSums(is.na(expression_data())) > 0), "\n\n")
    cat("====================\n")
    cat("Current missing filter setting: max allowed missing =", input$max_missing_fraction, "\n")
    cat("Filter mode:", detail$mode, "\n")
    cat("Predicted removal by missing filter:", detail$removed, "proteins\n")
    
    if (inf_removed > 0) {
      cat("\n--- Inf/Non-finite Value Filter ---\n")
      cat("Proteins removed due to non-finite (Inf/NaN) max intensity:", inf_removed, "\n")
      cat("List of removed protein IDs (first 30):\n")
      if (length(inf_ids) <= 30) {
        cat(paste(inf_ids, collapse = "\n"), "\n")
      } else {
        cat(paste(head(inf_ids, 30), collapse = "\n"))
        cat("\n... (total ", length(inf_ids), " proteins)\n", sep = "")
      }
      cat("----------------------------------\n")
    } else {
      cat("Inf/Non-finite filter: 0 proteins removed.\n")
    }
    cat("After missing/Inf filter:", after_inf, "proteins\n\n")
    
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
      threshold_val <- format(stats$thresholds[i], scientific = FALSE, big.mark = ",")
      filtered <- stats$filtered_counts[i]
      filtered_percent <- round(filtered / stats$total_proteins * 100, 1)
      retained <- stats$total_proteins - filtered
      cat(q_percent, "% quantile =", threshold_val,
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

# ========== 修复后的 Processed Data Table（自动修正行名） ==========
output$pre_processed_table <- DT::renderDT({
  message("[DEBUG] output$pre_processed_table called")
  req(processed_data())
  df <- processed_data()
  
  # 获取蛋白ID，如果行名是数字序号，则从原始数据中映射真实ID
  ids <- rownames(df)
  message("[DEBUG] pre_processed_table: first 5 rownames = ", paste(head(ids, 5), collapse = ", "))
  
  if (suppressWarnings(all(!is.na(as.numeric(ids))))) {
    message("[DEBUG] pre_processed_table: rownames are numeric, mapping to Master protein IDs")
    # 尝试从 rv$clean_data 获取原始ID列表
    if (!is.null(rv$clean_data) && "Master protein IDs" %in% colnames(rv$clean_data)) {
      original_ids <- rv$clean_data$`Master protein IDs`
      # 检查长度是否与 expression_data() 一致（未过滤前）
      # 实际上 processed_data() 的行是过滤后的，行名数字可能是相对于过滤后的矩阵的索引
      # 但更稳健的方式：直接使用 expression_data() 的原始行名，再通过匹配进行过滤？比较复杂。
      # 简单方案：如果行名全是数字，很可能就是索引，我们就直接用这些索引去 original_ids 取。
      idx <- as.integer(ids)
      if (max(idx, na.rm = TRUE) <= length(original_ids)) {
        ids <- original_ids[idx]
        message("[DEBUG] pre_processed_table: ID mapping done, first 5 IDs = ", paste(head(ids, 5), collapse = ", "))
      } else {
        message("[DEBUG] pre_processed_table: numeric rownames out of range, keeping original")
      }
    } else {
      message("[DEBUG] pre_processed_table: no clean_data or Master protein IDs column, cannot remap")
    }
  }
  
  df <- cbind(`Master Protein ID` = ids, df)
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

output$preprocessing_done <- reactive({ !is.null(processed_data()) })
outputOptions(output, "preprocessing_done", suspendWhenHidden = FALSE)

output$imputation_skipped <- reactive({
  !is.null(processed_data()) && preprocessing_params$imputation_method == "none"
})
outputOptions(output, "imputation_skipped", suspendWhenHidden = FALSE)

# ========== 新增：缺失数据摘要（Missing Value Filter 页卡底部） ==========
output$missing_data_info <- renderPrint({
  message("[DEBUG] output$missing_data_info called")
  tryCatch({
    req(expression_data())
    detail <- missing_filter_prediction_detail()
    data <- expression_data()
    mode <- input$missing_filter_mode
    threshold <- input$max_missing_fraction
    filtered_data <- apply_missing_filter(data, threshold, mode, rv$sample_info, rv$sample_names)
    max_int <- apply(filtered_data, 1, max, na.rm = TRUE)
    inf_removed <- sum(!is.finite(max_int))
    after_inf <- sum(is.finite(max_int))
    inf_ids <- rownames(filtered_data)[!is.finite(max_int)]
    
    cat("Raw data dimensions:", nrow(data), "proteins,", ncol(data), "samples\n")
    cat("Total missing values:", sum(is.na(data)), "\n")
    cat("Overall missing ratio:", round(sum(is.na(data))/(nrow(data)*ncol(data))*100, 2), "%\n")
    cat("Proteins with any missing:", sum(rowSums(is.na(data)) > 0), "\n\n")
    cat("====================\n")
    cat("Current missing filter setting: max allowed missing =", threshold, "\n")
    cat("Filter mode:", detail$mode, "\n")
    cat("Predicted removal by missing filter:", detail$removed, "proteins\n")
    
    if (inf_removed > 0) {
      cat("\n--- Inf/Non-finite Value Filter ---\n")
      cat("Proteins removed due to non-finite (Inf/NaN) max intensity:", inf_removed, "\n")
      cat("List of removed protein IDs (first 30):\n")
      if (length(inf_ids) <= 30) {
        cat(paste(inf_ids, collapse = "\n"), "\n")
      } else {
        cat(paste(head(inf_ids, 30), collapse = "\n"))
        cat("\n... (total ", length(inf_ids), " proteins)\n", sep = "")
      }
      cat("----------------------------------\n")
    } else {
      cat("Inf/Non-finite filter: 0 proteins removed.\n")
    }
    cat("After missing/Inf filter:", after_inf, "proteins\n")
  }, error = function(e) {
    cat("Error generating missing data info:\n", e$message, "\n")
  })
})
message("[DEBUG] output$missing_data_info defined")

# ========== 新增：强度分位数信息（Minimum Intensity Filter 折叠内） ==========
output$intensity_info <- renderPrint({
  message("[DEBUG] output$intensity_info called")
  tryCatch({
    req(expression_data())
    stats <- intensity_stats()
    
    cat("Intensity distribution (per protein max intensity):\n")
    qs <- c(0, 0.25, 0.5, 0.75, 1)
    vals <- quantile(stats$max_intensities, probs = qs, na.rm = TRUE)
    cat("Min:", format(vals[1], scientific = FALSE, big.mark = ","), "\n")
    cat("25%:", format(vals[2], scientific = FALSE, big.mark = ","), "\n")
    cat("Median:", format(vals[3], scientific = FALSE, big.mark = ","), "\n")
    cat("75%:", format(vals[4], scientific = FALSE, big.mark = ","), "\n")
    cat("Max:", format(vals[5], scientific = FALSE, big.mark = ","), "\n\n")
    
    cat("Different quantile thresholds and predicted removal:\n")
    for (i in seq_along(stats$thresholds)) {
      q_percent <- stats$quantiles[i] * 100
      threshold_val <- format(stats$thresholds[i], scientific = FALSE, big.mark = ",")
      filtered <- stats$filtered_counts[i]
      filtered_percent <- round(filtered / stats$total_proteins * 100, 1)
      retained <- stats$total_proteins - filtered
      cat(q_percent, "% quantile =", threshold_val,
          ": remove", filtered, "proteins (", filtered_percent, "%), keep", retained, "\n")
    }
    cat("\nRecommended minimum intensity threshold:", format(round(stats$recommended_threshold), scientific = FALSE, big.mark = ","), "\n")
    cat("(This threshold will filter out approximately 15% of low-quality proteins)\n")
  }, error = function(e) {
    cat("Error generating intensity info:\n", e$message, "\n")
  })
})
message("[DEBUG] output$intensity_info defined")

# ========== 新增：预处理步骤摘要（Processed Data Table 页卡） ==========
output$preprocessing_steps_summary <- renderPrint({
  message("[DEBUG] output$preprocessing_steps_summary called")
  if (is.null(processed_data())) {
    cat("Preprocessing has not been run yet. Please click 'Run Preprocessing' to see the processed data.\n")
    return()
  }
  cat("Preprocessing performed at:", format(preprocessing_params$last_run_time, "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  # 步骤1：缺失值过滤
  cat("1. Missing Value Filter:\n")
  cat("   Mode:", preprocessing_params$missing_filter_mode, "\n")
  cat("   Threshold:", input$max_missing_fraction, "\n")
  if (isTRUE(preprocessing_params$missing_filter_fallback)) {
    cat("   Fallback: global mode (", preprocessing_params$missing_filter_fallback_unmatched, " samples unmatched)\n")
  }
  
  # 步骤2：Inf过滤
  cat("2. Inf/Non-finite Filter:\n")
  cat("   Removed proteins:", preprocessing_params$inf_filtered_count, "\n")
  if (preprocessing_params$inf_filtered_count > 0) {
    cat("   First few removed IDs:", paste(head(preprocessing_params$inf_filtered_proteins, 5), collapse = ", "), "\n")
  }
  
  # 步骤3：强度过滤
  cat("3. Minimum Intensity Filter:\n")
  cat("   Threshold:", input$min_intensity, "\n")
  cat("   Min samples above threshold:", preprocessing_params$intensity_min_samples, "\n")
  
  # 步骤4：填充
  cat("4. Missing Value Imputation:\n")
  if (!is.null(preprocessing_params$imputation_method)) {
    method_display <- switch(preprocessing_params$imputation_method,
                             knn = "K-Nearest Neighbors",
                             ppca = "Probabilistic PCA",
                             minvalue = "Fixed Minimum Value",
                             quantile = "Quantile",
                             none = "None (skipped)")
    cat("   Method:", method_display, "\n")
    if (grepl("knn", preprocessing_params$imputation_method)) {
      cat("   k =", preprocessing_params$knn_k, "\n")
    } else if (preprocessing_params$imputation_method == "ppca") {
      cat("   nPcs = 2 (via pcaMethods::pca with log2 transform)\n")
    } else if (preprocessing_params$imputation_method == "minvalue") {
      cat("   Fixed value =", preprocessing_params$min_value, "\n")
    } else if (preprocessing_params$imputation_method == "quantile") {
      cat("   Quantile =", preprocessing_params$quantile_prob, "\n")
    }
  } else {
    cat("   No imputation performed yet.\n")
  }
  
  # 步骤5：批次校正
  cat("5. Batch Correction:\n")
  if (preprocessing_params$batch_performed) {
    cat("   ComBat applied to", length(preprocessing_params$batch_corrected_cols), "samples\n")
    if (!is.null(preprocessing_params$batch_uncorrected_cols) && length(preprocessing_params$batch_uncorrected_cols) > 0) {
      cat("   Uncorted samples:", paste(preprocessing_params$batch_uncorrected_cols, collapse = ", "), "\n")
    }
  } else {
    cat("   Not performed.\n")
  }
  
  cat("\nFinal data dimensions:", nrow(processed_data()), "proteins,", ncol(processed_data()), "samples\n")
  cat("Remaining missing values:", sum(is.na(processed_data())), "\n")
})
message("[DEBUG] output$preprocessing_steps_summary defined")