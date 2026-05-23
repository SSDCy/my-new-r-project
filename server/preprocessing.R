# server/preprocessing.R

# ---------- 动态填充批次列选项 ----------
observe({
  req(rv$sample_info)
  updateSelectInput(session, "batch_column", choices = colnames(rv$sample_info))
})

# ---------- 显示样本信息各列的概要 ----------
output$batch_info_preview <- renderPrint({
  req(rv$sample_info, input$batch_column)
  col_data <- rv$sample_info[[input$batch_column]]
  unique_vals <- unique(col_data)
  cat("样本信息表“", input$batch_column, "”列概览：\n", sep = "")
  cat("样本总数：", length(col_data), "\n")
  cat("唯一值数量：", length(unique_vals), "\n")
  cat("唯一值列表（前10个）：\n")
  if (length(unique_vals) <= 10) {
    cat(paste(unique_vals, collapse = "\n"), "\n")
  } else {
    cat(paste(head(unique_vals, 10), collapse = "\n"), "\n... 还有", length(unique_vals) - 10, "个\n")
  }
})

# ---------- 提取原始表达矩阵（仅 LFQ 列），使用 Master protein IDs 作为行名 ----------
expression_data <- reactive({
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
  df
})

# ---------- 动态缺失值过滤效果 ----------
output$missing_filter_effect <- renderPrint({
  req(expression_data())
  missing_per_protein <- rowMeans(is.na(expression_data()))
  filtered <- sum(missing_per_protein > input$max_missing_fraction)
  retained <- nrow(expression_data()) - filtered
  filtered_percent <- round(filtered / nrow(expression_data()) * 100, 1)
  cat("Predicted removal:", filtered, "proteins (", filtered_percent, "%)\n", sep = "")
  cat("Predicted retained:", retained, "proteins\n", sep = "")
})

# ---------- 强度过滤效果 ----------
output$intensity_filter_effect <- renderPrint({
  req(expression_data(), input$max_missing_fraction, input$min_intensity)
  data <- expression_data()
  missing_frac <- rowMeans(is.na(data))
  keep <- missing_frac <= input$max_missing_fraction
  data <- data[keep, , drop = FALSE]
  max_int <- apply(data, 1, max, na.rm = TRUE)
  keep_finite <- is.finite(max_int)
  data <- data[keep_finite, , drop = FALSE]
  max_int <- max_int[keep_finite]
  below_thresh <- max_int < input$min_intensity
  cat("After missing and Inf filter, ", sum(below_thresh), " proteins are below the intensity threshold and will be removed.\n")
})

# ---------- 强度统计（用于阈值推荐和原始数据概览） ----------
intensity_stats <- reactive({
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

# ---------- 缺失值填充函数 ----------
impute_missing_values <- function(data, method = "min", min_value = 1e-4) {
  if (method == "none") {
    return(data)
  }
  data_matrix <- as.matrix(data)
  if (method == "min") {
    data_matrix[is.na(data_matrix)] <- min_value
  } else if (method == "mean") {
    data_matrix <- apply(data_matrix, 2, function(x) {
      x[is.na(x)] <- mean(x, na.rm = TRUE); x
    })
  } else if (method == "median") {
    data_matrix <- apply(data_matrix, 2, function(x) {
      x[is.na(x)] <- median(x, na.rm = TRUE); x
    })
  } else if (method == "knn") {
    if (!requireNamespace("impute", quietly = TRUE))
      stop("impute package required for KNN. Run BiocManager::install('impute')")
    impute_result <- impute::impute.knn(data_matrix)
    data_matrix <- impute_result$data
  }
  result <- as.data.frame(data_matrix)
  rownames(result) <- rownames(data)
  return(result)
}

# ---------- ComBat 批次校正 ----------
combat_correction <- function(data, batch) {
  if (!requireNamespace("sva", quietly = TRUE))
    stop("sva package required. Run BiocManager::install('sva')")
  batch_factor <- as.factor(batch)
  if (length(levels(batch_factor)) < 2) {
    stop("批次校正需要至少两个不同的批次值。")
  }
  corrected <- sva::ComBat(dat = as.matrix(data), batch = batch_factor)
  result <- as.data.frame(corrected)
  rownames(result) <- rownames(data)
  return(result)
}

# ---------- 样本名标准化（用于匹配） ----------
normalize_sample_name <- function(x) {
  x <- gsub("[\\. _]+", "-", x)
  x <- gsub("-+", "-", x)
  x <- sub("^(LFQ-intensity-|Intensity-)", "", x)
  x <- sub("^(LFQ-intensity|Intensity)-", "", x)
  return(x)
}

# ---------- 预处理参数记录 ----------
preprocessing_params <- reactiveValues(
  imputation_method = NULL,
  last_run_time = NULL,
  inf_filtered_count = 0,
  inf_filtered_proteins = character(0),
  batch_performed = FALSE,
  batch_corrected_cols = NULL
)

# ============ 核心预处理反应式 ============
processed_data <- eventReactive(input$run_preprocessing, {
  showNotification("Running preprocessing...", type = "message", duration = NULL, id = "preprocess_notif")
  tryCatch({
    data <- expression_data()
    
    # 1. 缺失值过滤（如果阈值 < 1）
    if (input$max_missing_fraction < 1) {
      missing_frac <- rowMeans(is.na(data))
      data <- data[missing_frac <= input$max_missing_fraction, , drop = FALSE]
      if (nrow(data) == 0) stop("No proteins left after missing value filter. Relax the threshold.")
    }
    
    # 2. 强度过滤（始终移除 Inf）
    max_int <- apply(data, 1, max, na.rm = TRUE)
    keep_finite <- is.finite(max_int)
    preprocessing_params$inf_filtered_count <- sum(!keep_finite)
    preprocessing_params$inf_filtered_proteins <- rownames(data)[!keep_finite]
    data <- data[keep_finite, , drop = FALSE]
    if (nrow(data) == 0) stop("No proteins left after removing Inf values.")
    
    if (input$min_intensity > 0) {
      max_int_finite <- apply(data, 1, max, na.rm = TRUE)
      keep <- max_int_finite > input$min_intensity
      data <- data[keep, , drop = FALSE]
      if (nrow(data) == 0) stop("No proteins left after intensity filter. Lower the threshold.")
    }
    
    # 3. 缺失值填充
    preprocessing_params$imputation_method <- input$imputation_method
    data <- impute_missing_values(data, method = input$imputation_method, min_value = input$min_impute_value)
    
    # 4. 批次校正（可选）
    preprocessing_params$batch_performed <- FALSE
    preprocessing_params$batch_corrected_cols <- NULL
    if (input$perform_batch_correction && !is.null(input$batch_column) && input$batch_column != "") {
      expr_cols <- colnames(data)
      info_rows <- rownames(rv$sample_info)
      expr_norm <- normalize_sample_name(expr_cols)
      info_norm <- normalize_sample_name(info_rows)
      match_idx <- match(expr_norm, info_norm)
      if (all(is.na(match_idx))) {
        showNotification("批次校正跳过：无法将表达矩阵中的样本名与样本信息表匹配。", type = "warning")
      } else {
        keep_samples <- !is.na(match_idx)
        batch_info <- rv$sample_info[match_idx[keep_samples], input$batch_column]
        unique_vals <- unique(batch_info)
        if (length(unique_vals) < 2) {
          showNotification(
            paste0("批次校正被跳过：所选批次列“", input$batch_column, "”中只有一种批次值（", 
                   paste(unique_vals, collapse = ", "), "）。"),
            type = "warning", duration = 8
          )
        } else {
          data_to_correct <- data[, keep_samples, drop = FALSE]
          corrected_part <- combat_correction(data_to_correct, batch = batch_info)
          data[, keep_samples] <- corrected_part
          preprocessing_params$batch_performed <- TRUE
          preprocessing_params$batch_corrected_cols <- expr_cols[keep_samples]
        }
      }
    }
    
    preprocessing_params$last_run_time <- Sys.time()
    removeNotification("preprocess_notif")
    showNotification("Preprocessing completed! Redirecting to Analysis & Export...", type = "message", duration = 3)
    
    # 自动跳转到 Plots 页面（分析页面）
    updateNavbarPage(session, "main_navbar", selected = "plots")
    
    return(data)
  }, error = function(e) {
    removeNotification("preprocess_notif")
    showNotification(paste("Preprocessing failed:", e$message), type = "error", duration = 10)
    return(NULL)
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
    method_name <- switch(preprocessing_params$imputation_method,
                          min = "极小值填补法",
                          mean = "均值填补法",
                          median = "中位数填补法",
                          knn = "KNN填补法",
                          none = "无（跳过填充）")
    cat("缺失值处理方式：", method_name, "\n")
    cat("最后运行时间：", format(preprocessing_params$last_run_time, "%Y-%m-%d %H:%M:%S"), "\n")
    
    cat("\n预处理步骤顺序：\n")
    cat("1. 缺失值过滤 (max fraction = ", input$max_missing_fraction, ")\n", sep = "")
    cat("2. 异常值(Inf)过滤（始终执行）\n")
    cat("3. 强度过滤 (min intensity = ", input$min_intensity, ")\n", sep = "")
    cat("4. 缺失值填补 (", method_name, ")\n", sep = "")
    if (preprocessing_params$batch_performed) {
      cat("5. 批次校正\n")
      if (!is.null(preprocessing_params$batch_corrected_cols)) {
        cat("   已校正的样本列数量：", length(preprocessing_params$batch_corrected_cols), "\n")
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
    options = list(
      pageLength = 10,
      scrollX = TRUE,
      searchHighlight = TRUE,
      server = TRUE
    ),
    rownames = FALSE,
    filter = "top"
  )
  
  corrected <- preprocessing_params$batch_corrected_cols
  if (!is.null(corrected) && length(corrected) > 0) {
    for (col in corrected) {
      dt <- DT::formatStyle(dt, columns = col, backgroundColor = "#e6f2ff")
    }
  }
  dt
})