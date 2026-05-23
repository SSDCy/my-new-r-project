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

# ---------- 缺失值填充函数（含 PPCA 回退） ----------
impute_missing_values <- function(data, method = "knn", min_value = 1e-4) {
  if (method == "none") {
    return(data)
  }
  data_matrix <- as.matrix(data)
  if (method == "knn") {
    if (!requireNamespace("impute", quietly = TRUE))
      stop("impute package required for KNN. Run BiocManager::install('impute')")
    suppressMessages({
      impute_result <- impute::impute.knn(data_matrix, k = 10)
    })
    data_matrix <- impute_result$data
  } else if (method == "ppca") {
    if (!requireNamespace("pcaMethods", quietly = TRUE))
      stop("pcaMethods package required for PPCA. Run BiocManager::install('pcaMethods')")
    tryCatch({
      pc <- pcaMethods::ppca(data_matrix, nPcs = 2, scale = "uv", center = TRUE)
      data_matrix <- as.matrix(pcaMethods::completeObs(pc))
    }, error = function(e) {
      showNotification(
        paste("PPCA imputation failed:", e$message, "- Automatically switching to KNN imputation."),
        type = "warning", duration = 10
      )
      if (!requireNamespace("impute", quietly = TRUE))
        stop("impute package required for fallback KNN. Run BiocManager::install('impute')")
      suppressMessages({
        impute_result <- impute::impute.knn(data_matrix, k = 10)
      })
      data_matrix <<- impute_result$data
    })
  } else {
    stop("Unknown imputation method.")
  }
  rownames(data_matrix) <- rownames(data)
  colnames(data_matrix) <- colnames(data)
  result <- as.data.frame(data_matrix)
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
    
    # 1. 缺失值过滤
    if (input$max_missing_fraction < 1) {
      missing_frac <- rowMeans(is.na(data))
      data <- data[missing_frac <= input$max_missing_fraction, , drop = FALSE]
      if (nrow(data) == 0) stop("No proteins left after missing value filter. Relax the threshold.")
    }
    
    # 2. 强度过滤
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
    
    # 3. 缺失值填补
    preprocessing_params$imputation_method <- input$imputation_method
    data <- impute_missing_values(data, method = input$imputation_method)
    
    # 4. 批次校正
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
    method_display <- switch(preprocessing_params$imputation_method,
                             knn = "K近邻填补法",
                             ppca = "概率主成分分析填补法",
                             none = "无（跳过填充）")
    cat("缺失值处理方式：", method_display, "\n")
    cat("最后运行时间：", format(preprocessing_params$last_run_time, "%Y-%m-%d %H:%M:%S"), "\n")
    
    cat("\n预处理步骤顺序：\n")
    cat("1. 缺失值过滤 (max fraction = ", input$max_missing_fraction, ")\n", sep = "")
    cat("2. 异常值(Inf)过滤（始终执行）\n")
    cat("3. 强度过滤 (min intensity = ", input$min_intensity, ")\n", sep = "")
    cat("4. 缺失值填补 (", method_display, ")\n", sep = "")
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

# ============ 缺失值过滤前后对比 ============

output$preprocessing_done <- reactive({
  !is.null(processed_data())
})
outputOptions(output, "preprocessing_done", suspendWhenHidden = FALSE)

filter_comparison_data <- reactive({
  req(processed_data())
  before <- expression_data()
  after <- processed_data()
  
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
  pass_missing <- before_missing_rate <= input$max_missing_fraction
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
    Missing_Fraction_Threshold = input$max_missing_fraction,
    Intensity_Threshold = input$min_intensity,
    stringsAsFactors = FALSE
  )
  
  raw_for_export <- cbind(Protein_ID = protein_ids, before)
  
  before_impute <- impute_missing_values(before, method = "knn")
  after_impute <- after
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
  req(filter_comparison_data())
  dat <- filter_comparison_data()
  
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
  tryCatch(
    prcomp(mat, scale. = scale),
    error = function(e) NULL
  )
}

output$filter_pca_plot <- renderPlot({
  req(filter_comparison_data())
  dat <- filter_comparison_data()
  
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
  req(filter_comparison_data())
  detailed <- filter_comparison_data()$detailed
  DT::datatable(
    detailed,
    options = list(pageLength = 10, scrollX = TRUE),
    rownames = FALSE
  )
})

output$download_filter_table <- downloadHandler(
  filename = function() {
    paste0("Filter_Comparison_", Sys.Date(), ".xlsx")
  },
  content = function(file) {
    dat <- filter_comparison_data()
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
      Metric = c(
        "Total Proteins (Before Filtering)",
        "Inf Value Filter Removed",
        "Missing Rate Filter Removed",
        "Intensity Filter Removed",
        "Final Retained Proteins"
      ),
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
      "Experiment Name",
      "Number of Samples",
      "Analysis Time",
      "Analyst",
      "",
      "--- Filtering Procedure (applied in order) ---",
      "Step 1: Inf/Abnormal Value Filter",
      "  Rule",
      "  Filtered Protein Count",
      "Step 2: Missing Rate Filter",
      "  Formula (Missing Rate = missing samples / total samples)",
      "  Threshold (max allowed fraction)",
      "  Filtered Protein Count",
      "Step 3: Intensity Filter",
      "  Formula (Intensity = max expression value across all samples)",
      "  Threshold (min intensity)",
      "  Filtered Protein Count",
      "",
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
      uploaded_name,
      as.character(n_samples),
      format(now_time, "%Y-%m-%d %H:%M:%S"),
      "Not provided",
      "",
      "",
      "",
      "Proteins with non-finite (Inf/NaN) maximum intensity are removed",
      as.character(inf_removed),
      "",
      "",
      as.character(miss_threshold),
      as.character(missing_removed),
      "",
      "",
      as.character(inten_threshold),
      as.character(intensity_removed),
      "",
      "",
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
  if (input$max_missing_fraction < 1) {
    missing_frac <- rowMeans(is.na(data))
    data <- data[missing_frac <= input$max_missing_fraction, , drop = FALSE]
  }
  max_int <- apply(data, 1, max, na.rm = TRUE)
  keep_finite <- is.finite(max_int)
  data <- data[keep_finite, , drop = FALSE]
  if (nrow(data) > 0 && input$min_intensity > 0) {
    max_int_finite <- apply(data, 1, max, na.rm = TRUE)
    keep <- max_int_finite > input$min_intensity
    data <- data[keep, , drop = FALSE]
  }
  data
})

imputation_comparison_data <- reactive({
  req(processed_data(), pre_imputation_matrix())
  before_imp <- pre_imputation_matrix()
  after_imp <- processed_data()
  
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
  params <- if (method == "knn") "k = 10 (default), using 10 nearest proteins for imputation" else if (method == "ppca") "nPcs = 2, scale = 'uv', center = TRUE" else "None"
  
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
  req(imputation_comparison_data())
  dat <- imputation_comparison_data()
  cat("Total missing values before imputation:", dat$total_missing_before, "\n")
  cat("Missing values after imputation:", dat$missing_after, "\n")
  cat("Imputation method:", dat$method, "\n")
  cat("Parameters:", dat$params, "\n")
})

output$imputation_boxplot <- renderPlot({
  req(imputation_comparison_data())
  dat <- imputation_comparison_data()
  
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
  method_display <- switch(dat$method, knn = "KNN", ppca = "PPCA", none = "None")
  subtitle_text <- paste0(method_display, " Imputation | ", n_proteins, " proteins | ", n_samples, " samples")
  
  group_summary <- data.frame(
    Sample = sample_names,
    Group = group_levels,
    stringsAsFactors = FALSE
  )
  
  max_y <- max(plot_df$Intensity, na.rm = TRUE)
  y_tile <- max_y * 1.08
  tile_height <- max_y * 0.02
  
  group_colors <- c("L2" = "#E41A1C", "L4" = "#377EB8", "L6" = "#4DAF4A", "L8" = "#984EA3")
  
  ggplot(plot_df, aes(x = Sample, y = Intensity, fill = Stage)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
    geom_tile(data = group_summary, aes(x = Sample, y = y_tile, fill = Group),
              width = 0.9, height = tile_height, inherit.aes = FALSE) +
    scale_fill_manual(
      values = c("Before Imputation" = "#e67e22", "After Imputation" = "#9b59b6",
                 setNames(group_colors, names(group_colors))),
      breaks = c("Before Imputation", "After Imputation"),
      name = "Stage"
    ) +
    labs(title = "Intensity Distribution (log2) Before and After Imputation",
         subtitle = subtitle_text,
         y = "log2(Intensity)", x = "") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          legend.position = "bottom")
})

output$imputation_pca_plot <- renderPlot({
  req(imputation_comparison_data())
  dat <- imputation_comparison_data()
  
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
  req(imputation_comparison_data())
  dat <- imputation_comparison_data()
  
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
  req(imputation_comparison_data())
  DT::datatable(
    imputation_comparison_data()$detailed,
    options = list(pageLength = 10, scrollX = TRUE),
    rownames = FALSE
  )
})

output$download_imputation_table <- downloadHandler(
  filename = function() {
    paste0("Imputation_Comparison_", Sys.Date(), ".xlsx")
  },
  content = function(file) {
    dat <- imputation_comparison_data()
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