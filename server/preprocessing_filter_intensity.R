# server/preprocessing_filter_intensity.R
message("[DEBUG] preprocessing_filter_intensity.R loaded")

# ---------- 强度过滤效果 ----------
output$intensity_filter_effect <- renderPrint({
  message("[DEBUG] intensity_filter_effect called")
  req(expression_data(), input$max_missing_fraction, input$min_intensity)
  data <- expression_data()
  
  # 与缺失值过滤使用相同的模式
  mode <- input$missing_filter_mode
  threshold_missing <- input$max_missing_fraction
  data <- apply_missing_filter(data, threshold_missing, mode, rv$sample_info, rv$sample_names)
  
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
  mode <- input$missing_filter_mode
  data <- apply_missing_filter(data, input$max_missing_fraction, mode, rv$sample_info, rv$sample_names)
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

# ---------- 强度统计（使用当前缺失值过滤模式，与calc_result同步） ----------
intensity_stats <- reactive({
  message("[DEBUG] intensity_stats called, using expression_data from data_upload.R")
  req(expression_data(), input$max_missing_fraction)
  data <- expression_data()
  mode <- input$missing_filter_mode
  threshold_missing <- input$max_missing_fraction
  
  # 应用缺失值过滤
  data <- apply_missing_filter(data, threshold_missing, mode, rv$sample_info, rv$sample_names)
  # 去除Inf
  max_intensities <- apply(data, 1, max, na.rm = TRUE)
  keep_finite <- is.finite(max_intensities)
  if (sum(keep_finite) == 0) {
    validate(need(FALSE, "All max intensities are non-finite after filtering."))
  }
  max_intensities <- max_intensities[keep_finite]
  
  if (length(max_intensities) == 0) {
    validate(need(FALSE, "No proteins remain after missing/Inf filter."))
  }
  
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

# ============ 强度过滤输入数据（与强度过滤效果完全同步） ============
intensity_filter_input_data <- reactive({
  message("[DEBUG] intensity_filter_input_data: computing input for intensity filter")
  req(expression_data())
  
  data <- expression_data()
  mode <- input$missing_filter_mode
  threshold_missing <- input$max_missing_fraction
  
  # 1. 缺失值过滤
  filtered <- apply_missing_filter(data, threshold_missing, mode, rv$sample_info, rv$sample_names)
  
  # 2. Inf过滤
  max_int <- apply(filtered, 1, max, na.rm = TRUE)
  keep_finite <- is.finite(max_int)
  result <- filtered[keep_finite, , drop = FALSE]
  
  message("[DEBUG] intensity_filter_input_data: after missing+inf filter, dim = ", nrow(result), " x ", ncol(result))
  
  # 3. 获取蛋白ID（与缺失值过滤导出逻辑相同）
  row_names <- rownames(data)
  # 获取原始蛋白ID（所有蛋白）
  if (!is.null(rv$clean_data) && "Master protein IDs" %in% colnames(rv$clean_data)) {
    original_ids <- rv$clean_data$`Master protein IDs`
    if (length(original_ids) != nrow(data)) {
      # 若长度不匹配，可能数据已被处理，使用行名
      original_ids <- row_names
    }
  } else {
    original_ids <- row_names
  }
  
  # 4. 获取过滤后的蛋白ID：需要通过缺失值过滤和Inf过滤后的索引
  # 先得到缺失值过滤后的蛋白索引（基于original_ids）
  missing_keep <- which(rowMeans(is.na(data)) <= threshold_missing)
  # 但上面使用的是 apply_missing_filter，逻辑更复杂，所以我们利用 filtered 的行名来匹配
  filtered_row_names <- rownames(filtered)
  if (suppressWarnings(all(!is.na(as.numeric(filtered_row_names))))) {
    # 行名是数字，表示原始行名也是数字，直接使用数字作为索引
    filtered_indices <- as.integer(filtered_row_names)
    filtered_ids <- original_ids[filtered_indices]
  } else {
    # 行名就是ID
    filtered_ids <- filtered_row_names
  }
  
  # 再去除Inf，同步ID
  result_ids <- filtered_ids[keep_finite]
  
  message("[DEBUG] intensity_filter_input_data: first 5 IDs = ", paste(head(result_ids, 5), collapse = ", "))
  
  list(data = result, ids = result_ids)
})

# ============ 导出强度过滤结果 Excel ============
output$download_intensity_filter_excel <- downloadHandler(
  filename = function() {
    paste0("Intensity_Filter_Result_", Sys.Date(), ".xlsx")
  },
  content = function(file) {
    message("[DEBUG] download_intensity_filter_excel: starting export")
    req(intensity_filter_input_data())
    
    input_data <- intensity_filter_input_data()
    mat <- input_data$data
    protein_ids <- input_data$ids
    threshold <- input$min_intensity
    min_samples <- input$min_samples_above_intensity %||% 1
    
    message("[DEBUG] download_intensity_filter_excel: threshold = ", threshold, ", min_samples = ", min_samples)
    message("[DEBUG] download_intensity_filter_excel: total proteins = ", nrow(mat))
    
    # 计算过滤（与intensity_filter_effect一致）
    above_counts <- apply(mat, 1, function(x) sum(x > threshold, na.rm = TRUE))
    keep <- above_counts >= min_samples
    retained_ids <- protein_ids[keep]
    filtered_ids <- protein_ids[!keep]
    
    message("[DEBUG] download_intensity_filter_excel: retained = ", length(retained_ids), ", filtered = ", length(filtered_ids))
    
    # 创建 workbook
    wb <- openxlsx::createWorkbook()
    
    # ---- Info Sheet ----
    openxlsx::addWorksheet(wb, "Info")
    info_df <- data.frame(
      Info = c(
        "Intensity Filter Export",
        paste("Data source: After Missing Value Filter (mode:", input$missing_filter_mode, 
              ", threshold:", input$max_missing_fraction, ") and Inf/Non-finite Filter"),
        paste("Intensity Threshold:", threshold),
        paste("Minimum samples above threshold:", min_samples),
        paste("Total proteins after missing/Inf filter:", nrow(mat)),
        paste("Proteins retained (>= ", min_samples, " samples > ", threshold, "): ", length(retained_ids), sep = ""),
        paste("Proteins removed:", length(filtered_ids)),
        "",
        "Sheets:",
        "  Retained: proteins that passed intensity filter",
        "  Filtered_Out: proteins that failed intensity filter",
        "  Threshold_Details: for each protein, TRUE/FALSE if intensity > threshold in each sample, and final decision"
      )
    )
    openxlsx::writeData(wb, "Info", info_df)
    
    # ---- Retained Sheet ----
    retained_mat <- mat[keep, , drop = FALSE]
    retained_df <- cbind(ProteinID = retained_ids, retained_mat, stringsAsFactors = FALSE)
    openxlsx::addWorksheet(wb, "Retained")
    openxlsx::writeData(wb, "Retained", retained_df)
    message("[DEBUG] wrote Retained sheet")
    
    # ---- Filtered_Out Sheet ----
    filtered_mat <- mat[!keep, , drop = FALSE]
    filtered_df <- cbind(ProteinID = filtered_ids, filtered_mat, stringsAsFactors = FALSE)
    openxlsx::addWorksheet(wb, "Filtered_Out")
    openxlsx::writeData(wb, "Filtered_Out", filtered_df)
    message("[DEBUG] wrote Filtered_Out sheet")
    
    # ---- Threshold_Details Sheet ----
    above_threshold <- as.data.frame(mat > threshold)
    above_threshold[is.na(above_threshold)] <- FALSE
    above_threshold <- cbind(ProteinID = protein_ids, above_threshold, stringsAsFactors = FALSE)
    above_threshold$Samples_Above_Count <- above_counts
    above_threshold$Retained <- ifelse(keep, "Yes", "No")
    
    openxlsx::addWorksheet(wb, "Threshold_Details")
    openxlsx::writeData(wb, "Threshold_Details", above_threshold)
    message("[DEBUG] wrote Threshold_Details sheet")
    
    openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    message("[DEBUG] download_intensity_filter_excel: export completed")
  }
)