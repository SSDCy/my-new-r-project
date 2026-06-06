# server/preprocessing_filter_intensity.R
message("[DEBUG] preprocessing_filter_intensity.R loaded - intensity filter disabled, fixed group mode")

# 强度过滤效果说明（显示已禁用）
output$intensity_filter_effect <- renderPrint({
  cat("Minimum intensity filter has been disabled.\n")
})

# 强度分布图（不再调用强度过滤）
output$intensity_dist_plot <- renderPlot({
  req(expression_data())
  threshold <- 0   # 不再使用
  data <- expression_data()
  # 直接使用分组缺失值过滤
  mode <- "group"
  data <- apply_missing_filter(data, input$max_missing_fraction, mode, rv$sample_info, rv$sample_names)
  max_int <- apply(data, 1, max, na.rm = TRUE)
  max_int <- max_int[is.finite(max_int)]
  if (length(max_int) == 0) return(NULL)
  
  log_int <- log10(max_int + 1)
  df <- data.frame(log_int = log_int)
  
  total_count <- length(max_int)
  
  p <- ggplot(df, aes(x = log_int)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "steelblue", alpha = 0.6) +
    geom_density(color = "darkorange", linewidth = 1.2) +
    labs(title = "Protein Max Intensity Distribution (After Missing Filter)",
         subtitle = paste0("After missing & Inf filter: ", total_count, " proteins"),
         x = "log10(Max Intensity + 1)", y = "Density") +
    theme_bw() +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9),
      plot.margin = margin(t = 10, r = 10, b = 5, l = 10)
    )
  p
})

# 强度统计（用于缺失值过滤后的强度分布参考）
intensity_stats <- reactive({
  message("[DEBUG] intensity_stats called, using expression_data and group missing filter")
  req(expression_data(), input$max_missing_fraction)
  data <- expression_data()
  mode <- "group"   # 强制分组
  threshold_missing <- input$max_missing_fraction
  
  data <- apply_missing_filter(data, threshold_missing, mode, rv$sample_info, rv$sample_names)
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

# 阈值计算器（已不再需要，但保留）
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

# 强度过滤导出（已不可用）
output$download_intensity_filter_excel <- downloadHandler(
  filename = function() { "Intensity_Filter_Result_disabled.xlsx" },
  content = function(file) {
    showNotification("Intensity filter has been disabled.", type = "error")
  }
)

message("[DEBUG] preprocessing_filter_intensity.R: fixed to use group mode")