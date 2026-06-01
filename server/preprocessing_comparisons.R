# server/preprocessing_comparisons.R
message("[DEBUG] preprocessing_comparisons.R loaded")

filter_comparison_data <- reactive({
  req(processed_data())
  before <- expression_data()
  after <- processed_data()
  
  message("[DEBUG] filter_comparison_data: raw before dim = ", nrow(before), "x", ncol(before), 
          "; after dim = ", nrow(after), "x", ncol(after))
  
  filter_mode <- preprocessing_params$missing_filter_mode
  threshold <- input$max_missing_fraction
  message("[DEBUG] filter_comparison_data: mode = ", filter_mode, ", threshold = ", threshold)
  
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
  
  before_impute <- impute_missing_values(before, method = "knn")
  after_impute <- after
  
  message("[DEBUG] filter_comparison_data: before_impute dim = ", nrow(before_impute), "x", ncol(before_impute),
          "; after_impute dim = ", nrow(after_impute), "x", ncol(after_impute))
  
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