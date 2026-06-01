# server/preprocessing_imputation.R
message("[DEBUG] preprocessing_imputation.R loaded - academic English for export sheets")

# ---------- 填补前数据矩阵 ----------
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

# ---------- 填补前数据（带正确蛋白ID） ----------
pre_imputation_with_ids <- reactive({
  req(pre_imputation_matrix())
  mat <- pre_imputation_matrix()
  
  original_ids <- if (!is.null(rv$clean_data) && "Master protein IDs" %in% colnames(rv$clean_data)) {
    rv$clean_data$`Master protein IDs`
  } else {
    rownames(expression_data())
  }
  
  row_n <- rownames(mat)
  if (suppressWarnings(all(!is.na(as.numeric(row_n))))) {
    idx <- as.integer(row_n)
    message("[DEBUG] pre_imputation_with_ids: rownames are numeric")
  } else {
    message("[DEBUG] pre_imputation_with_ids: rownames are IDs")
    return(list(before = mat, ids = row_n))
  }
  
  protein_ids <- original_ids[idx]
  message("[DEBUG] pre_imputation_with_ids: first 5 IDs = ", paste(head(protein_ids, 5), collapse = ", "))
  
  list(before = mat, ids = protein_ids)
})

# ---------- 全局缺失率过滤后的矩阵（仅用于KNN邻居搜索） ----------
knn_input_data <- reactive({
  if (input$imputation_method != "knn" || is.null(pre_imputation_with_ids())) return(NULL)
  pre <- pre_imputation_with_ids()
  mat <- pre$before
  ids <- pre$ids
  thr <- input$max_missing_fraction
  missing_frac <- rowMeans(is.na(mat))
  keep <- missing_frac <= thr
  mat_filtered <- mat[keep, , drop = FALSE]
  ids_filtered <- ids[keep]
  message("[DEBUG] knn_input_data: after global filter, rows = ", nrow(mat_filtered), " (from ", length(ids), ")")
  list(mat = mat_filtered, ids = ids_filtered)
})

# ---------- 全局邻居信息（基于皮尔逊相关系数） ----------
knn_neighbors_data <- reactive({
  input_data <- knn_input_data()
  if (is.null(input_data)) {
    message("[DEBUG] knn_neighbors_data: no input data")
    return(NULL)
  }
  
  mat <- as.matrix(input_data$mat)
  protein_ids <- input_data$ids
  k_val <- min(input$knn_k, nrow(mat) - 1)
  
  message("[DEBUG] knn_neighbors_data: computing Pearson correlation for ", nrow(mat), " proteins")
  
  cor_mat <- cor(t(mat), use = "pairwise.complete.obs")
  diag(cor_mat) <- -Inf
  
  n_proteins <- nrow(mat)
  neighbors <- matrix(NA_integer_, nrow = n_proteins, ncol = k_val)
  cor_vals <- matrix(NA_real_, nrow = n_proteins, ncol = k_val)
  dist_vals <- matrix(NA_real_, nrow = n_proteins, ncol = k_val)
  
  for (i in 1:n_proteins) {
    ord <- order(cor_mat[i, ], decreasing = TRUE)[1:k_val]
    neighbors[i, ] <- ord
    cor_vals[i, ] <- cor_mat[i, ord]
    dist_vals[i, ] <- 1 - cor_mat[i, ord]
  }
  
  rownames(neighbors) <- protein_ids
  rownames(cor_vals) <- protein_ids
  rownames(dist_vals) <- protein_ids
  colnames(neighbors) <- paste0("N", 1:k_val)
  colnames(cor_vals) <- paste0("Cor", 1:k_val)
  colnames(dist_vals) <- paste0("Dist", 1:k_val)
  
  message("[DEBUG] knn_neighbors_data: neighbors dim = ", nrow(neighbors), " x ", ncol(neighbors))
  message("[DEBUG] knn_neighbors_data: first 5 IDs = ", paste(head(protein_ids, 5), collapse = ", "))
  
  list(neighbors = neighbors, correlations = cor_vals, distances = dist_vals, protein_ids = protein_ids, k = k_val)
})

# ---------- 下拉框更新 ----------
observe({
  data <- knn_neighbors_data()
  if (is.null(data)) {
    message("[DEBUG] updateSelectizeInput: data is NULL")
    return()
  }
  ids <- data$protein_ids
  message("[DEBUG] updateSelectizeInput: updating with ", length(ids), " choices")
  updateSelectizeInput(session, "knn_lookup_protein",
                       choices = ids,
                       server = TRUE,
                       selected = if (length(ids) > 0) ids[1] else NULL)
})

# ---------- 邻居表格输出（UI） ----------
output$knn_lookup_table <- DT::renderDT({
  data <- knn_neighbors_data()
  req(data, input$knn_lookup_protein)
  
  pid <- input$knn_lookup_protein
  idx <- which(data$protein_ids == pid)
  if (length(idx) == 0) {
    message("[DEBUG] knn_lookup_table: protein not found: ", pid)
    return(DT::datatable(data.frame(Message = "Protein not found in neighbor list."), options = list(dom = 't')))
  }
  
  idx <- idx[1]
  n_k <- ncol(data$neighbors)
  neighbors_idx <- data$neighbors[idx, ]
  correlations <- data$correlations[idx, ]
  distances <- data$distances[idx, ]
  
  valid <- !is.na(neighbors_idx)
  if (!any(valid)) {
    return(DT::datatable(data.frame(Message = "No valid neighbors for this protein."), options = list(dom = 't')))
  }
  
  df <- data.frame(
    Rank = seq_len(n_k)[valid],
    NeighborProteinID = data$protein_ids[neighbors_idx[valid]],
    Correlation = round(correlations[valid], 4),
    Distance = round(distances[valid], 4),
    stringsAsFactors = FALSE
  )
  message("[DEBUG] knn_lookup_table: showing ", nrow(df), " neighbors for ", pid)
  DT::datatable(df, options = list(dom = 't', pageLength = 10), rownames = FALSE) |>
    DT::formatStyle("Distance", color = "gray")
})

# ---------- 真正的PPCA填补函数（log2转换 + 反转换） ----------
true_ppca_impute <- function(data, nPcs = 2) {
  message("[DEBUG] true_ppca_impute: starting PPCA with log2 transformation")
  if (!requireNamespace("pcaMethods", quietly = TRUE))
    stop("pcaMethods package required for PPCA. Run BiocManager::install('pcaMethods')")
  
  data_matrix <- as.matrix(data)
  message("[DEBUG] true_ppca_impute: original data range = [", min(data_matrix, na.rm = TRUE), ", ", max(data_matrix, na.rm = TRUE), "]")
  
  log_data <- log2(data_matrix + 1)
  message("[DEBUG] true_ppca_impute: log2 transformed data range = [", min(log_data, na.rm = TRUE), ", ", max(log_data, na.rm = TRUE), "]")
  
  na_rows <- which(rowSums(is.na(log_data)) == ncol(log_data))
  na_cols <- which(colSums(is.na(log_data)) == nrow(log_data))
  if (length(na_rows) > 0 || length(na_cols) > 0) {
    message("[DEBUG] true_ppca_impute: removing all-NA rows/cols")
    if (length(na_rows) > 0) log_data <- log_data[-na_rows, , drop = FALSE]
    if (length(na_cols) > 0) log_data <- log_data[, -na_cols, drop = FALSE]
  }
  
  pc <- pcaMethods::pca(log_data, method = "ppca", nPcs = nPcs, scale = "uv", center = TRUE)
  imputed_log <- as.matrix(pcaMethods::completeObs(pc))
  message("[DEBUG] true_ppca_impute: PPCA on log2 data succeeded")
  
  imputed_original <- 2^imputed_log - 1
  imputed_original[imputed_original < 0] <- 0
  message("[DEBUG] true_ppca_impute: back-transformed range = [", min(imputed_original, na.rm = TRUE), ", ", max(imputed_original, na.rm = TRUE), "]")
  
  return(imputed_original)
}

# ---------- PPCA 可视化数据 ----------
ppca_visualization_data <- reactive({
  message("[DEBUG] ppca_visualization_data: triggered. Method = ", input$imputation_method, 
          ", processed_data exists = ", !is.null(processed_data()))
  if (input$imputation_method != "ppca" || is.null(processed_data())) {
    message("[DEBUG] ppca_visualization_data: not PPCA or no processed_data, returning NULL")
    return(NULL)
  }
  pre <- pre_imputation_with_ids()
  if (is.null(pre)) {
    message("[DEBUG] ppca_visualization_data: pre_imputation_with_ids is NULL")
    return(NULL)
  }
  mat <- as.matrix(pre$before)
  message("[DEBUG] ppca_visualization_data: matrix dim = ", nrow(mat), " x ", ncol(mat))
  
  if (!requireNamespace("pcaMethods", quietly = TRUE)) {
    message("[DEBUG] ppca_visualization_data: pcaMethods not installed")
    return(NULL)
  }
  
  tryCatch({
    log_mat <- log2(mat + 1)
    message("[DEBUG] ppca_visualization_data: starting pcaMethods::pca on log2 matrix")
    pc <- pcaMethods::pca(log_mat, method = "ppca", nPcs = 2, scale = "uv", center = TRUE)
    message("[DEBUG] ppca_visualization_data: PCA complete")
    
    scores <- as.data.frame(pcaMethods::scores(pc))
    scores$Sample <- rownames(scores)
    message("[DEBUG] ppca_visualization_data: scores dim = ", nrow(scores), " x ", ncol(scores))
    
    loadings <- as.data.frame(pcaMethods::loadings(pc))
    loadings$Protein <- rownames(loadings)
    
    imputed_log <- as.matrix(pcaMethods::completeObs(pc))
    imputed_original <- 2^imputed_log - 1
    imputed_original[imputed_original < 0] <- 0
    
    message("[DEBUG] ppca_visualization_data: imputed matrix dim = ", nrow(imputed_original), " x ", ncol(imputed_original))
    
    list(scores = scores, loadings = loadings, imputed = imputed_original, original = mat)
  }, error = function(e) {
    message("[DEBUG] ppca_visualization_data error: ", e$message)
    NULL
  })
})

# ---------- PPCA 得分图 ----------
output$ppca_score_plot <- renderPlot({
  message("[DEBUG] ppca_score_plot: rendering")
  data <- ppca_visualization_data()
  if (is.null(data)) {
    plot.new(); text(0.5, 0.5, "PPCA visualization not available. Please run preprocessing first with PPCA method.")
    message("[DEBUG] ppca_score_plot: data is NULL")
    return()
  }
  
  ggplot(data$scores, aes(x = PC1, y = PC2)) +
    geom_point(size = 3, color = "#3498db") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(title = "PPCA Score Plot (PC1 vs PC2) - on log2 scale",
         subtitle = "Samples plotted in the space of the first two principal components.",
         x = "Principal Component 1", y = "Principal Component 2") +
    annotate("text", x = max(data$scores$PC1)*0.7, y = max(data$scores$PC2)*0.7, 
             label = "Main trend\n→", color = "red", size = 5) +
    theme_bw()
})

# ---------- PPCA 填补值分布图 ----------
output$ppca_imputation_hist <- renderPlot({
  message("[DEBUG] ppca_imputation_hist: rendering")
  data <- ppca_visualization_data()
  if (is.null(data)) {
    plot.new(); text(0.5, 0.5, "Histogram not available.")
    message("[DEBUG] ppca_imputation_hist: data is NULL")
    return()
  }
  
  original_vals <- as.vector(data$original)
  imputed_vals <- as.vector(data$imputed)
  
  na_positions <- is.na(data$original)
  filled_vals <- imputed_vals[na_positions]
  non_missing <- original_vals[!is.na(original_vals)]
  
  df <- data.frame(
    Value = c(non_missing, filled_vals),
    Type = c(rep("Original (non-missing)", length(non_missing)),
             rep("Imputed", length(filled_vals)))
  )
  
  message("[DEBUG] ppca_imputation_hist: non-missing count = ", length(non_missing), 
          ", filled count = ", length(filled_vals))
  message("[DEBUG] ppca_imputation_hist: filled range = [", min(filled_vals), ", ", max(filled_vals), "]")
  message("[DEBUG] ppca_imputation_hist: original non-missing range = [", min(non_missing), ", ", max(non_missing), "]")
  
  ggplot(df, aes(x = Value, fill = Type)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    scale_fill_manual(values = c("Original (non-missing)" = "#3498db", "Imputed" = "#2ecc71")) +
    labs(title = "Distribution of Original vs Imputed Values (log2 PPCA)",
         subtitle = "Imputed values after log2 transformation and back-transformation.",
         x = "Expression Value", y = "Frequency") +
    theme_bw() + theme(legend.position = "bottom")
})

# ---------- 填补比较数据 ----------
imputation_comparison_data <- reactive({
  req(processed_data(), pre_imputation_matrix())
  before_imp <- pre_imputation_matrix()
  after_imp <- processed_data()
  
  if (nrow(before_imp) != nrow(after_imp) || ncol(before_imp) != ncol(after_imp)) {
    message("[DEBUG] imputation_comparison_data: dimension mismatch")
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
    "nPcs = 2, log2 transform applied"
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

# ---------- 填补前后图表 ----------
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

# ---------- 跳过填充时的专属图表 ----------
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

# ============ 导出填补结果 Excel（学术英文说明） ============
output$download_imputation_excel <- downloadHandler(
  filename = function() {
    paste0("Imputation_Result_", Sys.Date(), ".xlsx")
  },
  content = function(file) {
    message("[DEBUG] download_imputation_excel: starting export")
    if (is.null(processed_data())) {
      showNotification("Please run preprocessing first.", type = "error")
      return()
    }
    
    pre_data <- pre_imputation_with_ids()
    before_mat <- as.matrix(pre_data$before)
    protein_ids <- pre_data$ids
    method <- input$imputation_method
    thr <- input$max_missing_fraction
    k_val <- input$knn_k
    minval <- input$minvalue_fixed
    quant <- input$quantile_prob
    
    if (method == "knn") {
      missing_frac <- rowMeans(is.na(before_mat))
      keep <- missing_frac <= thr
      before_mat <- before_mat[keep, , drop = FALSE]
      protein_ids <- protein_ids[keep]
      message("[DEBUG] KNN mode: after global filter, rows = ", nrow(before_mat))
    } else {
      message("[DEBUG] Method ", method, ": keeping all rows = ", nrow(before_mat))
    }
    
    after_mat <- NULL
    neighbors <- NULL
    correlations <- NULL
    
    tryCatch({
      if (method == "knn") {
        message("[DEBUG] Running impute.knn with k = ", k_val)
        knn_result <- impute::impute.knn(before_mat, k = k_val)
        after_mat <- knn_result$data
        
        cor_mat <- cor(t(before_mat), use = "pairwise.complete.obs")
        diag(cor_mat) <- -Inf
        k_val_eff <- min(k_val, nrow(before_mat) - 1)
        neigh_mat <- matrix(NA_integer_, nrow = nrow(before_mat), ncol = k_val_eff)
        cor_mat_neighbors <- matrix(NA_real_, nrow = nrow(before_mat), ncol = k_val_eff)
        for (i in seq_len(nrow(before_mat))) {
          ord <- order(cor_mat[i, ], decreasing = TRUE)[1:k_val_eff]
          neigh_mat[i, ] <- ord
          cor_mat_neighbors[i, ] <- cor_mat[i, ord]
        }
        rownames(neigh_mat) <- protein_ids
        rownames(cor_mat_neighbors) <- protein_ids
        neighbors <- neigh_mat
        correlations <- cor_mat_neighbors
      } else if (method == "ppca") {
        message("[DEBUG] Running true PPCA imputation (with log2 transform)")
        after_mat <- true_ppca_impute(before_mat)
      } else if (method == "minvalue") {
        message("[DEBUG] Running minvalue imputation with value = ", minval)
        after_mat <- before_mat
        after_mat[is.na(after_mat)] <- minval
      } else if (method == "quantile") {
        message("[DEBUG] Running quantile imputation with prob = ", quant)
        after_mat <- as.matrix(impute_missing_values(before_mat, method = "quantile", quantile_prob = quant))
      } else {
        message("[DEBUG] No imputation (method = none)")
        after_mat <- before_mat
      }
    }, error = function(e) {
      message("[DEBUG] Imputation error: ", e$message)
      showNotification(paste("Imputation failed:", e$message), type = "error")
      return()
    })
    
    if (is.null(after_mat)) {
      showNotification("Imputation failed. Please check parameters.", type = "error")
      return()
    }
    
    after_mat <- as.matrix(after_mat)
    na_positions <- is.na(before_mat)
    total_imputed <- sum(na_positions)
    message("[DEBUG] total imputed cells = ", total_imputed)
    
    wb <- openxlsx::createWorkbook()
    
    # ---- Imputation_Info (English, academic style) ----
    openxlsx::addWorksheet(wb, "Imputation_Info")
    if (method == "knn") {
      desc_text <- paste0("K-Nearest Neighbors imputation (k = ", k_val, "). For each missing value, the algorithm identifies the k proteins with the most similar expression profiles (based on Pearson correlation across available samples). The missing value is replaced by the average of the k neighbors' values in the corresponding sample. Neighbors with missing data in that sample are excluded from the average.")
    } else if (method == "ppca") {
      desc_text <- "Probabilistic Principal Component Analysis (PPCA). The data are first log2(x+1) transformed to satisfy the normality assumption of the model. A PPCA model with 2 principal components is then fitted to the log2-transformed data using the Expectation-Maximization algorithm, which simultaneously estimates the principal components and the missing values. The imputed log2 matrix is back-transformed to the original scale (2^x - 1). This method accounts for the global covariance structure and does not rely on local neighbors."
    } else if (method == "minvalue") {
      desc_text <- paste0("Fixed minimum value imputation. All missing values are replaced with the constant value ", minval, ". This method is suitable when missingness is assumed to result from left-censoring at a detection limit.")
    } else if (method == "quantile") {
      desc_text <- paste0("Quantile imputation. For each sample column, missing values are replaced by the ", quant*100, "th percentile of the observed (non-missing) values in that column. This method assumes that missing values fall below the chosen quantile of the observed distribution.")
    } else {
      desc_text <- "No imputation was applied."
    }
    
    info_text <- c(
      "Missing Value Imputation Export",
      paste("Data source: after missing value filter (mode:", preprocessing_params$missing_filter_mode, 
            ", threshold:", input$max_missing_fraction, ")"),
      "After Inf/Non-finite filter",
      paste("After minimum intensity filter (threshold:", input$min_intensity,
            ", min samples:", input$min_samples_above_intensity, ")"),
      if (method == "knn") paste("Additional filtering: proteins with global missing rate >", thr, "were removed to ensure sufficient observations for KNN."),
      "",
      paste("Imputation method:", method),
      "Parameters:",
      if (method == "knn") paste("  k =", k_val),
      if (method == "ppca") "  nPcs = 2; data were log2(x+1) transformed before imputation and back-transformed afterwards.",
      if (method == "minvalue") paste("  fixed value =", minval),
      if (method == "quantile") paste("  quantile =", quant),
      "",
      "Method description:",
      desc_text,
      "",
      "Workbook sheets:",
      "- Imputation_Info: this information.",
      "- Before_Imputation: matrix before imputation (NA = missing).",
      "- After_Imputation: matrix after imputation. Cells that were imputed are highlighted in red.",
      if (method == "ppca") "- Imputation_Steps: step-by-step description of the PPCA algorithm including log2 transformation.",
      if (!is.null(neighbors)) "- KNN_Neighbors: protein-level list of nearest neighbors (Correlation and Distance).",
      if (!is.null(neighbors)) "- Missing_Imputation_Detail: detailed view of each imputed cell with neighbor values and validity flag."
    )
    info_df <- data.frame(Info = info_text, stringsAsFactors = FALSE)
    openxlsx::writeData(wb, "Imputation_Info", info_df)
    openxlsx::setColWidths(wb, "Imputation_Info", cols = 1, widths = 100)
    
    # ---- Imputation_Steps (English, for PPCA) ----
    if (method == "ppca") {
      openxlsx::addWorksheet(wb, "Imputation_Steps")
      steps <- c(
        "PPCA Imputation Procedure",
        "",
        "1. Log2 transformation: y = log2(x + 1) is applied to the original expression matrix to reduce skewness and approximate normality.",
        "2. Centering and scaling: columns (samples) are centered to mean = 0 and scaled to unit variance.",
        "3. Model fitting: a PPCA model with 2 principal components is fitted to the transformed data via the Expectation-Maximization (EM) algorithm.",
        "   - E-step: estimates the posterior distribution of the latent variables given the current parameter estimates.",
        "   - M-step: updates the model parameters (loadings, residual variance) by maximizing the expected complete-data log-likelihood.",
        "   - Missing values are treated as additional latent variables and estimated during the EM iterations.",
        "4. Imputation in log2 space: after convergence, the complete log2 matrix is reconstructed from the latent scores and loadings.",
        "5. Back-transformation: x = 2^y - 1 restores the original expression scale.",
        "6. Negative values arising from back-transformation are set to zero.",
        "",
        "Advantages:",
        "- Utilizes the global covariance structure across all proteins and samples.",
        "- Does not rely on a limited set of neighbors, reducing bias from outlier proteins.",
        "- Particularly suitable for data where missing values are MAR or MCAR.",
        "",
        "Reference: Tipping, M. E., & Bishop, C. M. (1999). Probabilistic principal component analysis. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 61(3), 611-622."
      )
      steps_df <- data.frame(Step = steps, stringsAsFactors = FALSE)
      openxlsx::writeData(wb, "Imputation_Steps", steps_df)
      openxlsx::setColWidths(wb, "Imputation_Steps", cols = 1, widths = 100)
    }
    
    # ---- Before_Imputation ----
    before_df <- cbind(ProteinID = protein_ids, before_mat, stringsAsFactors = FALSE)
    openxlsx::addWorksheet(wb, "Before_Imputation")
    openxlsx::writeData(wb, "Before_Imputation", before_df)
    
    # ---- After_Imputation (with red) ----
    after_df <- cbind(ProteinID = protein_ids, after_mat, stringsAsFactors = FALSE)
    openxlsx::addWorksheet(wb, "After_Imputation")
    openxlsx::writeData(wb, "After_Imputation", after_df)
    
    style_red <- openxlsx::createStyle(bgFill = "#FF9999")
    for (r in seq_len(nrow(na_positions))) {
      missing_cols <- which(na_positions[r, ])
      if (length(missing_cols) > 0) {
        openxlsx::addStyle(wb, "After_Imputation", style_red,
                           rows = r + 1, cols = missing_cols + 1, gridExpand = TRUE, stack = TRUE)
      }
    }
    message("[DEBUG] Applied red styles to ", total_imputed, " cells")
    
    # ---- KNN specific sheets ----
    if (method == "knn" && !is.null(neighbors)) {
      message("[DEBUG] Building KNN_Neighbors sheet")
      neighbor_list <- list()
      n_k <- ncol(neighbors)
      for (i in seq_len(nrow(neighbors))) {
        pid <- protein_ids[i]
        for (j in seq_len(n_k)) {
          neighbor_idx <- neighbors[i, j]
          if (is.na(neighbor_idx) || neighbor_idx < 1 || neighbor_idx > length(protein_ids)) next
          neighbor_id <- protein_ids[neighbor_idx]
          cor_val <- if (!is.null(correlations)) round(correlations[i, j], 4) else NA
          neighbor_list[[length(neighbor_list) + 1]] <- data.frame(
            ProteinID = pid,
            NeighborRank = j,
            NeighborProteinID = neighbor_id,
            Correlation = cor_val,
            Distance = if (!is.na(cor_val)) round(1 - cor_val, 4) else NA,
            stringsAsFactors = FALSE
          )
        }
      }
      if (length(neighbor_list) > 0) {
        neighbor_df <- do.call(rbind, neighbor_list)
        openxlsx::addWorksheet(wb, "KNN_Neighbors")
        openxlsx::writeData(wb, "KNN_Neighbors", neighbor_df)
        message("[DEBUG] KNN_Neighbors sheet written with ", nrow(neighbor_df), " rows")
      }
      
      # Missing_Imputation_Detail
      detail_list <- list()
      sample_names <- colnames(before_mat)
      for (i in seq_len(nrow(before_mat))) {
        pid <- protein_ids[i]
        missing_cols <- which(na_positions[i, ])
        if (length(missing_cols) == 0) next
        neigh_idx <- neighbors[i, ]
        neigh_cor <- correlations[i, ]
        for (col in missing_cols) {
          imputed_val <- after_mat[i, col]
          neighbor_values <- before_mat[neigh_idx, col]
          valid_neighbors <- !is.na(neighbor_values)
          for (rank in seq_along(neigh_idx)) {
            nv <- neighbor_values[rank]
            detail_list[[length(detail_list) + 1]] <- data.frame(
              ProteinID = pid,
              Sample = sample_names[col],
              ImputedValue = round(imputed_val, 4),
              NeighborRank = rank,
              NeighborProteinID = protein_ids[neigh_idx[rank]],
              NeighborValue = if (!is.na(nv)) round(nv, 4) else NA,
              Valid = valid_neighbors[rank],
              Correlation = round(neigh_cor[rank], 4),
              stringsAsFactors = FALSE
            )
          }
        }
      }
      if (length(detail_list) > 0) {
        detail_df <- do.call(rbind, detail_list)
        openxlsx::addWorksheet(wb, "Missing_Imputation_Detail")
        openxlsx::writeData(wb, "Missing_Imputation_Detail", detail_df)
        message("[DEBUG] Missing_Imputation_Detail sheet written with ", nrow(detail_df), " rows")
      }
    }
    
    openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    message("[DEBUG] download_imputation_excel: export completed")
  }
)