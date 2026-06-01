# server/preprocessing_imputation.R
message("[DEBUG] preprocessing_imputation.R loaded - imputation detail with valid neighbors and correct average")

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

# ---------- 全局缺失率过滤后的矩阵（用于邻居搜索，保证一致性） ----------
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

# ---------- 全局邻居信息（基于皮尔逊相关系数，快速稳定） ----------
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

# ============ 导出填补结果 Excel（含详细填补细节及 Valid 标记） ============
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
    thr <- input$max_missing_fraction
    
    missing_frac <- rowMeans(is.na(before_mat))
    keep <- missing_frac <= thr
    before_mat <- before_mat[keep, , drop = FALSE]
    protein_ids <- protein_ids[keep]
    message("[DEBUG] download_imputation_excel: after global filter, rows = ", nrow(before_mat))
    
    method <- input$imputation_method
    k_val <- input$knn_k
    minval <- input$minvalue_fixed
    quant <- input$quantile_prob
    
    after_mat <- NULL
    neighbors <- NULL
    correlations <- NULL
    
    tryCatch({
      if (method == "knn") {
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
        after_mat <- as.matrix(impute_missing_values(before_mat, method = "ppca"))
      } else if (method == "minvalue") {
        after_mat <- before_mat
        after_mat[is.na(after_mat)] <- minval
      } else if (method == "quantile") {
        after_mat <- as.matrix(impute_missing_values(before_mat, method = "quantile", quantile_prob = quant))
      } else {
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
    
    # ---------- 构建详细填补表（新增 Valid 列） ----------
    detail_list <- list()
    if (method == "knn" && !is.null(neighbors)) {
      sample_names <- colnames(before_mat)
      for (i in seq_len(nrow(before_mat))) {
        pid <- protein_ids[i]
        missing_cols <- which(na_positions[i, ])
        if (length(missing_cols) == 0) next
        neigh_idx <- neighbors[i, ]
        neigh_cor <- correlations[i, ]
        for (col in missing_cols) {
          neighbor_values <- before_mat[neigh_idx, col]
          valid_neighbors <- !is.na(neighbor_values)
          # 计算有效邻居的平均值
          if (any(valid_neighbors)) {
            imputed_value <- mean(neighbor_values[valid_neighbors])
          } else {
            imputed_value <- NA
          }
          for (rank in seq_along(neigh_idx)) {
            nv <- neighbor_values[rank]
            detail_list[[length(detail_list) + 1]] <- data.frame(
              ProteinID = pid,
              Sample = sample_names[col],
              ImputedValue = round(imputed_value, 4),
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
      message("[DEBUG] Built detail list with ", length(detail_list), " rows")
    }
    
    # 公式说明
    if (method == "knn") {
      formula_text <- paste0("KNN imputation (k=", k_val, "): each missing value is replaced by the average of the k nearest neighbors (based on Pearson correlation). Neighbors with missing values in the corresponding sample are excluded from the average (Valid = FALSE).")
    } else {
      formula_text <- paste0("Imputation method: ", method)
    }
    
    wb <- openxlsx::createWorkbook()
    
    # ---- Imputation_Info ----
    openxlsx::addWorksheet(wb, "Imputation_Info")
    info_text <- c(
      "Missing Value Imputation Export",
      paste("Data source: After Missing Value Filter (mode:", preprocessing_params$missing_filter_mode,
            ", threshold:", input$max_missing_fraction, ")"),
      "After Inf/Non-finite Filter",
      paste("After Minimum Intensity Filter (threshold:", input$min_intensity,
            ", min samples:", input$min_samples_above_intensity, ")"),
      paste("Note: Data further filtered to ensure global missing rate ≤", thr, "for KNN compatibility."),
      "",
      paste("Imputation method:", method),
      "Parameters:",
      if (method == "knn") paste("  k =", k_val),
      "",
      "Formula:",
      formula_text,
      "",
      "Sheets:",
      "- Imputation_Info: this information",
      "- Before_Imputation: matrix with NAs before imputation",
      "- After_Imputation: filled matrix (imputed cells highlighted in red)",
      "- KNN_Neighbors: protein-level neighbor list (Correlation & Distance)",
      if (!is.null(detail_list)) "- Missing_Imputation_Detail: each imputed cell with neighbor values and Valid flag",
      "",
      "Understanding Missing_Imputation_Detail:",
      "  - ProteinID & Sample: location of missing value",
      "  - ImputedValue: final filled value (average of valid neighbors)",
      "  - NeighborRank: rank of neighbor (1 = most correlated)",
      "  - NeighborProteinID: ID of the neighbor protein",
      "  - NeighborValue: expression of neighbor in this sample (NA if missing)",
      "  - Valid: TRUE if neighbor value was used in average; FALSE if neighbor was missing",
      "  - Correlation: Pearson correlation between target and neighbor protein"
    )
    info_df <- data.frame(Info = info_text, stringsAsFactors = FALSE)
    openxlsx::writeData(wb, "Imputation_Info", info_df)
    openxlsx::setColWidths(wb, "Imputation_Info", cols = 1, widths = 90)
    
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
    
    # ---- KNN_Neighbors (protein-level) ----
    if (!is.null(neighbors)) {
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
    }
    
    # ---- Missing_Imputation_Detail ----
    if (length(detail_list) > 0) {
      detail_df <- do.call(rbind, detail_list)
      openxlsx::addWorksheet(wb, "Missing_Imputation_Detail")
      openxlsx::writeData(wb, "Missing_Imputation_Detail", detail_df)
      message("[DEBUG] Missing_Imputation_Detail sheet written with ", nrow(detail_df), " rows")
    }
    
    openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    message("[DEBUG] download_imputation_excel: export completed")
  }
)