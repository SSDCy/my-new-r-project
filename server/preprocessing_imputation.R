# server/preprocessing_imputation.R
message("[DEBUG] preprocessing_imputation.R loaded - all imputation logic protected against missing input")

# 保护函数：如果 input$imputation_method 缺失，返回 NULL 或跳过
imputation_method_safe <- function() {
  if (is.null(input$imputation_method)) return("none")
  else input$imputation_method
}

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

# ---------- KNN 输入数据（仅在 method 为 knn 时计算） ----------
knn_input_data <- reactive({
  method <- imputation_method_safe()
  if (method != "knn" || is.null(pre_imputation_with_ids())) return(NULL)
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

# ---------- KNN 邻居信息 ----------
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
  
  list(neighbors = neighbors, correlations = cor_vals, distances = dist_vals, protein_ids = protein_ids, k = k_val)
})

# ---------- 下拉框更新（KNN 邻居查找） ----------
observe({
  data <- knn_neighbors_data()
  if (is.null(data)) return()
  ids <- data$protein_ids
  updateSelectizeInput(session, "knn_lookup_protein",
                       choices = ids,
                       server = TRUE,
                       selected = if (length(ids) > 0) ids[1] else NULL)
})

# ---------- 邻居表格输出 ----------
output$knn_lookup_table <- DT::renderDT({
  data <- knn_neighbors_data()
  req(data, input$knn_lookup_protein)
  
  pid <- input$knn_lookup_protein
  idx <- which(data$protein_ids == pid)
  if (length(idx) == 0) {
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
  DT::datatable(df, options = list(dom = 't', pageLength = 10), rownames = FALSE) |>
    DT::formatStyle("Distance", color = "gray")
})

# ---------- 真正的PPCA填补函数 ----------
true_ppca_impute <- function(data, nPcs = 2) {
  message("[DEBUG] true_ppca_impute: starting PPCA with log2 transformation")
  if (!requireNamespace("pcaMethods", quietly = TRUE))
    stop("pcaMethods package required for PPCA. Run BiocManager::install('pcaMethods')")
  
  data_matrix <- as.matrix(data)
  log_data <- log2(data_matrix + 1)
  
  na_rows <- which(rowSums(is.na(log_data)) == ncol(log_data))
  na_cols <- which(colSums(is.na(log_data)) == nrow(log_data))
  if (length(na_rows) > 0 || length(na_cols) > 0) {
    if (length(na_rows) > 0) log_data <- log_data[-na_rows, , drop = FALSE]
    if (length(na_cols) > 0) log_data <- log_data[, -na_cols, drop = FALSE]
  }
  
  pc <- pcaMethods::pca(log_data, method = "ppca", nPcs = nPcs, scale = "uv", center = TRUE)
  imputed_log <- as.matrix(pcaMethods::completeObs(pc))
  
  imputed_original <- 2^imputed_log - 1
  imputed_original[imputed_original < 0] <- 0
  return(imputed_original)
}

# ---------- PPCA 可视化数据 ----------
ppca_visualization_data <- reactive({
  method <- imputation_method_safe()
  if (method != "ppca" || is.null(processed_data())) return(NULL)
  pre <- pre_imputation_with_ids()
  if (is.null(pre)) return(NULL)
  mat <- as.matrix(pre$before)
  
  if (!requireNamespace("pcaMethods", quietly = TRUE)) return(NULL)
  
  tryCatch({
    log_mat <- log2(mat + 1)
    pc <- pcaMethods::pca(log_mat, method = "ppca", nPcs = 2, scale = "uv", center = TRUE)
    scores <- as.data.frame(pcaMethods::scores(pc))
    scores$Sample <- rownames(scores)
    loadings <- as.data.frame(pcaMethods::loadings(pc))
    loadings$Protein <- rownames(loadings)
    imputed_log <- as.matrix(pcaMethods::completeObs(pc))
    imputed_original <- 2^imputed_log - 1
    imputed_original[imputed_original < 0] <- 0
    list(scores = scores, loadings = loadings, imputed = imputed_original, original = mat)
  }, error = function(e) {
    NULL
  })
})

output$ppca_score_plot <- renderPlot({
  data <- ppca_visualization_data()
  if (is.null(data)) {
    plot.new(); text(0.5, 0.5, "PPCA visualization not available. Please run preprocessing first with PPCA method.")
    return()
  }
  ggplot(data$scores, aes(x = PC1, y = PC2)) +
    geom_point(size = 3, color = "#3498db") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(title = "PPCA Score Plot", x = "PC1", y = "PC2") +
    theme_bw()
})

output$ppca_imputation_hist <- renderPlot({
  data <- ppca_visualization_data()
  if (is.null(data)) {
    plot.new(); text(0.5, 0.5, "Histogram not available.")
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
  ggplot(df, aes(x = Value, fill = Type)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    scale_fill_manual(values = c("Original (non-missing)" = "#3498db", "Imputed" = "#2ecc71")) +
    labs(title = "Distribution of Original vs Imputed Values", x = "Expression Value", y = "Frequency") +
    theme_bw() + theme(legend.position = "bottom")
})

# ---------- 分位数填补可视化数据 ----------
quantile_imputation_data <- reactive({
  method <- imputation_method_safe()
  if (method != "quantile" || is.null(processed_data())) return(NULL)
  before_mat <- as.matrix(pre_imputation_matrix())
  q <- input$quantile_prob
  n_samples <- ncol(before_mat)
  thresholds <- sapply(1:n_samples, function(j) quantile(before_mat[, j], probs = q, na.rm = TRUE))
  missing_counts <- sapply(1:n_samples, function(j) sum(is.na(before_mat[, j])))
  
  sample_names_full <- colnames(before_mat)
  sample_names_short <- extract_sample_names(sample_names_full)
  
  groups <- rep("Unassigned", n_samples)
  if (!is.null(rv$sample_info) && "Group" %in% colnames(rv$sample_info)) {
    si <- rv$sample_info
    si$short <- extract_sample_names(rownames(si))
    idx <- match(sample_names_short, si$short)
    matched <- !is.na(idx)
    if (any(matched)) {
      groups[matched] <- si$Group[idx[matched]]
    }
  }
  
  list(thresholds = thresholds, missing_counts = missing_counts, 
       sample_names = sample_names_short, sample_names_full = sample_names_full, 
       q = q, groups = groups, mat = before_mat)
})

output$quantile_threshold_plot <- renderPlot({
  data <- quantile_imputation_data()
  if (is.null(data)) {
    plot.new(); text(0.5, 0.5, "Quantile visualization not available. Please run preprocessing first with Quantile method.")
    return()
  }
  df <- data.frame(
    Sample = data$sample_names,
    Threshold = data$thresholds,
    Missing = data$missing_counts,
    Group = data$groups,
    stringsAsFactors = FALSE
  )
  df <- df[order(df$Threshold, decreasing = TRUE), ]
  df$Sample <- factor(df$Sample, levels = df$Sample)
  df$SampleLabel <- paste0(df$Sample, " (", df$Missing, ")")
  
  groups <- unique(df$Group)
  if (length(groups) <= 8) {
    group_colors <- setNames(RColorBrewer::brewer.pal(length(groups), "Set1"), groups)
  } else {
    group_colors <- setNames(rainbow(length(groups)), groups)
  }
  
  y_max <- max(df$Threshold) * 1.12
  y_min <- min(df$Threshold) * 0.98
  
  ggplot(df, aes(x = Sample, y = Threshold, fill = Group)) +
    geom_col(alpha = 0.85, width = 0.7) +
    scale_fill_manual(values = group_colors, name = "Group") +
    scale_x_discrete(labels = setNames(df$SampleLabel, df$Sample)) +
    labs(title = paste0("Imputation Threshold (", data$q*100, "th Percentile) per Sample"),
         x = "Sample (missing count)", y = "Threshold Value") +
    coord_cartesian(ylim = c(y_min, y_max)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          legend.position = "bottom")
})

output$quantile_threshold_table <- renderTable({
  data <- quantile_imputation_data()
  req(data)
  df <- data.frame(
    Sample = data$sample_names_full,
    Threshold = format(round(data$thresholds, 2), big.mark = ",", scientific = FALSE, trim = TRUE),
    `Missing Count` = data$missing_counts,
    Group = data$groups,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  df
}, striped = TRUE, bordered = TRUE, width = "100%")

observe({
  data <- quantile_imputation_data()
  if (is.null(data)) return()
  updateSelectInput(session, "quantile_verify_sample", choices = data$sample_names_full)
})

output$quantile_distribution_plot <- renderPlot({
  data <- quantile_imputation_data()
  req(data, input$quantile_verify_sample)
  sample_idx <- which(data$sample_names_full == input$quantile_verify_sample)
  if (length(sample_idx) != 1) {
    plot.new(); text(0.5, 0.5, "Sample not found.")
    return()
  }
  values <- data$mat[, sample_idx]
  non_missing <- values[!is.na(values)]
  threshold <- data$thresholds[sample_idx]
  q <- data$q
  
  if (length(non_missing) == 0) {
    plot.new(); text(0.5, 0.5, "No non-missing values in this sample.")
    return()
  }
  
  actual_pct <- mean(non_missing <= threshold) * 100
  expected_pct <- q * 100
  
  df <- data.frame(Value = non_missing)
  ggplot(df, aes(x = Value)) +
    geom_histogram(bins = 50, fill = "#3498db", alpha = 0.7, boundary = 0) +
    geom_vline(xintercept = threshold, color = "red", linewidth = 1.5, linetype = "dashed") +
    annotate("text", x = threshold, y = Inf, 
             label = paste0("Threshold = ", format(round(threshold, 2), big.mark = ",", scientific = FALSE)),
             color = "red", vjust = 2, hjust = -0.1, size = 3.5) +
    labs(title = paste0("Distribution of Non‑Missing Values in ", input$quantile_verify_sample),
         subtitle = paste0(format(round(threshold, 2), big.mark = ",", scientific = FALSE),
                           " is the ", expected_pct,
                           "th percentile. Actual values ≤ threshold: ", round(actual_pct, 1), "% (expected ", expected_pct, "%)"),
         x = "Expression Value", y = "Frequency") +
    theme_bw()
})

output$quantile_threshold_position <- renderPrint({
  data <- quantile_imputation_data()
  req(data, input$quantile_verify_sample)
  sample_idx <- which(data$sample_names_full == input$quantile_verify_sample)
  if (length(sample_idx) != 1) {
    cat("Sample not found.\n")
    return()
  }
  values <- data$mat[, sample_idx]
  non_missing <- values[!is.na(values)]
  threshold <- data$thresholds[sample_idx]
  q <- data$q
  
  sorted_vals <- sort(non_missing)
  n <- length(sorted_vals)
  p <- (n - 1) * q + 1
  k <- floor(p)
  lambda <- p - k
  
  if (k >= 1 && k < n) {
    xk <- sorted_vals[k]
    xk1 <- sorted_vals[k + 1]
  } else if (k < 1) {
    xk <- sorted_vals[1]
    xk1 <- sorted_vals[1]
  } else {
    xk <- sorted_vals[n]
    xk1 <- sorted_vals[n]
  }
  
  fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
  
  cat("Threshold:", fmt(round(threshold, 2)), "\n")
  cat("Total non‑missing values:", n, "\n\n")
  cat("Formula (R default type=7):\n")
  cat("  Position p = (n-1) * q + 1\n")
  cat("  p = (", n, " - 1) * ", q, " + 1 = ", round(p, 4), "\n", sep = "")
  cat("  k = floor(p) = ", k, "\n", sep = "")
  cat("  λ = p - k = ", round(lambda, 4), "\n", sep = "")
  cat("  x[", k, "] = ", fmt(xk), "\n", sep = "")
  cat("  x[", k+1, "] = ", fmt(xk1), "\n", sep = "")
  cat("  Threshold = x[k] + λ * (x[k+1] - x[k])\n")
  cat("  Threshold = ", fmt(xk), " + ", round(lambda, 4), " * (", fmt(xk1), " - ", fmt(xk), ")\n", sep = "")
  cat("  Threshold = ", fmt(round(threshold, 2)), "\n\n", sep = "")
  
  pos <- findInterval(threshold, sorted_vals)
  if (pos == 0) {
    cat("Threshold is below the smallest value (", fmt(sorted_vals[1]), ").\n")
  } else if (pos == n) {
    cat("Threshold is above the largest value (", fmt(sorted_vals[n]), ").\n")
  } else {
    cat("Threshold lies between rank", pos, "(", fmt(sorted_vals[pos]), ")",
        "and rank", pos+1, "(", fmt(sorted_vals[pos+1]), ").\n")
  }
})

output$quantile_raw_data_table <- DT::renderDT({
  data <- quantile_imputation_data()
  req(data, input$quantile_verify_sample)
  sample_idx <- which(data$sample_names_full == input$quantile_verify_sample)
  if (length(sample_idx) != 1) return(DT::datatable(data.frame(Message = "Sample not found")))
  values <- data$mat[, sample_idx]
  non_missing <- values[!is.na(values)]
  threshold <- data$thresholds[sample_idx]
  
  if (length(non_missing) == 0) return(DT::datatable(data.frame(Message = "No non-missing values")))
  sorted_vals <- sort(non_missing)
  df <- data.frame(Rank = 1:length(sorted_vals), Value = sorted_vals, stringsAsFactors = FALSE)
  highlight_val <- sorted_vals[findInterval(threshold, sorted_vals) + 1]
  if (is.na(highlight_val) || highlight_val < threshold) {
    highlight_val <- sorted_vals[length(sorted_vals)]
  }
  DT::datatable(df, options = list(pageLength = 25, scrollY = "300px", dom = 'ftip'), rownames = FALSE) |>
    DT::formatStyle("Value", target = "row", backgroundColor = DT::styleEqual(highlight_val, "#FFFFCC"))
})

# ---------- 填补比较数据（当 imputation 为 none 时不再计算） ----------
imputation_comparison_data <- reactive({
  # 禁用填补，返回 NULL
  method <- imputation_method_safe()
  if (method == "none" || is.null(processed_data())) return(NULL)
  req(processed_data(), pre_imputation_matrix())
  # ... 原有逻辑，但不会被触发
})

output$imputation_stats_text <- renderPrint({
  if (is.null(processed_data())) {
    cat("Please run preprocessing first.\n")
    return()
  }
  cat("Imputation method: none (imputation disabled)\n")
})

output$imputation_boxplot <- renderPlot({
  plot.new()
  text(0.5, 0.5, "Imputation has been disabled.")
})

output$imputation_pca_plot <- renderPlot({
  plot.new()
  text(0.5, 0.5, "Imputation has been disabled.")
})

output$imputation_qq_plot <- renderPlot({
  plot.new()
  text(0.5, 0.5, "Imputation has been disabled.")
})

output$imputation_summary_table <- DT::renderDT({
  DT::datatable(data.frame(Message = "Imputation disabled."))
})

output$download_imputation_table <- downloadHandler(
  filename = function() "imputation_disabled.xlsx",
  content = function(file) showNotification("Imputation disabled.", type = "error")
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
                     main = "Missing Value Heatmap",
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

# ---------- 导出填补结果 Excel（禁用） ----------
output$download_imputation_excel <- downloadHandler(
  filename = function() "Imputation_Result_disabled.xlsx",
  content = function(file) {
    showNotification("Imputation has been disabled.", type = "error")
  }
)

message("[DEBUG] preprocessing_imputation.R: all imputation logic protected against missing input")