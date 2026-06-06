# server/sample_correlation_plot.R
message("[DEBUG] sample_correlation_plot.R loading... (log2, Z-score, enhanced debug)")

# ========== 原 Plots 相关（保留但 UI 已移除此选项卡） ==========
output$sample_cor_data_source_note <- renderUI({
  if (is.null(norm_data_full())) {
    return(div(style = "color: #e74c3c;", "Normalized data not available. Please run preprocessing."))
  }
  div(style = "margin-bottom: 10px; color: #27ae60; font-weight: bold;",
      icon("check-circle"), " Data source: Normalized expression data (total intensity normalization + log2 transformation)")
})

output$sample_cor_preprocess_steps <- renderUI({
  # 保留但不再使用
})

# ========== 数据质量页面使用的样本相关性（log2，Z-score，默认聚类） ==========
dq_sample_cor_data <- reactive({
  req(dq_expr_matrix())
  mat <- dq_expr_matrix()
  samples <- selected_samples()
  samples <- intersect(samples, colnames(mat))
  if (length(samples) < 2) {
    message("[DEBUG] dq_sample_cor_data: need at least 2 samples, currently ", length(samples))
    return(NULL)
  }
  mat <- mat[, samples, drop = FALSE]
  
  mat <- mat[rowSums(!is.na(mat)) > 1, , drop = FALSE]
  message("[DEBUG] dq_sample_cor_data: after removing all-NA rows: ", nrow(mat), " x ", ncol(mat))
  
  if (any(is.na(mat))) {
    message("[DEBUG] dq_sample_cor_data: imputing missing values with Quantile (1%)")
    mat <- as.matrix(impute_missing_values(as.data.frame(mat), method = "quantile", quantile_prob = 0.01))
  }
  
  # log2 转换
  log_expr <- log2(mat + 1)
  
  # Z-score 标准化
  z_expr <- t(scale(t(log_expr)))
  z_expr[!is.finite(z_expr)] <- 0
  
  # Top500 变异蛋白
  row_vars <- apply(z_expr, 1, var, na.rm = TRUE)
  n_keep <- min(500, nrow(z_expr))
  top_idx <- order(row_vars, decreasing = TRUE)[1:n_keep]
  z_expr <- z_expr[top_idx, , drop = FALSE]
  message("[DEBUG] dq_sample_cor_data: selected top ", n_keep, " variable proteins (Z-score)")
  
  # 计算相关性
  cor_mat <- cor(z_expr, use = "complete.obs")
  cor_mat[is.na(cor_mat)] <- 0
  
  # ========== 增强调试信息 ==========
  message("[DEBUG] dq_sample_cor_data: --- Correlation Matrix ---")
  # 打印完整矩阵（21x21 在控制台中可读）
  print(round(cor_mat, 4))
  
  message("[DEBUG] dq_sample_cor_data: --- Matrix Summary ---")
  cor_vals <- cor_mat[upper.tri(cor_mat)]
  message(sprintf("  min = %.4f", min(cor_vals)))
  message(sprintf("  25th percentile = %.4f", quantile(cor_vals, 0.25)))
  message(sprintf("  median = %.4f", median(cor_vals)))
  message(sprintf("  75th percentile = %.4f", quantile(cor_vals, 0.75)))
  message(sprintf("  max = %.4f", max(cor_vals)))
  message(sprintf("  fraction > 0.5 = %.2f", mean(cor_vals > 0.5)))
  
  # 每个样本的相关性概况
  sample_names <- colnames(cor_mat)
  sample_stats <- data.frame(
    Sample = sample_names,
    MinCor = apply(cor_mat, 1, min),
    MedianCor = apply(cor_mat, 1, median),
    MaxCor = apply(cor_mat, 1, max)
  )
  message("[DEBUG] dq_sample_cor_data: --- Per Sample Correlation Stats ---")
  print(sample_stats)
  
  # 注释信息
  ann_col <- data.frame(row.names = sample_names)
  ann_colors <- list()
  
  if (!is.null(rv$sample_info)) {
    si <- rv$sample_info
    si$ShortName <- extract_sample_names(rownames(si))
    idx <- match(sample_names, si$ShortName)
    
    if ("SubGroup" %in% colnames(si)) {
      groups <- si$SubGroup[idx]
      groups[is.na(groups)] <- "Unassigned"
      ann_col$SubGroup <- factor(groups)
      groups_uniq <- levels(ann_col$SubGroup)
      ann_colors$SubGroup <- get_group_colors(groups_uniq)
      message("[DEBUG] dq_sample_cor_data: SubGroup annotation added, levels: ", paste(groups_uniq, collapse=", "))
    } else if ("Group" %in% colnames(si)) {
      groups <- si$Group[idx]
      groups[is.na(groups)] <- "Unassigned"
      ann_col$Group <- factor(groups)
      groups_uniq <- levels(ann_col$Group)
      ann_colors$Group <- get_group_colors(groups_uniq)
      message("[DEBUG] dq_sample_cor_data: SubGroup not found, used Group annotation instead")
    }
  }
  
  if (ncol(ann_col) == 0) ann_col <- NULL
  if (length(ann_colors) == 0) ann_colors <- NULL
  
  message("[DEBUG] dq_sample_cor_data: correlation matrix generated, dim = ", nrow(cor_mat), "x", ncol(cor_mat))
  list(cor_mat = cor_mat, ann_col = ann_col, ann_colors = ann_colors, samples = sample_names)
})

output$dq_sample_cor_heatmap <- renderPlot({
  dat <- dq_sample_cor_data()
  if (is.null(dat)) {
    plot.new()
    text(0.5, 0.5, "Not enough samples selected")
    return()
  }
  cor_vals <- dat$cor_mat[upper.tri(dat$cor_mat)]
  min_cor <- min(cor_vals, na.rm = TRUE)
  max_cor <- max(cor_vals, na.rm = TRUE)
  if (min_cor == max_cor) {
    min_cor <- min_cor - 0.01
    max_cor <- max_cor + 0.01
  }
  palette <- colorRampPalette(c("blue", "white", "red"))(255)
  breaks <- seq(min_cor, max_cor, length.out = 256)
  
  message(sprintf("[DEBUG] dq_sample_cor_heatmap: color range: min=%.4f, max=%.4f", min_cor, max_cor))
  
  pheatmap::pheatmap(dat$cor_mat,
                     color = palette,
                     breaks = breaks,
                     legend_breaks = round(c(min_cor, (min_cor+max_cor)/2, max_cor), 2),
                     legend_labels = c(format(min_cor, digits=2), format((min_cor+max_cor)/2, digits=2), format(max_cor, digits=2)),
                     clustering_distance_rows = as.dist(1 - dat$cor_mat),
                     clustering_distance_cols = as.dist(1 - dat$cor_mat),
                     clustering_method = "ward.D2",
                     show_rownames = TRUE,
                     show_colnames = TRUE,
                     fontsize_row = 9,
                     fontsize_col = 9,
                     angle_col = 45,
                     annotation_col = dat$ann_col,
                     annotation_colors = dat$ann_colors,
                     annotation_legend = TRUE,
                     legend = TRUE,
                     main = "Sample Correlation (Z-score per protein, log2, Top500)")
})

# ========== 原 Plots 相关（保留但已从 UI 移除） ==========
sample_cor_data <- eventReactive(input$generate_sample_cor, {
  NULL
})

output$sample_cor_heatmap <- renderPlot({
  plot.new()
  text(0.5,0.5,"This plot has been moved to Data Quality Analysis.")
})

output$download_sample_cor_png <- downloadHandler(
  filename = function() "sample_correlation.png",
  content = function(file) {
    dat <- dq_sample_cor_data()
    if (!is.null(dat)) {
      png(file, width = 900, height = 700, res = 150)
      palette <- colorRampPalette(c("blue", "white", "red"))(255)
      cor_vals <- dat$cor_mat[upper.tri(dat$cor_mat)]
      min_cor <- min(cor_vals, na.rm = TRUE)
      max_cor <- max(cor_vals, na.rm = TRUE)
      if (min_cor == max_cor) { min_cor <- min_cor - 0.01; max_cor <- max_cor + 0.01 }
      breaks <- seq(min_cor, max_cor, length.out = 256)
      pheatmap::pheatmap(dat$cor_mat,
                         color = palette,
                         breaks = breaks,
                         legend_breaks = round(c(min_cor, (min_cor+max_cor)/2, max_cor), 2),
                         legend_labels = c(format(min_cor, digits=2), format((min_cor+max_cor)/2, digits=2), format(max_cor, digits=2)),
                         clustering_distance_rows = as.dist(1 - dat$cor_mat),
                         clustering_distance_cols = as.dist(1 - dat$cor_mat),
                         clustering_method = "ward.D2",
                         show_rownames = TRUE,
                         show_colnames = TRUE,
                         fontsize_row = 9,
                         fontsize_col = 9,
                         angle_col = 45,
                         annotation_col = dat$ann_col,
                         annotation_colors = dat$ann_colors,
                         annotation_legend = TRUE,
                         legend = TRUE,
                         main = "Sample Correlation (Z-score per protein, log2)")
      dev.off()
    }
  }
)

message("[DEBUG] sample_correlation_plot.R loaded successfully (enhanced debug)")