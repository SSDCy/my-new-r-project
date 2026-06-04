# server/sample_correlation_plot.R
message("[DEBUG] sample_correlation_plot.R loading...")

# ---------- 数据源说明 ----------
output$sample_cor_data_source_note <- renderUI({
  if (is.null(norm_data_full())) {
    return(div(style = "color: #e74c3c;", "Normalized data not available. Please run preprocessing."))
  }
  div(style = "margin-bottom: 10px; color: #27ae60; font-weight: bold;",
      icon("check-circle"), " Data source: Normalized expression data (total intensity normalization + log2 transformation)")
})

# ---------- 步骤指示器（可折叠） ----------
output$sample_cor_preprocess_steps <- renderUI({
  steps <- list()
  steps <- c(steps, paste0("Missing Value Filter: threshold = ", input$max_missing_fraction %||% 0.5,
                           ", mode = ", preprocessing_params$missing_filter_mode %||% "global"))
  steps <- c(steps, paste0("Minimum Intensity Filter: threshold = ", input$min_intensity,
                           ", min samples = ", input$min_samples_above_intensity %||% 1))
  imp <- preprocessing_params$imputation_method %||% "none"
  if (imp == "none") {
    steps <- c(steps, "Missing Value Imputation: none (Quantile 1% auto-applied if missing values present)")
  } else {
    steps <- c(steps, paste0("Missing Value Imputation: ", imp))
  }
  if (isTRUE(preprocessing_params$batch_performed)) {
    steps <- c(steps, "Batch Correction (ComBat): applied")
  } else {
    steps <- c(steps, "Batch Correction: not applied")
  }
  steps <- c(steps, "Normalization: Total intensity normalization (baseline sample)")
  steps <- c(steps, "Data source: Normalized expression data (Norm_LFQ intensity columns) + log2 transformation")
  steps <- c(steps, "Top 500 variable proteins selected; missing values imputed with Quantile 1% (if any)")
  
  step_tags <- lapply(seq_along(steps), function(i) {
    tagList(
      if (i > 1) tags$span(style = "font-size: 20px; color: #e67e22; margin: 0 8px;", "→"),
      tags$span(style = "background: #e8f0fe; padding: 6px 12px; border-radius: 15px; font-size: 13px;", steps[[i]])
    )
  })
  tags$details(
    tags$summary("Data preprocessing steps for Sample Correlation", style = "cursor: pointer; font-weight: bold; color: #2c3e50; margin-bottom: 10px;"),
    div(style = "display: flex; flex-wrap: wrap; align-items: center;", do.call(tagList, step_tags))
  )
})

# ---------- 计算样本相关性矩阵 ----------
sample_cor_data <- eventReactive(input$generate_sample_cor, {
  nd <- norm_data_full()
  if (is.null(nd)) {
    message("[DEBUG] sample_cor: norm_data_full() is NULL")
    return(NULL)
  }
  
  # 提取 Norm_ 开头的强度列
  norm_cols <- grep("^Norm_LFQ intensity ", colnames(nd), value = TRUE)
  if (length(norm_cols) == 0) {
    message("[DEBUG] sample_cor: no Norm_ columns found")
    return(NULL)
  }
  
  expr_mat <- as.matrix(nd[, norm_cols])
  # 移除全为缺失值的行
  expr_mat <- expr_mat[rowSums(!is.na(expr_mat)) > 1, , drop = FALSE]
  message("[DEBUG] sample_cor: matrix dim after removing all-NA rows: ", nrow(expr_mat), " x ", ncol(expr_mat))
  
  # 填补缺失值（若存在）——使用 1% 分位数填补
  if (any(is.na(expr_mat))) {
    message("[DEBUG] sample_cor: imputing missing values with Quantile (1%)")
    expr_mat <- as.matrix(impute_missing_values(
      as.data.frame(expr_mat), method = "quantile", quantile_prob = 0.01
    ))
  }
  
  # log2 转换
  log_expr <- log2(expr_mat + 1)
  
  # 选择 top 500 高变异蛋白
  row_vars <- apply(log_expr, 1, var, na.rm = TRUE)
  n_keep <- min(500, nrow(log_expr))
  top_idx <- order(row_vars, decreasing = TRUE)[1:n_keep]
  log_expr <- log_expr[top_idx, , drop = FALSE]
  message("[DEBUG] sample_cor: selected top ", n_keep, " variable proteins")
  
  # 计算相关性
  cor_mat <- cor(log_expr, use = "complete.obs")
  cor_mat[is.na(cor_mat)] <- 0
  
  # 获取样本短名
  colnames(cor_mat) <- gsub("^Norm_LFQ intensity ", "", colnames(cor_mat))
  rownames(cor_mat) <- colnames(cor_mat)
  
  # 生成注释信息（分组 + 批次）
  sample_short <- colnames(cor_mat)
  ann_col <- data.frame(row.names = sample_short)
  
  # 分组信息
  if (!is.null(rv$groups) && length(rv$groups) > 0) {
    groups <- rep("Unassigned", length(sample_short))
    for (i in seq_along(sample_short)) {
      sn <- sample_short[i]
      for (gn in names(rv$groups)) {
        if (sn %in% rv$groups[[gn]]) {
          groups[i] <- gn
          break
        }
      }
    }
    ann_col$Group <- groups
  }
  
  # 批次信息
  if (!is.null(rv$sample_info) && "Batch" %in% colnames(rv$sample_info)) {
    si <- rv$sample_info
    info_full <- rownames(si)
    info_short <- extract_sample_names(info_full)
    info_short_std <- standardize_sample_name(info_short)
    sample_std <- standardize_sample_name(sample_short)
    idx <- match(sample_std, info_short_std)
    batch <- si$Batch[idx]
    ann_col$Batch <- batch
    message("[DEBUG] sample_cor: matched ", sum(!is.na(idx)), " samples to Batch info")
  }
  
  # 如果 ann_col 没有任何列，则设为 NULL
  if (ncol(ann_col) == 0) {
    ann_col <- NULL
  } else {
    message("[DEBUG] sample_cor: annotation columns: ", paste(colnames(ann_col), collapse = ", "), 
            ", nrow = ", nrow(ann_col))
  }
  
  # 生成注释颜色
  ann_colors <- list()
  if (!is.null(ann_col)) {
    if ("Group" %in% colnames(ann_col)) {
      groups_uniq <- unique(ann_col$Group)
      ann_colors$Group <- get_group_colors(groups_uniq)
    }
    if ("Batch" %in% colnames(ann_col)) {
      batches_uniq <- unique(na.omit(ann_col$Batch))
      ann_colors$Batch <- setNames(rainbow(length(batches_uniq)), batches_uniq)
    }
  }
  if (length(ann_colors) == 0) ann_colors <- NULL
  
  message("[DEBUG] sample_cor: correlation matrix generated, dim = ", nrow(cor_mat), " x ", ncol(cor_mat))
  list(cor_mat = cor_mat, ann_col = ann_col, ann_colors = ann_colors)
})

# ---------- 热图绘制 ----------
output$sample_cor_heatmap <- renderPlot({
  dat <- sample_cor_data()
  if (is.null(dat)) {
    plot.new()
    text(0.5, 0.5, "No data available. Please click 'Generate Correlation Heatmap'.")
    return()
  }
  
  cor_mat <- dat$cor_mat
  ann_col <- dat$ann_col
  ann_colors <- dat$ann_colors
  
  pheatmap::pheatmap(cor_mat,
                     color = colorRampPalette(c("blue", "white", "red"))(255),
                     breaks = seq(-1, 1, length.out = 256),
                     legend_breaks = c(-1, -0.5, 0, 0.5, 1),
                     legend_labels = c("-1.0", "-0.5", "0.0", "0.5", "1.0"),
                     clustering_distance_rows = as.dist(1 - cor_mat),
                     clustering_distance_cols = as.dist(1 - cor_mat),
                     clustering_method = "ward.D2",
                     show_rownames = TRUE,
                     show_colnames = TRUE,
                     fontsize_row = 9,
                     fontsize_col = 9,
                     angle_col = 45,
                     annotation_col = ann_col,
                     annotation_colors = ann_colors,
                     main = "Sample Correlation Heatmap (based on top 500 variable proteins, log2 intensity)")
})

# ---------- 下载热图 PNG ----------
output$download_sample_cor_png <- downloadHandler(
  filename = function() paste0("Sample_Correlation_", Sys.Date(), ".png"),
  content = function(file) {
    dat <- sample_cor_data()
    if (is.null(dat)) return()
    png(file, width = 900, height = 700, res = 150)
    pheatmap::pheatmap(dat$cor_mat,
                       color = colorRampPalette(c("blue", "white", "red"))(255),
                       breaks = seq(-1, 1, length.out = 256),
                       legend_breaks = c(-1, -0.5, 0, 0.5, 1),
                       legend_labels = c("-1.0", "-0.5", "0.0", "0.5", "1.0"),
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
                       main = "Sample Correlation Heatmap")
    dev.off()
  }
)

# ---------- 下载相关性矩阵 CSV ----------
output$download_sample_cor_matrix <- downloadHandler(
  filename = function() paste0("Sample_Correlation_Matrix_", Sys.Date(), ".csv"),
  content = function(file) {
    dat <- sample_cor_data()
    if (is.null(dat)) return()
    write.csv(dat$cor_mat, file, row.names = TRUE)
  }
)

message("[DEBUG] sample_correlation_plot.R loaded successfully.")