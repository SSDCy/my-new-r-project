# server/pca_plot.R
message("[DEBUG] pca_plot.R loading...")

# ---------- 辅助函数 ----------
safe_pca <- function(mat, scale = TRUE) {
  mat <- as.matrix(mat)
  mat[!is.finite(mat)] <- NA
  mat <- mat[complete.cases(mat), , drop = FALSE]
  if (nrow(mat) < 2) return(NULL)
  row_vars <- apply(mat, 1, var, na.rm = TRUE)
  keep <- row_vars > 1e-12
  mat <- mat[keep, , drop = FALSE]
  if (nrow(mat) < 2) return(NULL)
  tryCatch(prcomp(mat, scale. = scale), error = function(e) NULL)
}

# ---------- 数据准备（基于归一化数据） ----------
pca_data <- reactive({
  nd <- norm_data_full()
  if (is.null(nd)) {
    message("[DEBUG] pca_plot: norm_data_full() is NULL")
    return(NULL)
  }
  
  # 提取 Norm_ 开头的强度列
  norm_cols <- grep("^Norm_LFQ intensity ", colnames(nd), value = TRUE)
  if (length(norm_cols) == 0) {
    message("[DEBUG] pca_plot: no Norm_ columns found")
    return(NULL)
  }
  
  expr_mat <- as.matrix(nd[, norm_cols])
  rownames(expr_mat) <- nd$`Master protein IDs`
  
  # 填补缺失值（如果存在）——使用 1% 分位数填补
  if (any(is.na(expr_mat))) {
    message("[DEBUG] pca_plot: imputing missing values with Quantile (1%)")
    expr_mat <- as.matrix(impute_missing_values(
      as.data.frame(expr_mat), method = "quantile", quantile_prob = 0.01
    ))
  }
  
  # log2 转换
  log_expr <- log2(expr_mat + 1)
  
  # 过滤低变异蛋白
  row_vars <- apply(log_expr, 1, var)
  keep <- row_vars > 1e-6
  log_expr <- log_expr[keep, , drop = FALSE]
  message("[DEBUG] pca_plot: after filtering low variance, ", nrow(log_expr), " proteins remain")
  
  if (nrow(log_expr) < 2) {
    message("[DEBUG] pca_plot: too few variable proteins")
    return(NULL)
  }
  
  # 转置：样本为行，蛋白为列
  mat_t <- t(log_expr)
  
  # PCA
  pca <- safe_pca(mat_t)
  if (is.null(pca)) {
    message("[DEBUG] pca_plot: PCA failed")
    return(NULL)
  }
  
  variance <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 1)
  scores <- as.data.frame(pca$x[, 1:2])
  scores$Sample <- rownames(scores)
  
  # 获取短样本名
  sample_short <- gsub("^Norm_LFQ intensity ", "", scores$Sample)
  
  # 附加分组信息
  groups <- rep("Unassigned", nrow(scores))
  batch <- rep(NA, nrow(scores))
  if (!is.null(rv$groups) && length(rv$groups) > 0) {
    for (i in seq_along(sample_short)) {
      sn <- sample_short[i]
      for (gn in names(rv$groups)) {
        if (sn %in% rv$groups[[gn]]) {
          groups[i] <- gn
          break
        }
      }
    }
  }
  
  # 修正 Batch 匹配：将样本信息表的行名转换为短名再匹配
  if (!is.null(rv$sample_info) && "Batch" %in% colnames(rv$sample_info)) {
    si <- rv$sample_info
    # 提取样本信息表的短名：去掉 "LFQ intensity " 前缀，再标准化
    info_full <- rownames(si)
    info_short <- extract_sample_names(info_full)   # 得到类似 "WT.1", "100.12.1" 等
    info_short_std <- standardize_sample_name(info_short)
    sample_std <- standardize_sample_name(sample_short)
    idx <- match(sample_std, info_short_std)
    batch <- si$Batch[idx]
    message("[DEBUG] pca_plot: matched ", sum(!is.na(idx)), " samples to Batch info, ",
            sum(is.na(idx)), " unmatched")
  } else {
    message("[DEBUG] pca_plot: no Batch info available in sample_info")
  }
  
  scores$Group <- groups
  scores$Batch <- batch
  scores$SampleShort <- sample_short
  
  # 离群样本检测
  z1 <- abs((scores$PC1 - mean(scores$PC1)) / sd(scores$PC1))
  z2 <- abs((scores$PC2 - mean(scores$PC2)) / sd(scores$PC2))
  scores$Outlier <- ifelse(z1 > 3 | z2 > 3, "Outlier", "Normal")
  
  message("[DEBUG] pca_plot: PCA completed. PC1=", variance[1], "%, PC2=", variance[2],
          "%, outliers=", sum(scores$Outlier == "Outlier"))
  
  list(scores = scores, variance = variance)
})

# ---------- 数据源提示 ----------
output$pca_data_source_note <- renderUI({
  if (is.null(norm_data_full())) {
    return(div(style = "color: #e74c3c;", "Normalized data not available. Please run preprocessing."))
  }
  div(style = "margin-bottom: 10px; color: #27ae60; font-weight: bold;",
      icon("check-circle"), " Data source: Normalized expression data (total intensity normalization)")
})

# ---------- PCA by Group ----------
output$pca_group_plot <- renderPlotly({
  dat <- pca_data()
  if (is.null(dat)) return(plotly::plot_ly() %>% layout(title = "No data"))
  scores <- dat$scores
  var <- dat$variance
  
  # 颜色配置
  groups <- unique(scores$Group)
  group_colors <- get_group_colors(groups)
  
  # 形状区分离群样本
  scores$shape <- ifelse(scores$Outlier == "Outlier", "circle-open", "circle")
  
  p <- plot_ly(
    data = scores,
    x = ~PC1, y = ~PC2,
    color = ~Group, colors = group_colors,
    symbol = ~Outlier, symbols = c("circle", "circle-open"),
    type = "scatter", mode = "markers",
    marker = list(size = 10),
    text = ~paste("Sample:", SampleShort, "<br>Group:", Group, "<br>Outlier:", Outlier),
    hoverinfo = "text"
  ) %>%
    layout(
      title = "PCA by Group",
      xaxis = list(title = paste0("PC1 (", var[1], "%)")),
      yaxis = list(title = paste0("PC2 (", var[2], "%)")),
      legend = list(title = list(text = "Group"))
    )
  p
})

# ---------- PCA by Batch ----------
output$pca_batch_plot <- renderPlotly({
  dat <- pca_data()
  if (is.null(dat)) return(plotly::plot_ly() %>% layout(title = "No data"))
  scores <- dat$scores
  var <- dat$variance
  
  # 检查是否有 Batch 信息
  if (all(is.na(scores$Batch))) {
    return(plotly::plot_ly() %>% layout(title = "No Batch information available"))
  }
  
  batches <- unique(na.omit(scores$Batch))
  batch_colors <- setNames(rainbow(length(batches)), batches)
  
  scores$shape <- ifelse(scores$Outlier == "Outlier", "circle-open", "circle")
  
  p <- plot_ly(
    data = scores,
    x = ~PC1, y = ~PC2,
    color = ~Batch, colors = batch_colors,
    symbol = ~Outlier, symbols = c("circle", "circle-open"),
    type = "scatter", mode = "markers",
    marker = list(size = 10),
    text = ~paste("Sample:", SampleShort, "<br>Batch:", Batch, "<br>Outlier:", Outlier),
    hoverinfo = "text"
  ) %>%
    layout(
      title = "PCA by Batch",
      xaxis = list(title = paste0("PC1 (", var[1], "%)")),
      yaxis = list(title = paste0("PC2 (", var[2], "%)")),
      legend = list(title = list(text = "Batch"))
    )
  p
})

# ---------- 离群样本信息 ----------
output$pca_outlier_info <- renderPrint({
  dat <- pca_data()
  if (is.null(dat)) {
    cat("PCA not available.\n")
    return()
  }
  outliers <- dat$scores[dat$scores$Outlier == "Outlier", "SampleShort"]
  if (length(outliers) == 0) {
    cat("No outlier samples detected (Z-score > 3 on PC1 or PC2).\n")
  } else {
    cat("Potential outlier samples (Z-score > 3):\n")
    cat(paste(outliers, collapse = ", "), "\n")
    cat("These samples may have abnormal expression profiles and could be investigated further.\n")
  }
})

# ---------- 下载 PCA 图 ----------
output$download_pca_group_png <- downloadHandler(
  filename = function() paste0("PCA_Group_", Sys.Date(), ".png"),
  content = function(file) {
    dat <- pca_data()
    if (is.null(dat)) return()
    scores <- dat$scores
    var <- dat$variance
    groups <- unique(scores$Group)
    group_colors <- get_group_colors(groups)
    p <- ggplot(scores, aes(x = PC1, y = PC2, color = Group, shape = Outlier)) +
      geom_point(size = 3) +
      scale_color_manual(values = group_colors) +
      scale_shape_manual(values = c("Normal" = 16, "Outlier" = 1)) +
      labs(title = "PCA by Group", x = paste0("PC1 (", var[1], "%)"), y = paste0("PC2 (", var[2], "%)")) +
      theme_bw()
    ggsave(file, plot = p, width = 8, height = 6, dpi = 150)
  }
)

output$download_pca_batch_png <- downloadHandler(
  filename = function() paste0("PCA_Batch_", Sys.Date(), ".png"),
  content = function(file) {
    dat <- pca_data()
    if (is.null(dat)) return()
    scores <- dat$scores
    var <- dat$variance
    if (all(is.na(scores$Batch))) return()
    p <- ggplot(scores, aes(x = PC1, y = PC2, color = Batch, shape = Outlier)) +
      geom_point(size = 3) +
      scale_shape_manual(values = c("Normal" = 16, "Outlier" = 1)) +
      labs(title = "PCA by Batch", x = paste0("PC1 (", var[1], "%)"), y = paste0("PC2 (", var[2], "%)")) +
      theme_bw()
    ggsave(file, plot = p, width = 8, height = 6, dpi = 150)
  }
)

# ---------- 步骤指示器（可折叠） ----------
output$pca_preprocess_steps <- renderUI({
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
  steps <- c(steps, "Missing values imputed with Quantile 1% (if any)")
  
  step_tags <- lapply(seq_along(steps), function(i) {
    tagList(
      if (i > 1) tags$span(style = "font-size: 20px; color: #27ae60; margin: 0 8px;", "→"),
      tags$span(style = "background: #e8f0fe; padding: 6px 12px; border-radius: 15px; font-size: 13px;", steps[[i]])
    )
  })
  tags$details(
    tags$summary("Data preprocessing steps for PCA", style = "cursor: pointer; font-weight: bold; color: #2c3e50; margin-bottom: 10px;"),
    div(style = "display: flex; flex-wrap: wrap; align-items: center;", do.call(tagList, step_tags))
  )
})

message("[DEBUG] pca_plot.R loaded successfully.")