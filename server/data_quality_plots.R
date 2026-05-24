# server/data_quality_plots.R

# ==================== 辅助函数 ====================
get_outlier_samples <- function(pca_result, z_threshold = 3) {
  if (is.null(pca_result)) return(character(0))
  scores <- pca_result$pca_df[, c("PC1", "PC2")]
  z1 <- abs((scores$PC1 - mean(scores$PC1)) / sd(scores$PC1))
  z2 <- abs((scores$PC2 - mean(scores$PC2)) / sd(scores$PC2))
  outlier_mask <- z1 > z_threshold | z2 > z_threshold
  pca_result$pca_df$Sample[outlier_mask]
}

# ==================== 数据质量评分 ====================
output$dq_score <- renderText({
  req(dq_expr_matrix())
  tryCatch(dq_quality_score()$score, error = function(e) paste("Error:", e$message))
})
output$dq_grade <- renderText({
  req(dq_expr_matrix())
  tryCatch(dq_quality_score()$grade, error = function(e) paste("Error:", e$message))
})
output$dq_missing_rate <- renderText({
  req(dq_expr_matrix())
  tryCatch(paste0(dq_quality_score()$details$missing_ratio, "%"), error = function(e) paste("Error:", e$message))
})
output$dq_missing_score_frac <- renderText({
  req(dq_expr_matrix())
  tryCatch(paste0(dq_quality_score()$details$missing_score, " / 30"), error = function(e) paste("Error:", e$message))
})
output$dq_consistency_score_frac <- renderText({
  req(dq_expr_matrix())
  tryCatch(paste0(dq_quality_score()$details$correlation_score, " / 20"), error = function(e) paste("Error:", e$message))
})
output$dq_protein_score_frac <- renderText({
  req(dq_expr_matrix())
  tryCatch(paste0(dq_quality_score()$details$protein_score, " / 20"), error = function(e) paste("Error:", e$message))
})
output$dq_missing_grade <- renderText({
  req(dq_expr_matrix())
  tryCatch({
    s <- dq_quality_score()$details$missing_score
    if (s >= 25) "Good" else if (s >= 15) "Fair" else "Poor"
  }, error = function(e) paste("Error:", e$message))
})
output$dq_consistency_grade <- renderText({
  req(dq_expr_matrix())
  tryCatch({
    if (dq_quality_score()$details$avg_correlation > 0.8) "Good"
    else if (dq_quality_score()$details$avg_correlation > 0.7) "Fair" else "Poor"
  }, error = function(e) paste("Error:", e$message))
})
output$dq_protein_grade <- renderText({
  req(dq_expr_matrix())
  tryCatch({
    if (dq_quality_score()$details$protein_valid_ratio > 80) "Good"
    else if (dq_quality_score()$details$protein_valid_ratio > 60) "Fair" else "Poor"
  }, error = function(e) paste("Error:", e$message))
})

# 保留旧版输出（兼容性）
output$dq_missing_score <- renderText({ req(dq_expr_matrix()); tryCatch(paste0("Score: ", dq_quality_score()$details$missing_score, "/30"), error = function(e) paste("Error:", e$message)) })
output$dq_consistency <- renderText({ req(dq_expr_matrix()); tryCatch({ if (dq_quality_score()$details$avg_correlation > 0.8) "Good" else if (dq_quality_score()$details$avg_correlation > 0.7) "Fair" else "Poor" }, error = function(e) paste("Error:", e$message)) })
output$dq_consistency_score <- renderText({ req(dq_expr_matrix()); tryCatch(paste0("Score: ", dq_quality_score()$details$correlation_score, "/20"), error = function(e) paste("Error:", e$message)) })
output$dq_protein_quality <- renderText({ req(dq_expr_matrix()); tryCatch({ if (dq_quality_score()$details$protein_valid_ratio > 80) "Good" else if (dq_quality_score()$details$protein_valid_ratio > 60) "Fair" else "Poor" }, error = function(e) paste("Error:", e$message)) })
output$dq_protein_score <- renderText({ req(dq_expr_matrix()); tryCatch(paste0("Score: ", dq_quality_score()$details$protein_score, "/20"), error = function(e) paste("Error:", e$message)) })

# ==================== 智能报告 ====================
output$dq_total_score <- renderText({ req(dq_expr_matrix()); tryCatch(dq_quality_score()$score, error = function(e) paste("Error:", e$message)) })
output$dq_total_grade <- renderText({ req(dq_expr_matrix()); tryCatch(dq_quality_score()$grade, error = function(e) paste("Error:", e$message)) })
output$dq_total_missing <- renderText({ req(dq_expr_matrix()); tryCatch(paste0(dq_quality_score()$details$missing_ratio, "%"), error = function(e) paste("Error:", e$message)) })
output$dq_total_consistency <- renderText({ req(dq_expr_matrix()); tryCatch({ if (dq_quality_score()$details$avg_correlation > 0.8) "Good" else if (dq_quality_score()$details$avg_correlation > 0.7) "Fair" else "Poor" }, error = function(e) paste("Error:", e$message)) })
output$dq_total_correlation <- renderText({ req(dq_expr_matrix()); tryCatch(dq_quality_score()$details$avg_correlation, error = function(e) paste("Error:", e$message)) })

output$dq_key_findings <- renderUI({
  req(dq_expr_matrix())
  tryCatch({
    req(dq_quality_score())
    report <- generate_quality_report(dq_quality_score(), dq_expr_matrix(), rv$sample_info)
    tagList(lapply(report$key_findings, render_key_finding))
  }, error = function(e) div(style = "color: red;", paste("Error:", e$message)))
})
output$dq_recommendations <- renderUI({
  req(dq_expr_matrix())
  tryCatch({
    req(dq_quality_score())
    report <- generate_quality_report(dq_quality_score(), dq_expr_matrix(), rv$sample_info)
    tagList(lapply(report$recommendations, render_recommendation))
  }, error = function(e) div(style = "color: red;", paste("Error:", e$message)))
})
output$dq_special_note <- renderUI({
  req(dq_expr_matrix())
  tryCatch({
    req(dq_quality_score())
    report <- generate_quality_report(dq_quality_score(), dq_expr_matrix(), rv$sample_info)
    if (report$special_note != "") {
      div(style = "background: #e3f2fd; padding: 15px; border-radius: 8px; margin-top: 15px;",
          h4(icon("info-circle"), " 特别说明", style = "margin: 0 0 10px 0; font-size: 16px; color: #1976d2;"),
          p(style = "margin: 0;", report$special_note))
    } else NULL
  }, error = function(e) div(style = "color: red;", paste("Error:", e$message)))
})

# ==================== 缺失值热图 ====================
dq_missing_heatmap_plot <- reactive({
  req(dq_expr_matrix())
  mat <- dq_expr_matrix()
  missing_mat <- is.na(mat) * 1
  if (nrow(missing_mat) > 1000) {
    set.seed(123)
    missing_mat <- missing_mat[sample(1:nrow(missing_mat), 1000), ]
  }
  df <- reshape2::melt(missing_mat)
  colnames(df) <- c("Protein", "Sample", "Missing")
  ggplot(df, aes(x = Sample, y = Protein, fill = factor(Missing))) +
    geom_tile() +
    scale_fill_manual(values = c("0" = "#3498db", "1" = "#e74c3c"),
                      labels = c("0" = "Detected", "1" = "Missing"), name = "Status") +
    labs(title = "Missing Value Heatmap (Blue = Detected, Red = Missing, Rows = Proteins, Columns = Samples)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none")
})

output$dq_missing_heatmap <- renderPlot({ dq_missing_heatmap_plot() })
output$download_missing_heatmap <- downloadHandler(
  filename = function() "missing_heatmap.png",
  content = function(file) ggsave(file, plot = dq_missing_heatmap_plot(), width = 10, height = 6, dpi = 150)
)
observeEvent(input$help_missing_heatmap, {
  showModal(modalDialog(
    title = "Missing Value Heatmap",
    "蓝色表示该蛋白在对应样本中被检测到，红色表示缺失。行代表蛋白，列代表样本。通过聚类可观察缺失模式是否与样本分组相关。",
    easyClose = TRUE, footer = modalButton("关闭")
  ))
})

# ==================== 有效值柱状图 ====================
dq_valid_plot <- reactive({
  req(dq_expr_matrix(), dq_missing_stats())
  stats <- dq_missing_stats()
  valid_percent <- (1 - stats$sample_missing) * 100
  df <- data.frame(Sample = factor(names(valid_percent), levels = names(valid_percent)), ValidPercent = valid_percent)
  threshold <- 70
  ggplot(df, aes(x = Sample, y = ValidPercent)) +
    geom_col(aes(fill = ValidPercent < threshold)) +
    geom_hline(yintercept = threshold, color = "red", linetype = "dashed", linewidth = 1) +
    scale_fill_manual(values = c("TRUE" = "#e74c3c", "FALSE" = "#2ecc71"), guide = "none") +
    labs(title = "Valid Values per Sample (Red dashed line = 70% quality threshold)", y = "Valid Values (%)") +
    ylim(0, 80) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
})

output$dq_valid_values_plot <- renderPlot({ dq_valid_plot() })
output$download_valid_bar <- downloadHandler(
  filename = function() "valid_values.png",
  content = function(file) ggsave(file, plot = dq_valid_plot(), width = 8, height = 5, dpi = 150)
)
observeEvent(input$help_valid_bar, {
  showModal(modalDialog(
    title = "Valid Values per Sample",
    "每个样本中有效（非缺失）蛋白的百分比。红色参考线为 70% 阈值，低于此线的样本可能需要进一步检查或移除。",
    easyClose = TRUE, footer = modalButton("关闭")
  ))
})

# ==================== 缺失值相关性热图 ====================
dq_missing_cor_plot_obj <- reactive({
  req(dq_expr_matrix())
  mat <- dq_expr_matrix()
  missing_mat <- is.na(mat) * 1
  if (ncol(missing_mat) < 2) return(NULL)
  cor_mat <- cor(missing_mat, use = "pairwise.complete.obs")
  cor_mat[is.na(cor_mat)] <- 0
  pheatmap::pheatmap(cor_mat, 
                     main = "Missing Value Correlation (Red = High correlation, Blue = Low correlation)",
                     color = colorRampPalette(c("blue", "white", "red"))(100),
                     show_rownames = TRUE, show_colnames = TRUE,
                     fontsize_row = 8, fontsize_col = 8,
                     angle_col = 45)
})

output$dq_missing_cor_plot <- renderPlot({
  obj <- dq_missing_cor_plot_obj()
  if (!is.null(obj)) {
    grid::grid.newpage()
    grid::grid.draw(obj$gtable)
  } else {
    plot.new(); text(0.5, 0.5, "Not enough samples")
  }
})
output$download_missing_cor <- downloadHandler(
  filename = function() "missing_correlation.png",
  content = function(file) {
    png(file, width = 800, height = 600, res = 150)
    obj <- dq_missing_cor_plot_obj()
    if (!is.null(obj)) grid::grid.draw(obj$gtable)
    dev.off()
  }
)
observeEvent(input$help_missing_cor, {
  showModal(modalDialog(
    title = "Missing Value Correlation",
    "展示样本间缺失模式的相关性。如果某些样本的缺失模式高度相关，提示可能存在批次效应或技术偏差。",
    easyClose = TRUE, footer = modalButton("关闭")
  ))
})

# ==================== 强度分布箱线图 ====================
dq_intensity_plot <- reactive({
  req(dq_expr_matrix())
  mat <- dq_expr_matrix()
  log_mat <- log2(mat + 1)
  df <- reshape2::melt(as.matrix(log_mat))
  colnames(df) <- c("Protein", "Sample", "Log2Intensity")
  ggplot(df, aes(x = Sample, y = Log2Intensity)) +
    geom_boxplot(fill = "#3498db", alpha = 0.7, outlier.size = 1) +
    labs(title = "Protein Intensity Distribution (log2-transformed)", y = "log2(Intensity)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
})

output$dq_intensity_dist_plot <- renderPlot({ dq_intensity_plot() })
output$download_intensity <- downloadHandler(
  filename = function() "intensity_distribution.png",
  content = function(file) ggsave(file, plot = dq_intensity_plot(), width = 10, height = 6, dpi = 150)
)
observeEvent(input$help_intensity, {
  showModal(modalDialog(
    title = "Protein Intensity Distribution",
    "箱线图展示每个样本中蛋白强度的分布（log2 转换）。异常样本通常表现出整体偏移或异常离散。",
    easyClose = TRUE, footer = modalButton("关闭")
  ))
})

# ==================== 样本相关性热图（标签旋转 + 调试） ====================
dq_cor_heatmap_plot_obj <- reactive({
  message("[DEBUG] dq_cor_heatmap_plot_obj: entering")
  req(dq_expr_matrix())
  message("[DEBUG] dq_cor_heatmap_plot_obj: dq_expr_matrix obtained")
  cor_mat <- calculate_sample_correlation(dq_expr_matrix())
  message("[DEBUG] dq_cor_heatmap_plot_obj: correlation matrix dimensions: ", if(is.null(cor_mat)) "NULL" else paste(nrow(cor_mat), "x", ncol(cor_mat)))
  if (is.null(cor_mat)) {
    message("[DEBUG] dq_cor_heatmap_plot_obj: cor_mat is NULL, returning NULL")
    return(NULL)
  }
  ann_col <- NULL; ann_colors <- NULL
  if (!is.null(rv$sample_info) && "Group" %in% colnames(rv$sample_info)) {
    message("[DEBUG] dq_cor_heatmap_plot_obj: sample_info found, extracting Group")
    sample_info_short <- rv$sample_info
    rownames(sample_info_short) <- extract_sample_names(rownames(sample_info_short))
    common_samples <- intersect(colnames(cor_mat), rownames(sample_info_short))
    message("[DEBUG] dq_cor_heatmap_plot_obj: common samples: ", length(common_samples))
    if (length(common_samples) > 0) {
      ann_col <- data.frame(Group = sample_info_short[common_samples, "Group"], row.names = common_samples)
      ann_col <- ann_col[!is.na(ann_col$Group) & ann_col$Group != "a", , drop = FALSE]
      if (nrow(ann_col) > 0) {
        cor_mat <- cor_mat[rownames(ann_col), rownames(ann_col)]
        groups <- unique(ann_col$Group)
        ann_colors <- list(Group = get_group_colors(groups))
        message("[DEBUG] dq_cor_heatmap_plot_obj: annotation prepared, groups: ", paste(groups, collapse = ", "))
      } else {
        message("[DEBUG] dq_cor_heatmap_plot_obj: ann_col is empty after filtering")
      }
    }
  } else {
    message("[DEBUG] dq_cor_heatmap_plot_obj: no sample_info or no Group column")
  }
  message("[DEBUG] dq_cor_heatmap_plot_obj: drawing pheatmap...")
  p <- pheatmap::pheatmap(cor_mat, 
                          main = "Sample Correlation Heatmap (Red = High correlation, Blue = Low correlation)",
                          color = colorRampPalette(c("blue", "white", "red"))(100),
                          breaks = seq(0.4, 1, length.out = 101),
                          show_rownames = TRUE, show_colnames = TRUE,
                          fontsize_row = 9, fontsize_col = 9,
                          angle_col = 45,
                          annotation_col = ann_col, annotation_colors = ann_colors,
                          silent = TRUE)
  message("[DEBUG] dq_cor_heatmap_plot_obj: pheatmap returned")
  return(p)
})

output$dq_cor_heatmap <- renderPlot({
  message("[DEBUG] dq_cor_heatmap render: called")
  obj <- dq_cor_heatmap_plot_obj()
  if (is.null(obj)) {
    plot.new(); text(0.5, 0.5, "Not enough data")
    message("[DEBUG] dq_cor_heatmap: obj is NULL, plotting text")
  } else {
    message("[DEBUG] dq_cor_heatmap: drawing object")
    grid::grid.newpage()
    grid::grid.draw(obj$gtable)
    message("[DEBUG] dq_cor_heatmap: drawing finished")
  }
})

output$download_cor_heatmap <- downloadHandler(
  filename = function() "sample_correlation.png",
  content = function(file) {
    png(file, width = 900, height = 700, res = 150)
    obj <- dq_cor_heatmap_plot_obj()
    if (!is.null(obj)) grid::grid.draw(obj$gtable)
    dev.off()
  }
)
observeEvent(input$help_cor_heatmap, {
  showModal(modalDialog(
    title = "Sample Correlation Heatmap",
    "基于高变异蛋白计算的样本间相关性。颜色越红表示相关性越高。理想情况下，同组样本应聚在一起。",
    easyClose = TRUE, footer = modalButton("关闭")
  ))
})

# ==================== PCA 图（双维度，图例精简） ====================
pca_group_plot_obj <- reactive({
  req(dq_expr_matrix())
  pca_result <- calculate_pca(dq_expr_matrix(), rv$sample_info)
  if (is.null(pca_result)) return(NULL)
  pca_df <- pca_result$pca_df
  if (!is.null(rv$sample_info) && "Group" %in% colnames(rv$sample_info)) {
    sample_info_short <- rv$sample_info
    rownames(sample_info_short) <- extract_sample_names(rownames(sample_info_short))
    pca_df$Group <- sample_info_short[extract_sample_names(pca_df$Sample), "Group"]
    pca_df <- pca_df[!is.na(pca_df$Group) & pca_df$Group != "a", ]
  } else {
    pca_df <- pca_df[!is.na(pca_df$Group) & pca_df$Group != "a", ]
  }
  outliers <- get_outlier_samples(pca_result)
  pca_df$Outlier <- ifelse(pca_df$Sample %in% outliers, "Outlier", "Normal")
  group_colors <- c("Control" = "#FF69B4", "Treatment" = "#00CED1")
  all_groups <- unique(pca_df$Group)
  missing_colors <- setdiff(all_groups, names(group_colors))
  if (length(missing_colors) > 0) {
    extra_colors <- rainbow(length(missing_colors))
    names(extra_colors) <- missing_colors
    group_colors <- c(group_colors, extra_colors)
  }
  # 图例中去掉 "(red)" 冗余描述，只保留 Normal/Outlier
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = Outlier)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_text(data = pca_df[pca_df$Outlier == "Outlier", ], aes(label = Sample), vjust = 1.5, size = 3, show.legend = FALSE) +
    scale_shape_manual(values = c("Normal" = 16, "Outlier" = 17),
                       labels = c("Normal" = "Normal", "Outlier" = "Outlier")) +
    scale_color_manual(values = group_colors) +
    labs(title = "PCA by Group (Outliers in red ▲)",
         x = paste0("PC1 (", pca_result$pc1_var, "%)"),
         y = paste0("PC2 (", pca_result$pc2_var, "%)")) +
    theme_bw() + theme(legend.position = "right",
                       axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
})

output$dq_pca_group_plot <- renderPlot({ req(pca_group_plot_obj()); pca_group_plot_obj() })
output$download_pca_group <- downloadHandler(
  filename = function() "pca_group.png",
  content = function(file) ggsave(file, plot = pca_group_plot_obj(), width = 8, height = 6, dpi = 150)
)
observeEvent(input$help_pca_group, {
  showModal(modalDialog(
    title = "PCA by Group",
    "按实验分组着色，用于观察组间分离程度。若组间分离明显，说明生物学差异是主要变异来源。异常样本已用红色三角形标注。",
    easyClose = TRUE, footer = modalButton("关闭")
  ))
})

pca_batch_plot_obj <- reactive({
  req(dq_expr_matrix())
  if (is.null(rv$sample_info) || !"Batch" %in% colnames(rv$sample_info)) return(NULL)
  pca_result <- calculate_pca(dq_expr_matrix(), rv$sample_info)
  if (is.null(pca_result)) return(NULL)
  pca_df <- pca_result$pca_df
  sample_info_short <- rv$sample_info
  rownames(sample_info_short) <- extract_sample_names(rownames(sample_info_short))
  pca_df$Batch <- sample_info_short[extract_sample_names(pca_df$Sample), "Batch"]
  pca_df <- pca_df[!is.na(pca_df$Batch) & pca_df$Batch != "a", ]
  outliers <- get_outlier_samples(pca_result)
  pca_df$Outlier <- ifelse(pca_df$Sample %in% outliers, "Outlier", "Normal")
  batch_colors <- c("Batch1" = "#E41A1C", "Batch2" = "#00CED1")
  all_batches <- unique(pca_df$Batch)
  if (length(all_batches) > 2) {
    batch_colors <- setNames(rainbow(length(all_batches)), all_batches)
  }
  # 图例精简
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch, shape = Outlier)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_text(data = pca_df[pca_df$Outlier == "Outlier", ], aes(label = Sample), vjust = 1.5, size = 3, show.legend = FALSE) +
    scale_shape_manual(values = c("Normal" = 16, "Outlier" = 17),
                       labels = c("Normal" = "Normal", "Outlier" = "Outlier")) +
    scale_color_manual(values = batch_colors) +
    labs(title = "PCA by Batch (Outliers in red ▲)",
         x = paste0("PC1 (", pca_result$pc1_var, "%)"),
         y = paste0("PC2 (", pca_result$pc2_var, "%)")) +
    theme_bw() + theme(legend.position = "right",
                       axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
})

output$dq_pca_batch_plot <- renderPlot({
  if (is.null(pca_batch_plot_obj())) {
    plot.new(); text(0.5, 0.5, "Batch information not available")
  } else {
    pca_batch_plot_obj()
  }
})
output$download_pca_batch <- downloadHandler(
  filename = function() "pca_batch.png",
  content = function(file) {
    p <- pca_batch_plot_obj()
    if (!is.null(p)) ggsave(file, plot = p, width = 8, height = 6, dpi = 150)
  }
)
observeEvent(input$help_pca_batch, {
  showModal(modalDialog(
    title = "PCA by Batch",
    "按实验批次着色，用于检测批次效应。若样本按批次聚集，说明批次是主要变异来源，需进行批次校正。异常样本已用红色三角形标注。",
    easyClose = TRUE, footer = modalButton("关闭")
  ))
})

# ==================== 预处理跳转按钮 ====================
observeEvent(input$goto_preprocessing, {
  message("Jumping to Data Preprocessing page...")
  updateNavbarPage(session, "main_navbar", selected = "preprocessing")
})