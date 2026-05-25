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
# 显示为 "分数/30"
output$dq_missing_score_frac <- renderText({
  req(dq_expr_matrix())
  tryCatch(paste0(dq_quality_score()$details$missing_score, " / 30"), error = function(e) paste("Error:", e$message))
})
# 显示为 "分数/40"
output$dq_consistency_score_frac <- renderText({
  req(dq_expr_matrix())
  tryCatch(paste0(dq_quality_score()$details$consistency_score, " / 40"), error = function(e) paste("Error:", e$message))
})
# 显示为 "分数/30"
output$dq_protein_score_frac <- renderText({
  req(dq_expr_matrix())
  tryCatch(paste0(dq_quality_score()$details$protein_score, " / 30"), error = function(e) paste("Error:", e$message))
})

# 分项评级：使用得分率映射到统一标准
get_grade_from_ratio <- function(score, max_score) {
  ratio <- score / max_score
  if (ratio >= 0.9) return("Excellent")
  if (ratio >= 0.8) return("Good")
  if (ratio >= 0.6) return("Fair")
  return("Poor")
}

output$dq_missing_grade <- renderText({
  req(dq_expr_matrix())
  tryCatch({
    s <- dq_quality_score()$details$missing_score
    get_grade_from_ratio(s, 30)
  }, error = function(e) paste("Error:", e$message))
})
output$dq_consistency_grade <- renderText({
  req(dq_expr_matrix())
  tryCatch({
    s <- dq_quality_score()$details$consistency_score
    get_grade_from_ratio(s, 40)
  }, error = function(e) paste("Error:", e$message))
})
output$dq_protein_grade <- renderText({
  req(dq_expr_matrix())
  tryCatch({
    s <- dq_quality_score()$details$protein_score
    get_grade_from_ratio(s, 30)
  }, error = function(e) paste("Error:", e$message))
})

# 保留旧版输出（兼容性）
output$dq_missing_score <- renderText({ req(dq_expr_matrix()); tryCatch(paste0("Score: ", dq_quality_score()$details$missing_score, "/30"), error = function(e) paste("Error:", e$message)) })
output$dq_consistency <- renderText({ req(dq_expr_matrix()); tryCatch({ if (dq_quality_score()$details$avg_correlation > 0.8) "Good" else if (dq_quality_score()$details$avg_correlation > 0.7) "Fair" else "Poor" }, error = function(e) paste("Error:", e$message)) })
output$dq_consistency_score <- renderText({ req(dq_expr_matrix()); tryCatch(paste0("Score: ", dq_quality_score()$details$consistency_score, "/40"), error = function(e) paste("Error:", e$message)) })
output$dq_protein_quality <- renderText({ req(dq_expr_matrix()); tryCatch({ if (dq_quality_score()$details$protein_valid_ratio > 80) "Good" else if (dq_quality_score()$details$protein_valid_ratio > 60) "Fair" else "Poor" }, error = function(e) paste("Error:", e$message)) })
output$dq_protein_score <- renderText({ req(dq_expr_matrix()); tryCatch(paste0("Score: ", dq_quality_score()$details$protein_score, "/30"), error = function(e) paste("Error:", e$message)) })

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

# 导出完整缺失矩阵 (Detected/Missing)
output$download_missing_matrix <- downloadHandler(
  filename = function() { paste0("missing_matrix_", Sys.Date(), ".csv") },
  content = function(file) {
    message("[DEBUG] download_missing_matrix started")
    req(dq_expr_matrix())
    mat <- dq_expr_matrix()
    miss_df <- as.data.frame(is.na(mat))
    miss_df[] <- lapply(miss_df, function(x) ifelse(x, "Missing", "Detected"))
    rownames(miss_df) <- rownames(mat)
    cat("# Missing value matrix. Each cell indicates Detected or Missing. Rows: proteins, Columns: samples.\n", file = file)
    write.csv(miss_df, file, append = TRUE, row.names = TRUE)
    message("[DEBUG] download_missing_matrix finished")
  }
)

# 导出样本缺失率统计表（增加 MissingLevel 和 KeepStatus）
output$download_sample_missing_stats <- downloadHandler(
  filename = function() { paste0("sample_missing_stats_", Sys.Date(), ".csv") },
  content = function(file) {
    message("[DEBUG] download_sample_missing_stats started")
    req(dq_expr_matrix())
    mat <- dq_expr_matrix()
    total_proteins <- nrow(mat)
    missing_per_sample <- colSums(is.na(mat))
    valid_per_sample <- total_proteins - missing_per_sample
    valid_pct <- round(valid_per_sample / total_proteins * 100, 2)
    missing_pct <- 100 - valid_pct
    
    MissingLevel <- ifelse(missing_pct < 30, "Low",
                           ifelse(missing_pct < 50, "Medium", "High"))
    KeepStatus <- ifelse(valid_pct >= 70, "Keep", "Filter")
    
    stats <- data.frame(
      SampleID = colnames(mat),
      TotalProteins = total_proteins,
      ValidProteins = valid_per_sample,
      MissingProteins = missing_per_sample,
      ValidPercentage = valid_pct,
      MissingPercentage = missing_pct,
      MissingLevel = MissingLevel,
      KeepStatus = KeepStatus
    )
    cat("# Sample missing statistics. ValidPercentage = ValidProteins / TotalProteins * 100. KeepStatus: Keep if ValidPercentage >= 70% (recommended threshold), else Filter.\n", file = file)
    write.csv(stats, file, append = TRUE, row.names = FALSE)
    message("[DEBUG] download_sample_missing_stats finished")
  }
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

# 导出缺失模式相关系数矩阵（带注释）
output$download_missing_cor_matrix <- downloadHandler(
  filename = function() { paste0("missing_correlation_matrix_", Sys.Date(), ".csv") },
  content = function(file) {
    message("[DEBUG] download_missing_cor_matrix started")
    req(dq_expr_matrix())
    mat <- dq_expr_matrix()
    missing_mat <- is.na(mat) * 1
    if (ncol(missing_mat) < 2) {
      write.csv(data.frame(Note = "Not enough samples"), file, row.names = FALSE)
      return()
    }
    cor_mat <- cor(missing_mat, use = "pairwise.complete.obs")
    cor_mat[is.na(cor_mat)] <- 0
    cat("# Missing pattern Pearson correlation matrix, values range 0-1. Higher values indicate more similar missing profiles between samples, suggesting potential systematic missing bias.\n", file = file)
    write.csv(cor_mat, file, append = TRUE, row.names = TRUE)
    message("[DEBUG] download_missing_cor_matrix finished")
  }
)

observeEvent(input$help_missing_cor, {
  showModal(modalDialog(
    title = "Missing Value Correlation",
    "展示样本间缺失模式的相关性。如果某些样本的缺失模式高度相关，提示可能存在批次效应或技术偏差。",
    easyClose = TRUE, footer = modalButton("关闭")
  ))
})

# ==================== 强度分布箱线图（图片 + Excel数据） ====================
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

# 强度数据 Excel（增加 IsOutlier 并重命名统计列）
output$download_intensity_data <- downloadHandler(
  filename = function() { paste0("Intensity_Data_", Sys.Date(), ".xlsx") },
  content = function(file) {
    message("[DEBUG] download_intensity_data started")
    req(dq_expr_matrix())
    mat <- dq_expr_matrix()
    log_mat <- log2(mat + 1)
    log_df <- as.data.frame(log_mat)
    sample_ids <- colnames(mat)
    pca_full <- dq_pca_full()
    if (!is.null(pca_full)) {
      outlier_samples <- pca_full$scores$Sample[pca_full$scores$Outlier == "Outlier"]
    } else {
      outlier_samples <- character(0)
    }
    IsOutlier <- ifelse(sample_ids %in% outlier_samples, "Yes", "No")
    
    stats_df <- data.frame(
      SampleID = sample_ids,
      MinIntensity = round(apply(log_mat, 2, min, na.rm = TRUE), 3),
      Q25Intensity = round(apply(log_mat, 2, quantile, 0.25, na.rm = TRUE), 3),
      MedianIntensity = round(apply(log_mat, 2, median, na.rm = TRUE), 3),
      Q75Intensity = round(apply(log_mat, 2, quantile, 0.75, na.rm = TRUE), 3),
      MaxIntensity = round(apply(log_mat, 2, max, na.rm = TRUE), 3),
      MeanIntensity = round(colMeans(log_mat, na.rm = TRUE), 3),
      StdDevIntensity = round(apply(log_mat, 2, sd, na.rm = TRUE), 3),
      IsOutlier = IsOutlier
    )
    
    wb <- openxlsx::createWorkbook()
    readme <- data.frame(Description = c(
      "Sheet 'log2_intensity_matrix': log2(Intensity+1) transformed intensity matrix. NA represents missing value.",
      "Sheet 'sample_intensity_stats': summary statistics of log2 intensities per sample. These values directly correspond to the boxplot.",
      "IsOutlier column indicates if the sample was flagged as an outlier based on PCA (z-score > 3)."
    ))
    openxlsx::addWorksheet(wb, "README")
    openxlsx::writeData(wb, "README", readme)
    openxlsx::addWorksheet(wb, "log2_intensity_matrix")
    openxlsx::writeData(wb, "log2_intensity_matrix", log_df)
    openxlsx::addWorksheet(wb, "sample_intensity_stats")
    openxlsx::writeData(wb, "sample_intensity_stats", stats_df)
    openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    message("[DEBUG] download_intensity_data finished")
  }
)

observeEvent(input$help_intensity, {
  showModal(modalDialog(
    title = "Protein Intensity Distribution",
    "箱线图展示每个样本中蛋白强度的分布（log2 转换）。异常样本通常表现出整体偏移或异常离散。",
    easyClose = TRUE, footer = modalButton("关闭")
  ))
})

# ==================== 样本相关性热图 ====================
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
  
  min_cor <- min(cor_mat, na.rm = TRUE)
  max_cor <- max(cor_mat, na.rm = TRUE)
  message("[DEBUG] dq_cor_heatmap_plot_obj: correlation range: [", min_cor, ", ", max_cor, "]")
  if (abs(max_cor - min_cor) < 1e-6) {
    min_cor <- min_cor - 0.01
    max_cor <- max_cor + 0.01
  }
  breaks <- seq(min_cor, max_cor, length.out = 101)
  
  message("[DEBUG] dq_cor_heatmap_plot_obj: drawing pheatmap with dynamic breaks...")
  p <- pheatmap::pheatmap(cor_mat, 
                          main = "Sample Correlation Heatmap (Red = High correlation, Blue = Low correlation)",
                          color = colorRampPalette(c("blue", "white", "red"))(100),
                          breaks = breaks,
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

# 导出样本表达相关矩阵 Excel（含组统计）
output$download_cor_matrix <- downloadHandler(
  filename = function() { paste0("correlation_matrix_", Sys.Date(), ".xlsx") },
  content = function(file) {
    message("[DEBUG] download_cor_matrix Excel started")
    cor_mat <- calculate_sample_correlation(dq_expr_matrix())
    if (is.null(cor_mat)) {
      showNotification("Correlation matrix not available", type = "error")
      return()
    }
    wb <- openxlsx::createWorkbook()
    # Sheet1: correlation matrix
    openxlsx::addWorksheet(wb, "correlation_matrix")
    openxlsx::writeData(wb, "correlation_matrix", cor_mat, rowNames = TRUE)
    
    # Sheet2: group statistics
    group_stats <- NULL
    if (!is.null(rv$sample_info) && "Group" %in% colnames(rv$sample_info)) {
      sample_info_short <- rv$sample_info
      rownames(sample_info_short) <- extract_sample_names(rownames(sample_info_short))
      common_samples <- intersect(colnames(cor_mat), rownames(sample_info_short))
      if (length(common_samples) > 1) {
        groups <- sample_info_short[common_samples, "Group"]
        cor_sub <- cor_mat[common_samples, common_samples]
        group_levels <- unique(groups)
        within_avg <- list()
        between_avg <- list()
        for (g in group_levels) {
          idx <- which(groups == g)
          if (length(idx) > 1) {
            within_avg[[g]] <- mean(cor_sub[idx, idx][lower.tri(cor_sub[idx, idx])], na.rm = TRUE)
          } else {
            within_avg[[g]] <- NA
          }
        }
        if (length(group_levels) >= 2) {
          for (i in 1:(length(group_levels)-1)) {
            for (j in (i+1):length(group_levels)) {
              g1 <- group_levels[i]; g2 <- group_levels[j]
              idx1 <- which(groups == g1); idx2 <- which(groups == g2)
              if (length(idx1) > 0 && length(idx2) > 0) {
                between_avg[[paste(g1, "vs", g2)]] <- mean(cor_sub[idx1, idx2], na.rm = TRUE)
              }
            }
          }
        }
        group_stats <- data.frame(
          Comparison = c(names(within_avg), names(between_avg)),
          Type = c(rep("Within-group", length(within_avg)), rep("Between-group", length(between_avg))),
          AverageCorrelation = unlist(c(within_avg, between_avg))
        )
        openxlsx::addWorksheet(wb, "correlation_group_stats")
        openxlsx::writeData(wb, "correlation_group_stats", group_stats)
      }
    }
    
    readme_text <- c(
      "Sheet 'correlation_matrix': Pearson correlation matrix of log2-transformed protein intensities (top 500 variable proteins), values range -1 to 1.",
      "Sheet 'correlation_group_stats': Average within-group and between-group correlations, if group information is available. This summary can be compared with the overall average correlation displayed in the heatmap caption."
    )
    openxlsx::addWorksheet(wb, "README")
    openxlsx::writeData(wb, "README", data.frame(Description = readme_text))
    
    openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    message("[DEBUG] download_cor_matrix Excel finished")
  }
)

observeEvent(input$help_cor_heatmap, {
  showModal(modalDialog(
    title = "Sample Correlation Heatmap",
    "基于高变异蛋白计算的样本间相关性。颜色越红表示相关性越高。理想情况下，同组样本应聚在一起。",
    easyClose = TRUE, footer = modalButton("关闭")
  ))
})

# ==================== PCA 公共反应式 ====================
dq_pca_full <- reactive({
  req(dq_expr_matrix())
  expr <- dq_expr_matrix()
  log_expr <- log2(expr + 1)
  log_expr[is.na(log_expr)] <- 0
  row_vars <- apply(log_expr, 1, var)
  n_keep <- min(500, nrow(log_expr))
  top_var <- order(row_vars, decreasing = TRUE)[1:n_keep]
  log_expr <- log_expr[top_var, , drop = FALSE]
  row_unique <- apply(log_expr, 1, function(x) length(unique(x)))
  log_expr <- log_expr[row_unique > 1, , drop = FALSE]
  if (nrow(log_expr) < 2) return(NULL)
  pca <- tryCatch(prcomp(t(log_expr), scale. = TRUE), error = function(e) NULL)
  if (is.null(pca)) return(NULL)
  variance <- pca$sdev^2 / sum(pca$sdev^2) * 100
  scores <- as.data.frame(pca$x[, 1:2])
  scores$Sample <- rownames(scores)
  if (!is.null(rv$sample_info)) {
    si_short <- rv$sample_info
    rownames(si_short) <- extract_sample_names(rownames(si_short))
    common <- intersect(scores$Sample, rownames(si_short))
    if (length(common) > 0) {
      scores$Group <- si_short[common, "Group"]
      scores$Batch <- if ("Batch" %in% colnames(si_short)) si_short[common, "Batch"] else NA
    } else {
      scores$Group <- "All"
      scores$Batch <- NA
    }
  } else {
    scores$Group <- "All"
    scores$Batch <- NA
  }
  outliers <- get_outlier_samples(list(pca_df = scores, pc1_var = variance[1], pc2_var = variance[2]))
  scores$Outlier <- ifelse(scores$Sample %in% outliers, "Outlier", "Normal")
  list(pca = pca, scores = scores, variance = variance, loadings = pca$rotation)
})

# ==================== PCA Group 图 ====================
pca_group_plot_obj <- reactive({
  pca_full <- dq_pca_full()
  if (is.null(pca_full)) return(NULL)
  scores <- pca_full$scores
  scores <- scores[!is.na(scores$Group) & scores$Group != "a", ]
  group_colors <- c("Control" = "#FF69B4", "Treatment" = "#00CED1")
  all_groups <- unique(scores$Group)
  missing_colors <- setdiff(all_groups, names(group_colors))
  if (length(missing_colors) > 0) {
    extra_colors <- rainbow(length(missing_colors))
    names(extra_colors) <- missing_colors
    group_colors <- c(group_colors, extra_colors)
  }
  ggplot(scores, aes(x = PC1, y = PC2, color = Group, shape = Outlier)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_text(data = scores[scores$Outlier == "Outlier", ], aes(label = Sample), vjust = 1.5, size = 3, show.legend = FALSE) +
    scale_shape_manual(values = c("Normal" = 16, "Outlier" = 17), labels = c("Normal" = "Normal", "Outlier" = "Outlier")) +
    scale_color_manual(values = group_colors) +
    labs(title = "PCA by Group (Outliers in red ▲)",
         x = paste0("PC1 (", round(pca_full$variance[1], 1), "%)"),
         y = paste0("PC2 (", round(pca_full$variance[2], 1), "%)")) +
    theme_bw() + theme(legend.position = "right",
                       axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
})

output$dq_pca_group_plot <- renderPlot({ req(pca_group_plot_obj()); pca_group_plot_obj() })
output$download_pca_group <- downloadHandler(
  filename = function() "pca_group.png",
  content = function(file) ggsave(file, plot = pca_group_plot_obj(), width = 8, height = 6, dpi = 150)
)
output$download_pca_group_data <- downloadHandler(
  filename = function() { paste0("PCA_Group_Data_", Sys.Date(), ".xlsx") },
  content = function(file) {
    message("[DEBUG] download_pca_group_data started")
    pca_full <- dq_pca_full()
    if (is.null(pca_full)) {
      showNotification("PCA not available", type = "error")
      return()
    }
    wb <- openxlsx::createWorkbook()
    scores_df <- pca_full$scores
    openxlsx::addWorksheet(wb, "pca_scores")
    openxlsx::writeData(wb, "pca_scores", scores_df)
    loadings_df <- as.data.frame(pca_full$loadings[, 1:2])
    loadings_df$ProteinID <- rownames(loadings_df)
    loadings_df <- loadings_df[, c("ProteinID", "PC1", "PC2")]
    openxlsx::addWorksheet(wb, "pca_loadings")
    openxlsx::writeData(wb, "pca_loadings", loadings_df)
    var_df <- data.frame(PC = seq_along(pca_full$variance), VarianceExplained = pca_full$variance)
    openxlsx::addWorksheet(wb, "pca_variance_explained")
    openxlsx::writeData(wb, "pca_variance_explained", var_df)
    readme <- data.frame(Description = c(
      "Sheet 'pca_scores': coordinates of samples in PC1/PC2 space. Columns: Sample, PC1, PC2, Group, Batch, Outlier (whether identified as potential outlier).",
      "Sheet 'pca_loadings': protein loadings for PC1 and PC2 (top 500 variable proteins).",
      "Sheet 'pca_variance_explained': percentage of variance explained by each principal component."
    ))
    openxlsx::addWorksheet(wb, "README")
    openxlsx::writeData(wb, "README", readme)
    openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    message("[DEBUG] download_pca_group_data finished")
  }
)
observeEvent(input$help_pca_group, {
  showModal(modalDialog(
    title = "PCA by Group",
    "按实验分组着色，用于观察组间分离程度。异常样本已用红色三角形标注。",
    easyClose = TRUE, footer = modalButton("关闭")
  ))
})

# ==================== PCA Batch 图 ====================
pca_batch_plot_obj <- reactive({
  pca_full <- dq_pca_full()
  if (is.null(pca_full)) return(NULL)
  scores <- pca_full$scores
  if (all(is.na(scores$Batch))) return(NULL)
  scores <- scores[!is.na(scores$Batch) & scores$Batch != "a", ]
  batch_colors <- c("Batch1" = "#E41A1C", "Batch2" = "#00CED1")
  all_batches <- unique(scores$Batch)
  if (length(all_batches) > 2) {
    batch_colors <- setNames(rainbow(length(all_batches)), all_batches)
  }
  ggplot(scores, aes(x = PC1, y = PC2, color = Batch, shape = Outlier)) +
    geom_point(size = 3, alpha = 0.8) +
    geom_text(data = scores[scores$Outlier == "Outlier", ], aes(label = Sample), vjust = 1.5, size = 3, show.legend = FALSE) +
    scale_shape_manual(values = c("Normal" = 16, "Outlier" = 17)) +
    scale_color_manual(values = batch_colors) +
    labs(title = "PCA by Batch (Outliers in red ▲)",
         x = paste0("PC1 (", round(pca_full$variance[1], 1), "%)"),
         y = paste0("PC2 (", round(pca_full$variance[2], 1), "%)")) +
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
output$download_pca_batch_data <- downloadHandler(
  filename = function() { paste0("PCA_Batch_Data_", Sys.Date(), ".xlsx") },
  content = function(file) {
    message("[DEBUG] download_pca_batch_data started")
    pca_full <- dq_pca_full()
    if (is.null(pca_full)) {
      showNotification("PCA not available", type = "error")
      return()
    }
    wb <- openxlsx::createWorkbook()
    scores_df <- pca_full$scores
    openxlsx::addWorksheet(wb, "pca_scores")
    openxlsx::writeData(wb, "pca_scores", scores_df)
    loadings_df <- as.data.frame(pca_full$loadings[, 1:2])
    loadings_df$ProteinID <- rownames(loadings_df)
    loadings_df <- loadings_df[, c("ProteinID", "PC1", "PC2")]
    openxlsx::addWorksheet(wb, "pca_loadings")
    openxlsx::writeData(wb, "pca_loadings", loadings_df)
    var_df <- data.frame(PC = seq_along(pca_full$variance), VarianceExplained = pca_full$variance)
    openxlsx::addWorksheet(wb, "pca_variance_explained")
    openxlsx::writeData(wb, "pca_variance_explained", var_df)
    readme <- data.frame(Description = c(
      "Sheet 'pca_scores': coordinates of samples in PC1/PC2 space. Columns: Sample, PC1, PC2, Group, Batch, Outlier.",
      "Sheet 'pca_loadings': protein loadings for PC1 and PC2.",
      "Sheet 'pca_variance_explained': percentage of variance explained by each PC."
    ))
    openxlsx::addWorksheet(wb, "README")
    openxlsx::writeData(wb, "README", readme)
    openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    message("[DEBUG] download_pca_batch_data finished")
  }
)
observeEvent(input$help_pca_batch, {
  showModal(modalDialog(
    title = "PCA by Batch",
    "按实验批次着色，用于检测批次效应。异常样本已用红色三角形标注。",
    easyClose = TRUE, footer = modalButton("关闭")
  ))
})

# ==================== 预处理跳转按钮 ====================
observeEvent(input$goto_preprocessing, {
  message("Jumping to Data Preprocessing page...")
  updateNavbarPage(session, "main_navbar", selected = "preprocessing")
})