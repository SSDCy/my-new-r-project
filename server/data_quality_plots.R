# server/data_quality_plots.R

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
output$dq_missing_score <- renderText({
  req(dq_expr_matrix())
  tryCatch(paste0("Score: ", dq_quality_score()$details$missing_score, "/30"), error = function(e) paste("Error:", e$message))
})
output$dq_consistency <- renderText({
  req(dq_expr_matrix())
  tryCatch({
    if (dq_quality_score()$details$avg_correlation > 0.8) "Good"
    else if (dq_quality_score()$details$avg_correlation > 0.7) "Fair" else "Poor"
  }, error = function(e) paste("Error:", e$message))
})
output$dq_consistency_score <- renderText({
  req(dq_expr_matrix())
  tryCatch(paste0("Score: ", dq_quality_score()$details$correlation_score, "/20"), error = function(e) paste("Error:", e$message))
})
output$dq_protein_quality <- renderText({
  req(dq_expr_matrix())
  tryCatch({
    if (dq_quality_score()$details$protein_valid_ratio > 80) "Good"
    else if (dq_quality_score()$details$protein_valid_ratio > 60) "Fair" else "Poor"
  }, error = function(e) paste("Error:", e$message))
})
output$dq_protein_score <- renderText({
  req(dq_expr_matrix())
  tryCatch(paste0("Score: ", dq_quality_score()$details$protein_score, "/20"), error = function(e) paste("Error:", e$message))
})

# ==================== 智能报告 ====================
output$dq_total_score <- renderText({
  req(dq_expr_matrix())
  tryCatch(dq_quality_score()$score, error = function(e) paste("Error:", e$message))
})
output$dq_total_grade <- renderText({
  req(dq_expr_matrix())
  tryCatch(dq_quality_score()$grade, error = function(e) paste("Error:", e$message))
})
output$dq_total_missing <- renderText({
  req(dq_expr_matrix())
  tryCatch(paste0(dq_quality_score()$details$missing_ratio, "%"), error = function(e) paste("Error:", e$message))
})
output$dq_total_consistency <- renderText({
  req(dq_expr_matrix())
  tryCatch({
    if (dq_quality_score()$details$avg_correlation > 0.8) "Good"
    else if (dq_quality_score()$details$avg_correlation > 0.7) "Fair" else "Poor"
  }, error = function(e) paste("Error:", e$message))
})
output$dq_total_correlation <- renderText({
  req(dq_expr_matrix())
  tryCatch(dq_quality_score()$details$avg_correlation, error = function(e) paste("Error:", e$message))
})

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
output$dq_missing_heatmap <- renderPlot({
  req(dq_expr_matrix())
  tryCatch({
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
                        labels = c("0" = "Present", "1" = "Missing"), name = "Status") +
      labs(title = "Missing Value Heatmap", x = "", y = "") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none")
  }, error = function(e) { plot.new(); text(0.5, 0.5, paste("Error:", e$message)) })
})

# ==================== 有效值柱状图 ====================
output$dq_valid_values_plot <- renderPlot({
  req(dq_expr_matrix())
  tryCatch({
    req(dq_missing_stats())
    stats <- dq_missing_stats()
    valid_counts <- round((1 - stats$sample_missing) * nrow(dq_expr_matrix()))
    valid_percent <- (1 - stats$sample_missing) * 100
    df <- data.frame(Sample = factor(names(valid_counts), levels = names(valid_counts)), ValidPercent = valid_percent)
    threshold <- 70
    ggplot(df, aes(x = Sample, y = ValidPercent)) +
      geom_col(aes(fill = ValidPercent < threshold)) +
      geom_hline(yintercept = threshold, color = "red", linetype = "dashed", linewidth = 1) +
      scale_fill_manual(values = c("TRUE" = "#e74c3c", "FALSE" = "#2ecc71"), guide = "none") +
      labs(title = "Valid Values per Sample", x = "", y = "Valid Values (%)") +
      ylim(0, 80) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
  }, error = function(e) { plot.new(); text(0.5, 0.5, paste("Error:", e$message)) })
})

# ==================== 缺失值相关性热图 ====================
output$debug_missing_cor_status <- renderText({
  if (!isTRUE(dq_expr_matrix())) return("No expression data uploaded.")
  return("Container is active, attempting plot...")
})

output$dq_missing_cor_plot <- renderPlot({
  req(dq_expr_matrix())
  tryCatch({
    mat <- dq_expr_matrix()
    missing_mat <- is.na(mat) * 1
    if (ncol(missing_mat) < 2) {
      plot.new(); text(0.5, 0.5, "Not enough samples")
      return()
    }
    cor_mat <- cor(missing_mat, use = "pairwise.complete.obs")
    cor_mat[is.na(cor_mat)] <- 0
    pheatmap(cor_mat, main = "Missing Value Correlation",
             color = colorRampPalette(c("blue", "white", "red"))(100),
             show_rownames = TRUE, show_colnames = TRUE,
             fontsize_row = 8, fontsize_col = 8)
  }, error = function(e) {
    plot.new(); text(0.5, 0.5, paste("Error:", e$message))
  })
})

# ==================== 强度分布箱线图 ====================
output$dq_intensity_dist_plot <- renderPlot({
  req(dq_expr_matrix())
  tryCatch({
    mat <- dq_expr_matrix()
    log_mat <- log2(mat + 1)
    df <- reshape2::melt(as.matrix(log_mat))
    colnames(df) <- c("Protein", "Sample", "Log2Intensity")
    ggplot(df, aes(x = Sample, y = Log2Intensity)) +
      geom_boxplot(fill = "#3498db", alpha = 0.7, outlier.size = 1) +
      labs(title = "Intensity Distribution (log2)", x = "", y = "log2(Intensity)") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
  }, error = function(e) { plot.new(); text(0.5, 0.5, paste("Error:", e$message)) })
})

# ==================== 样本相关性热图 ====================
output$debug_cor_heatmap_status <- renderText({
  if (!isTRUE(dq_expr_matrix())) return("No expression data uploaded.")
  return("Container is active, attempting plot...")
})

output$dq_cor_heatmap <- renderPlot({
  req(dq_expr_matrix())
  tryCatch({
    cor_mat <- calculate_sample_correlation(dq_expr_matrix())
    if (is.null(cor_mat)) {
      plot.new(); text(0.5, 0.5, "Not enough data")
      return()
    }
    ann_col <- NULL; ann_colors <- NULL
    if (!is.null(rv$sample_info) && "Group" %in% colnames(rv$sample_info)) {
      sample_info_short <- rv$sample_info
      rownames(sample_info_short) <- extract_sample_names(rownames(sample_info_short))
      common_samples <- intersect(colnames(cor_mat), rownames(sample_info_short))
      if (length(common_samples) > 0) {
        ann_col <- data.frame(Group = sample_info_short[common_samples, "Group"], row.names = common_samples)
        ann_col <- ann_col[!is.na(ann_col$Group) & ann_col$Group != "a", , drop = FALSE]
        if (nrow(ann_col) > 0) {
          cor_mat <- cor_mat[rownames(ann_col), rownames(ann_col)]
          groups <- unique(ann_col$Group)
          ann_colors <- list(Group = get_group_colors(groups))
        }
      }
    }
    pheatmap(cor_mat, main = "Sample Correlation Heatmap",
             color = colorRampPalette(c("blue", "white", "red"))(100),
             breaks = seq(0.4, 1, length.out = 101),
             show_rownames = TRUE, show_colnames = TRUE,
             fontsize_row = 9, fontsize_col = 9,
             annotation_col = ann_col, annotation_colors = ann_colors)
  }, error = function(e) {
    plot.new(); text(0.5, 0.5, paste("Error:", e$message))
  })
})

# ==================== PCA 图 ====================
output$dq_pca_plot <- renderPlot({
  req(dq_expr_matrix())
  tryCatch({
    pca_result <- calculate_pca(dq_expr_matrix(), rv$sample_info)
    if (is.null(pca_result)) {
      plot.new(); text(0.5, 0.5, "Not enough data for PCA")
      return()
    }
    pca_df <- pca_result$pca_df
    pca_df <- pca_df[!is.na(pca_df$Group) & pca_df$Group != "a", ]
    ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
      geom_point(size = 3, alpha = 0.8) +
      geom_text(vjust = 1.5, size = 3) +
      labs(title = "PCA Plot (Unsupervised Clustering)",
           x = paste0("PC1 (", pca_result$pc1_var, "%)"),
           y = paste0("PC2 (", pca_result$pc2_var, "%)")) +
      theme_bw() + theme(legend.position = "right")
  }, error = function(e) { plot.new(); text(0.5, 0.5, paste("Error:", e$message)) })
})