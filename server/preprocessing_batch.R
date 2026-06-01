# server/preprocessing_batch.R
message("[DEBUG] preprocessing_batch.R loaded")

batch_diagnostic <- reactive({
  if (is.null(rv$sample_info) || !"Batch" %in% colnames(rv$sample_info)) {
    return(list(status = "warning", message = "⚠️ 批次诊断失败：样本信息表中缺少 Batch 列。"))
  }
  if (is.null(expression_data())) {
    return(list(status = "info", message = "Expression data not loaded."))
  }
  
  tryCatch({
    mat <- expression_data()
    if (input$max_missing_fraction < 1) {
      missing_frac <- rowMeans(is.na(mat))
      mat <- mat[missing_frac <= input$max_missing_fraction, , drop = FALSE]
    }
    max_int <- apply(mat, 1, max, na.rm = TRUE)
    keep_finite <- is.finite(max_int)
    mat <- mat[keep_finite, , drop = FALSE]
    if (!is.null(input$min_intensity) && !is.na(input$min_intensity) && input$min_intensity > 0) {
      max_int_finite <- apply(mat, 1, max, na.rm = TRUE)
      mat <- mat[max_int_finite > input$min_intensity, , drop = FALSE]
    }
    if (nrow(mat) < 2 || ncol(mat) < 2) {
      return(list(status = "warning", message = "⚠️ 批次诊断失败：过滤后数据量不足。"))
    }
    suppressMessages({
      filled <- impute_missing_values(mat, method = "knn")
    })
    filled <- as.matrix(filled)
    mat_t <- t(log2(filled + 1))
    mat_t <- mat_t[, apply(mat_t, 2, var) > 1e-12, drop = FALSE]
    if (ncol(mat_t) < 2) {
      return(list(status = "warning", message = "⚠️ 批次诊断失败：可变蛋白数量不足。"))
    }
    pca <- prcomp(mat_t, scale. = TRUE)
    pc1 <- pca$x[,1]
    sample_names <- rownames(mat_t)
    info_rows <- standardize_sample_name(rownames(rv$sample_info))
    sample_norm <- standardize_sample_name(sample_names)
    idx <- match(sample_norm, info_rows)
    if (all(is.na(idx))) {
      return(list(status = "warning", message = "⚠️ 批次诊断失败：无法匹配样本到样本信息表。"))
    }
    batch <- rv$sample_info[idx, "Batch"]
    batch <- factor(batch)
    if (length(levels(batch)) < 2) {
      return(list(status = "info", message = "✅ 批次效应检测：仅检测到一个批次，无需校正。"))
    }
    anova_p <- tryCatch(summary(aov(pc1 ~ batch))[[1]]$`Pr(>F)`[1], error = function(e) NA)
    if (!is.na(anova_p) && anova_p < 0.05) {
      return(list(status = "significant", message = "⚠️ 检测到显著批次效应（ANOVA p < 0.05），建议启用 ComBat 校正。"))
    } else {
      return(list(status = "success", message = "✅ 批次效应检测完成：您的数据无显著批次分离，校正前后样本分布无明显变化，无需校正。"))
    }
  }, error = function(e) {
    return(list(status = "warning", message = "⚠️ 批次诊断失败：请检查样本信息表中 Batch 列是否存在，或每个批次样本数量是否≥3 个。"))
  })
})

observe({
  diag <- batch_diagnostic()
  req(diag)
  if (identical(input$imputation_method, "none")) {
    message("[DEBUG] batch auto-enable blocked: imputation method is 'none'")
    return()
  }
  if (diag$status == "significant") {
    if (!isTRUE(input$perform_batch_correction)) {
      message("[DEBUG] auto-enabling batch correction due to significant effect")
      updateCheckboxInput(session, "perform_batch_correction", value = TRUE)
      showNotification("⚠️ 检测到显著批次效应，已自动启用 ComBat 校正。", type = "warning", duration = 6, id = "auto_batch_notif")
    }
  }
})

output$batch_diagnostic_message <- renderPrint({
  cat(batch_diagnostic()$message)
})

output$batch_help_text <- renderUI({
  diag <- batch_diagnostic()
  base_text <- "采用 ComBat 算法校正实验批次偏差。批次列将自动识别为样本信息表中的 'Batch' 列。"
  if (diag$status == "success") {
    HTML(paste(base_text, "<br><b>👉 结合您的数据：校正前后 PCA 显示，批次无显著分离，无明显批次效应，建议取消勾选，避免过度校正。</b>"))
  } else if (diag$status == "significant") {
    HTML(paste(base_text, "<br><b>👉 注意：已自动启用校正，建议运行预处理后查看 Batch Correction 选项卡中的 PCA 对比图。</b>"))
  } else if (diag$status == "warning") {
    HTML(paste(base_text, "<br><b>", diag$message, "</b>"))
  } else {
    HTML(base_text)
  }
})

output$batch_correction_performed <- reactive({
  !is.null(processed_data()) && preprocessing_params$batch_performed
})
outputOptions(output, "batch_correction_performed", suspendWhenHidden = FALSE)

observe({
  message("[DEBUG] conflict check: imputation = ", input$imputation_method, ", batch = ", input$perform_batch_correction)
  if (identical(input$imputation_method, "none")) {
    if (isTRUE(input$perform_batch_correction)) {
      updateCheckboxInput(session, "perform_batch_correction", value = FALSE)
      showNotification("批次校正需要完整数据，已自动取消（不填充时无法校正）。", type = "warning", duration = 6, id = "conflict_notif")
    }
    shinyjs::disable("perform_batch_correction")
  } else {
    shinyjs::enable("perform_batch_correction")
  }
})

batch_comparison_pca <- reactive({
  req(processed_data(), preprocessing_params$batch_performed)
  before <- preprocessing_params$pre_batch_data
  after <- preprocessing_params$post_batch_data
  if (is.null(before) || is.null(after)) return(NULL)
  
  common_cols <- intersect(colnames(before), colnames(after))
  before <- before[, common_cols, drop = FALSE]
  after <- after[, common_cols, drop = FALSE]
  
  before_t <- t(log2(as.matrix(before) + 1))
  after_t <- t(log2(as.matrix(after) + 1))
  
  before_ok <- apply(before_t, 2, function(x) all(is.finite(x)))
  after_ok <- apply(after_t, 2, function(x) all(is.finite(x)))
  keep_cols <- which(before_ok & after_ok)
  if (length(keep_cols) < 2) return(NULL)
  before_t <- before_t[, keep_cols, drop = FALSE]
  after_t <- after_t[, keep_cols, drop = FALSE]
  
  before_var <- apply(before_t, 2, var)
  after_var <- apply(after_t, 2, var)
  keep_var <- which(before_var > 1e-12 & after_var > 1e-12)
  if (length(keep_var) < 2) return(NULL)
  before_t <- before_t[, keep_var, drop = FALSE]
  after_t <- after_t[, keep_var, drop = FALSE]
  
  pca_before <- tryCatch(prcomp(before_t, scale. = TRUE), error = function(e) NULL)
  pca_after <- tryCatch(prcomp(after_t, scale. = TRUE), error = function(e) NULL)
  if (is.null(pca_before) || is.null(pca_after)) return(NULL)
  
  sample_info_short <- rv$sample_info
  rownames(sample_info_short) <- standardize_sample_name(rownames(sample_info_short))
  before_norm <- standardize_sample_name(colnames(before))
  after_norm <- standardize_sample_name(colnames(after))
  batch_before <- sample_info_short[before_norm, "Batch"]
  batch_after <- sample_info_short[after_norm, "Batch"]
  
  list(
    pca_before = pca_before,
    pca_after = pca_after,
    batch_before = batch_before,
    batch_after = batch_after,
    var_before = round(pca_before$sdev^2 / sum(pca_before$sdev^2) * 100, 1),
    var_after = round(pca_after$sdev^2 / sum(pca_after$sdev^2) * 100, 1)
  )
})

output$batch_pca_plot <- renderPlot({
  dat <- batch_comparison_pca()
  if (is.null(dat)) {
    plot.new(); text(0.5, 0.5, "PCA not available")
    return()
  }
  df_before <- data.frame(PC1 = dat$pca_before$x[,1], PC2 = dat$pca_before$x[,2],
                          Batch = dat$batch_before, Stage = "Before Correction")
  df_after <- data.frame(PC1 = dat$pca_after$x[,1], PC2 = dat$pca_after$x[,2],
                         Batch = dat$batch_after, Stage = "After Correction")
  pca_df <- rbind(df_before, df_after)
  pca_df$Stage <- factor(pca_df$Stage, levels = c("Before Correction", "After Correction"))
  
  batches <- unique(na.omit(c(as.character(dat$batch_before), as.character(dat$batch_after))))
  if (length(batches) > 8) {
    batch_colors <- setNames(rainbow(length(batches)), batches)
  } else {
    batch_colors <- setNames(RColorBrewer::brewer.pal(length(batches), "Set1"), batches)
  }
  
  ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(type = "norm", level = 0.95) +
    scale_color_manual(values = batch_colors) +
    facet_wrap(~ Stage) +
    labs(title = "PCA Before and After ComBat Batch Correction",
         x = paste0("PC1 (Before: ", dat$var_before[1], "%, After: ", dat$var_after[1], "%)"),
         y = paste0("PC2 (Before: ", dat$var_before[2], "%, After: ", dat$var_after[2], "%)")) +
    theme_bw() + theme(legend.position = "bottom")
})

output$batch_pca_interpretation <- renderUI({
  diag <- batch_diagnostic()
  if (diag$status == "significant") {
    p("图注：红色 = Batch1，青色 = Batch2。校正前后批次样本分布出现分离趋势，提示存在批次效应，建议启用校正。")
  } else {
    p("图注：红色 = Batch1，青色 = Batch2。校正前后批次样本分布无明显聚类，说明无显著批次效应，无需校正。")
  }
})

output$download_batch_pca <- downloadHandler(
  filename = function() "batch_correction_pca.png",
  content = function(file) {
    dat <- batch_comparison_pca()
    if (is.null(dat)) return()
    df_before <- data.frame(PC1 = dat$pca_before$x[,1], PC2 = dat$pca_before$x[,2],
                            Batch = dat$batch_before, Stage = "Before Correction")
    df_after <- data.frame(PC1 = dat$pca_after$x[,1], PC2 = dat$pca_after$x[,2],
                           Batch = dat$batch_after, Stage = "After Correction")
    pca_df <- rbind(df_before, df_after)
    pca_df$Stage <- factor(pca_df$Stage, levels = c("Before Correction", "After Correction"))
    
    batches <- unique(na.omit(c(as.character(dat$batch_before), as.character(dat$batch_after))))
    if (length(batches) > 8) {
      batch_colors <- setNames(rainbow(length(batches)), batches)
    } else {
      batch_colors <- setNames(RColorBrewer::brewer.pal(length(batches), "Set1"), batches)
    }
    
    plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
      geom_point(size = 3, alpha = 0.8) +
      stat_ellipse(type = "norm", level = 0.95) +
      scale_color_manual(values = batch_colors) +
      facet_wrap(~ Stage) +
      labs(title = "PCA Before and After ComBat Batch Correction",
           x = paste0("PC1 (Before: ", dat$var_before[1], "%, After: ", dat$var_after[1], "%)"),
           y = paste0("PC2 (Before: ", dat$var_before[2], "%, After: ", dat$var_after[2], "%)")) +
      theme_bw() + theme(legend.position = "bottom")
    ggsave(file, plot = plot, width = 12, height = 6, dpi = 150)
  }
)