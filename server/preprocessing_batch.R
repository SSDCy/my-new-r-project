# server/preprocessing_batch.R
message("[DEBUG] preprocessing_batch.R loaded - step-by-step batch visualization")

# =========================================================================
# 注意：preprocessing_params 由 preprocessing_core.R 定义，
# 本文件不再重复定义，以免覆盖已有的值。
# =========================================================================

# ---------- 调试：在响应式上下文中确认 preprocessing_params 可访问 ----------
shiny::observe({
  # 使用 tryCatch 确保即使访问失败也不会导致闪退
  tryCatch({
    if (exists("preprocessing_params", envir = parent.env(environment()), inherits = TRUE)) {
      msg <- capture.output({
        cat("[DEBUG] preprocessing_batch.R: preprocessing_params is accessible")
        cat(", batch_performed =", preprocessing_params$batch_performed)
        cat(", batch_corrected_cols =", if (is.null(preprocessing_params$batch_corrected_cols)) "NULL" else paste(preprocessing_params$batch_corrected_cols, collapse=","))
      })
      message(msg)
    } else {
      message("[WARN] preprocessing_batch.R: preprocessing_params not found in parent environment.")
    }
  }, error = function(e) {
    message("[ERROR] preprocessing_batch.R: Could not access preprocessing_params: ", e$message)
  })
}, priority = 10)  # 较低优先级以便环境已建立

# ---------- 批次诊断 ----------
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

output$batch_diagnostic_ready <- reactive({
  diag <- batch_diagnostic()
  !is.null(rv$sample_info) && "Batch" %in% colnames(rv$sample_info) && length(unique(rv$sample_info$Batch)) >= 2
})
outputOptions(output, "batch_diagnostic_ready", suspendWhenHidden = FALSE)

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

# ---------- 更新批次选择下拉框 ----------
observe({
  if (!is.null(rv$sample_info) && "Batch" %in% colnames(rv$sample_info)) {
    batches <- unique(rv$sample_info$Batch)
    updateSelectInput(session, "batch_verification_batch1", choices = batches, selected = batches[1])
    updateSelectInput(session, "batch_verification_batch2", choices = batches, selected = if(length(batches)>1) batches[2] else batches[1])
    message("[DEBUG] batch verification choices updated: ", paste(batches, collapse = ", "))
  }
})

# ---------- 批次验证：准备数据（返回中间矩阵供可视化使用） ----------
batch_verification_data <- reactive({
  req(input$batch_verification_batch1, input$batch_verification_batch2)
  if (input$batch_verification_batch1 == input$batch_verification_batch2) {
    message("[DEBUG] batch_verification_data: same batch selected")
    return(NULL)
  }
  if (is.null(expression_data())) {
    message("[DEBUG] batch_verification_data: no expression data")
    return(NULL)
  }
  
  tryCatch({
    mat <- expression_data()
    original_proteins <- nrow(mat)
    if (input$max_missing_fraction < 1) {
      missing_frac <- rowMeans(is.na(mat))
      mat <- mat[missing_frac <= input$max_missing_fraction, , drop = FALSE]
    }
    after_missing <- nrow(mat)
    max_int <- apply(mat, 1, max, na.rm = TRUE)
    keep_finite <- is.finite(max_int)
    mat <- mat[keep_finite, , drop = FALSE]
    after_inf <- nrow(mat)
    if (!is.null(input$min_intensity) && !is.na(input$min_intensity) && input$min_intensity > 0) {
      max_int_finite <- apply(mat, 1, max, na.rm = TRUE)
      mat <- mat[max_int_finite > input$min_intensity, , drop = FALSE]
    }
    after_intensity <- nrow(mat)
    if (nrow(mat) < 2) {
      message("[DEBUG] batch_verification_data: not enough proteins after filtering")
      return(NULL)
    }
    
    # 保存过滤后的原始强度矩阵
    raw_filtered <- mat
    
    suppressMessages({
      filled <- impute_missing_values(mat, method = "knn")
    })
    filled <- as.matrix(filled)
    log_mat <- log2(filled + 1)
    
    # 转置，使样本为行，蛋白为列
    mat_t <- t(log_mat)
    before_var_filter <- ncol(mat_t)
    mat_t <- mat_t[, apply(mat_t, 2, var) > 1e-12, drop = FALSE]
    after_var_filter <- ncol(mat_t)
    if (ncol(mat_t) < 2) {
      message("[DEBUG] batch_verification_data: not enough variable proteins")
      return(NULL)
    }
    
    pca <- prcomp(mat_t, scale. = TRUE)
    pc1 <- pca$x[,1]
    pc2 <- pca$x[,2]
    sample_names_full <- rownames(mat_t)
    
    # 匹配批次信息
    si <- rv$sample_info
    info_norm <- standardize_sample_name(rownames(si))
    sample_norm <- standardize_sample_name(sample_names_full)
    idx <- match(sample_norm, info_norm)
    if (all(is.na(idx))) {
      message("[DEBUG] batch_verification_data: cannot match any sample to batch info")
      return(NULL)
    }
    
    batch <- si$Batch[idx]
    sample_display <- extract_sample_names(sample_names_full)
    
    b1 <- input$batch_verification_batch1
    b2 <- input$batch_verification_batch2
    
    pc1_b1 <- pc1[batch == b1]
    pc1_b2 <- pc1[batch == b2]
    names(pc1_b1) <- sample_display[batch == b1]
    names(pc1_b2) <- sample_display[batch == b2]
    
    # 为可视化准备 data.frame
    pca_scores <- data.frame(PC1 = pc1, PC2 = pc2, Batch = batch, Sample = sample_display,
                             stringsAsFactors = FALSE)
    
    message("[DEBUG] batch_verification_data: matched ", length(pc1_b1), " samples to ", b1,
            ", ", length(pc1_b2), " samples to ", b2)
    
    list(
      batch1 = b1, batch2 = b2,
      pc1_b1 = pc1_b1, pc1_b2 = pc1_b2,
      sample_names_b1 = names(pc1_b1),
      sample_names_b2 = names(pc1_b2),
      raw_filtered = raw_filtered,        # 原始过滤后强度
      log_mat = log_mat,                  # log2 转换后矩阵（蛋白×样本）
      pca_obj = pca,
      pca_scores = pca_scores,
      n_original_proteins = original_proteins,
      n_after_missing = after_missing,
      n_after_inf = after_inf,
      n_after_intensity = after_intensity,
      n_proteins_used = nrow(mat),
      n_variable_proteins = after_var_filter,
      n_removed_low_var = before_var_filter - after_var_filter
    )
  }, error = function(e) {
    message("[DEBUG] batch_verification_data error: ", e$message)
    NULL
  })
})

# ---------- 样本 PC1 值表格 ----------
output$batch_verification_table <- renderTable({
  data <- batch_verification_data()
  req(data)
  df1 <- data.frame(Sample = data$sample_names_b1, Batch = data$batch1, PC1 = round(data$pc1_b1, 4), stringsAsFactors = FALSE)
  df2 <- data.frame(Sample = data$sample_names_b2, Batch = data$batch2, PC1 = round(data$pc1_b2, 4), stringsAsFactors = FALSE)
  rbind(df1, df2)
}, striped = TRUE, bordered = TRUE, width = "100%")

# ---------- 小提琴图/箱线图 ----------
output$batch_verification_plot <- renderPlot({
  data <- batch_verification_data()
  req(data)
  
  df <- data.frame(
    PC1 = c(data$pc1_b1, data$pc1_b2),
    Batch = c(rep(data$batch1, length(data$pc1_b1)), rep(data$batch2, length(data$pc1_b2))),
    stringsAsFactors = FALSE
  )
  
  ggplot(df, aes(x = Batch, y = PC1, fill = Batch)) +
    geom_violin(alpha = 0.5, draw_quantiles = 0.5) +
    geom_jitter(width = 0.1, height = 0, size = 2, alpha = 0.8) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "red") +
    labs(title = "PC1 Distribution by Batch",
         subtitle = "Red diamond: mean. Horizontal line: median.",
         y = "PC1 score") +
    theme_bw() + theme(legend.position = "none")
})

# ---------- 详细计算过程（包含 PC1 计算说明） ----------
output$batch_verification_details <- renderPrint({
  data <- batch_verification_data()
  req(data)
  
  b1 <- data$batch1; b2 <- data$batch2
  x <- data$pc1_b1; y <- data$pc1_b2
  n1 <- length(x); n2 <- length(y)
  mean1 <- mean(x); mean2 <- mean(y)
  var1 <- var(x); var2 <- var(y)
  sp <- sqrt(((n1-1)*var1 + (n2-1)*var2) / (n1+n2-2))
  se <- sp * sqrt(1/n1 + 1/n2)
  t_stat <- (mean1 - mean2) / se
  df <- n1 + n2 - 2
  p_val <- 2 * pt(abs(t_stat), df, lower.tail = FALSE)
  
  cat("===========================================================================\n")
  cat("Step 1: Data preparation\n")
  cat("===========================================================================\n")
  cat("Original proteins: ", data$n_original_proteins, "\n", sep = "")
  cat("After missing value filter (threshold ", input$max_missing_fraction, "): ", data$n_after_missing, "\n", sep = "")
  cat("After Inf/non‑finite filter: ", data$n_after_inf, "\n", sep = "")
  if (!is.null(input$min_intensity) && input$min_intensity > 0)
    cat("After minimum intensity filter (threshold ", input$min_intensity, "): ", data$n_after_intensity, "\n", sep = "")
  cat("Proteins used: ", data$n_proteins_used, "\n", sep = "")
  cat("Samples: ", n1 + n2, "\n\n", sep = "")
  
  cat("===========================================================================\n")
  cat("Step 2: Log2 transformation\n")
  cat("===========================================================================\n")
  cat("Each expression value x is transformed: y = log2(x + 1)\n")
  cat("This reduces skewness and brings data closer to a normal distribution.\n\n")
  
  cat("===========================================================================\n")
  cat("Step 3: Remove low‑variance proteins\n")
  cat("===========================================================================\n")
  cat("Proteins with variance ≤ 1e-12 across samples are removed.\n")
  cat("Removed: ", data$n_removed_low_var, " low‑variance proteins.\n")
  cat("Variable proteins retained: ", data$n_variable_proteins, "\n\n", sep = "")
  
  cat("===========================================================================\n")
  cat("Step 4: Principal Component Analysis (PCA)\n")
  cat("===========================================================================\n")
  cat("PCA is performed on the samples × variable‑proteins matrix,\n")
  cat("with columns centered to mean = 0 and scaled to unit variance.\n")
  cat("The first principal component (PC1) captures the largest variance\n")
  cat("direction among the samples.\n")
  cat("Variance explained by PC1: ", round(summary(data$pca_obj)$importance[2,1]*100, 1), "%\n\n", sep = "")
  
  cat("===========================================================================\n")
  cat("Step 5: Extract PC1 scores\n")
  cat("===========================================================================\n")
  cat("Each sample receives a PC1 score, which is its projection onto\n")
  cat("the first principal axis. The scores are centered around 0.\n\n")
  
  cat("--- Sample PC1 values ---\n")
  cat(paste(b1, ": ", paste(sprintf("%.4f", x), collapse = ", "), "\n"))
  cat(paste(b2, ": ", paste(sprintf("%.4f", y), collapse = ", "), "\n\n"))
  
  cat("===========================================================================\n")
  cat("Step 6: Independent Samples t‑test on PC1\n")
  cat("===========================================================================\n")
  cat("Batch 1 (", b1, "): n1 = ", n1, "\n", sep = "")
  cat("Batch 2 (", b2, "): n2 = ", n2, "\n\n", sep = "")
  
  cat("--- Means ---\n")
  cat(sprintf("mean1 = %.4f\n", mean1))
  cat(sprintf("mean2 = %.4f\n\n", mean2))
  
  cat("--- Variances ---\n")
  cat(sprintf("var1 = %.4f\n", var1))
  cat(sprintf("var2 = %.4f\n\n", var2))
  
  cat("--- Pooled standard deviation ---\n")
  cat(sprintf("sp = sqrt[((n1-1)*var1 + (n2-1)*var2) / (n1+n2-2)]\n"))
  cat(sprintf("sp = sqrt[((%d-1)*%.4f + (%d-1)*%.4f) / (%d+%d-2)] = %.4f\n\n", n1, var1, n2, var2, n1, n2, sp))
  
  cat("--- Standard error of difference ---\n")
  cat(sprintf("SE = sp * sqrt(1/n1 + 1/n2) = %.4f * sqrt(1/%d + 1/%d) = %.4f\n\n", sp, n1, n2, se))
  
  cat("--- t‑statistic ---\n")
  cat(sprintf("t = (mean1 - mean2) / SE = (%.4f - %.4f) / %.4f = %.4f\n\n", mean1, mean2, se, t_stat))
  
  cat("--- Degrees of freedom ---\n")
  cat(sprintf("df = n1 + n2 - 2 = %d + %d - 2 = %d\n\n", n1, n2, df))
  
  cat("--- p‑value (two‑tailed) ---\n")
  cat(sprintf("p = 2 * P(T > |t|) = %.4f\n\n", p_val))
  
  if (p_val < 0.05) {
    cat("Conclusion: Significant batch effect detected (p < 0.05).\n")
    cat("Recommendation: Enable batch correction.\n")
  } else {
    cat("Conclusion: No significant batch effect detected (p ≥ 0.05).\n")
    cat("Recommendation: Batch correction is not necessary.\n")
  }
  
  message("[DEBUG] batch_verification_details: t=", round(t_stat,4), " p=", format.pval(p_val, digits=4))
})

# ============ 逐步可视化 ============

# 1. 原始强度分布（所有过滤后样本）
output$batch_viz_raw_hist <- renderPlot({
  data <- batch_verification_data()
  req(data)
  raw_vals <- as.vector(as.matrix(data$raw_filtered))
  raw_vals <- raw_vals[is.finite(raw_vals)]
  df <- data.frame(Value = raw_vals)
  ggplot(df, aes(x = Value)) +
    geom_histogram(bins = 50, fill = "#3498db", alpha = 0.7, boundary = 0) +
    labs(title = "Raw Intensity (after filtering)", x = "Intensity", y = "Frequency") +
    theme_bw()
})

# 2. log2 转换后分布
output$batch_viz_log_hist <- renderPlot({
  data <- batch_verification_data()
  req(data)
  log_vals <- as.vector(data$log_mat)
  log_vals <- log_vals[is.finite(log_vals)]
  df <- data.frame(Value = log_vals)
  ggplot(df, aes(x = Value)) +
    geom_histogram(bins = 50, fill = "#2ecc71", alpha = 0.7, boundary = 0) +
    labs(title = "log2(Intensity + 1)", x = "log2(Intensity)", y = "Frequency") +
    theme_bw()
})

# 3. PCA 得分图（PC1 vs PC2，按批次着色）
output$batch_viz_pca <- renderPlot({
  data <- batch_verification_data()
  req(data)
  pca_scores <- data$pca_scores
  var_explained <- round(summary(data$pca_obj)$importance[2, 1:2] * 100, 1)
  ggplot(pca_scores, aes(x = PC1, y = PC2, color = Batch)) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(type = "norm", level = 0.95) +
    labs(title = "PCA Score Plot",
         subtitle = paste0("PC1: ", var_explained[1], "%, PC2: ", var_explained[2], "%"),
         x = paste0("PC1 (", var_explained[1], "%)"),
         y = paste0("PC2 (", var_explained[2], "%)")) +
    theme_bw() + theme(legend.position = "bottom")
})

# 4. PC1 分组箱线图
output$batch_viz_pc1_box <- renderPlot({
  data <- batch_verification_data()
  req(data)
  pca_scores <- data$pca_scores
  ggplot(pca_scores, aes(x = Batch, y = PC1, fill = Batch)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.1, height = 0, size = 2, alpha = 0.6) +
    labs(title = "PC1 by Batch", y = "PC1 score") +
    theme_bw() + theme(legend.position = "none")
})

# ---------- 冲突检查 ----------
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

# ---------- 批次校正执行（使用 preprocessing_core.R 中定义的对象） ----------
output$batch_correction_performed <- reactive({
  # 安全检查：确保 preprocessing_params 存在且可访问
  tryCatch({
    !is.null(processed_data()) && isTRUE(preprocessing_params$batch_performed)
  }, error = function(e) {
    message("[WARN] batch_correction_performed: could not access preprocessing_params: ", e$message)
    FALSE
  })
})
outputOptions(output, "batch_correction_performed", suspendWhenHidden = FALSE)

batch_comparison_pca <- reactive({
  req(processed_data())
  tryCatch({
    if (!exists("preprocessing_params") || !isTRUE(preprocessing_params$batch_performed)) {
      message("[DEBUG] batch_comparison_pca: batch not performed or preprocessing_params missing")
      return(NULL)
    }
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
  }, error = function(e) {
    message("[ERROR] batch_comparison_pca: ", e$message)
    NULL
  })
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

message("[DEBUG] preprocessing_batch.R fully loaded (conflict resolved, all reactive access safe).")