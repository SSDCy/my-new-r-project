# server/data_quality_plots.R
message("[DEBUG] data_quality_plots.R loaded - optimized user-friendly messages")

# 辅助：安全的 validate
`%then%` <- function(a, b) { if (a) b else TRUE }
validate_condition <- function(condition, message) {
  if (!condition) {
    validate(message)
  }
}

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
  validate_condition(!is.null(dq_expr_matrix()), "Please upload expression data first.")
  tryCatch(dq_quality_score()$score, error = function(e) paste("Error:", e$message))
})

output$dq_missing_rate <- renderText({
  validate_condition(!is.null(dq_expr_matrix()), "Please upload expression data first.")
  tryCatch(paste0(dq_quality_score()$details$missing_ratio, "%"), error = function(e) paste("Error:", e$message))
})
output$dq_missing_score_frac <- renderText({
  validate_condition(!is.null(dq_expr_matrix()), "Please upload expression data first.")
  tryCatch(paste0(dq_quality_score()$details$missing_score, " / 30"), error = function(e) paste("Error:", e$message))
})
output$dq_consistency_score_frac <- renderText({
  validate_condition(!is.null(dq_expr_matrix()), "Please upload expression data first.")
  tryCatch(paste0(dq_quality_score()$details$consistency_score, " / 40"), error = function(e) paste("Error:", e$message))
})
output$dq_protein_score_frac <- renderText({
  validate_condition(!is.null(dq_expr_matrix()), "Please upload expression data first.")
  tryCatch(paste0(dq_quality_score()$details$protein_score, " / 30"), error = function(e) paste("Error:", e$message))
})

# ==================== 样本选择（缺失值热图） ====================
observeEvent(rv$sample_names, {
  message("[DEBUG] updating heatmap_sample_select choices")
  all_samples <- rv$sample_names
  if (length(all_samples) > 0) {
    updateSelectizeInput(session, "heatmap_sample_select",
                         choices = all_samples,
                         selected = all_samples,
                         server = TRUE)
    updateSelectizeInput(session, "missing_cor_sample_select",
                         choices = all_samples,
                         selected = all_samples,
                         server = TRUE)
  }
}, ignoreNULL = TRUE, once = FALSE)

output$heatmap_group_buttons_ui <- renderUI({
  req(rv$sample_info)
  if (!"Group" %in% colnames(rv$sample_info)) return(NULL)
  groups <- unique(rv$sample_info$Group)
  lapply(groups, function(g) {
    actionButton(inputId = paste0("heatmap_group_", g), label = g,
                 class = "btn-sm btn-outline-info", icon = icon("filter"))
  })
})

observeEvent(input$heatmap_select_all, {
  req(rv$sample_names)
  updateSelectizeInput(session, "heatmap_sample_select", selected = rv$sample_names)
})
observeEvent(input$heatmap_clear_all, {
  updateSelectizeInput(session, "heatmap_sample_select", selected = character(0))
})

observe({
  group_buttons <- grep("^heatmap_group_", names(input), value = TRUE)
  for (btn in group_buttons) {
    local({
      b <- btn
      observeEvent(input[[b]], {
        group <- sub("^heatmap_group_", "", b)
        req(rv$sample_info)
        si <- rv$sample_info
        si$SampleName <- rownames(si)
        si$ShortName <- extract_sample_names(si$SampleName)
        samples_in_group <- si$ShortName[si$Group == group]
        current_sel <- input$heatmap_sample_select
        new_sel <- union(current_sel, samples_in_group)
        updateSelectizeInput(session, "heatmap_sample_select", selected = new_sel)
      }, ignoreInit = TRUE)
    })
  }
})

selected_samples <- reactive({
  sel <- input$heatmap_sample_select
  if (is.null(sel) || length(sel) == 0) return(rv$sample_names)
  return(sel)
})

# ==================== 缺失值热图 ====================
dq_missing_heatmap_plot_obj <- reactive({
  req(dq_expr_matrix())
  mat <- dq_expr_matrix()
  sel_samples <- selected_samples()
  common_s <- intersect(sel_samples, colnames(mat))
  if (length(common_s) == 0) return(NULL)
  
  annot_df <- NULL
  if (!is.null(rv$sample_info) && "Group" %in% colnames(rv$sample_info)) {
    si <- rv$sample_info
    si$SampleName <- rownames(si)
    si$ShortName <- extract_sample_names(si$SampleName)
    idx <- match(common_s, si$ShortName)
    if (any(!is.na(idx))) {
      groups <- si$Group[idx]
      groups[is.na(idx)] <- "Unknown"
      annot_df <- data.frame(Group = factor(groups), row.names = common_s)
    }
  }
  
  if (!is.null(annot_df)) {
    common_s <- rownames(annot_df)[order(annot_df$Group, rownames(annot_df))]
    annot_df <- annot_df[common_s, , drop = FALSE]
  } else {
    common_s <- sort(common_s)
  }
  
  mat_sub <- mat[, common_s, drop = FALSE]
  missing_mat <- (is.na(mat_sub) * 1)
  n_prot <- nrow(missing_mat)
  n_samp <- ncol(missing_mat)
  
  cluster_rows <- if (is.null(input$heatmap_cluster_rows)) FALSE else input$heatmap_cluster_rows
  
  if (cluster_rows && n_prot >= 2) {
    dist_row <- dist(missing_mat, method = "binary")
    hc <- hclust(dist_row, method = "ward.D2")
    row_order <- hc$order
  } else {
    row_order <- seq_len(n_prot)
  }
  
  protein_levels <- rownames(missing_mat)[row_order]
  sample_levels <- common_s
  plot_df <- expand.grid(
    Protein = factor(protein_levels, levels = protein_levels),
    Sample = factor(sample_levels, levels = sample_levels),
    stringsAsFactors = FALSE
  )
  plot_df$Value <- as.vector(missing_mat[row_order, sample_levels])
  
  p <- ggplot(plot_df, aes(x = Sample, y = Protein, fill = factor(Value))) +
    geom_tile() +
    scale_fill_manual(
      values = c("0" = "#3498db", "1" = "#e74c3c"),
      labels = c("0" = "Detected", "1" = "Missing"),
      name = "Status"
    ) +
    labs(title = paste0("Missing Heatmap (", n_prot, " proteins)"), x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  return(p)
})

output$dq_missing_heatmap <- renderPlot({
  validate_condition(!is.null(dq_expr_matrix()), "Please upload expression data first.")
  p <- dq_missing_heatmap_plot_obj()
  if (is.null(p)) {
    plot.new(); text(0.5, 0.5, "Please select at least one sample")
  } else {
    print(p)
  }
})

# ==================== 蛋白缺失率分布直方图 ====================
dq_protein_missing_hist <- reactive({
  req(dq_expr_matrix())
  mat <- dq_expr_matrix()
  sel_s <- selected_samples()
  common_s <- intersect(sel_s, colnames(mat))
  if (length(common_s) == 0) return(NULL)
  mat_sub <- mat[, common_s, drop = FALSE]
  rates <- rowMeans(is.na(mat_sub))
  df <- data.frame(MissingRate = rates)
  ggplot(df, aes(x = MissingRate)) +
    geom_histogram(fill = "#3498db", bins = 30, alpha = 0.8, boundary = 0) +
    labs(title = paste0("Protein Missing Rate Distribution (", length(common_s), " samples)"),
         x = "Missing Rate", y = "Number of Proteins") +
    theme_bw()
})

output$dq_protein_missing_hist <- renderPlot({
  validate_condition(!is.null(dq_expr_matrix()), "Please upload expression data first.")
  req(dq_protein_missing_hist())
  dq_protein_missing_hist()
})
output$download_protein_missing_hist <- downloadHandler(
  filename = function() "protein_missing_hist.png",
  content = function(file) ggsave(file, plot = dq_protein_missing_hist(), width = 6, height = 4, dpi = 150)
)

# ==================== 样本缺失率条形图 ====================
dq_sample_missing_bar <- reactive({
  req(dq_expr_matrix())
  mat <- dq_expr_matrix()
  sel_s <- selected_samples()
  common_s <- intersect(sel_s, colnames(mat))
  if (length(common_s) == 0) return(NULL)
  mat_sub <- mat[, common_s, drop = FALSE]
  missing_pct <- colMeans(is.na(mat_sub)) * 100
  df <- data.frame(Sample = factor(common_s, levels = common_s),
                   MissingPct = missing_pct)
  ggplot(df, aes(x = Sample, y = MissingPct, fill = MissingPct)) +
    geom_col() +
    scale_fill_gradient(low = "#fcbba1", high = "#cb181d", guide = "none") +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    labs(title = paste0("Sample Missing Rate (", length(common_s), " samples)"),
         y = "Missing Percentage (%)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
          axis.title.x = element_blank())
})

output$dq_sample_missing_bar <- renderPlot({
  validate_condition(!is.null(dq_expr_matrix()), "Please upload expression data first.")
  req(dq_sample_missing_bar())
  dq_sample_missing_bar()
}, height = 400)
output$download_sample_missing_bar <- downloadHandler(
  filename = function() "sample_missing_rate.png",
  content = function(file) {
    p <- dq_sample_missing_bar()
    if (!is.null(p)) ggsave(file, plot = p, width = 12, height = 6, dpi = 150)
  }
)

# ==================== 缺失值相关性热图 ====================
missing_cor_subset <- reactiveVal(NULL)

observeEvent(rv$sample_names, {
  updateSelectizeInput(session, "missing_cor_sample_select",
                       choices = rv$sample_names,
                       selected = rv$sample_names,
                       server = TRUE)
})

observeEvent(input$missing_cor_subset_go, {
  sel <- input$missing_cor_sample_select
  if (is.null(sel) || length(sel) < 2) {
    showNotification("请至少选择两个样本以生成子热图", type = "warning")
    return()
  }
  missing_cor_subset(sel)
  showNotification(paste0("已切换至选中样本（", length(sel), " 个）的子热图"), type = "message")
})

observeEvent(input$missing_cor_reset, {
  missing_cor_subset(NULL)
  updateSelectizeInput(session, "missing_cor_sample_select", selected = rv$sample_names)
  showNotification("已恢复为全局热图", type = "message")
})

dq_missing_cor_plot_obj <- reactive({
  req(dq_expr_matrix())
  mat <- dq_expr_matrix()
  missing_mat <- is.na(mat) * 1
  
  subset_samples <- missing_cor_subset()
  if (is.null(subset_samples) || length(subset_samples) == 0) {
    common_s <- intersect(rv$sample_names, colnames(missing_mat))
  } else {
    common_s <- intersect(subset_samples, colnames(missing_mat))
  }
  if (length(common_s) < 2) return(NULL)
  
  missing_sub <- missing_mat[, common_s, drop = FALSE]
  cor_mat <- cor(missing_sub, use = "pairwise.complete.obs")
  cor_mat[is.na(cor_mat)] <- 0
  
  min_cor <- min(cor_mat, na.rm = TRUE)
  max_cor <- max(cor_mat, na.rm = TRUE)
  if (abs(max_cor - min_cor) < 1e-6) { min_cor <- min_cor - 0.01; max_cor <- max_cor + 0.01 }
  
  legend_breaks <- pretty(c(min_cor, max_cor), n = 5)
  legend_breaks <- legend_breaks[legend_breaks >= min_cor & legend_breaks <= max_cor]
  legend_labels <- sprintf("%.2f", legend_breaks)
  
  cor_dist <- as.dist(1 - cor_mat)
  
  annotation_col <- NULL
  annotation_colors <- NULL
  if (!is.null(rv$sample_info) && "Group" %in% colnames(rv$sample_info)) {
    si <- rv$sample_info
    si$ShortName <- extract_sample_names(rownames(si))
    common <- intersect(colnames(cor_mat), si$ShortName)
    if (length(common) > 0) {
      group_vec <- si$Group[match(common, si$ShortName)]
      annotation_col <- data.frame(Group = factor(group_vec), row.names = common)
      groups <- unique(group_vec)
      annotation_colors <- list(Group = get_group_colors(groups))
    }
  }
  
  main_title <- "缺失值相关性热图（Pearson 相关系数）"
  if (!is.null(subset_samples)) main_title <- paste0(main_title, " 子集")
  
  pheatmap::pheatmap(cor_mat,
                     main = main_title,
                     color = colorRampPalette(c("blue", "white", "red"))(100),
                     breaks = seq(min_cor, max_cor, length.out = 101),
                     legend_breaks = legend_breaks,
                     legend_labels = legend_labels,
                     clustering_distance_rows = cor_dist,
                     clustering_distance_cols = cor_dist,
                     clustering_method = "ward.D2",
                     show_rownames = TRUE,
                     show_colnames = TRUE,
                     fontsize_row = 10,
                     fontsize_col = 8,
                     angle_col = 45,
                     treeheight_row = 70,
                     treeheight_col = 70,
                     margins = c(12, 8),
                     annotation_col = annotation_col,
                     annotation_colors = annotation_colors,
                     legend = TRUE,
                     silent = TRUE)
})

output$dq_missing_cor_plot <- renderPlot({
  validate_condition(!is.null(dq_expr_matrix()), "Please upload expression data first.")
  obj <- dq_missing_cor_plot_obj()
  if (!is.null(obj)) {
    grid::grid.newpage()
    grid::grid.draw(obj$gtable)
  } else {
    plot.new(); text(0.5, 0.5, "样本数不足")
  }
})

output$download_missing_cor <- downloadHandler(
  filename = function() "missing_correlation.png",
  content = function(file) {
    png(file, width = 1600, height = 1400, res = 150)
    obj <- dq_missing_cor_plot_obj()
    if (!is.null(obj)) grid::grid.draw(obj$gtable)
    dev.off()
  }
)

output$download_missing_cor_matrix <- downloadHandler(
  filename = function() { paste0("missing_correlation_matrix_", Sys.Date(), ".csv") },
  content = function(file) {
    req(dq_expr_matrix())
    mat <- dq_expr_matrix()
    missing_mat <- is.na(mat) * 1
    subset_samples <- missing_cor_subset()
    if (is.null(subset_samples)) common_s <- intersect(rv$sample_names, colnames(missing_mat))
    else common_s <- intersect(subset_samples, colnames(missing_mat))
    if (length(common_s) < 2) { showNotification("Not enough samples", type = "error"); return() }
    missing_sub <- missing_mat[, common_s, drop = FALSE]
    cor_mat <- cor(missing_sub, use = "pairwise.complete.obs")
    cor_mat[is.na(cor_mat)] <- 0
    write.csv(cor_mat, file, row.names = TRUE)
  }
)

observeEvent(input$help_missing_cor, {
  showModal(modalDialog(
    title = "缺失值相关性",
    "基于样本缺失 0/1 矩阵的 Pearson 相关系数。聚类使用 (1 - 相关性) 作为距离。",
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
  
  group_vec <- NULL
  if (!is.null(rv$sample_info) && "Group" %in% colnames(rv$sample_info)) {
    si <- rv$sample_info
    si$ShortName <- extract_sample_names(rownames(si))
    idx <- match(df$Sample, si$ShortName)
    if (any(!is.na(idx))) {
      group_vec <- si$Group[idx]
      group_vec[is.na(idx)] <- "Unknown"
    }
  }
  
  if (!is.null(group_vec)) {
    df$Group <- group_vec
    p <- ggplot(df, aes(x = Sample, y = Log2Intensity, fill = Group)) +
      geom_boxplot(outlier.size = 1, alpha = 0.7) +
      labs(title = "Protein Intensity Distribution (log2-transformed)", y = "log2(Intensity)") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), legend.position = "right")
  } else {
    p <- ggplot(df, aes(x = Sample, y = Log2Intensity)) +
      geom_boxplot(fill = "#3498db", alpha = 0.7, outlier.size = 1) +
      labs(title = "Protein Intensity Distribution (log2-transformed)", y = "log2(Intensity)") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
  }
  return(p)
})

output$dq_intensity_dist_plot <- renderPlot({
  validate_condition(!is.null(dq_expr_matrix()), "Please upload expression data first.")
  dq_intensity_plot()
})
output$download_intensity <- downloadHandler(
  filename = function() "intensity_distribution.png",
  content = function(file) ggsave(file, plot = dq_intensity_plot(), width = 10, height = 6, dpi = 150)
)
output$download_intensity_data <- downloadHandler(
  filename = function() { paste0("Intensity_Data_", Sys.Date(), ".xlsx") },
  content = function(file) {
    req(dq_expr_matrix())
    mat <- dq_expr_matrix()
    log_mat <- log2(mat + 1)
    log_df <- as.data.frame(log_mat)
    sample_ids <- colnames(mat)
    pca_full <- dq_pca_full()
    if (!is.null(pca_full)) outlier_samples <- pca_full$scores$Sample[pca_full$scores$Outlier == "Outlier"]
    else outlier_samples <- character(0)
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
      "Sheet 'log2_intensity_matrix': log2(Intensity+1) transformed intensity matrix.",
      "Sheet 'sample_intensity_stats': summary statistics per sample."
    ))
    openxlsx::addWorksheet(wb, "README")
    openxlsx::writeData(wb, "README", readme)
    openxlsx::addWorksheet(wb, "log2_intensity_matrix")
    openxlsx::writeData(wb, "log2_intensity_matrix", log_df)
    openxlsx::addWorksheet(wb, "sample_intensity_stats")
    openxlsx::writeData(wb, "sample_intensity_stats", stats_df)
    openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
  }
)
observeEvent(input$help_intensity, {
  showModal(modalDialog(
    title = "Protein Intensity Distribution",
    "箱线图展示每个样本中蛋白强度的分布（log2 转换）。",
    easyClose = TRUE, footer = modalButton("关闭")
  ))
})

# ==================== 样本相关性热图 ====================
dq_cor_heatmap_plot_obj <- reactive({
  req(dq_expr_matrix())
  message("[DEBUG] dq_cor_heatmap_plot_obj: computing sample correlation")
  cor_mat <- calculate_sample_correlation(dq_expr_matrix())
  if (is.null(cor_mat)) {
    message("[DEBUG] dq_cor_heatmap_plot_obj: correlation matrix is NULL")
    return(NULL)
  }
  cor_mat[is.na(cor_mat)] <- 0
  
  ann_col <- NULL
  ann_colors <- NULL
  subcluster_info <- NULL
  
  if (!is.null(rv$sample_info) && "Group" %in% colnames(rv$sample_info)) {
    sample_info_short <- rv$sample_info
    rownames(sample_info_short) <- extract_sample_names(rownames(sample_info_short))
    common_samples <- intersect(colnames(cor_mat), rownames(sample_info_short))
    if (length(common_samples) > 0) {
      group_vec <- sample_info_short[common_samples, "Group"]
      ann_col <- data.frame(Group = group_vec, row.names = common_samples)
      
      if (ncol(cor_mat) >= 3) {
        cor_sub <- cor_mat[common_samples, common_samples, drop = FALSE]
        cor_dist <- as.dist(1 - cor_sub)
        hc_col <- hclust(cor_dist, method = "ward.D2")
        k <- ifelse(length(common_samples) >= 6, 3, 2)
        subcluster <- cutree(hc_col, k = k)
        while (any(table(subcluster) < 2) && k > 1) { k <- k - 1; subcluster <- cutree(hc_col, k = k) }
        subcluster_label <- paste0("SubClust", subcluster)
        ann_col$Subcluster <- factor(subcluster_label, levels = unique(subcluster_label))
        subcluster_info <- subcluster
      }
      
      groups <- unique(group_vec)
      group_colors <- get_group_colors(groups)
      ann_colors <- list(Group = group_colors)
      if (!is.null(subcluster_info)) {
        n_sub <- length(unique(subcluster_info))
        sub_colors <- RColorBrewer::brewer.pal(min(n_sub, 8), "Set2")[1:n_sub]
        names(sub_colors) <- unique(subcluster_label)
        ann_colors$Subcluster <- sub_colors
      }
      
      cor_mat <- cor_mat[common_samples, common_samples, drop = FALSE]
    }
  }
  
  min_cor <- min(cor_mat, na.rm = TRUE)
  max_cor <- max(cor_mat, na.rm = TRUE)
  message("[DEBUG] dq_cor_heatmap_plot_obj: cor range (actual): [", min_cor, ", ", max_cor, "]")
  
  limit <- 1
  my_colors <- colorRampPalette(c("blue", "white", "red"))(255)
  breaks <- seq(-limit, limit, length.out = 256)
  
  legend_breaks <- c(-1, -0.5, 0, 0.5, 1)
  legend_labels <- c("-1.0", "-0.5", "0.0", "0.5", "1.0")
  
  cor_dist <- as.dist(1 - cor_mat)
  
  pheatmap::pheatmap(cor_mat,
                     main = "样本相关性热图（基于 top 500 高变异蛋白的 log2 强度 Pearson 相关）",
                     color = my_colors,
                     breaks = breaks,
                     legend_breaks = legend_breaks,
                     legend_labels = legend_labels,
                     clustering_distance_rows = cor_dist,
                     clustering_distance_cols = cor_dist,
                     clustering_method = "ward.D2",
                     show_rownames = TRUE, show_colnames = TRUE,
                     fontsize_row = 9, fontsize_col = 9,
                     angle_col = 45,
                     annotation_col = ann_col, annotation_colors = ann_colors,
                     na_col = "grey",
                     silent = TRUE)
})

output$dq_cor_heatmap <- renderPlot({
  validate_condition(!is.null(dq_expr_matrix()), "Please upload expression data first.")
  obj <- dq_cor_heatmap_plot_obj()
  if (is.null(obj)) { plot.new(); text(0.5, 0.5, "Not enough data") }
  else { grid::grid.newpage(); grid::grid.draw(obj$gtable) }
})

output$download_cor_heatmap <- downloadHandler(
  filename = function() "sample_correlation.png",
  content = function(file) { png(file, width = 900, height = 700, res = 150); obj <- dq_cor_heatmap_plot_obj(); if (!is.null(obj)) grid::grid.draw(obj$gtable); dev.off() }
)

output$download_cor_matrix <- downloadHandler(
  filename = function() { paste0("correlation_matrix_", Sys.Date(), ".xlsx") },
  content = function(file) {
    cor_mat <- calculate_sample_correlation(dq_expr_matrix())
    if (is.null(cor_mat)) { showNotification("相关性矩阵不可用", type = "error"); return() }
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "相关性矩阵")
    openxlsx::writeData(wb, "相关性矩阵", cor_mat, rowNames = TRUE)
    openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
  }
)

observeEvent(input$help_cor_heatmap, {
  showModal(modalDialog(title = "样本相关性热图", "基于高变异蛋白计算的样本间 Pearson 相关性。颜色范围固定为 [-1, 1]。", easyClose = TRUE, footer = modalButton("关闭")))
})

# ==================== 缺失值类型诊断 ====================
dq_missing_type_data <- reactive({
  req(dq_expr_matrix())
  mat <- dq_expr_matrix()
  missing_ratio <- rowMeans(is.na(mat))
  mean_intensity <- rowMeans(mat, na.rm = TRUE)
  log2_mean <- log2(mean_intensity + 1)
  
  group <- ifelse(missing_ratio == 0, "Complete (0%)",
                  ifelse(missing_ratio <= 0.5, "Partial (0-50%)", "High (>50%)"))
  df <- data.frame(
    Protein = rownames(mat),
    MissingRatio = missing_ratio,
    Log2MeanIntensity = log2_mean,
    Group = factor(group, levels = c("Complete (0%)", "Partial (0-50%)", "High (>50%)")),
    stringsAsFactors = FALSE
  )
  df
})

dq_missing_type_plot_obj <- reactive({
  df <- dq_missing_type_data()
  req(df)
  ggplot(df, aes(x = Log2MeanIntensity, fill = Group, color = Group)) +
    geom_density(alpha = 0.4, linewidth = 1) +
    scale_fill_manual(values = c("Complete (0%)" = "#2ecc71", "Partial (0-50%)" = "#f1c40f", "High (>50%)" = "#e74c3c")) +
    scale_color_manual(values = c("Complete (0%)" = "#27ae60", "Partial (0-50%)" = "#f39c12", "High (>50%)" = "#c0392b")) +
    labs(title = "Missing Value Type Diagnosis: Intensity Distribution by Missing Rate",
         subtitle = "左移的红色曲线提示 MNAR（非随机缺失，与低表达相关）",
         x = "log2(Mean Intensity + 1)", y = "Density") +
    theme_bw() + theme(legend.position = "right")
})

output$dq_missing_type_plot <- renderPlot({
  validate_condition(!is.null(dq_expr_matrix()), "Please upload expression data first.")
  dq_missing_type_plot_obj()
})

output$download_missing_type <- downloadHandler(
  filename = function() "missing_type_diagnosis.png",
  content = function(file) ggsave(file, plot = dq_missing_type_plot_obj(), width = 8, height = 5, dpi = 150)
)

output$missing_type_summary <- renderText({
  df <- dq_missing_type_data()
  req(df)
  high_missing <- df[df$Group == "High (>50%)", ]
  complete <- df[df$Group == "Complete (0%)", ]
  if (nrow(high_missing) < 5 || nrow(complete) < 5) return("Insufficient data for diagnosis.")
  
  median_high <- median(high_missing$Log2MeanIntensity, na.rm = TRUE)
  median_complete <- median(complete$Log2MeanIntensity, na.rm = TRUE)
  shift <- median_complete - median_high
  if (shift > 1.5) {
    "诊断结果：**高度提示 MNAR（非随机缺失）**。建议使用最小值或分位数填补。"
  } else if (shift > 0.5) {
    "诊断结果：**中等提示 MNAR**。建议结合领域知识选择填补方法。"
  } else {
    "诊断结果：**提示 MAR（随机缺失）为主**。推荐使用 KNN 或 PPCA 填补。"
  }
})

observeEvent(input$help_missing_type, {
  showModal(modalDialog(
    title = "缺失值类型诊断（MNAR vs MAR）",
    "通过比较不同缺失率蛋白的强度分布，判断缺失值的主要类型。",
    easyClose = TRUE, footer = modalButton("关闭")
  ))
})

# ==================== 缺失值定量统计 ====================
dq_group_missing_data <- reactive({
  req(dq_expr_matrix())
  mat <- dq_expr_matrix()
  sample_missing_rate <- colMeans(is.na(mat)) * 100
  sample_names_short <- colnames(mat)
  group <- rep("Unassigned", length(sample_names_short))
  
  if (!is.null(rv$sample_info) && "Group" %in% colnames(rv$sample_info)) {
    si <- rv$sample_info
    si$ShortName <- extract_sample_names(rownames(si))
    idx <- match(sample_names_short, si$ShortName)
    if (any(!is.na(idx))) {
      group[!is.na(idx)] <- si$Group[idx[!is.na(idx)]]
    }
  }
  
  data.frame(Sample = sample_names_short, MissingRate = sample_missing_rate, Group = factor(group), stringsAsFactors = FALSE)
})

output$dq_group_missing_boxplot <- renderPlot({
  validate_condition(!is.null(dq_expr_matrix()), "Please upload expression data first.")
  df <- dq_group_missing_data()
  req(df)
  if (all(df$Group == "Unassigned")) {
    ggplot(df, aes(x = Group, y = MissingRate)) +
      geom_boxplot(fill = "#3498db", alpha = 0.7) +
      labs(title = "Sample Missing Rate (No group information)", y = "Missing Rate (%)") +
      theme_bw()
  } else {
    ggplot(df, aes(x = Group, y = MissingRate, fill = Group)) +
      geom_boxplot(alpha = 0.7) +
      labs(title = "Missing Rate Distribution by Group", y = "Missing Rate (%)") +
      theme_bw() + theme(legend.position = "none")
  }
})

output$dq_group_missing_table <- renderTable({
  req(dq_group_missing_data())
  df <- dq_group_missing_data()
  stats <- df %>%
    group_by(Group) %>%
    summarise(
      N = n(),
      Mean = round(mean(MissingRate), 2),
      SD = round(sd(MissingRate), 2),
      Median = round(median(MissingRate), 2),
      Min = round(min(MissingRate), 2),
      Max = round(max(MissingRate), 2),
      .groups = "drop"
    )
  stats
}, striped = TRUE, bordered = TRUE, width = "100%")

output$dq_group_missing_test <- renderPrint({
  df <- dq_group_missing_data()
  req(df)
  groups_valid <- unique(df$Group[df$Group != "Unassigned"])
  if (length(groups_valid) < 2) {
    cat("Not enough groups for statistical comparison.")
    return()
  }
  df_valid <- df[df$Group %in% groups_valid, ]
  df_valid$Group <- factor(df_valid$Group)
  test_result <- tryCatch(kruskal.test(MissingRate ~ Group, data = df_valid), error = function(e) NULL)
  if (is.null(test_result)) {
    cat("Statistical test failed.")
    return()
  }
  pval <- test_result$p.value
  cat("Kruskal-Wallis test p-value =", format.pval(pval, digits = 3), "\n")
  if (pval < 0.05) cat("Significant difference among groups.") else cat("No significant difference.")
})

# ==================== PCA 分析 ====================
dq_pca_full <- reactive({
  message("[DEBUG] dq_pca_full: starting PCA computation")
  
  expr <- tryCatch({
    mat <- expression_data()
    message("[DEBUG] dq_pca_full: expression_data() returned dim=", nrow(mat), "x", ncol(mat))
    mat
  }, error = function(e) {
    message("[DEBUG] dq_pca_full: expression_data() error - ", e$message)
    NULL
  })
  
  if (is.null(expr)) {
    message("[DEBUG] dq_pca_full: expression_data() is NULL, abort")
    return(NULL)
  }
  
  proc_mat <- tryCatch(get_analysis_matrix(), error = function(e) NULL)
  if (!is.null(proc_mat)) {
    message("[DEBUG] dq_pca_full: using preprocessed data (dim=", nrow(proc_mat), "x", ncol(proc_mat), ")")
    expr <- proc_mat
    data_source <- "Preprocessed"
  } else {
    message("[DEBUG] dq_pca_full: using raw expression_data")
    data_source <- "Raw"
  }
  
  sample_short <- extract_sample_names(colnames(expr))
  colnames(expr) <- sample_short
  message("[DEBUG] dq_pca_full: columns standardized, first 3: ", paste(head(sample_short, 3), collapse = ", "))
  
  filled <- tryCatch({
    if (requireNamespace("impute", quietly = TRUE)) {
      message("[DEBUG] dq_pca_full: running KNN imputation")
      impute_missing_values(expr, method = "knn")
    } else {
      expr[is.na(expr)] <- 1e-4
      expr
    }
  }, error = function(e) {
    message("[DEBUG] dq_pca_full: KNN failed, min fill - ", e$message)
    expr[is.na(expr)] <- 1e-4
    expr
  })
  
  log_expr <- log2(filled + 1)
  row_vars <- apply(log_expr, 1, var)
  log_expr <- log_expr[row_vars > 1e-12, , drop = FALSE]
  row_unique <- apply(log_expr, 1, function(x) length(unique(x)))
  log_expr <- log_expr[row_unique > 1, , drop = FALSE]
  
  message("[DEBUG] dq_pca_full: rows after filtering: ", nrow(log_expr))
  
  if (nrow(log_expr) < 2) {
    message("[DEBUG] dq_pca_full: too few variable rows")
    return(NULL)
  }
  
  pca <- tryCatch(prcomp(t(log_expr), scale. = TRUE), error = function(e) {
    message("[DEBUG] dq_pca_full: prcomp error - ", e$message)
    NULL
  })
  if (is.null(pca)) return(NULL)
  
  variance <- pca$sdev^2 / sum(pca$sdev^2) * 100
  cum_variance <- cumsum(variance)
  message(sprintf("[DEBUG] dq_pca_full: PC1=%.1f%%, PC2=%.1f%%, cumulative=%.1f%%", 
                  variance[1], variance[2], cum_variance[2]))
  
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
      scores$Group <- "All"; scores$Batch <- NA
    }
  } else {
    scores$Group <- "All"; scores$Batch <- NA
  }
  
  outliers <- get_outlier_samples(list(pca_df = scores, pc1_var = variance[1], pc2_var = variance[2]))
  scores$Outlier <- ifelse(scores$Sample %in% outliers, "Outlier", "Normal")
  
  message("[DEBUG] dq_pca_full: PCA completed, data_source = ", data_source)
  list(pca = pca, scores = scores, variance = variance, cum_variance = cum_variance, loadings = pca$rotation, data_source = data_source)
})

output$pca_data_source_note <- renderUI({
  pca_full <- dq_pca_full()
  if (is.null(pca_full)) {
    return(div(style = "margin-bottom: 10px; color: #e74c3c;", "PCA calculation failed. Please check your data or run preprocessing."))
  }
  if (pca_full$data_source == "Preprocessed") {
    div(style = "margin-bottom: 10px; color: #27ae60; font-weight: bold;",
        icon("check-circle"), " PCA data source: Preprocessed data")
  } else {
    div(style = "margin-bottom: 10px; color: #e67e22; font-weight: bold;",
        icon("exclamation-triangle"), " PCA data source: Raw data (KNN imputed). Running preprocessing will automatically switch to preprocessed data.")
  }
})

pca_group_plot_obj <- reactive({
  pca_full <- dq_pca_full()
  if (is.null(pca_full)) return(NULL)
  scores <- pca_full$scores
  scores <- scores[!is.na(scores$Group), ]
  if (nrow(scores) == 0) return(NULL)
  
  group_colors <- c("Control" = "#FF69B4", "Treatment" = "#00CED1")
  all_groups <- unique(scores$Group)
  missing_colors <- setdiff(all_groups, names(group_colors))
  if (length(missing_colors) > 0) {
    extra_colors <- rainbow(length(missing_colors))
    names(extra_colors) <- missing_colors
    group_colors <- c(group_colors, extra_colors)
  }
  pc1_label <- sprintf("PC1 (%.1f%%, cum. %.1f%%)", pca_full$variance[1], pca_full$cum_variance[1])
  pc2_label <- sprintf("PC2 (%.1f%%, cum. %.1f%%)", pca_full$variance[2], pca_full$cum_variance[2])
  
  ggplot(scores, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = group_colors) +
    labs(title = "PCA by Group",
         x = pc1_label,
         y = pc2_label) +
    theme_bw() + theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
})

output$dq_pca_group_plot <- renderPlot({
  validate_condition(!is.null(dq_expr_matrix()), "Please upload expression data first.")
  p <- pca_group_plot_obj()
  if (is.null(p)) {
    plot.new(); text(0.5, 0.5, "PCA not available")
  } else {
    print(p)
  }
})
output$download_pca_group <- downloadHandler(
  filename = function() "pca_group.png",
  content = function(file) {
    p <- pca_group_plot_obj()
    if (!is.null(p)) ggsave(file, plot = p, width = 8, height = 6, dpi = 150)
    else showNotification("PCA not available", type = "error")
  }
)

pca_batch_plot_obj <- reactive({
  pca_full <- dq_pca_full()
  if (is.null(pca_full)) return(NULL)
  scores <- pca_full$scores
  if (all(is.na(scores$Batch))) return(NULL)
  scores <- scores[!is.na(scores$Batch), ]
  if (nrow(scores) == 0) return(NULL)
  batch_colors <- c("Batch1" = "#E41A1C", "Batch2" = "#00CED1")
  all_batches <- unique(scores$Batch)
  if (length(all_batches) > 2) batch_colors <- setNames(rainbow(length(all_batches)), all_batches)
  
  pc1_label <- sprintf("PC1 (%.1f%%, cum. %.1f%%)", pca_full$variance[1], pca_full$cum_variance[1])
  pc2_label <- sprintf("PC2 (%.1f%%, cum. %.1f%%)", pca_full$variance[2], pca_full$cum_variance[2])
  
  ggplot(scores, aes(x = PC1, y = PC2, color = Batch)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = batch_colors) +
    labs(title = "PCA by Batch",
         x = pc1_label,
         y = pc2_label) +
    theme_bw() + theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
})

output$dq_pca_batch_plot <- renderPlot({
  validate_condition(!is.null(dq_expr_matrix()), "Please upload expression data first.")
  p <- pca_batch_plot_obj()
  if (is.null(p)) { plot.new(); text(0.5, 0.5, "Batch info not available") }
  else { print(p) }
})
output$download_pca_batch <- downloadHandler(
  filename = function() "pca_batch.png",
  content = function(file) {
    p <- pca_batch_plot_obj()
    if (!is.null(p)) ggsave(file, plot = p, width = 8, height = 6, dpi = 150)
    else showNotification("PCA batch not available", type = "error")
  }
)

observeEvent(input$help_pca_group, {
  showModal(modalDialog(title = "PCA by Group", "按实验分组着色，轴标签包含方差解释率和累计解释率。", easyClose = TRUE, footer = modalButton("关闭")))
})
observeEvent(input$help_pca_batch, {
  showModal(modalDialog(title = "PCA by Batch", "按实验批次着色，轴标签包含方差解释率和累计解释率。", easyClose = TRUE, footer = modalButton("关闭")))
})