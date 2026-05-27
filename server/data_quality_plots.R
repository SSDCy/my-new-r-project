# server/data_quality_plots.R
message("[DEBUG] data_quality_plots.R loaded - removed download_missing_cor_stats export")

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
  tryCatch(paste0(dq_quality_score()$details$consistency_score, " / 40"), error = function(e) paste("Error:", e$message))
})
output$dq_protein_score_frac <- renderText({
  req(dq_expr_matrix())
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

# ==================== 缺失值热图（ggplot2） ====================
dq_missing_heatmap_plot_obj <- reactive({
  req(dq_expr_matrix())
  mat <- dq_expr_matrix()
  sel_samples <- selected_samples()
  common_s <- intersect(sel_samples, colnames(mat))
  if (length(common_s) == 0) return(NULL)
  
  if (!is.null(rv$sample_info) && "Group" %in% colnames(rv$sample_info)) {
    si <- rv$sample_info
    si$SampleName <- rownames(si)
    si$ShortName <- extract_sample_names(si$SampleName)
    si <- si[si$ShortName %in% common_s, ]
    if (nrow(si) > 0) {
      si <- si[order(si$Group, si$ShortName), ]
      common_s <- si$ShortName
    }
  } else {
    common_s <- sort(common_s)
  }
  
  mat_sub <- mat[, common_s, drop = FALSE]
  missing_mat <- (is.na(mat_sub) * 1)
  n_prot <- nrow(missing_mat)
  n_samp <- ncol(missing_mat)
  
  message(sprintf("[DEBUG] heatmap FULL matrix: %d proteins x %d samples", n_prot, n_samp))
  
  cluster_rows <- if (is.null(input$heatmap_cluster_rows)) FALSE else input$heatmap_cluster_rows
  message("[DEBUG] cluster_rows: ", cluster_rows)
  
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
  
  message("[DEBUG] ggplot2 heatmap generated")
  return(p)
})

output$dq_missing_heatmap <- renderPlot({
  p <- dq_missing_heatmap_plot_obj()
  if (is.null(p)) {
    plot.new(); text(0.5, 0.5, "Please select at least one sample")
  } else {
    print(p)
  }
})

# ==================== 缺失值热图导出 ====================
output$download_missing_heatmap <- downloadHandler(
  filename = function() "missing_heatmap.png",
  content = function(file) {
    p <- dq_missing_heatmap_plot_obj()
    if (!is.null(p)) ggsave(file, plot = p, width = 12, height = 8, dpi = 150)
  }
)

output$download_missing_matrix <- downloadHandler(
  filename = function() { paste0("missing_matrix_", Sys.Date(), ".csv") },
  content = function(file) {
    req(dq_expr_matrix())
    mat <- dq_expr_matrix()
    sel_s <- selected_samples()
    common_s <- intersect(sel_s, colnames(mat))
    if (length(common_s) == 0) { showNotification("No samples selected", type = "error"); return() }
    mat_sub <- mat[, common_s, drop = FALSE]
    miss_df <- as.data.frame(is.na(mat_sub))
    miss_df[] <- lapply(miss_df, function(x) ifelse(x, "Missing", "Detected"))
    rownames(miss_df) <- rownames(mat)
    write.csv(miss_df, file, row.names = TRUE)
  }
)

output$download_sample_missing_stats <- downloadHandler(
  filename = function() { paste0("sample_missing_stats_", Sys.Date(), ".csv") },
  content = function(file) {
    req(dq_expr_matrix())
    mat <- dq_expr_matrix()
    sel_s <- selected_samples()
    common_s <- intersect(sel_s, colnames(mat))
    if (length(common_s) == 0) return()
    mat_sub <- mat[, common_s, drop = FALSE]
    total_prot <- nrow(mat_sub)
    missing_per_sample <- colSums(is.na(mat_sub))
    valid_per_sample <- total_prot - missing_per_sample
    stats <- data.frame(
      SampleID = common_s,
      TotalProteins = total_prot,
      ValidProteins = valid_per_sample,
      MissingProteins = missing_per_sample,
      ValidPercentage = round(valid_per_sample/total_prot*100,2),
      MissingPercentage = round(missing_per_sample/total_prot*100,2)
    )
    write.csv(stats, file, row.names = FALSE)
  }
)

observeEvent(input$help_missing_heatmap, {
  showModal(modalDialog(
    title = "Missing Value Heatmap",
    "蓝色 = Detected，红色 = Missing。样本按分组排序，可勾选“Cluster rows”进行行聚类。",
    easyClose = TRUE, footer = modalButton("关闭")
  ))
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
  n_sel <- length(common_s)
  df <- data.frame(MissingRate = rates)
  ggplot(df, aes(x = MissingRate)) +
    geom_histogram(fill = "#3498db", bins = 30, alpha = 0.8, boundary = 0) +
    labs(title = paste0("Protein Missing Rate Distribution (", n_sel, " samples)"),
         x = "Missing Rate", y = "Number of Proteins") +
    theme_bw()
})

output$dq_protein_missing_hist <- renderPlot({ req(dq_protein_missing_hist()); dq_protein_missing_hist() })
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

output$dq_sample_missing_bar <- renderPlot({ req(dq_sample_missing_bar()); dq_sample_missing_bar() }, height = 400)
output$download_sample_missing_bar <- downloadHandler(
  filename = function() "sample_missing_rate.png",
  content = function(file) {
    p <- dq_sample_missing_bar()
    if (!is.null(p)) ggsave(file, plot = p, width = 12, height = 6, dpi = 150)
  }
)

# ==================== 有效值柱状图 ====================
dq_valid_plot <- reactive({
  req(dq_expr_matrix(), dq_missing_stats())
  stats <- dq_missing_stats()
  valid_percent <- (1 - stats$sample_missing) * 100
  df <- data.frame(Sample = factor(names(valid_percent), levels = names(valid_percent)), ValidPercent = valid_percent)
  ggplot(df, aes(x = Sample, y = ValidPercent)) +
    geom_col(fill = "#3498db") +
    labs(title = "Valid Values per Sample", y = "Valid Values (%)") +
    ylim(0, max(valid_percent) * 1.1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
})

output$dq_valid_values_plot <- renderPlot({ dq_valid_plot() })
output$download_valid_bar <- downloadHandler(
  filename = function() "valid_values.png",
  content = function(file) ggsave(file, plot = dq_valid_plot(), width = 8, height = 5, dpi = 150)
)
observeEvent(input$help_valid_bar, {
  showModal(modalDialog(title = "Valid Values per Sample", "每个样本中有效（非缺失）蛋白的百分比。", easyClose = TRUE, footer = modalButton("关闭")))
})

# ==================== 缺失值相关性热图（保留子热图功能，移除统计导出） ====================
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
  cor_mat[is.na(cor_mat)] <- 0  # 缺失值相关性热图保留 0 替代 NA（用户已接受）
  
  n_samp <- ncol(cor_mat)
  
  min_cor <- min(cor_mat, na.rm = TRUE)
  max_cor <- max(cor_mat, na.rm = TRUE)
  if (abs(max_cor - min_cor) < 1e-6) { min_cor <- min_cor - 0.01; max_cor <- max_cor + 0.01 }
  
  n_ticks <- 7
  legend_breaks <- pretty(c(min_cor, max_cor), n = n_ticks)
  legend_breaks <- legend_breaks[legend_breaks >= min_cor & legend_breaks <= max_cor]
  if (length(legend_breaks) < 2) legend_breaks <- seq(min_cor, max_cor, length.out = 5)
  legend_labels <- sprintf("%.2f", legend_breaks)
  message("[DEBUG] Missing cor range: [", round(min_cor,4), ", ", round(max_cor,4), "] Ticks: ", paste(legend_breaks, collapse = ", "))
  
  dist_mat <- as.dist(1 - cor_mat)
  
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
  
  main_title <- "缺失值相关性热图（Pearson 相关系数，距离 = 1 - 相关性）"
  subtitle <- if (!is.null(subset_samples)) paste0(" 样本子集（", n_samp, " 个样本）") else paste0(" 全部样本（", n_samp, " 个样本）")
  main_title <- paste0(main_title, subtitle)
  
  pheatmap::pheatmap(cor_mat,
                     main = main_title,
                     color = colorRampPalette(c("blue", "white", "red"))(100),
                     breaks = seq(min_cor, max_cor, length.out = 101),
                     legend_breaks = legend_breaks,
                     legend_labels = legend_labels,
                     clustering_distance_rows = dist_mat,
                     clustering_distance_cols = dist_mat,
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
  showModal(modalDialog(title = "Protein Intensity Distribution", "箱线图展示每个样本中蛋白强度的分布（log2 转换）。", easyClose = TRUE, footer = modalButton("关闭")))
})

# ==================== 样本相关性热图（修复图例刻度，添加详细调试） ====================
dq_cor_heatmap_plot_obj <- reactive({
  req(dq_expr_matrix())
  message("[DEBUG] dq_cor_heatmap_plot_obj: computing sample correlation")
  cor_mat <- calculate_sample_correlation(dq_expr_matrix())
  if (is.null(cor_mat)) {
    message("[DEBUG] dq_cor_heatmap_plot_obj: correlation matrix is NULL")
    return(NULL)
  }
  
  message("[DEBUG] dq_cor_heatmap_plot_obj: cor matrix dim = ", nrow(cor_mat), "x", ncol(cor_mat))
  message("[DEBUG] dq_cor_heatmap_plot_obj: NA count in cor_mat = ", sum(is.na(cor_mat)))
  
  ann_col <- NULL
  ann_colors <- NULL
  subcluster_info <- NULL
  
  if (!is.null(rv$sample_info) && "Group" %in% colnames(rv$sample_info)) {
    sample_info_short <- rv$sample_info
    rownames(sample_info_short) <- extract_sample_names(rownames(sample_info_short))
    common_samples <- intersect(colnames(cor_mat), rownames(sample_info_short))
    message("[DEBUG] dq_cor_heatmap_plot_obj: common samples with group info = ", length(common_samples))
    if (length(common_samples) > 0) {
      group_vec <- sample_info_short[common_samples, "Group"]
      ann_col <- data.frame(Group = group_vec, row.names = common_samples)
      
      if (ncol(cor_mat) >= 3) {
        cor_sub <- cor_mat[common_samples, common_samples, drop = FALSE]
        # 为聚类创建无 NA 的副本
        cor_sub_clean <- cor_sub
        cor_sub_clean[is.na(cor_sub_clean)] <- 0
        cor_dist <- as.dist(1 - cor_sub_clean)
        hc_col <- hclust(cor_dist, method = "ward.D2")
        k <- 2
        if (length(common_samples) >= 6) k <- 3
        subcluster <- cutree(hc_col, k = k)
        while (any(table(subcluster) < 2) && k > 1) { k <- k - 1; subcluster <- cutree(hc_col, k = k) }
        subcluster_label <- paste0("SubClust", subcluster)
        ann_col$Subcluster <- factor(subcluster_label, levels = unique(subcluster_label))
        subcluster_info <- subcluster
        message("[DEBUG] dq_cor_heatmap_plot_obj: subclusters created, k = ", k)
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
  
  # 计算实际数据的 min/max，忽略 NA
  min_cor <- min(cor_mat, na.rm = TRUE)
  max_cor <- max(cor_mat, na.rm = TRUE)
  message("[DEBUG] dq_cor_heatmap_plot_obj: cor range (actual): [", min_cor, ", ", max_cor, "]")
  
  # 若范围过小，稍微扩展防止 breaks 出错
  if (abs(max_cor - min_cor) < 1e-6) { min_cor <- min_cor - 0.01; max_cor <- max_cor + 0.01 }
  
  # 生成初步的图例刻度（pretty 可能超出数据范围，我们将裁剪）
  n_ticks <- 7
  raw_breaks <- pretty(c(min_cor, max_cor), n = n_ticks)
  message("[DEBUG] dq_cor_heatmap_plot_obj: raw_breaks from pretty: ", paste(raw_breaks, collapse = ", "))
  
  # 裁剪到数据实际范围
  legend_breaks <- raw_breaks[raw_breaks >= min_cor & raw_breaks <= max_cor]
  # 如果裁剪后太少，则强制使用实际 min 和 max 构成的线性序列
  if (length(legend_breaks) < 2) {
    legend_breaks <- seq(min_cor, max_cor, length.out = 5)
  }
  legend_labels <- sprintf("%.2f", legend_breaks)
  message("[DEBUG] dq_cor_heatmap_plot_obj: final legend_breaks: ", paste(legend_breaks, collapse = ", "))
  
  # 颜色映射 breaks 基于实际范围
  color_breaks <- seq(min_cor, max_cor, length.out = 101)
  
  # 聚类距离矩阵（基于 1 - cor，NA 暂时用 0 替换以避免 dist 报错）
  cor_dist_mat <- cor_mat
  cor_dist_mat[is.na(cor_dist_mat)] <- 0
  cor_dist <- as.dist(1 - cor_dist_mat)
  message("[DEBUG] dq_cor_heatmap_plot_obj: clustering distance based on 1 - cor (NAs temporarily replaced by 0)")
  
  main_title <- "样本相关性热图（基于 top 500 高变异蛋白的 log2 强度 Pearson 相关）"
  
  pheatmap::pheatmap(cor_mat,
                     main = main_title,
                     color = colorRampPalette(c("blue", "white", "red"))(100),
                     breaks = color_breaks,
                     legend_breaks = legend_breaks,
                     legend_labels = legend_labels,
                     clustering_distance_rows = cor_dist,
                     clustering_distance_cols = cor_dist,
                     clustering_method = "ward.D2",
                     show_rownames = TRUE, show_colnames = TRUE,
                     fontsize_row = 9, fontsize_col = 9,
                     angle_col = 45,
                     annotation_col = ann_col, annotation_colors = ann_colors,
                     na_col = "grey",   # NA 显示为灰色
                     silent = TRUE)
})

output$dq_cor_heatmap <- renderPlot({
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
    
    group_stats <- NULL
    subcluster_stats <- NULL
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
          if (length(idx) > 1) within_avg[[g]] <- mean(cor_sub[idx, idx][lower.tri(cor_sub[idx, idx])], na.rm = TRUE)
          else within_avg[[g]] <- NA
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
          比较 = c(names(within_avg), names(between_avg)),
          类型 = c(rep("组内", length(within_avg)), rep("组间", length(between_avg))),
          平均相关系数 = round(unlist(c(within_avg, between_avg)), 4),
          stringsAsFactors = FALSE
        )
        
        if (ncol(cor_sub) >= 3) {
          cor_sub_clean <- cor_sub
          cor_sub_clean[is.na(cor_sub_clean)] <- 0
          cor_dist <- as.dist(1 - cor_sub_clean)
          hc_col <- hclust(cor_dist, method = "ward.D2")
          k <- 2
          if (length(common_samples) >= 6) k <- 3
          subcluster <- cutree(hc_col, k = k)
          while (any(table(subcluster) < 2) && k > 1) { k <- k - 1; subcluster <- cutree(hc_col, k = k) }
          subclusters <- unique(subcluster)
          sub_within <- list()
          sub_between <- list()
          for (sc in subclusters) {
            idx <- which(subcluster == sc)
            if (length(idx) > 1) sub_within[[paste0("亚簇", sc)]] <- mean(cor_sub[idx, idx][lower.tri(cor_sub[idx, idx])], na.rm = TRUE)
            else sub_within[[paste0("亚簇", sc)]] <- NA
          }
          if (length(subclusters) >= 2) {
            for (i in 1:(length(subclusters)-1)) {
              for (j in (i+1):length(subclusters)) {
                sc1 <- subclusters[i]; sc2 <- subclusters[j]
                idx1 <- which(subcluster == sc1); idx2 <- which(subcluster == sc2)
                if (length(idx1) > 0 && length(idx2) > 0) {
                  sub_between[[paste0("亚簇", sc1, " vs 亚簇", sc2)]] <- mean(cor_sub[idx1, idx2], na.rm = TRUE)
                }
              }
            }
          }
          subcluster_stats <- data.frame(
            比较 = c(names(sub_within), names(sub_between)),
            类型 = c(rep("亚簇内", length(sub_within)), rep("亚簇间", length(sub_between))),
            平均相关系数 = round(unlist(c(sub_within, sub_between)), 4),
            stringsAsFactors = FALSE
          )
        }
      }
    }
    
    if (!is.null(group_stats)) {
      openxlsx::addWorksheet(wb, "分组相关统计")
      openxlsx::writeData(wb, "分组相关统计", group_stats)
    }
    if (!is.null(subcluster_stats)) {
      openxlsx::addWorksheet(wb, "亚簇相关统计")
      openxlsx::writeData(wb, "亚簇相关统计", subcluster_stats)
    }
    
    readme_text <- c(
      "工作簿说明：",
      "1. 相关性矩阵：基于 top 500 高变异蛋白的 log2 强度计算 Pearson 相关系数，范围为 -1 到 1。",
      "2. 分组相关统计：根据样本信息中 Group 列划分的组内和组间平均相关系数。",
      "3. 亚簇相关统计：通过列聚类树自动识别的亚簇内和亚簇间平均相关系数，聚类距离为 1 - 相关系数。"
    )
    openxlsx::addWorksheet(wb, "说明")
    openxlsx::writeData(wb, "说明", data.frame(说明 = readme_text, stringsAsFactors = FALSE))
    
    openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
  }
)

observeEvent(input$help_cor_heatmap, {
  showModal(modalDialog(title = "样本相关性热图", "基于高变异蛋白计算的样本间 Pearson 相关性。热图顶部注释条包括分组和自动识别的亚簇。聚类基于 1 - 相关系数距离。", easyClose = TRUE, footer = modalButton("关闭")))
})

# ==================== 修复后的 PCA（缺失值用 KNN 填充，移除错误过滤） ====================
dq_pca_full <- reactive({
  req(dq_expr_matrix())
  message("[DEBUG] dq_pca_full: starting PCA computation")
  expr <- dq_expr_matrix()
  
  filled <- tryCatch({
    if (requireNamespace("impute", quietly = TRUE)) {
      message("[DEBUG] dq_pca_full: impute package available, running KNN imputation")
      impute_missing_values(expr, method = "knn")
    } else {
      message("[DEBUG] dq_pca_full: impute not installed, using minimal value fill")
      expr[is.na(expr)] <- 1e-4
      expr
    }
  }, error = function(e) {
    message("[DEBUG] dq_pca_full: KNN imputation failed: ", e$message, "; fallback to min fill")
    expr[is.na(expr)] <- 1e-4
    expr
  })
  
  log_expr <- log2(filled + 1)
  message("[DEBUG] dq_pca_full: log2 transformed, checking remaining NAs: ", sum(is.na(log_expr)))
  
  row_vars <- apply(log_expr, 1, var)
  log_expr <- log_expr[row_vars > 1e-12, , drop = FALSE]
  row_unique <- apply(log_expr, 1, function(x) length(unique(x)))
  log_expr <- log_expr[row_unique > 1, , drop = FALSE]
  
  if (nrow(log_expr) < 2) {
    message("[DEBUG] dq_pca_full: insufficient variable rows for PCA")
    return(NULL)
  }
  
  pca <- tryCatch(prcomp(t(log_expr), scale. = TRUE), error = function(e) {
    message("[DEBUG] dq_pca_full: prcomp error - ", e$message)
    NULL
  })
  if (is.null(pca)) return(NULL)
  
  variance <- pca$sdev^2 / sum(pca$sdev^2) * 100
  message("[DEBUG] dq_pca_full: PCA computed, PC1 var = ", variance[1], "%, PC2 = ", variance[2], "%")
  
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
  
  list(pca = pca, scores = scores, variance = variance, loadings = pca$rotation)
})

pca_group_plot_obj <- reactive({
  pca_full <- dq_pca_full()
  if (is.null(pca_full)) {
    message("[DEBUG] pca_group_plot_obj: no PCA data")
    return(NULL)
  }
  scores <- pca_full$scores
  scores <- scores[!is.na(scores$Group), ]
  message("[DEBUG] pca_group_plot_obj: number of samples with valid group: ", nrow(scores))
  
  group_colors <- c("Control" = "#FF69B4", "Treatment" = "#00CED1")
  all_groups <- unique(scores$Group)
  missing_colors <- setdiff(all_groups, names(group_colors))
  if (length(missing_colors) > 0) {
    extra_colors <- rainbow(length(missing_colors))
    names(extra_colors) <- missing_colors
    group_colors <- c(group_colors, extra_colors)
  }
  ggplot(scores, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = group_colors) +
    labs(title = "PCA by Group",
         x = paste0("PC1 (", round(pca_full$variance[1], 1), "%)"),
         y = paste0("PC2 (", round(pca_full$variance[2], 1), "%)")) +
    theme_bw() + theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
})

output$dq_pca_group_plot <- renderPlot({ req(pca_group_plot_obj()); pca_group_plot_obj() })
output$download_pca_group <- downloadHandler(
  filename = function() "pca_group.png",
  content = function(file) ggsave(file, plot = pca_group_plot_obj(), width = 8, height = 6, dpi = 150)
)
output$download_pca_group_data <- downloadHandler(
  filename = function() { paste0("PCA_Group_Data_", Sys.Date(), ".xlsx") },
  content = function(file) {
    pca_full <- dq_pca_full()
    if (is.null(pca_full)) { showNotification("PCA not available", type = "error"); return() }
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
    readme <- data.frame(Description = c("Sheet 'pca_scores': sample coordinates; Outlier column indicates potential outlier."))
    openxlsx::addWorksheet(wb, "README")
    openxlsx::writeData(wb, "README", readme)
    openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
  }
)
observeEvent(input$help_pca_group, {
  showModal(modalDialog(title = "PCA by Group", "按实验分组着色，用于观察组间分离程度。", easyClose = TRUE, footer = modalButton("关闭")))
})

pca_batch_plot_obj <- reactive({
  pca_full <- dq_pca_full()
  if (is.null(pca_full)) return(NULL)
  scores <- pca_full$scores
  if (all(is.na(scores$Batch))) return(NULL)
  scores <- scores[!is.na(scores$Batch), ]
  batch_colors <- c("Batch1" = "#E41A1C", "Batch2" = "#00CED1")
  all_batches <- unique(scores$Batch)
  if (length(all_batches) > 2) batch_colors <- setNames(rainbow(length(all_batches)), all_batches)
  ggplot(scores, aes(x = PC1, y = PC2, color = Batch)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = batch_colors) +
    labs(title = "PCA by Batch",
         x = paste0("PC1 (", round(pca_full$variance[1], 1), "%)"),
         y = paste0("PC2 (", round(pca_full$variance[2], 1), "%)")) +
    theme_bw() + theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
})

output$dq_pca_batch_plot <- renderPlot({
  if (is.null(pca_batch_plot_obj())) { plot.new(); text(0.5, 0.5, "Batch information not available") }
  else { pca_batch_plot_obj() }
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
    pca_full <- dq_pca_full()
    if (is.null(pca_full)) { showNotification("PCA not available", type = "error"); return() }
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
    readme <- data.frame(Description = c("Sheet 'pca_scores': sample coordinates; Outlier column indicates potential outlier."))
    openxlsx::addWorksheet(wb, "README")
    openxlsx::writeData(wb, "README", readme)
    openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
  }
)
observeEvent(input$help_pca_batch, {
  showModal(modalDialog(title = "PCA by Batch", "按实验批次着色，用于检测批次效应。", easyClose = TRUE, footer = modalButton("关闭")))
})