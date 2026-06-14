# server/data_quality_plots.R
message("[DEBUG] data_quality_plots.R loaded - init empty selection, unified Arial font, larger sizes")

`%then%` <- function(a, b) { if (a) b else TRUE }
validate_condition <- function(condition, message) {
  if (!condition) { validate(message) }
}

# ==================== 样本选择 ====================
observeEvent(rv$sample_names, {
  message("[DEBUG] Updating heatmap_sample_select choices, n=", length(rv$sample_names))
  # 关键修改：初始不选择任何样本，避免大型计算
  updateSelectizeInput(session, "heatmap_sample_select", 
                       choices = rv$sample_names, 
                       selected = character(0),  # 原为 rv$sample_names
                       server = TRUE)
}, ignoreNULL = TRUE, once = FALSE)

output$heatmap_group_buttons_ui <- renderUI({
  req(rv$sample_info)
  if (!"Group" %in% colnames(rv$sample_info)) return(NULL)
  groups <- unique(rv$sample_info$Group)
  lapply(groups, function(g) {
    actionButton(inputId = paste0("heatmap_group_", g), label = g, class = "btn-sm btn-outline-info", icon = icon("filter"))
  })
})

observeEvent(input$heatmap_select_all, {
  updateSelectizeInput(session, "heatmap_sample_select", selected = rv$sample_names)
  message("[DEBUG] heatmap_select_all")
})
observeEvent(input$heatmap_clear_all, {
  updateSelectizeInput(session, "heatmap_sample_select", selected = character(0))
  message("[DEBUG] heatmap_clear_all")
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
        new_sel <- union(input$heatmap_sample_select, samples_in_group)
        updateSelectizeInput(session, "heatmap_sample_select", selected = new_sel)
        message("[DEBUG] heatmap_group button added group ", group)
      }, ignoreInit = TRUE)
    })
  }
})

# ---------- 自然排序 ----------
natural_order <- function(x) {
  if (length(x) == 0) return(integer(0))
  parts <- strsplit(as.character(x), "_")
  max_len <- max(sapply(parts, length))
  num_mat <- t(sapply(parts, function(p) {
    nums <- suppressWarnings(as.numeric(p))
    nums[is.na(nums)] <- 0
    if (length(nums) < max_len) c(nums, rep(0, max_len - length(nums)))
    else nums
  }))
  num_df <- as.data.frame(num_mat)
  do.call(order, num_df)
}

ordered_samples <- reactive({
  sel <- input$heatmap_sample_select
  if (is.null(sel) || length(sel) == 0) sel <- rv$sample_names   # 展示所有样本，但热图选择为空时显示全部
  message("[DEBUG] ordered_samples: raw selection length = ", length(sel))
  if (length(sel) == 0) return(character(0))
  
  grp <- NULL
  if (!is.null(rv$sample_info) && "Group" %in% colnames(rv$sample_info)) {
    si <- rv$sample_info
    info_short <- extract_sample_names(rownames(si))
    info_std <- standardize_sample_name(info_short)
    sel_std <- standardize_sample_name(sel)
    idx <- match(sel_std, info_std)
    grp <- ifelse(is.na(idx), "Unassigned", si$Group[idx])
  } else {
    grp <- rep("All", length(sel))
  }
  
  is_ctrl <- grepl("control|wt|ck", grp, ignore.case = TRUE)
  ctrl_groups <- sort(unique(grp[is_ctrl]))
  other_groups <- sort(setdiff(unique(grp), ctrl_groups))
  group_order <- c(ctrl_groups, other_groups)
  
  ord <- integer(0)
  for (g in group_order) {
    idx_g <- which(grp == g)
    if (length(idx_g) > 0) {
      local_order <- natural_order(sel[idx_g])
      ord <- c(ord, idx_g[local_order])
    }
  }
  sel[ord]
})

selected_samples <- reactive({ ordered_samples() })

# ==================== 缺失值热图 ====================
dq_missing_heatmap_plot_obj <- reactive({
  req(dq_expr_matrix())
  mat <- dq_expr_matrix()
  common_s <- selected_samples()
  common_s <- intersect(common_s, colnames(mat))
  if (length(common_s) == 0) return(NULL)
  
  annot_df <- NULL
  if (!is.null(rv$sample_info) && "Group" %in% colnames(rv$sample_info)) {
    si <- rv$sample_info
    info_short <- extract_sample_names(rownames(si))
    info_std <- standardize_sample_name(info_short)
    common_std <- standardize_sample_name(common_s)
    idx <- match(common_std, info_std)
    groups <- ifelse(is.na(idx), "Unknown", si$Group[idx])
    annot_df <- data.frame(Group = factor(groups), row.names = common_s)
  }
  
  mat_sub <- mat[, common_s, drop = FALSE]
  missing_mat <- (is.na(mat_sub) * 1)
  n_prot <- nrow(missing_mat)
  n_samp <- ncol(missing_mat)
  
  cluster_rows <- isTRUE(input$heatmap_cluster_rows)
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
      values = c("0" = "#3498db", "1" = "white"),
      labels = c("0" = "Detected", "1" = "Missing"),
      name = "Status"
    ) +
    labs(title = paste0("Missing Heatmap (", n_prot, " proteins)"), x = NULL, y = NULL) +
    theme_bw(base_family = "Arial") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )
  message("[DEBUG] dq_missing_heatmap_plot_obj: ", n_prot, " proteins, ", n_samp, " samples (Arial large)")
  return(p)
})

output$dq_missing_heatmap <- renderPlot({
  validate_condition(!is.null(dq_expr_matrix()), "Please upload expression data first.")
  p <- dq_missing_heatmap_plot_obj()
  if (is.null(p)) { plot.new(); text(0.5,0.5,"Please select at least one sample") }
  else print(p)
})

output$download_missing_heatmap <- downloadHandler(
  filename = function() "missing_heatmap.png",
  content = function(file) {
    p <- dq_missing_heatmap_plot_obj()
    if (!is.null(p)) ggsave(file, plot = p, width = 10, height = 8, dpi = 150)
    else showNotification("No data to download", type = "error")
    message("[DEBUG] download_missing_heatmap")
  }
)
output$download_missing_matrix <- downloadHandler(
  filename = function() paste0("missing_matrix_", Sys.Date(), ".csv"),
  content = function(file) {
    mat <- dq_expr_matrix()
    if (is.null(mat)) { showNotification("No data", type="error"); return() }
    common_s <- selected_samples()
    common_s <- intersect(common_s, colnames(mat))
    if (length(common_s)==0) { showNotification("No samples", type="error"); return() }
    miss_mat <- is.na(mat[, common_s, drop=FALSE]) * 1
    write.csv(as.data.frame(miss_mat), file, row.names = TRUE)
    message("[DEBUG] download_missing_matrix")
  }
)
observeEvent(input$help_missing_heatmap, {
  showModal(modalDialog(title="Missing Value Heatmap","Blue=Detected, White=Missing.",easyClose=TRUE,footer=modalButton("Close")))
  message("[DEBUG] help_missing_heatmap")
})

# ==================== 非缺失值直方图（已更新字体） ====================
output$sample_nonmiss_hist <- renderPlot({
  req(dq_expr_matrix())
  mat <- dq_expr_matrix()
  common_s <- selected_samples()
  common_s <- intersect(common_s, colnames(mat))
  if (length(common_s) == 0) return(NULL)
  
  nonmiss <- colSums(!is.na(mat[, common_s, drop = FALSE]))
  
  extract_prefix <- function(sample_name) {
    pref <- sub("_[0-9]+(_[0-9]+)*$", "", sample_name)
    if (nchar(pref) == 0) sample_name else pref
  }
  prefixes <- sapply(common_s, extract_prefix)
  unique_prefixes <- unique(prefixes)
  
  color_map <- c(
    "WT" = "#FFB3B3",
    "100" = "#B3FFB3",
    "200" = "#B3D9FF"
  )
  extra_colors <- setNames(
    RColorBrewer::brewer.pal(min(length(unique_prefixes), 8), "Pastel1")[1:length(unique_prefixes)],
    unique_prefixes
  )
  fill_colors <- color_map[unique_prefixes]
  missing_pref <- unique_prefixes[!unique_prefixes %in% names(color_map)]
  if (length(missing_pref) > 0) {
    fill_colors[missing_pref] <- extra_colors[missing_pref]
  }
  fill_colors <- fill_colors[!is.na(fill_colors)]
  
  message("[DEBUG] sample_nonmiss_hist: prefixes found: ", paste(unique_prefixes, collapse = ", "),
          " colors: ", paste(fill_colors, collapse = ", "))
  
  df <- data.frame(Sample = factor(common_s, levels = common_s),
                   NonMissing = nonmiss,
                   Prefix = factor(prefixes, levels = unique_prefixes))
  
  ggplot(df, aes(x = Sample, y = NonMissing, fill = Prefix)) +
    geom_col() +
    geom_text(aes(label = NonMissing), vjust = -0.5, size = 4, color = "#333333", family = "Arial") +
    scale_fill_manual(values = fill_colors, name = "Treatment Group") +
    labs(title = "Number of Quantified Proteins per Sample",
         y = "Non-Missing Count", x = NULL) +
    theme_bw(base_family = "Arial") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.title = element_text(size = 16),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      axis.text.y = element_text(size = 14),
      legend.position = "right",
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      panel.grid = element_blank()
    )
})

output$download_sample_nonmiss_hist <- downloadHandler(
  filename = function() "valid_values_per_sample.png",
  content = function(file) {
    req(dq_expr_matrix())
    mat <- dq_expr_matrix()
    common_s <- selected_samples()
    common_s <- intersect(common_s, colnames(mat))
    if (length(common_s) == 0) {
      showNotification("No samples selected", type = "error")
      return()
    }
    nonmiss <- colSums(!is.na(mat[, common_s, drop = FALSE]))
    extract_prefix <- function(sample_name) {
      pref <- sub("_[0-9]+(_[0-9]+)*$", "", sample_name)
      if (nchar(pref) == 0) sample_name else pref
    }
    prefixes <- sapply(common_s, extract_prefix)
    unique_prefixes <- unique(prefixes)
    color_map <- c("WT"="#FFB3B3","100"="#B3FFB3","200"="#B3D9FF")
    extra_colors <- setNames(RColorBrewer::brewer.pal(min(length(unique_prefixes),8),"Pastel1")[1:length(unique_prefixes)], unique_prefixes)
    fill_colors <- color_map[unique_prefixes]
    missing_pref <- unique_prefixes[!unique_prefixes %in% names(color_map)]
    if (length(missing_pref)>0) fill_colors[missing_pref] <- extra_colors[missing_pref]
    fill_colors <- fill_colors[!is.na(fill_colors)]
    df <- data.frame(Sample = factor(common_s, levels=common_s), NonMissing=nonmiss, Prefix=factor(prefixes, levels=unique_prefixes))
    p <- ggplot(df, aes(x=Sample, y=NonMissing, fill=Prefix)) +
      geom_col() +
      geom_text(aes(label=NonMissing), vjust=-0.5, size=4, color="#333333", family="Arial") +
      scale_fill_manual(values=fill_colors, name="Treatment Group") +
      labs(title="Number of Quantified Proteins per Sample", y="Non-Missing Count", x=NULL) +
      theme_bw(base_family="Arial") +
      theme(
        plot.title = element_text(hjust=0.5, face="bold", size=18),
        axis.title = element_text(size=16),
        axis.text.x = element_text(angle=45, hjust=1, size=14),
        axis.text.y = element_text(size=14),
        legend.position = "right",
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        panel.grid = element_blank()
      )
    ggsave(file, plot=p, width=12, height=6, dpi=150)
    message("[DEBUG] download_sample_nonmiss_hist: saved with Arial large")
  }
)

# ==================== 原始数据 PCA ====================
dq_pca_data <- reactive({
  tryCatch({
    raw <- expression_data()
    message("[DEBUG] dq_pca_data: raw data dim = ", nrow(raw), " x ", ncol(raw))
    if (nrow(raw) < 2 || ncol(raw) < 2) return(NULL)
    
    message("[DEBUG] dq_pca_data: imputing missing values with 1% quantile method")
    filled <- suppressMessages(impute_missing_values(raw, method = "quantile", quantile_prob = 0.01))
    log_expr <- log2(as.matrix(filled) + 1)
    
    row_vars <- apply(log_expr, 1, var)
    log_expr <- log_expr[row_vars > 1e-6, , drop = FALSE]
    if (nrow(log_expr) < 2) {
      message("[DEBUG] dq_pca_data: too few variable proteins")
      return(NULL)
    }
    
    pca <- prcomp(t(log_expr), scale. = TRUE)
    var_explained <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 1)
    scores <- as.data.frame(pca$x[, 1:2])
    scores$Sample <- rownames(scores)
    
    if (!is.null(rv$sample_info)) {
      si <- rv$sample_info
      short_names <- extract_sample_names(scores$Sample)
      short_std <- standardize_sample_name(short_names)
      info_short <- extract_sample_names(rownames(si))
      info_std <- standardize_sample_name(info_short)
      idx <- match(short_std, info_std)
      matched_count <- sum(!is.na(idx))
      message("[DEBUG] dq_pca_data: matched ", matched_count, " out of ", nrow(scores), " PCA samples to sample info")
      
      color_col <- NULL
      if ("SubGroup" %in% colnames(si)) {
        color_col <- "SubGroup"
      } else if ("Group" %in% colnames(si)) {
        color_col <- "Group"
      }
      
      if (!is.null(color_col)) {
        raw_group <- ifelse(is.na(idx), "Unassigned", si[[color_col]][idx])
      } else {
        raw_group <- "All"
      }
    } else {
      raw_group <- "All"
    }
    
    custom_order <- function(groups) {
      uniq <- unique(groups)
      wt <- grep("^WT$", uniq, ignore.case = TRUE, value = TRUE)
      others <- setdiff(uniq, wt)
      if (length(others) > 0) {
        parsed <- strsplit(others, "-")
        valid <- sapply(parsed, function(x) length(x) == 2 && !any(is.na(suppressWarnings(as.numeric(x)))))
        if (any(valid)) {
          nums <- t(sapply(parsed[valid], as.numeric))
          ord_others <- others[valid][order(nums[,1], nums[,2])]
          invalid_others <- others[!valid]
          if (length(invalid_others) > 0) ord_others <- c(ord_others, sort(invalid_others))
        } else {
          ord_others <- sort(others)
        }
      } else {
        ord_others <- character(0)
      }
      c(wt, ord_others)
    }
    
    group_levels <- custom_order(raw_group)
    if ("Unassigned" %in% group_levels) {
      group_levels <- c(setdiff(group_levels, "Unassigned"), "Unassigned")
    }
    scores$Group <- factor(raw_group, levels = group_levels)
    message("[DEBUG] dq_pca_data: Group factor levels: ", paste(levels(scores$Group), collapse=", "))
    message("[DEBUG] dq_pca_data: PCA completed, PC1 = ", var_explained[1], "%, PC2 = ", var_explained[2], "%")
    list(scores = scores, var = var_explained)
  }, error = function(e) {
    message("[ERROR] dq_pca_data: ", e$message)
    NULL
  })
})

output$dq_pca_plot <- renderPlot({
  data <- dq_pca_data()
  if (is.null(data)) {
    plot.new()
    text(0.5, 0.5, "PCA not available (upload expression data first)")
    return()
  }
  
  scores <- data$scores
  var <- data$var
  groups <- levels(scores$Group)
  
  n_groups <- length(groups)
  if (n_groups <= 8) {
    group_colors <- setNames(
      RColorBrewer::brewer.pal(max(3, n_groups), "Dark2")[1:n_groups],
      groups
    )
  } else {
    group_colors <- setNames(
      colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(n_groups),
      groups
    )
  }
  message("[DEBUG] dq_pca_plot: using colors: ", paste(group_colors, collapse=", "))
  message("[DEBUG] dq_pca_plot: applying Arial large theme")
  
  ggplot(scores, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 4, alpha = 0.9) +
    scale_color_manual(values = group_colors, drop = FALSE) +
    labs(title = "PCA (Raw Data, 1% quantile imputation) by SubGroup",
         x = paste0("PC1 (", var[1], "%)"),
         y = paste0("PC2 (", var[2], "%)")) +
    theme_bw(base_family = "Arial") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )
})

# ==================== 样本相关性热图（增大行列标签） ====================
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
  log_expr <- log2(mat + 1)
  z_expr <- t(scale(t(log_expr)))
  z_expr[!is.finite(z_expr)] <- 0
  row_vars <- apply(z_expr, 1, var, na.rm = TRUE)
  n_keep <- min(500, nrow(z_expr))
  top_idx <- order(row_vars, decreasing = TRUE)[1:n_keep]
  z_expr <- z_expr[top_idx, , drop = FALSE]
  message("[DEBUG] dq_sample_cor_data: selected top ", n_keep, " variable proteins (Z-score)")
  cor_mat <- cor(z_expr, use = "complete.obs")
  cor_mat[is.na(cor_mat)] <- 0
  cor_vals <- cor_mat[upper.tri(cor_mat)]
  message(sprintf("[DEBUG] dq_sample_cor_data: cor matrix distribution - min=%.3f, max=%.3f, median=%.3f, mean=%.3f",
                  min(cor_vals), max(cor_vals), median(cor_vals), mean(cor_vals)))
  sample_names <- colnames(cor_mat)
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
    } else if ("Group" %in% colnames(si)) {
      groups <- si$Group[idx]
      groups[is.na(groups)] <- "Unassigned"
      ann_col$Group <- factor(groups)
      groups_uniq <- levels(ann_col$Group)
      ann_colors$Group <- get_group_colors(groups_uniq)
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
  if (min_cor == max_cor) { min_cor <- min_cor - 0.01; max_cor <- max_cor + 0.01 }
  palette <- colorRampPalette(c("blue", "white", "red"))(255)
  breaks <- seq(min_cor, max_cor, length.out = 256)
  
  message("[DEBUG] dq_sample_cor_heatmap: rendering with Arial, fontsize 14")
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
                     fontsize = 14,
                     fontsize_row = 12,
                     fontsize_col = 12,
                     angle_col = 45,
                     annotation_col = dat$ann_col,
                     annotation_colors = dat$ann_colors,
                     annotation_legend = TRUE,
                     legend = TRUE,
                     main = "Sample Correlation (Z-score per protein, log2, Top500)",
                     fontfamily = "Arial")
  message("[DEBUG] dq_sample_cor_heatmap: rendered")
})

output$download_dq_sample_cor_png <- downloadHandler(
  filename = function() "sample_correlation_heatmap.png",
  content = function(file) {
    dat <- dq_sample_cor_data()
    if (!is.null(dat)) {
      png(file, width = 1600, height = 1400, res = 150, family = "Arial")
      cor_vals <- dat$cor_mat[upper.tri(dat$cor_mat)]
      min_cor <- min(cor_vals, na.rm = TRUE)
      max_cor <- max(cor_vals, na.rm = TRUE)
      if (min_cor == max_cor) { min_cor <- min_cor - 0.01; max_cor <- max_cor + 0.01 }
      palette <- colorRampPalette(c("blue", "white", "red"))(255)
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
                         fontsize = 14,
                         fontsize_row = 12,
                         fontsize_col = 12,
                         angle_col = 45,
                         annotation_col = dat$ann_col,
                         annotation_colors = dat$ann_colors,
                         annotation_legend = TRUE,
                         legend = TRUE,
                         main = "Sample Correlation (Z-score per protein, log2, Top500)",
                         fontfamily = "Arial")
      dev.off()
      message("[DEBUG] download_dq_sample_cor_png: saved with Arial large")
    }
  }
)

output$download_dq_sample_cor_matrix <- downloadHandler(
  filename = function() "sample_correlation_matrix.csv",
  content = function(file) {
    dat <- dq_sample_cor_data()
    if (!is.null(dat)) {
      write.csv(dat$cor_mat, file, row.names = TRUE)
    } else {
      showNotification("No data to download", type = "error")
    }
  }
)

message("[DEBUG] data_quality_plots.R: all outputs defined (Arial large, empty init)")