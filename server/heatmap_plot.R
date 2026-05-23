# server/heatmap_plot.R

# ---------- 辅助：提取原始强度矩阵并进行 log2 和 Z-score ----------
prepare_expr_matrix <- function(data_src, samples, all_cols, protein_ids = NULL, top_n = NULL) {
  expr_matrix <- data.matrix(data_src[, all_cols, drop = FALSE])
  rownames(expr_matrix) <- rownames(data_src)
  colnames(expr_matrix) <- samples
  expr_matrix[is.infinite(expr_matrix)] <- NA
  
  keep_rows <- rowSums(!is.na(expr_matrix)) > 0
  expr_matrix <- expr_matrix[keep_rows, , drop = FALSE]
  
  if (nrow(expr_matrix) == 0) {
    return(list(error = "No expression data available after removing NA rows."))
  }
  
  if (!is.null(protein_ids)) {
    matched <- intersect(protein_ids, rownames(expr_matrix))
    if (length(matched) == 0) {
      return(list(error = "None of the provided protein IDs were found in the data."))
    }
    expr_matrix <- expr_matrix[matched, , drop = FALSE]
  }
  
  has_na <- any(is.na(expr_matrix))
  na_rows_removed <- 0
  if (has_na) {
    na_rows <- rowSums(is.na(expr_matrix)) > 0
    expr_matrix <- expr_matrix[!na_rows, , drop = FALSE]
    na_rows_removed <- sum(na_rows)
    if (nrow(expr_matrix) == 0) {
      return(list(error = "After removing rows with NA values, no proteins remain."))
    }
  }
  
  log_expr <- log2(expr_matrix + 1)
  
  if (!is.null(top_n)) {
    row_var <- apply(log_expr, 1, var, na.rm = TRUE)
    row_var[is.na(row_var)] <- 0
    n <- min(top_n, nrow(log_expr))
    top_idx <- order(row_var, decreasing = TRUE)[1:n]
    log_expr <- log_expr[top_idx, , drop = FALSE]
  }
  
  if (nrow(log_expr) == 0) {
    return(list(error = "No proteins left after filtering."))
  }
  
  expr_z <- t(apply(log_expr, 1, function(x) {
    s <- sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(0, length(x)))
    (x - mean(x, na.rm = TRUE)) / s
  }))
  colnames(expr_z) <- colnames(log_expr)
  expr_z[!is.finite(expr_z)] <- 0
  
  list(mat = expr_z, has_na = has_na, na_rows_removed = na_rows_removed)
}

# ---------- 原始强度样本名提取（用于 Raw Intensity 模式） ----------
heatmap_raw_sample_names <- reactive({
  req(rv$raw_data)
  int_cols <- grep("^Intensity ", colnames(rv$raw_data), value = TRUE)
  if (length(int_cols) == 0) return(character(0))
  gsub("^Intensity ", "", int_cols)
})

observe({
  req(input$heatmap_data_source == "Intensity")
  samples <- heatmap_raw_sample_names()
  if (length(samples) == 0) {
    updateSelectInput(session, "heatmap_group_level", choices = list("No Intensity columns" = ""))
    return()
  }
  levels <- parse_sample_levels(samples)
  if (length(levels) == 0) {
    updateSelectInput(session, "heatmap_group_level", choices = list("Default (prefix)" = "default"))
    return()
  }
  display_texts <- sapply(names(levels), function(level_num) {
    l <- levels[[level_num]]
    paste0("Level ", l$level, " (e.g., ", l$example, ")")
  })
  return_values <- names(levels)
  choices <- setNames(return_values, display_texts)
  choices <- c(choices, "Default (prefix)" = "default")
  updateSelectInput(session, "heatmap_group_level", choices = choices, selected = "default")
})

heatmap_raw_groups <- reactiveVal(NULL)

observeEvent(input$heatmap_apply_grouping, {
  req(input$heatmap_data_source == "Intensity")
  samples <- heatmap_raw_sample_names()
  if (length(samples) == 0) {
    showNotification("No Intensity columns found.", type = "error")
    return()
  }
  level <- input$heatmap_group_level
  if (is.null(level) || level == "") {
    showNotification("Please select a grouping level.", type = "warning")
    return()
  }
  levels_info <- parse_sample_levels(samples)
  separator <- if (level != "default" && !is.null(levels_info[[level]])) {
    levels_info[[level]]$separator
  } else {
    "-"
  }
  group_assign <- sapply(samples, function(s) extract_group_prefix(s, level, separator))
  group_assign[is.na(group_assign) | group_assign == ""] <- "Other"
  groups <- split(samples, group_assign)
  heatmap_raw_groups(groups)
  showNotification("Grouping applied!", type = "message")
})

output$heatmap_group_selection_ui <- renderUI({
  groups <- heatmap_raw_groups()
  if (is.null(groups)) return(p("Click 'Apply Grouping' to generate groups."))
  group_names <- names(groups)
  checkboxGroupInput("heatmap_selected_groups", "Select Groups to Include",
                     choices = group_names, selected = group_names, inline = TRUE)
})

# ---------- 热图数据（修正判断条件） ----------
heatmap_data <- eventReactive(input$generate_heatmap, {
  result <- list(error = NULL, mat = NULL, has_na = FALSE, na_rows_removed = 0)
  
  if (is.null(rv$raw_data)) {
    result$error <- "No data uploaded. Please upload a file first."
    return(result)
  }
  
  data_src <- get_analysis_matrix()
  if (is.null(data_src)) {
    result$error <- "No expression data available. Please upload data or run preprocessing."
    return(result)
  }
  
  all_lfq_cols <- rv$lfq_cols
  all_samples  <- rv$sample_names
  
  if (length(all_lfq_cols) == 0) {
    result$error <- "No intensity columns found. Please check your Data Upload settings."
    return(result)
  }
  
  mode <- input$heatmap_protein_mode
  if (is.null(mode)) mode <- "top_n"
  
  # 修正：LFQ 或 Intensity
  if (input$heatmap_data_source == "LFQ") {
    if (is.null(input$heatmap_groups) || length(input$heatmap_groups) == 0) {
      result$error <- "Please select at least one group."
      return(result)
    }
    groups_sel <- input$heatmap_groups
    selected_samples <- unlist(rv$groups[groups_sel])
    if (length(selected_samples) == 0) {
      result$error <- "Selected groups contain no samples."
      return(result)
    }
    keep <- all_samples %in% selected_samples
    if (!any(keep)) {
      result$error <- "None of the selected group samples match the columns in the data."
      return(result)
    }
    samples <- all_samples[keep]
    col_positions <- which(keep)
    all_cols <- colnames(data_src)[col_positions]
    group_vec <- rep("Unassigned", length(samples))
    for (g in groups_sel) group_vec[samples %in% rv$groups[[g]]] <- g
  } else {   # "Intensity"
    if (mode == "custom") {
      selected_samples <- all_samples
    } else {
      groups <- heatmap_raw_groups()
      if (is.null(groups) || is.null(input$heatmap_selected_groups)) {
        selected_samples <- all_samples
      } else {
        selected_groups <- input$heatmap_selected_groups
        if (length(selected_groups) == 0) {
          result$error <- "No groups selected."
          return(result)
        }
        selected_samples <- unlist(groups[selected_groups])
        if (length(selected_samples) == 0) {
          result$error <- "Selected groups contain no samples."
          return(result)
        }
      }
    }
    keep <- all_samples %in% selected_samples
    if (!any(keep)) {
      result$error <- "None of the selected samples match the columns."
      return(result)
    }
    samples <- all_samples[keep]
    col_positions <- which(keep)
    all_cols <- colnames(data_src)[col_positions]
    group_vec <- rep("All", length(samples))
    if (mode != "custom" && !is.null(groups)) {
      group_map <- setNames(rep(names(groups), lengths(groups)), unlist(groups))
      group_vec <- group_map[samples]
      group_vec[is.na(group_vec)] <- "Other"
    }
  }
  
  protein_ids <- NULL
  top_n <- NULL
  if (mode == "custom") {
    custom_text <- input$heatmap_custom_ids
    if (is.null(custom_text) || trimws(custom_text) == "") {
      result$error <- "Please enter at least one Master Protein ID."
      return(result)
    }
    ids <- trimws(unlist(strsplit(custom_text, "[,\n;]+")))
    ids <- ids[ids != ""]
    if (length(ids) == 0) {
      result$error <- "Invalid protein ID input."
      return(result)
    }
    protein_ids <- ids
  } else {
    top_n <- input$heatmap_top_n
    if (is.null(top_n) || is.na(top_n)) top_n <- 20
  }
  
  res_mat <- prepare_expr_matrix(data_src, samples, all_cols, protein_ids = protein_ids, top_n = top_n)
  if (!is.null(res_mat$error)) {
    result$error <- res_mat$error
    return(result)
  }
  
  result$has_na <- isTRUE(res_mat$has_na)
  result$na_rows_removed <- res_mat$na_rows_removed %||% 0
  
  ann_col <- data.frame(Group = factor(group_vec), row.names = samples)
  ann_colors <- list(Group = get_group_colors(levels(ann_col$Group)))
  
  result$mat <- res_mat$mat
  result$annotation_col <- ann_col
  result$annotation_colors <- ann_colors
  result$show_sample_names <- input$heatmap_show_sample_names
  result
})

make_heatmap_breaks <- function(mat, n_colors = 100) {
  rng <- range(mat, na.rm = TRUE)
  if (rng[2] - rng[1] < 1e-6) {
    rng <- c(rng[1] - 0.5, rng[2] + 0.5)
  }
  seq(rng[1], rng[2], length.out = n_colors + 1)
}

output$heatmap_plot <- renderPlot({
  dat <- heatmap_data()
  if (is.null(dat) || (!is.null(dat$error) && dat$error != "")) {
    plot.new()
    text(0.5, 0.5, ifelse(is.null(dat), "No data returned.", dat$error), cex = 1.2)
    return()
  }
  if (is.null(dat$mat) || nrow(dat$mat) == 0 || ncol(dat$mat) == 0) {
    plot.new()
    text(0.5, 0.5, "No data to display")
    return()
  }
  
  if (isTRUE(dat$has_na)) {
    showNotification(
      paste0("Note: ", dat$na_rows_removed, " proteins with missing values were excluded from the heatmap. Consider using imputation in preprocessing for a more complete view."),
      type = "message", duration = 10, id = "heatmap_na_note"
    )
  }
  
  brks <- make_heatmap_breaks(dat$mat)
  ann <- dat$annotation_col
  if (!is.null(ann) && nrow(ann) == 0) ann <- NULL
  pheatmap(dat$mat,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           breaks = brks,
           scale = "none",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = TRUE,
           show_colnames = dat$show_sample_names,
           annotation_col = ann,
           annotation_colors = dat$annotation_colors,
           main = "Expression Heatmap",
           fontsize_row = 8,
           fontsize_col = 8)
}, height = 700)

output$download_heatmap_png <- downloadHandler(
  filename = function() paste0("Heatmap_", Sys.Date(), ".png"),
  content = function(file) {
    dat <- heatmap_data()
    if (!is.null(dat$error) && dat$error != "") {
      png(file, width = 1200, height = 1000, res = 150)
      plot.new(); text(0.5, 0.5, dat$error, cex = 1.5)
      dev.off()
      return()
    }
    if (is.null(dat$mat) || nrow(dat$mat) == 0 || ncol(dat$mat) == 0) {
      png(file, width = 1200, height = 1000, res = 150)
      plot.new(); text(0.5, 0.5, "No data to display", cex = 1.5)
      dev.off()
      return()
    }
    png(file, width = 1200, height = 1000, res = 150)
    brks <- make_heatmap_breaks(dat$mat)
    ann <- dat$annotation_col
    if (!is.null(ann) && nrow(ann) == 0) ann <- NULL
    p <- pheatmap(dat$mat,
                  color = colorRampPalette(c("blue", "white", "red"))(100),
                  breaks = brks,
                  scale = "none",
                  cluster_rows = TRUE,
                  cluster_cols = TRUE,
                  show_rownames = TRUE,
                  show_colnames = dat$show_sample_names,
                  annotation_col = ann,
                  annotation_colors = dat$annotation_colors,
                  main = "Expression Heatmap")
    print(p)
    dev.off()
  }
)