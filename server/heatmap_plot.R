# server/heatmap_plot.R
# ============================================================
# 热图数据准备与绘制模块
# ============================================================

# ---------- 辅助函数 ----------
prepare_expr_matrix <- function(data_src, samples, all_cols, protein_ids = NULL, top_n = NULL, id_map = NULL) {
  message("[DEBUG] prepare_expr_matrix: starting with samples=", length(samples), ", all_cols=", length(all_cols))
  missing_cols <- setdiff(all_cols, colnames(data_src))
  if (length(missing_cols) > 0) {
    msg <- paste("The following columns are missing from the data source:",
                 paste(missing_cols, collapse = ", "))
    message("[DEBUG] prepare_expr_matrix: ERROR - ", msg)
    return(list(error = msg))
  }
  expr_matrix <- data.matrix(data_src[, all_cols, drop = FALSE])
  original_rownames <- rownames(data_src)
  message("[DEBUG] prepare_expr_matrix: original rownames first 5 = ", paste(head(original_rownames, 5), collapse = ", "))
  
  # 如果行名全是数字，尝试用 id_map 映射（后备方案）
  if (suppressWarnings(all(!is.na(as.numeric(original_rownames))))) {
    message("[DEBUG] prepare_expr_matrix: rownames are numeric, attempting ID remapping")
    if (!is.null(id_map)) {
      idx <- as.integer(original_rownames)
      if (max(idx) <= length(id_map)) {
        new_ids <- id_map[idx]
        rownames(expr_matrix) <- new_ids
        message("[DEBUG] prepare_expr_matrix: remapped to IDs, first 5 = ", paste(head(new_ids, 5), collapse = ", "))
      } else {
        message("[DEBUG] prepare_expr_matrix: numeric indices out of range, using original rownames")
        rownames(expr_matrix) <- original_rownames
      }
    } else {
      rownames(expr_matrix) <- original_rownames
    }
  } else {
    rownames(expr_matrix) <- original_rownames
  }
  
  colnames(expr_matrix) <- samples
  expr_matrix[is.infinite(expr_matrix)] <- NA
  
  keep_rows <- rowSums(!is.na(expr_matrix)) > 0
  expr_matrix <- expr_matrix[keep_rows, , drop = FALSE]
  message("[DEBUG] prepare_expr_matrix: after removing empty rows, dim = ", nrow(expr_matrix), "x", ncol(expr_matrix))
  
  if (nrow(expr_matrix) == 0) {
    return(list(error = "No expression data available after removing NA rows."))
  }
  
  if (!is.null(protein_ids)) {
    matched <- intersect(protein_ids, rownames(expr_matrix))
    if (length(matched) == 0) {
      return(list(error = "None of the provided protein IDs were found in the data."))
    }
    expr_matrix <- expr_matrix[matched, , drop = FALSE]
    message("[DEBUG] prepare_expr_matrix: after custom protein filtering, dim = ", nrow(expr_matrix), "x", ncol(expr_matrix))
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
    message("[DEBUG] prepare_expr_matrix: removed ", na_rows_removed, " rows with NA")
  }
  
  log_expr <- log2(expr_matrix + 1)
  
  if (!is.null(top_n)) {
    row_var <- apply(log_expr, 1, var, na.rm = TRUE)
    row_var[is.na(row_var)] <- 0
    n <- min(top_n, nrow(log_expr))
    top_idx <- order(row_var, decreasing = TRUE)[1:n]
    log_expr <- log_expr[top_idx, , drop = FALSE]
    message("[DEBUG] prepare_expr_matrix: top_n selection, now ", nrow(log_expr), " rows")
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
  
  message("[DEBUG] prepare_expr_matrix: final matrix dimensions=", nrow(expr_z), "x", ncol(expr_z))
  list(mat = expr_z, has_na = has_na, na_rows_removed = na_rows_removed)
}

# ---------- 原始强度样本名提取（用于 Intensity 模式） ----------
heatmap_raw_sample_names <- reactive({
  req(rv$raw_data)
  int_cols <- grep("^Intensity ", colnames(rv$raw_data), value = TRUE)
  if (length(int_cols) == 0) {
    message("[DEBUG] heatmap_raw_sample_names: no Intensity columns found")
    return(character(0))
  }
  samples <- gsub("^Intensity ", "", int_cols)
  message("[DEBUG] heatmap_raw_sample_names: found ", length(samples), " Intensity samples, first 3: ", paste(head(samples, 3), collapse = ", "))
  samples
})

# ---------- Intensity 模式数据矩阵（已修复行名为蛋白ID） ----------
heatmap_raw_data <- reactive({
  req(rv$raw_data)
  int_cols <- grep("^Intensity ", colnames(rv$raw_data), value = TRUE)
  if (length(int_cols) == 0) {
    message("[DEBUG] heatmap_raw_data: no Intensity columns, returning NULL")
    return(NULL)
  }
  mat <- rv$raw_data
  
  # 获取蛋白ID：优先用 raw_data 的 Master protein IDs，否则用 clean_data
  if ("Master protein IDs" %in% colnames(mat)) {
    ids <- as.character(mat[["Master protein IDs"]])
    message("[DEBUG] heatmap_raw_data: got IDs from raw_data$`Master protein IDs`")
  } else if (!is.null(rv$clean_data) && "Master protein IDs" %in% colnames(rv$clean_data)) {
    ids <- as.character(rv$clean_data[["Master protein IDs"]])
    message("[DEBUG] heatmap_raw_data: got IDs from clean_data$`Master protein IDs`")
  } else {
    ids <- rownames(mat)
    message("[DEBUG] heatmap_raw_data: using rownames as IDs")
  }
  message("[DEBUG] heatmap_raw_data: first 5 IDs = ", paste(head(ids, 5), collapse = ", "))
  
  mat <- mat[, int_cols, drop = FALSE]
  mat <- suppressWarnings(as.data.frame(lapply(mat, as.numeric)))
  mat[mat == 0] <- NA
  rownames(mat) <- ids
  colnames(mat) <- gsub("^Intensity ", "", int_cols)
  message("[DEBUG] heatmap_raw_data: created matrix with dim ", nrow(mat), "x", ncol(mat))
  mat
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

# ---------- 核心热图数据准备（按钮触发） ----------
heatmap_data <- eventReactive(input$generate_heatmap, {
  message("[DEBUG] heatmap_data triggered by Generate Heatmap button (click=", input$generate_heatmap, ")")
  result <- list(error = NULL, mat = NULL, has_na = FALSE, na_rows_removed = 0)
  
  if (is.null(rv$raw_data)) {
    result$error <- "Please upload an expression matrix file first."
    return(result)
  }
  
  mode <- input$heatmap_protein_mode
  if (is.null(mode)) mode <- "top_n"
  
  # 全局蛋白ID映射向量（后备）
  global_id_map <- NULL
  if (!is.null(rv$clean_data) && "Master protein IDs" %in% colnames(rv$clean_data)) {
    global_id_map <- as.character(rv$clean_data$`Master protein IDs`)
    message("[DEBUG] heatmap_data: global_id_map length = ", length(global_id_map))
  } else {
    message("[DEBUG] heatmap_data: no global_id_map available")
  }
  
  if (input$heatmap_data_source == "LFQ") {
    data_src <- get_analysis_matrix()
    if (is.null(data_src)) {
      result$error <- "Preprocessing has not been run or data is not available for the current intensity type. Please run preprocessing first."
      return(result)
    }
    src_colnames <- colnames(data_src)
    src_short <- extract_sample_names(src_colnames)
    message("[DEBUG] heatmap_data LFQ: src_colnames first 3: ", paste(head(src_colnames, 3), collapse = ", "))
    message("[DEBUG] heatmap_data LFQ: src_short first 3: ", paste(head(src_short, 3), collapse = ", "))
    
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
    keep <- src_short %in% selected_samples
    message("[DEBUG] heatmap_data LFQ: selected_samples count=", length(selected_samples),
            ", matched in data source=", sum(keep))
    if (!any(keep)) {
      result$error <- "None of the selected group samples match the columns in the data source."
      return(result)
    }
    samples <- src_short[keep]
    all_cols <- src_colnames[keep]
    group_vec <- rep("Unassigned", length(samples))
    for (g in groups_sel) group_vec[samples %in% rv$groups[[g]]] <- g
  } else {   # "Intensity" 模式
    data_src <- heatmap_raw_data()
    if (is.null(data_src)) {
      result$error <- "No Intensity columns found in the uploaded data."
      return(result)
    }
    src_short <- colnames(data_src)
    message("[DEBUG] heatmap_data Intensity: src_short first 3: ", paste(head(src_short, 3), collapse = ", "))
    message("[DEBUG] heatmap_data Intensity: first 5 rownames = ", paste(head(rownames(data_src), 5), collapse = ", "))
    
    if (mode == "custom") {
      selected_samples <- src_short
    } else {
      groups <- heatmap_raw_groups()
      if (is.null(groups) || is.null(input$heatmap_selected_groups)) {
        selected_samples <- src_short
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
    keep <- src_short %in% selected_samples
    message("[DEBUG] heatmap_data Intensity: selected_samples count=", length(selected_samples),
            ", matched=", sum(keep))
    if (!any(keep)) {
      result$error <- "None of the selected samples match the columns."
      return(result)
    }
    samples <- src_short[keep]
    all_cols <- samples
    group_vec <- rep("All", length(samples))
    if (mode != "custom" && !is.null(groups)) {
      group_map <- setNames(rep(names(groups), lengths(groups)), unlist(groups))
      group_vec <- group_map[samples]
      group_vec[is.na(group_vec)] <- "Other"
    }
  }
  
  message("[DEBUG] heatmap_data: final samples count=", length(samples),
          ", first 3 samples: ", paste(head(samples, 3), collapse = ", "))
  message("[DEBUG] heatmap_data: final all_cols first 3: ", paste(head(all_cols, 3), collapse = ", "))
  
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
  
  res_mat <- prepare_expr_matrix(data_src, samples, all_cols, 
                                 protein_ids = protein_ids, top_n = top_n,
                                 id_map = global_id_map)
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
  
  heatmap_generated_version(data_changed_trigger())
  message("[DEBUG] heatmap_data: successfully generated, heatmap_generated_version set to ", data_changed_trigger())
  result
})

# ---------- 显示数据反应式（整合版本判断） ----------
heatmap_display_data <- reactive({
  current_trigger <- data_changed_trigger()
  dat <- heatmap_data()
  
  if (is.null(dat)) {
    return(list(error = "Click 'Generate Heatmap' to start"))
  }
  
  if (current_trigger != heatmap_generated_version()) {
    return(list(error = "Data has changed. Please click 'Generate Heatmap' to update."))
  }
  
  dat
})

make_heatmap_breaks <- function(mat, n_colors = 100) {
  rng <- range(mat, na.rm = TRUE)
  if (rng[2] - rng[1] < 1e-6) {
    rng <- c(rng[1] - 0.5, rng[2] + 0.5)
  }
  seq(rng[1], rng[2], length.out = n_colors + 1)
}

output$heatmap_plot <- renderPlot({
  dat <- heatmap_display_data()
  
  if (!is.null(dat$error)) {
    plot.new()
    text(0.5, 0.5, dat$error, cex = 1.2)
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
    dat <- heatmap_display_data()
    if (!is.null(dat$error)) {
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

# ========== 热图数据源信息 ==========
output$heatmap_data_source_info <- renderPrint({
  message("[DEBUG] output$heatmap_data_source_info called")
  src <- input$heatmap_data_source
  cat("Data Source:", if (src == "LFQ") "LFQ Intensity (per-row Z-score)" else "Intensity (per-row Z-score)", "\n")
  if (src == "LFQ") {
    cat("Derived from: Preprocessed data (after all steps: filtering, imputation, batch correction if applied)\n")
    if (is.null(processed_data())) {
      cat("Status: Preprocessing has NOT been run. Please click 'Run Preprocessing' in Data Preprocessing tab before generating heatmaps.\n")
    } else {
      cat("Preprocessing was performed at:", format(preprocessing_params$last_run_time, "%Y-%m-%d %H:%M:%S"), "\n")
    }
    message("[DEBUG] heatmap_data_source_info: LFQ mode, preprocessing done = ", !is.null(processed_data()))
  } else {
    cat("Derived from: Raw uploaded data (Intensity columns, no preprocessing applied).\n")
    cat("Note: Missing values are removed row-wise; imputation is NOT used for this mode.\n")
    cat("This data source reflects the original Intensity values without any filtering or imputation.\n")
    message("[DEBUG] heatmap_data_source_info: Intensity mode (raw data)")
  }
})

# ========== 热图预处理步骤指示器 ==========
output$heatmap_preprocess_steps <- renderUI({
  steps <- list()
  src <- input$heatmap_data_source
  if (src == "LFQ") {
    # 读取预处理参数
    steps <- c(steps, "Data source: Preprocessed data (after filtering, imputation, batch correction)")
    if (!is.null(preprocessing_params$last_run_time)) {
      steps <- c(steps, paste0("Last preprocessing: ", format(preprocessing_params$last_run_time, "%Y-%m-%d %H:%M")))
    }
    steps <- c(steps, "log2(Intensity + 1) transformation applied")
    steps <- c(steps, "Per-row Z-score normalization (scale)")
  } else {
    steps <- c(steps, "Data source: Raw Intensity columns")
    steps <- c(steps, "log2(Intensity + 1) transformation applied")
    steps <- c(steps, "Per-row Z-score normalization (scale)")
    steps <- c(steps, "Note: Missing values are removed row-wise; no imputation is performed.")
  }
  
  step_tags <- lapply(seq_along(steps), function(i) {
    tagList(
      if (i > 1) tags$span(style = "font-size: 20px; color: #e67e22; margin: 0 8px;", "→"),
      tags$span(style = "background: #e8f0fe; padding: 6px 12px; border-radius: 15px; font-size: 13px;", steps[[i]])
    )
  })
  div(
    style = "margin-bottom: 15px; padding: 10px; background: #f8f9fa; border-radius: 8px; border: 1px solid #dee2e6;",
    p(strong(icon("info-circle"), " Data processing steps for heatmap:")),
    div(style = "display: flex; flex-wrap: wrap; align-items: center;", do.call(tagList, step_tags))
  )
})

message("[DEBUG] heatmap_plot.R fully loaded")