# server/heatmap_plot.R
message("[DEBUG] heatmap_plot.R loading... (Arial font in pheatmap)")

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
  
  if (suppressWarnings(all(!is.na(as.numeric(original_rownames))))) {
    if (!is.null(id_map)) {
      idx <- as.integer(original_rownames)
      if (max(idx) <= length(id_map)) {
        new_ids <- id_map[idx]
        rownames(expr_matrix) <- new_ids
        message("[DEBUG] prepare_expr_matrix: remapped to IDs, first 5 = ", paste(head(new_ids, 5), collapse = ", "))
      } else {
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
    message("[DEBUG] prepare_expr_matrix: removed ", na_rows_removed, " rows with NA")
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
  if (length(int_cols) == 0) return(character(0))
  samples <- gsub("^Intensity ", "", int_cols)
  samples
})

# ---------- Intensity 模式数据矩阵 ----------
heatmap_raw_data <- reactive({
  req(rv$raw_data)
  int_cols <- grep("^Intensity ", colnames(rv$raw_data), value = TRUE)
  if (length(int_cols) == 0) return(NULL)
  mat <- rv$raw_data
  
  if ("Master protein IDs" %in% colnames(mat)) {
    ids <- as.character(mat[["Master protein IDs"]])
  } else if (!is.null(rv$clean_data) && "Master protein IDs" %in% colnames(rv$clean_data)) {
    ids <- as.character(rv$clean_data[["Master protein IDs"]])
  } else {
    ids <- rownames(mat)
  }
  mat <- mat[, int_cols, drop = FALSE]
  mat <- suppressWarnings(as.data.frame(lapply(mat, as.numeric)))
  mat[mat == 0] <- NA
  rownames(mat) <- ids
  colnames(mat) <- gsub("^Intensity ", "", int_cols)
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
  choices <- c(choices, "Use SubGroup from Sample Info" = "subgroup")
  updateSelectInput(session, "heatmap_group_level", choices = choices, selected = "subgroup")
})

heatmap_raw_groups <- reactiveVal(NULL)

observeEvent(input$heatmap_apply_grouping, {
  req(input$heatmap_data_source == "Intensity")
  samples <- heatmap_raw_sample_names()
  if (length(samples) == 0) {
    showNotification("No Intensity columns found.", type = "error")
    return()
  }
  
  # 优先使用样本信息表的 SubGroup
  if (!is.null(rv$sample_info) && "SubGroup" %in% colnames(rv$sample_info)) {
    message("[DEBUG] heatmap_apply_grouping: using SubGroup from sample info")
    si <- rv$sample_info
    # 提取样本信息表的短名，并标准化
    info_raw <- rownames(si)
    info_std <- standardize_sample_name(info_raw)
    
    # 标准化当前样本名
    samples_std <- standardize_sample_name(samples)
    
    # 匹配
    idx <- match(samples_std, info_std)
    group_assign <- rep("Unassigned", length(samples))
    matched <- !is.na(idx)
    if (any(matched)) {
      group_assign[matched] <- si$SubGroup[idx[matched]]
      message("[DEBUG] heatmap_apply_grouping: matched ", sum(matched), " samples to SubGroup")
    } else {
      message("[DEBUG] heatmap_apply_grouping: no samples matched to SubGroup")
    }
    groups <- split(samples, group_assign)
    heatmap_raw_groups(groups)
    showNotification("Groups created based on SubGroup.", type = "message")
    message("[DEBUG] heatmap_apply_grouping: groups = ", paste(names(groups), collapse = ", "))
    return()
  }
  
  # 回退到原有的按层级分组
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
  message("[DEBUG] heatmap_apply_grouping: groups = ", paste(names(groups), collapse = ", "))
})

output$heatmap_group_selection_ui <- renderUI({
  groups <- heatmap_raw_groups()
  if (is.null(groups)) return(p("Click 'Apply Grouping' to generate groups."))
  group_names <- names(groups)
  checkboxGroupInput("heatmap_selected_groups", "Select Groups to Include",
                     choices = group_names, selected = group_names, inline = TRUE)
})

# ---------- 核心热图数据准备 ----------
heatmap_data <- eventReactive(input$generate_heatmap, {
  message("[DEBUG] heatmap_data triggered by Generate Heatmap button")
  result <- list(error = NULL, mat = NULL, has_na = FALSE, na_rows_removed = 0)
  
  if (is.null(rv$raw_data)) {
    result$error <- "Please upload an expression matrix file first."
    return(result)
  }
  
  mode <- input$heatmap_protein_mode
  if (is.null(mode)) mode <- "top_n"
  
  global_id_map <- NULL
  if (!is.null(rv$clean_data) && "Master protein IDs" %in% colnames(rv$clean_data)) {
    global_id_map <- as.character(rv$clean_data$`Master protein IDs`)
  }
  
  data_source_type <- input$heatmap_data_source
  use_normalized <- FALSE
  if (data_source_type == "LFQ" && !is.null(input$heatmap_normalization)) {
    use_normalized <- (input$heatmap_normalization == "total")
  }
  
  if (data_source_type == "LFQ") {
    if (use_normalized) {
      nd <- norm_data_full()
      if (is.null(nd)) {
        result$error <- "Normalized data not available. Please run preprocessing and set a baseline sample."
        return(result)
      }
      norm_cols <- grep("^Norm_LFQ intensity ", colnames(nd), value = TRUE)
      if (length(norm_cols) == 0) {
        result$error <- "No normalized LFQ intensity columns found."
        return(result)
      }
      data_src <- nd
      src_colnames <- norm_cols
      src_short <- gsub("^Norm_LFQ intensity ", "", norm_cols)
    } else {
      data_src <- get_analysis_matrix()
      if (is.null(data_src)) {
        result$error <- "Preprocessing has not been run or data is not available for the current intensity type. Please run preprocessing first."
        return(result)
      }
      src_colnames <- colnames(data_src)
      src_short <- extract_sample_names(src_colnames)
    }
    
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
    if (!any(keep)) {
      result$error <- "None of the selected group samples match the columns in the data source."
      return(result)
    }
    samples <- src_short[keep]
    all_cols <- src_colnames[keep]
    group_vec <- rep("Unassigned", length(samples))
    for (g in groups_sel) group_vec[samples %in% rv$groups[[g]]] <- g
  } else {
    data_src <- heatmap_raw_data()
    if (is.null(data_src)) {
      result$error <- "No Intensity columns found in the uploaded data."
      return(result)
    }
    src_short <- colnames(data_src)
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
    if (is.null(top_n) || is.na(top_n)) top_n <- 50
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
  result$use_normalized <- use_normalized
  
  heatmap_generated_version(data_changed_trigger())
  message("[DEBUG] heatmap_data: successfully generated, heatmap_generated_version set to ", data_changed_trigger())
  result
})

heatmap_display_data <- reactive({
  dat <- heatmap_data()
  if (is.null(dat)) return(list(error = "Click 'Generate Heatmap' to start"))
  return(dat)
})

make_heatmap_breaks <- function(mat, n_colors = 100) {
  rng <- range(mat, na.rm = TRUE)
  if (rng[2] - rng[1] < 1e-6) rng <- c(rng[1] - 0.5, rng[2] + 0.5)
  seq(rng[1], rng[2], length.out = n_colors + 1)
}

output$heatmap_plot <- renderPlot({
  message("[DEBUG] heatmap_plot renderPlot called")
  dat <- heatmap_display_data()
  if (!is.null(dat$error)) {
    message("[DEBUG] heatmap_plot error: ", dat$error)
    plot.new(); text(0.5, 0.5, dat$error, cex = 1.2); return()
  }
  if (is.null(dat$mat) || nrow(dat$mat) == 0 || ncol(dat$mat) == 0) {
    message("[DEBUG] heatmap_plot: empty matrix")
    plot.new(); text(0.5, 0.5, "No data to display"); return()
  }
  if (isTRUE(dat$has_na)) {
    showNotification(
      paste0("Note: ", dat$na_rows_removed, " proteins with missing values were excluded from the heatmap."),
      type = "message", duration = 10, id = "heatmap_na_note"
    )
  }
  
  brks <- make_heatmap_breaks(dat$mat)
  ann <- dat$annotation_col
  if (!is.null(ann) && nrow(ann) == 0) ann <- NULL
  
  tryCatch({
    p <- pheatmap::pheatmap(dat$mat,
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
                            fontsize_col = 8,
                            fontfamily = "Arial",   # 设置字体
                            silent = FALSE)
    message("[DEBUG] heatmap_plot: pheatmap rendered with Arial font")
  }, error = function(e) {
    message("[ERROR] heatmap_plot: ", e$message)
    plot.new()
    text(0.5, 0.5, paste("Heatmap error:", e$message), cex = 1.2)
  })
}, height = 700)

output$download_heatmap_png <- downloadHandler(
  filename = function() paste0("Heatmap_", Sys.Date(), ".png"),
  content = function(file) {
    dat <- heatmap_display_data()
    if (!is.null(dat$error)) {
      png(file, width = 1200, height = 1000, res = 150)
      plot.new(); text(0.5, 0.5, dat$error, cex = 1.5); dev.off()
      return()
    }
    if (is.null(dat$mat) || nrow(dat$mat) == 0 || ncol(dat$mat) == 0) {
      png(file, width = 1200, height = 1000, res = 150)
      plot.new(); text(0.5, 0.5, "No data to display", cex = 1.5); dev.off()
      return()
    }
    png(file, width = 1200, height = 1000, res = 150, family = "Arial")
    brks <- make_heatmap_breaks(dat$mat)
    ann <- dat$annotation_col
    if (!is.null(ann) && nrow(ann) == 0) ann <- NULL
    tryCatch({
      pheatmap::pheatmap(dat$mat,
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
                         fontfamily = "Arial")
      dev.off()
    }, error = function(e) {
      message("[ERROR] download_heatmap_png: ", e$message)
    })
  }
)

output$heatmap_data_source_info <- renderPrint({
  message("[DEBUG] output$heatmap_data_source_info called")
  src <- input$heatmap_data_source
  norm <- input$heatmap_normalization %||% "none"
  cat("Data Source:", if (src == "LFQ") "LFQ Intensity" else "Intensity", "\n")
  if (src == "LFQ") {
    if (norm == "total") {
      cat("Normalization: Total Intensity (baseline sample) applied\n")
    } else {
      cat("Normalization: None (preprocessed raw intensity)\n")
    }
    if (is.null(processed_data())) {
      cat("Status: Preprocessing has NOT been run.\n")
    } else {
      cat("Preprocessing was performed at:", format(preprocessing_params$last_run_time, "%Y-%m-%d %H:%M:%S"), "\n")
    }
  } else {
    cat("Derived from: Raw uploaded data (Intensity columns, no preprocessing applied).\n")
    cat("Note: Missing values are removed row-wise; imputation is NOT used for this mode.\n")
  }
})

output$heatmap_preprocess_steps <- renderUI({
  steps <- list()
  steps <- c(steps, paste0("Missing Value Filter: threshold = ", input$max_missing_fraction %||% 0.5,
                           ", mode = ", preprocessing_params$missing_filter_mode %||% "global"))
  min_int <- input$min_intensity
  if (!is.null(min_int) && !is.na(min_int) && min_int > 0) {
    steps <- c(steps, paste0("Minimum Intensity Filter: threshold = ", min_int,
                             ", min samples = ", input$min_samples_above_intensity %||% 1))
  } else {
    steps <- c(steps, "Minimum Intensity Filter: disabled")
  }
  imp <- preprocessing_params$imputation_method %||% "none"
  if (imp == "none") {
    steps <- c(steps, "Missing Value Imputation: none (rows with missing values will be removed)")
  } else {
    steps <- c(steps, paste0("Missing Value Imputation: ", imp))
  }
  if (isTRUE(preprocessing_params$batch_performed)) {
    steps <- c(steps, "Batch Correction (ComBat): applied")
  } else {
    steps <- c(steps, "Batch Correction: not applied")
  }
  
  norm_choice <- input$heatmap_normalization %||% "none"
  if (norm_choice == "total") {
    steps <- c(steps, "Normalization: Total intensity normalization (baseline sample) applied")
    steps <- c(steps, "Data source: Normalized expression data (Norm_LFQ intensity columns)")
  } else {
    steps <- c(steps, "Normalization: No total intensity normalization applied (uses preprocessed raw intensity)")
    steps <- c(steps, "Data source: Preprocessed data (LFQ/Intensity columns)")
  }
  steps <- c(steps, "log2(Intensity + 1) transformation applied")
  steps <- c(steps, "Per-row Z-score normalization (scale)")
  
  step_tags <- lapply(seq_along(steps), function(i) {
    tagList(
      if (i > 1) tags$span(style = "font-size: 20px; color: #e67e22; margin: 0 8px;", "→"),
      tags$span(style = "background: #e8f0fe; padding: 6px 12px; border-radius: 15px; font-size: 13px;", steps[[i]])
    )
  })
  tags$details(
    tags$summary("Data processing steps for heatmap", style = "cursor: pointer; font-weight: bold; color: #2c3e50; margin-bottom: 10px;"),
    div(style = "display: flex; flex-wrap: wrap; align-items: center;", do.call(tagList, step_tags))
  )
})

message("[DEBUG] heatmap_plot.R fully loaded (Arial font)")