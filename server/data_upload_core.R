# server/data_upload_core.R
message("[DEBUG] data_upload_core.R loading...")

# ================== 全局辅助函数（本文件专用） ==================
extract_sample_names <- function(cols) {
  short <- sub("^(LFQ intensity |Intensity )", "", cols, ignore.case = TRUE)
  still_full <- which(short == cols)
  if (length(still_full) > 0) {
    short[still_full] <- sub("^(LFQ[._]?intensity[._]?|Intensity[._]?)", "", cols[still_full], ignore.case = TRUE)
  }
  short <- sub("^[^[:alnum:]]+", "", short)
  short <- sub("[^[:alnum:]]+$", "", short)
  short <- gsub("-", ".", short)
  short
}

standardize_sample_name <- function(x) {
  if (is.null(x) || length(x) == 0) return(character(0))
  original <- x
  x <- as.character(x)
  x <- gsub("[-_]+", ".", x)
  x <- gsub("\\s+", ".", x)
  x <- gsub("\\.+", ".", x)
  x <- gsub("^\\.", "", x)
  x <- gsub("\\.$", "", x)
  if (length(x) > 0) {
    message("[DEBUG] standardize_sample_name: transformed first 3: ",
            paste(head(original, 3), collapse = ", "), " -> ",
            paste(head(x, 3), collapse = ", "))
  }
  return(x)
}

get_raw_prefix <- function(type = input$intensity_type) {
  if (type == "LFQ") "LFQ intensity " else "Intensity "
}

get_norm_prefix <- function(type = input$intensity_type) {
  if (type == "LFQ") "Norm_LFQ intensity " else "Norm_Intensity "
}

extract_group_name <- function(x) {
  x <- sub("^Norm_([A-Za-z]+ intensity) ", "", x)
  sapply(strsplit(x, "-"), `[`, 1)
}

get_group_colors <- function(groups) {
  if (length(groups) == 0) return(character(0))
  pal <- c("#FFB6C1","#90EE90","#87CEEB","#DDA0DD","#FFD700",
           "#FFA07A","#98FB98","#B0C4DE","#FFB347","#C9A0DC")
  if (length(groups) > length(pal)) {
    pal <- colorRampPalette(pal)(length(groups))
  }
  setNames(pal[1:length(groups)], groups)
}

# ================== 数据上传与预处理 ==================
output$download_sample_template <- downloadHandler(
  filename = function() { "sample_info_template.xlsx" },
  content = function(file) {
    prefix <- get_raw_prefix()
    template <- data.frame(
      SampleName = paste0(prefix, c("L2.1.1", "L2.1.2", "L2.1.3", "L2.2.1", "L2.2.2", "L2.2.3")),
      Group = c("Control", "Control", "Control", "Treatment", "Treatment", "Treatment"),
      Batch = c("Batch1", "Batch1", "Batch2", "Batch1", "Batch2", "Batch2"),
      Note = c("", "", "", "", "", "")
    )
    writexl::write_xlsx(template, file)
  }
)

read_sample_info <- function(file_path) {
  ext <- tools::file_ext(file_path)
  if (ext %in% c("csv", "txt")) {
    df <- read.csv(file_path, header = TRUE, row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
  } else if (ext %in% c("xlsx", "xls")) {
    df <- as.data.frame(readxl::read_excel(file_path, col_names = TRUE))
  } else {
    stop("Unsupported file format for sample info.")
  }
  rownames(df) <- as.character(df[[1]])
  df <- df[, -1, drop = FALSE]
  return(df)
}

cached_sample_info <- reactiveValues(LFQ = NULL, Intensity = NULL)

observeEvent(input$intensity_type, {
  if (!is.null(rv$sample_info)) {
    if (input$intensity_type == "LFQ") cached_sample_info$Intensity <- rv$sample_info
    else cached_sample_info$LFQ <- rv$sample_info
  }
  target_cache <- if (input$intensity_type == "LFQ") cached_sample_info$LFQ else cached_sample_info$Intensity
  if (!is.null(target_cache)) rv$sample_info <- target_cache
  else { rv$sample_info <- NULL; reset("sample_info_file") }
  if (!is.null(rv$raw_data)) {
    new_cols <- grep(paste0("^", get_raw_prefix()), colnames(rv$raw_data), value = TRUE)
    if (length(new_cols) > 0) {
      rv$lfq_cols <- new_cols
      rv$sample_names <- extract_sample_names(new_cols)
      updateSelectInput(session, "baseline_sample", choices = c("Auto", rv$sample_names), selected = "Auto")
    } else {
      showNotification("No matching intensity columns found for the selected type.", type = "error", duration = 5)
    }
  }
  
  updateRadioButtons(session, "heatmap_data_source", selected = "LFQ")
  updateRadioButtons(session, "heatmap_protein_mode", selected = "top_n")
  updateNumericInput(session, "heatmap_top_n", value = 20)
  updateTextAreaInput(session, "heatmap_custom_ids", value = "")
  heatmap_raw_groups(NULL)
  data_changed_trigger(data_changed_trigger() + 1)
  message("[DEBUG] intensity type changed, data_changed_trigger increased to ", data_changed_trigger())
}, ignoreInit = TRUE)

observeEvent(input$expression_file, {
  req(input$expression_file)
  rv$raw_data <- NULL; rv$clean_data <- NULL; rv$lfq_cols <- NULL
  rv$sample_names <- NULL; rv$groups <- list(); rv$comparisons <- list()
  rv$analysis_results <- NULL
  rv$comp_id_counter <- 0; rv$group_id_counter <- 0
  rv$group_id_map <- list(); rv$current_profile_protein <- NULL
  rv$batch_vector <- NULL; rv$sample_info <- NULL
  manual_sort_active(FALSE)
  updateSelectInput(session, "comp_treat", choices = character(0))
  updateSelectInput(session, "comp_ctrl", choices = character(0))
  updateSelectInput(session, "selected_comparison", choices = character(0))
  updateSelectInput(session, "baseline_sample", choices = c("Auto"), selected = "Auto")
  updateSelectInput(session, "batch_ref_group", choices = character(0))
  updateSelectizeInput(session, "venn_comparisons_select", choices = character(0), selected = character(0))
  updateCheckboxGroupInput(session, "venn_comparisons_checkbox", choices = character(0), selected = character(0))
  
  updateSelectInput(session, "missing_filter_mode", selected = "global")
  updateSliderInput(session, "max_missing_fraction", value = 0.5)
  updateNumericInput(session, "min_intensity", value = 1e5)
  updateNumericInput(session, "min_samples_above_intensity", value = 1)
  updateSelectInput(session, "imputation_method", selected = "quantile")   # 默认 quantile
  updateNumericInput(session, "knn_k", value = 10)
  updateCheckboxInput(session, "perform_batch_correction", value = FALSE)
  updateNumericInput(session, "fc_up", value = 1.2)
  updateNumericInput(session, "fc_down", value = 0.84)
  updateSelectInput(session, "p_cut", selected = "0.05")
  updateNumericInput(session, "min_treat_valid", value = 2)
  updateNumericInput(session, "min_ctrl_valid", value = 2)
  updateNumericInput(session, "min_rep_ttest", value = 2)
  updateNumericInput(session, "min_rep_inc", value = 2)
  updateNumericInput(session, "min_rep_dec", value = 2)
  updateNumericInput(session, "min_unique_pep", value = 2)
  
  updateRadioButtons(session, "heatmap_data_source", selected = "LFQ")
  updateRadioButtons(session, "heatmap_protein_mode", selected = "top_n")
  updateNumericInput(session, "heatmap_top_n", value = 20)
  updateTextAreaInput(session, "heatmap_custom_ids", value = "")
  heatmap_raw_groups(NULL)
  data_changed_trigger(data_changed_trigger() + 1)
  message("[DEBUG] expression file uploaded, data_changed_trigger increased to ", data_changed_trigger())
  
  tryCatch({
    file_path <- input$expression_file$datapath
    data <- fread(file_path, sep = "\t", stringsAsFactors = FALSE, data.table = FALSE, check.names = FALSE, colClasses = "character")
    for (cn in names(data)) {
      if (grepl("^(LFQ intensity |Intensity )", cn)) data[[cn]] <- as.numeric(data[[cn]])
    }
    lfq_cols <- grep(paste0("^", get_raw_prefix()), colnames(data), value = TRUE)
    if (length(lfq_cols) == 0) {
      showNotification("No matching intensity columns found. Please check intensity type.", type = "error", duration = 5)
      return()
    }
    sample_names <- extract_sample_names(lfq_cols)
    clean_data <- data
    if ("Reverse" %in% colnames(clean_data)) clean_data <- filter(clean_data, is.na(Reverse) | Reverse != "+")
    if ("Potential contaminant" %in% colnames(clean_data)) clean_data <- filter(clean_data, is.na(`Potential contaminant`) | `Potential contaminant` != "+")
    if ("Protein IDs" %in% colnames(clean_data)) clean_data <- filter(clean_data, !grepl("^CON_", `Protein IDs`))
    if ("Protein IDs" %in% colnames(clean_data)) clean_data <- mutate(clean_data, `Master protein IDs` = sub(";.*", "", `Protein IDs`), .after = `Majority protein IDs`)
    rv$raw_data <- data
    if ("Protein IDs" %in% colnames(rv$raw_data)) rv$raw_data$`Master protein IDs` <- sub(";.*", "", rv$raw_data$`Protein IDs`)
    rv$clean_data <- clean_data
    rv$lfq_cols <- lfq_cols
    rv$sample_names <- sample_names
    updateSelectInput(session, "baseline_sample", choices = c("Auto", sample_names), selected = "Auto")
    showNotification("Expression matrix uploaded successfully!", type = "message", duration = 3)
  }, error = function(e) {
    showNotification(paste("Error reading expression file:", e$message), type = "error", duration = 5)
  })
})

observeEvent(input$sample_info_file, {
  req(input$sample_info_file)
  tryCatch({
    df <- read_sample_info(input$sample_info_file$datapath)
    rv$sample_info <- df
    if (input$intensity_type == "LFQ") cached_sample_info$LFQ <- df
    else cached_sample_info$Intensity <- df
    showNotification("Sample info uploaded successfully!", type = "message", duration = 3)
  }, error = function(e) {
    showNotification(paste("Error reading sample info:", e$message), type = "error", duration = 5)
  })
})

sample_match_validation <- reactive({
  if (is.null(rv$lfq_cols) || is.null(rv$sample_info)) {
    return(list(status = "waiting", message = "Please upload both expression matrix and sample information.",
                matched = character(0), unmatched_info = character(0), unmatched_expr = character(0)))
  }
  expr_col_full <- rv$lfq_cols
  info_names_full <- rownames(rv$sample_info)
  expr_std <- standardize_sample_name(expr_col_full)
  info_std <- standardize_sample_name(info_names_full)
  matched_expr <- expr_col_full[expr_std %in% info_std]
  unmatched_info_full <- info_names_full[!info_std %in% expr_std]
  unmatched_expr_full <- expr_col_full[!expr_std %in% info_std]
  matched <- extract_sample_names(matched_expr)
  unmatched_info <- extract_sample_names(unmatched_info_full)
  unmatched_expr <- extract_sample_names(unmatched_expr_full)
  if (length(unmatched_info) == 0 && length(unmatched_expr) == 0) {
    return(list(status = "success", message = paste0("All ", length(matched), " samples are successfully matched!"),
                matched = matched, unmatched_info = character(0), unmatched_expr = character(0)))
  } else {
    return(list(status = "warning", message = paste0(length(matched), " samples matched. ",
                                                     length(unmatched_info), " sample(s) in info but not in expression; ",
                                                     length(unmatched_expr), " sample(s) in expression but not in info."),
                matched = matched, unmatched_info = unmatched_info, unmatched_expr = unmatched_expr))
  }
}) %>% bindCache(rv$lfq_cols, rv$sample_info)

output$sample_match_hint <- renderUI({
  req(rv$sample_info)
  div(style = "margin-top: 8px; padding: 8px 12px; background: #fff3cd; border-radius: 6px; color: #856404; font-weight: bold;",
      icon("info-circle"), " Green highlighted samples are matched with the uploaded sample info. Samples without fill color are not matched.")
})

output$upload_preview <- DT::renderDataTable({
  message("[DEBUG] upload_preview: rv$lfq_cols length = ", length(rv$lfq_cols))
  req(rv$clean_data, rv$lfq_cols)
  df <- rv$clean_data[, rv$lfq_cols, drop = FALSE]
  DT::datatable(df,
                options = list(pageLength = 10, scrollX = TRUE),
                rownames = FALSE)
})

output$sample_info_preview <- DT::renderDataTable({
  message("[DEBUG] sample_info_preview triggered")
  req(rv$sample_info)
  df <- rv$sample_info
  df_display <- data.frame(SampleName = rownames(df), df, check.names = FALSE, stringsAsFactors = FALSE)
  DT::datatable(df_display,
                options = list(pageLength = 10, scrollX = TRUE),
                rownames = FALSE)
})

output$data_summary_ui <- renderUI({
  req(rv$raw_data)
  type_label <- if (input$intensity_type == "LFQ") "LFQ intensity" else "Intensity"
  div(
    p(strong("Dimensions:"), sprintf("%d rows × %d columns", nrow(rv$raw_data), ncol(rv$raw_data))),
    p(strong(paste(type_label, "columns:", sep = " ")), length(rv$lfq_cols)),
    p(strong("Samples:"), length(rv$sample_names))
  )
})

output$detected_samples_ui <- renderUI({
  req(rv$sample_names)
  type_label <- if (input$intensity_type == "LFQ") "LFQ intensity" else "Intensity"
  samples <- rv$sample_names
  validation <- sample_match_validation()
  matched <- if (validation$status %in% c("success", "warning")) validation$matched else character(0)
  tagList(
    h4(icon("vial"), paste(" Detected Samples (", type_label, ")", sep = "")),
    div(style = "max-height: 200px; overflow-y: auto;",
        lapply(samples, function(s) {
          if (s %in% matched) {
            div(class = "sample-item", style = "background: #d4edda; border-color: #c3e6cb;", icon("vial"), " ", s)
          } else {
            div(class = "sample-item", icon("vial"), " ", s)
          }
        })
    )
  )
})

observeEvent(input$reset_all, {
  raw <- rv$raw_data; clean <- rv$clean_data; lfq <- rv$lfq_cols; sn <- rv$sample_names; si <- rv$sample_info
  rv$groups <- list(); rv$comparisons <- list(); rv$analysis_results <- NULL
  rv$pending_duplicate <- NULL; rv$reset_counter <- rv$reset_counter + 1
  rv$comp_id_counter <- 0; rv$group_id_counter <- 0; rv$group_id_map <- list()
  rv$current_profile_protein <- NULL; rv$batch_vector <- NULL; rv$sample_info <- NULL
  manual_sort_active(FALSE)
  cached_sample_info$LFQ <- NULL; cached_sample_info$Intensity <- NULL
  for (name in names(subplot_old_values)) subplot_old_values[[name]] <- NULL
  updateSelectInput(session, "comp_treat", choices = character(0))
  updateSelectInput(session, "comp_ctrl", choices = character(0))
  updateSelectInput(session, "selected_comparison", choices = character(0))
  updateSelectInput(session, "batch_ref_group", choices = character(0))
  updateSelectizeInput(session, "venn_comparisons_select", choices = character(0), selected = character(0))
  updateCheckboxGroupInput(session, "venn_comparisons_checkbox", choices = character(0), selected = character(0))
  updateNumericInput(session, "fc_up", value = 1.2)
  updateNumericInput(session, "fc_down", value = 0.84)
  updateSelectInput(session, "p_cut", selected = "0.05")
  updateNumericInput(session, "min_treat_valid", value = 2)
  updateNumericInput(session, "min_ctrl_valid", value = 2)
  updateNumericInput(session, "min_rep_ttest", value = 2)
  updateNumericInput(session, "min_rep_inc", value = 2)
  updateNumericInput(session, "min_rep_dec", value = 2)
  updateNumericInput(session, "min_unique_pep", value = 2)
  updateNumericInput(session, "point_size", value = 4)
  cols <- default_colors()
  colourpicker::updateColourInput(session, "color_up", value = cols$Up)
  colourpicker::updateColourInput(session, "color_down", value = cols$Down)
  colourpicker::updateColourInput(session, "color_increase", value = cols$Increase)
  colourpicker::updateColourInput(session, "color_decrease", value = cols$Decrease)
  colourpicker::updateColourInput(session, "color_ns", value = cols$NS)
  updateRadioButtons(session, "stat_method", selected = "t-test")
  updateNumericInput(session, "replicate_fill_all", value = 2)
  updateTextInput(session, "download_single_title", value = "")
  updateTextInput(session, "combined_plot_title", value = "Combined Volcano Plots")
  updateTextInput(session, "single_plot_title", value = "")
  shinyjs::reset("plot_format")
  updateTextInput(session, "plot_width", value = "10")
  updateTextInput(session, "plot_height", value = "8")
  sub_ids <- grep("^subplot_title_", names(input), value = TRUE)
  for (sid in sub_ids) updateTextInput(session, sid, value = "")
  updateRadioButtons(session, "heatmap_data_source", selected = "LFQ")
  heatmap_raw_groups(NULL)
  # 重置填补方法为 quantile
  updateSelectInput(session, "imputation_method", selected = "quantile")
  rv$raw_data <- raw; rv$clean_data <- clean; rv$lfq_cols <- lfq; rv$sample_names <- sn; rv$sample_info <- si
  if (!is.null(sn) && length(sn) > 0) {
    updateSelectInput(session, "baseline_sample", choices = c("Auto", sn), selected = "Auto")
  } else {
    updateSelectInput(session, "baseline_sample", choices = c("Auto"), selected = "Auto")
  }
  showNotification("All settings reset.", type = "message", duration = 2)
})

get_base_sample <- function() {
  user_sel <- input$baseline_sample
  if (!is.null(user_sel) && user_sel != "Auto" && user_sel %in% rv$sample_names) return(user_sel)
  groups <- rv$groups
  if (length(groups) > 0) {
    wt_names <- grep("WT|Control|CK", names(groups), value = TRUE, ignore.case = TRUE)
    if (length(wt_names) > 0 && length(groups[[wt_names[1]]]) > 0) return(groups[[wt_names[1]]][1])
  }
  comps <- rv$comparisons
  if (length(comps) > 0) {
    ctrl_group <- comps[[1]]$ctrl
    if (!is.null(ctrl_group) && ctrl_group %in% names(groups) && length(groups[[ctrl_group]]) > 0)
      return(groups[[ctrl_group]][1])
  }
  if (length(rv$sample_names) > 0) return(rv$sample_names[1])
  NULL
}

current_baseline <- reactive({ get_base_sample() })

raw_totals <- reactive({
  req(rv$raw_data, rv$lfq_cols)
  totals <- sapply(rv$lfq_cols, function(col) sum(as.numeric(rv$raw_data[[col]]), na.rm = TRUE))
  names(totals) <- extract_sample_names(rv$lfq_cols)
  totals
})

get_analysis_matrix <- reactive({
  if (is.null(preprocessing_params$intensity_type_used) || 
      preprocessing_params$intensity_type_used != input$intensity_type) {
    message("[DEBUG] get_analysis_matrix: processed data unavailable or intensity type mismatch")
    # 不再显示右下角通知，避免在 Data Quality 页面弹出警告
    return(NULL)
  }
  
  proc <- tryCatch(processed_data(), error = function(e) NULL)
  if (!is.null(proc)) {
    message("[DEBUG] get_analysis_matrix: returning processed data")
    return(proc)
  } else {
    message("[DEBUG] get_analysis_matrix: processed_data is NULL")
    return(NULL)
  }
})

# ---------- norm_data_before_batch ----------
norm_data_before_batch <- reactive({
  mat <- get_analysis_matrix()
  if (is.null(mat)) {
    showNotification("No expression data available.", type = "error")
    return(NULL)
  }
  base_sample <- current_baseline()
  if (is.null(base_sample)) {
    showNotification("Unable to determine baseline sample.", type = "error")
    return(NULL)
  }
  sample_short <- extract_sample_names(colnames(mat))
  base_idx <- which(sample_short == base_sample)
  if (length(base_idx) == 0) {
    showNotification(paste0("Baseline sample '", base_sample, "' not found."), type = "error")
    return(NULL)
  }
  base_sum <- sum(mat[, base_idx], na.rm = TRUE)
  if (base_sum <= 0) {
    showNotification("Baseline sample total intensity is zero.", type = "error")
    return(NULL)
  }
  norm_mat <- mat
  for (i in seq_len(ncol(mat))) {
    s <- sum(mat[, i], na.rm = TRUE)
    if (s > 0) norm_mat[, i] <- mat[, i] * base_sum / s
    else norm_mat[, i] <- mat[, i]
  }
  norm_prefix <- get_norm_prefix()
  colnames(norm_mat) <- paste0(norm_prefix, sample_short)
  
  # 获取原始强度列
  orig_prefix <- paste0("Original_", get_raw_prefix())   # e.g. "Original_LFQ intensity "
  orig_mat <- mat
  colnames(orig_mat) <- paste0(orig_prefix, sample_short)
  
  rn <- rownames(norm_mat)
  protein_ids <- rn
  if (suppressWarnings(all(!is.na(as.numeric(rn))))) {
    clean_ids <- rv$clean_data$`Master protein IDs`
    if (!is.null(clean_ids) && length(clean_ids) >= max(as.integer(rn))) {
      protein_ids <- clean_ids[as.integer(rn)]
    }
  }
  
  norm_df <- as.data.frame(norm_mat)
  rownames(norm_df) <- protein_ids
  norm_df$`Master protein IDs` <- protein_ids
  
  clean <- rv$clean_data
  extra_cols <- intersect(c("Protein IDs", "Majority protein IDs", "Unique peptides"), colnames(clean))
  if (length(extra_cols) > 0) {
    idx <- match(protein_ids, clean$`Master protein IDs`)
    for (col in extra_cols) {
      norm_df[[col]] <- clean[[col]][idx]
    }
  }
  
  orig_df <- as.data.frame(orig_mat)
  for (cn in colnames(orig_df)) {
    norm_df[[cn]] <- orig_df[[cn]]
  }
  
  if ("Fasta headers" %in% colnames(norm_df)) {
    norm_df <- norm_df[, setdiff(colnames(norm_df), "Fasta headers"), drop = FALSE]
  }
  
  desired_annotation <- c("Protein IDs", "Majority protein IDs", "Master protein IDs", "Unique peptides")
  existing_annotation <- intersect(desired_annotation, colnames(norm_df))
  orig_cols <- grep(paste0("^", orig_prefix), colnames(norm_df), value = TRUE)
  norm_cols_all <- grep(paste0("^", norm_prefix), colnames(norm_df), value = TRUE)
  other_cols <- setdiff(colnames(norm_df), c(existing_annotation, orig_cols, norm_cols_all))
  new_order <- c(existing_annotation, orig_cols, norm_cols_all, other_cols)
  norm_df <- norm_df[, new_order, drop = FALSE]
  
  message("[DEBUG] norm_data_before_batch: final column order (first 8): ",
          paste(head(colnames(norm_df), 8), collapse = ", "))
  norm_df
})

norm_totals <- reactive({
  nd <- norm_data_before_batch()
  if (is.null(nd)) return(NULL)
  norm_prefix <- get_norm_prefix()
  norm_cols <- grep(paste0("^", norm_prefix), colnames(nd), value = TRUE)
  if (length(norm_cols) == 0) return(NULL)
  totals <- sapply(norm_cols, function(col) sum(as.numeric(nd[[col]]), na.rm = TRUE))
  names(totals) <- gsub(paste0("^", norm_prefix), "", norm_cols)
  totals
})

norm_data_full <- reactive({
  norm_data_before_batch()
})

output$norm_comparison_plot <- renderPlotly({
  req(raw_totals())
  raw <- raw_totals()
  nrt <- tryCatch(norm_totals(), error = function(e) NULL)
  samples <- names(raw)
  baseline <- current_baseline()
  p <- plot_ly()
  p <- add_bars(p, x = samples, y = raw[samples], name = "Raw",
                marker = list(color = "steelblue"),
                hovertemplate = paste0("Sample: %{x}<br>Raw: %{y:.0f}<extra></extra>"))
  if (!is.null(nrt)) {
    p <- add_bars(p, x = samples, y = nrt[samples], name = "Normalized",
                  marker = list(color = "darkorange"),
                  hovertemplate = paste0("Sample: %{x}<br>Normalized: %{y:.0f}<extra></extra>"))
    title_text <- paste0("Total Intensity: Raw vs Normalized (Baseline: ", baseline, ")")
  } else {
    title_text <- paste0("Total Intensity (Raw only) (Baseline: ", baseline %||% "N/A", ")")
  }
  p %>% layout(title = title_text, yaxis = list(title = "Total Intensity"),
               xaxis = list(tickangle = -45), legend = list(title = list(text = "Type")), barmode = "group")
})

output$upload_status_ui <- renderUI({
  if (is.null(rv$raw_data)) {
    div(class = "status-badge status-warning", icon("exclamation-triangle"), " No expression file uploaded")
  } else {
    div(class = "status-badge status-success", icon("check-circle"), " Expression matrix uploaded!")
  }
})

observeEvent(input$reset_color, {
  cols <- default_colors()
  colourpicker::updateColourInput(session, "color_up", value = cols$Up)
  colourpicker::updateColourInput(session, "color_down", value = cols$Down)
  colourpicker::updateColourInput(session, "color_increase", value = cols$Increase)
  colourpicker::updateColourInput(session, "color_decrease", value = cols$Decrease)
  colourpicker::updateColourInput(session, "color_ns", value = cols$NS)
  showNotification("Colors reset to defaults.", type = "message", duration = 2)
})

expression_data <- reactive({
  message("[DEBUG] expression_data (from data_upload_core.R) triggered")
  req(rv$clean_data)
  
  if (is.null(rv$lfq_cols) || length(rv$lfq_cols) == 0) {
    message("[DEBUG] expression_data: no lfq_cols found, will validate")
    validate(need(FALSE, "No intensity columns found. Please upload data first."))
  }
  
  df <- rv$clean_data
  if (!"Master protein IDs" %in% colnames(df)) {
    message("[DEBUG] expression_data: Master protein IDs column missing")
    validate(need(FALSE, "Master protein IDs column not found in cleaned data."))
  }
  
  rownames(df) <- as.character(df$`Master protein IDs`)
  df <- df[, rv$lfq_cols, drop = FALSE]
  df <- suppressWarnings(as.data.frame(lapply(df, as.numeric)))
  df[df == 0] <- NA
  
  if (ncol(df) == 0) {
    message("[DEBUG] expression_data: zero columns after subsetting")
    validate(need(FALSE, "No intensity columns found."))
  }
  if (nrow(df) == 0) {
    message("[DEBUG] expression_data: zero rows")
    validate(need(FALSE, "No protein rows found."))
  }
  
  message(sprintf("[DEBUG] expression_data (data_upload_core.R): returning %d proteins, %d samples", nrow(df), ncol(df)))
  df
})

message("[DEBUG] data_upload_core.R loaded successfully.")