# server/preprocessing_core.R
message("[DEBUG] preprocessing_core.R loaded - minvalue sample shows full underscored column name")

preprocessing_params <- reactiveValues(
  imputation_method = NULL,
  last_run_time = NULL,
  inf_filtered_count = 0,
  inf_filtered_proteins = character(0),
  batch_performed = FALSE,
  batch_corrected_cols = NULL,
  batch_uncorrected_cols = NULL,
  batch_match_summary = NULL,
  pre_batch_data = NULL,
  post_batch_data = NULL,
  missing_filter_mode = "group",
  missing_filtered_by_group = 0,
  knn_k = 10,
  min_value = 1e-4,
  quantile_prob = 0.01,
  intensity_min_samples = 1,
  missing_filter_fallback = FALSE,
  missing_filter_fallback_unmatched = 0,
  intensity_type_used = NULL,
  global_min_intensity = NULL,
  global_min_protein = NULL,
  global_min_sample = NULL,
  fill_factor = NULL,
  actual_fill_value = NULL
)

processed_data <- eventReactive(input$run_preprocessing, {
  showNotification("Running preprocessing...", type = "message", duration = NULL, id = "preprocess_notif")
  tryCatch({
    data <- expression_data()
    message("[DEBUG] processed_data: starting preprocessing")
    
    preprocessing_params$intensity_type_used <- input$intensity_type
    preprocessing_params$missing_filter_mode <- "group"
    preprocessing_params$missing_filter_fallback <- FALSE
    preprocessing_params$missing_filter_fallback_unmatched <- 0
    
    # 1. 缺失值过滤
    data <- apply_missing_filter(
      data, 
      threshold = input$max_missing_fraction,
      mode = "group",
      sample_info = rv$sample_info,
      sample_names_short = rv$sample_names
    )
    if (isTRUE(attr(data, "fallback_triggered"))) {
      preprocessing_params$missing_filter_fallback <- TRUE
      preprocessing_params$missing_filter_fallback_unmatched <- attr(data, "unmatched_count")
    }
    if (nrow(data) == 0) stop("No proteins left after missing filter.")
    
    # 2. Inf 过滤
    max_int <- apply(data, 1, max, na.rm = TRUE)
    keep_finite <- is.finite(max_int)
    preprocessing_params$inf_filtered_count <- sum(!keep_finite)
    preprocessing_params$inf_filtered_proteins <- rownames(data)[!keep_finite]
    data <- data[keep_finite, , drop = FALSE]
    if (nrow(data) == 0) stop("No proteins left after Inf filter.")
    
    # 3. 缺失值填补
    imp_method <- input$imputation_method
    if (is.null(imp_method) || imp_method == "none") {
      preprocessing_params$imputation_method <- "none"
      preprocessing_params$global_min_intensity <- NULL
      preprocessing_params$global_min_protein <- NULL
      preprocessing_params$global_min_sample <- NULL
      preprocessing_params$fill_factor <- NULL
      preprocessing_params$actual_fill_value <- NULL
    } else if (imp_method == "minvalue") {
      mat <- as.matrix(data)
      pos_indices <- which(!is.na(mat) & mat > 0, arr.ind = TRUE)
      if (nrow(pos_indices) == 0) {
        global_min <- 1e-4
        min_protein <- "N/A"
        min_sample <- "N/A"
      } else {
        vals <- mat[pos_indices]
        min_idx <- which.min(vals)
        min_row <- pos_indices[min_idx, "row"]
        min_col <- pos_indices[min_idx, "col"]
        global_min <- vals[min_idx]
        
        # 蛋白 ID
        protein_id_raw <- rownames(data)[min_row]
        if (!is.null(rv$clean_data) && "Master protein IDs" %in% colnames(rv$clean_data)) {
          if (suppressWarnings(!is.na(as.numeric(protein_id_raw)))) {
            idx <- as.integer(protein_id_raw)
            clean_ids <- rv$clean_data$`Master protein IDs`
            if (idx >= 1 && idx <= length(clean_ids)) {
              protein_id_raw <- clean_ids[idx]
            }
          }
        }
        min_protein <- protein_id_raw
        
        # 样本完整列名（已经是下划线格式）
        sample_long <- colnames(data)[min_col]
        min_sample <- sample_long   # 直接使用，如 "LFQ intensity WT_3"
      }
      factor_val <- input$minvalue_fixed
      actual_fill <- global_min * factor_val
      
      message(sprintf("[DEBUG] processed_data: minvalue gmin=%g from protein '%s', sample '%s'", 
                      global_min, min_protein, min_sample))
      
      data_matrix <- as.matrix(data)
      data_matrix[is.na(data_matrix)] <- actual_fill
      data <- as.data.frame(data_matrix)
      
      preprocessing_params$imputation_method <- "minvalue"
      preprocessing_params$global_min_intensity <- global_min
      preprocessing_params$global_min_protein <- min_protein
      preprocessing_params$global_min_sample <- min_sample
      preprocessing_params$fill_factor <- factor_val
      preprocessing_params$actual_fill_value <- actual_fill
    } else {
      preprocessing_params$knn_k <- input$knn_k
      preprocessing_params$min_value <- input$minvalue_fixed
      preprocessing_params$quantile_prob <- input$quantile_prob
      message("[DEBUG] processed_data: running imputation with method = ", imp_method)
      data <- impute_missing_values(data, method = imp_method,
                                    k = input$knn_k,
                                    min_value = input$minvalue_fixed,
                                    quantile_prob = input$quantile_prob)
      actual_method <- attr(data, "actual_method")
      if (is.null(actual_method)) actual_method <- imp_method
      preprocessing_params$imputation_method <- actual_method
      preprocessing_params$global_min_intensity <- NULL
      preprocessing_params$global_min_protein <- NULL
      preprocessing_params$global_min_sample <- NULL
      preprocessing_params$fill_factor <- NULL
      preprocessing_params$actual_fill_value <- NULL
    }
    
    preprocessing_params$batch_performed <- FALSE
    preprocessing_params$batch_corrected_cols <- NULL
    preprocessing_params$last_run_time <- Sys.time()
    removeNotification("preprocess_notif")
    showNotification("Preprocessing completed!", type = "message", duration = 3)
    return(data)
  }, error = function(e) {
    removeNotification("preprocess_notif")
    showNotification(paste("Preprocessing failed:", e$message), type = "error", duration = 10)
    return(NULL)
  })
})

observeEvent(input$expression_file, { preprocessing_params$batch_performed <- FALSE })
observeEvent(input$intensity_type, { preprocessing_params$batch_performed <- FALSE })

output$pre_processed_table <- DT::renderDT({
  req(processed_data())
  df <- processed_data()
  ids <- rownames(df)
  if (suppressWarnings(all(!is.na(as.numeric(ids))))) {
    if (!is.null(rv$clean_data) && "Master protein IDs" %in% colnames(rv$clean_data)) {
      original_ids <- rv$clean_data$`Master protein IDs`
      idx <- as.integer(ids)
      if (max(idx, na.rm = TRUE) <= length(original_ids)) ids <- original_ids[idx]
    }
  }
  df <- cbind(`Master Protein ID` = ids, df)
  rownames(df) <- NULL
  DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE, searchHighlight = TRUE, server = TRUE),
                rownames = FALSE, filter = "top")
})

output$preprocessing_done <- reactive({ !is.null(processed_data()) })
outputOptions(output, "preprocessing_done", suspendWhenHidden = FALSE)

output$imputation_skipped <- reactive({
  is.null(preprocessing_params$imputation_method) || preprocessing_params$imputation_method == "none"
})
outputOptions(output, "imputation_skipped", suspendWhenHidden = FALSE)

output$intensity_info <- renderPrint({ cat("Minimum intensity filter disabled.\n") })

output$preprocessing_steps_summary <- renderPrint({
  proc_df <- processed_data()
  if (is.null(proc_df)) {
    cat("Please run preprocessing first.\n")
    return()
  }
  cat("Preprocessing performed at:", format(preprocessing_params$last_run_time, "%Y-%m-%d %H:%M:%S"), "\n\n")
  cat("1. Missing Filter (within groups): Threshold", input$max_missing_fraction, "\n")
  cat("2. Inf Filter: Removed", preprocessing_params$inf_filtered_count, "\n")
  cat("3. Imputation:\n")
  method <- preprocessing_params$imputation_method
  if (!is.null(method) && method != "none") {
    cat("   Method:", method, "\n")
    if (method == "minvalue") {
      cat("   Factor:", format(preprocessing_params$fill_factor), "\n")
      cat("   Global min intensity:", format(preprocessing_params$global_min_intensity), "\n")
      cat("     (Protein:", preprocessing_params$global_min_protein, 
          ", Sample:", preprocessing_params$global_min_sample, ")\n")
      cat("   Actual fill value:", format(preprocessing_params$actual_fill_value), "\n")
      message("[DEBUG] steps_summary: displayed minvalue params with full underscored column name")
    }
  } else cat("   Skipped.\n")
  cat("Final data dimensions:", nrow(proc_df), "proteins,", ncol(proc_df), "samples\n")
  cat("Remaining missing values:", sum(is.na(proc_df)), "\n")
})

message("[DEBUG] preprocessing_core.R: minvalue sample name now uses full column name (underscored)")