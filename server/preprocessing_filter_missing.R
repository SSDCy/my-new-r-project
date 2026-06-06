# server/preprocessing_filter_missing.R
message("[DEBUG] preprocessing_filter_missing.R loaded (group mode default)")

# 样本信息Batch列预览
output$batch_info_preview <- renderPrint({
  req(rv$sample_info)
  col_data <- rv$sample_info[["Batch"]]
  unique_vals <- unique(col_data)
  cat("样本信息表“Batch”列概览：\n")
  cat("样本总数：", length(col_data), "\n")
  cat("唯一值数量：", length(unique_vals), "\n")
  cat("唯一值列表（前10个）：\n")
  if (length(unique_vals) <= 10) {
    cat(paste(unique_vals, collapse = "\n"), "\n")
  } else {
    cat(paste(head(unique_vals, 10), collapse = "\n"), "\n... 还有", length(unique_vals) - 10, "个\n")
  }
})

# 快速预设按钮
observeEvent(input$preset_missing_0.3, {
  updateSliderInput(session, "max_missing_fraction", value = 0.3)
})
observeEvent(input$preset_missing_0.5, {
  updateSliderInput(session, "max_missing_fraction", value = 0.5)
})
observeEvent(input$preset_missing_0.7, {
  updateSliderInput(session, "max_missing_fraction", value = 0.7)
})

# 分组匹配状态 UI（已简化，只显示分组信息）
output$filter_mode_group_match_ui <- renderUI({
  req(rv$sample_names)
  
  if (is.null(rv$sample_info) || !"Group" %in% colnames(rv$sample_info)) {
    return(div(style = "color: #d9534f;", 
               icon("times-circle"), " No sample info with Group column uploaded. Filtering will fall back to global mode."))
  }
  
  si <- rv$sample_info
  si$short <- extract_sample_names(rownames(si))
  matched <- sum(rv$sample_names %in% si$short)
  unmatched <- length(rv$sample_names) - matched
  
  if (unmatched > 0) {
    unmatched_samples <- setdiff(rv$sample_names, si$short)
    sample_list <- paste(head(unmatched_samples, 5), collapse = ", ")
    if (length(unmatched_samples) > 5) sample_list <- paste0(sample_list, " ... and ", length(unmatched_samples)-5, " more")
    return(div(style = "color: #f0ad4e;",
               icon("exclamation-triangle"), 
               paste0(" Warning: ", unmatched, " sample(s) not matched to sample info. ",
                      "Filtering will fall back to Global mode."),
               br(),
               tags$small("Unmatched samples: ", sample_list)))
  }
  
  if (!is.null(rv$groups) && length(rv$groups) > 0) {
    groups <- names(rv$groups)
    assigned_count <- sum(lengths(rv$groups))
    return(div(style = "color: #5cb85c;",
               icon("check-circle"),
               paste0(" Using user-defined groups: ", length(groups), " groups with ", assigned_count, " assigned samples.")))
  } else {
    groups <- unique(si$Group[match(rv$sample_names, si$short)])
    return(div(style = "color: #5cb85c;",
               icon("check-circle"),
               paste0(" Using groups from sample info: ", length(groups), " groups.")))
  }
})

# 缺失值过滤效果（只显示分组模式）
missing_filter_prediction_detail <- reactive({
  req(expression_data(), input$max_missing_fraction)
  data <- expression_data()
  threshold <- input$max_missing_fraction
  # 强制使用分组模式
  mode <- "group"
  
  detail <- list(
    mode = "Within Groups",
    total_proteins = nrow(data),
    total_samples = ncol(data),
    removed = NA,
    retained = NA,
    groups_detail = NULL
  )
  
  if (is.null(rv$sample_info) || !"Group" %in% colnames(rv$sample_info)) {
    detail$warning <- "No sample info with Group column available. Filtering will fall back to Global mode."
    return(detail)
  }
  
  si <- rv$sample_info
  si$short <- extract_sample_names(rownames(si))
  sample_short <- rv$sample_names
  idx <- match(sample_short, si$short)
  group_vec <- si$Group[idx]
  group_vec[is.na(idx)] <- "Unknown"
  
  keep <- rep(FALSE, nrow(data))
  groups_list <- unique(group_vec)
  groups_detail <- list()
  
  per_group_pass <- matrix(FALSE, nrow = nrow(data), ncol = length(groups_list))
  colnames(per_group_pass) <- groups_list
  
  for (i in seq_along(groups_list)) {
    g <- groups_list[i]
    cols_in_group <- which(group_vec == g)
    if (length(cols_in_group) == 0) next
    missing_frac_group <- rowMeans(is.na(data[, cols_in_group, drop = FALSE]))
    keep_g <- missing_frac_group <= threshold
    keep <- keep | keep_g
    per_group_pass[, i] <- keep_g
    groups_detail[[g]] <- list(
      n_samples = length(cols_in_group),
      proteins_kept_in_group = sum(keep_g)
    )
  }
  
  detail$retained <- sum(keep)
  detail$removed <- nrow(data) - detail$retained
  detail$groups_detail <- groups_detail
  detail$keep <- keep
  detail$per_group_pass <- per_group_pass
  detail$missing_rate <- rowMeans(is.na(data))
  
  if (sum(is.na(idx)) > 0) {
    detail$warning <- paste0(sum(is.na(idx)), " samples not matched, treated as group 'Unknown'.")
  }
  
  return(detail)
})

output$missing_filter_effect <- renderPrint({
  detail <- missing_filter_prediction_detail()
  req(detail)
  
  cat("Current filter mode: Within Groups\n")
  if (!is.null(detail$warning)) {
    cat("Warning: ", detail$warning, "\n")
  }
  cat("Evaluating missing rate within each group.\n")
  cat("Groups and sample counts:\n")
  if (!is.null(detail$groups_detail)) {
    for (g in names(detail$groups_detail)) {
      info <- detail$groups_detail[[g]]
      cat(sprintf("  %s: %d samples, %d proteins pass threshold within group\n", g, info$n_samples, info$proteins_kept_in_group))
    }
  }
  cat(sprintf("\nOverall predicted removal: %d proteins\n", detail$removed))
  cat(sprintf("Overall predicted retained: %d proteins\n", detail$retained))
  cat("\nNote: This prediction only considers missing rate, without prior Inf filter.\n")
})

# 导出缺失值过滤结果 Excel（分组模式）
output$download_missing_filter_excel <- downloadHandler(
  filename = function() {
    paste0("Missing_Filter_Result_", Sys.Date(), ".xlsx")
  },
  content = function(file) {
    shiny::withProgress(message = 'Exporting Missing Filter Excel...', value = 0, {
      incProgress(0.1, detail = "Preparing filter details...")
      message("[DEBUG] download_missing_filter_excel: starting export")
      req(expression_data())
      
      data <- expression_data()
      threshold <- input$max_missing_fraction
      mode <- "group"   # 强制分组
      
      detail <- missing_filter_prediction_detail()
      if (is.null(detail$keep)) {
        showNotification("Unable to compute filter details.", type = "error")
        return()
      }
      
      keep <- detail$keep
      row_names <- rownames(data)
      all_numeric <- suppressWarnings(all(!is.na(as.numeric(row_names))))
      if (all_numeric) {
        if (!is.null(rv$clean_data) && "Master protein IDs" %in% colnames(rv$clean_data)) {
          clean_ids <- rv$clean_data$`Master protein IDs`
          if (length(clean_ids) == nrow(data)) {
            protein_ids <- clean_ids
          } else {
            protein_ids <- row_names
          }
        } else {
          protein_ids <- row_names
        }
      } else {
        protein_ids <- row_names
      }
      
      retained_ids <- protein_ids[keep]
      filtered_ids <- protein_ids[!keep]
      
      group_extra_ids <- character(0)
      if (mode == "group") {
        missing_global <- rowMeans(is.na(data))
        global_keep <- missing_global <= threshold
        group_extra_ids <- protein_ids[keep & !global_keep]
      }
      
      incProgress(0.3, detail = "Creating workbook...")
      wb <- openxlsx::createWorkbook()
      
      retained_df <- data[which(keep), , drop = FALSE]
      retained_df <- cbind(ProteinID = retained_ids, retained_df, stringsAsFactors = FALSE)
      if (length(group_extra_ids) > 0) {
        retained_df$GroupExtra <- ifelse(retained_df$ProteinID %in% group_extra_ids, "Yes (Group saved, Global would remove)", "")
      }
      openxlsx::addWorksheet(wb, "Retained")
      openxlsx::writeData(wb, "Retained", retained_df)
      
      filtered_df <- data[which(!keep), , drop = FALSE]
      filtered_df <- cbind(ProteinID = filtered_ids, filtered_df, stringsAsFactors = FALSE)
      openxlsx::addWorksheet(wb, "Filtered_Out")
      openxlsx::writeData(wb, "Filtered_Out", filtered_df)
      
      incProgress(0.3, detail = "Writing group details...")
      if (mode == "group" && !is.null(detail$per_group_pass)) {
        per_group <- detail$per_group_pass
        group_names <- colnames(per_group)
        pass_matrix <- per_group
        rownames(pass_matrix) <- protein_ids
        
        retained_pass <- pass_matrix[retained_ids, , drop = FALSE]
        group_assign <- apply(retained_pass, 1, function(x) {
          g <- group_names[x]
          if (length(g) == 0) return("None")
          paste(g, collapse = ";")
        })
        
        group_detail_df <- data.frame(ProteinID = retained_ids, GroupPass = group_assign, stringsAsFactors = FALSE)
        openxlsx::addWorksheet(wb, "Group_Details")
        openxlsx::writeData(wb, "Group_Details", group_detail_df)
        
        openxlsx::addWorksheet(wb, "Legend")
        legend_text <- c(
          "Retained sheet: 'GroupExtra' column indicates proteins retained by Group mode but would be filtered out by Global mode.",
          "Group_Details sheet: 'GroupPass' column shows which group(s) the protein passes the threshold in."
        )
        openxlsx::writeData(wb, "Legend", data.frame(Info = legend_text))
      }
      
      incProgress(0.2, detail = "Saving workbook...")
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
      incProgress(0.1, detail = "Done")
      message("[DEBUG] download_missing_filter_excel: export completed")
    })
  }
)