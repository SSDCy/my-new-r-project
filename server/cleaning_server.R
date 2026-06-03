# server/cleaning_server.R
message("[DEBUG] cleaning_server.R loaded")

# 清洗统计（基于 rv$raw_data 和 rv$clean_data 的对比）
output$cleaning_summary <- renderPrint({
  message("[DEBUG] cleaning_summary called")
  req(rv$raw_data)
  cat("Data cleaning summary:\n")
  cat("=====================\n")
  
  raw <- rv$raw_data
  clean <- rv$clean_data
  
  total_raw <- nrow(raw)
  total_clean <- nrow(clean)
  
  # 计算被移除的各类蛋白
  reverse_removed <- if ("Reverse" %in% colnames(raw)) {
    raw$Reverse == "+" & !is.na(raw$Reverse)
  } else rep(FALSE, total_raw)
  contaminant_removed <- if ("Potential contaminant" %in% colnames(raw)) {
    raw$`Potential contaminant` == "+" & !is.na(raw$`Potential contaminant`)
  } else rep(FALSE, total_raw)
  con_removed <- if ("Protein IDs" %in% colnames(raw)) {
    grepl("^CON_", raw$`Protein IDs`)
  } else rep(FALSE, total_raw)
  
  # 实际删除的行：raw 中不在 clean 中的行（通过 Master protein IDs 对比）
  if ("Master protein IDs" %in% colnames(raw) && "Master protein IDs" %in% colnames(clean)) {
    removed_ids <- setdiff(raw$`Master protein IDs`, clean$`Master protein IDs`)
  } else if ("Protein IDs" %in% colnames(raw) && "Protein IDs" %in% colnames(clean)) {
    removed_ids <- setdiff(raw$`Protein IDs`, clean$`Protein IDs`)
  } else {
    removed_ids <- setdiff(rownames(raw), rownames(clean))
  }
  
  n_reverse <- sum(reverse_removed)
  n_contaminant <- sum(contaminant_removed)
  n_con <- sum(con_removed)
  n_total_removed <- length(removed_ids)
  
  cat("Proteins in uploaded file:", total_raw, "\n")
  cat("Proteins retained after cleaning:", total_clean, "\n")
  cat("Total proteins removed:", n_total_removed, "\n\n")
  
  cat("Breakdown:\n")
  cat("  - Reverse hits:", n_reverse, "\n")
  cat("  - Potential contaminants:", n_contaminant, "\n")
  cat("  - CON_ contaminants:", n_con, "\n")
  cat("\nNote: Some proteins may be removed by multiple filters.\n")
  
  message("[DEBUG] cleaning_summary: total_raw=", total_raw, ", total_clean=", total_clean, ", removed=", n_total_removed)
})

# 被移除的蛋白ID列表
cleaning_removed_ids <- reactive({
  req(rv$raw_data, rv$clean_data)
  raw <- rv$raw_data
  clean <- rv$clean_data
  if ("Master protein IDs" %in% colnames(raw) && "Master protein IDs" %in% colnames(clean)) {
    list(
      all_removed = setdiff(raw$`Master protein IDs`, clean$`Master protein IDs`),
      reverse = if ("Reverse" %in% colnames(raw)) {
        raw$`Master protein IDs`[raw$Reverse == "+" & !is.na(raw$Reverse)]
      } else character(0),
      contaminant = if ("Potential contaminant" %in% colnames(raw)) {
        raw$`Master protein IDs`[raw$`Potential contaminant` == "+" & !is.na(raw$`Potential contaminant`)]
      } else character(0),
      con = if ("Protein IDs" %in% colnames(raw)) {
        raw$`Master protein IDs`[grepl("^CON_", raw$`Protein IDs`)]
      } else character(0)
    )
  } else {
    list(all_removed = character(0), reverse = character(0), contaminant = character(0), con = character(0))
  }
})

output$cleaning_reverse_ids <- renderPrint({
  ids <- cleaning_removed_ids()$reverse
  if (length(ids) == 0) cat("No reverse hits found.\n")
  else cat(paste(ids, collapse = "\n"), "\n")
})

output$cleaning_contaminant_ids <- renderPrint({
  ids <- cleaning_removed_ids()$contaminant
  if (length(ids) == 0) cat("No potential contaminants found.\n")
  else cat(paste(ids, collapse = "\n"), "\n")
})

output$cleaning_con_ids <- renderPrint({
  ids <- cleaning_removed_ids()$con
  if (length(ids) == 0) cat("No CON_ contaminants found.\n")
  else cat(paste(ids, collapse = "\n"), "\n")
})

message("[DEBUG] cleaning_server.R fully loaded")