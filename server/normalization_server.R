# server/normalization_server.R
message("[DEBUG] normalization_server.R loaded")

# 更新基线样本选择器的选项
observeEvent(rv$sample_names, {
  message("[DEBUG] normalization_server: updating baseline_sample choices")
  samples <- rv$sample_names
  if (length(samples) > 0) {
    # 默认选中第一个包含 "Control" 或 "WT" 的样本
    default <- samples[1]
    ctrl <- grep("WT|Control|CK", samples, value = TRUE, ignore.case = TRUE)
    if (length(ctrl) > 0) default <- ctrl[1]
    updateSelectInput(session, "baseline_sample", choices = samples, selected = default)
    message("[DEBUG] normalization_server: baseline_sample choices set, default = ", default)
  }
})

# 显示当前基线信息
output$norm_baseline_info <- renderPrint({
  message("[DEBUG] norm_baseline_info called")
  req(input$baseline_sample)
  baseline <- input$baseline_sample
  totals <- raw_totals()
  if (!is.null(totals) && baseline %in% names(totals)) {
    cat("Baseline sample:", baseline, "\n")
    cat("Raw total intensity:", format(totals[baseline], big.mark = ",", scientific = FALSE), "\n")
  } else {
    cat("Baseline sample not found in expression data.\n")
  }
  message("[DEBUG] norm_baseline_info: baseline = ", baseline)
})

# 归一化前后总强度表格
output$norm_totals_table <- DT::renderDT({
  message("[DEBUG] norm_totals_table called")
  req(raw_totals())
  raw <- raw_totals()
  norm <- tryCatch(norm_totals(), error = function(e) NULL)
  
  df <- data.frame(
    Sample = names(raw),
    Raw_Total = raw,
    Normalized_Total = if (!is.null(norm)) norm[names(raw)] else NA_real_,
    stringsAsFactors = FALSE
  )
  df$Normalized_Total <- ifelse(is.na(df$Normalized_Total), "N/A", format(df$Normalized_Total, big.mark = ",", scientific = FALSE))
  df$Raw_Total <- format(df$Raw_Total, big.mark = ",", scientific = FALSE)
  
  message("[DEBUG] norm_totals_table: returning table with ", nrow(df), " rows")
  DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
})

message("[DEBUG] normalization_server.R fully loaded")