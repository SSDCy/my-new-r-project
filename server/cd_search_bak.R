# server/cd_search.R
message("[DEBUG] cd_search.R loading... (httr first, then curl, manual fallback)")

# ---------- 读取上传的 FASTA ----------
cd_fasta <- reactive({
  req(input$cd_search_fasta)
  file_path <- input$cd_search_fasta$datapath
  message("[DEBUG] cd_search: reading FASTA from ", file_path)
  if (!requireNamespace("seqinr", quietly = TRUE)) {
    showNotification("Package 'seqinr' is required.", type = "error")
    return(NULL)
  }
  tryCatch({
    seqs <- seqinr::read.fasta(file_path, seqtype = "AA", as.string = TRUE)
    ids <- names(seqs)
    sequences <- unlist(seqs)
    message("[DEBUG] cd_search: loaded ", length(ids), " sequences")
    data.frame(ID = ids, Sequence = sequences, stringsAsFactors = FALSE)
  }, error = function(e) {
    message("[ERROR] cd_search: ", conditionMessage(e))
    showNotification(paste("Failed to read FASTA:", conditionMessage(e)), type = "error")
    NULL
  })
})

# ---------- 自动搜索按钮 ----------
observeEvent(input$run_cd_search, {
  message("[DEBUG] cd_search: button clicked")
  fasta_df <- cd_fasta()
  if (is.null(fasta_df) || nrow(fasta_df) == 0) {
    showNotification("Please upload a FASTA file first.", type = "error")
    return()
  }
  fasta_lines <- paste0(">", fasta_df$ID, "\n", fasta_df$Sequence)
  fasta_text <- paste(fasta_lines, collapse = "\n")
  temp_fasta <- tempfile(fileext = ".fasta")
  writeLines(fasta_text, temp_fasta)
  message("[DEBUG] cd_search: temp FASTA at ", temp_fasta)
  
  progress <- shiny::Progress$new()
  progress$set(message = "Submitting to NCBI CD-Search...", value = 0)
  on.exit({
    if (file.exists(temp_fasta)) file.remove(temp_fasta)
    progress$close()
  })
  
  # ---- 方法一：使用 httr ----
  tryCatch({
    message("[DEBUG] cd_search: trying httr POST")
    res <- httr::POST(
      url = "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi",
      config = httr::config(ssl_verifyhost = FALSE, ssl_verifypeer = FALSE),
      body = list(
        queries = httr::upload_file(temp_fasta, type = "text/plain"),
        db = "cdd_delta",
        smode = "auto",
        useid1 = "true",
        format = "tsv"
      ),
      encode = "multipart"
    )
    httr::stop_for_status(res)
    text_content <- httr::content(res, as = "text", encoding = "UTF-8")
    con <- textConnection(text_content)
    cd_result <- read.csv(con, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    close(con)
    message("[DEBUG] cd_search: parsed result, nrow = ", nrow(cd_result))
    rv$cd_auto_result <- cd_result
    showNotification("CD-Search completed!", type = "message")
    progress$set(value = 1, detail = "Done")
    return()
  }, error = function(e) {
    message("[ERROR] cd_search (httr): ", conditionMessage(e))
  })
  
  # ---- 方法二：回退到系统 curl ----
  message("[DEBUG] cd_search: trying system curl")
  system_curl <- tryCatch({
    system("curl --version", ignore.stdout = TRUE, ignore.stderr = TRUE)
    TRUE
  }, error = function(e) FALSE)
  
  if (!system_curl) {
    showNotification("curl not found. Please use manual upload.", type = "error")
    return()
  }
  
  tryCatch({
    cmd <- sprintf('curl -k -s -X POST "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi" -F "queries=@%s;type=text/plain" -F "db=cdd_delta" -F "smode=auto" -F "useid1=true" -F "format=tsv"',
                   shQuote(temp_fasta))
    message("[DEBUG] cd_search: running curl command")
    res <- system(cmd, intern = TRUE)
    if (any(grepl("Error|error|failed|SSL", res))) {
      stop("SSL connect error.")
    }
    text_content <- paste(res, collapse = "\n")
    con <- textConnection(text_content)
    cd_result <- read.csv(con, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    close(con)
    if (ncol(cd_result) < 2) {
      stop("Invalid result format.")
    }
    message("[DEBUG] cd_search: parsed result (curl), nrow = ", nrow(cd_result))
    rv$cd_auto_result <- cd_result
    showNotification("CD-Search completed via curl!", type = "message")
    progress$set(value = 1, detail = "Done")
  }, error = function(e) {
    message("[ERROR] cd_search (curl): ", conditionMessage(e))
    showNotification(
      "Batch CD-Search failed. Please use the manual upload method:\n1. Visit https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi\n2. Upload the FASTA file, select 'CDD--62456 PSSMs'.\n3. Download the result TSV file.\n4. Upload it below.",
      type = "error", duration = 15
    )
  })
})

# ---------- 手动上传结果 ----------
cd_manual <- reactive({
  req(input$cd_manual_file)
  file_path <- input$cd_manual_file$datapath
  message("[DEBUG] cd_manual: reading file: ", file_path)
  tryCatch({
    line <- readLines(file_path, n = 1, warn = FALSE)
    sep <- if (grepl("\t", line)) "\t" else ","
    df <- read.csv(file_path, sep = sep, stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)
    if (ncol(df) < 2) {
      showNotification("The uploaded CD-Search result has too few columns.", type = "error")
      return(NULL)
    }
    message("[DEBUG] cd_manual: parsed, nrow = ", nrow(df), ", ncol = ", ncol(df))
    df
  }, error = function(e) {
    message("[ERROR] cd_manual: ", e$message)
    showNotification("Failed to read CD-Search result file.", type = "error")
    NULL
  })
})

# ---------- 最终使用的 CD 结果（自动优先，否则手动） ----------
cd_result <- reactive({
  if (!is.null(rv$cd_auto_result)) {
    message("[DEBUG] cd_result: using auto result")
    return(rv$cd_auto_result)
  }
  if (!is.null(cd_manual())) {
    message("[DEBUG] cd_result: using manual result")
    return(cd_manual())
  }
  NULL
})

# ---------- UI 输出 ----------
output$cd_loaded <- reactive({
  !is.null(cd_result())
})
outputOptions(output, "cd_loaded", suspendWhenHidden = FALSE)

output$cd_preview <- DT::renderDT({
  req(cd_result())
  DT::datatable(cd_result(), options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
})

output$cd_search_status_ui <- renderUI({
  if (!is.null(cd_result())) {
    div(style = "color: #5cb85c;", icon("check-circle"), " CD-Search result loaded (", nrow(cd_result()), " rows)")
  } else {
    div(style = "color: #999;", "No CD-Search result yet.")
  }
})

output$cd_status_ui <- renderUI({
  if (is.null(cd_result())) {
    return(div(style = "margin-top: 10px; color: #999;", "No CD-Search annotation loaded."))
  }
  div(style = "margin-top: 10px; color: #5cb85c;",
      icon("check-circle"), " CD-Search annotation loaded: ", nrow(cd_result()), " entries")
})

# ---------- 通用合并函数（用于蛋白表格） ----------
add_cd_to_table <- function(protein_table, id_col_name = "Protein") {
  cds <- cd_result()
  if (is.null(cds) || nrow(cds) == 0) {
    message("[DEBUG] add_cd_to_table: no CD result, returning original")
    return(protein_table)
  }
  if (!(id_col_name %in% colnames(protein_table))) {
    message("[DEBUG] add_cd_to_table: column '", id_col_name, "' not found in protein table")
    return(protein_table)
  }
  # CD 结果的第一列通常是查询 ID
  id_col_cd <- colnames(cds)[1]
  names(cds)[names(cds) == id_col_cd] <- id_col_name
  merged <- merge(protein_table, cds, by = id_col_name, all.x = TRUE, sort = FALSE)
  message("[DEBUG] add_cd_to_table: merged CD annotation, new cols: ", 
          paste(setdiff(colnames(merged), colnames(protein_table)), collapse = ", "))
  return(merged)
}