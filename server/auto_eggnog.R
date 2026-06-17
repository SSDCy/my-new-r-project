# server/auto_eggnog.R
message("[DEBUG] auto_eggnog.R loading... (httr first, then curl, manual fallback)")

# ---------- 读取参考 FASTA ----------
reference_sequences <- reactive({
  req(input$reference_fasta)
  file_path <- input$reference_fasta$datapath
  message("[DEBUG] auto_eggnog: reference_fasta uploaded: ", file_path)
  if (!requireNamespace("seqinr", quietly = TRUE)) {
    showNotification("Package 'seqinr' is required. Please install it.", type = "error")
    return(NULL)
  }
  tryCatch({
    seqs <- seqinr::read.fasta(file_path, seqtype = "AA", as.string = TRUE)
    ids <- names(seqs)
    sequences <- unlist(seqs)
    message("[DEBUG] auto_eggnog: loaded ", length(ids), " sequences")
    if (length(ids) == 0) {
      showNotification("The uploaded FASTA file appears to be empty.", type = "error")
      return(NULL)
    }
    data.frame(ID = ids, Sequence = sequences, stringsAsFactors = FALSE)
  }, error = function(e) {
    message("[ERROR] auto_eggnog: failed to read FASTA: ", conditionMessage(e))
    showNotification(paste("Failed to read reference FASTA:", conditionMessage(e)), type = "error")
    NULL
  })
})

# ---------- 辅助函数：使用 httr 提交 eggNOG 作业 ----------
eggnog_post_httr <- function(fasta_path) {
  res <- httr::POST(
    url = "https://eggnog-mapper.embl.de/api/job",
    config = httr::config(ssl_verifyhost = FALSE, ssl_verifypeer = FALSE),
    body = list(
      data = httr::upload_file(fasta_path, type = "text/plain"),
      tax_scope = "auto",
      target_orthologs = "all",
      go_evidence = "non-electronic",
      pfam_realign = "false"
    ),
    encode = "multipart"
  )
  httr::stop_for_status(res)
  httr::content(res, as = "parsed")
}

# ---------- 辅助函数：使用 httr 获取作业状态 ----------
eggnog_get_status_httr <- function(job_id) {
  res <- httr::GET(
    paste0("https://eggnog-mapper.embl.de/api/job/", job_id),
    config = httr::config(ssl_verifyhost = FALSE, ssl_verifypeer = FALSE)
  )
  httr::stop_for_status(res)
  httr::content(res, as = "parsed")
}

# ---------- 辅助函数：使用 httr 下载结果 ----------
eggnog_download_httr <- function(job_id) {
  res <- httr::GET(
    paste0("https://eggnog-mapper.embl.de/api/job/", job_id, "/results/annotations"),
    config = httr::config(ssl_verifyhost = FALSE, ssl_verifypeer = FALSE)
  )
  httr::stop_for_status(res)
  httr::content(res, as = "text", encoding = "UTF-8")
}

# ---------- 自动注释按钮 ----------
observeEvent(input$auto_eggnog, {
  message("[DEBUG] auto_eggnog: button clicked")
  
  if (is.null(input$reference_fasta)) {
    showNotification("Please upload a reference proteome FASTA file first.", type = "error")
    return()
  }
  ref <- reference_sequences()
  if (is.null(ref) || nrow(ref) == 0) {
    showNotification("Failed to read sequences from the uploaded FASTA file.", type = "error")
    return()
  }
  if (is.null(rv$clean_data)) {
    showNotification("No expression data loaded. Please upload MaxQuant data first.", type = "error")
    return()
  }
  protein_ids <- rv$clean_data$`Master protein IDs`
  protein_ids <- unique(protein_ids[!is.na(protein_ids) & protein_ids != ""])
  if (length(protein_ids) == 0) {
    showNotification("No protein IDs found in data.", type = "error")
    return()
  }
  
  matched <- ref[ref$ID %in% protein_ids, ]
  if (nrow(matched) == 0) {
    showNotification("None of the protein IDs could be matched in the reference FASTA. Check ID consistency.", type = "error")
    return()
  }
  
  fasta_lines <- paste0(">", matched$ID, "\n", matched$Sequence)
  fasta_text <- paste(fasta_lines, collapse = "\n")
  temp_fasta <- tempfile(fileext = ".fasta")
  writeLines(fasta_text, temp_fasta)
  message("[DEBUG] auto_eggnog: wrote FASTA to temp file: ", temp_fasta)
  
  progress <- shiny::Progress$new()
  progress$set(message = "Submitting to eggNOG...", value = 0)
  on.exit({
    if (file.exists(temp_fasta)) file.remove(temp_fasta)
    progress$close()
  })
  
  # ---- 方法一：使用 httr ----
  tryCatch({
    message("[DEBUG] auto_eggnog: trying httr POST")
    job <- eggnog_post_httr(temp_fasta)
    job_id <- job$id
    message("[DEBUG] auto_eggnog: job submitted (httr), id = ", job_id)
    progress$set(detail = paste("Job", job_id, "submitted"))
    
    # 轮询状态
    status <- ""
    while (!status %in% c("done", "failed")) {
      Sys.sleep(10)
      job2 <- eggnog_get_status_httr(job_id)
      status <- job2$status
      message("[DEBUG] auto_eggnog: job status = ", status)
      progress$set(value = 0.3, detail = paste("Status:", status))
    }
    if (status == "failed") {
      showNotification("eggNOG annotation job failed.", type = "error")
      return()
    }
    
    # 下载结果
    text_content <- eggnog_download_httr(job_id)
    message("[DEBUG] auto_eggnog: downloaded annotation, length = ", nchar(text_content), " chars")
    
    con <- textConnection(text_content)
    eggnog_annot <- read.csv(con, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    close(con)
    message("[DEBUG] auto_eggnog: parsed annotation, nrow = ", nrow(eggnog_annot))
    
    id_col <- colnames(eggnog_annot)[1]
    names(eggnog_annot)[names(eggnog_annot) == id_col] <- "Master protein IDs"
    eggnog_annot <- eggnog_annot[!duplicated(eggnog_annot$`Master protein IDs`), ]
    rv$auto_eggnog_result <- eggnog_annot
    showNotification("eggNOG annotation completed successfully!", type = "message", duration = 10)
    progress$set(value = 1, detail = "Done")
    return()
  }, error = function(e) {
    message("[ERROR] auto_eggnog (httr): ", conditionMessage(e))
  })
  
  # ---- 方法二：回退到系统 curl ----
  message("[DEBUG] auto_eggnog: trying system curl")
  system_curl <- tryCatch({
    system("curl --version", ignore.stdout = TRUE, ignore.stderr = TRUE)
    TRUE
  }, error = function(e) FALSE)
  
  if (!system_curl) {
    showNotification(
      "Auto annotation unavailable. Please use manual upload:\n1. Download protein FASTA from 'Download Protein Sequences' section.\n2. Submit it to http://eggnog-mapper.embl.de\n3. Download the result file.\n4. Upload it using 'Upload eggNOG Annotation (Manual)'.",
      type = "error", duration = 15
    )
    return()
  }
  
  tryCatch({
    cmd <- sprintf('curl -k -s -X POST "https://eggnog-mapper.embl.de/api/job" -F "data=@%s;type=text/plain" -F "tax_scope=auto" -F "target_orthologs=all" -F "go_evidence=non-electronic" -F "pfam_realign=false"',
                   shQuote(temp_fasta))
    message("[DEBUG] auto_eggnog: running curl command")
    res <- system(cmd, intern = TRUE)
    if (any(grepl("SEC_E_ILLEGAL_MESSAGE|SSL connect error|failed", res))) {
      stop("SSL connect error.")
    }
    job <- jsonlite::fromJSON(paste(res, collapse = "\n"))
    job_id <- job$id
    message("[DEBUG] auto_eggnog: job submitted (curl), id = ", job_id)
    progress$set(detail = paste("Job", job_id, "submitted"))
    
    # 轮询状态
    status <- ""
    while (!status %in% c("done", "failed")) {
      Sys.sleep(10)
      res2 <- system(sprintf('curl -k -s "https://eggnog-mapper.embl.de/api/job/%s"', job_id), intern = TRUE)
      job2 <- jsonlite::fromJSON(paste(res2, collapse = "\n"))
      status <- job2$status
      message("[DEBUG] auto_eggnog: job status = ", status)
      progress$set(value = 0.3, detail = paste("Status:", status))
    }
    if (status == "failed") {
      showNotification("eggNOG annotation job failed.", type = "error")
      return()
    }
    
    # 下载结果
    res3 <- system(sprintf('curl -k -s "https://eggnog-mapper.embl.de/api/job/%s/results/annotations"', job_id), intern = TRUE)
    text_content <- paste(res3, collapse = "\n")
    con <- textConnection(text_content)
    eggnog_annot <- read.csv(con, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    close(con)
    message("[DEBUG] auto_eggnog: parsed annotation (curl), nrow = ", nrow(eggnog_annot))
    
    id_col <- colnames(eggnog_annot)[1]
    names(eggnog_annot)[names(eggnog_annot) == id_col] <- "Master protein IDs"
    eggnog_annot <- eggnog_annot[!duplicated(eggnog_annot$`Master protein IDs`), ]
    rv$auto_eggnog_result <- eggnog_annot
    showNotification("eggNOG annotation completed via curl!", type = "message", duration = 10)
    progress$set(value = 1, detail = "Done")
  }, error = function(e) {
    message("[ERROR] auto_eggnog (curl): ", conditionMessage(e))
    showNotification(
      "Auto annotation failed due to network/SSL. Please use the manual upload method:\n1. Download protein FASTA from 'Download Protein Sequences' section.\n2. Submit it to http://eggnog-mapper.embl.de\n3. Download the result file.\n4. Upload it using 'Upload eggNOG Annotation (Manual)'.",
      type = "error", duration = 15
    )
  })
})

# ---------- UI 状态显示 ----------
output$auto_eggnog_status_ui <- renderUI({
  if (!is.null(rv$auto_eggnog_result)) {
    div(style = "color: #5cb85c;", icon("check-circle"), " Auto annotation ready (",
        nrow(rv$auto_eggnog_result), " entries)")
  } else {
    div(style = "color: #999;", "Auto annotation may not work. Use manual upload.")
  }
})

output$reference_fasta_hint <- renderUI({
  if (is.null(input$reference_fasta)) {
    div(style = "color: #d9534f; font-weight: bold;", 
        icon("exclamation-triangle"), " Please upload the reference FASTA file to enable auto annotation.")
  } else {
    div(style = "color: #5cb85c;", 
        icon("check-circle"), " FASTA file uploaded (auto may fail due to SSL).")
  }
})