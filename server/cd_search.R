# server/cd_search.R
# ============================================
# 最终版：NCBI Batch CD-Search Web API（全程验证通过）
#   submit(tdata 触发文本) -> #cdsid -> poll(#status) -> fetch -> 解析(洗 Query ID)
#   批大小 950；当前工作目录下 cd_cache/<文件名>/ 存 chunk rds + progress.txt + 自动导出 TSV
#   断点续跑；跨平台 curl；保留一蛋白多行命中
# ============================================

message("[DEBUG] cd_search.R loading... (Final + auto TSV export)")

# ---------- 配置 ----------
CD_USE_PROXY <- FALSE                          # ← 测试直连设 FALSE；走代理设 TRUE
CD_PROXY     <- "socks5h://127.0.0.1:10808"
CD_PROXY_ARG <- if (CD_USE_PROXY) paste("-x", CD_PROXY) else ""
CD_ENDPOINT <- "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi"
CD_CHUNK    <- 950                           # 1000 会被 NCBI 判超限，必须 <1000
CD_CURL     <- if (.Platform$OS.type == "windows") "curl.exe" else "curl"

# 自动：当前工作目录下建 cd_cache/ 总缓存目录
CD_CACHE_DIR <- normalizePath(file.path(getwd(), "cd_cache"), mustWork = FALSE)
dir.create(CD_CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
message("[DEBUG] cd_search: cache root = ", CD_CACHE_DIR)

# 用上传文件名作为任务子目录名（去扩展名 + 清理非法字符）
make_task_id <- function(fname) {
  base <- tools::file_path_sans_ext(basename(fname))   # Ptenuiflora_genome_V1.prot.fasta -> Ptenuiflora_genome_V1.prot
  gsub('[\\\\/:*?"<>|]', "_", base)
}

# ---------- 读取 FASTA ----------
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
    ids  <- names(seqs)
    sequences <- unlist(seqs)
    message("[DEBUG] cd_search: loaded ", length(ids), " sequences")
    data.frame(ID = ids, Sequence = sequences, stringsAsFactors = FALSE)
  }, error = function(e) {
    message("[ERROR] cd_search: ", conditionMessage(e))
    showNotification(paste("Failed to read FASTA:", conditionMessage(e)), type = "error")
    NULL
  })
})

# ---------- 1) 提交，返回 cdsid（tdata 触发文本；超限硬失败，临时错误退避重试） ----------
cd_submit <- function(temp_fasta, db = "cdd", smode = "auto", retries = 3) {
  for (attempt in seq_len(retries)) {
    cmd <- sprintf(
      paste0('%s %s -s -X POST "%s" ',
             '--data-urlencode "queries@%s" ',
             '--data "db=%s" --data "smode=%s" ',
             '--data "tdata=hits" --data "dmode=rep" ',
             '--data "cddefl=false" --data "qdefl=false" --data "clonly=false" ',
             '--data "useid1=true" --max-time 120'),
      CD_CURL, CD_PROXY_ARG, CD_ENDPOINT, temp_fasta, db, smode
    )
    out <- tryCatch(system(cmd, intern = TRUE, ignore.stderr = FALSE),
                    error = function(e) character(0))
    txt <- paste(out, collapse = "\n")
    
    m <- regmatches(txt, regexpr("QM3-qcdsearch-[A-Za-z0-9-]+", txt))
    if (length(m) > 0) return(list(success = TRUE, cdsid = m[1]))
    
    if (grepl("too many|1000 or less", txt, ignore.case = TRUE)) {
      return(list(success = FALSE,
                  error = "数量超限（NCBI 拒绝），请把 CD_CHUNK 降到更小（如 900）"))
    }
    transient <- nchar(txt) == 0 ||
      grepl("temporarily|try again|service unavailable|busy", txt, ignore.case = TRUE)
    if (transient && attempt < retries) { Sys.sleep(30 * attempt); next }
    return(list(success = FALSE, error = paste("提交被拒，预览:", substr(txt, 1, 300))))
  }
  list(success = FALSE, error = "重试多次仍失败")
}

# ---------- 2) 轮询 + 取结果（用 #status 判断） ----------
cd_poll_fetch <- function(cdsid, max_wait = 900, interval = 20) {
  elapsed <- 0
  repeat {
    cmd <- sprintf(
      paste0('%s %s -s -X POST "%s" --data "cdsid=%s" --data "tdata=hits" ',
             '--data "cddefl=true" --data "qdefl=true" --data "dmode=rep" ',
             '--data "clonly=false" --max-time 120'),
      CD_CURL, CD_PROXY_ARG, CD_ENDPOINT, cdsid
    )
    out <- tryCatch(system(cmd, intern = TRUE, ignore.stderr = FALSE),
                    error = function(e) character(0))
    if (length(out) == 0) {
      Sys.sleep(interval); elapsed <- elapsed + interval
      if (elapsed >= max_wait) return(list(success = FALSE, error = "轮询无响应超时"))
      next
    }
    
    status_line <- grep("^#status", out, value = TRUE)
    status <- if (length(status_line)) {
      as.integer(sub(".*#status[^0-9]*([0-9]+).*", "\\1", status_line[1]))
    } else NA_integer_
    
    if (!is.na(status) && status == 0) {
      return(list(success = TRUE, content = paste(out, collapse = "\n")))
    }
    if (!is.na(status) && status == 3) {
      message("[DEBUG] running... (", elapsed, "s) cdsid=", cdsid)
      Sys.sleep(interval); elapsed <- elapsed + interval
      if (elapsed >= max_wait) return(list(success = FALSE, error = "等待结果超时"))
      next
    }
    if (is.na(status)) {
      Sys.sleep(interval); elapsed <- elapsed + interval
      if (elapsed >= max_wait) return(list(success = FALSE, error = "未读到 #status"))
      next
    }
    return(list(success = FALSE, error = paste("NCBI status =", status)))
  }
}

# ---------- TSV 文本 -> data.frame（去 # 注释行；洗 Query 列 ID；空命中返回 0 行） ----------
cd_parse_tsv <- function(content) {
  lines <- strsplit(content, "\n")[[1]]
  data_lines <- lines[!grepl("^#", lines) & nchar(trimws(lines)) > 0]
  if (length(data_lines) < 1) return(NULL)
  df <- tryCatch(
    read.csv(text = paste(data_lines, collapse = "\n"),
             sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(df) || ncol(df) < 2) return(NULL)
  qcol <- names(df)[1]
  df[[qcol]] <- trimws(sub("^Q#\\d+\\s*-\\s*>?", "", df[[qcol]]))   # Q#1 - >Pt_Chr01 -> Pt_Chr01
  df
}

# ============================================
# 主按钮：自动分批 + 断点续跑 + 进度文件 + 自动导出 TSV
# ============================================
observeEvent(input$run_cd_search, {
  message("[DEBUG] cd_search: button clicked")
  fasta_df <- cd_fasta()
  if (is.null(fasta_df) || nrow(fasta_df) == 0) {
    showNotification("Please upload a FASTA file first.", type = "error")
    return()
  }
  
  n_total <- nrow(fasta_df)
  chunks  <- split(fasta_df, ceiling(seq_len(n_total) / CD_CHUNK))
  n_chunk <- length(chunks)
  
  # 每个任务独立子目录（按上传文件名）
  task_id <- make_task_id(input$cd_search_fasta$name)
  out_dir <- file.path(CD_CACHE_DIR, task_id)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  message("[DEBUG] cd_search: task dir = ", out_dir)
  
  # PowerShell 可实时跟踪的纯文本进度文件
  progress_file <- file.path(out_dir, "progress.txt")
  log_progress <- function(msg) {
    cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg),
        file = progress_file, append = TRUE)
  }
  log_progress(sprintf("==== START: %d 序列 / %d 批 (chunk=%d) ====",
                       n_total, n_chunk, CD_CHUNK))
  
  progress <- shiny::Progress$new()
  progress$set(message = sprintf("CD-Search: %d 序列 / %d 批", n_total, n_chunk), value = 0)
  on.exit(progress$close())
  
  collected <- vector("list", n_chunk)
  failed    <- integer(0)
  
  for (i in seq_len(n_chunk)) {
    progress$set(value = (i - 1) / n_chunk, detail = sprintf("Chunk %d/%d ...", i, n_chunk))
    
    cache_file <- file.path(out_dir, sprintf("chunk_%03d.rds", i))
    if (file.exists(cache_file)) {
      collected[[i]] <- readRDS(cache_file)
      message("[DEBUG] chunk ", i, " 命中缓存，跳过")
      log_progress(sprintf("chunk %d/%d  SKIP (已缓存)", i, n_chunk))
      next
    }
    
    ck <- chunks[[i]]
    tf <- tempfile(fileext = ".fasta")
    writeLines(paste0(">", ck$ID, "\n", ck$Sequence), tf)
    
    sub <- cd_submit(tf)
    if (file.exists(tf)) file.remove(tf)
    if (!sub$success) {
      message("[ERROR] chunk ", i, " 提交失败: ", sub$error)
      log_progress(sprintf("chunk %d/%d  FAIL (submit): %s", i, n_chunk, sub$error))
      failed <- c(failed, i); next
    }
    message("[DEBUG] chunk ", i, " cdsid = ", sub$cdsid)
    log_progress(sprintf("chunk %d/%d  submitted cdsid=%s", i, n_chunk, sub$cdsid))
    
    res <- cd_poll_fetch(sub$cdsid)
    if (!res$success) {
      message("[ERROR] chunk ", i, " 取结果失败: ", res$error)
      log_progress(sprintf("chunk %d/%d  FAIL (fetch): %s", i, n_chunk, res$error))
      failed <- c(failed, i); next
    }
    
    df <- cd_parse_tsv(res$content)
    if (is.null(df)) {
      message("[ERROR] chunk ", i, " 解析为空")
      log_progress(sprintf("chunk %d/%d  FAIL (parse)", i, n_chunk))
      failed <- c(failed, i); next
    }
    
    saveRDS(df, cache_file)
    collected[[i]] <- df
    log_progress(sprintf("chunk %d/%d  DONE  rows=%d", i, n_chunk, nrow(df)))
    Sys.sleep(10)
  }
  
  collected <- collected[!vapply(collected, is.null, logical(1))]
  if (length(collected) == 0) {
    log_progress("==== END: 全部失败 ====")
    showNotification("❌ 所有批次都失败了，请检查 proxy / captive portal。",
                     type = "error", duration = 15)
    return()
  }
  
  cd_all <- tryCatch(do.call(rbind, collected),
                     error = function(e) {
                       cols <- Reduce(intersect, lapply(collected, names))
                       do.call(rbind, lapply(collected, function(d) d[, cols, drop = FALSE]))
                     })
  
  rv$cd_auto_result <- cd_all
  progress$set(value = 1, detail = "Done!")
  
  # 跑完自动导出合并 TSV 到任务目录
  tsv_path <- file.path(out_dir, "cd_search_all.tsv")
  tryCatch({
    write.table(cd_all, tsv_path, sep = "\t", row.names = FALSE,
                quote = FALSE, fileEncoding = "UTF-8")
    message("[DEBUG] 已导出 TSV: ", tsv_path)
    log_progress(sprintf("==== END: %d 行，成功 %d/%d 批；TSV -> %s ====",
                         nrow(cd_all), length(collected), n_chunk, tsv_path))
  }, error = function(e) {
    message("[ERROR] TSV 导出失败: ", conditionMessage(e))
    log_progress(sprintf("==== END: %d 行，成功 %d/%d 批；TSV 导出失败: %s ====",
                         nrow(cd_all), length(collected), n_chunk, conditionMessage(e)))
  })
  
  msg <- sprintf("✅ CD-Search 完成：%d 行，成功 %d/%d 批；已导出 cd_search_all.tsv",
                 nrow(cd_all), length(collected), n_chunk)
  if (length(failed) > 0) {
    msg <- paste0(msg, "；失败批次: ", paste(failed, collapse = ", "),
                  "（缓存已保留，再次点击可续跑）")
  }
  showNotification(msg, type = if (length(failed) > 0) "warning" else "message", duration = 12)
})

# ---------- 手动上传（保持不变） ----------
cd_manual <- reactive({
  req(input$cd_manual_file)
  file_path <- input$cd_manual_file$datapath
  message("[DEBUG] cd_manual: reading file: ", file_path)
  tryCatch({
    line <- readLines(file_path, n = 1, warn = FALSE)
    sep  <- if (grepl("\t", line)) "\t" else ","
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

output$cd_loaded <- reactive({ !is.null(cd_result()) })
outputOptions(output, "cd_loaded", suspendWhenHidden = FALSE)

output$cd_preview <- DT::renderDT({
  req(cd_result())
  DT::datatable(cd_result(), options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
})

output$cd_search_status_ui <- renderUI({
  if (!is.null(cd_result())) {
    div(style = "color: #5cb85c;", icon("check-circle"),
        " CD-Search result loaded (", nrow(cd_result()), " rows)")
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

add_cd_to_table <- function(protein_table, id_col_name = "Protein") {
  cds <- cd_result()
  if (is.null(cds) || nrow(cds) == 0) return(protein_table)
  if (!(id_col_name %in% colnames(protein_table))) return(protein_table)
  id_col_cd <- colnames(cds)[1]
  names(cds)[names(cds) == id_col_cd] <- id_col_name
  merge(protein_table, cds, by = id_col_name, all.x = TRUE, sort = FALSE)
}

output$download_cd_fasta <- downloadHandler(
  filename = function() paste0("cd_search_sequences_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".fasta"),
  content = function(file) {
    fasta_df <- cd_fasta()
    if (is.null(fasta_df) || nrow(fasta_df) == 0) {
      writeLines(">No_sequences\n", file)
      return()
    }
    fasta_lines <- paste0(">", fasta_df$ID, "\n", fasta_df$Sequence)
    writeLines(fasta_lines, file)
  }
)

message("[DEBUG] cd_search.R loaded successfully (Final + auto TSV export)")