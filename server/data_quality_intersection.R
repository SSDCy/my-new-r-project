# server/data_quality_intersection.R
message("[DEBUG] data_quality_intersection.R loading... (barplot style, square shape, no cut-off)")

# ---------- 更新可选样本列表 ----------
observe({
  req(rv$sample_names)
  updateSelectizeInput(session, "intersection_samples", choices = rv$sample_names, server = TRUE)
  message("[DEBUG] intersection_samples choices updated: n = ", length(rv$sample_names))
})

# ---------- 监控用户实际选择的样本 ----------
observe({
  sel <- input$intersection_samples
  message("[DEBUG] intersection_samples current selection: ", 
          if (is.null(sel)) "NULL" else paste0("(", length(sel), ") ", paste(sel, collapse = ", ")))
})

# ---------- 蛋白存在矩阵生成函数 ----------
get_presence_matrix <- function(sample_names) {
  mat <- dq_expr_matrix()
  if (is.null(mat) || nrow(mat) == 0) {
    message("[DEBUG] get_presence_matrix: dq_expr_matrix is empty")
    return(NULL)
  }
  message("[DEBUG] get_presence_matrix: dq_expr_matrix columns = ", paste(head(colnames(mat), 10), collapse = ", "), "...")
  
  available_samples <- intersect(sample_names, colnames(mat))
  if (length(available_samples) == 0) {
    message("[DEBUG] get_presence_matrix: no selected samples found in expression matrix")
    return(NULL)
  }
  message("[DEBUG] get_presence_matrix: selected samples found: ", paste(available_samples, collapse = ", "))
  
  sub_mat <- mat[, available_samples, drop = FALSE]
  presence <- ifelse(!is.na(sub_mat) & sub_mat > 0, 1, 0)
  presence <- as.data.frame(presence)
  presence$Protein <- rownames(presence)
  presence$Sum <- rowSums(presence[, available_samples, drop = FALSE])
  
  message("[DEBUG] get_presence_matrix: presence matrix created, nrow = ", nrow(presence), 
          ", ncol = ", ncol(presence) - 2)
  message("[DEBUG] get_presence_matrix: first 3 rows:")
  print(head(presence, 3))
  return(presence)
}

# ---------- 基于样本选择计算数据 ----------
intersection_data <- eventReactive(input$generate_intersection, {
  samples <- input$intersection_samples
  message("[DEBUG] intersection_data: GENERATE button clicked. Current selection: ", 
          if (is.null(samples)) "NULL" else paste0("(", length(samples), ") ", paste(samples, collapse = ", ")))
  
  if (length(samples) < 2 || length(samples) > 15) {
    showNotification("Please select 2–15 samples.", type = "error")
    message("[DEBUG] intersection_data: INVALID selection count (", length(samples), ")")
    return(NULL)
  }
  
  presence_matrix <- get_presence_matrix(samples)
  if (is.null(presence_matrix)) {
    showNotification("No data available for selected samples.", type = "error")
    return(NULL)
  }
  
  sample_cols <- setdiff(colnames(presence_matrix), c("Protein", "Sum"))
  n_total <- nrow(presence_matrix)
  shared_all <- sum(presence_matrix$Sum == length(sample_cols))
  unique_counts <- sapply(sample_cols, function(s) {
    sum(presence_matrix[[s]] == 1 & presence_matrix$Sum == 1)
  })
  unique_total <- sum(unique_counts)
  
  detected_proteins <- presence_matrix$Protein[presence_matrix$Sum > 0]
  message("[DEBUG] intersection_data: detected proteins with Sum>0 = ", length(detected_proteins))
  
  message("[DEBUG] intersection_data: total proteins = ", n_total, 
          ", shared_all = ", shared_all, 
          ", unique_per_sample = ", paste(unique_counts, collapse = ", "))
  
  list(
    presence_matrix = presence_matrix,
    samples = sample_cols,
    total_proteins = n_total,
    shared_all = shared_all,
    unique_counts = unique_counts,
    unique_total = unique_total,
    detected_proteins = detected_proteins
  )
})

# ---------- 摘要 ----------
output$intersection_summary <- renderText({
  data <- intersection_data()
  if (is.null(data)) {
    msg <- "No data to display. Please select 2-15 samples and click 'Generate Table'."
    return(msg)
  }
  
  lines <- c()
  lines <- c(lines, paste("Selected samples:", paste(data$samples, collapse = ", ")))
  lines <- c(lines, paste("Total proteins evaluated:", data$total_proteins))
  lines <- c(lines, paste("Proteins detected in all selected samples (shared):", data$shared_all))
  lines <- c(lines, paste("Proteins detected in at least one sample:", length(data$detected_proteins)))
  lines <- c(lines, "")
  lines <- c(lines, "Unique proteins per sample:")
  for (s in names(data$unique_counts)) {
    lines <- c(lines, sprintf("  %s: %d", s, data$unique_counts[s]))
  }
  lines <- c(lines, "")
  lines <- c(lines, paste("Total unique proteins (detected in exactly one sample):", data$unique_total))
  
  paste(lines, collapse = "\n")
})

# ---------- 蛋白存在矩阵表格 ----------
output$intersection_protein_table <- DT::renderDataTable({
  data <- intersection_data()
  if (is.null(data)) {
    message("[DEBUG] intersection_protein_table: data is NULL, placeholder")
    return(DT::datatable(data.frame(Message = "No data")))
  }
  
  df <- data$presence_matrix
  message("[DEBUG] intersection_protein_table: rendering datatable, nrow = ", nrow(df),
          ", ncol = ", ncol(df))
  
  if (nrow(df) == 0) {
    return(DT::datatable(data.frame(Message = "No rows")))
  }
  
  dt <- DT::datatable(
    df,
    options = list(pageLength = 10),
    rownames = FALSE
  )
  
  if ("Sum" %in% colnames(df)) {
    dt <- DT::formatStyle(dt, columns = "Sum", fontWeight = "bold", color = "#2c3e50")
  }
  for (s in data$samples) {
    if (s %in% colnames(df)) {
      dt <- DT::formatStyle(dt, columns = s,
                            backgroundColor = DT::styleEqual(c(0, 1), c("#f2f2f2", "#c6efce")))
    }
  }
  dt
})

# ---------- 下载蛋白表 CSV ----------
output$download_intersection_proteins <- downloadHandler(
  filename = function() paste0("protein_presence_matrix_", Sys.Date(), ".csv"),
  content = function(file) {
    data <- intersection_data()
    if (is.null(data)) return(write.csv(data.frame(Error = "No data"), file))
    write.csv(data$presence_matrix, file, row.names = FALSE)
    message("[DEBUG] download_intersection_proteins: saved")
  }
)

# ---------- 肽段表格 ----------
output$intersection_peptide_table <- DT::renderDataTable({
  data <- intersection_data()
  req(data)
  
  mode <- input$peptide_display_mode
  if (is.null(mode)) mode <- "merged"
  message("[DEBUG] intersection_peptide_table: display mode = ", mode)
  
  target_proteins <- data$detected_proteins
  message("[DEBUG] intersection_peptide_table: using ", length(target_proteins), " detected proteins")
  
  if (length(target_proteins) == 0) {
    return(DT::datatable(data.frame(Message = "No proteins detected in selected samples.")))
  }
  
  clean <- rv$clean_data
  if (is.null(clean) || !"Master protein IDs" %in% colnames(clean)) {
    df <- data.frame(Protein = target_proteins, Peptide_Sequence = "Not available")
    return(DT::datatable(df, options = list(pageLength = 15), rownames = FALSE))
  }
  
  idx <- match(target_proteins, clean$`Master protein IDs`)
  if (all(is.na(idx))) {
    df <- data.frame(Protein = target_proteins, Peptide_Sequence = "No match")
    return(DT::datatable(df, options = list(pageLength = 15), rownames = FALSE))
  }
  
  if (!"Peptide sequences" %in% colnames(clean)) {
    sub <- clean[idx, , drop = FALSE]
    display_cols <- c("Master protein IDs", "Unique peptides", "Protein IDs")
    available <- intersect(display_cols, colnames(sub))
    if (length(available) == 0) available <- "Master protein IDs"
    sub <- sub[, available, drop = FALSE]
    message("[DEBUG] intersection_peptide_table: no Peptide sequences column")
    return(DT::datatable(sub, options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE))
  }
  
  peptide_seq <- clean$`Peptide sequences`[idx]
  names(peptide_seq) <- clean$`Master protein IDs`[idx]
  
  presence <- data$presence_matrix
  rownames(presence) <- presence$Protein
  sample_cols <- data$samples
  
  get_detected_samples <- function(protein_id) {
    if (!protein_id %in% rownames(presence)) return(NA_character_)
    row <- presence[protein_id, sample_cols, drop = FALSE]
    detected <- sample_cols[row == 1]
    if (length(detected) == 0) return("None")
    paste(detected, collapse = ", ")
  }
  
  if (mode == "merged") {
    merged_list <- lapply(names(peptide_seq), function(pid) {
      seqs <- trimws(unlist(strsplit(peptide_seq[pid], ";")))
      seqs <- seqs[seqs != ""]
      if (length(seqs) == 0) return(NULL)
      data.frame(`Master protein IDs` = pid, 
                 `Peptide Sequences` = paste(seqs, collapse = "; "),
                 stringsAsFactors = FALSE, check.names = FALSE)
    })
    df_merged <- do.call(rbind, merged_list)
    if (is.null(df_merged) || nrow(df_merged) == 0) {
      return(DT::datatable(data.frame(Message = "No peptide sequences found.")))
    }
    message("[DEBUG] intersection_peptide_table: merged mode, nrow = ", nrow(df_merged))
    message("[DEBUG] merged table first rows:")
    print(head(df_merged, 3))
    dt <- DT::datatable(df_merged, 
                        options = list(pageLength = 25, scrollX = TRUE),
                        rownames = FALSE,
                        selection = 'single')
    message("[DEBUG] intersection_peptide_table: merged mode with single selection configured")
    return(dt)
  } else {
    seq_list <- lapply(seq_along(peptide_seq), function(i) {
      pid <- names(peptide_seq)[i]
      seqs <- trimws(unlist(strsplit(peptide_seq[i], ";")))
      seqs <- seqs[seqs != ""]
      if (length(seqs) == 0) return(NULL)
      detected_samples <- get_detected_samples(pid)
      data.frame(`Master protein IDs` = pid, 
                 `Peptide Sequence` = seqs,
                 `Detected in Samples` = detected_samples,
                 stringsAsFactors = FALSE, check.names = FALSE)
    })
    df_seq <- do.call(rbind, seq_list)
    if (is.null(df_seq) || nrow(df_seq) == 0) {
      return(DT::datatable(data.frame(Message = "No peptide sequences found.")))
    }
    message("[DEBUG] intersection_peptide_table: expanded mode, nrow = ", nrow(df_seq))
    return(DT::datatable(df_seq, options = list(pageLength = 25, scrollX = TRUE), rownames = FALSE))
  }
})

# ---------- 下载肽段 CSV ----------
output$download_intersection_peptides <- downloadHandler(
  filename = function() "shared_unique_peptides.csv",
  content = function(file) {
    data <- intersection_data()
    req(data)
    mode <- input$peptide_display_mode
    target_proteins <- data$detected_proteins
    clean <- rv$clean_data
    if (is.null(clean) || !"Master protein IDs" %in% colnames(clean)) {
      return(write.csv(data.frame(Error = "No clean data"), file))
    }
    idx <- match(target_proteins, clean$`Master protein IDs`)
    sub <- clean[idx, , drop = FALSE]
    if ("Peptide sequences" %in% colnames(sub)) {
      peptide_seq <- sub$`Peptide sequences`
      names(peptide_seq) <- sub$`Master protein IDs`
      
      if (mode == "merged") {
        merged_list <- lapply(names(peptide_seq), function(pid) {
          seqs <- trimws(unlist(strsplit(peptide_seq[pid], ";")))
          seqs <- seqs[seqs != ""]
          if (length(seqs) == 0) return(NULL)
          data.frame(`Master protein IDs` = pid, 
                     `Peptide Sequences` = paste(seqs, collapse = "; "),
                     stringsAsFactors = FALSE, check.names = FALSE)
        })
        df_out <- do.call(rbind, merged_list)
      } else {
        presence <- data$presence_matrix
        rownames(presence) <- presence$Protein
        sample_cols <- data$samples
        get_detected_samples <- function(protein_id) {
          if (!protein_id %in% rownames(presence)) return(NA_character_)
          row <- presence[protein_id, sample_cols, drop = FALSE]
          detected <- sample_cols[row == 1]
          if (length(detected) == 0) return("None")
          paste(detected, collapse = ", ")
        }
        seq_list <- lapply(seq_along(peptide_seq), function(i) {
          pid <- names(peptide_seq)[i]
          seqs <- trimws(unlist(strsplit(peptide_seq[i], ";")))
          seqs <- seqs[seqs != ""]
          if (length(seqs) == 0) return(NULL)
          data.frame(`Master protein IDs` = pid, 
                     `Peptide Sequence` = seqs,
                     `Detected in Samples` = get_detected_samples(pid),
                     stringsAsFactors = FALSE, check.names = FALSE)
        })
        df_out <- do.call(rbind, seq_list)
      }
      if (!is.null(df_out)) write.csv(df_out, file, row.names = FALSE)
    } else {
      write.csv(sub, file, row.names = FALSE)
    }
    message("[DEBUG] download_intersection_peptides: saved (mode = ", mode, ")")
  }
)

# ========== 选中蛋白驱动直方图 ==========
selected_protein_for_hist <- reactive({
  if (input$peptide_display_mode != "merged") {
    return(NULL)
  }
  sel <- input$intersection_peptide_table_rows_selected
  message("[DEBUG] selected_protein_for_hist: raw selection input = ", 
          if (is.null(sel)) "NULL" else paste(sel, collapse = ", "))
  if (is.null(sel) || length(sel) == 0) {
    message("[DEBUG] selected_protein_for_hist: no row selected")
    return(NULL)
  }
  
  data <- intersection_data()
  if (is.null(data)) return(NULL)
  clean <- rv$clean_data
  if (is.null(clean) || !"Master protein IDs" %in% colnames(clean)) return(NULL)
  
  target_proteins <- data$detected_proteins
  idx <- match(target_proteins, clean$`Master protein IDs`)
  if (all(is.na(idx))) return(NULL)
  peptide_seq <- clean$`Peptide sequences`[idx]
  names(peptide_seq) <- clean$`Master protein IDs`[idx]
  
  merged_list <- lapply(names(peptide_seq), function(pid) {
    seqs <- trimws(unlist(strsplit(peptide_seq[pid], ";")))
    seqs <- seqs[seqs != ""]
    if (length(seqs) == 0) return(NULL)
    data.frame(ID = pid, stringsAsFactors = FALSE)
  })
  df_merged <- do.call(rbind, merged_list)
  if (is.null(df_merged) || nrow(df_merged) == 0) return(NULL)
  
  message("[DEBUG] selected_protein_for_hist: merged df nrow = ", nrow(df_merged), 
          ", selected row index = ", sel[1])
  message("[DEBUG] merged df first column name (colnames): ", colnames(df_merged)[1])
  message("[DEBUG] merged df first few IDs: ", paste(head(df_merged[[1]], 5), collapse = ", "))
  
  selected_row <- sel[1]
  if (selected_row < 1 || selected_row > nrow(df_merged)) return(NULL)
  
  pid <- df_merged[selected_row, 1]
  message("[DEBUG] selected_protein_for_hist: extracted protein ID = ", pid)
  return(as.character(pid))
})

# ========== 肽段长度数据 ==========
peptide_lengths <- reactive({
  data <- intersection_data()
  if (is.null(data)) return(NULL)
  clean <- rv$clean_data
  if (is.null(clean) || !"Peptide sequences" %in% colnames(clean)) return(NULL)
  
  target_proteins <- data$detected_proteins
  selected <- selected_protein_for_hist()
  
  if (!is.null(selected)) {
    if (!selected %in% target_proteins) {
      message("[DEBUG] peptide_lengths: selected protein not in detected proteins, fallback to all")
      selected <- NULL
    }
  }
  
  if (!is.null(selected)) {
    idx <- match(selected, clean$`Master protein IDs`)
    if (is.na(idx)) {
      message("[DEBUG] peptide_lengths: protein ID not found in clean data")
      return(numeric(0))
    }
    pep_seq <- clean$`Peptide sequences`[idx]
    seqs <- trimws(unlist(strsplit(pep_seq, ";")))
    seqs <- seqs[seqs != ""]
    lengths <- nchar(seqs)
    message("[DEBUG] peptide_lengths: using selected protein ", selected, 
            " with ", length(lengths), " peptides")
    message("[DEBUG] peptide_lengths: lengths = ", paste(lengths, collapse = ", "))
  } else {
    idx <- match(target_proteins, clean$`Master protein IDs`)
    peptide_seq <- clean$`Peptide sequences`[idx]
    names(peptide_seq) <- clean$`Master protein IDs`[idx]
    lengths <- c()
    for (pid in names(peptide_seq)) {
      seqs <- trimws(unlist(strsplit(peptide_seq[pid], ";")))
      seqs <- seqs[seqs != ""]
      lengths <- c(lengths, nchar(seqs))
    }
    message("[DEBUG] peptide_lengths: using all detected proteins, total peptides = ", length(lengths))
  }
  lengths
})

# ---------- 条形图输出（因子化x轴，无空隙，防止标签截断） ----------
output$intersection_peptide_length_hist <- renderPlot({
  lengths <- peptide_lengths()
  if (is.null(lengths) || length(lengths) == 0) {
    plot.new()
    text(0.5, 0.5, "Peptide sequence data not available", family = "Arial", cex = 1.5)
    return()
  }
  
  selected <- selected_protein_for_hist()
  data <- intersection_data()
  if (is.null(data)) return()
  
  if (!is.null(selected)) {
    n_peptides <- length(lengths)
    plot_title <- paste0("Peptide Length Distribution for Protein: ", selected,
                         "\n(", n_peptides, " peptides)")
    message("[DEBUG] rendering barplot for protein ", selected, ", ", n_peptides, " peptides")
  } else {
    n_proteins <- length(data$detected_proteins)
    n_peptides <- length(lengths)
    plot_title <- paste0("Peptide Length Distribution\n(",
                         n_proteins, " proteins, ", n_peptides, " peptides)")
    message("[DEBUG] rendering barplot for all detected proteins: ", n_proteins, " proteins, ", n_peptides, " peptides")
  }
  
  # 统计频数
  freq_table <- as.data.frame(table(Length = lengths), stringsAsFactors = FALSE)
  freq_table$Length <- as.integer(as.character(freq_table$Length))
  freq_table <- freq_table[order(freq_table$Length), ]
  message("[DEBUG] barplot freq_table: ")
  print(freq_table)
  
  # x轴因子，按顺序
  freq_table$LengthFactor <- factor(freq_table$Length, levels = sort(unique(freq_table$Length)))
  
  p <- ggplot(freq_table, aes(x = LengthFactor, y = Freq)) +
    geom_col(fill = "#3498db", alpha = 0.8, width = 0.7) +
    labs(title = plot_title,
         x = "Peptide Length (characters)", y = "Frequency") +
    theme_bw(base_family = "Arial") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 14)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))  # 顶部留白15%，防止标签被截断
  
  # 在柱子上方添加计数标签（当长度种类 ≤30 时）
  if (nrow(freq_table) <= 30) {
    p <- p + geom_text(aes(label = Freq), vjust = -0.5, size = 6, color = "black", family = "Arial")
  }
  
  print(p)
}, height = 500)  # 设置高度，使图形更方正

# ---------- 下载条形图（同样使用条形图） ----------
output$download_peptide_length_hist <- downloadHandler(
  filename = function() {
    selected <- selected_protein_for_hist()
    if (!is.null(selected)) {
      paste0("peptide_length_", selected, "_", Sys.Date(), ".png")
    } else {
      paste0("peptide_length_all_", Sys.Date(), ".png")
    }
  },
  content = function(file) {
    lengths <- peptide_lengths()
    if (is.null(lengths) || length(lengths) == 0) {
      png(file, family = "Arial")
      plot.new(); text(0.5,0.5,"No data", family = "Arial", cex = 1.5); dev.off()
      return()
    }
    selected <- selected_protein_for_hist()
    data <- intersection_data()
    if (is.null(data)) return()
    
    if (!is.null(selected)) {
      n_peptides <- length(lengths)
      plot_title <- paste0("Peptide Length Distribution for Protein: ", selected,
                           " (", n_peptides, " peptides)")
    } else {
      n_proteins <- length(data$detected_proteins)
      n_peptides <- length(lengths)
      plot_title <- paste0("Peptide Length Distribution (",
                           n_proteins, " proteins, ", n_peptides, " peptides)")
    }
    
    freq_table <- as.data.frame(table(Length = lengths), stringsAsFactors = FALSE)
    freq_table$Length <- as.integer(as.character(freq_table$Length))
    freq_table <- freq_table[order(freq_table$Length), ]
    freq_table$LengthFactor <- factor(freq_table$Length, levels = sort(unique(freq_table$Length)))
    
    p <- ggplot(freq_table, aes(x = LengthFactor, y = Freq)) +
      geom_col(fill = "#3498db", alpha = 0.8, width = 0.7) +
      labs(title = plot_title,
           x = "Peptide Length (characters)", y = "Frequency") +
      theme_bw(base_family = "Arial") +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14)
      ) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
    
    if (nrow(freq_table) <= 30) {
      p <- p + geom_text(aes(label = Freq), vjust = -0.5, size = 6, color = "black", family = "Arial")
    }
    
    ggsave(file, plot = p, width = 8, height = 6, dpi = 150)
    message("[DEBUG] download_peptide_length_hist: saved with title ", plot_title)
  }
)

message("[DEBUG] data_quality_intersection.R loaded successfully (barplot, square, no truncation)")