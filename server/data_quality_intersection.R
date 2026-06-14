# server/data_quality_intersection.R
message("[DEBUG] data_quality_intersection.R loading... (UpSetR, order by list, height 800)")

# ---------- 监控热图选中的样本 ----------
observe({
  sel <- selected_samples()
  message("[DEBUG] heatmap sample selection: ", 
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

# ---------- 自然排序辅助函数 ----------
natural_order <- function(x) {
  if (length(x) == 0) return(integer(0))
  parts <- strsplit(as.character(x), "_")
  max_len <- max(sapply(parts, length))
  num_mat <- t(sapply(parts, function(p) {
    nums <- suppressWarnings(as.numeric(p))
    nums[is.na(nums)] <- 0
    if (length(nums) < max_len) c(nums, rep(0, max_len - length(nums)))
    else nums
  }))
  num_df <- as.data.frame(num_mat)
  do.call(order, num_df)
}

# ---------- 基于热图选择自动计算数据 ----------
intersection_data <- reactive({
  samples <- selected_samples()
  message("[DEBUG] intersection_data: reactive triggered. Heatmap selection: ", 
          if (is.null(samples)) "NULL" else paste0("(", length(samples), ") ", paste(samples, collapse = ", ")))
  
  if (is.null(samples) || length(samples) < 2) {
    message("[DEBUG] intersection_data: insufficient samples (<2), returning NULL")
    return(NULL)
  }
  
  presence_matrix <- get_presence_matrix(samples)
  if (is.null(presence_matrix)) {
    message("[DEBUG] intersection_data: presence matrix is NULL, returning NULL")
    return(NULL)
  }
  
  sample_cols <- setdiff(colnames(presence_matrix), c("Protein", "Sum"))
  n_total <- nrow(presence_matrix)
  
  # 计算各样本蛋白检测数目（未排序）
  per_sample_counts_unsorted <- colSums(presence_matrix[, sample_cols, drop = FALSE])
  
  # 构建蛋白集合列表
  protein_sets <- lapply(sample_cols, function(s) {
    presence_matrix$Protein[presence_matrix[[s]] == 1]
  })
  names(protein_sets) <- sample_cols
  
  # 自然排序：WT_1, WT_2, ..., 100_6_1, 100_6_2, ...
  ordered_names <- names(protein_sets)[natural_order(names(protein_sets))]
  protein_sets <- protein_sets[ordered_names]          # 列表顺序固定
  
  # per_sample_counts 也按相同顺序
  per_sample_counts <- per_sample_counts_unsorted[ordered_names]
  
  # 最终样本列表
  sorted_sample_cols <- ordered_names
  
  message("[DEBUG] intersection_data: protein_sets reordered to: ", paste(names(protein_sets), collapse = ", "))
  message("[DEBUG] intersection_data: per_sample_counts reordered accordingly")
  
  shared_all <- sum(presence_matrix$Sum == length(sample_cols))
  
  unique_counts <- sapply(sorted_sample_cols, function(s) {
    sum(presence_matrix[[s]] == 1 & presence_matrix$Sum == 1)
  })
  unique_total <- sum(unique_counts)
  
  upset_df <- presence_matrix[, sample_cols, drop = FALSE]
  upset_df[] <- lapply(upset_df, as.logical)
  upset_df <- upset_df[, ordered_names, drop = FALSE]
  
  message("[DEBUG] intersection_data: built protein sets for UpSet, sample count = ", length(protein_sets))
  message("[DEBUG] intersection_data: set sizes = ", paste(sapply(protein_sets, length), collapse = ", "))
  
  detected_proteins <- presence_matrix$Protein[presence_matrix$Sum > 0]
  message("[DEBUG] intersection_data: detected proteins with Sum>0 = ", length(detected_proteins))
  
  message("[DEBUG] intersection_data: total proteins = ", n_total, 
          ", shared_all = ", shared_all, 
          ", unique_per_sample = ", paste(unique_counts, collapse = ", "))
  
  list(
    presence_matrix = presence_matrix,
    samples = sorted_sample_cols,          # 排序后的样本列表
    total_proteins = n_total,
    shared_all = shared_all,
    unique_counts = unique_counts,
    unique_total = unique_total,
    per_sample_counts = per_sample_counts,  # 已排序
    detected_proteins = detected_proteins,
    protein_sets = protein_sets,
    upset_df = upset_df
  )
})

# ---------- 摘要 ----------
output$intersection_summary <- renderText({
  data <- intersection_data()
  if (is.null(data)) {
    return("Select at least 2 samples in the Missing Heatmap above. The table and plots will update automatically.")
  }
  
  lines <- c()
  lines <- c(lines, paste("Selected samples:", paste(data$samples, collapse = ", ")))
  lines <- c(lines, paste("Total proteins evaluated:", data$total_proteins))
  lines <- c(lines, paste("Proteins detected in all selected samples (shared):", data$shared_all))
  lines <- c(lines, paste("Proteins detected in at least one sample:", length(data$detected_proteins)))
  lines <- c(lines, "")
  lines <- c(lines, "Proteins detected per sample:")
  for (s in names(data$per_sample_counts)) {
    lines <- c(lines, sprintf("  %s: %d", s, data$per_sample_counts[s]))
  }
  lines <- c(lines, "")
  lines <- c(lines, "Unique proteins per sample (detected in exactly one sample):")
  if (data$unique_total == 0) {
    lines <- c(lines, "  None (all proteins appear in at least 2 samples)")
  } else {
    for (s in names(data$unique_counts)) {
      lines <- c(lines, sprintf("  %s: %d", s, data$unique_counts[s]))
    }
    lines <- c(lines, paste("Total unique proteins:", data$unique_total))
  }
  
  paste(lines, collapse = "\n")
})

# ---------- 蛋白存在矩阵表格 ----------
output$intersection_protein_table <- DT::renderDataTable({
  data <- intersection_data()
  if (is.null(data)) {
    return(DT::datatable(data.frame(Message = "Select at least 2 samples in the Missing Heatmap."), 
                         options = list(dom = 't'), rownames = FALSE))
  }
  
  df <- data$presence_matrix
  message("[DEBUG] intersection_protein_table: rendering datatable, nrow = ", nrow(df),
          ", ncol = ", ncol(df))
  
  if (nrow(df) == 0) {
    return(DT::datatable(data.frame(Message = "No rows")))
  }
  
  dt <- DT::datatable(
    df,
    options = list(pageLength = 10, server = FALSE),
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
}, server = FALSE)

# ---------- 下载蛋白表 CSV ----------
output$download_intersection_proteins <- downloadHandler(
  filename = function() paste0("protein_presence_matrix_", Sys.Date(), ".csv"),
  content = function(file) {
    data <- intersection_data()
    if (is.null(data)) {
      write.csv(data.frame(Error = "Not enough samples selected"), file, row.names = FALSE)
      return()
    }
    write.csv(data$presence_matrix, file, row.names = FALSE)
    message("[DEBUG] download_intersection_proteins: saved")
  }
)

# ========== UpSet 图（使用已排序列表，不传 sets，高度800） ==========
output$intersection_upset_plot <- renderPlot({
  start_time <- Sys.time()
  data <- intersection_data()
  if (is.null(data)) {
    plot.new()
    text(0.5, 0.5, "Select at least 2 samples in the Missing Heatmap.")
    output$intersection_upset_time <- renderText({ "请先选择样本" })
    return()
  }
  
  sets <- data$protein_sets
  n_sets <- length(sets)
  set_order <- names(sets)   # 已按自然顺序
  message("[DEBUG] drawing UpSetR plot for ", n_sets, " samples, list order: ", paste(set_order, collapse = ", "))
  
  # 选择颜色
  if (n_sets <= 8) {
    set_colors <- RColorBrewer::brewer.pal(max(3, n_sets), "Set2")[1:n_sets]
  } else {
    set_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_sets)
  }
  
  # 动态调整文字大小
  if (n_sets > 10) {
    text_scale <- c(1.8, 1.6, 1.6, 1.4, 1.0, 1.4)
    message("[DEBUG] n_sets > 10: text_scale[5] set to 1.0")
  } else {
    text_scale <- c(1.8, 1.6, 1.6, 1.4, 2.0, 1.4)
    message("[DEBUG] n_sets <= 10: default text_scale")
  }
  
  # Arial 字体
  old_par <- par(family = "Arial", cex.main = 1.8, cex.lab = 1.6, cex.axis = 1.4)
  on.exit(par(old_par), add = TRUE)
  
  tryCatch({
    # 注意：这里不传 sets 参数，完全依赖 fromList(sets) 中列表的顺序
    print(
      UpSetR::upset(
        UpSetR::fromList(sets),
        nsets = n_sets,
        order.by = "freq",
        decreasing = TRUE,
        main.bar.color = "#3498db",
        sets.bar.color = set_colors,
        matrix.color = "#34495e",
        mainbar.y.label = "Intersection Size",
        sets.x.label = "Set Size",
        set_size.angles = 45,
        text.scale = text_scale
      )
    )
    message("[DEBUG] intersection_upset_plot: UpSetR plot printed, order from list")
    elapsed <- Sys.time() - start_time
    output$intersection_upset_time <- renderText({
      sprintf("%.1f 秒", as.numeric(elapsed))
    })
    message(sprintf("[TIMING] UpSetR: %.2f seconds", as.numeric(elapsed)))
  }, error = function(e) {
    message("[ERROR] UpSetR: ", e$message)
    plot.new()
    text(0.5, 0.5, paste("UpSet error:", e$message))
  })
}, height = 800)

# ---------- 显示耗时 ----------
output$intersection_upset_time <- renderText({
  "计算中..."
})

# ---------- 下载 UpSet 图 ----------
output$download_intersection_upset <- downloadHandler(
  filename = function() paste0("UpSet_", Sys.Date(), ".png"),
  content = function(file) {
    data <- intersection_data()
    if (is.null(data)) return()
    sets <- data$protein_sets
    n_sets <- length(sets)
    if (n_sets <= 8) {
      set_colors <- RColorBrewer::brewer.pal(max(3, n_sets), "Set2")[1:n_sets]
    } else {
      set_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_sets)
    }
    
    if (n_sets > 10) {
      text_scale <- c(1.8, 1.6, 1.6, 1.4, 1.0, 1.4)
      png_height <- 1000
    } else {
      text_scale <- c(1.8, 1.6, 1.6, 1.4, 2.0, 1.4)
      png_height <- 800
    }
    
    png(file, width = 1200, height = png_height, res = 150)
    old_par <- par(family = "Arial", cex.main = 1.8, cex.lab = 1.6, cex.axis = 1.4)
    on.exit(par(old_par), add = TRUE)
    
    tryCatch({
      print(
        UpSetR::upset(
          UpSetR::fromList(sets),
          nsets = n_sets,
          order.by = "freq",
          decreasing = TRUE,
          main.bar.color = "#3498db",
          sets.bar.color = set_colors,
          matrix.color = "#34495e",
          mainbar.y.label = "Intersection Size",
          sets.x.label = "Set Size",
          set_size.angles = 45,
          text.scale = text_scale
        )
      )
    }, error = function(e) {
      message("[ERROR] download UpSetR: ", e$message)
    })
    dev.off()
    message("[DEBUG] download_intersection_upset: saved (order from list)")
  }
)

# ---------- 肽段表格 ----------
output$intersection_peptide_table <- DT::renderDataTable({
  data <- intersection_data()
  if (is.null(data)) {
    return(DT::datatable(data.frame(Message = "Select at least 2 samples to view peptide sequences.")))
  }
  
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
    dt <- DT::datatable(df_merged, 
                        options = list(pageLength = 25, scrollX = TRUE),
                        rownames = FALSE,
                        selection = 'single')
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
    if (is.null(data)) {
      write.csv(data.frame(Error = "No data"), file)
      return()
    }
    mode <- input$peptide_display_mode
    target_proteins <- data$detected_proteins
    clean <- rv$clean_data
    if (is.null(clean) || !"Master protein IDs" %in% colnames(clean)) {
      write.csv(data.frame(Error = "No clean data"), file)
      return()
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

# ========== 选中蛋白驱动条形图 ==========
selected_protein_for_hist <- reactive({
  if (input$peptide_display_mode != "merged") return(NULL)
  sel <- input$intersection_peptide_table_rows_selected
  if (is.null(sel) || length(sel) == 0) return(NULL)
  
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
  
  selected_row <- sel[1]
  if (selected_row < 1 || selected_row > nrow(df_merged)) return(NULL)
  pid <- df_merged[selected_row, 1]
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
  
  if (!is.null(selected) && !selected %in% target_proteins) selected <- NULL
  
  if (!is.null(selected)) {
    idx <- match(selected, clean$`Master protein IDs`)
    if (is.na(idx)) return(numeric(0))
    pep_seq <- clean$`Peptide sequences`[idx]
    seqs <- trimws(unlist(strsplit(pep_seq, ";")))
    seqs <- seqs[seqs != ""]
    lengths <- nchar(seqs)
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
  }
  lengths
})

# ---------- 条形图输出 ----------
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
  } else {
    n_proteins <- length(data$detected_proteins)
    n_peptides <- length(lengths)
    plot_title <- paste0("Peptide Length Distribution\n(",
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
  
  print(p)
}, height = 500)

# ---------- 下载条形图 ----------
output$download_peptide_length_hist <- downloadHandler(
  filename = function() {
    selected <- selected_protein_for_hist()
    if (!is.null(selected)) paste0("peptide_length_", selected, "_", Sys.Date(), ".png")
    else paste0("peptide_length_all_", Sys.Date(), ".png")
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
  }
)

message("[DEBUG] data_quality_intersection.R loaded successfully (order from list, height 800)")