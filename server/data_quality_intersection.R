# server/data_quality_intersection.R
message("[DEBUG] data_quality_intersection.R loading... (eggNOG removed, CD annotation retained, fixed extract_group_prefix conflict, added FASTA sequence to table and download)")

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

# ---------- 从样本名提取分组前缀（已删除局部定义，复用 global.R 中的 extract_group_prefix） ----------

# ---------- 格式化显示名称：100_6 -> "100mM 6h"，200_12 -> "200mM 12h"，WT -> "WT" ----------
format_group_display_name <- function(prefix) {
  if (grepl("^WT$", prefix, ignore.case = TRUE)) return("WT")
  parts <- strsplit(prefix, "_")[[1]]
  concentration <- parts[1]
  time_val <- if (length(parts) >= 2) parts[2] else ""
  if (concentration %in% c("100", "200")) {
    return(paste0(concentration, "mM ", time_val, "h"))
  } else {
    return(prefix)
  }
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
  
  per_sample_counts_unsorted <- colSums(presence_matrix[, sample_cols, drop = FALSE])
  
  protein_sets <- lapply(sample_cols, function(s) {
    presence_matrix$Protein[presence_matrix[[s]] == 1]
  })
  names(protein_sets) <- sample_cols
  
  ordered_names <- sample_cols
  per_sample_counts <- per_sample_counts_unsorted[ordered_names]
  sorted_sample_cols <- ordered_names
  
  message("[DEBUG] intersection_data: protein_sets order kept as selected")
  
  shared_all <- sum(presence_matrix$Sum == length(sample_cols))
  
  unique_counts <- sapply(sorted_sample_cols, function(s) {
    sum(presence_matrix[[s]] == 1 & presence_matrix$Sum == 1)
  })
  unique_total <- sum(unique_counts)
  
  upset_df <- presence_matrix[, sample_cols, drop = FALSE]
  upset_df[] <- lapply(upset_df, as.logical)
  upset_df <- upset_df[, ordered_names, drop = FALSE]
  
  message("[DEBUG] intersection_data: built protein sets for Venn, sample count = ", length(protein_sets))
  message("[DEBUG] intersection_data: set sizes = ", paste(sapply(protein_sets, length), collapse = ", "))
  
  detected_proteins <- presence_matrix$Protein[presence_matrix$Sum > 0]
  message("[DEBUG] intersection_data: detected proteins with Sum>0 = ", length(detected_proteins))
  
  message("[DEBUG] intersection_data: total proteins = ", n_total, 
          ", shared_all = ", shared_all, 
          ", unique_per_sample = ", paste(unique_counts, collapse = ", "))
  
  list(
    presence_matrix = presence_matrix,
    samples = sorted_sample_cols,
    total_proteins = n_total,
    shared_all = shared_all,
    unique_counts = unique_counts,
    unique_total = unique_total,
    per_sample_counts = per_sample_counts,
    detected_proteins = detected_proteins,
    protein_sets = protein_sets,
    upset_df = upset_df
  )
})

# ---------- 按前缀分组（用于韦恩图） ----------
grouped_sample_sets <- reactive({
  data <- intersection_data()
  if (is.null(data)) return(NULL)
  samples <- data$samples
  groups <- sapply(samples, extract_group_prefix)   # 使用全局函数
  unique_groups <- unique(groups)
  group_list <- lapply(unique_groups, function(g) {
    members <- samples[groups == g]
    sets <- data$protein_sets[members]
    names(sets) <- members
    sets
  })
  names(group_list) <- unique_groups
  message("[DEBUG] grouped_sample_sets: groups = ", paste(unique_groups, collapse = ", "))
  group_list
})

# ---------- 蛋白存在矩阵表格（添加 CD 注释，移除 eggNOG，尝试加入 FASTA 序列） ----------
output$intersection_protein_table <- DT::renderDataTable({
  data <- intersection_data()
  if (is.null(data)) {
    return(DT::datatable(data.frame(Message = "Select at least 2 samples in the Missing Heatmap."), 
                         options = list(dom = 't'), rownames = FALSE))
  }
  df <- data$presence_matrix
  message("[DEBUG] intersection_protein_table: original dims = ", nrow(df), " x ", ncol(df))
  if (nrow(df) == 0) {
    return(DT::datatable(data.frame(Message = "No rows")))
  }
  # eggNOG 注释已移除
  message("[DEBUG] intersection_protein_table: eggNOG annotation merge removed, keeping CD only.")
  # 合并 CD-Search 注释
  if (exists("add_cd_to_table", mode = "function")) {
    tryCatch({
      df <- add_cd_to_table(df, id_col_name = "Protein")
      message("[DEBUG] intersection_protein_table: after CD merge, dims = ", nrow(df), " x ", ncol(df))
    }, error = function(e) {
      message("[ERROR] add_cd_to_table failed: ", conditionMessage(e))
    })
  }
  
  # ---------- 尝试添加 FASTA 序列到表格显示 ----------
  message("[DEBUG] intersection_protein_table: checking for FASTA sequence data for display")
  fasta_df <- NULL
  tryCatch({
    fasta_df <- cd_fasta()   # 调用 cd_search.R 中的反应式
  }, error = function(e) {
    message("[DEBUG] intersection_protein_table: cd_fasta() call failed - ", conditionMessage(e))
  })
  
  if (!is.null(fasta_df) && nrow(fasta_df) > 0) {
    message("[DEBUG] intersection_protein_table: FASTA data found with ", nrow(fasta_df), " sequences")
    if ("ID" %in% colnames(fasta_df) && "Sequence" %in% colnames(fasta_df)) {
      fasta_df$Protein <- fasta_df$ID
      fasta_sub <- fasta_df[, c("Protein", "Sequence"), drop = FALSE]
      df <- merge(df, fasta_sub, by = "Protein", all.x = TRUE, sort = FALSE)
      message("[DEBUG] intersection_protein_table: after adding FASTA sequence for display, dims = ", nrow(df), " x ", ncol(df))
      n_matched <- sum(!is.na(df$Sequence))
      message("[DEBUG] intersection_protein_table: sequences matched for ", n_matched, " proteins")
    } else {
      message("[DEBUG] intersection_protein_table: FASTA data format unexpected, columns: ", paste(colnames(fasta_df), collapse = ", "))
    }
  } else {
    message("[DEBUG] intersection_protein_table: no FASTA data available, adding placeholder")
    df$Sequence <- "Sequence not available (upload FASTA for CD-Search)"
    message("[DEBUG] intersection_protein_table: added placeholder Sequence column")
  }
  
  dt <- DT::datatable(df, options = list(pageLength = 10, server = FALSE, scrollX = TRUE), rownames = FALSE)
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

# ---------- 下载蛋白表 CSV（合并注释，移除 eggNOG，尝试加入 FASTA 序列） ----------
output$download_intersection_proteins <- downloadHandler(
  filename = function() paste0("protein_presence_matrix_", Sys.Date(), ".csv"),
  content = function(file) {
    data <- intersection_data()
    if (is.null(data)) {
      write.csv(data.frame(Error = "Not enough samples selected"), file, row.names = FALSE)
      return()
    }
    df <- data$presence_matrix
    message("[DEBUG] download_intersection_proteins: original dims = ", nrow(df), " x ", ncol(df))
    
    # 合并 CD 注释
    if (exists("add_cd_to_table", mode = "function")) {
      df <- tryCatch(add_cd_to_table(df, id_col_name = "Protein"), error = function(e) df)
      message("[DEBUG] download_intersection_proteins: after CD merge, dims = ", nrow(df), " x ", ncol(df))
    }
    
    # ---------- 尝试添加 FASTA 序列 ----------
    message("[DEBUG] download_intersection_proteins: checking for FASTA sequence data")
    fasta_df <- NULL
    tryCatch({
      fasta_df <- cd_fasta()   # 尝试调用 cd_search.R 中的反应式
    }, error = function(e) {
      message("[DEBUG] download_intersection_proteins: cd_fasta() call failed - ", conditionMessage(e))
    })
    
    if (!is.null(fasta_df) && nrow(fasta_df) > 0) {
      message("[DEBUG] download_intersection_proteins: FASTA data found with ", nrow(fasta_df), " sequences")
      # 确保有 Protein 列用于匹配（cd_fasta 返回 ID 和 Sequence）
      if ("ID" %in% colnames(fasta_df) && "Sequence" %in% colnames(fasta_df)) {
        fasta_df$Protein <- fasta_df$ID
        # 只保留 Protein 和 Sequence 列，避免重复
        fasta_sub <- fasta_df[, c("Protein", "Sequence"), drop = FALSE]
        # 左连接，保留所有蛋白，匹配不到则 Sequence 为 NA
        df <- merge(df, fasta_sub, by = "Protein", all.x = TRUE, sort = FALSE)
        message("[DEBUG] download_intersection_proteins: after adding FASTA sequence, dims = ", nrow(df), " x ", ncol(df))
        # 检查匹配情况
        n_matched <- sum(!is.na(df$Sequence))
        message("[DEBUG] download_intersection_proteins: sequences matched for ", n_matched, " proteins")
      } else {
        message("[DEBUG] download_intersection_proteins: FASTA data format unexpected, columns: ", paste(colnames(fasta_df), collapse = ", "))
      }
    } else {
      message("[DEBUG] download_intersection_proteins: no FASTA data available (did you upload a FASTA for CD-Search?)")
      # 如果用户没有上传 FASTA，添加一个提示列
      df$Sequence <- "Sequence not available (upload FASTA for CD-Search)"
      message("[DEBUG] download_intersection_proteins: added placeholder Sequence column")
    }
    
    # 写入 CSV
    write.csv(df, file, row.names = FALSE)
    message("[DEBUG] download_intersection_proteins: saved CSV with ", nrow(df), " rows and ", ncol(df), " columns")
  }
)

# ========== 韦恩图 UI 生成 ==========
output$group_venn_plots_ui <- renderUI({
  groups <- grouped_sample_sets()
  if (is.null(groups) || length(groups) == 0) return(div("No groups to display."))
  plot_ids <- paste0("venn_", names(groups))
  plot_list <- lapply(seq_along(groups), function(i) {
    gname <- names(groups)[i]
    sets <- groups[[gname]]
    total_proteins <- length(Reduce(union, sets))
    display_name <- format_group_display_name(gname)
    column(width = 4,
           div(
             plotOutput(plot_ids[i], height = "280px"),
             div(style = "text-align: center; font-weight: bold; font-size: 14px; margin-top: 5px;",
                 paste0(display_name, " (", total_proteins, ")"))
           )
    )
  })
  do.call(tagList, list(fluidRow(plot_list)))
})

# ========== 绘制每个分组的韦恩图（正圆，小尺寸，内部标签） ==========
observe({
  groups <- grouped_sample_sets()
  if (is.null(groups)) return()
  for (gname in names(groups)) {
    local({
      local_gname <- gname
      output[[paste0("venn_", local_gname)]] <- renderPlot({
        sets <- groups[[local_gname]]
        n <- length(sets)
        message("[DEBUG] venn_", local_gname, ": plotting ", n, " sets")
        if (n < 2 || n > 5) {
          plot.new()
          text(0.5, 0.5, if (n < 2) "Need at least 2 samples" else "Too many samples for Venn")
          return()
        }
        par(pty = "s")
        border_colors <- if (n <= 8) {
          RColorBrewer::brewer.pal(max(3, n), "Dark2")[1:n]
        } else {
          rainbow(n)
        }
        labels <- paste0("R", 1:n)
        if (n == 2) {
          VennDiagram::draw.pairwise.venn(
            area1 = length(sets[[1]]), area2 = length(sets[[2]]),
            cross.area = length(intersect(sets[[1]], sets[[2]])),
            category = labels,
            col = border_colors, lty = "solid", lwd = 2.5,
            fill = NA, alpha = 0,
            cat.cex = 1.2, cex = 1.2,
            cat.pos = c(-30, 30), cat.dist = 0.03,
            margin = 0.05
          )
        } else if (n == 3) {
          VennDiagram::draw.triple.venn(
            area1 = length(sets[[1]]), area2 = length(sets[[2]]), area3 = length(sets[[3]]),
            n12 = length(intersect(sets[[1]], sets[[2]])),
            n23 = length(intersect(sets[[2]], sets[[3]])),
            n13 = length(intersect(sets[[1]], sets[[3]])),
            n123 = length(Reduce(intersect, sets)),
            category = labels,
            col = border_colors, lty = "solid", lwd = 2.5,
            fill = rep(NA, 3), alpha = rep(0, 3),
            cat.cex = 1.2, cex = 1.2,
            cat.pos = c(-30, 30, 180), cat.dist = 0.05,
            margin = 0.05
          )
        } else if (n == 4) {
          VennDiagram::draw.quad.venn(
            area1 = length(sets[[1]]), area2 = length(sets[[2]]), area3 = length(sets[[3]]), area4 = length(sets[[4]]),
            n12 = length(intersect(sets[[1]], sets[[2]])),
            n13 = length(intersect(sets[[1]], sets[[3]])),
            n14 = length(intersect(sets[[1]], sets[[4]])),
            n23 = length(intersect(sets[[2]], sets[[3]])),
            n24 = length(intersect(sets[[2]], sets[[4]])),
            n34 = length(intersect(sets[[3]], sets[[4]])),
            n123 = length(Reduce(intersect, sets[1:3])),
            n124 = length(Reduce(intersect, sets[c(1,2,4)])),
            n134 = length(Reduce(intersect, sets[c(1,3,4)])),
            n234 = length(Reduce(intersect, sets[2:4])),
            n1234 = length(Reduce(intersect, sets)),
            col = border_colors, lty = "solid", lwd = 2.5,
            fill = rep(NA, 4), alpha = rep(0, 4),
            category = labels,
            cat.cex = 1.0, cex = 1.0,
            cat.pos = c(-15, 15, 180, 180), cat.dist = 0.05,
            margin = 0.05
          )
        } else if (n == 5) {
          VennDiagram::draw.quintuple.venn(
            area1 = length(sets[[1]]), area2 = length(sets[[2]]), area3 = length(sets[[3]]), area4 = length(sets[[4]]), area5 = length(sets[[5]]),
            n12 = length(intersect(sets[[1]], sets[[2]])),
            n13 = length(intersect(sets[[1]], sets[[3]])),
            n14 = length(intersect(sets[[1]], sets[[4]])),
            n15 = length(intersect(sets[[1]], sets[[5]])),
            n23 = length(intersect(sets[[2]], sets[[3]])),
            n24 = length(intersect(sets[[2]], sets[[4]])),
            n25 = length(intersect(sets[[2]], sets[[5]])),
            n34 = length(intersect(sets[[3]], sets[[4]])),
            n35 = length(intersect(sets[[3]], sets[[5]])),
            n45 = length(intersect(sets[[4]], sets[[5]])),
            n123 = length(Reduce(intersect, sets[1:3])),
            n124 = length(Reduce(intersect, sets[c(1,2,4)])),
            n125 = length(Reduce(intersect, sets[c(1,2,5)])),
            n134 = length(Reduce(intersect, sets[c(1,3,4)])),
            n135 = length(Reduce(intersect, sets[c(1,3,5)])),
            n145 = length(Reduce(intersect, sets[c(1,4,5)])),
            n234 = length(Reduce(intersect, sets[2:4])),
            n235 = length(Reduce(intersect, sets[c(2,3,5)])),
            n245 = length(Reduce(intersect, sets[c(2,4,5)])),
            n345 = length(Reduce(intersect, sets[3:5])),
            n1234 = length(Reduce(intersect, sets[1:4])),
            n1235 = length(Reduce(intersect, sets[c(1,2,3,5)])),
            n1245 = length(Reduce(intersect, sets[c(1,2,4,5)])),
            n1345 = length(Reduce(intersect, sets[c(1,3,4,5)])),
            n2345 = length(Reduce(intersect, sets[2:5])),
            n12345 = length(Reduce(intersect, sets)),
            col = border_colors, lty = "solid", lwd = 2.5,
            fill = rep(NA, 5), alpha = rep(0, 5),
            category = labels,
            cat.cex = 0.9, cex = 0.9,
            cat.pos = c(-20, 20, 180, 0, 0), cat.dist = 0.04,
            margin = 0.05
          )
        }
        message("[DEBUG] venn_", local_gname, ": done")
      })
    })
  }
})

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
      data.frame(`Master protein IDs` = pid, `Peptide Sequences` = paste(seqs, collapse = "; "), stringsAsFactors = FALSE, check.names = FALSE)
    })
    df_merged <- do.call(rbind, merged_list)
    if (is.null(df_merged) || nrow(df_merged) == 0) return(DT::datatable(data.frame(Message = "No peptide sequences found.")))
    dt <- DT::datatable(df_merged, options = list(pageLength = 25, scrollX = TRUE), rownames = FALSE)
    return(dt)
  } else {
    seq_list <- lapply(seq_along(peptide_seq), function(i) {
      pid <- names(peptide_seq)[i]
      seqs <- trimws(unlist(strsplit(peptide_seq[i], ";")))
      seqs <- seqs[seqs != ""]
      if (length(seqs) == 0) return(NULL)
      detected_samples <- get_detected_samples(pid)
      data.frame(`Master protein IDs` = pid, `Peptide Sequence` = seqs, `Detected in Samples` = detected_samples, stringsAsFactors = FALSE, check.names = FALSE)
    })
    df_seq <- do.call(rbind, seq_list)
    if (is.null(df_seq) || nrow(df_seq) == 0) return(DT::datatable(data.frame(Message = "No peptide sequences found.")))
    return(DT::datatable(df_seq, options = list(pageLength = 25, scrollX = TRUE), rownames = FALSE))
  }
})

# ---------- 下载肽段 CSV ----------
output$download_intersection_peptides <- downloadHandler(
  filename = function() "shared_unique_peptides.csv",
  content = function(file) {
    data <- intersection_data()
    if (is.null(data)) return()
    mode <- input$peptide_display_mode
    target_proteins <- data$detected_proteins
    clean <- rv$clean_data
    if (is.null(clean) || !"Master protein IDs" %in% colnames(clean)) return()
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
          data.frame(`Master protein IDs` = pid, `Peptide Sequences` = paste(seqs, collapse = "; "), stringsAsFactors = FALSE, check.names = FALSE)
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
          detected_samples <- get_detected_samples(pid)
          data.frame(`Master protein IDs` = pid, `Peptide Sequence` = seqs, `Detected in Samples` = detected_samples, stringsAsFactors = FALSE, check.names = FALSE)
        })
        df_out <- do.call(rbind, seq_list)
      }
      write.csv(df_out, file, row.names = FALSE)
    } else write.csv(sub, file, row.names = FALSE)
  }
)

# ========== 肽段长度数据（所有检测蛋白） ==========
peptide_lengths <- reactive({
  data <- intersection_data()
  if (is.null(data)) return(NULL)
  clean <- rv$clean_data
  if (is.null(clean) || !"Peptide sequences" %in% colnames(clean)) return(NULL)
  seq_col <- grep("^Peptide sequences$", colnames(clean), ignore.case = TRUE, value = TRUE)
  if (length(seq_col) == 0) {
    message("[DEBUG] peptide_lengths: 'Peptide sequences' column not found")
    return(NULL)
  }
  seq_col <- seq_col[1]
  target_proteins <- data$detected_proteins
  idx <- match(target_proteins, clean$`Master protein IDs`)
  peptide_seq <- clean[[seq_col]][idx]
  names(peptide_seq) <- clean$`Master protein IDs`[idx]
  lengths <- c()
  for (pid in names(peptide_seq)) {
    seqs <- trimws(unlist(strsplit(peptide_seq[pid], ";")))
    seqs <- seqs[seqs != ""]
    lengths <- c(lengths, nchar(seqs))
  }
  message("[DEBUG] peptide_lengths: total peptides = ", length(lengths))
  if (length(lengths) == 0) return(NULL)
  lengths
})

# ---------- 条形图输出（所有蛋白，无网格线） ----------
output$intersection_peptide_length_hist <- renderPlot({
  lengths <- peptide_lengths()
  if (is.null(lengths) || length(lengths) == 0) {
    plot.new()
    text(0.5, 0.5, "Peptide sequence data not available\n(No 'Peptide sequences' column in MaxQuant output)", family = "Arial", cex = 1.2)
    return()
  }
  data <- intersection_data()
  if (is.null(data)) return()
  n_proteins <- length(data$detected_proteins)
  n_peptides <- length(lengths)
  plot_title <- paste0("Peptide Length Distribution\n(", n_proteins, " proteins, ", n_peptides, " peptides)")
  freq_table <- as.data.frame(table(Length = lengths), stringsAsFactors = FALSE)
  freq_table$Length <- as.integer(as.character(freq_table$Length))
  freq_table <- freq_table[order(freq_table$Length), ]
  freq_table$LengthFactor <- factor(freq_table$Length, levels = sort(unique(freq_table$Length)))
  p <- ggplot(freq_table, aes(x = LengthFactor, y = Freq)) +
    geom_col(fill = "#3498db", alpha = 0.8, width = 0.7) +
    labs(title = plot_title, x = "Peptide Length (characters)", y = "Frequency") +
    theme_bw(base_family = "Arial") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 14),
      panel.grid = element_blank()
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  if (nrow(freq_table) <= 30) {
    p <- p + geom_text(aes(label = Freq), vjust = -0.5, size = 6, color = "black", family = "Arial")
  }
  message("[DEBUG] intersection_peptide_length_hist: rendered for all proteins (no grid)")
  print(p)
}, height = 500)

# ---------- 下载条形图 ----------
output$download_peptide_length_hist <- downloadHandler(
  filename = function() paste0("peptide_length_all_", Sys.Date(), ".png"),
  content = function(file) {
    lengths <- peptide_lengths()
    if (is.null(lengths) || length(lengths) == 0) {
      png(file, family = "Arial")
      plot.new(); text(0.5,0.5,"No data", family = "Arial", cex = 1.5); dev.off()
      return()
    }
    data <- intersection_data()
    if (is.null(data)) return()
    n_proteins <- length(data$detected_proteins)
    n_peptides <- length(lengths)
    plot_title <- paste0("Peptide Length Distribution\n(", n_proteins, " proteins, ", n_peptides, " peptides)")
    freq_table <- as.data.frame(table(Length = lengths), stringsAsFactors = FALSE)
    freq_table$Length <- as.integer(as.character(freq_table$Length))
    freq_table <- freq_table[order(freq_table$Length), ]; freq_table$LengthFactor <- factor(freq_table$Length, levels = sort(unique(freq_table$Length)))
    p <- ggplot(freq_table, aes(x = LengthFactor, y = Freq)) +
      geom_col(fill = "#3498db", alpha = 0.8, width = 0.7) +
      labs(title = plot_title, x = "Peptide Length (characters)", y = "Frequency") +
      theme_bw(base_family = "Arial") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
            axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
            axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
            axis.text.y = element_text(size = 14),
            panel.grid = element_blank()) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
    if (nrow(freq_table) <= 30) p <- p + geom_text(aes(label = Freq), vjust = -0.5, size = 6, color = "black", family = "Arial")
    ggsave(file, plot = p, width = 8, height = 6, dpi = 150)
  }
)

message("[DEBUG] data_quality_intersection.R loaded successfully (extract_group_prefix conflict resolved, table and download both support FASTA sequence)")