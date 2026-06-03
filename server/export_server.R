# server/export_server.R
message("[DEBUG] export_server.R loading...")

# ---------- 辅助函数：组合图绘制（移植自旧平台） ----------
get_optimal_layout <- function(n_plots) {
  if (n_plots <= 0) return(list(ncol = 1, nrow = 1))
  if (n_plots == 1) return(list(ncol = 1, nrow = 1))
  if (n_plots == 2) return(list(ncol = 2, nrow = 1))
  if (n_plots == 3) return(list(ncol = 3, nrow = 1))
  if (n_plots == 4) return(list(ncol = 2, nrow = 2))
  if (n_plots == 5) return(list(ncol = 3, nrow = 2))
  if (n_plots == 6) return(list(ncol = 3, nrow = 2))
  if (n_plots <= 8) return(list(ncol = 4, nrow = 2))
  if (n_plots <= 9) return(list(ncol = 3, nrow = 3))
  if (n_plots <= 12) return(list(ncol = 4, nrow = 3))
  return(list(ncol = ceiling(sqrt(n_plots)), nrow = ceiling(n_plots / ceiling(sqrt(n_plots)))))
}

plot_volcano_core_combined <- function(df, fc_up, fc_down, p_cut, cols, point_size = 1.8) {
  req(nrow(df) > 0)
  cnt <- list(Up = sum(df$regulation == "Up", na.rm = TRUE),
              Down = sum(df$regulation == "Down", na.rm = TRUE),
              Increase = sum(df$regulation == "Increase", na.rm = TRUE),
              Decrease = sum(df$regulation == "Decrease", na.rm = TRUE))
  x1 <- -8.5; x2 <- 8.5; fixed_y_max <- 8
  p <- ggplot(df, aes(log2FC, log10P, color = regulation)) + 
    geom_point(size = point_size, alpha = 0.6) +
    scale_color_manual(values = cols) + 
    geom_vline(xintercept = log2(c(fc_down, fc_up)), lty = 2, color = "gray40", linewidth = 0.4) +
    geom_hline(yintercept = -log10(p_cut), lty = 2, color = "gray40", linewidth = 0.4) +
    coord_cartesian(xlim = c(x1, x2), ylim = c(0, fixed_y_max)) + 
    volcano_theme() +
    theme(legend.position = "none", plot.title = element_blank(), axis.title = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
          panel.border = element_blank(), plot.margin = margin(5, 5, 5, 5))
  p <- p + annotate("text", x = -4.5, y = 8.0, label = sprintf("%d", cnt$Down), color = cols$Down, fontface = "bold", size = 5, hjust = 1) +
    annotate("text", x = -4.5 + 0.1, y = 8.0, label = sprintf("(%d)", cnt$Decrease), color = cols$Decrease, fontface = "bold", size = 5, hjust = 0) +
    annotate("text", x = 4.5, y = 8.0, label = sprintf("%d", cnt$Up), color = cols$Up, fontface = "bold", size = 5, hjust = 1) +
    annotate("text", x = 4.5 + 0.1, y = 8.0, label = sprintf("(%d)", cnt$Increase), color = cols$Increase, fontface = "bold", size = 5, hjust = 0)
  return(p)
}

build_combined_plot <- function(results_list, fcu, fcd, pc, cols, main_title, sub_titles = NULL, point_size = 1.8) {
  if (length(results_list) == 0) return(NULL)
  n_plots <- length(results_list)
  layout_info <- get_optimal_layout(n_plots)
  ncol <- layout_info$ncol; nrow <- layout_info$nrow
  
  volcano_plots <- lapply(results_list, function(r) {
    plot_volcano_core_combined(r$data, fcu, fcd, pc, cols, point_size)
  })
  
  base_height_per_plot <- 4.0
  title_height <- 0.5
  title_space <- 1.0
  bottom_space <- 0.8
  top_title_space <- 0.5
  bottom_title_space <- 0.5
  
  blank <- ggdraw()
  row_plots <- list()
  row_heights <- c()
  top_titles <- vector("list", ncol)
  bottom_titles <- vector("list", ncol)
  
  for (r in seq_len(nrow)) {
    start_idx <- (r-1)*ncol + 1
    end_idx <- min(r*ncol, n_plots)
    cols_in_row <- end_idx - start_idx + 1
    col_units <- vector("list", cols_in_row)
    for (j in seq_len(cols_in_row)) {
      i <- start_idx + j - 1
      vp <- volcano_plots[[i]]
      st <- if (!is.null(sub_titles) && i <= length(sub_titles)) trimws(sub_titles[i]) else ""
      if (r == 1) {
        col_units[[j]] <- vp
        if (st != "") top_titles[[j]] <- st
      } else if (r == nrow) {
        col_units[[j]] <- vp
        if (st != "") bottom_titles[[j]] <- st
      } else {
        if (st != "") {
          tg <- ggdraw() + draw_label(st, size = 12, fontface = "bold")
          col_units[[j]] <- plot_grid(vp, tg, ncol = 1, rel_heights = c(base_height_per_plot, title_height))
        } else {
          col_units[[j]] <- plot_grid(vp, blank, ncol = 1, rel_heights = c(base_height_per_plot, title_height))
        }
      }
    }
    if (length(col_units) > 1) {
      row_grob <- plot_grid(plotlist = col_units, ncol = length(col_units), align = "hv", axis = "none")
    } else {
      row_grob <- col_units[[1]]
    }
    row_plots[[r]] <- row_grob
    row_heights[r] <- if (r == 1 || r == nrow) base_height_per_plot else base_height_per_plot + title_height
  }
  
  grid_height <- sum(row_heights)
  plot_height <- grid_height + title_space + bottom_space + top_title_space + bottom_title_space
  
  grid_plot <- plot_grid(plotlist = row_plots, ncol = 1, rel_heights = row_heights, align = "v", axis = "none")
  
  final_plot <- ggdraw() +
    draw_plot(grid_plot,
              x = 0.12, y = (bottom_space + bottom_title_space) / plot_height,
              width = 0.76, height = grid_height / plot_height)
  
  for (col_i in seq_len(ncol)) {
    txt <- top_titles[[col_i]]
    if (!is.null(txt) && txt != "") {
      x_pos <- 0.12 + (col_i - 0.5) * (0.76 / ncol)
      y_pos <- (bottom_space + bottom_title_space + grid_height + top_title_space/2) / plot_height
      final_plot <- final_plot + draw_label(txt, x = x_pos, y = y_pos, size = 12, fontface = "bold")
    }
  }
  for (col_i in seq_len(ncol)) {
    txt <- bottom_titles[[col_i]]
    if (!is.null(txt) && txt != "") {
      x_pos <- 0.12 + (col_i - 0.5) * (0.76 / ncol)
      y_pos <- (bottom_space + bottom_title_space/2) / plot_height
      final_plot <- final_plot + draw_label(txt, x = x_pos, y = y_pos, size = 12, fontface = "bold")
    }
  }
  
  final_plot <- final_plot + draw_label(main_title, x = 0.5, y = 1 - (title_space/2)/plot_height, size = 24, fontface = "bold")
  y_axis_y <- (bottom_space + bottom_title_space + grid_height/2) / plot_height
  final_plot <- final_plot + draw_label(expression(bold(-Log[10]~(P-Value))), x = 0.06, y = y_axis_y, angle = 90, size = 16, fontface = "bold")
  x_axis_y <- (bottom_space + bottom_title_space)/2 / plot_height
  final_plot <- final_plot + draw_label(expression(bold(Log[2]~(Fold~Change))), x = 0.5, y = x_axis_y, size = 16, fontface = "bold")
  
  final_plot <- final_plot + draw_grob(rectGrob(gp = gpar(lwd = 2, col = "black", fill = NA)),
                                       x = 0.12, y = (bottom_space + bottom_title_space) / plot_height,
                                       width = 0.76, height = grid_height / plot_height)
  
  if (ncol > 1) {
    for (j in 1:(ncol - 1)) {
      x_line <- 0.12 + j * (0.76 / ncol)
      final_plot <- final_plot + draw_grob(linesGrob(x = c(x_line, x_line),
                                                     y = c((bottom_space + bottom_title_space) / plot_height,
                                                           (bottom_space + bottom_title_space + grid_height) / plot_height),
                                                     gp = gpar(lwd = 1, col = "black")))
    }
  }
  if (nrow > 1) {
    y_cumulative <- (bottom_space + bottom_title_space) / plot_height
    for (k in 1:(nrow - 1)) {
      y_cumulative <- y_cumulative + row_heights[k] / plot_height
      final_plot <- final_plot + draw_grob(linesGrob(x = c(0.12, 0.88),
                                                     y = c(y_cumulative, y_cumulative),
                                                     gp = gpar(lwd = 1, col = "black")))
    }
  }
  
  return(final_plot)
}

# ---------- 子标题UI ----------
output$subplot_titles_ui <- renderUI({
  comps <- sorted_comps()
  if (length(comps) == 0) return(NULL)
  title_inputs <- lapply(seq_along(comps), function(i) {
    comp <- comps[[i]]
    textInputMax(paste0("subplot_title_", i), label = paste0("Sub-plot ", i, " (", comp$name, ")"), value = "", placeholder = paste0("Default: ", comp$name), maxlength = 25, allowed_pattern = "[^a-zA-Z0-9 _-]")
  })
  do.call(tagList, title_inputs)
})

get_subplot_title <- function(index, default_display_name) {
  title_input_id <- paste0("subplot_title_", index)
  custom_title <- input[[title_input_id]]
  if (is.null(custom_title) || trimws(custom_title) == "") default_display_name else trimws(custom_title)
}

# ---------- 所有比较的差异分析结果 ----------
all_analysis_results <- reactive({
  message("[DEBUG] export: computing all_analysis_results...")
  
  nd <- norm_data_full()
  if (is.null(nd)) {
    message("[DEBUG] export: norm_data_full() is NULL, cannot compute DE")
    return(NULL)
  }
  
  # 已经移除了 Fasta headers，无需再次移除
  
  fcu <- input$fc_up; fcd <- input$fc_down; pc <- as.numeric(input$p_cut)
  stat_method <- input$stat_method
  message("[DEBUG] export: params: FC_up=", fcu, ", FC_down=", fcd, ", p_cut=", pc, ", method=", stat_method)
  
  # 提取归一化后的强度列（以 Norm_ 开头）
  norm_cols <- grep("^Norm_LFQ intensity ", colnames(nd), value = TRUE)
  if (length(norm_cols) == 0) {
    message("[DEBUG] export: no Norm_LFQ intensity columns found")
    return(NULL)
  }
  
  # 可选：根据 unique peptides 过滤
  unique_col <- grep("^Unique peptides$", colnames(nd), value = TRUE)[1]
  if (!is.na(unique_col) && input$min_unique_pep > 1) {
    nd[[unique_col]] <- as.numeric(nd[[unique_col]])
    nd <- nd[nd[[unique_col]] >= input$min_unique_pep, ]
    message("[DEBUG] export: after unique peptide filter (>= ", input$min_unique_pep, "), nrow = ", nrow(nd))
  }
  
  results <- list()
  sorted_comp <- sorted_comps()
  message("[DEBUG] export: number of comparisons = ", length(sorted_comp))
  
  for (comp in sorted_comp) {
    treat_group <- comp$treat; ctrl_group <- comp$ctrl; comp_name <- comp$name
    treat_samples <- rv$groups[[treat_group]]; ctrl_samples <- rv$groups[[ctrl_group]]
    
    treat_cols <- paste0("Norm_LFQ intensity ", treat_samples)
    ctrl_cols <- paste0("Norm_LFQ intensity ", ctrl_samples)
    treat_cols <- treat_cols[treat_cols %in% colnames(nd)]
    ctrl_cols <- ctrl_cols[ctrl_cols %in% colnames(nd)]
    
    if (length(treat_cols) == 0 || length(ctrl_cols) == 0) {
      message("[DEBUG] export: skipping ", comp_name, " - missing Norm columns")
      next
    }
    
    # 保留注释列（已调整顺序）
    annotation_cols <- intersect(c("Protein IDs", "Majority protein IDs", "Master protein IDs", "Unique peptides"), colnames(nd))
    select_cols <- c(annotation_cols, treat_cols, ctrl_cols)
    sub_df <- nd[, select_cols, drop = FALSE]
    
    message("[DEBUG] export: running DE for ", comp_name, " with ", length(treat_cols), " treat and ", length(ctrl_cols), " ctrl columns")
    
    res <- tryCatch({
      run_de_analysis(
        data_subset = sub_df,
        treat_cols = treat_cols,
        ctrl_cols = ctrl_cols,
        fc_up = fcu,
        fc_down = fcd,
        p_cut = pc,
        stat_method = stat_method
      )
    }, error = function(e) {
      message("[ERROR] export: DE failed for ", comp_name, ": ", e$message)
      NULL
    })
    
    if (!is.null(res) && nrow(res) > 0) {
      results[[comp_name]] <- list(data = res, treat = treat_group, ctrl = ctrl_group, name = comp_name)
      message("[DEBUG] export: DE completed for ", comp_name, ", nrow=", nrow(res),
              ", Up=", sum(res$regulation == "Up"), ", Down=", sum(res$regulation == "Down"),
              ", Increase=", sum(res$regulation == "Increase"), ", Decrease=", sum(res$regulation == "Decrease"))
    } else {
      message("[DEBUG] export: DE result empty for ", comp_name)
    }
  }
  
  # 准备 Raw/Clean 数据（过滤掉多余的列）
  raw_data <- rv$raw_data
  clean_data <- rv$clean_data
  cols_to_remove <- c("Fasta headers", "Number of proteins", "Peptide counts (total)", 
                      "Peptide counts (razor+unique)", "Peptide counts (unique)")
  raw_data <- raw_data[, setdiff(colnames(raw_data), cols_to_remove), drop = FALSE]
  clean_data <- clean_data[, setdiff(colnames(clean_data), cols_to_remove), drop = FALSE]
  
  list(
    raw = raw_data,
    clean = clean_data,
    norm = nd,
    filtered = nd,
    unique_col = if (!is.na(unique_col)) unique_col else NA_character_,
    results = results
  )
})

# 更新导出比较选择器
observe({
  comp_names <- sapply(sorted_comps(), function(c) c$name)
  updateSelectInput(session, "export_comparisons", choices = comp_names, selected = comp_names)
})

# ---------- 单图下载 ----------
output$download_plot <- downloadHandler(
  filename = function() {
    s <- gsub("[^a-zA-Z0-9_]", "_", input$selected_comparison %||% "plot")
    paste0("Volcano_", s, "_", Sys.Date(), ".", input$plot_format)
  },
  content = function(f) {
    res <- volcano_de_result()
    if (is.null(res)) {
      png(f, 800, 600); plot(0,0,type="n"); text(0.5,0.5,"No data"); dev.off()
      return()
    }
    df <- res$data; comp_name <- res$name
    cols <- color_mapping_vector()
    width_val <- as.integer(input$plot_width); height_val <- as.integer(input$plot_height)
    if (is.na(width_val) || width_val < 5 || width_val > 30) width_val <- 10
    if (is.na(height_val) || height_val < 5 || height_val > 30) height_val <- 8
    point_size <- input$point_size
    download_title <- if (!is.null(input$download_single_title) && input$download_single_title != "") input$download_single_title else comp_name
    
    p <- ggplot(df, aes(log2FC, log10P, color = regulation)) + geom_point(size = point_size, alpha = 0.6) + scale_color_manual(values = cols) +
      geom_vline(xintercept = log2(c(input$fc_down, input$fc_up)), lty = 2, color = "gray40") +
      geom_hline(yintercept = -log10(as.numeric(input$p_cut)), lty = 2, color = "gray40") +
      labs(title = download_title, x = expression(Log[2]~"(Fold Change)"), y = expression(-Log[10]~"(P-Value)")) +
      volcano_theme() + guides(color = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 1))
    ggsave(f, plot = p, width = width_val, height = height_val, dpi = 300, device = input$plot_format)
  }
)

# ---------- 组合图下载（强制输出 PNG） ----------
output$download_combined_plot <- downloadHandler(
  filename = function() {
    title_clean <- gsub("[^a-zA-Z0-9_]", "_", input$combined_plot_title)
    if (title_clean == "" || is.null(title_clean)) title_clean <- "Combined_Volcano_Plots"
    paste0(title_clean, "_", Sys.Date(), ".png")
  },
  content = function(f) {
    progress <- shiny::Progress$new()
    progress$set(message = "Generating Combined Volcano Plot", value = 0)
    on.exit(progress$close())
    
    req(all_analysis_results())
    all_res <- all_analysis_results()
    results_list <- all_res$results
    if (length(results_list) == 0) {
      png(f, 800, 600); plot(0,0,type="n"); text(0.5,0.5,"No comparisons available"); dev.off()
      return()
    }
    
    fcu <- input$fc_up; fcd <- input$fc_down; pc <- as.numeric(input$p_cut)
    cols <- color_mapping_vector()
    custom_title <- if (is.null(input$combined_plot_title) || trimws(input$combined_plot_title) == "") {
      "Combined Volcano Plots"
    } else {
      trimws(input$combined_plot_title)
    }
    
    comp_names <- names(results_list)
    default_display <- sapply(results_list, function(r) r$name)
    sub_titles <- sapply(seq_along(results_list), function(i) get_subplot_title(i, default_display[i]))
    point_size <- input$point_size
    
    progress$set(detail = "Building combined plot...", value = 0.5)
    final_plot <- build_combined_plot(results_list, fcu, fcd, pc, cols, custom_title, sub_titles, point_size)
    if (is.null(final_plot)) {
      png(f, 800, 600); plot(0,0,type="n"); text(0.5,0.5,"No plots to combine"); dev.off()
      return()
    }
    
    n_plots <- length(results_list)
    layout_info <- get_optimal_layout(n_plots)
    ncol <- layout_info$ncol; nrow <- layout_info$nrow
    base_height_per_plot <- 4.0
    title_height <- 0.5
    title_space <- 1.0
    bottom_space <- 0.8
    top_title_space <- 0.5
    bottom_title_space <- 0.5
    row_heights <- c()
    for (r in seq_len(nrow)) {
      if (r == 1 || r == nrow) row_heights[r] <- base_height_per_plot
      else row_heights[r] <- base_height_per_plot + title_height
    }
    grid_height <- sum(row_heights)
    plot_height <- grid_height + title_space + bottom_space + top_title_space + bottom_title_space
    plot_width <- 4.5 * ncol
    
    ggsave(f, plot = final_plot, device = "png", width = plot_width, height = plot_height, dpi = 300, bg = "white")
    progress$set(detail = "Complete!", value = 1.0)
  }
)

# ---------- Excel 导出 ----------
output$download_excel <- downloadHandler(
  filename = function() paste0("Proteomics_Report_", Sys.Date(), ".xlsx"),
  content = function(f) {
    showNotification("The data is large. The Excel report may take 1–2 minutes. Please wait...", type = "warning", duration = 10)
    
    req(all_analysis_results(), input$export_comparisons)
    all_res <- all_analysis_results()
    selected_comps <- input$export_comparisons
    if (length(selected_comps) == 0) {
      showNotification("Please select at least one comparison.", type = "error")
      return()
    }
    message("[DEBUG] export_excel: selected comparisons: ", paste(selected_comps, collapse = ", "))
    
    progress <- shiny::Progress$new()
    progress$set(message = "Generating Excel report...", value = 0)
    on.exit(progress$close())
    
    wb <- createWorkbook()
    sanitize_sheet <- function(name) sanitize_name(name, max_len = 31)
    
    addWorksheet(wb, "Raw Data")
    addWorksheet(wb, "Clean Data")
    addWorksheet(wb, "Normalized")
    addWorksheet(wb, "Unique Filtered")
    addWorksheet(wb, "DE Summary")
    
    progress$set(detail = "Writing main data...", value = 0.1)
    writeData(wb, "Raw Data", all_res$raw)
    writeData(wb, "Clean Data", all_res$clean)
    writeData(wb, "Normalized", all_res$norm)
    writeData(wb, "Unique Filtered", all_res$filtered)
    
    results_list <- all_res$results
    stats_list <- lapply(selected_comps, function(nm) {
      res <- results_list[[nm]]
      if (is.null(res)) return(data.frame(Comparison = nm, Up = 0, Down = 0, Increase = 0, Decrease = 0, NS = 0, Total = 0))
      cnt <- attr(res$data, "counts")
      data.frame(Comparison = nm, Up = cnt$Up, Down = cnt$Down, Increase = cnt$Increase, Decrease = cnt$Decrease, NS = cnt$NS,
                 Total = cnt$Up + cnt$Down + cnt$Increase + cnt$Decrease, stringsAsFactors = FALSE)
    })
    stats_df <- do.call(rbind, stats_list)
    writeData(wb, "DE Summary", stats_df)
    
    # 样式设置（安全处理）
    all_norm_cols <- grep("^Norm_LFQ intensity", colnames(all_res$norm), value = TRUE)
    all_groups <- unique(sapply(all_norm_cols, extract_group_name))
    group_colors <- get_group_colors(all_groups)
    message("[DEBUG] export_excel: groups for coloring: ", paste(names(group_colors), collapse = ", "))
    
    style_red <- createStyle(fontColour = "#C00000")
    style_bold <- createStyle(textDecoration = "bold")
    style_note <- createStyle(fontColour = "#2c3e50", wrapText = TRUE)
    
    for (sheet_name in c("Normalized", "Unique Filtered")) {
      df_sheet <- if (sheet_name == "Normalized") all_res$norm else all_res$filtered
      if (is.null(df_sheet) || nrow(df_sheet) == 0) next
      norm_cols <- grep("^Norm_LFQ intensity", colnames(df_sheet), value = TRUE)
      if (length(group_colors) > 0) {
        for (col_name in norm_cols) {
          col_idx <- which(colnames(df_sheet) == col_name)
          if (length(col_idx) == 0) next
          group_name <- extract_group_name(col_name)
          if (group_name %in% names(group_colors)) {
            addStyle(wb, sheet_name, createStyle(fgFill = group_colors[[group_name]], halign = "center"),
                     rows = 2:(nrow(df_sheet)+1), cols = col_idx, gridExpand = TRUE, stack = TRUE)
          }
        }
      }
      if ("Master protein IDs" %in% colnames(df_sheet)) {
        addStyle(wb, sheet_name, style_red, rows = 2:(nrow(df_sheet)+1),
                 cols = which(colnames(df_sheet) == "Master protein IDs"), gridExpand = TRUE)
      }
      if (!is.null(all_res$unique_col) && all_res$unique_col %in% colnames(df_sheet)) {
        addStyle(wb, sheet_name, style_red, rows = 2:(nrow(df_sheet)+1),
                 cols = which(colnames(df_sheet) == all_res$unique_col), gridExpand = TRUE)
      }
      freezePane(wb, sheet_name, firstRow = TRUE)
    }
    
    for (sn in c("Raw Data", "Clean Data")) {
      df_sn <- switch(sn, "Raw Data" = all_res$raw, "Clean Data" = all_res$clean)
      if (!is.null(df_sn) && "Master protein IDs" %in% colnames(df_sn)) {
        addStyle(wb, sn, style_red, rows = 2:(nrow(df_sn)+1),
                 cols = which(colnames(df_sn) == "Master protein IDs"), gridExpand = TRUE)
      }
      freezePane(wb, sn, firstRow = TRUE)
    }
    
    cols_map <- color_mapping_vector()
    style_up <- createStyle(fontColour = cols_map["Up"])
    style_down <- createStyle(fontColour = cols_map["Down"])
    style_inc <- createStyle(fontColour = cols_map["Increase"])
    style_dec <- createStyle(fontColour = cols_map["Decrease"])
    style_ns <- createStyle(fontColour = cols_map["NS"])
    
    n_comps <- length(selected_comps)
    for (i in seq_along(selected_comps)) {
      comp_name <- selected_comps[i]
      sheet_name <- sanitize_sheet(comp_name)
      addWorksheet(wb, sheet_name)
      res_item <- results_list[[comp_name]]
      if (is.null(res_item)) {
        writeData(wb, sheet_name, data.frame(Message = paste("Analysis for", comp_name, "not available.")))
        next
      }
      df_diff <- res_item$data
      writeData(wb, sheet_name, df_diff)
      
      # 差异表也应用组颜色（仅对 Norm_ 列）
      norm_cols_diff <- grep("^Norm_LFQ intensity", colnames(df_diff), value = TRUE)
      if (length(group_colors) > 0 && length(norm_cols_diff) > 0) {
        for (col_name in norm_cols_diff) {
          col_idx <- which(colnames(df_diff) == col_name)
          if (length(col_idx) == 0) next
          group_name <- extract_group_name(col_name)
          if (group_name %in% names(group_colors)) {
            addStyle(wb, sheet_name, createStyle(fgFill = group_colors[[group_name]], halign = "center"),
                     rows = 2:(nrow(df_diff)+1), cols = col_idx, gridExpand = TRUE, stack = TRUE)
          }
        }
      }
      target_cols_red <- c("Master protein IDs", "n_treat", "n_control", "Unique peptides")
      for (tcol in target_cols_red) {
        if (tcol %in% colnames(df_diff)) {
          addStyle(wb, sheet_name, style_red, rows = 2:(nrow(df_diff)+1),
                   cols = which(colnames(df_diff) == tcol), gridExpand = TRUE)
        }
      }
      if ("regulation" %in% colnames(df_diff)) {
        reg_col <- which(colnames(df_diff) == "regulation")
        addStyle(wb, sheet_name, style_bold, rows = 1, cols = reg_col)
        for (j in seq_len(nrow(df_diff))) {
          val <- df_diff$regulation[j]
          if (!is.na(val) && val %in% names(cols_map)) {
            style <- switch(val, Up = style_up, Down = style_down, Increase = style_inc, Decrease = style_dec, NS = style_ns)
            addStyle(wb, sheet_name, style, rows = j+1, cols = reg_col)
          }
        }
      }
      if ("regulation_note" %in% colnames(df_diff)) {
        note_col_idx <- which(colnames(df_diff) == "regulation_note")
        addStyle(wb, sheet_name, style_note, rows = 2:(nrow(df_diff)+1),
                 cols = note_col_idx, gridExpand = TRUE)
        setColWidths(wb, sheet_name, cols = note_col_idx, widths = 50)
      }
      freezePane(wb, sheet_name, firstRow = TRUE)
      progress$inc(0.4 / n_comps, detail = paste0("Writing ", comp_name, "..."))
    }
    
    progress$set(detail = "Saving...", value = 0.95)
    saveWorkbook(wb, f, overwrite = TRUE)
    progress$set(detail = "Complete!", value = 1.0)
    message("[DEBUG] export_excel: report saved to ", f)
  }
)

# ---------- PDF 导出 ----------
output$download_pdf_report <- downloadHandler(
  filename = function() paste0("Proteomics_Analysis_Report_", Sys.Date(), ".pdf"),
  content = function(file) {
    progress <- shiny::Progress$new()
    progress$set(message = "Generating PDF Report...", value = 0)
    on.exit(progress$close())
    
    req(all_analysis_results(), input$export_comparisons)
    all_res <- all_analysis_results()
    selected_comps <- input$export_comparisons
    if (length(selected_comps) == 0) {
      showNotification("Please select at least one comparison.", type = "error")
      return()
    }
    
    results_list <- all_res$results[selected_comps]
    fcu <- input$fc_up; fcd <- input$fc_down; pc <- as.numeric(input$p_cut)
    cols <- color_mapping_vector()
    point_size <- input$point_size
    
    custom_title <- if (is.null(input$combined_plot_title) || trimws(input$combined_plot_title) == "") {
      "Combined Volcano Plots"
    } else {
      trimws(input$combined_plot_title)
    }
    default_display <- sapply(results_list, function(r) r$name)
    sub_titles <- sapply(seq_along(results_list), function(i) get_subplot_title(i, default_display[i]))
    
    final_plot <- build_combined_plot(results_list, fcu, fcd, pc, cols, custom_title, sub_titles, point_size)
    
    stats_list <- lapply(selected_comps, function(nm) {
      res <- results_list[[nm]]
      cnt <- attr(res$data, "counts")
      data.frame(Comparison = nm, Up = cnt$Up, Down = cnt$Down, Increase = cnt$Increase, Decrease = cnt$Decrease, NS = cnt$NS,
                 Total = cnt$Up + cnt$Down + cnt$Increase + cnt$Decrease, stringsAsFactors = FALSE)
    })
    stats_df <- do.call(rbind, stats_list)
    
    param_text <- paste0(
      "Analysis Parameters:\n",
      "Fold Change up > ", fcu, "\n",
      "Fold Change down < ", fcd, "\n",
      "P-value threshold: ", pc, "\n",
      "Treatment min valid replicates: ", input$min_treat_valid, "\n",
      "Control min valid replicates: ", input$min_ctrl_valid, "\n",
      "Min replicates for t-test: ", input$min_rep_ttest, "\n",
      "Min replicates for Increase: ", input$min_rep_inc, "\n",
      "Min replicates for Decrease: ", input$min_rep_dec, "\n",
      "Min Unique Peptides filter: ", input$min_unique_pep
    )
    
    pdf(file, width = 10, height = 8, onefile = TRUE)
    if (!is.null(final_plot)) {
      print(final_plot)
    } else {
      grid.text("No comparisons available.", x = 0.5, y = 0.5)
    }
    
    grid.newpage()
    vp_table <- viewport(x = 0.5, y = 0.6, width = 0.9, height = 0.6, just = "center")
    pushViewport(vp_table)
    grid.draw(tableGrob(stats_df, rows = NULL, theme = ttheme_minimal()))
    popViewport()
    
    vp_text <- viewport(x = 0.5, y = 0.15, width = 0.9, height = 0.3, just = "center")
    pushViewport(vp_text)
    grid.draw(textGrob(param_text, x = 0, y = 1, just = c("left", "top"), gp = gpar(fontsize = 9, lineheight = 1.2)))
    popViewport()
    
    dev.off()
    progress$set(detail = "PDF report created.", value = 1.0)
  }
)

observeEvent(input$reset_plot_size, {
  updateTextInput(session, "plot_width", value = "10")
  updateTextInput(session, "plot_height", value = "8")
  shinyjs::runjs("$('#plot_width_warning').text(''); $('#plot_height_warning').text('');")
  shinyjs::runjs("$('#plot_width').removeClass('is-invalid'); $('#plot_height').removeClass('is-invalid');")
  showNotification("Plot size reset to default (10x8 inches)", type = "message", duration = 2)
})

observeEvent(input$goto_comparisons, {
  updateNavbarPage(session, "main_navbar", selected = "plots")
})

message("[DEBUG] export_server.R loaded successfully.")