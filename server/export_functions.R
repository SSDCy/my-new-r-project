# server/export_functions.R

# ---------- 辅助函数：标准化导出数据框的列名并移除不需要的列 ----------
normalize_export_data <- function(df) {
  if (is.null(df)) return(NULL)
  colnames(df) <- gsub("\\.", " ", colnames(df))
  colnames(df) <- gsub("\\s+", " ", colnames(df))
  if ("Fasta headers" %in% colnames(df)) {
    df <- df[, colnames(df) != "Fasta headers", drop = FALSE]
  }
  df
}

# 组合图核心绘制
plot_volcano_core_combined <- function(df, fc_up, fc_down, p_cut, cols, point_size = 1.8) {
  req(nrow(df) > 0)
  cnt <- list(Up = sum(df$regulation == "Up", na.rm = TRUE),
              Down = sum(df$regulation == "Down", na.rm = TRUE),
              Increase = sum(df$regulation == "Increase", na.rm = TRUE),
              Decrease = sum(df$regulation == "Decrease", na.rm = TRUE))
  x1 <- -8.5; x2 <- 8.5
  y_max <- max(df$log10P, na.rm = TRUE, 2) * 1.15
  if (!is.finite(y_max) || y_max <= 1) y_max <- 5
  p <- ggplot(df, aes(log2FC, log10P, color = regulation)) +
    geom_point(size = point_size, alpha = 0.6) +
    scale_color_manual(values = cols) +
    geom_vline(xintercept = log2(c(fc_down, fc_up)), lty = 2, color = "gray40", linewidth = 0.4) +
    geom_hline(yintercept = -log10(p_cut), lty = 2, color = "gray40", linewidth = 0.4) +
    coord_cartesian(xlim = c(x1, x2), ylim = c(0, y_max)) +
    volcano_theme() +
    theme(legend.position = "none", plot.title = element_blank(), axis.title = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
          panel.border = element_blank(), plot.margin = margin(5, 5, 5, 5))
  p <- p + annotate("text", x = -4.5, y = y_max*0.95, label = sprintf("%d", cnt$Down), color = cols$Down, fontface = "bold", size = 5, hjust = 1) +
    annotate("text", x = -4.5 + 0.1, y = y_max*0.95, label = sprintf("(%d)", cnt$Decrease), color = cols$Decrease, fontface = "bold", size = 5, hjust = 0) +
    annotate("text", x = 4.5, y = y_max*0.95, label = sprintf("%d", cnt$Up), color = cols$Up, fontface = "bold", size = 5, hjust = 1) +
    annotate("text", x = 4.5 + 0.1, y = y_max*0.95, label = sprintf("(%d)", cnt$Increase), color = cols$Increase, fontface = "bold", size = 5, hjust = 0)
  return(p)
}

build_combined_plot <- function(results_list, fcu, fcd, pc, cols, main_title, sub_titles = NULL, point_size = 1.8) {
  if (length(results_list) == 0) return(NULL)
  n_plots <- length(results_list)
  layout_info <- get_optimal_layout(n_plots)
  ncol <- layout_info$ncol; nrow <- layout_info$nrow
  
  volcano_plots <- lapply(results_list, function(r) plot_volcano_core_combined(r$data, fcu, fcd, pc, cols, point_size))
  
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
    if (length(col_units) > 1) row_grob <- plot_grid(plotlist = col_units, ncol = length(col_units), align = "hv", axis = "none")
    else row_grob <- col_units[[1]]
    row_plots[[r]] <- row_grob
    row_heights[r] <- if (r == 1 || r == nrow) base_height_per_plot else base_height_per_plot + title_height
  }
  
  grid_height <- sum(row_heights)
  plot_height <- grid_height + title_space + bottom_space + top_title_space + bottom_title_space
  
  grid_plot <- plot_grid(plotlist = row_plots, ncol = 1, rel_heights = row_heights, align = "v", axis = "none")
  
  final_plot <- ggdraw() +
    draw_plot(grid_plot, x = 0.12, y = (bottom_space + bottom_title_space) / plot_height, width = 0.76, height = grid_height / plot_height)
  
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

get_subplot_title <- function(index, default_display_name) {
  title_input_id <- paste0("subplot_title_", index)
  custom_title <- input[[title_input_id]]
  if (is.null(custom_title) || trimws(custom_title) == "") default_display_name else trimws(custom_title)
}

output$subplot_titles_ui <- renderUI({
  comps <- sorted_comps()
  if (length(comps) == 0) return(NULL)
  title_inputs <- lapply(seq_along(comps), function(i) {
    comp <- comps[[i]]
    textInputMax(paste0("subplot_title_", i), label = paste0("Sub-plot ", i, " (", comp$name, ")"),
                 value = "", placeholder = paste0("Default: ", comp$name), maxlength = 25,
                 allowed_pattern = "[^a-zA-Z0-9 _-]")
  })
  do.call(tagList, title_inputs)
})

# 子标题重复检查
observe({
  all_ids <- grep("^subplot_title_", names(input), value = TRUE)
  if (length(all_ids) == 0) return()
  if (!is.null(rv$pending_duplicate)) return()
  
  current_vals <- sapply(all_ids, function(id) input[[id]])
  for (id in all_ids) {
    old_val <- isolate(subplot_old_values[[id]])
    if (is.null(old_val)) old_val <- ""
    new_val <- current_vals[[id]]
    if (!identical(new_val, old_val)) {
      other_ids <- setdiff(all_ids, id)
      other_vals <- current_vals[other_ids]
      if (new_val != "" && new_val %in% other_vals) {
        rv$pending_duplicate <- list(id = id, new_val = new_val, old_val = old_val)
        showModal(modalDialog(
          title = "Warning",
          paste0("A sub-plot titled \"", new_val, "\" already exists. Do you wish to continue?"),
          footer = tagList(
            actionButton("confirm_duplicate_title", "Continue", class = "btn btn-secondary"),
            actionButton("cancel_duplicate_title", "Cancel", class = "btn btn-primary")
          ),
          easyClose = FALSE
        ))
        break
      } else {
        subplot_old_values[[id]] <- new_val
      }
    }
  }
})

observeEvent(input$confirm_duplicate_title, {
  req(rv$pending_duplicate)
  info <- rv$pending_duplicate
  subplot_old_values[[info$id]] <- info$new_val
  rv$pending_duplicate <- NULL
  removeModal()
})

observeEvent(input$cancel_duplicate_title, {
  req(rv$pending_duplicate)
  info <- rv$pending_duplicate
  updateTextInput(session, info$id, value = info$old_val)
  subplot_old_values[[info$id]] <- info$old_val
  rv$pending_duplicate <- NULL
  removeModal()
})

# 单图下载
output$download_plot <- downloadHandler(
  filename = function() {
    s <- gsub("[^a-zA-Z0-9_]", "_", input$selected_comparison %||% "plot")
    paste0("Volcano_", s, "_", Sys.Date(), ".", input$plot_format)
  },
  content = function(f) {
    res <- current_result(); req(res); df <- res$data; comp_name <- res$name; cols <- color_mapping()
    if (nrow(df) == 0) { png(f, 800, 600); plot(0,0,type="n"); text(0.5,0.5,"No data"); dev.off(); return() }
    cnt <- attr(df, "counts")
    y_max <- max(df$log10P, na.rm = TRUE, 5) * 1.15
    x1 <- -8.5; x2 <- 8.5
    download_title <- if (!is.null(input$download_single_title) && input$download_single_title != "") input$download_single_title else comp_name
    if (is.null(download_title) || download_title == "") download_title <- comp_name
    y_annot <- y_max * 0.95
    width_val <- as.integer(input$plot_width)
    height_val <- as.integer(input$plot_height)
    if (is.na(width_val) || width_val < 5 || width_val > 30) width_val <- 10
    if (is.na(height_val) || height_val < 5 || height_val > 30) height_val <- 8
    point_size <- input$point_size
    p <- ggplot(df, aes(log2FC, log10P, color = regulation)) + geom_point(size = point_size, alpha = 0.6) +
      scale_color_manual(values = cols) +
      geom_vline(xintercept = log2(c(input$fc_down, input$fc_up)), lty = 2, color = "gray40") +
      geom_hline(yintercept = -log10(as.numeric(input$p_cut)), lty = 2, color = "gray40") +
      coord_cartesian(xlim = c(x1, x2), ylim = c(0, y_max)) +
      labs(title = download_title, x = expression(Log[2]~"(Fold Change)"), y = expression(-Log[10]~"(P-Value)")) +
      volcano_theme() +
      theme(axis.ticks = element_blank(), axis.text = element_blank()) +
      guides(color = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 1))
    p <- p + annotate("text", x = -4.5, y = y_annot, label = sprintf("%d", cnt$Down), color = cols$Down, fontface = "bold", size = 5, hjust = 1) +
      annotate("text", x = -4.5 + 0.1, y = y_annot, label = sprintf("(%d)", cnt$Decrease), color = cols$Decrease, fontface = "bold", size = 5, hjust = 0) +
      annotate("text", x = 4.5, y = y_annot, label = sprintf("%d", cnt$Up), color = cols$Up, fontface = "bold", size = 5, hjust = 1) +
      annotate("text", x = 4.5 + 0.1, y = y_annot, label = sprintf("(%d)", cnt$Increase), color = cols$Increase, fontface = "bold", size = 5, hjust = 0)
    ggsave(f, plot = p, width = width_val, height = height_val, dpi = 300, device = input$plot_format)
  }
)

# 组合图下载
output$download_combined_plot <- downloadHandler(
  filename = function() {
    title_clean <- gsub("[^a-zA-Z0-9_]", "_", input$combined_plot_title)
    if (title_clean == "" || is.null(title_clean)) title_clean <- "Combined_Volcano_Plots"
    paste0(title_clean, "_", Sys.Date(), ".png")
  },
  content = function(f) {
    progress <- shiny::Progress$new(); on.exit(progress$close())
    progress$set(message = "Generating Combined Volcano Plot", value = 0)
    progress$set(detail = "Loading comparison data...", value = 0.1)
    all_res <- all_analysis_results()
    if (is.null(all_res)) { png(f, 800, 600); plot(0,0,type="n"); text(0.5,0.5,"No analysis results"); dev.off(); return() }
    results_list <- all_res$results
    if (length(results_list) == 0) {
      png(f, 800, 600); plot(0,0,type="n"); text(0.5,0.5,"No comparisons available"); dev.off()
      progress$set(detail = "Complete!", value = 1.0); return()
    }
    fcu <- input$fc_up; fcd <- input$fc_down; pc <- as.numeric(input$p_cut)
    cols <- color_mapping()
    custom_title <- if (is.null(input$combined_plot_title) || trimws(input$combined_plot_title) == "") "Combined Volcano Plots" else trimws(input$combined_plot_title)
    comp_names <- names(results_list)
    default_display <- sapply(results_list, function(r) r$name)
    sub_titles <- sapply(seq_along(results_list), function(i) get_subplot_title(i, default_display[i]))
    point_size <- input$point_size
    progress$set(detail = "Building combined plot...", value = 0.5)
    final_plot <- build_combined_plot(results_list, fcu, fcd, pc, cols, custom_title, sub_titles, point_size)
    if (is.null(final_plot)) { png(f, 800, 600); plot(0,0,type="n"); text(0.5,0.5,"No plots to combine"); dev.off(); return() }
    progress$set(detail = "Saving PNG file...", value = 0.9)
    n_plots <- length(results_list)
    layout_info <- get_optimal_layout(n_plots); ncol <- layout_info$ncol; nrow <- layout_info$nrow
    base_height_per_plot <- 4.0; title_height <- 0.5; title_space <- 1.0; bottom_space <- 0.8
    top_title_space <- 0.5; bottom_title_space <- 0.5
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

# Excel 报告下载
output$download_excel <- downloadHandler(
  filename = function() paste0("Proteomics_Report_", Sys.Date(), ".xlsx"),
  content = function(f) {
    progress <- shiny::Progress$new(); on.exit(progress$close())
    progress$set(message = "Generating Excel report...", value = 0)
    
    all_res <- all_analysis_results()
    if (is.null(all_res) || is.null(all_res$raw)) {
      showNotification("No analysis results available. Please complete comparisons.", type = "error")
      return()
    }
    req(all_res, input$export_comparisons)
    selected_comps <- input$export_comparisons
    if (length(selected_comps) == 0) {
      showNotification("Please select at least one comparison.", type = "error"); return()
    }
    if (length(selected_comps) > 10) {
      showNotification("Excel export is limited to 10 comparisons. Only the first 10 will be used.", type = "warning")
      selected_comps <- selected_comps[1:10]
    }
    
    progress$set(detail = "Creating workbook...", value = 0.05)
    wb <- createWorkbook()
    sanitize_sheet <- function(name) sanitize_name(name, max_len = 31)
    
    addWorksheet(wb, "Raw Data")
    addWorksheet(wb, "Clean Data")
    addWorksheet(wb, "Normalized")
    addWorksheet(wb, "Unique Filtered")
    addWorksheet(wb, "DE Summary")
    
    raw_df <- normalize_export_data(all_res$raw)
    clean_df <- normalize_export_data(all_res$clean)
    norm_df <- normalize_export_data(all_res$norm)
    filtered_df <- normalize_export_data(all_res$filtered)
    
    progress$set(detail = "Writing main data...", value = 0.1)
    writeData(wb, "Raw Data", raw_df)
    writeData(wb, "Clean Data", clean_df)
    writeData(wb, "Normalized", norm_df)
    writeData(wb, "Unique Filtered", filtered_df)
    
    results_list <- all_res$results[selected_comps]
    stats_list <- lapply(selected_comps, function(nm) {
      res <- results_list[[nm]]
      cnt <- attr(res$data, "counts")
      data.frame(Comparison = nm, Up = cnt$Up, Down = cnt$Down, Increase = cnt$Increase, Decrease = cnt$Decrease, NS = cnt$NS,
                 Total = cnt$Up + cnt$Down + cnt$Increase + cnt$Decrease, stringsAsFactors = FALSE)
    })
    stats_df <- do.call(rbind, stats_list)
    writeData(wb, "DE Summary", stats_df)
    
    progress$set(detail = "Applying styles...", value = 0.3)
    all_norm_cols <- grep("^Norm_LFQ intensity", colnames(norm_df), value = TRUE)
    all_groups <- unique(sapply(all_norm_cols, extract_group_name))
    group_colors <- get_group_colors(all_groups)
    style_red <- createStyle(fontColour = "#C00000")
    style_bold <- createStyle(textDecoration = "bold")
    style_note <- createStyle(fontColour = "#2c3e50", wrapText = TRUE)
    
    unique_col <- grep("^Unique peptides$", colnames(filtered_df), value = TRUE)[1]
    
    for (sheet_name in c("Normalized", "Unique Filtered")) {
      df_sheet <- if (sheet_name == "Normalized") norm_df else filtered_df
      norm_cols <- grep("^Norm_LFQ intensity", colnames(df_sheet), value = TRUE)
      for (col_name in norm_cols) {
        col_idx <- which(colnames(df_sheet) == col_name)
        group_name <- extract_group_name(col_name)
        if (group_name %in% names(group_colors)) {
          addStyle(wb, sheet_name, createStyle(fgFill = group_colors[[group_name]], halign = "center"),
                   rows = 2:(nrow(df_sheet)+1), cols = col_idx, gridExpand = TRUE, stack = TRUE)
        }
      }
      if ("Master protein IDs" %in% colnames(df_sheet)) {
        addStyle(wb, sheet_name, style_red, rows = 2:(nrow(df_sheet)+1),
                 cols = which(colnames(df_sheet) == "Master protein IDs"), gridExpand = TRUE)
      }
      if (!is.null(unique_col) && unique_col %in% colnames(df_sheet)) {
        addStyle(wb, sheet_name, style_red, rows = 2:(nrow(df_sheet)+1),
                 cols = which(colnames(df_sheet) == unique_col), gridExpand = TRUE)
      }
      freezePane(wb, sheet_name, firstRow = TRUE)
    }
    for (sn in c("Raw Data", "Clean Data")) {
      df_sn <- if (sn == "Raw Data") raw_df else clean_df
      if ("Master protein IDs" %in% colnames(df_sn)) {
        addStyle(wb, sn, style_red, rows = 2:(nrow(df_sn)+1),
                 cols = which(colnames(df_sn) == "Master protein IDs"), gridExpand = TRUE)
      }
      freezePane(wb, sn, firstRow = TRUE)
    }
    
    progress$set(detail = "Writing comparison sheets...", value = 0.5)
    cols_map <- color_mapping()
    style_up <- createStyle(fontColour = cols_map$Up)
    style_down <- createStyle(fontColour = cols_map$Down)
    style_inc <- createStyle(fontColour = cols_map$Increase)
    style_dec <- createStyle(fontColour = cols_map$Decrease)
    style_ns <- createStyle(fontColour = cols_map$NS)
    
    n_comps <- length(selected_comps)
    for (i in seq_along(selected_comps)) {
      comp_name <- selected_comps[i]
      sheet_name <- sanitize_sheet(comp_name)
      counter <- 1
      while (sheet_name %in% sheets(wb)) {
        sheet_name <- sanitize_sheet(paste0(comp_name, "_", counter))
        counter <- counter + 1
      }
      addWorksheet(wb, sheet_name)
      df_diff <- results_list[[comp_name]]$data
      df_diff <- normalize_export_data(df_diff)
      writeData(wb, sheet_name, df_diff)
      
      norm_cols_diff <- grep("^Norm_LFQ intensity", colnames(df_diff), value = TRUE)
      for (col_name in norm_cols_diff) {
        col_idx <- which(colnames(df_diff) == col_name)
        group_name <- extract_group_name(col_name)
        if (group_name %in% names(group_colors)) {
          addStyle(wb, sheet_name, createStyle(fgFill = group_colors[[group_name]], halign = "center"),
                   rows = 2:(nrow(df_diff)+1), cols = col_idx, gridExpand = TRUE, stack = TRUE)
        }
      }
      target_cols_red <- c("Master protein IDs", "n_treat", "n_control")
      if (!is.null(unique_col)) target_cols_red <- c(target_cols_red, unique_col)
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
  }
)

# PDF 报告下载
output$download_pdf_report <- downloadHandler(
  filename = function() paste0("Proteomics_Analysis_Report_", Sys.Date(), ".pdf"),
  content = function(file) {
    progress <- shiny::Progress$new(); on.exit(progress$close())
    progress$set(message = "Generating PDF Report...", value = 0)
    
    all_res <- all_analysis_results()
    if (is.null(all_res) || is.null(all_res$raw)) {
      showNotification("No analysis results available. Please complete comparisons.", type = "error")
      return()
    }
    req(all_res, input$export_comparisons)
    selected_comps <- input$export_comparisons
    if (length(selected_comps) == 0) {
      showNotification("Please select at least one comparison.", type = "error"); return()
    }
    
    results_list <- all_res$results[selected_comps]
    fcu <- input$fc_up; fcd <- input$fc_down; pc <- as.numeric(input$p_cut)
    cols <- color_mapping(); point_size <- input$point_size
    custom_title <- if (is.null(input$combined_plot_title) || trimws(input$combined_plot_title) == "") "Combined Volcano Plots" else trimws(input$combined_plot_title)
    default_display <- sapply(results_list, function(r) r$name)
    sub_titles <- sapply(seq_along(results_list), function(i) get_subplot_title(i, default_display[i]))
    
    progress$set(detail = "Building combined plot...", value = 0.3)
    final_plot <- build_combined_plot(results_list, fcu, fcd, pc, cols, custom_title, sub_titles, point_size)
    
    stats_list <- lapply(selected_comps, function(nm) {
      res <- results_list[[nm]]; cnt <- attr(res$data, "counts")
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
    
    progress$set(detail = "Writing PDF...", value = 0.6)
    pdf(file, width = 10, height = 8, onefile = TRUE)
    
    if (!is.null(final_plot)) { print(final_plot) } else { grid.text("No comparisons available.", x = 0.5, y = 0.5) }
    
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