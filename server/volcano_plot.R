# server/volcano_plot.R
message("[DEBUG] volcano_plot.R loading... (safe stat_method and min_unique_pep)")

# ---------- 默认颜色 ----------
default_colors <- reactive(list(Up="#FF0000", Down="#0000FF", Increase="#C00000", Decrease="#0945A5", NS="#7f7e83"))

output$val_up <- renderText({ input$color_up })
output$val_down <- renderText({ input$color_down })
output$val_inc <- renderText({ input$color_increase })
output$val_dec <- renderText({ input$color_decrease })
output$val_ns <- renderText({ input$color_ns })

output$color_preview <- renderUI({
  cols <- list(Up = input$color_up, Down = input$color_down, Increase = input$color_increase, Decrease = input$color_decrease, NS = input$color_ns)
  items <- lapply(names(cols), function(nm) div(class = "color-preview-item", div(class = "color-preview-swatch", style = paste0("background-color:", cols[[nm]], ";")), div(class = "color-preview-label", nm)))
  div(class = "color-preview-container", do.call(tagList, items))
})

observeEvent(input$reset_color, {
  cols <- default_colors()
  colourpicker::updateColourInput(session, "color_up", value = cols$Up)
  colourpicker::updateColourInput(session, "color_down", value = cols$Down)
  colourpicker::updateColourInput(session, "color_increase", value = cols$Increase)
  colourpicker::updateColourInput(session, "color_decrease", value = cols$Decrease)
  colourpicker::updateColourInput(session, "color_ns", value = cols$NS)
  updateNumericInput(session, "point_size", value = 4)
  showNotification("Colors and point size reset", type = "message", duration = 2)
})

color_mapping_vector <- reactive({
  cols <- list(Up = input$color_up, Down = input$color_down, Increase = input$color_increase, 
               Decrease = input$color_decrease, NS = input$color_ns)
  vec <- unlist(cols)
  message("[DEBUG] volcano: color_mapping_vector = ", paste(names(vec), vec, sep = ":", collapse = ", "))
  vec
})

color_mapping <- reactive(list(Up=input$color_up, Down=input$color_down, Increase=input$color_increase, Decrease=input$color_decrease, NS=input$color_ns))

safe_get_analysis_matrix <- function() {
  tryCatch({
    mat <- get_analysis_matrix()
    if (is.null(mat)) { message("[DEBUG] volcano: get_analysis_matrix() returned NULL"); return(NULL) }
    message("[DEBUG] volcano: got analysis matrix, dim = ", nrow(mat), " x ", ncol(mat))
    mat
  }, error = function(e) { message("[ERROR] volcano: failed to get analysis matrix: ", e$message); NULL })
}

selected_volcano_comparison <- reactive({
  comp_name <- input$selected_comparison
  if (is.null(comp_name) || comp_name == "") return(NULL)
  comps <- rv$comparisons
  if (length(comps) == 0) return(NULL)
  comp <- Find(function(c) c$name == comp_name, comps)
  if (is.null(comp)) { message("[DEBUG] volcano: comparison not found: ", comp_name); return(NULL) }
  message("[DEBUG] volcano: selected comparison = ", comp_name, ", treat = ", comp$treat, ", ctrl = ", comp$ctrl)
  comp
})

observe({
  comps <- rv$comparisons
  if (length(comps) > 0) {
    choices <- sapply(comps, `[[`, "name")
    updateSelectInput(session, "selected_comparison", choices = choices, selected = choices[1])
    message("[DEBUG] volcano: updated comparison choices, n = ", length(choices))
  } else {
    updateSelectInput(session, "selected_comparison", choices = character(0))
  }
})

volcano_de_result <- reactive({
  comp <- selected_volcano_comparison(); req(comp)
  nd <- norm_data_full()
  if (is.null(nd)) { message("[DEBUG] volcano: norm_data_full() is NULL"); return(NULL) }
  
  # 应用 unique peptide 过滤（安全处理）
  unique_col <- grep("^Unique peptides$", colnames(nd), value = TRUE)[1]
  if (!is.na(unique_col) && !is.null(input$min_unique_pep) && input$min_unique_pep > 1) {
    nd[[unique_col]] <- as.numeric(nd[[unique_col]])
    nd <- nd[nd[[unique_col]] >= input$min_unique_pep, ]
    message("[DEBUG] volcano: after unique peptide filter (>= ", input$min_unique_pep, "), nrow = ", nrow(nd))
  }
  
  treat_group <- comp$treat; ctrl_group <- comp$ctrl
  treat_samples <- rv$groups[[treat_group]]; ctrl_samples <- rv$groups[[ctrl_group]]
  if (length(treat_samples) == 0 || length(ctrl_samples) == 0) {
    showNotification("One of the groups has no samples.", type = "error"); return(NULL)
  }
  
  treat_cols <- paste0("Norm_LFQ intensity ", treat_samples)
  ctrl_cols <- paste0("Norm_LFQ intensity ", ctrl_samples)
  treat_cols <- treat_cols[treat_cols %in% colnames(nd)]
  ctrl_cols <- ctrl_cols[ctrl_cols %in% colnames(nd)]
  if (length(treat_cols) == 0 || length(ctrl_cols) == 0) {
    showNotification("Sample columns mismatch. Please check group assignments.", type = "error"); return(NULL)
  }
  
  message("[DEBUG] volcano: matched treat = ", length(treat_cols), " columns, ctrl = ", length(ctrl_cols), " columns")
  
  sub_df <- nd[, c("Master protein IDs", treat_cols, ctrl_cols), drop = FALSE]
  
  fc_up <- input$fc_up; fc_down <- input$fc_down; p_cut <- as.numeric(input$p_cut)
  stat_method <- input$stat_method %||% "t-test"   # 安全默认
  message("[DEBUG] volcano: running DE with method=", stat_method, ", FC_up=", fc_up, ", FC_down=", fc_down, ", p_cut=", p_cut)
  
  res <- tryCatch({
    run_de_analysis(data_subset = sub_df, treat_cols = treat_cols,
                    ctrl_cols = ctrl_cols,
                    fc_up = fc_up, fc_down = fc_down, p_cut = p_cut, stat_method = stat_method)
  }, error = function(e) { message("[ERROR] volcano: DE failed: ", e$message); NULL })
  
  if (!is.null(res)) {
    for (col in c("log2FC","log10P","regulation")) {
      if (!(col %in% colnames(res))) res[[col]] <- NA_real_
    }
  }
  if (is.null(res) || nrow(res) == 0) { message("[DEBUG] volcano: DE result empty"); return(NULL) }
  message("[DEBUG] volcano: DE completed, nrow = ", nrow(res), ", Up = ", sum(res$regulation == "Up"),
          ", Down = ", sum(res$regulation == "Down"), ", Increase = ", sum(res$regulation == "Increase"),
          ", Decrease = ", sum(res$regulation == "Decrease"))
  res
})

volcano_counts <- reactive({
  res <- volcano_de_result()
  if (is.null(res)) return(list(Up=0, Down=0, Increase=0, Decrease=0, Total=0))
  list(Up=sum(res$regulation=="Up"), Down=sum(res$regulation=="Down"),
       Increase=sum(res$regulation=="Increase"), Decrease=sum(res$regulation=="Decrease"),
       Total=nrow(res))
})

output$plot_info_ui <- renderUI({
  res <- volcano_de_result(); req(res); cnt <- volcano_counts(); cols <- color_mapping()
  comp_name <- selected_volcano_comparison()$name
  div(
    fluidRow(
      column(3, div(style = "text-align: center; padding:20px; background:#fdecec; border-radius:12px;",
                    h4(icon("arrow-up"), style=paste0("color:", cols$Up, "; margin:0; font-size:2em;")),
                    h3(cnt$Up, style=paste0("color:", cols$Up, "; margin:0; font-size:2.5em; font-weight:bold;")),
                    p("Up-regulated", style="margin:0; color:#666; font-size:1.1em;"))),
      column(3, div(style = "text-align: center; padding:20px; background:#ebf5fb; border-radius:12px;",
                    h4(icon("arrow-down"), style=paste0("color:", cols$Down, "; margin:0; font-size:2em;")),
                    h3(cnt$Down, style=paste0("color:", cols$Down, "; margin:0; font-size:2.5em; font-weight:bold;")),
                    p("Down-regulated", style="margin:0; color:#666; font-size:1.1em;"))),
      column(3, div(style = "background:#fdecec; padding:20px; border-radius:12px; text-align:center;",
                    h4(icon("plus-circle"), style=paste0("color:", cols$Increase, "; margin:0; font-size:2em;")),
                    h3(cnt$Increase, style=paste0("color:", cols$Increase, "; margin:0; font-size:2.5em; font-weight:bold;")),
                    p("Increase", style="margin:0; color:#666; font-size:1.1em;"))),
      column(3, div(style = "background:#ebf5fb; padding:20px; border-radius:12px; text-align:center;",
                    h4(icon("minus-circle"), style=paste0("color:", cols$Decrease, "; margin:0; font-size:2em;")),
                    h3(cnt$Decrease, style=paste0("color:", cols$Decrease, "; margin:0; font-size:2.5em; font-weight:bold;")),
                    p("Decrease", style="margin:0; color:#666; font-size:1.1em;")))
    ), br(),
    p(strong("Comparison:"), comp_name),
    p(strong("Parameters:"), sprintf("FC > %.2f or < %.2f, P-value < %s", input$fc_up, input$fc_down, input$p_cut))
  )
})

output$volcano_plot <- renderPlotly({
  tryCatch({
    res <- volcano_de_result()
    if (is.null(res)) return(plotly::plot_ly() %>% layout(title = "No DE result"))
    df <- res; if (!"Master protein IDs" %in% colnames(df)) df$`Master protein IDs` <- rownames(df)
    color_vec <- color_mapping_vector(); point_size <- input$point_size
    comp_name <- selected_volcano_comparison()$name
    plot_title <- if (!is.null(input$single_plot_title) && input$single_plot_title != "") input$single_plot_title else comp_name
    fc_up <- input$fc_up; fc_down <- input$fc_down; p_cut <- as.numeric(input$p_cut)
    x_up <- log2(fc_up); x_down <- log2(fc_down); y_cut <- -log10(p_cut)
    df$hover_text <- paste("ID:", df$`Master protein IDs`, "<br>log2FC:", round(df$log2FC,3), "<br>-log10(P):", round(df$log10P,3), "<br>Regulation:", df$regulation)
    df$customdata <- df$`Master protein IDs`
    p <- plot_ly(data = df, x = ~log2FC, y = ~log10P, text = ~hover_text, hoverinfo = "text",
                 customdata = ~customdata, color = ~regulation, colors = color_vec,
                 type = "scatter", mode = "markers", marker = list(size = point_size, opacity = 0.6), source = "volcano_plot") %>%
      event_register("plotly_click") %>%
      layout(title = plot_title, xaxis = list(title = "log2(Fold Change)"), yaxis = list(title = "-log10(P-value)"),
             legend = list(title = list(text = "Regulation")),
             shapes = list(
               list(type = "line", x0 = x_up, x1 = x_up, y0 = 0, y1 = 1, yref = "paper", line = list(dash = "dash", color = "gray40")),
               list(type = "line", x0 = x_down, x1 = x_down, y0 = 0, y1 = 1, yref = "paper", line = list(dash = "dash", color = "gray40")),
               list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = y_cut, y1 = y_cut, line = list(dash = "dash", color = "gray40"))
             ))
    message("[DEBUG] volcano: plotly object created with threshold lines")
    p
  }, error = function(e) { plotly::plot_ly() %>% layout(title = paste("Error:", e$message)) })
})

observeEvent(event_data("plotly_click", source = "volcano_plot"), {
  cd <- event_data("plotly_click", source = "volcano_plot")$customdata
  if (!is.null(cd) && cd != "") { message("[DEBUG] volcano: clicked protein = ", cd); clicked_protein(list(id = cd, ts = Sys.time())) }
})

observeEvent(clicked_protein(), {
  req(clicked_protein()); pro_id <- clicked_protein()$id
  showModal(modalDialog(title = div(icon("chart-line"), paste("Expression Profile:", pro_id)), size = "l", easyClose = TRUE, footer = modalButton("Close"), plotlyOutput("protein_profile_plot", height = "450px")))
})

output$protein_profile_plot <- renderPlotly({
  req(clicked_protein()); pro_id <- clicked_protein()$id
  nd <- tryCatch(norm_data_full(), error = function(e) NULL)
  if (is.null(nd)) { nd <- get_analysis_matrix(); if (is.null(nd)) return(plotly::plot_ly() %>% layout(title = "Normalized data not available.")) }
  idx <- which(nd[["Master protein IDs"]] == pro_id)
  if (length(idx) == 0) return(plotly::plot_ly() %>% layout(title = paste("Protein", pro_id, "not found.")))
  if (length(idx) > 1) idx <- idx[1]
  mat <- get_analysis_matrix()
  if (!is.null(mat) && pro_id %in% rownames(mat)) {
    intensities <- as.numeric(mat[pro_id, ]); samples <- extract_sample_names(colnames(mat))
  } else {
    norm_cols <- grep("^Norm_LFQ intensity ", colnames(nd), value = TRUE)
    if (length(norm_cols) == 0) norm_cols <- grep("^LFQ intensity ", colnames(nd), value = TRUE)
    intensities <- as.numeric(nd[idx, norm_cols]); samples <- gsub("^(Norm_)?LFQ intensity ", "", norm_cols)
  }
  all_groups <- rv$groups; group_map <- setNames(rep(names(all_groups), lengths(all_groups)), unlist(all_groups))
  groups_vec <- group_map[samples]; groups_vec[is.na(groups_vec)] <- "Unassigned"
  df_plot <- data.frame(Sample = samples, Group = groups_vec, Intensity = intensities, stringsAsFactors = FALSE)
  df_plot$Intensity[is.na(df_plot$Intensity)] <- 0
  p <- ggplot(df_plot, aes(x = Group, y = Intensity, color = Group, text = paste("Sample:", Sample, "<br>Intensity:", round(Intensity,2)))) +
    geom_jitter(width = 0.15, size = 3, alpha = 0.7) + stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "black", fatten = 2) +
    labs(title = paste("Expression profile of", pro_id), y = "Intensity", x = "") + theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  ggplotly(p, tooltip = "text") %>% layout(plot_bgcolor = 'white', paper_bgcolor = 'white', margin = list(b = 80))
})

# ---------- 预处理步骤指示器（可折叠） ----------
output$volcano_preprocess_steps <- renderUI({
  steps <- list()
  steps <- c(steps, paste0("Missing Value Filter: threshold = ", input$max_missing_fraction %||% 0.5,
                           ", mode = ", preprocessing_params$missing_filter_mode %||% "global"))
  steps <- c(steps, paste0("Minimum Intensity Filter: threshold = ", input$min_intensity,
                           ", min samples = ", input$min_samples_above_intensity %||% 1))
  imp <- preprocessing_params$imputation_method %||% "none"
  if (imp == "none") {
    steps <- c(steps, "Missing Value Imputation: none (some proteins may have missing values)")
  } else {
    steps <- c(steps, paste0("Missing Value Imputation: ", imp))
  }
  if (isTRUE(preprocessing_params$batch_performed)) {
    steps <- c(steps, "Batch Correction (ComBat): applied")
  } else {
    steps <- c(steps, "Batch Correction: not applied")
  }
  steps <- c(steps, "Normalization: Total intensity normalization (baseline sample) + unique peptide filter")
  steps <- c(steps, "Data source: Normalized expression data (Norm_LFQ intensity columns)")
  
  step_tags <- lapply(seq_along(steps), function(i) {
    tagList(
      if (i > 1) tags$span(style = "font-size: 20px; color: #3498db; margin: 0 8px;", "→"),
      tags$span(style = "background: #e8f0fe; padding: 6px 12px; border-radius: 15px; font-size: 13px;", steps[[i]])
    )
  })
  tags$details(
    tags$summary("Data preprocessing steps applied before this volcano plot", style = "cursor: pointer; font-weight: bold; color: #2c3e50; margin-bottom: 10px;"),
    div(style = "display: flex; flex-wrap: wrap; align-items: center;", do.call(tagList, step_tags))
  )
})

# ---------- 下载单图 ----------
output$download_volcano_png <- downloadHandler(
  filename = function() paste0("Volcano_", Sys.Date(), ".png"),
  content = function(file) {
    tryCatch({
      res <- volcano_de_result()
      if (is.null(res)) { png(file); plot.new(); text(0.5,0.5,"No data"); dev.off(); return() }
      df <- res; if (!"Master protein IDs" %in% colnames(df)) df$`Master protein IDs` <- rownames(df)
      required <- c("log2FC", "log10P", "regulation")
      if (!all(required %in% colnames(df))) { png(file); plot.new(); text(0.5,0.5,"Missing columns"); dev.off(); return() }
      cols <- color_mapping_vector()
      point_size <- input$point_size
      comp_name <- selected_volcano_comparison()$name
      download_title <- if (!is.null(input$single_plot_title) && input$single_plot_title != "") input$single_plot_title else comp_name
      cnt <- volcano_counts()
      y_annot <- max(df$log10P, na.rm = TRUE) * 1.05
      if (is.infinite(y_annot) || is.na(y_annot)) y_annot <- 8
      p <- ggplot(df, aes(log2FC, log10P, color = regulation)) +
        geom_point(size = point_size, alpha = 0.6) + scale_color_manual(values = cols) +
        geom_vline(xintercept = log2(c(input$fc_down, input$fc_up)), lty = 2, color = "gray40") +
        geom_hline(yintercept = -log10(as.numeric(input$p_cut)), lty = 2, color = "gray40") +
        labs(title = download_title, x = expression(Log[2]~"(Fold Change)"), y = expression(-Log[10]~"(P-Value)")) +
        volcano_theme() +
        annotate("text", x = -4.5, y = y_annot, label = sprintf("%d", cnt$Down), color = cols["Down"], fontface = "bold", size = 5, hjust = 1) +
        annotate("text", x = -4.5 + 0.1, y = y_annot, label = sprintf("(%d)", cnt$Decrease), color = cols["Decrease"], fontface = "bold", size = 5, hjust = 0) +
        annotate("text", x = 4.5, y = y_annot, label = sprintf("%d", cnt$Up), color = cols["Up"], fontface = "bold", size = 5, hjust = 1) +
        annotate("text", x = 4.5 + 0.1, y = y_annot, label = sprintf("(%d)", cnt$Increase), color = cols["Increase"], fontface = "bold", size = 5, hjust = 0)
      ggsave(file, plot = p, width = 8, height = 6, dpi = 150)
    }, error = function(e) {
      message("[ERROR] download_volcano_png: ", e$message)
      png(file); plot.new(); text(0.5,0.5, paste("Error:", e$message)); dev.off()
    })
  }
)

message("[DEBUG] volcano_plot.R loaded successfully (safe).")