# server/volcano_plot.R

color_mapping <- reactive(list(
  Up = input$color_up,
  Down = input$color_down,
  Increase = input$color_increase,
  Decrease = input$color_decrease,
  NS = input$color_ns
))

current_result <- reactive({
  req(all_analysis_results(), input$selected_comparison)
  all_res <- all_analysis_results()
  all_res$results[[input$selected_comparison]]
})

volcano_plot_reactive <- reactive({
  res <- current_result()
  req(res)
  df <- res$data
  comp_name <- res$name
  cols <- color_mapping()
  
  if (nrow(df) == 0) {
    return(plotly::plot_ly() %>%
             add_text(x = 0.5, y = 0.5, text = "No data") %>%
             layout(xaxis = list(showticklabels = FALSE), yaxis = list(showticklabels = FALSE)))
  }
  
  if (!"Master protein IDs" %in% colnames(df)) {
    df$`Master protein IDs` <- "Unknown"
  }
  df$`Master protein IDs` <- as.character(df$`Master protein IDs`)
  
  cnt <- attr(df, "counts")
  plot_title <- if (!is.null(input$single_plot_title) && nzchar(input$single_plot_title))
    input$single_plot_title else comp_name
  
  x1 <- if (any(df$regulation %in% c("Increase", "Decrease"))) -8.5 else -8
  x2 <- if (any(df$regulation %in% c("Increase", "Decrease"))) 8.5 else 8
  y_max <- max(df$log10P, na.rm = TRUE, 5) * 1.15
  point_size <- input$point_size
  
  p <- ggplot(df, aes(log2FC, log10P, color = regulation,
                      customdata = `Master protein IDs`,
                      text = paste0("Protein: ", `Master protein IDs`,
                                    "\nlog2FC: ", round(log2FC, 3),
                                    "\n-log10P: ", round(log10P, 3)))) +
    geom_point(size = point_size, alpha = 0.6) +
    scale_color_manual(values = cols) +
    geom_vline(xintercept = log2(c(input$fc_down, input$fc_up)), lty = 2, color = "gray40") +
    geom_hline(yintercept = -log10(as.numeric(input$p_cut)), lty = 2, color = "gray40") +
    coord_cartesian(xlim = c(x1, x2), ylim = c(0, y_max)) +
    labs(title = plot_title, x = expression(Log[2]~FC), y = expression(-Log[10]~P)) +
    volcano_theme() +
    guides(color = guide_legend(override.aes = list(size = 4), ncol = 1))
  
  ggplotly(p, tooltip = "text", source = "volcano_plot") %>%
    event_register("plotly_click") %>%
    config(displayModeBar = FALSE) %>%
    layout(showlegend = TRUE, margin = list(t = 80, b = 50, l = 60, r = 60))
})

output$volcano_plot <- renderPlotly({
  tryCatch(
    volcano_plot_reactive(),
    error = function(e) {
      plotly::plot_ly() %>%
        add_text(x = 0.5, y = 0.5, text = paste("Error:\n", e$message)) %>%
        layout(xaxis = list(showticklabels = FALSE), yaxis = list(showticklabels = FALSE))
    }
  )
})

observeEvent(event_data("plotly_click", source = "volcano_plot"), {
  cd <- event_data("plotly_click", source = "volcano_plot")$customdata
  if (!is.null(cd) && cd != "") {
    if (!is.null(rv$current_profile_protein)) removeModal()
    rv$current_profile_protein <- cd
    clicked_protein(list(id = cd, ts = Sys.time()))
  }
})

observeEvent(clicked_protein(), {
  req(clicked_protein())
  pro_id <- clicked_protein()$id
  showModal(modalDialog(
    title = div(icon("chart-line"), paste("Expression Profile:", pro_id)),
    size = "l",
    easyClose = TRUE,
    footer = modalButton("Close"),
    plotlyOutput("protein_profile_plot", height = "450px")
  ))
  observeEvent(input$`modalButton`, { rv$current_profile_protein <- NULL }, once = TRUE)
})

output$protein_profile_plot <- renderPlotly({
  req(clicked_protein())
  pro_id <- clicked_protein()$id
  nd <- norm_data_full()
  if (is.null(nd)) {
    return(plotly_empty() %>% layout(
      title = "Normalized data not available. Please ensure you have selected a baseline sample and run preprocessing if needed."
    ))
  }
  idx <- which(nd[["Master protein IDs"]] == pro_id)
  if (length(idx) == 0) return(plotly_empty() %>% layout(title = paste("Protein", pro_id, "not found in normalized data.")))
  if (length(idx) > 1) idx <- idx[1]
  norm_cols <- grep("^Norm_LFQ intensity ", colnames(nd), value = TRUE)
  intensities <- as.numeric(nd[idx, norm_cols])
  samples <- gsub("^Norm_LFQ intensity ", "", norm_cols)
  all_groups <- rv$groups
  group_map <- setNames(rep(names(all_groups), lengths(all_groups)), unlist(all_groups))
  groups_vec <- group_map[samples]
  groups_vec[is.na(groups_vec)] <- "Unassigned"
  df_plot <- data.frame(
    Sample = as.character(samples),
    Group = as.character(groups_vec),
    Intensity = as.numeric(intensities),
    stringsAsFactors = FALSE
  )
  df_plot$Intensity[is.na(df_plot$Intensity)] <- 0
  df_plot$Group <- factor(df_plot$Group, levels = c(unique(groups_vec[groups_vec != "Unassigned"]), "Unassigned"))
  
  p <- ggplot(df_plot, aes(x = Group, y = Intensity, color = Group,
                           text = paste0("Sample: ", Sample,
                                         "\nGroup: ", Group,
                                         "\nIntensity: ", format(Intensity, scientific = FALSE)))) +
    geom_jitter(width = 0.15, size = 3, alpha = 0.7) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "black", fatten = 2) +
    labs(title = paste("Expression profile of", pro_id), y = "Normalized Intensity", x = "") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggplotly(p, tooltip = "text") %>%
    layout(plot_bgcolor = 'white', paper_bgcolor = 'white', margin = list(b = 80))
})

output$plot_info_ui <- renderUI({
  res <- current_result()
  req(res)
  df <- res$data
  cnt <- attr(df, "counts")
  cols <- color_mapping()
  div(
    fluidRow(
      column(3, div(style = "text-align: center; padding:20px; background:#fdecec; border-radius:12px;",
                    h4(icon("arrow-up"), style = paste0("color:", cols$Up, "; margin:0; font-size:2em;")),
                    h3(cnt$Up, style = paste0("color:", cols$Up, "; margin:0; font-size:2.5em; font-weight:bold;")),
                    p("Up-regulated", style = "margin:0; color:#666; font-size:1.1em;"))),
      column(3, div(style = "text-align: center; padding:20px; background:#ebf5fb; border-radius:12px;",
                    h4(icon("arrow-down"), style = paste0("color:", cols$Down, "; margin:0; font-size:2em;")),
                    h3(cnt$Down, style = paste0("color:", cols$Down, "; margin:0; font-size:2.5em; font-weight:bold;")),
                    p("Down-regulated", style = "margin:0; color:#666; font-size:1.1em;"))),
      column(3, div(style = "background:#fdecec; padding:20px; border-radius:12px; text-align:center;",
                    h4(icon("plus-circle"), style = paste0("color:", cols$Increase, "; margin:0; font-size:2em;")),
                    h3(cnt$Increase, style = paste0("color:", cols$Increase, "; margin:0; font-size:2.5em; font-weight:bold;")),
                    p("Increase", style = "margin:0; color:#666; font-size:1.1em;"))),
      column(3, div(style = "background:#ebf5fb; padding:20px; border-radius:12px; text-align:center;",
                    h4(icon("minus-circle"), style = paste0("color:", cols$Decrease, "; margin:0; font-size:2em;")),
                    h3(cnt$Decrease, style = paste0("color:", cols$Decrease, "; margin:0; font-size:2.5em; font-weight:bold;")),
                    p("Decrease", style = "margin:0; color:#666; font-size:1.1em;")))
    ), br(),
    p(strong("Comparison:"), res$name),
    p(strong("Parameters:"), sprintf("FC > %.2f or < %.2f, P-value < %s", input$fc_up, input$fc_down, input$p_cut))
  )
})