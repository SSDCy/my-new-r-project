# server/volcano_plot.R
# 火山图模块：颜色、点大小、交互、蛋白表达谱、计数卡片、阈值线
message("[DEBUG] volcano_plot.R loading...")

# ---------- 默认颜色 ----------
default_colors <- reactive(list(Up="#FF0000", Down="#0000FF", Increase="#C00000", Decrease="#0945A5", NS="#7f7e83"))

# ---------- 颜色值输出（用于 UI 显示） ----------
output$val_up <- renderText({ input$color_up })
output$val_down <- renderText({ input$color_down })
output$val_inc <- renderText({ input$color_increase })
output$val_dec <- renderText({ input$color_decrease })
output$val_ns <- renderText({ input$color_ns })

# ---------- 颜色预览 ----------
output$color_preview <- renderUI({
  cols <- list(Up = input$color_up, Down = input$color_down, Increase = input$color_increase, Decrease = input$color_decrease, NS = input$color_ns)
  items <- lapply(names(cols), function(nm) div(class = "color-preview-item", div(class = "color-preview-swatch", style = paste0("background-color:", cols[[nm]], ";")), div(class = "color-preview-label", nm)))
  div(class = "color-preview-container", do.call(tagList, items))
})

# ---------- 重置颜色 ----------
observeEvent(input$reset_color, {
  cols <- default_colors()
  colourpicker::updateColourInput(session, "color_up", value = cols$Up)
  colourpicker::updateColourInput(session, "color_down", value = cols$Down)
  colourpicker::updateColourInput(session, "color_increase", value = cols$Increase)
  colourpicker::updateColourInput(session, "color_decrease", value = cols$Decrease)
  colourpicker::updateColourInput(session, "color_ns", value = cols$NS)
  updateNumericInput(session, "point_size", value = 1.8)
  showNotification("Colors and point size reset", type = "message", duration = 2)
})

# ---------- 颜色映射（转换为命名向量供 plotly 使用） ----------
color_mapping_vector <- reactive({
  cols <- list(Up = input$color_up, Down = input$color_down, Increase = input$color_increase, 
               Decrease = input$color_decrease, NS = input$color_ns)
  vec <- unlist(cols)
  message("[DEBUG] volcano: color_mapping_vector = ", paste(names(vec), vec, sep = ":", collapse = ", "))
  vec
})

# 保留原列表形式供 plot_info_ui 使用
color_mapping <- reactive(list(Up=input$color_up, Down=input$color_down, Increase=input$color_increase, Decrease=input$color_decrease, NS=input$color_ns))

# ---------- 分析矩阵（复用预处理数据） ----------
safe_get_analysis_matrix <- function() {
  tryCatch({
    mat <- get_analysis_matrix()
    if (is.null(mat)) {
      message("[DEBUG] volcano: get_analysis_matrix() returned NULL")
      return(NULL)
    }
    message("[DEBUG] volcano: got analysis matrix, dim = ", nrow(mat), " x ", ncol(mat))
    mat
  }, error = function(e) {
    message("[ERROR] volcano: failed to get analysis matrix: ", e$message)
    NULL
  })
}

# ---------- 当前选中的比较 ----------
selected_volcano_comparison <- reactive({
  comp_name <- input$selected_comparison
  if (is.null(comp_name) || comp_name == "") return(NULL)
  comps <- rv$comparisons
  if (length(comps) == 0) return(NULL)
  comp <- Find(function(c) c$name == comp_name, comps)
  if (is.null(comp)) {
    message("[DEBUG] volcano: comparison not found: ", comp_name)
    return(NULL)
  }
  message("[DEBUG] volcano: selected comparison = ", comp_name,
          ", treat = ", comp$treat, ", ctrl = ", comp$ctrl)
  comp
})

# ---------- 更新比较选择器 ----------
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

# ---------- 差异分析结果（针对所选比较） ----------
volcano_de_result <- reactive({
  comp <- selected_volcano_comparison()
  req(comp)
  
  mat <- safe_get_analysis_matrix()
  req(mat)
  
  treat_group <- comp$treat
  ctrl_group <- comp$ctrl
  treat_samples <- rv$groups[[treat_group]]
  ctrl_samples <- rv$groups[[ctrl_group]]
  if (length(treat_samples) == 0 || length(ctrl_samples) == 0) {
    message("[DEBUG] volcano: one group has no samples")
    showNotification("One of the groups has no samples.", type = "error")
    return(NULL)
  }
  
  mat_colnames <- colnames(mat)
  mat_short <- extract_sample_names(mat_colnames)
  
  treat_idx <- which(mat_short %in% treat_samples)
  ctrl_idx <- which(mat_short %in% ctrl_samples)
  
  if (length(treat_idx) == 0 || length(ctrl_idx) == 0) {
    message("[DEBUG] volcano: could not match sample names in matrix")
    showNotification("Sample names mismatch. Please check group assignments.", type = "error")
    return(NULL)
  }
  
  message("[DEBUG] volcano: matched treat = ", length(treat_idx), " samples, ctrl = ", length(ctrl_idx), " samples")
  
  sub_mat <- mat[, c(treat_idx, ctrl_idx), drop = FALSE]
  protein_ids <- rownames(sub_mat)
  sub_df <- as.data.frame(sub_mat)
  sub_df$`Master protein IDs` <- protein_ids
  
  fc_up <- input$fc_up
  fc_down <- input$fc_down
  p_cut <- as.numeric(input$p_cut)
  stat_method <- input$stat_method
  message("[DEBUG] volcano: running DE with method=", stat_method, ", FC_up=", fc_up, ", FC_down=", fc_down, ", p_cut=", p_cut)
  
  res <- tryCatch({
    run_de_analysis(
      data_subset = sub_df,
      treat_cols = colnames(sub_mat)[1:length(treat_idx)],
      ctrl_cols = colnames(sub_mat)[(length(treat_idx)+1):ncol(sub_mat)],
      fc_up = fc_up,
      fc_down = fc_down,
      p_cut = p_cut,
      stat_method = stat_method
    )
  }, error = function(e) {
    message("[ERROR] volcano: DE analysis failed: ", e$message)
    showNotification(paste("DE analysis error:", e$message), type = "error")
    NULL
  })
  
  if (is.null(res) || nrow(res) == 0) {
    message("[DEBUG] volcano: DE result empty")
    return(NULL)
  }
  
  message("[DEBUG] volcano: DE completed, nrow = ", nrow(res), ", Up = ", sum(res$regulation == "Up"),
          ", Down = ", sum(res$regulation == "Down"), ", Increase = ", sum(res$regulation == "Increase"),
          ", Decrease = ", sum(res$regulation == "Decrease"))
  res
})

# ---------- 计数卡片 ----------
volcano_counts <- reactive({
  res <- volcano_de_result()
  if (is.null(res)) return(list(Up = 0, Down = 0, Increase = 0, Decrease = 0, Total = 0))
  list(
    Up = sum(res$regulation == "Up", na.rm = TRUE),
    Down = sum(res$regulation == "Down", na.rm = TRUE),
    Increase = sum(res$regulation == "Increase", na.rm = TRUE),
    Decrease = sum(res$regulation == "Decrease", na.rm = TRUE),
    Total = nrow(res)
  )
})

# ---------- 绘图信息卡片 ----------
output$plot_info_ui <- renderUI({
  res <- volcano_de_result()
  req(res)
  cnt <- volcano_counts()
  cols <- color_mapping()
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

# ---------- 交互式火山图（含阈值线） ----------
output$volcano_plot <- renderPlotly({
  tryCatch({
    res <- volcano_de_result()
    if (is.null(res)) {
      return(plotly::plot_ly() %>% layout(title = "No DE result"))
    }
    df <- res
    if (!"Master protein IDs" %in% colnames(df)) df$`Master protein IDs` <- rownames(df)
    
    color_vec <- color_mapping_vector()
    point_size <- input$point_size
    comp_name <- selected_volcano_comparison()$name
    plot_title <- if (!is.null(input$single_plot_title) && input$single_plot_title != "") input$single_plot_title else comp_name
    
    fc_up <- input$fc_up
    fc_down <- input$fc_down
    p_cut <- as.numeric(input$p_cut)
    x_up <- log2(fc_up)
    x_down <- log2(fc_down)
    y_cut <- -log10(p_cut)
    
    # 计算坐标范围，为阈值线留出空间
    x_range <- range(df$log2FC, na.rm = TRUE)
    y_range <- range(df$log10P, na.rm = TRUE)
    x_margin <- diff(x_range) * 0.1
    y_margin <- diff(y_range) * 0.1
    x_lim <- c(x_range[1] - x_margin, x_range[2] + x_margin)
    y_lim <- c(0, y_range[2] + y_margin)
    
    df$hover_text <- paste(
      "ID:", df$`Master protein IDs`,
      "<br>log2FC:", round(df$log2FC, 3),
      "<br>-log10(P):", round(df$log10P, 3),
      "<br>Regulation:", df$regulation
    )
    df$customdata <- df$`Master protein IDs`
    
    p <- plot_ly(
      data = df,
      x = ~log2FC,
      y = ~log10P,
      text = ~hover_text,
      hoverinfo = "text",
      customdata = ~customdata,
      color = ~regulation,
      colors = color_vec,
      type = "scatter",
      mode = "markers",
      marker = list(size = point_size, opacity = 0.6),
      source = "volcano_plot"
    ) %>%
      event_register("plotly_click") %>%
      layout(
        title = plot_title,
        xaxis = list(title = "log2(Fold Change)", range = x_lim),
        yaxis = list(title = "-log10(P-value)", range = y_lim),
        legend = list(title = list(text = "Regulation")),
        # 添加阈值线
        shapes = list(
          list(type = "line", x0 = x_up, x1 = x_up, y0 = y_lim[1], y1 = y_lim[2],
               line = list(dash = "dash", color = "gray40", width = 1)),
          list(type = "line", x0 = x_down, x1 = x_down, y0 = y_lim[1], y1 = y_lim[2],
               line = list(dash = "dash", color = "gray40", width = 1)),
          list(type = "line", x0 = x_lim[1], x1 = x_lim[2], y0 = y_cut, y1 = y_cut,
               line = list(dash = "dash", color = "gray40", width = 1))
        )
      )
    message("[DEBUG] volcano: plotly object created with threshold lines")
    p
  }, error = function(e) {
    message("[ERROR] volcano: plotly rendering failed: ", e$message)
    plotly::plot_ly() %>% layout(title = paste("Error:", e$message))
  })
})

# ---------- 点击事件：显示蛋白表达谱 ----------
observeEvent(event_data("plotly_click", source = "volcano_plot"), {
  cd <- event_data("plotly_click", source = "volcano_plot")$customdata
  if (!is.null(cd) && cd != "") {
    message("[DEBUG] volcano: clicked protein = ", cd)
    clicked_protein(list(id = cd, ts = Sys.time()))
  }
})

observeEvent(clicked_protein(), {
  req(clicked_protein())
  pro_id <- clicked_protein()$id
  message("[DEBUG] volcano: showing expression profile for ", pro_id)
  
  showModal(modalDialog(
    title = div(icon("chart-line"), paste("Expression Profile:", pro_id)),
    size = "l",
    easyClose = TRUE,
    footer = modalButton("Close"),
    plotlyOutput("protein_profile_plot", height = "450px")
  ))
  
  observeEvent(input$`modalButton`, {
    message("[DEBUG] volcano: modal closed")
  }, once = TRUE)
})

# ---------- 蛋白表达谱绘图 ----------
output$protein_profile_plot <- renderPlotly({
  req(clicked_protein())
  pro_id <- clicked_protein()$id
  
  nd <- tryCatch({
    norm_data_full()
  }, error = function(e) {
    NULL
  })
  
  if (is.null(nd)) {
    nd <- get_analysis_matrix()
    if (is.null(nd)) return(plotly::plot_ly() %>% layout(title = "Normalized data not available."))
  }
  
  idx <- which(nd[["Master protein IDs"]] == pro_id)
  if (length(idx) == 0) {
    return(plotly::plot_ly() %>% layout(title = paste("Protein", pro_id, "not found.")))
  }
  if (length(idx) > 1) idx <- idx[1]
  
  mat <- get_analysis_matrix()
  if (!is.null(mat) && pro_id %in% rownames(mat)) {
    intensities <- as.numeric(mat[pro_id, ])
    samples <- extract_sample_names(colnames(mat))
  } else {
    norm_cols <- grep("^Norm_LFQ intensity ", colnames(nd), value = TRUE)
    if (length(norm_cols) == 0) norm_cols <- grep("^LFQ intensity ", colnames(nd), value = TRUE)
    intensities <- as.numeric(nd[idx, norm_cols])
    samples <- gsub("^(Norm_)?LFQ intensity ", "", norm_cols)
  }
  
  all_groups <- rv$groups
  group_map <- setNames(rep(names(all_groups), lengths(all_groups)), unlist(all_groups))
  groups_vec <- group_map[samples]
  groups_vec[is.na(groups_vec)] <- "Unassigned"
  
  df_plot <- data.frame(
    Sample = samples,
    Group = groups_vec,
    Intensity = intensities,
    stringsAsFactors = FALSE
  )
  df_plot$Intensity[is.na(df_plot$Intensity)] <- 0
  
  p <- ggplot(df_plot, aes(x = Group, y = Intensity, color = Group,
                           text = paste("Sample:", Sample, "<br>Intensity:", round(Intensity, 2)))) +
    geom_jitter(width = 0.15, size = 3, alpha = 0.7) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "black", fatten = 2) +
    labs(title = paste("Expression profile of", pro_id), y = "Intensity", x = "") +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1),
          plot.background = element_rect(fill = "white", color = NA))
  
  ggplotly(p, tooltip = "text") %>% 
    layout(plot_bgcolor = 'white', paper_bgcolor = 'white', margin = list(b = 80))
})

# ---------- 下载火山图 PNG ----------
output$download_volcano_png <- downloadHandler(
  filename = function() { paste0("Volcano_", Sys.Date(), ".png") },
  content = function(file) {
    res <- volcano_de_result()
    if (is.null(res)) return()
    df <- res
    if (!"Master protein IDs" %in% colnames(df)) df$`Master protein IDs` <- rownames(df)
    cols <- color_mapping_vector()
    p <- ggplot(df, aes(x = log2FC, y = log10P, color = regulation)) +
      geom_point(alpha = 0.6, size = 2) +
      scale_color_manual(values = cols) +
      labs(title = paste("Volcano Plot:", selected_volcano_comparison()$name),
           x = "log2(Fold Change)", y = "-log10(P-value)") +
      theme_bw()
    ggsave(file, plot = p, width = 8, height = 6, dpi = 150)
  }
)

message("[DEBUG] volcano_plot.R loaded successfully.")