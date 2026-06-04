# server/plots_venn_upset.R
message("[DEBUG] plots_venn_upset.R loading...")

# ---------- 韦恩图/UpSet图 辅助函数 ----------
plot_venn_diagram <- function(sets, fill_colors = NULL) {
  n <- length(sets)
  if (n < 2 || n > 5) {
    message("[DEBUG] plot_venn_diagram: number of sets out of range (2-5): ", n)
    return()
  }
  if (is.null(fill_colors)) fill_colors <- RColorBrewer::brewer.pal(5, "Pastel1")[1:n]
  if (n == 2) {
    VennDiagram::draw.pairwise.venn(area1=length(sets[[1]]), area2=length(sets[[2]]), cross.area=length(intersect(sets[[1]],sets[[2]])), category=names(sets)[1:2], col="black", fill=fill_colors[1:2], alpha=0.5, cat.cex=1.2, cex=1.2)
  } else if (n == 3) {
    VennDiagram::draw.triple.venn(area1=length(sets[[1]]), area2=length(sets[[2]]), area3=length(sets[[3]]), n12=length(intersect(sets[[1]],sets[[2]])), n23=length(intersect(sets[[2]],sets[[3]])), n13=length(intersect(sets[[1]],sets[[3]])), n123=length(Reduce(intersect, sets)), category=names(sets)[1:3], col="black", fill=fill_colors[1:3], alpha=0.5, cat.cex=1.2, cex=1.2)
  } else if (n == 4) {
    VennDiagram::draw.quad.venn(area1=length(sets[[1]]), area2=length(sets[[2]]), area3=length(sets[[3]]), area4=length(sets[[4]]), n12=length(intersect(sets[[1]],sets[[2]])), n13=length(intersect(sets[[1]],sets[[3]])), n14=length(intersect(sets[[1]],sets[[4]])), n23=length(intersect(sets[[2]],sets[[3]])), n24=length(intersect(sets[[2]],sets[[4]])), n34=length(intersect(sets[[3]],sets[[4]])), n123=length(Reduce(intersect, sets[1:3])), n124=length(Reduce(intersect, sets[c(1,2,4)])), n134=length(Reduce(intersect, sets[c(1,3,4)])), n234=length(Reduce(intersect, sets[2:4])), n1234=length(Reduce(intersect, sets)), category=names(sets)[1:4], col="black", fill=fill_colors[1:4], alpha=0.5, cat.cex=1.2, cex=1.2)
  } else if (n == 5) {
    VennDiagram::draw.quintuple.venn(area1=length(sets[[1]]), area2=length(sets[[2]]), area3=length(sets[[3]]), area4=length(sets[[4]]), area5=length(sets[[5]]), n12=length(intersect(sets[[1]],sets[[2]])), n13=length(intersect(sets[[1]],sets[[3]])), n14=length(intersect(sets[[1]],sets[[4]])), n15=length(intersect(sets[[1]],sets[[5]])), n23=length(intersect(sets[[2]],sets[[3]])), n24=length(intersect(sets[[2]],sets[[4]])), n25=length(intersect(sets[[2]],sets[[5]])), n34=length(intersect(sets[[3]],sets[[4]])), n35=length(intersect(sets[[3]],sets[[5]])), n45=length(intersect(sets[[4]],sets[[5]])), n123=length(Reduce(intersect, sets[1:3])), n124=length(Reduce(intersect, sets[c(1,2,4)])), n125=length(Reduce(intersect, sets[c(1,2,5)])), n134=length(Reduce(intersect, sets[c(1,3,4)])), n135=length(Reduce(intersect, sets[c(1,3,5)])), n145=length(Reduce(intersect, sets[c(1,4,5)])), n234=length(Reduce(intersect, sets[2:4])), n235=length(Reduce(intersect, sets[c(2,3,5)])), n245=length(Reduce(intersect, sets[c(2,4,5)])), n345=length(Reduce(intersect, sets[3:5])), n1234=length(Reduce(intersect, sets[1:4])), n1235=length(Reduce(intersect, sets[c(1,2,3,5)])), n1245=length(Reduce(intersect, sets[c(1,2,4,5)])), n1345=length(Reduce(intersect, sets[c(1,3,4,5)])), n2345=length(Reduce(intersect, sets[2:5])), n12345=length(Reduce(intersect, sets)), category=names(sets)[1:5], col="black", fill=fill_colors[1:5], alpha=0.5, cat.cex=1.2, cex=1.2)
  }
}

# ---------- 颜色选择器 ----------
venn_colors_for_plot <- reactive({
  n <- length(input$venn_comparisons)
  if (n == 0) return(character(0))
  if (n <= 8) {
    RColorBrewer::brewer.pal(max(3, n), "Set2")[1:n]
  } else {
    colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n)
  }
})

# ---------- 更新比较选择器 ----------
observe({
  comp_names <- sapply(sorted_comps(), function(c) c$name)
  updateSelectizeInput(session, "venn_comparisons", choices = comp_names, server = TRUE)
  message("[DEBUG] plots_venn_upset: updated venn_comparisons choices, n = ", length(comp_names))
})

# ---------- 消息提示 ----------
output$venn_message_ui <- renderUI({
  if (is.null(input$venn_comparisons) || length(input$venn_comparisons) < 2)
    div(class = "input-warning", "Please select at least 2 comparisons.")
  else if (input$venn_upset_method == "venn" && length(input$venn_comparisons) > 5)
    div(class = "input-warning", "Venn diagram supports at most 5 comparisons. Switch to UpSet or reduce selection.")
  else
    NULL
})

# ---------- 韦恩图/UpSet 数据 ----------
venn_data <- eventReactive(input$generate_venn, {
  req(all_analysis_results(), input$venn_comparisons, input$venn_type)
  req(length(input$venn_comparisons) >= 2)
  
  all_res <- all_analysis_results()
  comps <- input$venn_comparisons
  reg_type <- input$venn_type
  message("[DEBUG] plots_venn_upset: generating data for ", reg_type, " in comparisons: ", paste(comps, collapse = ", "))
  
  protein_sets <- lapply(comps, function(comp_name) {
    res <- all_res$results[[comp_name]]
    if (is.null(res)) {
      message("[DEBUG] plots_venn_upset: comparison ", comp_name, " not found in results")
      return(character(0))
    }
    df <- res$data
    ids <- df[df$regulation == reg_type, "Master protein IDs"]
    ids <- ids[!is.na(ids) & ids != ""]
    message("[DEBUG] plots_venn_upset: ", comp_name, " - ", length(ids), " proteins")
    unique(ids)
  })
  names(protein_sets) <- comps
  
  non_empty <- lengths(protein_sets) > 0
  if (!all(non_empty)) {
    message("[DEBUG] plots_venn_upset: empty sets removed: ", paste(names(protein_sets)[!non_empty], collapse = ", "))
    protein_sets <- protein_sets[non_empty]
  }
  if (length(protein_sets) < 2) {
    message("[DEBUG] plots_venn_upset: not enough non-empty sets to generate venn/upset")
    return(list(sets = list(), regions = list()))
  }
  
  n <- length(protein_sets)
  combn_ids <- expand.grid(rep(list(c(FALSE, TRUE)), n))
  combn_ids <- combn_ids[-1, ]
  region_list <- lapply(1:nrow(combn_ids), function(i) {
    flags <- as.logical(combn_ids[i, ])
    region_name <- paste(names(protein_sets)[flags], collapse = " & ")
    sets_to_intersect <- protein_sets[flags]
    others <- protein_sets[!flags]
    intersect_set <- Reduce(intersect, sets_to_intersect)
    if (length(others) > 0) {
      union_others <- Reduce(union, others)
      result <- setdiff(intersect_set, union_others)
    } else result <- intersect_set
    list(name = region_name, proteins = result)
  })
  names(region_list) <- sapply(region_list, `[[`, "name")
  message("[DEBUG] plots_venn_upset: generated ", length(region_list), " regions")
  list(sets = protein_sets, regions = region_list)
})

# ---------- 韦恩图 ----------
output$venn_plot <- renderPlot({
  req(venn_data())
  sets <- venn_data()$sets
  if (length(sets) == 0) {
    plot.new()
    text(0.5, 0.5, "No data to display")
    return()
  }
  colors <- venn_colors_for_plot()
  if (input$venn_upset_method == "venn" && length(sets) > 5) {
    plot.new()
    text(0.5, 0.5, "Venn diagram limited to 5 comparisons.\nUse UpSet plot for more.", cex = 1.5)
  } else {
    plot_venn_diagram(sets, fill_colors = colors)
  }
})

# ---------- UpSet图 ----------
output$upset_plot <- renderPlot({
  req(venn_data())
  sets <- venn_data()$sets
  if (length(sets) == 0) {
    plot.new()
    text(0.5, 0.5, "No data to display")
    return()
  }
  colors <- venn_colors_for_plot()
  
  upset_list <- UpSetR::fromList(sets)
  colnames(upset_list) <- names(sets)
  message("[DEBUG] UpSet columns: ", paste(colnames(upset_list), collapse = ", "))
  
  tryCatch({
    UpSetR::upset(
      upset_list,
      nsets = length(sets),
      order.by = "freq",
      decreasing = TRUE,
      main.bar.color = "#3498db",
      sets.bar.color = colors,
      matrix.color = "#34495e",
      mainbar.y.label = "Intersection Size",
      sets.x.label = "Set Size",
      text.scale = c(1.5, 1.5, 1.2, 1.2, 1.5, 1.2)
    )
  }, error = function(e) {
    message("[ERROR] UpSet plot: ", e$message)
    plot.new()
    text(0.5, 0.5, paste("Error generating UpSet plot:\n", e$message), cex = 1.2)
  })
})

# ---------- 下载韦恩图 PNG ----------
output$download_venn_png <- downloadHandler(
  filename = function() paste0("Venn_", input$venn_type, "_", Sys.Date(), ".png"),
  content = function(file) {
    sets <- venn_data()$sets
    if (length(sets) == 0) {
      png(file, 800, 600); plot.new(); text(0.5,0.5,"No data"); dev.off()
      return()
    }
    colors <- venn_colors_for_plot()
    png(file, width = 900, height = 900, res = 120)
    if (input$venn_upset_method == "venn" && length(sets) > 5) {
      plot.new(); text(0.5,0.5,"Too many sets for Venn")
    } else {
      plot_venn_diagram(sets, fill_colors = colors)
    }
    dev.off()
  }
)

# ---------- 下载UpSet图 PNG ----------
output$download_upset_png <- downloadHandler(
  filename = function() paste0("UpSet_", input$venn_type, "_", Sys.Date(), ".png"),
  content = function(file) {
    sets <- venn_data()$sets
    if (length(sets) == 0) {
      png(file, 800, 600); plot.new(); text(0.5,0.5,"No data"); dev.off()
      return()
    }
    colors <- venn_colors_for_plot()
    upset_list <- UpSetR::fromList(sets)
    colnames(upset_list) <- names(sets)
    png(file, width = 1200, height = 800, res = 120)
    tryCatch({
      UpSetR::upset(upset_list, nsets = length(sets), order.by = "freq", decreasing = TRUE,
                    main.bar.color = "#3498db", sets.bar.color = colors, matrix.color = "#34495e",
                    text.scale = c(1.5, 1.5, 1.2, 1.2, 1.5, 1.2))
    }, error = function(e) {
      plot.new(); text(0.5,0.5, paste("Error:", e$message))
    })
    dev.off()
  }
)

# ---------- 区域选择器 ----------
output$venn_region_select_ui <- renderUI({
  req(venn_data())
  regions <- venn_data()$regions
  if (length(regions) == 0) return(NULL)
  region_names <- names(regions)
  counts <- sapply(regions, function(x) length(x$proteins))
  region_choices <- paste0(region_names, " (", counts, " proteins)")
  selectInput("venn_region", "Select Region", choices = setNames(region_names, region_choices), selected = region_names[1])
})

# ---------- 区域蛋白表格 ----------
output$venn_region_table <- renderDT({
  req(venn_data(), input$venn_region)
  proteins <- venn_data()$regions[[input$venn_region]]$proteins
  if (length(proteins) == 0) df <- data.frame(Protein = character(0))
  else df <- data.frame(`Master Protein IDs` = proteins, check.names = FALSE, stringsAsFactors = FALSE)
  datatable(df, rownames = FALSE, options = list(pageLength = 10, dom = 'Bfrtip'))
})

# ---------- 下载区域蛋白 CSV ----------
output$download_venn_region <- downloadHandler(
  filename = function() paste0("Region_", input$venn_type, "_", gsub(" ", "_", input$venn_region), "_", Sys.Date(), ".csv"),
  content = function(file) {
    write.csv(data.frame(`Master Protein IDs` = venn_data()$regions[[input$venn_region]]$proteins, check.names = FALSE), file, row.names = FALSE)
  }
)

# ---------- 步骤指示器（可折叠） ----------
output$venn_preprocess_steps <- renderUI({
  steps <- list()
  steps <- c(steps, paste0("Missing Value Filter: threshold = ", input$max_missing_fraction %||% 0.5,
                           ", mode = ", preprocessing_params$missing_filter_mode %||% "global"))
  steps <- c(steps, paste0("Minimum Intensity Filter: threshold = ", input$min_intensity,
                           ", min samples = ", input$min_samples_above_intensity %||% 1))
  imp <- preprocessing_params$imputation_method %||% "none"
  if (imp == "none") {
    steps <- c(steps, "Missing Value Imputation: none")
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
    tags$summary("Data preprocessing steps for Venn/UpSet", style = "cursor: pointer; font-weight: bold; color: #2c3e50; margin-bottom: 10px;"),
    div(style = "display: flex; flex-wrap: wrap; align-items: center;", do.call(tagList, step_tags))
  )
})

message("[DEBUG] plots_venn_upset.R loaded successfully.")