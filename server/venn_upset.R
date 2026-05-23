# server/venn_upset.R

selected_venn_comps <- reactive({
  if (input$venn_upset_method == "venn") input$venn_comparisons_select else input$venn_comparisons_checkbox
})

observe({
  comps <- sapply(sorted_comps(), `[[`, "name")
  if (length(comps) > 0) {
    updateSelectizeInput(session, "venn_comparisons_select", choices = comps, server = TRUE)
  } else {
    updateSelectizeInput(session, "venn_comparisons_select", choices = character(0), selected = character(0), server = TRUE)
  }
})

observeEvent(input$venn_upset_method, {
  updateSelectizeInput(session, "venn_comparisons_select", selected = character(0))
  updateCheckboxGroupInput(session, "venn_comparisons_checkbox", selected = character(0))
})

output$venn_message_ui <- renderUI({
  if (is.null(selected_venn_comps()) || length(selected_venn_comps()) < 2)
    div(class = "input-warning", "Please select at least 2 comparisons.")
  else if (input$venn_upset_method == "venn" && length(selected_venn_comps()) > 5)
    div(class = "input-warning", "Venn diagram supports at most 5 comparisons. Switch to UpSet or reduce selection.")
  else if (input$venn_upset_method == "upset" && length(selected_venn_comps()) > 15)
    div(class = "input-warning", "UpSet is limited to 15 comparisons due to memory constraints. Please reduce your selection.")
  else NULL
})

venn_colors_for_plot <- reactive({
  n <- length(selected_venn_comps())
  if (n == 0) return(character(0))
  if (n <= 8) RColorBrewer::brewer.pal(max(3, n), "Set2")[1:n]
  else colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n)
})

venn_data <- eventReactive(input$generate_venn, {
  all_res <- all_analysis_results()
  comps <- selected_venn_comps()
  if (is.null(comps) || length(comps) < 2) return(list(error = "Please select at least 2 comparisons."))
  if (input$venn_upset_method == "upset" && length(comps) > 15) {
    showNotification("UpSet is limited to 15 comparisons. Only the first 15 will be used.", type = "warning")
    comps <- comps[1:15]
  }
  reg_type <- input$venn_type
  protein_sets <- lapply(comps, function(comp_name) {
    df <- all_res$results[[comp_name]]$data
    ids <- df[df$regulation == reg_type, "Master protein IDs"]
    ids <- ids[!is.na(ids) & ids != ""]
    unique(ids)
  })
  names(protein_sets) <- comps
  
  n <- length(comps)
  combn_ids <- expand.grid(rep(list(c(FALSE, TRUE)), n))
  combn_ids <- combn_ids[-1, ]
  region_list <- lapply(1:nrow(combn_ids), function(i) {
    flags <- as.logical(combn_ids[i, ])
    region_name <- paste(comps[flags], collapse = " & ")
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
  list(sets = protein_sets, regions = region_list)
})

plot_venn_diagram <- function(sets, fill_colors = NULL, margin_val = 0.15, cat_cex_val = 1.5, cex_val = 1.5) {
  n <- length(sets)
  if (n < 2 || n > 5) return()
  if (is.null(fill_colors)) fill_colors <- RColorBrewer::brewer.pal(5, "Pastel1")[1:n]
  if (n == 2) {
    VennDiagram::draw.pairwise.venn(area1 = length(sets[[1]]), area2 = length(sets[[2]]),
                                    cross.area = length(intersect(sets[[1]], sets[[2]])),
                                    category = names(sets)[1:2], col = "black", fill = fill_colors[1:2],
                                    alpha = 0.5, margin = margin_val, cat.cex = cat_cex_val, cex = cex_val)
  } else if (n == 3) {
    VennDiagram::draw.triple.venn(area1 = length(sets[[1]]), area2 = length(sets[[2]]), area3 = length(sets[[3]]),
                                  n12 = length(intersect(sets[[1]], sets[[2]])),
                                  n23 = length(intersect(sets[[2]], sets[[3]])),
                                  n13 = length(intersect(sets[[1]], sets[[3]])),
                                  n123 = length(Reduce(intersect, sets)),
                                  category = names(sets)[1:3], col = "black", fill = fill_colors[1:3],
                                  alpha = 0.5, margin = margin_val, cat.cex = cat_cex_val, cex = cex_val)
  } else if (n == 4) {
    VennDiagram::draw.quad.venn(area1 = length(sets[[1]]), area2 = length(sets[[2]]), area3 = length(sets[[3]]), area4 = length(sets[[4]]),
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
                                category = names(sets)[1:4], col = "black", fill = fill_colors[1:4],
                                alpha = 0.5, margin = margin_val, cat.cex = cat_cex_val, cex = cex_val)
  } else if (n == 5) {
    VennDiagram::draw.quintuple.venn(area1 = length(sets[[1]]), area2 = length(sets[[2]]), area3 = length(sets[[3]]), area4 = length(sets[[4]]), area5 = length(sets[[5]]),
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
                                     category = names(sets)[1:5], col = "black", fill = fill_colors[1:5],
                                     alpha = 0.5, margin = margin_val, cat.cex = cat_cex_val, cex = cex_val)
  }
}

output$venn_plot <- renderPlot({
  vd <- venn_data(); req(vd)
  if (!is.null(vd$error)) { plot.new(); text(0.5, 0.5, vd$error, cex = 1.5); return() }
  sets <- vd$sets; colors <- venn_colors_for_plot()
  if (input$venn_upset_method == "venn" && length(sets) > 5) {
    plot.new(); text(0.5, 0.5, "Venn diagram limited to 5 comparisons.\nUse UpSet plot for more.", cex = 1.5)
  } else {
    plot_venn_diagram(sets, fill_colors = colors, margin_val = 0.2, cat_cex_val = 1.2, cex_val = 1.2)
  }
}, width = 800, height = 600)

output$upset_plot <- renderPlot({
  vd <- venn_data(); req(vd)
  if (!is.null(vd$error)) { plot.new(); text(0.5, 0.5, vd$error, cex = 1.5); return() }
  sets <- vd$sets
  if (length(sets) < 2) { plot.new(); text(0.5, 0.5, "Not enough comparisons selected."); return() }
  tryCatch({
    UpSetR::upset(UpSetR::fromList(sets), nsets = length(sets), order.by = "freq", decreasing = TRUE,
                  main.bar.color = "#3498db", sets.bar.color = venn_colors_for_plot(),
                  matrix.color = "#34495e", mainbar.y.label = "Intersection Size", sets.x.label = "Set Size")
  }, error = function(e) { plot.new(); text(0.5, 0.5, paste("Error generating UpSet plot:", e$message)) })
}, width = 900, height = 600)

output$download_venn_png <- downloadHandler(
  filename = function() paste0("Venn_", input$venn_type, "_", Sys.Date(), ".png"),
  content = function(file) {
    vd <- venn_data(); req(vd)
    png(file, width = 1200, height = 900, res = 150)
    if (!is.null(vd$error)) { plot.new(); text(0.5, 0.5, vd$error) }
    else if (input$venn_upset_method == "venn" && length(vd$sets) > 5) { plot.new(); text(0.5, 0.5, "Too many sets for Venn") }
    else plot_venn_diagram(vd$sets, fill_colors = venn_colors_for_plot(), margin_val = 0.5, cat_cex_val = 1.8, cex_val = 1.3)
    dev.off()
  }
)

output$download_upset_png <- downloadHandler(
  filename = function() paste0("UpSet_", input$venn_type, "_", Sys.Date(), ".png"),
  content = function(file) {
    vd <- venn_data(); req(vd)
    png(file, width = 1200, height = 800, res = 150)
    if (!is.null(vd$error)) { plot.new(); text(0.5, 0.5, vd$error) }
    else if (length(vd$sets) < 2) { plot.new(); text(0.5, 0.5, "Not enough comparisons selected.") }
    else {
      tryCatch({
        p <- UpSetR::upset(UpSetR::fromList(vd$sets), nsets = length(vd$sets), order.by = "freq",
                           decreasing = TRUE, main.bar.color = "#3498db",
                           sets.bar.color = venn_colors_for_plot(), matrix.color = "#34495e")
        print(p)
      }, error = function(e) { plot.new(); text(0.5, 0.5, paste("Error:", e$message)) })
    }
    dev.off()
  }
)

output$venn_region_select_ui <- renderUI({
  req(venn_data()); vd <- venn_data()
  if (!is.null(vd$error)) return(NULL)
  regions <- vd$regions; region_names <- names(regions)
  counts <- sapply(regions, function(x) length(x$proteins))
  region_choices <- paste0(region_names, " (", counts, " proteins)")
  selectInput("venn_region", "Select Region", choices = setNames(region_names, region_choices), selected = region_names[1])
})

output$venn_region_table <- renderDT({
  req(venn_data(), input$venn_region)
  vd <- venn_data()
  if (!is.null(vd$error)) return(datatable(data.frame(Message = vd$error)))
  proteins <- vd$regions[[input$venn_region]]$proteins
  if (length(proteins) == 0) df <- data.frame(Protein = character(0))
  else df <- data.frame(`Master Protein IDs` = proteins, check.names = FALSE, stringsAsFactors = FALSE)
  datatable(df, rownames = FALSE, options = list(pageLength = 10, dom = 'Bfrtip'))
})

output$download_venn_region <- downloadHandler(
  filename = function() paste0("Region_", input$venn_type, "_", gsub(" ", "_", input$venn_region), "_", Sys.Date(), ".csv"),
  content = function(file) {
    vd <- venn_data(); req(vd)
    if (!is.null(vd$error)) return()
    write.csv(data.frame(`Master Protein IDs` = vd$regions[[input$venn_region]]$proteins, check.names = FALSE), file, row.names = FALSE)
  }
)

output$stat_table <- renderDT({
  req(all_analysis_results())
  all_res <- all_analysis_results()
  results_list <- all_res$results
  if (length(results_list) == 0) return(datatable(data.frame(Message="No results")))
  stats_list <- lapply(names(results_list), function(nm) {
    res <- results_list[[nm]]; cnt <- attr(res$data, "counts")
    data.frame(Comparison = nm, Up = cnt$Up, Down = cnt$Down, Increase = cnt$Increase, Decrease = cnt$Decrease, NS = cnt$NS,
               Total = cnt$Up + cnt$Down + cnt$Increase + cnt$Decrease, stringsAsFactors = FALSE)
  })
  stats_df <- do.call(rbind, stats_list)
  datatable(stats_df, rownames = FALSE, selection = 'single',
            options = list(pageLength = 100, dom = 'Bt'),
            class = 'display compact stripe hover')
})