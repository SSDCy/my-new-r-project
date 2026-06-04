# ui/plots_ui.R

plots_ui <- function() {
  tabPanel(
    title = div(icon("chart-bar"), "Plots"),
    value = "plots",
    br(),
    tabsetPanel(
      id = "plots_subnav",
      # ---- 热图选项卡 ----
      tabPanel(
        title = "Heatmap",
        value = "heatmap_sub",
        fluidRow(
          column(12,
                 div(class = "card-modern",
                     div(class = "card-header-modern", icon("th"), " Differential Protein Heatmap"),
                     div(style = "padding: 20px;",
                         p("Expression patterns of proteins. Per-protein log2 + Z-score normalization is applied, followed by hierarchical clustering."),
                         fluidRow(
                           column(4,
                                  radioButtons("heatmap_data_source", "Data Source",
                                               choices = c("LFQ Intensity (per-row Z-score)" = "LFQ",
                                                           "Intensity (per-row Z-score)" = "Intensity"),
                                               selected = "LFQ"),
                                  verbatimTextOutput("heatmap_data_source_info"),
                                  hr(),
                                  radioButtons("heatmap_protein_mode", "Protein Selection Mode",
                                               choices = c("Top N proteins (by variance)" = "top_n",
                                                           "Custom protein list" = "custom"),
                                               selected = "top_n"),
                                  conditionalPanel(
                                    condition = "input.heatmap_protein_mode == 'top_n'",
                                    numericInput("heatmap_top_n", "Top N proteins", value = 20, min = 5, max = 200, step = 5)
                                  ),
                                  conditionalPanel(
                                    condition = "input.heatmap_protein_mode == 'custom'",
                                    textAreaInput("heatmap_custom_ids", "Enter Master Protein IDs", rows = 5, placeholder = "P12345\nP67890")
                                  ),
                                  hr(),
                                  conditionalPanel(
                                    condition = "input.heatmap_protein_mode == 'top_n'",
                                    selectInput("heatmap_groups", "Select Groups", choices = NULL, multiple = TRUE)
                                  ),
                                  conditionalPanel(
                                    condition = "input.heatmap_data_source == 'Intensity'",
                                    selectInput("heatmap_group_level", "Grouping Level", choices = NULL),
                                    actionButton("heatmap_apply_grouping", "Apply Grouping"),
                                    uiOutput("heatmap_group_selection_ui")
                                  ),
                                  checkboxInput("heatmap_show_sample_names", "Show sample names", value = TRUE),
                                  actionButton("generate_heatmap", "Generate Heatmap", class = "btn btn-primary btn-block"),
                                  hr(),
                                  downloadButton("download_heatmap_png", "Download Heatmap PNG", class = "btn btn-sm btn-outline-success")
                           ),
                           column(8,
                                  shinycssloaders::withSpinner(plotOutput("heatmap_plot", height = "700px"), type = 4, color = "#e67e22")
                           )
                         )
                     )
                 )
          )
        )
      ),
      # ---- 火山图选项卡（左对齐） ----
      tabPanel(
        title = "Volcano Plot",
        value = "volcano_sub",
        fluidRow(
          column(12,
                 div(class = "card-modern",
                     div(class = "card-header-modern",
                         icon("chart-line"),
                         " Volcano Plot Visualization ",
                         span(style = "font-size: 14px; font-weight: normal;",
                              "(Click on any dot to view the protein's expression profile across all groups)")),
                     div(style = "padding: 20px;",
                         # 颜色选择器行
                         div(class = "color-palette-row",
                             div(class = "color-card",
                                 div(class = "color-card-label", "Up"),
                                 colourpicker::colourInput("color_up", NULL, value = "#FF0000",
                                                           showColour = "background",
                                                           allowTransparent = FALSE),
                                 div(class = "color-card-value", textOutput("val_up", inline = TRUE))
                             ),
                             div(class = "color-card",
                                 div(class = "color-card-label", "Down"),
                                 colourpicker::colourInput("color_down", NULL, value = "#0000FF",
                                                           showColour = "background",
                                                           allowTransparent = FALSE),
                                 div(class = "color-card-value", textOutput("val_down", inline = TRUE))
                             ),
                             div(class = "color-card",
                                 div(class = "color-card-label", "Increase"),
                                 colourpicker::colourInput("color_increase", NULL, value = "#C00000",
                                                           showColour = "background",
                                                           allowTransparent = FALSE),
                                 div(class = "color-card-value", textOutput("val_inc", inline = TRUE))
                             ),
                             div(class = "color-card",
                                 div(class = "color-card-label", "Decrease"),
                                 colourpicker::colourInput("color_decrease", NULL, value = "#0945A5",
                                                           showColour = "background",
                                                           allowTransparent = FALSE),
                                 div(class = "color-card-value", textOutput("val_dec", inline = TRUE))
                             ),
                             div(class = "color-card",
                                 div(class = "color-card-label", "NS"),
                                 colourpicker::colourInput("color_ns", NULL, value = "#7f7e83",
                                                           showColour = "background",
                                                           allowTransparent = FALSE),
                                 div(class = "color-card-value", textOutput("val_ns", inline = TRUE))
                             ),
                             div(class = "reset-btn-wrapper",
                                 actionButton("reset_color", icon("undo"), class = "btn-circle-modern")
                             )
                         ),
                         uiOutput("color_preview"),
                         fluidRow(
                           column(4, numericInput("point_size", "Point Size", value = 4, min = 0.5, max = 10, step = 0.1))
                         ),
                         hr(),
                         fluidRow(
                           column(8,
                                  selectInput("selected_comparison", "Select Comparison", choices = NULL, width = "100%")
                           ),
                           column(4,
                                  textInputMax("single_plot_title", "Custom Title", value = "", placeholder = "Auto-generated", maxlength = 50, width = "100%")
                           )
                         ),
                         fluidRow(column(12, uiOutput("plot_info_ui"))),
                         fluidRow(
                           column(12,
                                  shinycssloaders::withSpinner(plotlyOutput("volcano_plot", height = "700px"), type = 4, color = "#3498db")
                           )
                         )
                     )
                 )
          )
        )
      )
    )
  )
}