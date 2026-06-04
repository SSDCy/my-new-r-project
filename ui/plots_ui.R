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
                         uiOutput("heatmap_preprocess_steps"),
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
      # ---- 火山图选项卡 ----
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
                         uiOutput("volcano_preprocess_steps"),
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
      ),
      # ---- 韦恩图 / UpSet 选项卡 ----
      tabPanel(
        title = "Venn / UpSet",
        value = "venn_upset_sub",
        fluidRow(
          column(12,
                 div(class = "card-modern",
                     div(class = "card-header-modern", icon("chart-pie"), " Shared & Unique Proteins"),
                     div(style = "padding: 20px;",
                         uiOutput("venn_preprocess_steps"),
                         p("Select regulation type and at least 2 comparisons. Venn diagram works best with 2-5; UpSet can handle all."),
                         radioButtons("venn_upset_method", "Visualization Method",
                                      choices = c("Venn Diagram (max 5)" = "venn", "UpSet Plot (unlimited)" = "upset"),
                                      selected = "venn", inline = TRUE),
                         fluidRow(
                           column(4, selectInput("venn_type", "Regulation Type", choices = c("Up", "Down", "Increase", "Decrease"), selected = "Up")),
                           column(6, selectizeInput("venn_comparisons", "Comparisons (min. 2)", choices = NULL, multiple = TRUE, 
                                                    options = list(placeholder = 'Select at least 2 comparisons'))),
                           column(2, br(), actionButton("generate_venn", "Generate", class = "btn btn-primary btn-block"))
                         ),
                         uiOutput("venn_message_ui"),
                         conditionalPanel(
                           condition = "input.venn_upset_method == 'venn'",
                           shinycssloaders::withSpinner(plotOutput("venn_plot", height = "500px"), type = 4, color = "#3498db"),
                           fluidRow(
                             column(6, downloadButton("download_venn_png", "Download Venn PNG", class = "btn btn-sm btn-outline-success"))
                           )
                         ),
                         conditionalPanel(
                           condition = "input.venn_upset_method == 'upset'",
                           shinycssloaders::withSpinner(plotOutput("upset_plot", height = "600px"), type = 4, color = "#3498db"),
                           fluidRow(
                             column(6, downloadButton("download_upset_png", "Download UpSet PNG", class = "btn btn-sm btn-outline-success"))
                           )
                         ),
                         uiOutput("venn_region_select_ui"),
                         shinycssloaders::withSpinner(DTOutput("venn_region_table"), type = 4, color = "#3498db"),
                         downloadButton("download_venn_region", "Download Region Proteins", class = "btn btn-sm btn-outline-secondary")
                     )
                 )
          )
        )
      ),
      # ---- PCA 选项卡 ----
      tabPanel(
        title = "PCA",
        value = "pca_sub",
        fluidRow(
          column(12,
                 div(class = "card-modern",
                     div(class = "card-header-modern", icon("project-diagram"), " Principal Component Analysis"),
                     div(style = "padding: 20px;",
                         uiOutput("pca_preprocess_steps"),
                         p("PCA is performed on normalized expression data after log2 transformation. Outliers are detected based on a Z-score > 3 on PC1 or PC2."),
                         uiOutput("pca_data_source_note"),
                         fluidRow(
                           column(6,
                                  h4(icon("users"), " PCA by Group"),
                                  shinycssloaders::withSpinner(plotlyOutput("pca_group_plot", height = "500px"), type = 4, color = "#3498db"),
                                  downloadButton("download_pca_group_png", "Download Group PCA", class = "btn btn-sm btn-outline-success")
                           ),
                           column(6,
                                  h4(icon("layer-group"), " PCA by Batch"),
                                  shinycssloaders::withSpinner(plotlyOutput("pca_batch_plot", height = "500px"), type = 4, color = "#3498db"),
                                  downloadButton("download_pca_batch_png", "Download Batch PCA", class = "btn btn-sm btn-outline-success")
                           )
                         ),
                         hr(),
                         h4(icon("exclamation-triangle"), " Outlier Detection"),
                         verbatimTextOutput("pca_outlier_info")
                     )
                 )
          )
        )
      ),
      # ---- 样本相关性热图 ----
      tabPanel(
        title = "Sample Correlation",
        value = "sample_cor_sub",
        fluidRow(
          column(12,
                 div(class = "card-modern",
                     div(class = "card-header-modern", icon("th"), " Sample Correlation Heatmap"),
                     div(style = "padding: 20px;",
                         uiOutput("sample_cor_preprocess_steps"),
                         p("Pearson correlation between samples based on the normalized expression data (log2 transformed). The heatmap uses the top 500 most variable proteins."),
                         uiOutput("sample_cor_data_source_note"),
                         fluidRow(
                           column(4,
                                  actionButton("generate_sample_cor", "Generate Correlation Heatmap", class = "btn btn-primary btn-block"),
                                  hr(),
                                  downloadButton("download_sample_cor_png", "Download Heatmap PNG", class = "btn btn-sm btn-outline-success"),
                                  downloadButton("download_sample_cor_matrix", "Download Correlation Matrix CSV", class = "btn btn-sm btn-outline-secondary")
                           ),
                           column(8,
                                  shinycssloaders::withSpinner(plotOutput("sample_cor_heatmap", height = "600px"), type = 4, color = "#e67e22")
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