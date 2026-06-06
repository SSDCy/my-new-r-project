# ui/plots_ui.R
message("[DEBUG] plots_ui.R loaded - export formats SVG/TIFF, sample names use underscores")

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
                                  conditionalPanel(
                                    condition = "input.heatmap_data_source == 'LFQ'",
                                    radioButtons("heatmap_normalization", "Normalization",
                                                 choices = c("None (preprocessed raw intensity)" = "none",
                                                             "Total Intensity (baseline sample)" = "total"),
                                                 selected = "none")
                                  ),
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
      # ---- 火山图选项卡（整合 Parameters 和 Export） ----
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
                         
                         # Parameters
                         tags$details(
                           tags$summary("Parameters", style = "cursor: pointer; font-weight: bold; color: #2c3e50; margin-top: 10px; margin-bottom: 10px;"),
                           div(style = "margin-top: 10px;",
                               fluidRow(
                                 column(6,
                                        div(style = "background: #f8f9fa; padding: 15px; border-radius: 10px; margin-bottom: 15px;",
                                            h5(icon("balance-scale"), " Fold Change Thresholds"),
                                            numericInput("fc_up", "FC up >", value = 1.2, min = 1, step = 0.1),
                                            uiOutput("fc_up_warning"),
                                            p(class = "param-hint", "Must be ≥ 1.0. Up-regulated proteins fold change threshold."),
                                            numericInput("fc_down", "FC down <", value = 0.84, min = 0, max = 1, step = 0.01),
                                            uiOutput("fc_down_warning"),
                                            p(class = "param-hint", "Must be between 0.0 and 1.0. Down-regulated proteins fold change threshold."),
                                            hr(),
                                            h5(icon("chart-line"), " Statistical Significance"),
                                            selectInput("p_cut", "P-value threshold", choices = c("0.05", "0.01"), selected = "0.05"),
                                            p(class = "param-hint", "Select significance level (P-value threshold)."),
                                            h5(icon("flask"), " Statistical Method"),
                                            selectInput("stat_method", "Method", choices = c("t-test", "wilcoxon", "limma"), selected = "t-test"),
                                            p(class = "param-hint", "Choose the statistical test for differential analysis.")
                                        )
                                 ),
                                 column(6,
                                        div(style = "background: #f8f9fa; padding: 15px; border-radius: 10px; margin-bottom: 15px;",
                                            h5(icon("check-circle"), " Valid Replicates"),
                                            div(style = "display: flex; align-items: center; gap: 10px; margin-bottom: 10px;",
                                                numericInput("replicate_fill_all", "Set All Thresholds", value = 2, min = 1, max = 10, step = 1, width = "100px"),
                                                actionButton("apply_replicate_fill", "Apply to All", class = "btn-sm btn-info")
                                            ),
                                            numericInput("min_treat_valid", "Treatment group min valid replicates", value = 2, min = 1, max = 20),
                                            uiOutput("min_treat_valid_warning"),
                                            numericInput("min_ctrl_valid", "Control group min valid replicates", value = 2, min = 1, max = 20),
                                            uiOutput("min_ctrl_valid_warning"),
                                            numericInput("min_rep_ttest", "Min replicates for t-test", value = 2, min = 1, max = 10),
                                            uiOutput("min_rep_ttest_warning"),
                                            numericInput("min_rep_inc", "Min replicates for 'Increase'", value = 2, min = 1, max = 10),
                                            uiOutput("min_rep_inc_warning"),
                                            numericInput("min_rep_dec", "Min replicates for 'Decrease'", value = 2, min = 1, max = 10),
                                            uiOutput("min_rep_dec_warning"),
                                            p(class = "param-hint", "Use the 'Set All' field above to fill all replicate thresholds at once, or adjust individually.")
                                        )
                                 )
                               ),
                               fluidRow(
                                 column(12,
                                        div(style = "background: #e3f2fd; padding: 15px; border-radius: 10px;",
                                            h5(icon("filter"), " Protein Filtering (Unique Peptides)"),
                                            numericInput("min_unique_pep", "Minimum Unique Peptides", value = 2, min = 1, max = 20, step = 1),
                                            uiOutput("min_unique_pep_warning"),
                                            p(class = "param-hint", "Filter proteins with unique peptides ≥ this value (integer ≥ 1). Common values: 1, 2, 3, 6.")
                                        )
                                 )
                               )
                           )
                         ),
                         
                         # Color Palette & Point Size
                         tags$details(
                           tags$summary("Color Palette & Point Size", style = "cursor: pointer; font-weight: bold; color: #2c3e50; margin-bottom: 10px;"),
                           div(style = "margin-top: 10px;",
                               div(class = "color-palette-row",
                                   div(class = "color-card",
                                       div(class = "color-card-label", "Up"),
                                       colourpicker::colourInput("color_up", NULL, value = "#FF0000", showColour = "background", allowTransparent = FALSE),
                                       div(class = "color-card-value", textOutput("val_up", inline = TRUE))
                                   ),
                                   div(class = "color-card",
                                       div(class = "color-card-label", "Down"),
                                       colourpicker::colourInput("color_down", NULL, value = "#0000FF", showColour = "background", allowTransparent = FALSE),
                                       div(class = "color-card-value", textOutput("val_down", inline = TRUE))
                                   ),
                                   div(class = "color-card",
                                       div(class = "color-card-label", "Increase"),
                                       colourpicker::colourInput("color_increase", NULL, value = "#C00000", showColour = "background", allowTransparent = FALSE),
                                       div(class = "color-card-value", textOutput("val_inc", inline = TRUE))
                                   ),
                                   div(class = "color-card",
                                       div(class = "color-card-label", "Decrease"),
                                       colourpicker::colourInput("color_decrease", NULL, value = "#0945A5", showColour = "background", allowTransparent = FALSE),
                                       div(class = "color-card-value", textOutput("val_dec", inline = TRUE))
                                   ),
                                   div(class = "color-card",
                                       div(class = "color-card-label", "NS"),
                                       colourpicker::colourInput("color_ns", NULL, value = "#7f7e83", showColour = "background", allowTransparent = FALSE),
                                       div(class = "color-card-value", textOutput("val_ns", inline = TRUE))
                                   ),
                                   div(class = "reset-btn-wrapper",
                                       actionButton("reset_color", icon("undo"), class = "btn-circle-modern")
                                   )
                               ),
                               uiOutput("color_preview"),
                               numericInput("point_size", "Point Size", value = 4, min = 0.5, max = 10, step = 0.1)
                           )
                         ),
                         
                         hr(),
                         fluidRow(
                           column(8, selectInput("selected_comparison", "Select Comparison", choices = NULL, width = "100%")),
                           column(4, textInputMax("single_plot_title", "Custom Title", value = "", placeholder = "Auto-generated", maxlength = 50, width = "100%"))
                         ),
                         fluidRow(column(12, uiOutput("plot_info_ui"))),
                         fluidRow(
                           column(12,
                                  shinycssloaders::withSpinner(plotlyOutput("volcano_plot", height = "700px"), type = 4, color = "#3498db")
                           )
                         ),
                         
                         # Export
                         tags$details(
                           tags$summary("Export", style = "cursor: pointer; font-weight: bold; color: #2c3e50; margin-top: 10px; margin-bottom: 10px;"),
                           div(style = "margin-top: 10px;",
                               fluidRow(
                                 column(6,
                                        div(style = "background: #f8f9fa; padding: 20px; border-radius: 10px;",
                                            h5(icon("image"), " Single Plot Export"),
                                            selectInput("plot_format", "Image Format", choices = c("SVG" = "svg", "TIFF" = "tiff"), selected = "svg"),
                                            fluidRow(
                                              column(12,
                                                     div(class="input-row-with-reset",
                                                         div(class="form-group shiny-input-container",
                                                             tags$label(class="control-label", "Width (inch)"),
                                                             tags$input(type="text", id="plot_width", class="form-control", value="10", placeholder="Enter width (5-30)")
                                                         ),
                                                         actionButton("reset_plot_size", icon("undo"), class = "btn btn-outline-secondary btn-reset-small", title = "Reset to default size")
                                                     ),
                                                     div(id="plot_width_warning", class="input-warning")
                                              )
                                            ),
                                            fluidRow(
                                              column(12,
                                                     div(class="form-group shiny-input-container",
                                                         tags$label(class="control-label", "Height (inch)"),
                                                         tags$input(type="text", id="plot_height", class="form-control", value="8", placeholder="Enter height (5-30)")
                                                     ),
                                                     div(id="plot_height_warning", class="input-warning")
                                              )
                                            ),
                                            br(),
                                            textInputMax("download_single_title", "Custom Title", value = "", placeholder = "Use plot page title", maxlength = 25),
                                            downloadButton("download_plot", "Download Single Plot",
                                                           style = "width: 100%; background: #3498db; color: white; padding: 10px; margin-top: 10px;")
                                        )
                                 ),
                                 column(6,
                                        div(style = "background: #f8f9fa; padding: 20px; border-radius: 10px;",
                                            h5(icon("th"), " Combined Plot Export"),
                                            textInputMax("combined_plot_title", "Main Title", value = "Combined Volcano Plots", maxlength = 30, width = "100%"),
                                            hr(),
                                            h5(icon("edit"), " Sub-plot Titles"),
                                            p(class = "param-hint", "Customize titles for each sub-plot according to the order of comparisons. Leave blank for defaults."),
                                            uiOutput("subplot_titles_ui"),
                                            div(style = "margin: 10px 0;",
                                                actionButton("goto_comparisons", "Go to Set Comparisons", icon = icon("arrow-left"), class = "btn-sm btn-outline-danger"),
                                                span(class = "red-text", "Reorder comparisons in 'Set Comparisons' tab if needed.")
                                            ),
                                            downloadButton("download_combined_plot", "Download Combined Plots",
                                                           style = "width: 100%; background: #9b59b6; color: white; padding: 10px; margin-top: 10px;")
                                        )
                                 )
                               ),
                               hr(),
                               fluidRow(
                                 column(12,
                                        div(style = "background: #f8f9fa; padding: 20px; border-radius: 10px;",
                                            h5(icon("file-excel"), " Excel Export (select comparisons)"),
                                            selectInput("export_comparisons", "Select comparisons to export", choices = NULL, multiple = TRUE),
                                            downloadButton("download_excel", "Download Selected Excel Report",
                                                           style = "width: 100%; background: #27ae60; color: white; padding: 10px;")
                                        )
                                 )
                               )
                           )
                         )
                     )
                 )
          )
        )
      ),
      # ---- UpSet 选项卡 ----
      tabPanel(
        title = "UpSet",
        value = "upset_sub",
        fluidRow(
          column(12,
                 div(class = "card-modern",
                     div(class = "card-header-modern", icon("chart-bar"), " UpSet Plot"),
                     div(style = "padding: 20px;",
                         uiOutput("venn_preprocess_steps"),
                         p("Select regulation type and at least 2 comparisons to display shared and unique proteins using UpSet plot."),
                         fluidRow(
                           column(4, selectInput("venn_type", "Regulation Type", choices = c("Up", "Down", "Increase", "Decrease"), selected = "Up")),
                           column(6, selectizeInput("venn_comparisons", "Comparisons (min. 2)", choices = NULL, multiple = TRUE, 
                                                    options = list(placeholder = 'Select at least 2 comparisons'))),
                           column(2, br(), actionButton("generate_venn", "Generate", class = "btn btn-primary btn-block"))
                         ),
                         shinycssloaders::withSpinner(plotOutput("upset_plot", height = "600px"), type = 4, color = "#3498db"),
                         fluidRow(
                           column(6, downloadButton("download_upset_png", "Download UpSet PNG", class = "btn btn-sm btn-outline-success"))
                         ),
                         uiOutput("venn_region_select_ui"),
                         shinycssloaders::withSpinner(DTOutput("venn_region_table"), type = 4, color = "#3498db"),
                         downloadButton("download_venn_region", "Download Region Proteins", class = "btn btn-sm btn-outline-secondary")
                     )
                 )
          )
        )
      )
    )
  )
}
message("[DEBUG] plots_ui.R fully defined - export formats updated")