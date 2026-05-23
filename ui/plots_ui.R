# ui/plots_ui.R

plots_ui <- function() {
  tabPanel(
    title = div(icon("chart-bar"), "Plots"),
    value = "plots",
    fluidRow(
      column(12,
             step_indicator(c("Upload Data", "Data Preprocessing", "Analyze & Export"), 3)
      )
    ),
    br(),
    # 顶部三个折叠面板：Define Groups, Set Comparisons, Parameters
    tags$details(
      tags$summary(icon("users"), " Define Groups", style = "font-size: 16px; font-weight: bold; color: #2c3e50;"),
      style = "margin-bottom: 15px;",
      fluidRow(
        column(4,
               div(style = "margin-bottom: 20px;",
                   h5(icon("plus-circle"), " Group Operations"),
                   div(class = "batch-group-row",
                       selectInput("group_level", "Grouping Level", choices = NULL, width = "180px"),
                       actionButton("batch_create_groups", "Batch Create Groups", icon = icon("cubes"), class = "btn-warning"),
                       actionButton("reset_groups", "Reset Groups", icon = icon("refresh"), class = "btn-danger")
                   ),
                   div(style = "display: flex; gap: 10px; margin-bottom: 15px;",
                       textInputMax("new_group_name", NULL, value = "", placeholder = "Enter group name", maxlength = 31, width = "200px", allowed_pattern = "[^a-zA-Z0-9 _-]"),
                       actionButton("add_group", "Add Group", icon = icon("plus"), class = "btn-primary")
                   ),
                   div(style = "display: flex; gap: 10px; margin-bottom: 15px;",
                       actionButton("auto_assign", "Auto-Assign Samples", icon = icon("magic"), class = "btn-info")
                   ),
                   hr(),
                   h5(icon("archive"), " Unassigned Samples"),
                   div(class = "sample-pool",
                       p(class = "param-hint", "Ctrl+Click / Shift+Click to multi-select, then drag into groups."),
                       uiOutput("unassigned_samples_ui")
                   )
               )
        ),
        column(8,
               h5(icon("folder"), " Groups & Samples"),
               div(id = "groups_container", style = "max-height: 65vh; overflow-y: auto;",
                   uiOutput("groups_ui")
               )
        )
      )
    ),
    tags$details(
      tags$summary(icon("exchange-alt"), " Set Comparisons", style = "font-size: 16px; font-weight: bold; color: #2c3e50;"),
      style = "margin-bottom: 15px;",
      fluidRow(
        column(12,
               p("Quickly add all pairwise comparisons using the automatic method, or manually define a custom comparison below."),
               fluidRow(
                 column(12,
                        div(style = "margin-bottom: 20px; padding: 15px; background: #f0f8ff; border-radius: 10px;",
                            h5(icon("layer-group"), " Auto Pairwise Comparisons"),
                            div(style = "display: flex; gap: 10px; align-items: flex-end; flex-wrap: wrap;",
                                selectInput("batch_ref_group", "Select Reference Control Group", choices = NULL, width = "250px"),
                                actionButton("batch_add_pairwise", "Add All Pairwise vs Selected Control", icon = icon("plus-circle"), class = "btn-info")
                            )
                        )
                 )
               ),
               div(style = "display: flex; gap: 10px; align-items: flex-end; margin-bottom: 20px; padding: 15px; background: #f0f8ff; border-radius: 10px; flex-wrap: wrap;",
                   selectInput("comp_treat", "Treatment Group", choices = NULL, width = "200px"),
                   div(style = "font-size: 20px; font-weight: bold; color: #666; padding-bottom: 10px;", "vs"),
                   selectInput("comp_ctrl", "Control Group", choices = NULL, width = "200px"),
                   textInputMax("comp_name", "Comparison Name (optional)", value = "", placeholder = "e.g., Mutant vs WT", maxlength = 50, width = "250px"),
                   actionButton("add_comparison", "Add Comparison", icon = icon("plus"), class = "btn-primary")
               ),
               hr(),
               div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;",
                   h5(icon("list"), uiOutput("comparisons_count_text")),
                   div(style = "display: flex; align-items: center; gap: 6px;",
                       actionButton("auto_sort_comparisons", "Auto-Sort", icon = icon("sort-alpha-down"), class = "btn-sm btn-outline-info"),
                       span(class = "param-hint", style = "margin: 0; color: #d9534f;", "Click to sort by natural order in group names")
                   )
               ),
               uiOutput("comparisons_list_ui"),
               div(style = "margin-top: 20px;",
                   actionButton("clear_comparisons", "Clear All Comparisons", class = "btn-sm btn-outline-danger")
               )
        )
      )
    ),
    tags$details(
      tags$summary(icon("sliders-h"), " Parameters", style = "font-size: 16px; font-weight: bold; color: #2c3e50;"),
      style = "margin-bottom: 15px;",
      fluidRow(
        column(6,
               div(style = "background: #f8f9fa; padding: 15px; border-radius: 10px; margin-bottom: 15px;",
                   h5(icon("balance-scale"), " Statistical Method"),
                   radioButtons("stat_method", "Choose statistical test for two-group comparisons:",
                                choices = c("t-test" = "t-test",
                                            "Wilcoxon rank-sum" = "wilcoxon",
                                            "limma (moderated t-test)" = "limma"),
                                selected = "t-test")
               )
        ),
        column(6,
               # Batch Correction 已移除，仅保留空列保持布局
               div(style = "background: #f8f9fa; padding: 15px; border-radius: 10px; margin-bottom: 15px;",
                   h5(icon("info-circle"), " Note"),
                   p("Batch correction can be performed in the Data Preprocessing step if needed.")
               )
        )
      ),
      fluidRow(
        column(6,
               div(style = "background: #f8f9fa; padding: 15px; border-radius: 10px; margin-bottom: 15px;",
                   h5(icon("balance-scale"), " Fold Change Thresholds"),
                   numericInput("fc_up", "FC up >", value = 1.2, min = 1, step = 0.1),
                   numericInput("fc_down", "FC down <", value = 0.84, min = 0, max = 1, step = 0.01),
                   hr(),
                   h5(icon("chart-line"), " Statistical Significance"),
                   selectInput("p_cut", "P-value threshold", choices = c("0.05", "0.1"), selected = "0.05")
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
                   numericInput("min_ctrl_valid", "Control group min valid replicates", value = 2, min = 1, max = 20),
                   numericInput("min_rep_ttest", "Min replicates for t-test", value = 2, min = 1, max = 10),
                   numericInput("min_rep_inc", "Min replicates for 'Increase'", value = 2, min = 1, max = 10),
                   numericInput("min_rep_dec", "Min replicates for 'Decrease'", value = 2, min = 1, max = 10)
               )
        )
      ),
      fluidRow(
        column(12,
               div(style = "background: #e3f2fd; padding: 15px; border-radius: 10px;",
                   h5(icon("filter"), " Protein Filtering (Unique Peptides)"),
                   numericInput("min_unique_pep", "Minimum Unique Peptides", value = 2, min = 1, max = 20, step = 1)
               )
        )
      )
    ),
    hr(),
    # 图表区域（Volcano / Venn-Upset / Heatmap 保持不变）
    tabsetPanel(
      id = "plots_subnav",
      tabPanel(
        title = "Volcano Plot",
        value = "volcano_sub",
        fluidRow(
          column(12,
                 div(class = "card-modern",
                     div(class = "card-header-modern", icon("balance-scale"), " Normalization Baseline & Effect"),
                     div(style = "padding: 20px;",
                         fluidRow(
                           column(4,
                                  selectInput("baseline_sample", "Baseline Sample", choices = c("Auto"), selected = "Auto"),
                                  p(class = "param-hint", "Select reference sample for global normalization.")
                           ),
                           column(8,
                                  shinycssloaders::withSpinner(
                                    plotlyOutput("norm_comparison_plot", height = "400px"),
                                    type = 4, color = "#3498db"
                                  )
                           )
                         )
                     )
                 ),
                 br(),
                 div(class = "card-modern",
                     div(class = "card-header-modern", 
                         icon("chart-line"), 
                         " Volcano Plot ",
                         span(style = "font-size: 14px; font-weight: normal;",
                              "(Click dot to view expression profile)")),
                     div(style = "padding: 20px;",
                         tags$details(
                           tags$summary("Color & Point Size Settings"),
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
                           fluidRow(column(4, numericInput("point_size", "Point Size", value = 1.8, min = 0.5, max = 10, step = 0.1)))
                         ),
                         hr(),
                         fluidRow(
                           column(8, selectInput("selected_comparison", "Select Comparison", choices = NULL, width = "100%")),
                           column(4, textInputMax("single_plot_title", "Custom Title", value = "", placeholder = "Auto-generated", maxlength = 50, width = "100%"))
                         ),
                         fluidRow(column(12, uiOutput("plot_info_ui"))),
                         shinycssloaders::withSpinner(plotlyOutput("volcano_plot", height = "700px"), type = 4, color = "#3498db"),
                         
                         tags$div(class = "export-details",
                                  tags$details(
                                    tags$summary(icon("download"), " Export Options"),
                                    div(style = "padding: 10px;",
                                        fluidRow(
                                          column(6,
                                                 h5("Single Plot Export"),
                                                 selectInput("plot_format", "Image Format", choices = c("PNG" = "png", "JPG" = "jpg"), selected = "png"),
                                                 fluidRow(
                                                   column(12, div(class="input-row-with-reset",
                                                                  div(class="form-group shiny-input-container",
                                                                      tags$label("Width (inch)"),
                                                                      tags$input(type="text", id="plot_width", class="form-control", value="10", placeholder="5-30")),
                                                                  actionButton("reset_plot_size", icon("undo"), class = "btn btn-outline-secondary btn-reset-small")
                                                   )),
                                                   column(12, div(class="form-group shiny-input-container",
                                                                  tags$label("Height (inch)"),
                                                                  tags$input(type="text", id="plot_height", class="form-control", value="8", placeholder="5-30")))
                                                 ),
                                                 textInputMax("download_single_title", "Custom Title", value = "", placeholder = "Use plot page title", maxlength = 25),
                                                 downloadButton("download_plot", "Download Single Plot", style = "width:100%; margin-top:10px;")
                                          ),
                                          column(6,
                                                 h5("Combined Plot Export"),
                                                 textInputMax("combined_plot_title", "Main Title", value = "Combined Volcano Plots", maxlength = 30),
                                                 h5("Sub-plot Titles"),
                                                 uiOutput("subplot_titles_ui"),
                                                 downloadButton("download_combined_plot", "Download Combined Plots", style = "width:100%; margin-top:10px;")
                                          )
                                        ),
                                        hr(),
                                        fluidRow(
                                          column(12,
                                                 h5("Excel & PDF Report Export"),
                                                 uiOutput("export_comparisons_ui"),
                                                 div(style = "margin-top: 5px;",
                                                     actionButton("select_all_export", "Select All", class = "btn-sm btn-outline-secondary"),
                                                     actionButton("deselect_all_export", "Deselect All", class = "btn-sm btn-outline-secondary")
                                                 ),
                                                 div(style = "display: flex; gap: 10px; margin-top: 10px;",
                                                     downloadButton("download_excel", "Download Excel Report", style = "background: #27ae60; color: white; flex: 1;"),
                                                     downloadButton("download_pdf_report", "Download PDF Report", style = "background: #e74c3c; color: white; flex: 1;")
                                                 )
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
      tabPanel(
        title = "Venn / UpSet",
        value = "venn_upset_sub",
        fluidRow(
          column(12,
                 div(class = "card-modern",
                     div(class = "card-header-modern", icon("chart-pie"), " Differential Proteins Summary"),
                     div(style = "padding: 20px;",
                         tags$details(
                           tags$summary("Comparison Statistics"),
                           shinycssloaders::withSpinner(DTOutput("stat_table"), type = 4, color = "#3498db")
                         ),
                         hr(),
                         h4(icon("venus-mars"), " Shared & Unique Proteins"),
                         radioButtons("venn_upset_method", "Visualization Method",
                                      choices = c("Venn Diagram (max 5)" = "venn", "UpSet Plot (max 15)" = "upset"),
                                      selected = "venn", inline = TRUE),
                         fluidRow(
                           column(4, selectInput("venn_type", "Regulation Type", choices = c("Up", "Down", "Increase", "Decrease"), selected = "Up")),
                           column(6, 
                                  conditionalPanel(
                                    condition = "input.venn_upset_method == 'venn'",
                                    selectizeInput("venn_comparisons_select", "Comparisons (min. 2, max 5)", choices = NULL, multiple = TRUE,
                                                   options = list(maxItems = 5, placeholder = 'Select 2-5 comparisons'))
                                  ),
                                  conditionalPanel(
                                    condition = "input.venn_upset_method == 'upset'",
                                    uiOutput("upset_comparisons_checkbox_ui")
                                  )
                           ),
                           column(2, br(), actionButton("generate_venn", "Generate", class = "btn btn-primary btn-block"))
                         ),
                         uiOutput("venn_message_ui"),
                         conditionalPanel(
                           condition = "input.venn_upset_method == 'venn'",
                           shinycssloaders::withSpinner(plotOutput("venn_plot", height = "600px"), type = 4, color = "#3498db"),
                           fluidRow(column(6, downloadButton("download_venn_png", "Download Venn PNG")))
                         ),
                         conditionalPanel(
                           condition = "input.venn_upset_method == 'upset'",
                           shinycssloaders::withSpinner(plotOutput("upset_plot", height = "600px"), type = 4, color = "#3498db"),
                           fluidRow(column(6, downloadButton("download_upset_png", "Download UpSet PNG")))
                         ),
                         uiOutput("venn_region_select_ui"),
                         shinycssloaders::withSpinner(DTOutput("venn_region_table"), type = 4, color = "#3498db"),
                         downloadButton("download_venn_region", "Download Region Proteins", class = "btn btn-sm btn-outline-secondary")
                     )
                 )
          )
        )
      ),
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
      )
    )
  )
}