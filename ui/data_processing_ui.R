# ui/data_processing_ui.R
message("[DEBUG] data_processing_ui.R loaded - removed Post-Processed Overview")

make_collapsible_comparison <- function(id, title, content_ui) {
  tags$details(
    tags$summary(
      title,
      style = "cursor: pointer; font-weight: bold; color: #2c3e50; margin-bottom: 10px; pointer-events: auto;",
      tabindex = "0"
    ),
    div(style = "margin-top: 10px; pointer-events: auto;", content_ui)
  )
}

data_processing_ui <- function() {
  tabPanel(
    title = div(icon("cogs"), "Data Processing"),
    value = "data_processing",
    br(),
    tabsetPanel(
      id = "data_processing_tabs",
      # ==================== Data Cleaning ====================
      tabPanel(
        title = "Data Cleaning",
        div(class = "card-modern",
            div(class = "card-header-modern", icon("table"), " Cleaning Summary"),
            div(style = "padding: 20px;",
                verbatimTextOutput("cleaning_summary")
            )
        ),
        div(class = "card-modern",
            div(class = "card-header-modern", icon("list"), " Removed Protein Details"),
            div(style = "padding: 20px;",
                p("The following proteins were removed (click to expand):"),
                tags$details(
                  tags$summary("Reverse hits", style = "cursor: pointer; font-weight: bold; margin-bottom: 5px;"),
                  verbatimTextOutput("cleaning_reverse_ids")
                ),
                tags$details(
                  tags$summary("CON_ contaminants", style = "cursor: pointer; font-weight: bold; margin-bottom: 5px;"),
                  verbatimTextOutput("cleaning_con_ids")
                )
            )
        )
      ),
      # ==================== Data Preprocessing ====================
      tabPanel(
        title = "Data Preprocessing",
        br(),
        tags$style(HTML("
          .sidebar-menu-item {
            display: block; padding: 10px 15px; margin: 4px 0;
            color: #2c3e50; background-color: transparent; border-radius: 8px;
            text-decoration: none; font-size: 15px; font-weight: 500;
            transition: background-color 0.2s ease, color 0.2s ease;
            cursor: pointer; border: none; width: 100%; text-align: left;
          }
          .sidebar-menu-item:hover, .sidebar-menu-item:focus {
            background-color: #e8f0fe; color: #1a73e8; outline: none;
          }
          .sidebar-menu-item:active {
            background-color: #d2e3fc;
          }
          .sidebar-menu-item.action-button {
            background: none;
          }
          .nav-link-custom {
            display: block; padding: 8px 12px; color: #2c3e50;
            text-decoration: none; border-radius: 6px; margin-bottom: 4px;
            transition: background 0.2s;
          }
          .nav-link-custom:hover {
            background: #e9ecef; color: #1a252f;
          }
          details {
            pointer-events: auto !important;
          }
          details summary {
            pointer-events: auto !important; cursor: pointer;
          }
          details summary::-webkit-details-marker {
            display: inline-block;
          }
          .data-source-note {
            background: #f0f8ff; border-left: 4px solid #3498db;
            padding: 10px; margin-bottom: 15px; border-radius: 4px;
            font-size: 14px; color: #2c3e50;
          }
        ")),
        bslib::page_sidebar(
          sidebar = bslib::sidebar(
            open = "open", width = 260,
            tags$h5("Preprocessing Steps", style = "color: #337ab7; font-weight: bold; margin-top: 0;"),
            tags$hr(),
            tags$div(
              actionLink("nav_missing_filter", "Missing Value Filter", class = "sidebar-menu-item"),
              actionLink("nav_imputation", "Missing Value Imputation", class = "sidebar-menu-item"),
              actionLink("nav_data_table", "Processed Data Table", class = "sidebar-menu-item")
            ),
            tags$hr(),
            div(
              style = "position: absolute; bottom: 20px; left: 15px; right: 15px;",
              actionButton("run_preprocessing", "Run Preprocessing", class = "btn-primary btn-block", style = "width: 100%;"),
              helpText("Click to execute missing value filter and imputation.")
            )
          ),
          div(
            id = "main_preprocessing_content",
            tabsetPanel(
              id = "main_nav_tabs",
              type = "hidden",
              # 1. Missing Value Filter
              tabPanel(
                title = "missing_filter", value = "missing_filter",
                div(class = "data-source-note",
                    icon("info-circle"), " Data source: Original expression matrix (raw data). Filtering is performed within each sample group."
                ),
                sliderInput("max_missing_fraction", "Max allowed missing fraction (0-1)",
                            min = 0, max = 1, value = 0.5, step = 0.05, ticks = TRUE),
                div(style = "margin-bottom: 10px;",
                    actionButton("preset_missing_0.3", "0.3", class = "btn-xs btn-outline-secondary"),
                    actionButton("preset_missing_0.5", "0.5", class = "btn-xs btn-outline-secondary"),
                    actionButton("preset_missing_0.7", "0.7", class = "btn-xs btn-outline-secondary")
                ),
                verbatimTextOutput("missing_filter_effect", placeholder = TRUE),
                div(style = "margin-top: 15px;",
                    downloadButton("download_missing_filter_excel", "Export Missing Filter Results (Excel)", class = "btn-success btn-block")
                ),
                hr(),
                h5("Data Summary (missing‑related)", style = "color: #2c3e50;"),
                verbatimTextOutput("missing_data_info")
              ),
              # 2. Missing Value Imputation
              tabPanel(
                title = "imputation", value = "imputation",
                div(class = "data-source-note",
                    icon("info-circle"), " Data source: After Missing Value Filter & Inf Filter.",
                    br(),
                    span("Imputation is applied to the filtered matrix. Download of imputation results requires preprocessing to be run.")
                ),
                selectInput("imputation_method", "Imputation method",
                            choices = c("k-Nearest Neighbors (KNN)" = "knn", "Probabilistic PCA (PPCA)" = "ppca",
                                        "Minimum Value (Fixed)" = "minvalue", "Quantile (e.g. 1% quantile)" = "quantile",
                                        "None (skip imputation)" = "none"),
                            selected = "quantile"),
                conditionalPanel(condition = "input.imputation_method == 'knn'", numericInput("knn_k", "KNN: number of neighbors (k)", value = 10, min = 1, max = 50, step = 1)),
                conditionalPanel(condition = "input.imputation_method == 'minvalue'", numericInput("minvalue_fixed", "Fixed minimum value", value = 1e-4, min = 0, step = 1e-5)),
                conditionalPanel(condition = "input.imputation_method == 'quantile'", numericInput("quantile_prob", "Quantile (e.g. 0.01 for 1%)", value = 0.01, min = 0.001, max = 0.5, step = 0.01)),
                conditionalPanel(condition = "input.imputation_method == 'knn' && output.preprocessing_done == true",
                                 hr(), h5("KNN Neighbor Lookup"), selectizeInput("knn_lookup_protein", "Select or type a protein ID", choices = NULL, multiple = FALSE, width = "100%"), DT::dataTableOutput("knn_lookup_table")),
                conditionalPanel(condition = "input.imputation_method == 'ppca' && output.preprocessing_done == true",
                                 hr(), h5("PPCA Visualization"),
                                 tags$div(style = "background: #f9f9f9; border-radius: 8px; padding: 10px; margin-bottom: 10px;", h6(icon("chart-line"), " Score Plot (PC1 vs PC2)"), plotOutput("ppca_score_plot", height = "300px")),
                                 tags$div(style = "background: #f9f9f9; border-radius: 8px; padding: 10px;", h6(icon("chart-bar"), " Distribution of Original vs. Imputed Values"), plotOutput("ppca_imputation_hist", height = "250px"))
                ),
                conditionalPanel(condition = "input.imputation_method == 'quantile' && output.preprocessing_done == true",
                                 hr(), h5("Quantile Imputation Visualization"),
                                 tags$div(style = "background: #f9f9f9; border-radius: 8px; padding: 10px; margin-bottom: 10px;", h6(icon("info-circle"), " How It Works"), p("For each sample column, the chosen quantile of non‑missing values replaces all missing values.")),
                                 h6(icon("bar-chart"), " Threshold per Sample"), plotOutput("quantile_threshold_plot", height = "400px"),
                                 selectInput("quantile_verify_sample", "Sample to view", choices = NULL, width = "100%"), plotOutput("quantile_distribution_plot", height = "300px"),
                                 h6(icon("info-circle"), " Threshold Position"), verbatimTextOutput("quantile_threshold_position"),
                                 h6(icon("table"), " Threshold Summary Table"), tableOutput("quantile_threshold_table")
                ),
                conditionalPanel(
                  condition = "output.preprocessing_done == true",
                  div(style = "margin-top: 15px;", downloadButton("download_imputation_excel", "Export Imputation Results (Excel)", class = "btn-success btn-block"))
                ),
                hr(),
                make_collapsible_comparison(
                  "imputation_comparison", "Imputation Comparison (Before/After)",
                  tagList(
                    conditionalPanel(
                      condition = "output.preprocessing_done == false",
                      div(style = "margin-top: 20px; color: #999; text-align: center;", icon("exclamation-triangle", "fa-3x"), h4("Please run preprocessing first to see imputation comparison."))
                    ),
                    conditionalPanel(
                      condition = "output.preprocessing_done == true",
                      conditionalPanel(
                        condition = "output.imputation_skipped == true",
                        div(style = "background: #fff3cd; border: 1px solid #ffeeba; padding: 15px; border-radius: 8px; margin-bottom: 20px;", h4(icon("exclamation-triangle"), " Imputation Skipped"), p("Missing value imputation was not performed. The data still contains NA values.")),
                        h4("Missing Value Heatmap"), shinycssloaders::withSpinner(plotOutput("missing_heatmap_skipped", height = "550px"), type = 4, color = "#3498db"),
                        hr(), h4("Valid Values per Sample"), shinycssloaders::withSpinner(plotOutput("valid_barplot_skipped", height = "350px"), type = 4, color = "#3498db"),
                        hr(), h4("Missing Value Summary"), tableOutput("missing_summary_table_skipped"),
                        hr(), h4("Downstream Analysis Risk"), div(style = "background: #e7f3ff; border-left: 4px solid #3498db; padding: 15px; border-radius: 4px;", tags$ul(tags$li("PCA, clustering require complete data matrix."), tags$li("Some tests may ignore proteins with missing values.")))
                      ),
                      conditionalPanel(
                        condition = "output.imputation_skipped == false",
                        fluidRow(column(12, h4("Imputation Statistics"), verbatimTextOutput("imputation_stats_text"), hr())),
                        fluidRow(column(12, h4("Boxplot: Before vs After Imputation"), shinycssloaders::withSpinner(plotOutput("imputation_boxplot", height = "600px"), type = 4, color = "#3498db"))),
                        hr(), fluidRow(column(12, h4("PCA: Before vs After Imputation"), shinycssloaders::withSpinner(plotOutput("imputation_pca_plot", height = "500px"), type = 4, color = "#3498db"))),
                        hr(), fluidRow(column(12, h4("Q-Q Plot: Before vs After Imputation"), shinycssloaders::withSpinner(plotOutput("imputation_qq_plot", height = "500px"), type = 4, color = "#3498db"))),
                        hr(), fluidRow(column(12, h4("Summary Statistics"), DT::dataTableOutput("imputation_summary_table"), br(), downloadButton("download_imputation_table", "Download Comparison Table", class = "btn btn-sm btn-outline-success")))
                      )
                    )
                  )
                )
              ),
              # 3. Processed Data Table
              tabPanel(
                title = "data_table", value = "data_table",
                h4("Preprocessing Steps Summary"),
                verbatimTextOutput("preprocessing_steps_summary"),
                hr(),
                h4("Processed Data Table"),
                DT::dataTableOutput("pre_processed_table")
              )
            )
          )
        )
      ),
      # ==================== Data Normalization ====================
      tabPanel(
        title = "Data Normalization",
        fluidRow(
          column(12,
                 div(style = "background: #f8f9fa; border-radius: 10px; padding: 15px; margin-bottom: 20px;",
                     h4(icon("calculator"), " Total Intensity Normalization", style = "margin-top: 0;"),
                     p("To correct for differences in total protein amounts between samples, we normalize each sample's total intensity to that of a chosen baseline sample. This step is performed automatically after preprocessing, but you can select the baseline here and preview the result."),
                     p("This normalization is applied to the preprocessed data and is used for all downstream analyses including volcano plots and heatmaps.")
                 )
          )
        ),
        fluidRow(
          column(12,
                 div(class = "card-modern",
                     div(class = "card-header-modern", icon("sliders-h"), " Normalization Settings"),
                     div(style = "padding: 20px;",
                         fluidRow(
                           column(6,
                                  selectInput("baseline_sample", "Select Baseline Sample",
                                              choices = NULL,
                                              selected = NULL,
                                              width = "100%"),
                                  helpText("The total intensity of all other samples will be scaled to match this sample.")
                           ),
                           column(6,
                                  verbatimTextOutput("norm_baseline_info")
                           )
                         )
                     )
                 )
          )
        ),
        fluidRow(
          column(12,
                 div(class = "card-modern",
                     div(class = "card-header-modern", icon("chart-bar"), " Normalization Comparison"),
                     div(style = "padding: 20px;",
                         plotlyOutput("norm_comparison_plot", height = "500px"),
                         hr(),
                         h5("Raw and Normalized Totals"),
                         DT::dataTableOutput("norm_totals_table")
                     )
                 )
          )
        )
      )
    )
  )
}
message("[DEBUG] data_processing_ui.R fully defined – Post-Processed Overview removed")