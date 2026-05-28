# ui/preprocessing_ui.R

preprocessing_ui <- function() {
  tabPanel(
    title = div(icon("filter"), "Data Preprocessing"),
    value = "preprocessing",
    fluidRow(
      column(12,
             step_indicator(c("Upload Data", "Data Preprocessing", "Analyze & Export"), 2)
      )
    ),
    br(),
    div(class = "card-modern",
        div(class = "card-header-modern",
            div(icon("filter"), " Data Preprocessing")
        ),
        div(style = "padding: 20px;",
            sidebarLayout(
              sidebarPanel(
                width = 3,
                h4("Preprocessing Steps (in order)", style = "margin-top: 0; color: #337ab7; font-weight: bold;"),
                hr(),
                h4("1. Missing Value Filter", style = "color: #337ab7;"),
                selectInput("missing_filter_mode", "Filter Mode",
                            choices = c("Global (all samples)" = "global",
                                        "Within Groups" = "group"),
                            selected = "global"),
                sliderInput("max_missing_fraction", "Max allowed missing fraction (0-1)",
                            min = 0, max = 1, value = 0.5, step = 0.05, ticks = TRUE),
                div(style = "margin-bottom: 10px;",
                    actionButton("preset_missing_0.3", "0.3", class = "btn-xs btn-outline-secondary"),
                    actionButton("preset_missing_0.5", "0.5", class = "btn-xs btn-outline-secondary"),
                    actionButton("preset_missing_0.7", "0.7", class = "btn-xs btn-outline-secondary")
                ),
                verbatimTextOutput("missing_filter_effect", placeholder = TRUE),
                helpText("Proteins with missing value proportion > this threshold will be removed. In 'Within Groups' mode, a protein is retained if at least one group has missing rate ≤ threshold."),
                hr(),
                h4("2. Minimum Intensity Filter", style = "color: #337ab7;"),
                numericInput("min_intensity", "Minimum intensity threshold",
                             value = 1e5, step = 1e4, min = 0),
                numericInput("min_samples_above_intensity", "At least N samples above threshold",
                             value = 1, min = 1, max = 100, step = 1),
                helpText("A protein is retained if at least this many samples have intensity above the minimum intensity threshold. Set to 1 for original behavior (max value)."),
                plotOutput("intensity_dist_plot", height = "250px"),
                verbatimTextOutput("intensity_filter_effect", placeholder = TRUE),
                helpText("Set to 0 to skip this filter. Recommended range: 1,000–10,000 to remove low-intensity background noise while preserving valid signals. Use the distribution plot above to identify a suitable cutoff."),
                hr(),
                h4("3. Missing Value Imputation", style = "color: #337ab7;"),
                selectInput("imputation_method", "Imputation method",
                            choices = c("k-Nearest Neighbors (KNN)" = "knn",
                                        "Probabilistic PCA (PPCA)" = "ppca",
                                        "Minimum Value (Fixed)" = "minvalue",
                                        "Quantile (e.g. 1% quantile)" = "quantile",
                                        "None (skip imputation)" = "none"),
                            selected = "knn"),
                conditionalPanel(
                  condition = "input.imputation_method == 'knn'",
                  numericInput("knn_k", "KNN: number of neighbors (k)", value = 10, min = 1, max = 50, step = 1),
                  helpText(style = "color: orange; font-weight: bold;",
                           "Note: KNN imputation requires the 'impute' package.\n",
                           "Please run: BiocManager::install('impute') if not installed.\n",
                           "Default k = 10. Adjust based on sample size (e.g., 3-5 for small datasets).")
                ),
                conditionalPanel(
                  condition = "input.imputation_method == 'ppca'",
                  helpText(style = "color: orange; font-weight: bold;",
                           "Note: PPCA imputation requires the 'pcaMethods' package.\n",
                           "Please run: BiocManager::install('pcaMethods') if not installed.\n",
                           "Default parameters: nPcs = 2.")
                ),
                conditionalPanel(
                  condition = "input.imputation_method == 'minvalue'",
                  numericInput("minvalue_fixed", "Fixed minimum value", value = 1e-4, min = 0, step = 1e-5),
                  helpText("Replace missing values with this constant value. Commonly used: 1e-4 or 1e-3. Suitable for MNAR (left-censored missing data).")
                ),
                conditionalPanel(
                  condition = "input.imputation_method == 'quantile'",
                  numericInput("quantile_prob", "Quantile (e.g. 0.01 for 1%)", value = 0.01, min = 0.001, max = 0.5, step = 0.01),
                  helpText("Value below which data are considered low-abundance. Commonly used: 0.01 (1%) or 0.05 (5%).")
                ),
                conditionalPanel(
                  condition = "input.imputation_method == 'none'",
                  helpText("No imputation will be performed. Missing values will remain as NA (may affect downstream analysis).")
                ),
                hr(),
                h4("4. Batch Correction (Optional)", style = "color: #337ab7;"),
                checkboxInput("perform_batch_correction", "启用 ComBat 批次校正", value = FALSE),
                verbatimTextOutput("batch_diagnostic_message", placeholder = TRUE),
                uiOutput("batch_help_text"),
                hr(),
                actionButton("run_preprocessing", "Run Preprocessing",
                             class = "btn-primary btn-lg btn-block"),
                helpText("Click to execute all steps in the above order.")
              ),
              mainPanel(
                width = 9,
                tabsetPanel(
                  id = "preprocessing_tabs",
                  tabPanel("Pre-Raw Overview", value = "pre_raw_overview",
                           h4("Data Basic Information"),
                           verbatimTextOutput("pre_raw_summary"),
                           hr(),
                           h4("Interactive Threshold Calculator"),
                           fluidRow(
                             column(6, numericInput("calc_threshold", "Enter threshold:", value = 1e5, step = 1e4)),
                             column(6, verbatimTextOutput("calc_result"))
                           ),
                           hr(),
                           h4("Missing Value Distribution"),
                           plotOutput("pre_raw_missing_plot", height = "400px")
                  ),
                  tabPanel("Post-Processed Overview", value = "pre_processed_overview",
                           h4("Data Basic Information"),
                           verbatimTextOutput("pre_processed_summary"),
                           hr(),
                           h4("Missing Value Distribution"),
                           plotOutput("pre_processed_missing_plot", height = "400px")
                  ),
                  tabPanel("Processed Data Table", value = "pre_processed_table",
                           DT::dataTableOutput("pre_processed_table")
                  ),
                  tabPanel("Filter Comparison", value = "filter_comparison",
                           conditionalPanel(
                             condition = "output.preprocessing_done == false",
                             div(style = "margin-top: 20px; color: #999; text-align: center;",
                                 icon("exclamation-triangle", "fa-3x"),
                                 h4("Please run preprocessing first to see filter comparison.")
                             )
                           ),
                           conditionalPanel(
                             condition = "output.preprocessing_done == true",
                             fluidRow(
                               column(12,
                                      h4("Boxplot: Before vs After Filtering"),
                                      shinycssloaders::withSpinner(
                                        plotOutput("filter_boxplot", height = "500px"),
                                        type = 4, color = "#3498db"
                                      )
                               )
                             ),
                             hr(),
                             fluidRow(
                               column(12,
                                      h4("PCA: Before vs After Filtering"),
                                      shinycssloaders::withSpinner(
                                        plotOutput("filter_pca_plot", height = "500px"),
                                        type = 4, color = "#3498db"
                                      )
                               )
                             ),
                             hr(),
                             fluidRow(
                               column(12,
                                      h4("Summary Statistics"),
                                      DT::dataTableOutput("filter_summary_table"),
                                      br(),
                                      downloadButton("download_filter_table", "Download Comparison Table", class = "btn btn-sm btn-outline-success")
                               )
                             )
                           )
                  ),
                  tabPanel("Imputation Comparison", value = "imputation_comparison",
                           conditionalPanel(
                             condition = "output.preprocessing_done == false",
                             div(style = "margin-top: 20px; color: #999; text-align: center;",
                                 icon("exclamation-triangle", "fa-3x"),
                                 h4("Please run preprocessing first to see imputation comparison.")
                             )
                           ),
                           conditionalPanel(
                             condition = "output.preprocessing_done == true",
                             conditionalPanel(
                               condition = "output.imputation_skipped == true",
                               div(style = "background: #fff3cd; border: 1px solid #ffeeba; padding: 15px; border-radius: 8px; margin-bottom: 20px;",
                                   h4(icon("exclamation-triangle"), " Imputation Skipped", style = "color: #856404; margin-top: 0;"),
                                   p("Missing value imputation was not performed. The data still contains NA values, and all downstream analyses will be based on the original data with missing values. Consider re-running preprocessing with an imputation method for more reliable results.")
                               ),
                               h4("Missing Value Heatmap"),
                               shinycssloaders::withSpinner(plotOutput("missing_heatmap_skipped", height = "550px"), type = 4, color = "#3498db"),
                               hr(),
                               h4("Valid Values per Sample"),
                               shinycssloaders::withSpinner(plotOutput("valid_barplot_skipped", height = "350px"), type = 4, color = "#3498db"),
                               hr(),
                               h4("Missing Value Summary"),
                               tableOutput("missing_summary_table_skipped"),
                               hr(),
                               h4("Downstream Analysis Risk"),
                               div(style = "background: #e7f3ff; border-left: 4px solid #3498db; padding: 15px; border-radius: 4px;",
                                   tags$ul(
                                     tags$li("PCA, clustering, and other multivariate methods require a complete data matrix without missing values."),
                                     tags$li("Some statistical tests (e.g., t-test) may ignore proteins with missing values, reducing statistical power."),
                                     tags$li("We strongly recommend applying a missing value filter and imputation before proceeding to downstream analysis.")
                                   )
                               )
                             ),
                             conditionalPanel(
                               condition = "output.imputation_skipped == false",
                               fluidRow(
                                 column(12,
                                        h4("Imputation Statistics"),
                                        verbatimTextOutput("imputation_stats_text"),
                                        hr()
                                 )
                               ),
                               fluidRow(
                                 column(12,
                                        h4("Boxplot: Before vs After Imputation"),
                                        shinycssloaders::withSpinner(
                                          plotOutput("imputation_boxplot", height = "600px"),
                                          type = 4, color = "#3498db"
                                        )
                                 )
                               ),
                               hr(),
                               fluidRow(
                                 column(12,
                                        h4("PCA: Before vs After Imputation"),
                                        shinycssloaders::withSpinner(
                                          plotOutput("imputation_pca_plot", height = "500px"),
                                          type = 4, color = "#3498db"
                                        )
                                 )
                               ),
                               hr(),
                               fluidRow(
                                 column(12,
                                        h4("Q-Q Plot: Before vs After Imputation"),
                                        shinycssloaders::withSpinner(
                                          plotOutput("imputation_qq_plot", height = "500px"),
                                          type = 4, color = "#3498db"
                                        )
                                 )
                               ),
                               hr(),
                               fluidRow(
                                 column(12,
                                        h4("Summary Statistics"),
                                        DT::dataTableOutput("imputation_summary_table"),
                                        br(),
                                        downloadButton("download_imputation_table", "Download Comparison Table", class = "btn btn-sm btn-outline-success")
                                 )
                               )
                             )
                           )
                  ),
                  tabPanel("Batch Correction", value = "batch_correction",
                           conditionalPanel(
                             condition = "output.preprocessing_done == false",
                             div(style = "margin-top: 20px; color: #999; text-align: center;",
                                 icon("exclamation-triangle", "fa-3x"),
                                 h4("Please run preprocessing first to see batch correction results.")
                             )
                           ),
                           conditionalPanel(
                             condition = "output.preprocessing_done == true",
                             conditionalPanel(
                               condition = "output.batch_correction_performed == false",
                               div(style = "background: #fff3cd; border: 1px solid #ffeeba; padding: 15px; border-radius: 8px; margin-bottom: 20px;",
                                   h4(icon("info-circle"), " Batch Correction Not Performed", style = "color: #856404; margin-top: 0;"),
                                   p("Batch correction was not applied in the last preprocessing run.")
                               )
                             ),
                             conditionalPanel(
                               condition = "output.batch_correction_performed == true",
                               fluidRow(
                                 column(12,
                                        h4("PCA: Before vs After Batch Correction"),
                                        shinycssloaders::withSpinner(
                                          plotOutput("batch_pca_plot", height = "600px"),
                                          type = 4, color = "#3498db"
                                        ),
                                        uiOutput("batch_pca_interpretation"),
                                        downloadButton("download_batch_pca", "Download PCA Comparison", class = "btn btn-sm btn-outline-success")
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
  )
}