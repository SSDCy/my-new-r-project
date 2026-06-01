# ui/preprocessing_ui.R

message("[DEBUG] preprocessing_ui.R loaded - added PPCA visualization panel with explanations")

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
                width = 4,
                style = "word-break: break-word; overflow-wrap: break-word;",
                h4("Preprocessing Steps (in order)", style = "margin-top: 0; color: #337ab7; font-weight: bold;"),
                hr(),
                # ========== 1. Missing Value Filter ==========
                h4("1. Missing Value Filter", style = "color: #337ab7;"),
                selectInput("missing_filter_mode", "Filter Mode",
                            choices = c("Global (all samples)" = "global",
                                        "Within Groups" = "group"),
                            selected = "global"),
                tags$details(
                  tags$summary("Mode Help", style = "cursor: pointer; color: #2c3e50; font-weight: bold;"),
                  div(style = "background: #f8f9fa; padding: 10px; border-radius: 5px; margin-top: 5px;",
                      p(strong("Global (all samples):"), " Calculate missing rate across all samples. A protein is removed if it is missing in more than the allowed fraction of ALL samples."),
                      p(strong("Within Groups:"), " Calculate missing rate within each group separately. A protein is kept if it satisfies the missing rate threshold in AT LEAST ONE group. This is useful when some groups have systematically more missing values (e.g., a treatment that reduces protein detection)."),
                      p(icon("exclamation-triangle"), strong("Note:"), " Within Groups mode requires that your sample info file contains a 'Group' column and that all samples are matched. Unmatched samples will trigger automatic fallback to Global mode."),
                      uiOutput("filter_mode_group_match_ui")
                  )
                ),
                div(style = "margin-top: 10px;",
                    verbatimTextOutput("missing_filter_comparison", placeholder = TRUE)
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
                    downloadButton("download_missing_filter_excel", "Export Missing Filter Results (Excel)",
                                   class = "btn-success btn-block")
                ),
                helpText("Proteins with missing value proportion > this threshold will be removed. In 'Within Groups' mode, a protein is retained if at least one group has missing rate ≤ threshold."),
                hr(),
                # ========== 2. Minimum Intensity Filter ==========
                h4("2. Minimum Intensity Filter", style = "color: #337ab7;"),
                numericInput("min_intensity", "Minimum intensity threshold",
                             value = 1e5, step = 1e4, min = 0),
                numericInput("min_samples_above_intensity", "At least N samples above threshold",
                             value = 1, min = 1, max = 100, step = 1),
                helpText("A protein is retained if at least this many samples have intensity above the minimum intensity threshold. Set to 1 for original behavior (max value)."),
                plotOutput("intensity_dist_plot", height = "250px"),
                verbatimTextOutput("intensity_filter_effect", placeholder = TRUE),
                div(style = "margin-top: 15px;",
                    downloadButton("download_intensity_filter_excel", "Export Intensity Filter Results (Excel)",
                                   class = "btn-success btn-block")
                ),
                helpText("Set to 0 to skip this filter. Recommended range: 1,000–10,000 to remove low-intensity background noise while preserving valid signals. Use the distribution plot above to identify a suitable cutoff."),
                hr(),
                # ========== 3. Missing Value Imputation ==========
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
                           "Default parameters: nPcs = 2, with automatic log2 transformation to improve accuracy.")
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
                # ---- KNN Neighbor Lookup (仅KNN模式且预处理完成) ----
                conditionalPanel(
                  condition = "input.imputation_method == 'knn' && output.preprocessing_done == true",
                  hr(),
                  h5("KNN Neighbor Lookup"),
                  selectizeInput("knn_lookup_protein", "Select or type a protein ID to view its neighbors",
                                 choices = NULL, multiple = FALSE, width = "100%"),
                  DT::dataTableOutput("knn_lookup_table")
                ),
                # ---- PPCA Visualization (仅PPCA模式且预处理完成) ----
                conditionalPanel(
                  condition = "input.imputation_method == 'ppca' && output.preprocessing_done == true",
                  hr(),
                  h5("PPCA Visualization", style = "color: #2c3e50;"),
                  tags$div(
                    style = "background: #f9f9f9; border-radius: 8px; padding: 10px; margin-bottom: 15px;",
                    h6(icon("chart-line"), " Score Plot (PC1 vs PC2)"),
                    p("Each dot is a sample. The distance between dots reflects how similar their overall protein expression patterns are. The red arrow points in the direction of the largest variation in the data – think of it as the main “trend” that distinguishes your samples."),
                    p("If dots of the same color (same experimental group) cluster together, it means the biological differences are stronger than random noise.", style = "font-size: 12px;"),
                    plotOutput("ppca_score_plot", height = "300px")
                  ),
                  tags$div(
                    style = "background: #f9f9f9; border-radius: 8px; padding: 10px;",
                    h6(icon("chart-bar"), " Distribution of Original vs. Imputed Values"),
                    p("Blue bars show the values that were originally present. Green bars show the values that PPCA filled in. When the two distributions overlap well (as they do now after log2‑transformation), it means the imputed values are realistic and consistent with the measured data."),
                    p("If the green bars were only on the far left (near zero), the imputation would be poor – it would indicate the algorithm couldn't learn the real data pattern.", style = "font-size: 12px;"),
                    plotOutput("ppca_imputation_hist", height = "250px")
                  )
                ),
                # ---- 导出填补结果 ----
                conditionalPanel(
                  condition = "output.preprocessing_done == true",
                  div(style = "margin-top: 15px;",
                      downloadButton("download_imputation_excel", "Export Imputation Results (Excel)",
                                     class = "btn-success btn-block")
                  )
                ),
                conditionalPanel(
                  condition = "output.preprocessing_done == false",
                  div(style = "margin-top: 15px; color: #999; font-style: italic;",
                      "Run preprocessing to enable imputation export.")
                ),
                hr(),
                # ========== 4. Batch Correction ==========
                h4("4. Batch Correction (Optional)", style = "color: #337ab7;"),
                checkboxInput("perform_batch_correction", "Enable ComBat Batch Correction", value = FALSE),
                verbatimTextOutput("batch_diagnostic_message", placeholder = TRUE),
                uiOutput("batch_help_text"),
                hr(),
                actionButton("run_preprocessing", "Run Preprocessing",
                             class = "btn-primary btn-lg btn-block"),
                helpText("Click to execute all steps in the above order.")
              ),
              mainPanel(
                width = 8,
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