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
                sliderInput("max_missing_fraction", "Max allowed missing fraction (0-1)",
                            min = 0, max = 1, value = 0.5, step = 0.05),
                verbatimTextOutput("missing_filter_effect", placeholder = TRUE),
                helpText("Proteins with missing value proportion above this threshold will be removed. Set to 1 to keep all proteins."),
                hr(),
                h4("2. Minimum Intensity Filter", style = "color: #337ab7;"),
                numericInput("min_intensity", "Minimum intensity threshold",
                             value = 1e5, step = 1e4, min = 0),
                verbatimTextOutput("intensity_filter_effect", placeholder = TRUE),
                helpText("Proteins whose maximum intensity across all samples is below this value (or Inf) will be removed. Set to 0 to skip this filter."),
                hr(),
                h4("3. Missing Value Imputation", style = "color: #337ab7;"),
                selectInput("imputation_method", "Imputation method",
                            choices = c("Minimum value" = "min",
                                        "Mean" = "mean",
                                        "Median" = "median",
                                        "KNN" = "knn",
                                        "None (skip imputation)" = "none")),
                conditionalPanel(
                  condition = "input.imputation_method == 'min'",
                  numericInput("min_impute_value", "Minimum value size",
                               value = 1e-4, step = 1e-5),
                  helpText("Recommended: 1e-4 (0.0001) as a small constant.")
                ),
                conditionalPanel(
                  condition = "input.imputation_method == 'knn'",
                  helpText(style = "color: orange; font-weight: bold;",
                           "Note: KNN imputation requires the 'impute' package.\n",
                           "Please run: BiocManager::install('impute') if not installed.")
                ),
                conditionalPanel(
                  condition = "input.imputation_method == 'none'",
                  helpText("No imputation will be performed. Missing values will remain as NA (may affect downstream analysis).")
                ),
                hr(),
                h4("4. Batch Correction (optional)", style = "color: #337ab7;"),
                checkboxInput("perform_batch_correction", "Perform batch correction", value = FALSE),
                conditionalPanel(
                  condition = "input.perform_batch_correction == true",
                  selectInput("batch_column", "Select batch column", choices = NULL),
                  helpText("Select the column from sample info that contains batch identifiers."),
                  br(),
                  verbatimTextOutput("batch_info_preview")
                ),
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
                  )
                )
              )
            )
        )
    )
  )
}