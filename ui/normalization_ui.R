# ui/normalization_ui.R
message("[DEBUG] normalization_ui.R loaded")

normalization_ui <- function() {
  tabPanel(
    title = div(icon("balance-scale"), "Data Normalization"),
    value = "normalization",
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
}
message("[DEBUG] normalization_ui.R fully defined")