# ui/venn_upset_ui.R
message("[DEBUG] venn_upset_ui.R loaded")

venn_upset_ui <- function() {
  tabPanel(
    title = div(icon("venus-mars"), "Venn / UpSet"),
    value = "venn_upset",
    fluidRow(
      column(12,
             div(class = "card-modern",
                 div(class = "card-header-modern", icon("chart-pie"), " Shared & Unique Proteins"),
                 div(style = "padding: 20px;",
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
  )
}
message("[DEBUG] venn_upset_ui.R fully defined")