# ui/cleaning_ui.R
message("[DEBUG] cleaning_ui.R loaded")

cleaning_ui <- function() {
  tabPanel(
    title = div(icon("broom"), "Data Cleaning"),
    value = "cleaning",
    fluidRow(
      column(12,
             div(style = "background: #f8f9fa; border-radius: 10px; padding: 15px; margin-bottom: 20px;",
                 h4(icon("filter"), " Automatic Data Cleaning", style = "margin-top: 0;"),
                 p("Upon uploading, the following filters are applied automatically to remove unreliable identifications:",
                   tags$ul(
                     tags$li(tags$b("Reverse hits:"), " proteins matching the reverse database (column 'Reverse' equals '+')."),
                     tags$li(tags$b("CON_ contaminants:"), " proteins with IDs starting with 'CON_' (common contaminants).")
                   )
                 )
             )
      )
    ),
    fluidRow(
      column(12,
             div(class = "card-modern",
                 div(class = "card-header-modern", icon("table"), " Cleaning Summary"),
                 div(style = "padding: 20px;",
                     verbatimTextOutput("cleaning_summary")
                 )
             )
      )
    ),
    fluidRow(
      column(12,
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
      )
    )
  )
}

message("[DEBUG] cleaning_ui.R fully defined")