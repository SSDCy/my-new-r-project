# ui/sample_assignment_ui.R
message("[DEBUG] sample_assignment_ui.R loaded")

sample_assignment_ui <- function() {
  tabPanel(
    title = div(icon("tags"), "Sample Assignment"),
    value = "sample_assignment",
    fluidRow(
      column(12,
             div(class = "card-modern",
                 div(class = "card-header-modern", icon("tasks"), " Assign Groups to Samples"),
                 div(style = "padding: 20px;",
                     p("Set each sample's Group (Control/Treatment) and SubGroup (e.g., 100-6, 200-12). SubGroup suggestions are extracted from sample names."),
                     actionButton("save_assignment", "Save Assignment", icon = icon("save"), class = "btn-primary"),
                     hr(),
                     DT::dataTableOutput("assignment_table")
                 )
             )
      )
    )
  )
}
message("[DEBUG] sample_assignment_ui.R fully defined")