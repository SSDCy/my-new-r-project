# ui/grouping_ui.R

grouping_ui <- function() {
  tabPanel(
    title = div(icon("users"), "Define Groups"), 
    value = "grouping",
    fluidRow(
      column(12, step_indicator(c("Upload Data", "Data Preprocessing", "Define Groups", "Set Comparisons", "Set Parameters", "Analyze & Export"), 3)),
      column(4,
             div(class = "sticky-panel",
                 div(style = "margin-bottom: 20px;",
                     h5(icon("plus-circle"), " Group Operations"),
                     div(class = "batch-group-row",
                         selectInput("group_level", "Grouping Level", choices = NULL, width = "180px"),
                         actionButton("batch_create_groups", "Batch Create Groups", icon = icon("cubes"), class = "btn-warning"),
                         actionButton("reset_groups", "Reset Groups", icon = icon("refresh"), class = "btn-danger")
                     ),
                     tags$details(
                       tags$summary("Group Management"),
                       div(style = "display: flex; gap: 10px; margin-bottom: 15px;",
                           textInputMax("new_group_name", NULL, value = "", placeholder = "Enter group name", maxlength = 31, width = "200px", allowed_pattern = "[^a-zA-Z0-9 _-]"),
                           actionButton("add_group", "Add Group", icon = icon("plus"), class = "btn-primary")
                       ),
                       div(style = "display: flex; gap: 10px; margin-bottom: 15px;",
                           actionButton("auto_assign", "Auto-Assign Samples", icon = icon("magic"), class = "btn-info")
                       ),
                       div(style = "display: flex; gap: 10px;",
                           actionButton("confirm_groups", "Confirm & Go to Comparisons", icon = icon("check"), class = "btn-success")
                       )
                     )
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
             div(id = "groups_container", style = "max-height: 70vh; overflow-y: auto;",
                 uiOutput("groups_ui")
             )
      )
    )
  )
}