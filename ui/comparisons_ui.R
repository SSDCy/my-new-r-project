# ui/comparisons_ui.R

comparisons_ui <- function() {
  tabPanel(
    title = div(icon("exchange-alt"), "Set Comparisons"), 
    value = "comparisons",
    fluidRow(
      column(12, 
             step_indicator(c("Upload Data", "Data Preprocessing", "Define Groups", "Set Comparisons", "Set Parameters", "Analyze & Export"), 4),
             div(class = "card-modern",
                 div(class = "card-header-modern", icon("balance-scale"), " Define Comparisons"),
                 div(style = "padding: 20px;",
                     p("You can quickly add all pairwise comparisons using the automatic method, or manually define a custom comparison below."),
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
                     tags$details(
                       tags$summary("Manual Comparison Entry"),
                       fluidRow(
                         column(12,
                                div(style = "display: flex; gap: 10px; align-items: flex-end; margin-bottom: 20px; padding: 15px; background: #f0f8ff; border-radius: 10px; flex-wrap: wrap;",
                                    selectInput("comp_treat", "Treatment Group", choices = NULL, width = "200px"),
                                    div(style = "font-size: 20px; font-weight: bold; color: #666; padding-bottom: 10px;", "vs"),
                                    selectInput("comp_ctrl", "Control Group", choices = NULL, width = "200px"),
                                    textInputMax("comp_name", "Comparison Name (optional)", value = "", placeholder = "e.g., Mutant vs WT", maxlength = 50, width = "250px"),
                                    actionButton("add_comparison", "Add Comparison", icon = icon("plus"), class = "btn-primary")
                                )
                         )
                       )
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
      )
    )
  )
}