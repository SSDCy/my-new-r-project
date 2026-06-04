# ui/parameters_ui.R

parameters_ui <- function() {
  tabPanel(
    title = div(icon("sliders-h"), "Parameters"), 
    value = "setting",
    fluidRow(
      column(12, 
             div(class = "card-modern",
                 div(class = "card-header-modern", 
                     style = "display: flex; justify-content: space-between; align-items: center;",
                     div(icon("cog"), " Analysis Parameters")
                 ),
                 div(style = "padding: 20px;",
                     fluidRow(
                       column(6, 
                              div(style = "background: #f8f9fa; padding: 15px; border-radius: 10px; margin-bottom: 15px;",
                                  h5(icon("balance-scale"), " Fold Change Thresholds"),
                                  numericInput("fc_up", "FC up >", value = 1.2, min = 1, step = 0.1),
                                  uiOutput("fc_up_warning"),
                                  p(class = "param-hint", "Must be ≥ 1.0. Up-regulated proteins fold change threshold."),
                                  numericInput("fc_down", "FC down <", value = 0.84, min = 0, max = 1, step = 0.01),
                                  uiOutput("fc_down_warning"),
                                  p(class = "param-hint", "Must be between 0.0 and 1.0. Down-regulated proteins fold change threshold."),
                                  hr(), 
                                  h5(icon("chart-line"), " Statistical Significance"),
                                  selectInput("p_cut", "P-value threshold", choices = c("0.05", "0.01"), selected = "0.05"),
                                  p(class = "param-hint", "Select significance level (P-value threshold).")
                              )
                       ),
                       column(6, 
                              div(style = "background: #f8f9fa; padding: 15px; border-radius: 10px; margin-bottom: 15px;",
                                  h5(icon("check-circle"), " Valid Replicates"),
                                  div(style = "display: flex; align-items: center; gap: 10px; margin-bottom: 10px;",
                                      numericInput("replicate_fill_all", "Set All Thresholds", value = 2, min = 1, max = 10, step = 1, width = "100px"),
                                      actionButton("apply_replicate_fill", "Apply to All", class = "btn-sm btn-info")
                                  ),
                                  numericInput("min_treat_valid", "Treatment group min valid replicates", value = 2, min = 1, max = 20),
                                  uiOutput("min_treat_valid_warning"),
                                  numericInput("min_ctrl_valid", "Control group min valid replicates", value = 2, min = 1, max = 20),
                                  uiOutput("min_ctrl_valid_warning"),
                                  numericInput("min_rep_ttest", "Min replicates for t-test", value = 2, min = 1, max = 10),
                                  uiOutput("min_rep_ttest_warning"),
                                  numericInput("min_rep_inc", "Min replicates for 'Increase'", value = 2, min = 1, max = 10),
                                  uiOutput("min_rep_inc_warning"),
                                  numericInput("min_rep_dec", "Min replicates for 'Decrease'", value = 2, min = 1, max = 10),
                                  uiOutput("min_rep_dec_warning"),
                                  p(class = "param-hint", "Use the 'Set All' field above to fill all replicate thresholds at once, or adjust individually.")
                              )
                       )
                     ),
                     fluidRow(
                       column(12, 
                              div(style = "background: #e3f2fd; padding: 15px; border-radius: 10px;",
                                  h5(icon("filter"), " Protein Filtering (Unique Peptides)"),
                                  numericInput("min_unique_pep", "Minimum Unique Peptides", value = 2, min = 1, max = 20, step = 1),
                                  uiOutput("min_unique_pep_warning"),
                                  p(class = "param-hint", "Filter proteins with unique peptides ≥ this value (integer ≥ 1). Common values: 1, 2, 3, 6.")
                              )
                       )
                     )
                 )
             )
      )
    )
  )
}