# server/input_validation.R

output$fc_up_warning <- renderUI({
  val <- input$fc_up
  if (is.null(val) || is.na(val)) return(div(class = "input-warning", "Cannot be empty."))
  if (val < 1.0) return(div(class = "input-warning", "Must be >= 1.0."))
  NULL
})
output$fc_down_warning <- renderUI({
  val <- input$fc_down
  if (is.null(val) || is.na(val)) return(div(class = "input-warning", "Cannot be empty."))
  if (val <= 0 || val >= 1) return(div(class = "input-warning", "Must be 0 < x < 1."))
  NULL
})
output$min_treat_valid_warning <- renderUI({
  val <- input$min_treat_valid
  if (is.null(val) || is.na(val)) return(div(class = "input-warning", "Cannot be empty."))
  if (val < 1 || val != round(val)) return(div(class = "input-warning", "Integer >= 1."))
  NULL
})
output$min_ctrl_valid_warning <- renderUI({
  val <- input$min_ctrl_valid
  if (is.null(val) || is.na(val)) return(div(class = "input-warning", "Cannot be empty."))
  if (val < 1 || val != round(val)) return(div(class = "input-warning", "Integer >= 1."))
  NULL
})
output$min_rep_ttest_warning <- renderUI({
  val <- input$min_rep_ttest
  if (is.null(val) || is.na(val)) return(div(class = "input-warning", "Cannot be empty."))
  if (val < 1 || val != round(val)) return(div(class = "input-warning", "Integer >= 1."))
  NULL
})
output$min_rep_inc_warning <- renderUI({
  val <- input$min_rep_inc
  if (is.null(val) || is.na(val)) return(div(class = "input-warning", "Cannot be empty."))
  if (val < 1 || val != round(val)) return(div(class = "input-warning", "Integer >= 1."))
  NULL
})
output$min_rep_dec_warning <- renderUI({
  val <- input$min_rep_dec
  if (is.null(val) || is.na(val)) return(div(class = "input-warning", "Cannot be empty."))
  if (val < 1 || val != round(val)) return(div(class = "input-warning", "Integer >= 1."))
  NULL
})
output$min_unique_pep_warning <- renderUI({
  val <- input$min_unique_pep
  if (is.null(val) || is.na(val)) return(div(class = "input-warning", "Cannot be empty."))
  if (val < 1 || val != round(val)) return(div(class = "input-warning", "Integer >= 1."))
  if (val > 10) return(div(class = "input-warning", "<= 10 recommended."))
  NULL
})