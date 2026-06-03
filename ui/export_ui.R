# ui/export_ui.R
message("[DEBUG] export_ui.R loaded")

export_ui <- function() {
  tabPanel(
    title = div(icon("download"), "Export"),
    value = "export",
    fluidRow(
      column(12,
             div(class = "card-modern",
                 div(class = "card-header-modern", icon("file-download"), " Results Export"),
                 div(style = "padding: 20px;",
                     fluidRow(
                       column(6,
                              div(style = "background: #f8f9fa; padding: 20px; border-radius: 10px;",
                                  h5(icon("image"), " Single Plot Export"),
                                  selectInput("plot_format", "Image Format", choices = c("PNG" = "png", "JPG" = "jpg"), selected = "png"),
                                  fluidRow(
                                    column(12,
                                           div(class="input-row-with-reset",
                                               div(class="form-group shiny-input-container",
                                                   tags$label(class="control-label", "Width (inch)"),
                                                   tags$input(type="text", id="plot_width", class="form-control", value="10", placeholder="Enter width (5-30)")
                                               ),
                                               actionButton("reset_plot_size", icon("undo"), class = "btn btn-outline-secondary btn-reset-small", title = "Reset to default size")
                                           ),
                                           div(id="plot_width_warning", class="input-warning")
                                    )
                                  ),
                                  fluidRow(
                                    column(12,
                                           div(class="form-group shiny-input-container",
                                               tags$label(class="control-label", "Height (inch)"),
                                               tags$input(type="text", id="plot_height", class="form-control", value="8", placeholder="Enter height (5-30)")
                                           ),
                                           div(id="plot_height_warning", class="input-warning")
                                    )
                                  ),
                                  br(),
                                  textInputMax("download_single_title", "Custom Title", value = "", placeholder = "Use plot page title", maxlength = 25),
                                  downloadButton("download_plot", "Download Single Plot",
                                                 style = "width: 100%; background: #3498db; color: white; padding: 10px; margin-top: 10px;")
                              )
                       ),
                       column(6,
                              div(style = "background: #f8f9fa; padding: 20px; border-radius: 10px;",
                                  h5(icon("th"), " Combined Plot Export"),
                                  textInputMax("combined_plot_title", "Main Title", value = "Combined Volcano Plots", maxlength = 30, width = "100%"),
                                  hr(),
                                  h5(icon("edit"), " Sub-plot Titles"),
                                  p(class = "param-hint", "Customize titles for each sub-plot according to the order of comparisons. Leave blank for defaults."),
                                  uiOutput("subplot_titles_ui"),
                                  div(style = "margin: 10px 0;",
                                      actionButton("goto_comparisons", "Go to Set Comparisons", icon = icon("arrow-left"), class = "btn-sm btn-outline-danger"),
                                      span(class = "red-text", "Reorder comparisons in 'Set Comparisons' tab if needed.")
                                  ),
                                  downloadButton("download_combined_plot", "Download Combined Plots",
                                                 style = "width: 100%; background: #9b59b6; color: white; padding: 10px; margin-top: 10px;")
                              )
                       )
                     ),
                     hr(),
                     fluidRow(
                       column(12,
                              div(style = "background: #f8f9fa; padding: 20px; border-radius: 10px;",
                                  h5(icon("file-excel"), " Excel Export (select comparisons)"),
                                  selectInput("export_comparisons", "Select comparisons to export", choices = NULL, multiple = TRUE),
                                  downloadButton("download_excel", "Download Selected Excel Report",
                                                 style = "width: 100%; background: #27ae60; color: white; padding: 10px;"),
                                  hr(),
                                  h5(icon("file-pdf"), " PDF Report (select comparisons)"),
                                  p("Uses the same comparison selection above."),
                                  downloadButton("download_pdf_report", "Download PDF Report",
                                                 style = "width: 100%; background: #e74c3c; color: white; padding: 10px;")
                              )
                       )
                     )
                 )
             )
      )
    )
  )
}
message("[DEBUG] export_ui.R fully defined")