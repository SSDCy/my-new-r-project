# ui/upload_ui.R

upload_ui <- function() {
  tabPanel(
    title = div(icon("upload"), "Data Upload"), 
    value = "upload",
    fluidRow(
      column(12, 
             step_indicator(c("Upload Data", "Data Preprocessing", "Analyze & Export"), 1)
      )
    ),
    br(),
    tabsetPanel(
      id = "upload_tabs",
      # ---- 上传与预览 ----
      tabPanel(
        title = "Upload & Preview",
        div(class = "card-modern",
            div(class = "card-header-modern", 
                style = "display: flex; justify-content: space-between; align-items: center;",
                div(icon("database"), " Data Import & Preview"),
                actionButton("reset_all", "Reset All", icon = icon("power-off"), class = "btn-sm btn-light")
            ),
            div(style = "padding: 20px;",
                fluidRow(
                  column(6, 
                         div(style = "background: #f0f8ff; padding: 20px; border-radius: 10px; margin-bottom: 15px;",
                             h4(icon("file-upload"), " Upload Expression Matrix"),
                             fileInput("expression_file", "Choose MaxQuant proteinGroups.txt",
                                       accept = c(".txt"), 
                                       buttonLabel = "Browse", 
                                       placeholder = "No file selected"),
                             radioButtons("intensity_type", "Intensity Type",
                                          choices = c("LFQ intensity" = "LFQ", "Intensity" = "Intensity"),
                                          selected = "LFQ", inline = TRUE),
                             p(style = "color: #666; font-size: 12px;",
                               icon("info-circle"), " Select which type of intensity columns to extract from the MaxQuant output. ",
                               "If your file contains both LFQ and Intensity columns, this choice will determine which set is used for analysis. ",
                               "LFQ intensity is recommended for label-free quantification.")
                         ),
                         div(style = "background: #f8f9fa; padding: 20px; border-radius: 10px;",
                             h4(icon("tags"), " Upload Sample Information"),
                             downloadButton("download_sample_template", "Download Sample Template", 
                                            class = "btn-success btn-block", style = "margin-bottom: 10px;"),
                             fileInput("sample_info_file", "Choose Sample Info File (CSV/TXT/Excel)",
                                       accept = c(".csv", ".txt", ".xlsx", ".xls"), 
                                       buttonLabel = "Browse", 
                                       placeholder = "No file selected"),
                             p(style = "color: #666; font-size: 12px;",
                               icon("exclamation-triangle"), " First column must be sample names, matching the sample names derived from intensity column headers (without the 'LFQ intensity ' or 'Intensity ' prefix).")
                         )
                  ),
                  column(6, 
                         div(style = "background: #f8f9fa; padding: 20px; border-radius: 10px; margin-bottom: 15px;",
                             h4(icon("chart-line"), " Data Summary"), 
                             uiOutput("data_summary_ui"),
                             br(),
                             uiOutput("detected_samples_ui"),
                             uiOutput("sample_match_hint")
                         )
                  )
                ),
                hr(),
                tags$details(
                  tags$summary(icon("table"), " Expression Matrix Preview (Sample Columns Only)", 
                               style = "cursor: pointer; font-weight: bold; color: #2c3e50; margin-bottom: 10px;"),
                  DT::dataTableOutput("upload_preview")
                ),
                hr(),
                tags$details(
                  tags$summary(icon("table"), " Sample Information Preview", 
                               style = "cursor: pointer; font-weight: bold; color: #2c3e50; margin-bottom: 10px;"),
                  DT::dataTableOutput("sample_info_preview")
                )
            )
        )
      ),
      # ---- 数据质量分析 ----
      tabPanel(
        title = "Data Quality Analysis",
        div(style = "padding: 10px;",
            div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 12px; margin-bottom: 20px;",
                h3(icon("star"), " Data Quality Score", style = "margin: 0 0 15px 0; font-size: 18px;"),
                fluidRow(
                  column(3, align = "center",
                         div(style = "font-size: 48px; font-weight: bold; margin-bottom: 5px;", textOutput("dq_score", inline = TRUE)),
                         div(style = "font-size: 24px;", textOutput("dq_grade", inline = TRUE))
                  ),
                  column(3, align = "center",
                         div(class = "card-modern", style = "background: rgba(255,255,255,0.1); padding: 15px; border-radius: 8px;",
                             h4("Missing Rate", style = "margin: 0 0 10px 0; font-size: 16px;"),
                             div(style = "font-size: 24px; font-weight: bold;", textOutput("dq_missing_rate", inline = TRUE)),
                             div(style = "font-size: 12px; opacity: 0.8;", textOutput("dq_missing_score", inline = TRUE))
                         )
                  ),
                  column(3, align = "center",
                         div(class = "card-modern", style = "background: rgba(255,255,255,0.1); padding: 15px; border-radius: 8px;",
                             h4("Sample Consistency", style = "margin: 0 0 10px 0; font-size: 16px;"),
                             div(style = "font-size: 24px; font-weight: bold;", textOutput("dq_consistency", inline = TRUE)),
                             div(style = "font-size: 12px; opacity: 0.8;", textOutput("dq_consistency_score", inline = TRUE))
                         )
                  ),
                  column(3, align = "center",
                         div(class = "card-modern", style = "background: rgba(255,255,255,0.1); padding: 15px; border-radius: 8px;",
                             h4("Protein Quality", style = "margin: 0 0 10px 0; font-size: 16px;"),
                             div(style = "font-size: 24px; font-weight: bold;", textOutput("dq_protein_quality", inline = TRUE)),
                             div(style = "font-size: 12px; opacity: 0.8;", textOutput("dq_protein_score", inline = TRUE))
                         )
                  )
                )
            ),
            
            div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 12px; margin-bottom: 20px;",
                h3(icon("lightbulb"), " 智能数据分析报告", style = "margin: 0 0 15px 0; font-size: 18px;"),
                div(style = "background: white; color: #333; padding: 15px; border-radius: 8px;",
                    div(style = "background: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 15px;",
                        h4(icon("bar-chart"), " 数据质量综合分析", style = "margin: 0 0 10px 0; font-size: 16px;"),
                        fluidRow(
                          column(3, align = "center",
                                 div(style = "font-size: 32px; font-weight: bold;", textOutput("dq_total_score", inline = TRUE)),
                                 div(style = "font-size: 18px;", textOutput("dq_total_grade", inline = TRUE))
                          ),
                          column(3, align = "center",
                                 div(style = "font-size: 14px; color: #666; margin-bottom: 5px;", "缺失率:"),
                                 div(style = "font-size: 16px; font-weight: bold;", textOutput("dq_total_missing", inline = TRUE))
                          ),
                          column(3, align = "center",
                                 div(style = "font-size: 14px; color: #666; margin-bottom: 5px;", "样本一致性:"),
                                 div(style = "font-size: 16px; font-weight: bold;", textOutput("dq_total_consistency", inline = TRUE))
                          ),
                          column(3, align = "center",
                                 div(style = "font-size: 14px; color: #666; margin-bottom: 5px;", "样本相关性:"),
                                 div(style = "font-size: 16px; font-weight: bold;", textOutput("dq_total_correlation", inline = TRUE))
                          )
                        )
                    ),
                    h4(icon("search"), " 关键发现", style = "margin: 0 0 15px 0; font-size: 16px;"),
                    uiOutput("dq_key_findings"),
                    hr(),
                    h4(icon("arrow-right"), " 下一步操作建议", style = "margin: 0 0 15px 0; font-size: 16px; color: #e67e22;"),
                    uiOutput("dq_recommendations"),
                    uiOutput("dq_special_note")
                )
            ),
            
            hr(),
            h4(icon("exclamation-triangle"), " Missing Value Analysis"),
            div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 12px; margin-bottom: 20px;",
                h3(icon("database"), " Missing Value Analysis", style = "margin: 0 0 15px 0; font-size: 18px;"),
                div(style = "background: white; color: #333; padding: 15px; border-radius: 8px;",
                    h5("Missing Value Heatmap (Blue = Present, Red = Missing)"),
                    plotOutput("dq_missing_heatmap", height = "400px")
                )
            ),
            fluidRow(
              column(6,
                     div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 12px; margin-bottom: 20px;",
                         h3(icon("bar-chart"), " Valid Values per Sample", style = "margin: 0 0 15px 0; font-size: 18px;"),
                         div(style = "background: white; color: #333; padding: 15px; border-radius: 8px;",
                             plotOutput("dq_valid_values_plot", height = "300px")
                         )
                     )
              ),
              column(6,
                     div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 12px; margin-bottom: 20px;",
                         h3(icon("th"), " Missing Value Correlation", style = "margin: 0 0 15px 0; font-size: 18px;"),
                         div(style = "background: white; color: #333; padding: 15px; border-radius: 8px;",
                             verbatimTextOutput("debug_missing_cor_status"),
                             plotOutput("dq_missing_cor_plot", height = "300px")
                         )
                     )
              )
            ),
            hr(),
            fluidRow(
              column(6,
                     div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 12px; margin-bottom: 20px;",
                         h3(icon("chart-bar"), " Intensity Distribution", style = "margin: 0 0 15px 0; font-size: 18px;"),
                         div(style = "background: white; color: #333; padding: 15px; border-radius: 8px;",
                             plotOutput("dq_intensity_dist_plot", height = "400px")
                         )
                     )
              ),
              column(6,
                     div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 12px; margin-bottom: 20px;",
                         h3(icon("th"), " Sample Correlation Heatmap", style = "margin: 0 0 15px 0; font-size: 18px;"),
                         div(style = "background: white; color: #333; padding: 15px; border-radius: 8px;",
                             verbatimTextOutput("debug_cor_heatmap_status"),
                             plotOutput("dq_cor_heatmap", height = "500px")
                         )
                     )
              )
            ),
            hr(),
            div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 12px; margin-bottom: 20px;",
                h3(icon("project-diagram"), " PCA Plot (Unsupervised Clustering)", style = "margin: 0 0 15px 0; font-size: 18px;"),
                div(style = "background: white; color: #333; padding: 15px; border-radius: 8px;",
                    plotOutput("dq_pca_plot", height = "500px")
                )
            )
        )
      )
    )
  )
}