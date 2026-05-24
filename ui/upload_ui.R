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
            # ========== 顶部总评分卡片 + 扣分点 ==========
            div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 20px; border-radius: 12px; margin-bottom: 20px;",
                fluidRow(
                  column(4, align = "center",
                         div(style = "font-size: 60px; font-weight: bold; margin-bottom: 5px;", textOutput("dq_score", inline = TRUE)),
                         div(style = "font-size: 28px;", textOutput("dq_grade", inline = TRUE)),
                         div(style = "margin-top: 10px; font-size: 14px; opacity: 0.9;", "Data Quality Score")
                  ),
                  column(8,
                         h4("Score Breakdown", style = "margin-top: 0;"),
                         fluidRow(
                           column(4, align = "center",
                                  div(style = "background: rgba(255,255,255,0.15); border-radius: 8px; padding: 10px;",
                                      div(style = "font-weight: bold; font-size: 16px;", "Missing Rate"),
                                      div(style = "font-size: 22px; font-weight: bold; margin-top: 5px;", textOutput("dq_missing_score_frac", inline = TRUE)),
                                      div(style = "font-size: 14px; margin-top: 3px;", textOutput("dq_missing_rate", inline = TRUE)),
                                      div(style = "font-size: 13px; margin-top: 5px; font-weight: bold;", textOutput("dq_missing_grade", inline = TRUE))
                                  )
                           ),
                           column(4, align = "center",
                                  div(style = "background: rgba(255,255,255,0.15); border-radius: 8px; padding: 10px;",
                                      div(style = "font-weight: bold; font-size: 16px;", "Sample Consistency"),
                                      div(style = "font-size: 22px; font-weight: bold; margin-top: 5px;", textOutput("dq_consistency_score_frac", inline = TRUE)),
                                      div(style = "font-size: 13px; margin-top: 5px; font-weight: bold;", textOutput("dq_consistency_grade", inline = TRUE))
                                  )
                           ),
                           column(4, align = "center",
                                  div(style = "background: rgba(255,255,255,0.15); border-radius: 8px; padding: 10px;",
                                      div(style = "font-weight: bold; font-size: 16px;", "Protein Quality"),
                                      div(style = "font-size: 22px; font-weight: bold; margin-top: 5px;", textOutput("dq_protein_score_frac", inline = TRUE)),
                                      div(style = "font-size: 13px; margin-top: 5px; font-weight: bold;", textOutput("dq_protein_grade", inline = TRUE))
                                  )
                           )
                         )
                  )
                ),
                div(style = "margin-top: 15px; padding: 10px; background: rgba(255,255,255,0.1); border-radius: 8px;",
                    p(style = "margin: 0; font-size: 13px;",
                      icon("info-circle"), " Industry benchmark for proteomics data: ≥80 = Excellent, 60–79 = Fair, <60 = Poor")
                ),
                div(style = "margin-top: 10px; padding: 10px; background: rgba(255,255,255,0.1); border-radius: 8px;",
                    p(style = "margin: 0; font-size: 13px;",
                      "✅ After completing the recommended preprocessing workflow (Missing Value Filtering → Minimum Intensity Filter → Missing Value Imputation → Batch Correction), the expected score can reach 85+ (Excellent level)")
                )
            ),
            # ========== 关键发现 ==========
            div(style = "margin-bottom: 20px;",
                h4(icon("search"), " Key Findings", style = "color: #2c3e50; margin-bottom: 15px;"),
                uiOutput("dq_key_findings")
            ),
            # ========== 预处理引导按钮（文案优化） ==========
            div(style = "margin-bottom: 30px; text-align: center;",
                actionButton("goto_preprocessing", "立即进入数据预处理 →", 
                             icon = icon("arrow-right"), 
                             class = "btn-warning btn-lg",
                             style = "border-radius: 30px; padding: 12px 40px; font-size: 18px; font-weight: bold;")
            ),
            # ========== 智能数据分析报告 ==========
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
                    hr(),
                    h4(icon("arrow-right"), " 下一步操作建议", style = "margin: 0 0 15px 0; font-size: 16px; color: #e67e22;"),
                    uiOutput("dq_recommendations"),
                    uiOutput("dq_special_note")
                )
            ),
            # ========== 缺失值分析图表 ==========
            hr(),
            h4(icon("exclamation-triangle"), " Missing Value Analysis"),
            div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 12px; margin-bottom: 20px;",
                h3(icon("database"), " Missing Value Analysis", style = "margin: 0 0 15px 0; font-size: 18px;"),
                div(style = "background: white; color: #333; padding: 15px; border-radius: 8px; position: relative;",
                    div(style = "position: absolute; top: 10px; right: 10px; display: flex; gap: 5px;",
                        downloadButton("download_missing_heatmap", "", icon = icon("download"), class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "Download high-quality image"),
                        actionButton("help_missing_heatmap", "", icon = icon("question-circle"), class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "View chart interpretation")
                    ),
                    plotOutput("dq_missing_heatmap", height = "400px"),
                    p(class = "text-muted", style = "font-size: 12px; margin-top: 8px;",
                      "图中可见大量缺失值；完成缺失值过滤 + 填充后，热图将变为全蓝。")
                )
            ),
            fluidRow(
              column(6,
                     div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 12px; margin-bottom: 20px;",
                         h3(icon("bar-chart"), " Valid Values per Sample", style = "margin: 0 0 15px 0; font-size: 18px;"),
                         div(style = "background: white; color: #333; padding: 15px; border-radius: 8px; position: relative;",
                             div(style = "position: absolute; top: 10px; right: 10px; display: flex; gap: 5px;",
                                 downloadButton("download_valid_bar", "", icon = icon("download"), class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "Download high-quality image"),
                                 actionButton("help_valid_bar", "", icon = icon("question-circle"), class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "View chart interpretation")
                             ),
                             plotOutput("dq_valid_values_plot", height = "300px"),
                             p(class = "text-muted", style = "font-size: 12px; margin-top: 8px;",
                               "所有样本有效值占比均低于 70% 阈值；完成缺失值过滤后，有效值占比预期提升至 95% 以上。")
                         )
                     )
              ),
              column(6,
                     div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 12px; margin-bottom: 20px;",
                         h3(icon("th"), " Missing Value Correlation", style = "margin: 0 0 15px 0; font-size: 18px;"),
                         div(style = "background: white; color: #333; padding: 15px; border-radius: 8px; position: relative;",
                             div(style = "position: absolute; top: 10px; right: 10px; display: flex; gap: 5px;",
                                 downloadButton("download_missing_cor", "", icon = icon("download"), class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "Download high-quality image"),
                                 actionButton("help_missing_cor", "", icon = icon("question-circle"), class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "View chart interpretation")
                             ),
                             plotOutput("dq_missing_cor_plot", height = "300px"),
                             p(class = "text-muted", style = "font-size: 12px; margin-top: 8px;",
                               "图中可见明显的缺失模式聚类；完成缺失值填充后，缺失模式相关性将被消除。")
                         )
                     )
              )
            ),
            hr(),
            fluidRow(
              column(6,
                     div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 12px; margin-bottom: 20px;",
                         h3(icon("chart-bar"), " Protein Intensity Distribution", style = "margin: 0 0 15px 0; font-size: 18px;"),
                         div(style = "background: white; color: #333; padding: 15px; border-radius: 8px; position: relative;",
                             div(style = "position: absolute; top: 10px; right: 10px; display: flex; gap: 5px;",
                                 downloadButton("download_intensity", "", icon = icon("download"), class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "Download high-quality image"),
                                 actionButton("help_intensity", "", icon = icon("question-circle"), class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "View chart interpretation")
                             ),
                             plotOutput("dq_intensity_dist_plot", height = "400px"),
                             p(class = "text-muted", style = "font-size: 12px; margin-top: 8px;",
                               "样本间强度分布中位数一致，数据重复性良好；完成最小强度过滤后，分布将更稳定。")
                         )
                     )
              ),
              column(6,
                     div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 12px; margin-bottom: 20px;",
                         h3(icon("th"), " Sample Correlation Heatmap", style = "margin: 0 0 15px 0; font-size: 18px;"),
                         div(style = "background: white; color: #333; padding: 15px; border-radius: 8px; position: relative;",
                             div(style = "position: absolute; top: 10px; right: 10px; display: flex; gap: 5px;",
                                 downloadButton("download_cor_heatmap", "", icon = icon("download"), class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "Download high-quality image"),
                                 actionButton("help_cor_heatmap", "", icon = icon("question-circle"), class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "View chart interpretation")
                             ),
                             plotOutput("dq_cor_heatmap", height = "500px"),
                             p(class = "text-muted", style = "font-size: 12px; margin-top: 8px;",
                               "样本间平均相关性仅为 0.555；完成缺失值填充 + 异常样本移除后，平均相关性预期提升至 0.8 以上。")
                         )
                     )
              )
            ),
            hr(),
            # ========== PCA 分析（双维度，标题精简） ==========
            h4(icon("project-diagram"), " PCA Analysis (Dual-view)"),
            fluidRow(
              column(6,
                     div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 12px; margin-bottom: 20px;",
                         h3(icon("users"), " PCA by Group", style = "margin: 0 0 15px 0; font-size: 18px;"),
                         div(style = "background: white; color: #333; padding: 15px; border-radius: 8px; position: relative;",
                             div(style = "position: absolute; top: 10px; right: 10px; display: flex; gap: 5px;",
                                 downloadButton("download_pca_group", "", icon = icon("download"), class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "Download high-quality image"),
                                 actionButton("help_pca_group", "", icon = icon("question-circle"), class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "View chart interpretation")
                             ),
                             plotOutput("dq_pca_group_plot", height = "450px"),
                             p(class = "text-muted", style = "font-size: 12px; margin-top: 8px;",
                               "异常样本已在图中标注；完成缺失值过滤 + 填充 + 异常样本移除后，组间分离将更清晰。")
                         )
                     )
              ),
              column(6,
                     div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 12px; margin-bottom: 20px;",
                         h3(icon("layer-group"), " PCA by Batch", style = "margin: 0 0 15px 0; font-size: 18px;"),
                         div(style = "background: white; color: #333; padding: 15px; border-radius: 8px; position: relative;",
                             div(style = "position: absolute; top: 10px; right: 10px; display: flex; gap: 5px;",
                                 downloadButton("download_pca_batch", "", icon = icon("download"), class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "Download high-quality image"),
                                 actionButton("help_pca_batch", "", icon = icon("question-circle"), class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "View chart interpretation")
                             ),
                             plotOutput("dq_pca_batch_plot", height = "450px"),
                             p(class = "text-muted", style = "font-size: 12px; margin-top: 8px;",
                               "异常样本已在图中标注；若聚类与批次对应，提示存在批次效应，需进行批次校正。")
                         )
                     )
              )
            )
        )
      )
    )
  )
}