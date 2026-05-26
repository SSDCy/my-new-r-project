# ui/upload_ui.R

message("[DEBUG] upload_ui.R loaded - removed '统计' download button from Missing Value Correlation")

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
            # ========== 顶部总评分卡片 ==========
            div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 20px; border-radius: 12px; margin-bottom: 20px;",
                fluidRow(
                  column(4, align = "center",
                         div(style = "font-size: 60px; font-weight: bold; margin-bottom: 5px;", textOutput("dq_score", inline = TRUE)),
                         div(style = "margin-top: 10px; font-size: 14px; opacity: 0.9;", "Data Quality Score")
                  ),
                  column(8,
                         h4("Score Breakdown", style = "margin-top: 0;"),
                         fluidRow(
                           column(4, align = "center",
                                  div(style = "background: rgba(255,255,255,0.15); border-radius: 8px; padding: 10px;",
                                      div(style = "font-weight: bold; font-size: 16px;", "Missing Rate"),
                                      div(style = "font-size: 22px; font-weight: bold; margin-top: 5px;", textOutput("dq_missing_score_frac", inline = TRUE)),
                                      div(style = "font-size: 14px; margin-top: 3px;", textOutput("dq_missing_rate", inline = TRUE))
                                  )
                           ),
                           column(4, align = "center",
                                  div(style = "background: rgba(255,255,255,0.15); border-radius: 8px; padding: 10px;",
                                      div(style = "font-weight: bold; font-size: 16px;", "Sample Consistency"),
                                      div(style = "font-size: 22px; font-weight: bold; margin-top: 5px;", textOutput("dq_consistency_score_frac", inline = TRUE))
                                  )
                           ),
                           column(4, align = "center",
                                  div(style = "background: rgba(255,255,255,0.15); border-radius: 8px; padding: 10px;",
                                      div(style = "font-weight: bold; font-size: 16px;", "Protein Quality"),
                                      div(style = "font-size: 22px; font-weight: bold; margin-top: 5px;", textOutput("dq_protein_score_frac", inline = TRUE))
                                  )
                           )
                         )
                  )
                )
            ),
            # ========== 第一行：缺失值分析 ==========
            hr(),
            h4(icon("exclamation-triangle"), " Missing Value Analysis"),
            div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 12px; margin-bottom: 20px;",
                h3(icon("database"), " Missing Value Analysis", style = "margin: 0 0 15px 0; font-size: 18px;"),
                div(style = "background: white; color: #333; padding: 15px; border-radius: 8px; position: relative;",
                    div(style = "position: absolute; top: 10px; right: 10px; display: flex; gap: 5px;",
                        downloadButton("download_missing_heatmap", "", icon = icon("download"), class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "Download high-quality image"),
                        downloadButton("download_missing_matrix", "Matrix", class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "Download missing matrix CSV (0/1)"),
                        actionButton("help_missing_heatmap", "", icon = icon("question-circle"), class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "View chart interpretation")
                    ),
                    div(style = "margin-top: 10px;",
                        div(style = "display: flex; align-items: center; gap: 10px; margin-bottom: 10px;",
                            actionButton("heatmap_select_all", "Select All", class = "btn-sm btn-primary", icon = icon("check-square")),
                            actionButton("heatmap_clear_all", "Clear All", class = "btn-sm btn-secondary", icon = icon("square")),
                            uiOutput("heatmap_group_buttons_ui")
                        ),
                        selectizeInput("heatmap_sample_select", label = NULL,
                                       choices = character(0),
                                       multiple = TRUE,
                                       width = "100%",
                                       options = list(plugins = list('remove_button'),
                                                      placeholder = 'Click to select samples...')
                        ),
                        div(style = "margin-top: 5px;",
                            checkboxInput("heatmap_cluster_rows", "Cluster rows", value = FALSE)
                        )
                    ),
                    plotOutput("dq_missing_heatmap", height = "400px")
                )
            ),
            # ========== 第二行：蛋白缺失率分布 + 样本缺失率 ==========
            fluidRow(
              column(6,
                     div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 12px; margin-bottom: 20px;",
                         h3(icon("chart-bar"), " Protein Missing Rate Distribution", style = "margin: 0 0 15px 0; font-size: 18px;"),
                         div(style = "background: white; color: #333; padding: 15px; border-radius: 8px; position: relative;",
                             div(style = "position: absolute; top: 10px; right: 10px; display: flex; gap: 5px;",
                                 downloadButton("download_protein_missing_hist", "", icon = icon("download"), class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "Download image")
                             ),
                             plotOutput("dq_protein_missing_hist", height = "400px")
                         )
                     )
              ),
              column(6,
                     div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 12px; margin-bottom: 20px;",
                         h3(icon("bar-chart"), " Sample Missing Rate", style = "margin: 0 0 15px 0; font-size: 18px;"),
                         div(style = "background: white; color: #333; padding: 15px; border-radius: 8px; position: relative;",
                             div(style = "position: absolute; top: 10px; right: 10px; display: flex; gap: 5px;",
                                 downloadButton("download_sample_missing_bar", "", icon = icon("download"), class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "Download image")
                             ),
                             plotOutput("dq_sample_missing_bar", height = "400px")
                         )
                     )
              )
            ),
            # ========== 第三行：缺失值相关性热图（统计按钮已删除） ==========
            fluidRow(
              column(12,
                     div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 12px; margin-bottom: 20px;",
                         h3(icon("th"), " Missing Value Correlation", style = "margin: 0 0 15px 0; font-size: 18px;"),
                         div(style = "background: white; color: #333; padding: 15px; border-radius: 8px; position: relative;",
                             div(style = "position: absolute; top: 10px; right: 10px; display: flex; gap: 5px;",
                                 downloadButton("download_missing_cor", "", icon = icon("download"), class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "Download high-quality image"),
                                 downloadButton("download_missing_cor_matrix", "CSV", class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "Download missing pattern correlation matrix CSV"),
                                 actionButton("help_missing_cor", "", icon = icon("question-circle"), class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "View chart interpretation")
                             ),
                             # 子热图控件
                             div(style = "margin-top: 10px; display: flex; align-items: center; gap: 10px;",
                                 selectizeInput("missing_cor_sample_select", "子热图样本选择",
                                                choices = character(0),
                                                multiple = TRUE,
                                                width = "60%",
                                                options = list(plugins = list('remove_button'),
                                                               placeholder = '按住 Ctrl/Shift 多选样本...')
                                 ),
                                 actionButton("missing_cor_subset_go", "生成子热图", class = "btn-sm btn-primary", icon = icon("play")),
                                 actionButton("missing_cor_reset", "重置为全局热图", class = "btn-sm btn-secondary", icon = icon("refresh"))
                             ),
                             plotOutput("dq_missing_cor_plot", height = "600px")
                         )
                     )
              )
            ),
            # ========== 第四行：样本相关性热图 ==========
            fluidRow(
              column(12,
                     div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 12px; margin-bottom: 20px;",
                         h3(icon("th"), " Sample Correlation Heatmap", style = "margin: 0 0 15px 0; font-size: 18px;"),
                         div(style = "background: white; color: #333; padding: 15px; border-radius: 8px; position: relative;",
                             div(style = "position: absolute; top: 10px; right: 10px; display: flex; gap: 5px;",
                                 downloadButton("download_cor_heatmap", "", icon = icon("download"), class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "Download high-quality image"),
                                 downloadButton("download_cor_matrix", "CSV", class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "Download correlation matrix CSV"),
                                 actionButton("help_cor_heatmap", "", icon = icon("question-circle"), class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "View chart interpretation")
                             ),
                             plotOutput("dq_cor_heatmap", height = "600px")
                         )
                     )
              )
            ),
            # ========== 第五行：蛋白强度分布 ==========
            fluidRow(
              column(12,
                     div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 12px; margin-bottom: 20px;",
                         h3(icon("chart-bar"), " Protein Intensity Distribution", style = "margin: 0 0 15px 0; font-size: 18px;"),
                         div(style = "background: white; color: #333; padding: 15px; border-radius: 8px; position: relative;",
                             div(style = "position: absolute; top: 10px; right: 10px; display: flex; gap: 5px;",
                                 downloadButton("download_intensity", "", icon = icon("download"), class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "Download high-quality image"),
                                 downloadButton("download_intensity_data", "Excel", class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "Download intensity data (multi-sheet Excel)"),
                                 actionButton("help_intensity", "", icon = icon("question-circle"), class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "View chart interpretation")
                             ),
                             plotOutput("dq_intensity_dist_plot", height = "400px")
                         )
                     )
              )
            ),
            # ========== 第六行：PCA 分析 ==========
            h4(icon("project-diagram"), " PCA Analysis (Dual-view)"),
            fluidRow(
              column(6,
                     div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 12px; margin-bottom: 20px;",
                         h3(icon("users"), " PCA by Group", style = "margin: 0 0 15px 0; font-size: 18px;"),
                         div(style = "background: white; color: #333; padding: 15px; border-radius: 8px; position: relative;",
                             div(style = "position: absolute; top: 10px; right: 10px; display: flex; gap: 5px;",
                                 downloadButton("download_pca_group", "", icon = icon("download"), class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "Download high-quality image"),
                                 downloadButton("download_pca_group_data", "Excel", class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "Download PCA Group data (multi-sheet Excel)"),
                                 actionButton("help_pca_group", "", icon = icon("question-circle"), class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "View chart interpretation")
                             ),
                             plotOutput("dq_pca_group_plot", height = "400px")
                         )
                     )
              ),
              column(6,
                     div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 12px; margin-bottom: 20px;",
                         h3(icon("layer-group"), " PCA by Batch", style = "margin: 0 0 15px 0; font-size: 18px;"),
                         div(style = "background: white; color: #333; padding: 15px; border-radius: 8px; position: relative;",
                             div(style = "position: absolute; top: 10px; right: 10px; display: flex; gap: 5px;",
                                 downloadButton("download_pca_batch", "", icon = icon("download"), class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "Download high-quality image"),
                                 downloadButton("download_pca_batch_data", "Excel", class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "Download PCA Batch data (multi-sheet Excel)"),
                                 actionButton("help_pca_batch", "", icon = icon("question-circle"), class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "View chart interpretation")
                             ),
                             plotOutput("dq_pca_batch_plot", height = "400px")
                         )
                     )
              )
            )
        )
      )
    )
  )
}