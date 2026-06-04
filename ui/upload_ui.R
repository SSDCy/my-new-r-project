# ui/upload_ui.R
message("[DEBUG] upload_ui.R loaded - removed baseline_sample selector (moved to Normalization)")

upload_ui <- function() {
  tabPanel(
    title = div(icon("upload"), "Data Upload"), 
    value = "upload",
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
            # ========== 缺失值分析 ==========
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
            # ========== 蛋白缺失率分布 + 样本缺失率 ==========
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
            # ========== 缺失值相关性热图 ==========
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
            # ========== 缺失值定量统计 ==========
            fluidRow(
              column(12,
                     div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 12px; margin-bottom: 20px;",
                         h3(icon("calculator"), " Missing Value Quantitative Analysis", style = "margin: 0 0 15px 0; font-size: 18px;"),
                         div(style = "background: white; color: #333; padding: 15px; border-radius: 8px;",
                             fluidRow(
                               column(6,
                                      h4("Missing Rate by Group"),
                                      plotOutput("dq_group_missing_boxplot", height = "300px")
                               ),
                               column(6,
                                      h4("Group Statistics"),
                                      tableOutput("dq_group_missing_table"),
                                      h4("Statistical Test"),
                                      verbatimTextOutput("dq_group_missing_test")
                               )
                             )
                         )
                     )
              )
            ),
            # ========== 缺失值类型诊断 ==========
            fluidRow(
              column(12,
                     div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 12px; margin-bottom: 20px;",
                         h3(icon("diagnoses"), " Missing Value Type Diagnosis (MNAR vs MAR)", style = "margin: 0 0 15px 0; font-size: 18px;"),
                         div(style = "background: white; color: #333; padding: 15px; border-radius: 8px; position: relative;",
                             div(style = "position: absolute; top: 10px; right: 10px; display: flex; gap: 5px;",
                                 downloadButton("download_missing_type", "", icon = icon("download"), class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "Download image"),
                                 actionButton("help_missing_type", "", icon = icon("question-circle"), class = "btn-sm btn-outline-secondary") %>% tagAppendAttributes(title = "View interpretation")
                             ),
                             plotOutput("dq_missing_type_plot", height = "400px"),
                             br(),
                             h5("诊断结果"),
                             textOutput("missing_type_summary")
                         )
                     )
              )
            ),
            # ========== 蛋白强度分布 ==========
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
            # ========== PCA 分析（含数据源提示） ==========
            h4(icon("project-diagram"), " PCA Analysis (Dual-view)"),
            uiOutput("pca_data_source_note"),
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