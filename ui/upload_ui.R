# ui/upload_ui.R
message("[DEBUG] upload_ui.R loaded - added Sample Correlation heatmap in Data Quality")

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
                    plotOutput("dq_missing_heatmap", height = "600px"),
                    h5("Valid Values per Sample"),
                    plotOutput("sample_nonmiss_hist", height = "300px")
                )
            ),
            # ========== PCA 分析 ==========
            div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 12px; margin-bottom: 20px;",
                h3(icon("project-diagram"), " PCA Analysis (Raw Data by Group)", style = "margin: 0 0 15px 0; font-size: 18px;"),
                div(style = "background: white; color: #333; padding: 15px; border-radius: 8px;",
                    p("PCA is performed on the raw expression data after log2 transformation. Missing values are imputed by 1% quantile method. Samples are colored by their SubGroup from the uploaded sample info."),
                    plotOutput("dq_pca_plot", height = "500px")
                )
            ),
            # ========== 样本相关性热图（新增） ==========
            div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 12px; margin-bottom: 20px;",
                h3(icon("th"), " Sample Correlation Heatmap", style = "margin: 0 0 15px 0; font-size: 18px;"),
                div(style = "background: white; color: #333; padding: 15px; border-radius: 8px;",
                    p("Pearson correlation between selected samples based on raw expression data (1% quantile imputation + log2 transformation). Top 500 most variable proteins are used. Samples are colored by Group and Batch if available."),
                    plotOutput("dq_sample_cor_heatmap", height = "600px")
                )
            )
        )
      )
    )
  )
}
message("[DEBUG] upload_ui.R fully defined - Sample Correlation added")