# ui/upload_ui.R
message("[DEBUG] upload_ui.R loaded - eggNOG and Sequence Download removed")

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
                                       accept = c(".txt"), buttonLabel = "Browse", placeholder = "No file selected"),
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
                                       accept = c(".csv", ".txt", ".xlsx", ".xls"), buttonLabel = "Browse", placeholder = "No file selected"),
                             p(style = "color: #666; font-size: 12px;",
                               icon("exclamation-triangle"), " First column must be sample names, matching the sample names derived from intensity column headers (without the 'LFQ intensity ' or 'Intensity ' prefix).")
                         ),
                         # ========= Batch CD-Search =========
                         div(style = "background: #f8f9fa; padding: 20px; border-radius: 10px; margin-top: 15px;",
                             h4(icon("search"), " Batch CD-Search (NCBI Conserved Domains)"),
                             p(style = "color: #666; font-size: 12px;",
                               icon("info-circle"), " Search protein sequences against the NCBI Conserved Domain Database (CDD). You can either run the search directly (if network allows) or manually submit a FASTA file to the NCBI website and upload the result."),
                             fileInput("cd_search_fasta", "Upload FASTA for CD-Search (optional, can use the same reference FASTA)",
                                       accept = c(".fasta", ".fa", ".txt"), buttonLabel = "Browse", placeholder = "Select FASTA file"),
                             actionButton("run_cd_search", "Run Batch CD-Search", icon = icon("play"), class = "btn-success"),
                             br(), br(),
                             uiOutput("cd_search_status_ui"),
                             tags$hr(),
                             p(strong("Or upload CD-Search result file manually:")),
                             fileInput("cd_manual_file", "Choose CD-Search result file (TSV/TXT)",
                                       accept = c(".txt", ".tsv", ".csv"), buttonLabel = "Browse", placeholder = "No file selected"),
                             p(style = "color: #666; font-size: 12px;",
                               icon("info-circle"), " If automatic search fails, visit https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi, upload the FASTA, select 'CDD--62456 PSSMs', run the search, and download the result as a TSV file. Then upload it here.")
                         )
                  ),
                  column(6,
                         div(style = "background: #f8f9fa; padding: 20px; border-radius: 10px; margin-bottom: 15px;",
                             h4(icon("chart-line"), " Data Summary"),
                             uiOutput("data_summary_ui"),
                             br(),
                             uiOutput("detected_samples_ui"),
                             uiOutput("sample_match_hint")
                         ),
                         uiOutput("cd_status_ui")   # 只保留 CD-Search 的状态显示
                  )
                ),
                hr(),
                tags$details(tags$summary(icon("table"), " Expression Matrix Preview (Sample Columns Only)", style = "cursor: pointer; font-weight: bold; color: #2c3e50; margin-bottom: 10px;"), DT::dataTableOutput("upload_preview")),
                hr(),
                tags$details(tags$summary(icon("table"), " Sample Information Preview", style = "cursor: pointer; font-weight: bold; color: #2c3e50; margin-bottom: 10px;"), DT::dataTableOutput("sample_info_preview")),
                # eggNOG 预览已删除，只保留 CD 预览
                conditionalPanel(
                  condition = "output.cd_loaded",
                  tags$details(tags$summary(icon("table"), " CD-Search Result Preview", style = "cursor: pointer; font-weight: bold; color: #2c3e50; margin-bottom: 10px;"), DT::dataTableOutput("cd_preview"))
                )
            )
        )
      ),
      # ---- 数据质量分析 ----
      tabPanel(
        title = "Data Quality Analysis",
        div(style = "padding: 10px;",
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
                    downloadButton("download_sample_nonmiss_hist", "Download Histogram", class = "btn-sm btn-outline-success", style = "margin-bottom: 5px;"),
                    plotOutput("sample_nonmiss_hist", height = "500px")
                )
            ),
            # 共有/独有蛋白与肽段
            div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 12px; margin-bottom: 20px;",
                h3(icon("venus-mars"), " Shared & Unique Proteins (Table View)", style = "margin: 0 0 15px 0; font-size: 18px;"),
                div(style = "background: white; color: #333; padding: 15px; border-radius: 8px;",
                    p("Select samples in the Missing Heatmap above (minimum 2). The table and plots will update automatically."),
                    
                    h4("Venn Diagrams by Treatment Group"),
                    p("For each treatment group (based on sample name prefix), a Venn diagram showing overlap of detected proteins among replicates is displayed."),
                    uiOutput("group_venn_plots_ui"),
                    
                    h4("Protein Presence Matrix"),
                    p("Shows 1 if protein was detected (non-missing) in the sample, 0 otherwise. 'Sum' column = number of samples where detected. CD-Search annotation columns (if uploaded) are appended."),
                    div(style = "margin-bottom: 10px;",
                        downloadButton("download_intersection_proteins", "Download Protein Table CSV", class = "btn-sm btn-outline-success")
                    ),
                    DT::dataTableOutput("intersection_protein_table"),
                    
                    hr(),
                    tags$details(
                      tags$summary("Peptide Sequences (for all proteins above)", style = "cursor: pointer; font-weight: bold; color: #2c3e50; margin-bottom: 10px;"),
                      div(style = "margin-top: 10px;",
                          p("If 'Peptide sequences' column is present in the original data, the peptide sequences associated with each Master protein ID are displayed."),
                          radioButtons("peptide_display_mode", "Display mode:",
                                       choices = c("Merged (one protein per row)" = "merged",
                                                   "Expanded (one peptide per row)" = "expanded"),
                                       selected = "merged", inline = TRUE),
                          downloadButton("download_intersection_peptides", "Download Peptide Sequences CSV", class = "btn-sm btn-outline-secondary"),
                          DT::dataTableOutput("intersection_peptide_table")
                      )
                    ),
                    
                    hr(),
                    h4("Peptide Length Distribution"),
                    p("Distribution of peptide lengths for all detected proteins in the selected samples."),
                    downloadButton("download_peptide_length_hist", "Download Histogram PNG", class = "btn-sm btn-outline-success", style = "margin-bottom: 5px;"),
                    plotOutput("intersection_peptide_length_hist", height = "500px")
                )
            ),
            # 样本相关性热图
            div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 12px; margin-bottom: 20px;",
                h3(icon("th"), " Sample Correlation Heatmap", style = "margin: 0 0 15px 0; font-size: 18px;"),
                div(style = "background: white; color: #333; padding: 15px; border-radius: 8px;",
                    p("Pearson correlation between selected samples based on raw expression data (1% quantile imputation + log2 transformation). Top 500 most variable proteins are used. Samples are colored by SubGroup if available."),
                    downloadButton("download_dq_sample_cor_png", "Download Heatmap PNG", class = "btn-sm btn-outline-success", style = "margin-right: 5px;"),
                    downloadButton("download_dq_sample_cor_matrix", "Download Correlation Matrix CSV", class = "btn-sm btn-outline-secondary"),
                    plotOutput("dq_sample_cor_heatmap", height = "600px")
                )
            ),
            # PCA 分析
            div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 12px; margin-bottom: 20px;",
                h3(icon("project-diagram"), " PCA Analysis (Raw Data by Group)", style = "margin: 0 0 15px 0; font-size: 18px;"),
                div(style = "background: white; color: #333; padding: 15px; border-radius: 8px;",
                    p("PCA is performed on the raw expression data after log2 transformation. Missing values are imputed by 1% quantile method. Samples are colored by their SubGroup from the uploaded sample info."),
                    plotOutput("dq_pca_plot", height = "500px")
                )
            )
        )
      )
    )
  )
}
message("[DEBUG] upload_ui.R fully defined (eggNOG and Sequence Download removed)")