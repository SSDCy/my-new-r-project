# ui/home_ui.R
message("[DEBUG] home_ui.R loaded")

home_ui <- function() {
  current_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[DEBUG] home_ui: DG Lab footer updated, timestamp = %s", current_time))
  message("[DEBUG] home_ui: 'Universal' removed, stats sample cor removed, export 'High-quality' removed, QC expanded, preprocessing renamed")
  
  tabPanel(
    title = div(icon("home"), "Home"),
    value = "home",
    # ---- 顶部大标题区域（去掉 Universal） ----
    div(style = "background: linear-gradient(135deg, #2c3e50 0%, #4ca1af 100%); padding: 40px 20px; color: white; text-align: center; border-radius: 0 0 20px 20px; margin-bottom: 30px;",
        h1("Proteomics Differential Analysis Platform", style = "font-weight: bold; font-size: 2.5em;"),
        p("from raw MaxQuant output to publication‑ready figures, annotations and biological insights", style = "font-size: 1.2em; opacity: 0.9;"),
        br(),
        actionButton("home_start_btn", "Click here to start", class = "btn btn-warning btn-lg", style = "font-weight: bold; border-radius: 50px; padding: 10px 30px;")
    ),
    
    fluidRow(
      column(12,
             # ---- 平台概述 ----
             div(class = "card-modern",
                 div(class = "card-header-modern", icon("info-circle"), " Overview"),
                 div(style = "padding: 20px;",
                     p("This platform offers a complete workflow for label‑free proteomics differential analysis. ",
                       "Starting from MaxQuant proteinGroups.txt, you can perform data cleaning, missing value filtering, ",
                       "normalization, statistical testing, and generate interactive volcano plots, heatmaps, ",
                       "and sample correlation heatmaps. In addition, it supports Batch CD‑Search for conserved domain annotation ",
                       "and comprehensive Excel report export. All steps are visualized with real‑time feedback.")
                 )
             ),
             
             # ---- 功能模块网格（新顺序：数据上传与功能注释 → 质量控制 → 数据预处理 → 统计分析 → 导出） ----
             div(class = "card-modern",
                 div(class = "card-header-modern", icon("cubes"), " Key Modules"),
                 div(style = "padding: 20px;",
                     fluidRow(
                       # 数据上传与功能注释
                       column(4,
                              div(style = "background: #f8f9fa; border-radius: 12px; padding: 20px; height: 220px;",
                                  div(style = "font-size: 36px; color: #3498db; text-align: center;", icon("upload")),
                                  h5(style = "text-align: center; font-weight: bold;", "Data Upload & Annotation"),
                                  p(style = "font-size: 14px; text-align: center;", "Import MaxQuant proteinGroups.txt, sample info, and FASTA for CD‑Search annotation.")
                              )
                       ),
                       # 质量控制（补充详细内容）
                       column(4,
                              div(style = "background: #f8f9fa; border-radius: 12px; padding: 20px; height: 220px;",
                                  div(style = "font-size: 36px; color: #f39c12; text-align: center;", icon("search")),
                                  h5(style = "text-align: center; font-weight: bold;", "Quality Control"),
                                  p(style = "font-size: 14px; text-align: center;",
                                    "Missing value heatmaps · Valid Values per Sample · ",
                                    "Venn Diagrams by Treatment Group · Protein Presence Matrix · ",
                                    "Peptide Length Distribution")
                              )
                       ),
                       # 数据预处理（改名为数据预处理，包含清洗）
                       column(4,
                              div(style = "background: #f8f9fa; border-radius: 12px; padding: 20px; height: 220px;",
                                  div(style = "font-size: 36px; color: #e67e22; text-align: center;", icon("filter")),
                                  h5(style = "text-align: center; font-weight: bold;", "Data Preprocessing"),
                                  p(style = "font-size: 14px; text-align: center;",
                                    "Data cleaning · Missing value filtering · ",
                                    "Imputation (KNN / PPCA / Quantile) · Total intensity normalization")
                              )
                       )
                     ),
                     br(),
                     fluidRow(
                       # 统计分析（删除样本相关性）
                       column(4,
                              div(style = "background: #f8f9fa; border-radius: 12px; padding: 20px; height: 180px;",
                                  div(style = "font-size: 36px; color: #2ecc71; text-align: center;", icon("chart-bar")),
                                  h5(style = "text-align: center; font-weight: bold;", "Statistical Analysis"),
                                  p(style = "font-size: 14px; text-align: center;", "t‑test · Volcano plots · Heatmaps · Sample correlation")
                              )
                       ),
                       # 导出（删除 High‑quality）
                       column(4,
                              div(style = "background: #f8f9fa; border-radius: 12px; padding: 20px; height: 180px;",
                                  div(style = "font-size: 36px; color: #e74c3c; text-align: center;", icon("file-export")),
                                  h5(style = "text-align: center; font-weight: bold;", "Export & Reporting"),
                                  p(style = "font-size: 14px; text-align: center;", "Figures (SVG / TIFF) · Formatted Excel reports with all DE results and annotations")
                              )
                       ),
                       column(4)
                     )
                 )
             ),
             
             # ---- 工作流程（新顺序） ----
             div(class = "card-modern",
                 div(class = "card-header-modern", icon("route"), " Analysis Workflow"),
                 div(style = "padding: 30px; text-align: center;",
                     div(style = "display: flex; justify-content: center; align-items: center; flex-wrap: wrap; gap: 20px;",
                         # 1. 数据上传与功能注释
                         div(style = "width: 150px;",
                             div(style = "background: #3498db; color: white; border-radius: 50%; width: 80px; height: 80px; display: flex; align-items: center; justify-content: center; margin: 0 auto;",
                                 icon("upload", "fa-2x")),
                             h5("1. Upload & Annotation", style = "margin-top: 10px; font-weight: bold;"),
                             p("Import & CD‑Search", style = "font-size: 12px;")
                         ),
                         div(style = "font-size: 30px; color: #aaa;", icon("arrow-right")),
                         # 2. 质量控制
                         div(style = "width: 150px;",
                             div(style = "background: #f39c12; color: white; border-radius: 50%; width: 80px; height: 80px; display: flex; align-items: center; justify-content: center; margin: 0 auto;",
                                 icon("search", "fa-2x")),
                             h5("2. Quality Control", style = "margin-top: 10px; font-weight: bold;"),
                             p("Missing & Cor", style = "font-size: 12px;")
                         ),
                         div(style = "font-size: 30px; color: #aaa;", icon("arrow-right")),
                         # 3. 数据预处理
                         div(style = "width: 150px;",
                             div(style = "background: #e67e22; color: white; border-radius: 50%; width: 80px; height: 80px; display: flex; align-items: center; justify-content: center; margin: 0 auto;",
                                 icon("filter", "fa-2x")),
                             h5("3. Preprocessing", style = "margin-top: 10px; font-weight: bold;"),
                             p("Filter & Impute", style = "font-size: 12px;")
                         ),
                         div(style = "font-size: 30px; color: #aaa;", icon("arrow-right")),
                         # 4. 统计分析
                         div(style = "width: 150px;",
                             div(style = "background: #2ecc71; color: white; border-radius: 50%; width: 80px; height: 80px; display: flex; align-items: center; justify-content: center; margin: 0 auto;",
                                 icon("chart-bar", "fa-2x")),
                             h5("4. Statistical Analysis", style = "margin-top: 10px; font-weight: bold;"),
                             p("Volcano & Heatmap", style = "font-size: 12px;")
                         ),
                         div(style = "font-size: 30px; color: #aaa;", icon("arrow-right")),
                         # 5. 导出
                         div(style = "width: 150px;",
                             div(style = "background: #e74c3c; color: white; border-radius: 50%; width: 80px; height: 80px; display: flex; align-items: center; justify-content: center; margin: 0 auto;",
                                 icon("file-export", "fa-2x")),
                             h5("5. Export", style = "margin-top: 10px; font-weight: bold;"),
                             p("Excel & Figures", style = "font-size: 12px;")
                         )
                     )
                 )
             ),
             
             # ---- 底部信息 ----
             div(style = "margin-top: 40px; padding: 20px; background: #f8f9fa; border-radius: 12px; text-align: center;",
                 p(paste0("DG Lab · Last updated ", current_time), style = "color: #666;")
             )
      )
    )
  )
}

message("[DEBUG] home_ui.R fully defined")