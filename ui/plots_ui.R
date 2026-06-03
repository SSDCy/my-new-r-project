# ui/plots_ui.R

plots_ui <- function() {
  tabPanel(
    title = div(icon("chart-bar"), "Plots"),
    value = "plots",
    fluidRow(
      column(12,
             step_indicator(c("Upload Data", "Data Preprocessing", "Analyze & Export"), 3)
      )
    ),
    br(),
    # 整合后的分析设置面板
    div(class = "card-modern",
        div(class = "card-header-modern", icon("cogs"), " Analysis Setup"),
        div(style = "padding: 20px;",
            fluidRow(
              column(6,
                     # ---------- 分组区域 ----------
                     h5(icon("users"), " Define Groups", style = "margin-top: 0; color: #2c3e50;"),
                     div(class = "batch-group-row",
                         selectInput("group_level", "Grouping Level", choices = NULL, width = "180px"),
                         actionButton("batch_create_groups", "Batch Create", icon = icon("cubes"), class = "btn-warning"),
                         actionButton("reset_groups", "Reset", icon = icon("refresh"), class = "btn-danger")
                     ),
                     # 折叠：添加分组和自动分配
                     tags$details(
                       tags$summary("Add / Auto-Assign", style = "cursor: pointer; font-weight: bold; color: #2c3e50; margin-bottom: 10px;"),
                       div(style = "display: flex; gap: 10px; margin-bottom: 15px;",
                           textInputMax("new_group_name", NULL, value = "", placeholder = "Enter group name", maxlength = 31, width = "200px", allowed_pattern = "[^a-zA-Z0-9 _-]"),
                           actionButton("add_group", "Add Group", icon = icon("plus"), class = "btn-primary")
                       ),
                       div(style = "margin-bottom: 15px;",
                           actionButton("auto_assign", "Auto-Assign Samples", icon = icon("magic"), class = "btn-info")
                       )
                     ),
                     hr(),
                     h5(icon("archive"), " Unassigned Samples"),
                     div(class = "sample-pool",
                         p(class = "param-hint", "Ctrl+Click / Shift+Click to multi-select, then drag into groups."),
                         uiOutput("unassigned_samples_ui")
                     ),
                     br(),
                     h5(icon("folder"), " Groups & Samples"),
                     div(id = "groups_container", style = "max-height: 50vh; overflow-y: auto;",
                         uiOutput("groups_ui")
                     )
              ),
              column(6,
                     # ---------- 比较区域 ----------
                     h5(icon("exchange-alt"), " Set Comparisons", style = "margin-top: 0; color: #2c3e50;"),
                     div(style = "margin-bottom: 20px; padding: 15px; background: #f0f8ff; border-radius: 10px;",
                         h5(icon("layer-group"), " Auto Pairwise Comparisons"),
                         div(style = "display: flex; gap: 10px; align-items: flex-end; flex-wrap: wrap;",
                             selectInput("batch_ref_group", "Reference Control Group", choices = NULL, width = "200px"),
                             actionButton("batch_add_pairwise", "Add All Pairwise", icon = icon("plus-circle"), class = "btn-info")
                         )
                     ),
                     # 折叠：手动添加比较
                     tags$details(
                       tags$summary("Manual Comparison Entry", style = "cursor: pointer; font-weight: bold; color: #2c3e50; margin-bottom: 10px;"),
                       div(style = "margin-bottom: 20px; padding: 15px; background: #f0f8ff; border-radius: 10px;",
                           div(style = "display: flex; gap: 10px; align-items: flex-end; flex-wrap: wrap; margin-bottom: 10px;",
                               selectInput("comp_treat", "Treatment Group", choices = NULL, width = "150px"),
                               div(style = "font-size: 20px; font-weight: bold; color: #666; padding-bottom: 10px;", "vs"),
                               selectInput("comp_ctrl", "Control Group", choices = NULL, width = "150px")
                           ),
                           textInputMax("comp_name", "Comparison Name (optional)", value = "", placeholder = "e.g., Mutant vs WT", maxlength = 50, width = "100%"),
                           div(style = "margin-top: 10px;",
                               actionButton("add_comparison", "Add Comparison", icon = icon("plus"), class = "btn-primary")
                           )
                       )
                     ),
                     hr(),
                     div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;",
                         h5(icon("list"), uiOutput("comparisons_count_text"), style = "margin: 0;"),
                         div(style = "display: flex; align-items: center; gap: 6px;",
                             actionButton("auto_sort_comparisons", "Auto-Sort", icon = icon("sort-alpha-down"), class = "btn-sm btn-outline-info"),
                             actionButton("clear_comparisons", "Clear All", class = "btn-sm btn-outline-danger")
                         )
                     ),
                     # 带滚动条的对比列表容器
                     div(style = "max-height: 300px; overflow-y: auto; border: 1px solid #dee2e6; border-radius: 8px; padding: 10px; background: #fff;",
                         uiOutput("comparisons_list_ui")
                     )
              )
            ),
            hr(),
            # ---------- 参数区域（折叠） ----------
            fluidRow(
              column(12,
                     tags$details(
                       tags$summary(icon("sliders-h"), " Statistical Method, Fold Change & Replicates", style = "cursor: pointer; font-weight: bold; font-size: 16px; color: #2c3e50; margin-bottom: 15px;"),
                       fluidRow(
                         column(4,
                                div(style = "background: #f8f9fa; padding: 15px; border-radius: 10px; margin-bottom: 15px;",
                                    h5(icon("balance-scale"), " Statistical Method"),
                                    radioButtons("stat_method", NULL,
                                                 choices = c("t-test" = "t-test",
                                                             "Wilcoxon rank-sum" = "wilcoxon",
                                                             "limma (moderated t-test)" = "limma"),
                                                 selected = "t-test")
                                )
                         ),
                         column(4,
                                div(style = "background: #f8f9fa; padding: 15px; border-radius: 10px; margin-bottom: 15px;",
                                    h5(icon("balance-scale"), " Fold Change & Significance"),
                                    numericInput("fc_up", "FC up >", value = 1.2, min = 1, step = 0.1),
                                    numericInput("fc_down", "FC down <", value = 0.84, min = 0, max = 1, step = 0.01),
                                    selectInput("p_cut", "P-value threshold", choices = c("0.05", "0.1"), selected = "0.05")
                                )
                         ),
                         column(4,
                                div(style = "background: #f8f9fa; padding: 15px; border-radius: 10px; margin-bottom: 15px;",
                                    h5(icon("check-circle"), " Valid Replicates"),
                                    div(style = "display: flex; align-items: center; gap: 10px; margin-bottom: 10px;",
                                        numericInput("replicate_fill_all", "Set All", value = 2, min = 1, max = 10, step = 1, width = "80px"),
                                        actionButton("apply_replicate_fill", "Apply to All", class = "btn-sm btn-info")
                                    ),
                                    fluidRow(
                                      column(6, numericInput("min_treat_valid", "Treat min", value = 2, min = 1, max = 20)),
                                      column(6, numericInput("min_ctrl_valid", "Ctrl min", value = 2, min = 1, max = 20))
                                    ),
                                    fluidRow(
                                      column(6, numericInput("min_rep_ttest", "t-test min", value = 2, min = 1, max = 10)),
                                      column(6, numericInput("min_rep_inc", "Increase min", value = 2, min = 1, max = 10))
                                    ),
                                    numericInput("min_rep_dec", "Decrease min", value = 2, min = 1, max = 10, width = "100%")
                                )
                         )
                       )
                     )
              )
            ),
            fluidRow(
              column(12,
                     div(style = "background: #e3f2fd; padding: 15px; border-radius: 10px;",
                         h5(icon("filter"), " Protein Filtering (Unique Peptides)"),
                         numericInput("min_unique_pep", "Minimum Unique Peptides", value = 2, min = 1, max = 20, step = 1)
                     )
              )
            )
        )
    ),
    hr(),
    # 图表区域（仅保留 Heatmap）
    tabsetPanel(
      id = "plots_subnav",
      tabPanel(
        title = "Heatmap",
        value = "heatmap_sub",
        fluidRow(
          column(12,
                 div(class = "card-modern",
                     div(class = "card-header-modern", icon("th"), " Differential Protein Heatmap"),
                     div(style = "padding: 20px;",
                         p("Expression patterns of proteins. Per-protein log2 + Z-score normalization is applied, followed by hierarchical clustering."),
                         fluidRow(
                           column(4,
                                  radioButtons("heatmap_data_source", "Data Source",
                                               choices = c("LFQ Intensity (per-row Z-score)" = "LFQ",
                                                           "Intensity (per-row Z-score)" = "Intensity"),
                                               selected = "LFQ"),
                                  # 新增：数据源信息显示
                                  verbatimTextOutput("heatmap_data_source_info"),
                                  hr(),
                                  radioButtons("heatmap_protein_mode", "Protein Selection Mode",
                                               choices = c("Top N proteins (by variance)" = "top_n",
                                                           "Custom protein list" = "custom"),
                                               selected = "top_n"),
                                  conditionalPanel(
                                    condition = "input.heatmap_protein_mode == 'top_n'",
                                    numericInput("heatmap_top_n", "Top N proteins", value = 20, min = 5, max = 200, step = 5)
                                  ),
                                  conditionalPanel(
                                    condition = "input.heatmap_protein_mode == 'custom'",
                                    textAreaInput("heatmap_custom_ids", "Enter Master Protein IDs", rows = 5, placeholder = "P12345\nP67890")
                                  ),
                                  hr(),
                                  conditionalPanel(
                                    condition = "input.heatmap_protein_mode == 'top_n'",
                                    selectInput("heatmap_groups", "Select Groups", choices = NULL, multiple = TRUE)
                                  ),
                                  conditionalPanel(
                                    condition = "input.heatmap_data_source == 'Intensity'",
                                    selectInput("heatmap_group_level", "Grouping Level", choices = NULL),
                                    actionButton("heatmap_apply_grouping", "Apply Grouping"),
                                    uiOutput("heatmap_group_selection_ui")
                                  ),
                                  checkboxInput("heatmap_show_sample_names", "Show sample names", value = TRUE),
                                  actionButton("generate_heatmap", "Generate Heatmap", class = "btn btn-primary btn-block"),
                                  hr(),
                                  downloadButton("download_heatmap_png", "Download Heatmap PNG", class = "btn btn-sm btn-outline-success")
                           ),
                           column(8,
                                  shinycssloaders::withSpinner(plotOutput("heatmap_plot", height = "700px"), type = 4, color = "#e67e22")
                           )
                         )
                     )
                 )
          )
        )
      )
    )
  )
}