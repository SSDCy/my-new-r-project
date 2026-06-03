# server/preprocessing_nav.R
message("[DEBUG] preprocessing_nav.R loaded - handling sidebar navigation")

# 验证 session 是否可用
shiny::observe({
  message("[DEBUG] preprocessing_nav: initial session check - session exists")
})

# 监听每个侧边栏菜单项的点击
observeEvent(input$nav_missing_filter, {
  message("[DEBUG] nav: Missing Value Filter clicked")
  updateTabsetPanel(session, "main_nav_tabs", selected = "missing_filter")
})
observeEvent(input$nav_intensity_filter, {
  message("[DEBUG] nav: Minimum Intensity Filter clicked")
  updateTabsetPanel(session, "main_nav_tabs", selected = "intensity_filter")
})
observeEvent(input$nav_imputation, {
  message("[DEBUG] nav: Missing Value Imputation clicked")
  updateTabsetPanel(session, "main_nav_tabs", selected = "imputation")
})
observeEvent(input$nav_batch, {
  message("[DEBUG] nav: Batch Correction clicked")
  updateTabsetPanel(session, "main_nav_tabs", selected = "batch")
})
observeEvent(input$nav_processed_overview, {
  message("[DEBUG] nav: Post-Processed Overview clicked")
  updateTabsetPanel(session, "main_nav_tabs", selected = "processed_overview")
})
observeEvent(input$nav_data_table, {
  message("[DEBUG] nav: Processed Data Table clicked")
  updateTabsetPanel(session, "main_nav_tabs", selected = "data_table")
})

# 默认选中第一个菜单项（Missing Value Filter）
shiny::observe({
  message("[DEBUG] preprocessing_nav: setting default tab to missing_filter")
  updateTabsetPanel(session, "main_nav_tabs", selected = "missing_filter")
})

message("[DEBUG] preprocessing_nav.R loaded successfully")