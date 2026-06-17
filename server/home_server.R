# server/home_server.R
message("[DEBUG] home_server.R loaded")

# 处理首页“开始”按钮，跳转到 Data Upload 页面
observeEvent(input$home_start_btn, {
  message("[DEBUG] home_start_btn clicked - navigating to upload tab")
  updateNavbarPage(session, "main_navbar", selected = "upload")
})

message("[DEBUG] home_server.R fully loaded")