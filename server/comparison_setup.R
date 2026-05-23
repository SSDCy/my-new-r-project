# server/comparison_setup.R

# 更新选择器
observe({
  group_names <- names(rv$groups)
  if (length(group_names) > 0) {
    updateSelectInput(session, "comp_treat", choices = group_names, selected = group_names[1])
    updateSelectInput(session, "comp_ctrl", choices = group_names, selected = if(length(group_names)>1) group_names[2] else group_names[1])
    updateSelectInput(session, "batch_ref_group", choices = group_names, selected = group_names[1])
  } else {
    updateSelectInput(session, "comp_treat", choices = character(0))
    updateSelectInput(session, "comp_ctrl", choices = character(0))
    updateSelectInput(session, "batch_ref_group", choices = character(0))
  }
})

# 添加单个比较
observeEvent(input$add_comparison, {
  req(input$comp_treat, input$comp_ctrl)
  treat <- input$comp_treat; ctrl <- input$comp_ctrl; custom_name <- input$comp_name
  if (treat == ctrl) { showNotification("Treatment and control groups must be different.", type = "error", duration = 2); return() }
  if (length(rv$groups[[treat]]) == 0) { showNotification(paste("Group", treat, "has no samples."), type = "error", duration = 2); return() }
  if (length(rv$groups[[ctrl]]) == 0) { showNotification(paste("Group", ctrl, "has no samples."), type = "error", duration = 2); return() }
  comp_name <- if (is.null(custom_name) || custom_name == "") paste(treat, "vs", ctrl) else custom_name
  if (any(sapply(rv$comparisons, function(c) c$name == comp_name))) {
    showNotification("A comparison with this name already exists.", type = "error", duration = 3)
    return()
  }
  if (any(sapply(rv$comparisons, function(c) c$treat == treat && c$ctrl == ctrl))) {
    showNotification("This comparison already exists.", type = "warning", duration = 2); return()
  }
  rv$comp_id_counter <- rv$comp_id_counter + 1
  comp_id <- as.character(rv$comp_id_counter)
  rv$comparisons[[length(rv$comparisons) + 1]] <- list(id = comp_id, treat = treat, ctrl = ctrl, name = comp_name)
  updateTextInput(session, "comp_name", value = "")
  manual_sort_active(FALSE)
  showNotification("Comparison added.", type = "message", duration = 2)
})

# 批量添加配对比较
observeEvent(input$batch_add_pairwise, {
  req(input$batch_ref_group)
  ref_group <- input$batch_ref_group
  all_groups <- names(rv$groups)
  other_groups <- setdiff(all_groups, ref_group)
  if (length(other_groups) == 0) {
    showNotification("No other groups to compare against.", type = "warning", id = "comp_msg")
    return()
  }
  added <- 0
  for (other in other_groups) {
    comp_name <- paste(other, "vs", ref_group)
    if (any(sapply(rv$comparisons, function(c) c$name == comp_name))) next
    if (any(sapply(rv$comparisons, function(c) c$treat == other && c$ctrl == ref_group))) next
    if (length(rv$groups[[other]]) == 0 || length(rv$groups[[ref_group]]) == 0) next
    rv$comp_id_counter <- rv$comp_id_counter + 1
    comp_id <- as.character(rv$comp_id_counter)
    rv$comparisons[[length(rv$comparisons) + 1]] <- list(id = comp_id, treat = other, ctrl = ref_group, name = comp_name)
    added <- added + 1
  }
  manual_sort_active(FALSE)
  showNotification(paste("Added", added, "pairwise comparisons."), type = "message", duration = 3, id = "comp_msg")
})

# 排序函数
sort_comparisons <- function(comps) {
  if (length(comps) == 0) return(comps)
  comps_clean <- lapply(comps, function(c) { c$ctrl <- trimws(c$ctrl); c$treat <- trimws(c$treat); c })
  all_names <- unique(c(sapply(comps_clean, `[[`, "ctrl"), sapply(comps_clean, `[[`, "treat")))
  sorted_names <- all_names[mixedorder(all_names)]
  name_rank <- setNames(seq_along(sorted_names), sorted_names)
  ord <- order(name_rank[sapply(comps_clean, `[[`, "ctrl")], name_rank[sapply(comps_clean, `[[`, "treat")])
  comps_clean[ord]
}

sorted_comps <- reactive({
  if (manual_sort_active()) rv$comparisons else sort_comparisons(rv$comparisons)
})

# 比较列表UI
output$comparisons_list_ui <- renderUI({
  sorted <- sorted_comps()
  if (length(sorted) == 0) {
    div(style = "color: #999; text-align: center; padding: 20px;", "No comparisons defined yet.")
  } else {
    lapply(sorted, function(comp) {
      comp_id <- comp$id
      remove_btn_id <- paste0("remove_comp_", comp_id)
      div(class = "comparison-item", 
          draggable = "true",
          `data-comp-id` = comp_id,
          span(class = "drag-handle", icon("grip-vertical")),
          div(class = "comparison-label", comp$name),
          div(style = "color: #666; font-size: 12px;", paste0("[", comp$treat, " vs ", comp$ctrl, "]")),
          div(style = "flex-grow: 1;"),
          actionButton(remove_btn_id, "", icon = icon("times"), class = "btn-sm btn-outline-danger")
      )
    })
  }
})

output$comparisons_count_text <- renderUI({
  count <- length(rv$comparisons)
  HTML(paste0(" Defined Comparisons (", count, ")"))
})

# 拖放排序
observeEvent(input$comparison_drag, {
  req(input$comparison_drag)
  drag_info <- input$comparison_drag
  dragged_id <- drag_info$dragged; target_id <- drag_info$target
  if (is.null(dragged_id) || is.null(target_id) || dragged_id == target_id) return()
  comps <- isolate(rv$comparisons)
  ids <- sapply(comps, `[[`, "id")
  i <- which(ids == dragged_id); j <- which(ids == target_id)
  if (length(i) != 1 || length(j) != 1) return()
  elem <- comps[[i]]
  comps <- comps[-i]
  if (i < j) insert_pos <- j - 1 else insert_pos <- j
  if (insert_pos <= 1) comps <- c(list(elem), comps)
  else if (insert_pos > length(comps)) comps <- c(comps, list(elem))
  else comps <- c(comps[1:(insert_pos-1)], list(elem), comps[insert_pos:length(comps)])
  rv$comparisons <- comps
  manual_sort_active(TRUE)
  showNotification("Comparison order updated.", type = "message", duration = 1)
})

# 自动排序
observeEvent(input$auto_sort_comparisons, {
  manual_sort_active(FALSE)
  showNotification("Comparisons re-sorted naturally.", type = "message", duration = 2)
})

# 删除单个比较
observe({
  req(rv$comparisons)
  for (comp in rv$comparisons) {
    local({
      comp_id <- comp$id
      remove_btn <- paste0("remove_comp_", comp_id)
      observeEvent(input[[remove_btn]], {
        rv$comparisons <- Filter(function(x) x$id != comp_id, rv$comparisons)
        manual_sort_active(FALSE)
        showNotification("Comparison removed.", type = "message", duration = 2)
      }, ignoreInit = TRUE, once = TRUE)
    })
  }
})

# 清除所有比较
observeEvent(input$clear_comparisons, {
  rv$comparisons <- list()
  manual_sort_active(FALSE)
  showNotification("All comparisons cleared.", type = "message", duration = 2)
})

# 一键填充重复数
observeEvent(input$apply_replicate_fill, {
  val <- input$replicate_fill_all
  if (!is.na(val) && val >= 1 && val <= 10) {
    updateNumericInput(session, "min_treat_valid", value = val)
    updateNumericInput(session, "min_ctrl_valid", value = val)
    updateNumericInput(session, "min_rep_ttest", value = val)
    updateNumericInput(session, "min_rep_inc", value = val)
    updateNumericInput(session, "min_rep_dec", value = val)
  }
})

# 跳转到比较页面
observeEvent(input$goto_comparisons, {
  updateNavbarPage(session, "main_navbar", selected = "comparisons")
})

# ---------- 导出相关 UI ----------
output$export_comparisons_ui <- renderUI({
  comps <- sapply(sorted_comps(), `[[`, "name")
  if (length(comps) == 0) return(p("No comparisons defined yet."))
  checkboxGroupInput("export_comparisons", NULL, choices = comps, selected = character(0), inline = TRUE)
})

output$export_comps_count_text <- renderText({
  comps <- sapply(sorted_comps(), `[[`, "name")
  paste0("(", length(comps), " comparisons)")
})

observeEvent(input$select_all_export, {
  comps <- sapply(sorted_comps(), `[[`, "name")
  if (length(comps) > 10) {
    showNotification("Excel export is limited to 10 comparisons. Only the first 10 will be selected.", type = "warning", duration = 5)
    comps <- comps[1:10]
  }
  updateCheckboxGroupInput(session, "export_comparisons", selected = comps)
})

observeEvent(input$deselect_all_export, {
  updateCheckboxGroupInput(session, "export_comparisons", selected = character(0))
})

# ---------- UpSet 复选框 UI ----------
output$upset_comparisons_checkbox_ui <- renderUI({
  comps <- sapply(sorted_comps(), `[[`, "name")
  if (length(comps) == 0) return(p("No comparisons defined."))
  tagList(
    div(style = "display: flex; align-items: center; gap: 10px; margin-bottom: 5px;",
        actionButton("select_toggle_upset", "Select All / Deselect All", class = "btn-sm btn-outline-secondary")
    ),
    checkboxGroupInput("venn_comparisons_checkbox", NULL, choices = comps, selected = NULL, inline = TRUE)
  )
})

observeEvent(input$select_toggle_upset, {
  comps <- sapply(sorted_comps(), `[[`, "name")
  current_selected <- input$venn_comparisons_checkbox
  if (length(current_selected) == length(comps)) {
    updateCheckboxGroupInput(session, "venn_comparisons_checkbox", selected = character(0))
  } else {
    if (length(comps) > 15) {
      showNotification("UpSet is limited to 15 comparisons. Selecting the first 15.", type = "warning")
      updateCheckboxGroupInput(session, "venn_comparisons_checkbox", selected = comps[1:15])
    } else {
      updateCheckboxGroupInput(session, "venn_comparisons_checkbox", selected = comps)
    }
  }
})