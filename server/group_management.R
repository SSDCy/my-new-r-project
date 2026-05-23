# server/group_management.R

# 未分配样本UI
output$unassigned_samples_ui <- renderUI({
  req(rv$sample_names)
  assigned_samples <- unlist(rv$groups)
  unassigned <- setdiff(rv$sample_names, assigned_samples)
  if (length(unassigned) == 0) {
    div(style = "color: #999; text-align: center; padding: 10px;", "All samples assigned!")
  } else {
    div(
      lapply(unassigned, function(s) {
        div(class = "sample-item-draggable", draggable = "true", `data-sample` = s,
            icon("vial"), " ", s)
      })
    )
  }
})

# 分组层级选择更新
observe({
  req(rv$sample_names)
  assigned <- unlist(rv$groups)
  unassigned <- setdiff(rv$sample_names, assigned)
  if (length(unassigned) == 0) {
    updateSelectInput(session, "group_level", choices = list("No unassigned samples" = ""))
    return()
  }
  levels <- parse_sample_levels(unassigned)
  if (length(levels) == 0) {
    updateSelectInput(session, "group_level", choices = list("Default (prefix)" = "default"))
    return()
  }
  display_texts <- sapply(names(levels), function(level_num) {
    l <- levels[[level_num]]
    paste0("Level ", l$level, " (e.g., ", l$example, ")")
  })
  return_values <- names(levels)
  choices <- setNames(return_values, display_texts)
  choices <- c(choices, "Default (prefix)" = "default")
  updateSelectInput(session, "group_level", choices = choices, selected = "default")
})

# 添加分组
observeEvent(input$add_group, {
  req(input$new_group_name)
  new_name <- trimws(input$new_group_name)
  if (new_name == "") { showNotification("Please enter a group name.", type = "error", duration = 2); return() }
  if (new_name %in% names(rv$groups)) { showNotification("Group name already exists.", type = "error", duration = 2); return() }
  rv$group_id_counter <- rv$group_id_counter + 1
  gid <- as.character(rv$group_id_counter)
  rv$group_id_map[[gid]] <- new_name
  rv$groups[[new_name]] <- character(0)
  updateTextInput(session, "new_group_name", value = "")
  showNotification(paste("Group", new_name, "created."), type = "message", duration = 2, id = "group_msg")
})

# 重置分组
observeEvent(input$reset_groups, {
  rv$groups <- list()
  rv$group_id_counter <- 0
  rv$group_id_map <- list()
  rv$comparisons <- list()
  manual_sort_active(FALSE)
  
  updateSelectInput(session, "comp_treat", choices = character(0))
  updateSelectInput(session, "comp_ctrl", choices = character(0))
  updateSelectInput(session, "selected_comparison", choices = character(0))
  updateSelectInput(session, "batch_ref_group", choices = character(0))
  
  showNotification("All groups and comparisons have been reset.", type = "message", duration = 3, id = "group_msg")
})

# 删除分组
observeEvent(input$remove_group, {
  gid <- input$remove_group
  gname <- rv$group_id_map[[gid]]
  if (is.null(gname)) return()
  rv$groups[[gname]] <- NULL
  rv$group_id_map[[gid]] <- NULL
  rv$comparisons <- Filter(function(comp) !(comp$treat == gname | comp$ctrl == gname), rv$comparisons)
  manual_sort_active(FALSE)
  showNotification(paste("Group", gname, "removed."), type = "message", duration = 2, id = "group_msg")
})

# 清空分组
observeEvent(input$empty_group, {
  gid <- input$empty_group
  gname <- rv$group_id_map[[gid]]
  if (is.null(gname)) return()
  rv$groups[[gname]] <- character(0)
  showNotification(paste("Group", gname, "emptied."), type = "message", duration = 2, id = "group_msg")
})

# 重命名分组请求
observeEvent(input$rename_group_request, {
  gid <- input$rename_group_request
  gname <- rv$group_id_map[[gid]]
  if (is.null(gname)) return()
  showModal(modalDialog(
    title = paste("Rename Group:", gname),
    textInput("rename_group_new_name", "New group name", value = gname),
    footer = tagList(
      actionButton("confirm_rename_group", "Rename", class = "btn btn-secondary"),
      tags$button("Cancel", type = "button", class = "btn btn-primary", `data-dismiss` = "modal", `data-bs-dismiss` = "modal")
    ),
    easyClose = TRUE
  ))
  rv$pending_rename_gid <- gid
})

# 确认重命名
observeEvent(input$confirm_rename_group, {
  new_name <- trimws(input$rename_group_new_name)
  gid <- rv$pending_rename_gid
  if (is.null(gid) || new_name == "") { removeModal(); return() }
  old_name <- rv$group_id_map[[gid]]
  if (is.null(old_name) || new_name == old_name) { removeModal(); return() }
  if (new_name %in% names(rv$groups)) {
    showNotification("A group with this name already exists.", type = "error", duration = 2); return()
  }
  names(rv$groups)[names(rv$groups) == old_name] <- new_name
  rv$group_id_map[[gid]] <- new_name
  for (i in seq_along(rv$comparisons)) {
    if (rv$comparisons[[i]]$treat == old_name) rv$comparisons[[i]]$treat <- new_name
    if (rv$comparisons[[i]]$ctrl == old_name) rv$comparisons[[i]]$ctrl <- new_name
  }
  removeModal()
  rv$pending_rename_gid <- NULL
  showNotification(paste("Group renamed from", old_name, "to", new_name), type = "message", duration = 2, id = "group_msg")
})

# 从分组中移除选中样本
observeEvent(input$unassign_sel_group, {
  gid <- input$unassign_sel_group
  gname <- rv$group_id_map[[gid]]
  if (is.null(gname)) return()
  selected <- input[[paste0("group_samples_", gid)]]
  if (!is.null(selected) && length(selected) > 0) {
    rv$groups[[gname]] <- setdiff(rv$groups[[gname]], selected)
    showNotification(paste("Removed", length(selected), "samples from", gname), type = "message", duration = 2, id = "group_msg")
  } else {
    showNotification("No samples selected to remove.", type = "warning", duration = 2, id = "group_msg")
  }
})

# 分组UI
output$groups_ui <- renderUI({
  req(rv$groups)
  if (length(rv$groups) == 0) {
    div(style = "color: #999; text-align: center; padding: 20px;", "No groups defined.")
  } else {
    group_names <- names(rv$groups)
    lapply(group_names, function(gname) {
      gid <- names(rv$group_id_map)[rv$group_id_map == gname]
      if (length(gid) == 0) return(NULL)
      samples <- rv$groups[[gname]]
      div(class = "group-box",
          div(class = "group-header",
              h5(style = "margin: 0;", gname),
              span(style = "color: #666; font-size: 12px;", paste0("(", length(samples), " samples)")),
              div(style = "flex-grow: 1;"),
              tags$button("Rename", class = "btn btn-sm btn-outline-info",
                          onclick = sprintf("Shiny.setInputValue('rename_group_request', '%s', {priority: 'event'})", gid)),
              tags$button("Remove Group", class = "btn btn-sm btn-danger",
                          onclick = sprintf("Shiny.setInputValue('remove_group', '%s', {priority: 'event'})", gid)),
              tags$button("Empty Group", class = "btn btn-sm btn-outline-secondary",
                          onclick = sprintf("Shiny.setInputValue('empty_group', '%s', {priority: 'event'})", gid))
          ),
          div(class = "group-drop-zone", `data-group` = gname,
              icon("arrow-down"), " Drop samples here"
          ),
          if (length(samples) > 0) {
            div(style = "max-height: 150px; overflow-y: auto;",
                p(class = "param-hint", "Select samples to remove from this group:"),
                checkboxGroupInput(paste0("group_samples_", gid), NULL, choices = samples, selected = NULL, inline = TRUE),
                tags$button("Remove Selected", class = "btn btn-sm btn-outline-warning",
                            onclick = sprintf("Shiny.setInputValue('unassign_sel_group', '%s', {priority: 'event'})", gid))
            )
          } else {
            div(style = "color: #999; font-style: italic;", "No samples in this group yet.")
          }
      )
    })
  }
})

# 拖放分配
observeEvent(input$drag_assign, {
  req(input$drag_assign)
  data <- jsonlite::fromJSON(input$drag_assign)
  group <- data$group
  samples <- data$samples
  if (!is.null(samples) && length(samples) > 0 && !is.null(group) && group %in% names(rv$groups)) {
    assigned_samples <- unlist(rv$groups)
    valid_samples <- samples[samples %in% rv$sample_names & !samples %in% assigned_samples]
    if (length(valid_samples) > 0) {
      rv$groups[[group]] <- unique(c(rv$groups[[group]], valid_samples))
      showNotification(paste("Assigned", length(valid_samples), "samples to", group), type = "message", duration = 2, id = "group_msg")
    }
  }
})

# 自动分配
observeEvent(input$auto_assign, {
  req(rv$sample_names)
  groups <- rv$groups
  if (length(groups) == 0) {
    showNotification("No groups defined. Please create groups first.", type = "warning")
    return()
  }
  assigned <- unlist(groups)
  unassigned <- setdiff(rv$sample_names, assigned)
  if (length(unassigned) == 0) {
    showNotification("All samples are already assigned.", type = "message")
    return()
  }
  group_names <- names(groups)
  assigned_count <- 0
  for (samp in unassigned) {
    matched_group <- NULL
    matches <- group_names[sapply(group_names, function(gn) grepl(tolower(gn), tolower(samp), fixed = TRUE))]
    if (length(matches) > 0) {
      matched_group <- matches[which.max(nchar(matches))]
    } else {
      dists <- adist(samp, group_names, ignore.case = TRUE)
      best_idx <- which.min(dists)
      if (dists[best_idx] <= 2) matched_group <- group_names[best_idx]
    }
    if (!is.null(matched_group)) {
      rv$groups[[matched_group]] <- unique(c(rv$groups[[matched_group]], samp))
      assigned_count <- assigned_count + 1
    }
  }
  showNotification(paste("Auto-assigned", assigned_count, "out of", length(unassigned), "unassigned samples."), type = "message", duration = 5, id = "group_msg")
})

# 批量创建分组
observeEvent(input$batch_create_groups, {
  req(rv$sample_names, input$group_level)
  assigned <- unlist(rv$groups)
  unassigned <- setdiff(rv$sample_names, assigned)
  if (length(unassigned) == 0) {
    showNotification("All samples are already assigned.", type = "message")
    return()
  }
  levels <- parse_sample_levels(unassigned)
  selected_level <- input$group_level
  separator <- if (selected_level != "default" && !is.null(levels[[selected_level]])) {
    levels[[selected_level]]$separator
  } else {
    "-"
  }
  prefix_map <- sapply(unassigned, function(s) extract_group_prefix(s, selected_level, separator))
  unique_prefixes <- unique(prefix_map)
  valid_prefixes <- unique_prefixes[nchar(unique_prefixes) >= 2 & !is.na(unique_prefixes)]
  if (length(valid_prefixes) == 0) {
    showNotification("Cannot extract valid group prefixes.", type = "warning")
    return()
  }
  new_groups_created <- 0
  for (pref in valid_prefixes) {
    if (!pref %in% names(rv$groups)) {
      base_name <- pref
      if (base_name %in% names(rv$groups)) base_name <- paste0(pref, "_group")
      rv$group_id_counter <- rv$group_id_counter + 1
      gid <- as.character(rv$group_id_counter)
      rv$group_id_map[[gid]] <- base_name
      rv$groups[[base_name]] <- character(0)
      new_groups_created <- new_groups_created + 1
    }
  }
  for (samp in unassigned) {
    pref <- extract_group_prefix(samp, selected_level, separator)
    if (!is.na(pref) && pref %in% names(rv$groups)) {
      rv$groups[[pref]] <- unique(c(rv$groups[[pref]], samp))
    }
  }
  showNotification(paste("Batch groups created:", new_groups_created, "new groups with auto-assigned samples."), type = "message", duration = 5, id = "group_msg")
})

observeEvent(input$confirm_groups, {
  updateNavbarPage(session, "main_navbar", selected = "comparisons")
})