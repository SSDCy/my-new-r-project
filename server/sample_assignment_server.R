# server/sample_assignment_server.R
message("[DEBUG] sample_assignment_server.R loaded (custom SubGroup with prompt, suspendWhenHidden = FALSE)")

assignment_data <- reactiveVal(data.frame())

extract_subgroup_suggestions <- function(sample_names) {
  suggestions <- sapply(sample_names, function(s) {
    m <- regmatches(s, gregexpr("\\d+", s))[[1]]
    if (length(m) >= 2) {
      paste0(m[1], "-", m[2])
    } else {
      s
    }
  })
  unique(suggestions)
}

observeEvent(rv$sample_names, {
  message("[DEBUG] sample_assignment_server: sample_names changed")
  samples <- rv$sample_names
  current_info <- rv$sample_info
  
  df <- data.frame(
    Sample = samples,
    Group = rep("Treatment", length(samples)),
    SubGroup = sapply(samples, function(s) {
      m <- regmatches(s, gregexpr("\\d+", s))[[1]]
      if (length(m) >= 2) paste0(m[1], "-", m[2]) else s
    }),
    stringsAsFactors = FALSE
  )
  assignment_data(df)
  message("[DEBUG] assignment table initialized with ", nrow(df), " rows")
})

subgroup_options <- reactive({
  req(rv$sample_names)
  extract_subgroup_suggestions(rv$sample_names)
})

output$assignment_table <- DT::renderDataTable({
  df <- assignment_data()
  req(nrow(df) > 0)
  
  sub_opts <- subgroup_options()
  sub_js_array <- paste0("[", paste(shQuote(sub_opts), collapse = ","), "]")
  
  DT::datatable(
    df,
    editable = list(target = "cell", disable = list(columns = 0)),
    options = list(
      pageLength = 50,
      autoWidth = TRUE,
      columnDefs = list(
        list(
          targets = 1,
          render = DT::JS(
            "function(data, type, row, meta) {
              var select = '<select class=\"form-select form-select-sm\" onchange=\"updateGroup(this, '\" + row + \"', '\" + meta.col + \"')\">' +
                           '<option value=\"Control\"' + (data == 'Control' ? ' selected' : '') + '>Control</option>' +
                           '<option value=\"Treatment\"' + (data == 'Treatment' ? ' selected' : '') + '>Treatment</option>' +
                           '</select>';
              return select;
            }"
          )
        ),
        list(
          targets = 2,
          render = DT::JS(
            sprintf("function(data, type, row, meta) {
              var subgroups = %s;
              var current = data;
              var inList = subgroups.includes(current);
              var select = '<select class=\"form-select form-select-sm\" onchange=\"handleSubgroupChange(this, '\" + row + \"', '\" + meta.col + \"')\">';
              subgroups.forEach(function(sub) {
                select += '<option value=\"' + sub + '\"' + (current == sub ? ' selected' : '') + '>' + sub + '</option>';
              });
              if (!inList && current !== '') {
                select += '<option value=\"__custom__\" selected>自定义 (' + current + ')</option>';
              } else {
                select += '<option value=\"__custom__\">自定义</option>';
              }
              select += '</select>';
              return select;
            }", sub_js_array)
          )
        )
      ),
      initComplete = DT::JS(
        "function(settings, json) {",
        "  $(window).trigger('resize');",  # 触发窗口调整，使表格正确计算宽度
        "}"
      )
    ),
    selection = "none",
    rownames = FALSE,
    escape = FALSE,
    callback = DT::JS("
      window.updateGroup = function(selectElement, row, col) {
        var newValue = selectElement.value;
        var table = $('#assignment_table').DataTable();
        table.cell(row, col).data(newValue).draw(false);
        Shiny.setInputValue('assignment_table_cell_edit', {row: row + 1, col: col, value: newValue}, {priority: 'event'});
      };
      
      window.handleSubgroupChange = function(selectElement, row, col) {
        var value = selectElement.value;
        if (value === '__custom__') {
          var custom = prompt('请输入自定义 SubGroup：', '');
          if (custom !== null && custom.trim() !== '') {
            custom = custom.trim();
            var option = document.createElement('option');
            option.value = custom;
            option.text = custom;
            option.selected = true;
            var customOption = selectElement.querySelector('option[value=\"__custom__\"]');
            if (customOption) customOption.remove();
            selectElement.add(option);
            var table = $('#assignment_table').DataTable();
            table.cell(row, col).data(custom).draw(false);
            Shiny.setInputValue('assignment_table_cell_edit', {row: row + 1, col: col, value: custom}, {priority: 'event'});
          } else {
            var table = $('#assignment_table').DataTable();
            var currentData = table.cell(row, col).data();
            selectElement.value = currentData;
          }
        } else {
          var table = $('#assignment_table').DataTable();
          table.cell(row, col).data(value).draw(false);
          Shiny.setInputValue('assignment_table_cell_edit', {row: row + 1, col: col, value: value}, {priority: 'event'});
        }
      };
    ")
  )
})

# 确保表格在隐藏状态下也进行渲染
outputOptions(output, "assignment_table", suspendWhenHidden = FALSE)

proxy <- DT::dataTableProxy("assignment_table")
observeEvent(input$assignment_table_cell_edit, {
  info <- input$assignment_table_cell_edit
  df <- assignment_data()
  row <- info$row
  col <- info$col + 1
  value <- info$value
  if (col == 2 || col == 3) {
    df[row, col] <- value
    assignment_data(df)
    message(sprintf("[DEBUG] cell edited: row %d, col %d, value %s", row, col, value))
  }
})

observeEvent(input$save_assignment, {
  df <- assignment_data()
  if (nrow(df) == 0) return()
  message("[DEBUG] save_assignment: saving assignment")
  
  full_cols <- rv$lfq_cols
  short_names <- extract_sample_names(full_cols)
  short_to_full <- setNames(full_cols, short_names)
  
  full_sample_names <- short_to_full[df$Sample]
  
  new_info <- data.frame(
    Group = df$Group,
    SubGroup = df$SubGroup,
    row.names = full_sample_names
  )
  rv$sample_info <- new_info
  message("[DEBUG] sample_info updated with full column names")
  
  subgroups <- unique(df$SubGroup)
  new_groups <- list()
  for (sg in subgroups) {
    samples_in_group <- df$Sample[df$SubGroup == sg]
    new_groups[[sg]] <- samples_in_group
  }
  rv$groups <- new_groups
  message("[DEBUG] groups created: ", paste(names(new_groups), collapse = ", "))
  
  rv$comparisons <- list()
  manual_sort_active(FALSE)
  
  showNotification("Sample assignment saved! Groups and SubGroups updated.", type = "message")
})

message("[DEBUG] sample_assignment_server.R loaded successfully")