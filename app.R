# app.R

library(shiny)

source("global.R", local = TRUE)
source("ui/upload_ui.R", local = TRUE)
source("ui/cleaning_ui.R", local = TRUE)
source("ui/normalization_ui.R", local = TRUE)   # 归一化 UI
source("ui/preprocessing_ui.R", local = TRUE)
source("ui/grouping_ui.R", local = TRUE)         # 分组页面
source("ui/comparisons_ui.R", local = TRUE)      # 比较页面
source("ui/parameters_ui.R", local = TRUE)       # 参数页面
source("ui/plots_ui.R", local = TRUE)
source("ui/export_ui.R", local = TRUE)

build_ui <- function() {
  fluidPage(
    title = "Universal Proteomics Differential Analysis Platform",
    theme = bslib::bs_theme(
      version = 5, bootswatch = "flatly",
      primary = "#2c3e50", secondary = "#3498db",
      success = "#18bc9c", info = "#3498db",
      warning = "#f39c12", danger = "#e74c3c"
    ),
    shinyjs::useShinyjs(),
    tags$head(
      tags$style(HTML("
        .navbar-brand {font-weight: bold; font-size: 18px;}
        .card-modern {border-radius: 12px; box-shadow: 0 2px 10px rgba(0,0,0,0.08); margin-bottom: 20px; border: none;}
        .card-header-modern {background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; border-radius: 12px 12px 0 0; padding: 15px 20px; font-weight: bold;}
        .status-badge {display: inline-block; padding: 8px 16px; border-radius: 20px; font-size: 14px; font-weight: bold;}
        .status-success {background: #d4edda; color: #155724; border: 1px solid #c3e6cb;}
        .status-warning {background: #fff3cd; color: #856404; border: 1px solid #ffeeba;}
        .sample-item {padding: 6px 10px; margin: 2px 0; background: white; border: 1px solid #ddd; border-radius: 4px; font-family: monospace; font-size: 11px; display: inline-block;}
        .sample-item-draggable { padding: 8px 12px; margin: 4px 0; background: white; border: 1px solid #bbb; border-radius: 6px; font-family: monospace; font-size: 12px; cursor: pointer; display: inline-block; transition: background 0.2s, box-shadow 0.2s; user-select: none; }
        .sample-item-draggable:hover { background: #e3f2fd; box-shadow: 0 2px 6px rgba(0,0,0,0.1); }
        .sample-item-draggable.selected { background: #cce5ff; border-color: #3498db; box-shadow: 0 0 0 2px rgba(52,152,219,0.3); }
        .group-box {border: 2px solid #667eea; border-radius: 10px; padding: 15px; margin-bottom: 15px; background: #f8f9fa;}
        .group-header {display: flex; align-items: center; gap: 10px; margin-bottom: 10px;}
        .group-drop-zone { min-height: 40px; border: 2px dashed #aaa; border-radius: 8px; padding: 10px; text-align: center; color: #888; font-size: 12px; margin-bottom: 10px; transition: background 0.2s, border-color 0.2s; }
        .group-drop-zone.drag-over { background: #eef; border-color: #667eea; }
        .comparison-item {display: flex; align-items: center; gap: 10px; padding: 10px; background: #f8f9fa; border-radius: 8px; margin-bottom: 8px; border: 2px solid transparent; transition: background 0.2s, border-color 0.2s;}
        .comparison-label {font-weight: bold; color: #2c3e50;}
        .comparison-item[draggable=true] { cursor: grab; }
        .comparison-item[draggable=true]:active { cursor: grabbing; }
        .comparison-item.dragging { opacity: 0.4; background: #fce4ec; border-color: #e57373; }
        .comparison-item.drag-over { border: 2px dashed #667eea; background: #eef; }
        .comparison-item .drag-handle { color: #aaa; font-size: 16px; margin-right: 4px; cursor: grab; }
        .color-palette-row { display: flex; justify-content: center; gap: 20px; flex-wrap: wrap; margin: 20px 0; align-items: flex-start; }
        .color-card { background: white; border-radius: 16px; box-shadow: 0 4px 14px rgba(0,0,0,0.06); padding: 15px 12px; width: 130px; text-align: center; transition: transform 0.2s ease, box-shadow 0.2s ease; }
        .color-card:hover { transform: translateY(-4px); box-shadow: 0 10px 24px rgba(0,0,0,0.12); }
        .color-card-label { font-size: 14px; font-weight: 700; color: #2c3e50; margin-bottom: 12px; letter-spacing: 0.5px; }
        .colourpicker-container { width: 60px !important; height: 120px !important; margin: 0 auto !important; cursor: pointer !important; position: relative !important; display: block !important; }
        .colourpicker-container .form-control { width: 60px !important; height: 120px !important; border-radius: 12px !important; border: 3px solid #e0e0e0 !important; cursor: pointer !important; opacity: 1 !important; display: block !important; visibility: visible !important; padding: 0 !important; margin: 0 !important; text-indent: -9999px !important; font-size: 0 !important; line-height: 0 !important; background-clip: padding-box !important; transition: all 0.3s ease !important; }
        .colourpicker-container .form-control:hover { border-color: #3498db !important; box-shadow: 0 0 0 5px rgba(52,152,219,0.15) !important; }
        .color-card-value { margin-top: 12px; font-size: 12px; font-weight: 600; color: #555; background: #f5f7fa; padding: 3px 10px; border-radius: 20px; display: inline-block; font-family: 'Courier New', monospace; }
        .reset-btn-wrapper { display: flex; align-items: center; justify-content: center; padding-top: 40px; }
        .btn-circle-modern { background-color: #3498db; color: white; font-size: 20px; width: 50px; height: 50px; border-radius: 50%; border: none; box-shadow: 0 4px 12px rgba(52,152,219,0.3); transition: all 0.2s; display: flex; align-items: center; justify-content: center; cursor: pointer; }
        .btn-circle-modern:hover { background-color: #2c3e50; box-shadow: 0 6px 16px rgba(44,62,80,0.4); }
        .color-preview-container { display: flex; justify-content: center; align-items: center; gap: 25px; margin: 10px 0 15px 0; padding: 12px; background: #f8f9fa; border-radius: 12px; border: 1px solid #dee2e6; }
        .color-preview-item {text-align: center;}
        .color-preview-swatch { width: 50px; height: 20px; border-radius: 6px; border: 1px solid #ccc; margin: 0 auto 4px auto; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
        .color-preview-label {font-size: 11px; font-weight: 700; color: #333;}
        .title-setting-panel {background: #f8f9fa; padding: 15px; border-radius: 10px; margin-bottom: 15px; border-left: 4px solid #667eea;}
        .param-hint {font-size: 12px; color: #7f8c8d; margin-top: 2px; margin-bottom: 10px;}
        .input-warning {color: #d9534f; font-size: 12px; margin-top: 4px; font-weight: bold;}
        .sample-pool {background: #fff3cd; border: 2px dashed #f0ad4e; border-radius: 10px; padding: 15px; min-height: 100px;}
        .sticky-panel { position: sticky; top: 80px; z-index: 10; }
        .floating-confirm-btn { position: fixed; bottom: 20px; right: 30px; z-index: 999; box-shadow: 0 6px 20px rgba(0,0,0,0.25); border-radius: 50px; padding: 12px 28px; font-size: 18px; font-weight: bold; }
        .red-text { color: #d9534f; font-size: 13px; font-weight: bold; }
        .input-row-with-reset { display: flex; align-items: flex-end; gap: 10px; }
        .input-row-with-reset .form-group { margin-bottom: 0; flex-grow: 1; }
        .btn-reset-small { width: 38px; height: 38px; border-radius: 6px; display: flex; align-items: center; justify-content: center; margin-bottom: 1rem; }
        .batch-group-row { display: flex; gap: 10px; align-items: flex-end; margin-bottom: 15px; flex-wrap: wrap; }
        details { margin-bottom: 15px; }
        details summary { font-weight: bold; cursor: pointer; color: #2c3e50; padding: 5px 0; outline: none; }
        .checkbox-inline-group { display: flex; flex-wrap: wrap; gap: 10px; }
        .scrollable-box { max-height: 300px; overflow-y: auto; border: 1px solid #dee2e6; border-radius: 8px; padding: 10px; background: #fff; }
        .export-details { margin-top: 30px; border: 1px solid #ddd; border-radius: 8px; padding: 15px; background: #fafafa; }
      ")),
      tags$script(HTML("
        var lastSelectedIndex = -1;
        $(document).on('click', '.sample-item-draggable', function(e) {
          var items = $('.sample-item-draggable');
          var idx = items.index(this);
          if (e.ctrlKey || e.metaKey) { $(this).toggleClass('selected'); lastSelectedIndex = idx; }
          else if (e.shiftKey && lastSelectedIndex >= 0) {
            var start = Math.min(lastSelectedIndex, idx);
            var end = Math.max(lastSelectedIndex, idx);
            items.slice(start, end + 1).addClass('selected');
          } else { items.removeClass('selected'); $(this).addClass('selected'); lastSelectedIndex = idx; }
        });
        $(document).on('dragstart', '.sample-item-draggable', function(e) {
          var selected = [];
          $('.sample-item-draggable.selected').each(function() { selected.push($(this).attr('data-sample')); });
          if (selected.length === 0) { selected.push($(this).attr('data-sample')); }
          e.originalEvent.dataTransfer.setData('text/plain', JSON.stringify(selected));
        });
        $(document).on('dragover', '.group-drop-zone', function(e) { e.preventDefault(); $(this).addClass('drag-over'); });
        $(document).on('dragleave', '.group-drop-zone', function(e) { $(this).removeClass('drag-over'); });
        $(document).on('drop', '.group-drop-zone', function(e) {
          e.preventDefault(); $(this).removeClass('drag-over');
          var data = e.originalEvent.dataTransfer.getData('text/plain');
          var group = $(this).attr('data-group');
          if (data && group) { Shiny.setInputValue('drag_assign', JSON.stringify({samples: JSON.parse(data), group: group}), {priority: 'event'}); $('.sample-item-draggable.selected').removeClass('selected'); }
        });
        $(document).on('dragstart', '.comparison-item[draggable=true]', function(e) {
          var compId = $(this).data('comp-id');
          e.originalEvent.dataTransfer.setData('text/plain', compId);
          $(this).addClass('dragging');
          e.originalEvent.dataTransfer.effectAllowed = 'move';
        });
        $(document).on('dragend', '.comparison-item', function(e) { $(this).removeClass('dragging'); $('.comparison-item').removeClass('drag-over'); });
        $(document).on('dragover', '.comparison-item', function(e) {
          e.preventDefault();
          if (!$(this).hasClass('dragging')) { $(this).addClass('drag-over'); }
          e.originalEvent.dataTransfer.dropEffect = 'move';
        });
        $(document).on('dragleave', '.comparison-item', function(e) { $(this).removeClass('drag-over'); });
        $(document).on('drop', '.comparison-item', function(e) {
          e.preventDefault(); e.stopPropagation(); $(this).removeClass('drag-over');
          var draggedId = e.originalEvent.dataTransfer.getData('text/plain');
          var targetId = $(this).data('comp-id');
          if (draggedId && targetId && draggedId !== targetId) {
            Shiny.setInputValue('comparison_drag', {dragged: draggedId, target: targetId}, {priority: 'event'});
          }
        });
        $(document).on('input', '#plot_width, #plot_height', function(e) {
          this.value = this.value.replace(/[^0-9]/g, '');
          var val = parseInt(this.value);
          var id = this.id;
          var warningId = id + '_warning';
          if (isNaN(val) || val === 0) { $('#' + warningId).text('Please enter a valid positive integer'); $(this).addClass('is-invalid'); }
          else if (val > 30) { $('#' + warningId).text('Value too large, max recommended is 30 inches'); $(this).addClass('is-invalid'); }
          else { $('#' + warningId).text(''); $(this).removeClass('is-invalid'); }
        });
      "))
    ),
    navbarPage(
      title = div(icon("dna", style = "margin-right: 8px;"), "Universal Proteomics Platform"),
      id = "main_navbar", collapsible = TRUE, windowTitle = "Proteomics Analysis",
      upload_ui(),
      cleaning_ui(),
      preprocessing_ui(),
      normalization_ui(),
      grouping_ui(),          # 独立的分组页面
      comparisons_ui(),       # 独立的比较页面
      parameters_ui(),        # 独立的参数页面
      plots_ui(),
      export_ui()
    )
  )
}

ui <- build_ui()

server <- function(input, output, session) {
  source("server/reactive_values.R", local = TRUE)
  source("server/data_upload.R", local = TRUE)
  source("server/cleaning_server.R", local = TRUE)
  source("server/normalization_server.R", local = TRUE)
  source("server/preprocessing_helpers.R", local = TRUE)
  source("server/preprocessing_filter_missing.R", local = TRUE)
  source("server/preprocessing_filter_intensity.R", local = TRUE)
  source("server/preprocessing_core.R", local = TRUE)
  source("server/preprocessing_imputation.R", local = TRUE)
  source("server/preprocessing_batch.R", local = TRUE)
  source("server/preprocessing_comparisons.R", local = TRUE)
  source("server/group_management.R", local = TRUE)
  source("server/comparison_setup.R", local = TRUE)
  source("server/de_analysis.R", local = TRUE)
  source("server/heatmap_plot.R", local = TRUE)
  source("server/volcano_plot.R", local = TRUE)
  source("server/export_server.R", local = TRUE)
  source("server/input_validation.R", local = TRUE)
  source("server/data_quality_plots.R", local = TRUE)
  source("server/preprocessing_nav.R", local = TRUE)
}

shinyApp(ui = ui, server = server)