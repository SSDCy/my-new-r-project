# server/reactive_values.R
rv <- reactiveValues(
  raw_data = NULL,
  clean_data = NULL,
  lfq_cols = NULL,
  sample_names = NULL,
  sample_info = NULL,
  groups = list(),
  comparisons = list(),
  analysis_results = NULL,
  pending_duplicate = NULL,
  reset_counter = 0,
  comp_id_counter = 0,
  group_id_counter = 0,
  group_id_map = list(),
  current_profile_protein = NULL,
  batch_vector = NULL
)

subplot_old_values <- reactiveValues()
manual_sort_active <- reactiveVal(FALSE)
clicked_protein <- reactiveVal(NULL)

# 热图数据变化触发器
data_changed_trigger <- reactiveVal(0)
# 最近一次成功生成热图时的版本号
heatmap_generated_version <- reactiveVal(-1)