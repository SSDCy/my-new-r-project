# server/reactive_values.R
message("[DEBUG] reactive_values.R: loading...")

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
  # auto_eggnog_result removed
)

message("[DEBUG] reactive_values.R: auto_eggnog_result field removed.")

subplot_old_values <- reactiveValues()

manual_sort_active <- reactiveVal(FALSE)
message("[DEBUG] reactive_values.R: manual_sort_active initialized to FALSE")

clicked_protein <- reactiveVal(NULL)

data_changed_trigger <- reactiveVal(0)
message("[DEBUG] reactive_values.R: data_changed_trigger initialized to 0")

heatmap_generated_version <- reactiveVal(-1)
message("[DEBUG] reactive_values.R: heatmap_generated_version initialized to -1")

message("[DEBUG] reactive_values.R: all reactive values loaded successfully.")