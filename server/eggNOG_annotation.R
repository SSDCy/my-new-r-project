# server/eggNOG_annotation.R
message("[DEBUG] eggNOG_annotation.R loading... (manual upload only, automatic branch removed)")

# 仅使用手动上传文件
eggnog_raw <- reactive({
  req(input$eggnog_file)
  file_path <- input$eggnog_file$datapath
  message("[DEBUG] eggnog_raw: reading manual file: ", file_path)
  tryCatch({
    first_lines <- readLines(file_path, n = 3, warn = FALSE)
    if (any(grepl("^>", first_lines))) {
      showNotification("The uploaded file appears to be a FASTA file, not an eggNOG annotation. Please upload the .emapper.annotations file.", type = "error", duration = 10)
      return(NULL)
    }
    sep <- if (grepl("\t", first_lines[1])) "\t" else ","
    df <- read.csv(file_path, sep = sep, stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)
    if (ncol(df) < 2) {
      showNotification("The uploaded file has only one column. Please upload a valid eggNOG annotation file (TSV/CSV).", type = "error", duration = 10)
      return(NULL)
    }
    message("[DEBUG] eggNOG raw dimensions: ", nrow(df), " rows, ", ncol(df), " cols")
    df
  }, error = function(e) {
    message("[ERROR] eggNOG file read failed: ", e$message)
    showNotification("Failed to read eggNOG annotation file. Ensure it is a valid TSV/CSV.", type = "error")
    NULL
  })
})

eggnog_id_col <- reactive({
  df <- eggnog_raw()
  if (is.null(df)) return(NULL)
  possible <- c("#query", "query", "protein_id", "Protein_ID")
  idx <- match(possible, colnames(df))
  idx <- which(!is.na(idx))[1]
  if (!is.na(idx)) {
    id_col <- colnames(df)[idx]
    message("[DEBUG] eggnog_id_col detected: ", id_col)
  } else {
    id_col <- colnames(df)[1]
    message("[DEBUG] eggnog_id_col defaulting to first column: ", id_col)
  }
  id_col
})

eggnog_clean <- reactive({
  df <- eggnog_raw()
  id_col <- eggnog_id_col()
  if (is.null(df) || is.null(id_col)) return(NULL)
  names(df)[names(df) == id_col] <- "Master protein IDs"
  df <- df[!duplicated(df$`Master protein IDs`), ]
  message("[DEBUG] eggnog_clean: after dedup, ", nrow(df), " rows")
  df
})

eggnog_annot_cols <- reactive({
  df <- eggnog_clean()
  if (is.null(df)) return(character(0))
  cols <- setdiff(colnames(df), "Master protein IDs")
  message("[DEBUG] eggnog_annot_cols: ", paste(cols, collapse = ", "))
  cols
})

output$eggnog_loaded <- reactive({
  !is.null(eggnog_clean())
})
outputOptions(output, "eggnog_loaded", suspendWhenHidden = FALSE)

output$eggnog_status_ui <- renderUI({
  if (is.null(eggnog_clean())) {
    return(div(style = "margin-top: 10px; color: #999;", "No eggNOG annotation loaded. Upload manually."))
  }
  div(style = "margin-top: 10px; color: #5cb85c;",
      icon("check-circle"), " eggNOG annotation loaded: ", nrow(eggnog_clean()), " entries")
})

output$eggnog_preview <- DT::renderDT({
  req(eggnog_clean())
  DT::datatable(eggnog_clean(), options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
})

add_eggnog_to_table <- function(protein_table, id_col_name = "Master protein IDs") {
  annot <- eggnog_clean()
  if (is.null(annot) || nrow(annot) == 0) {
    message("[DEBUG] add_eggnog_to_table: no annotation, returning original")
    return(protein_table)
  }
  if (!(id_col_name %in% colnames(protein_table))) {
    message("[DEBUG] add_eggnog_to_table: column '", id_col_name, "' not found in protein table, returning original")
    return(protein_table)
  }
  annot_copy <- annot
  names(annot_copy)[names(annot_copy) == "Master protein IDs"] <- id_col_name
  merged <- merge(protein_table, annot_copy, by = id_col_name, all.x = TRUE, sort = FALSE)
  message("[DEBUG] add_eggnog_to_table: merged annotation using column '", id_col_name, "', new cols: ", 
          paste(setdiff(colnames(merged), colnames(protein_table)), collapse = ", "))
  return(merged)
}