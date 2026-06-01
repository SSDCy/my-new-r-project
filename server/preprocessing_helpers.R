# server/preprocessing_helpers.R
message("[DEBUG] preprocessing_helpers.R loaded")

combat_correction <- function(data, batch) {
  if (!requireNamespace("sva", quietly = TRUE))
    stop("sva package required. Run BiocManager::install('sva')")
  batch_factor <- as.factor(batch)
  if (length(levels(batch_factor)) < 2) {
    stop("Batch correction requires at least two distinct batch values.")
  }
  data_matrix <- as.matrix(data)
  message("[DEBUG] combat_correction: input dim = ", nrow(data_matrix), " x ", ncol(data_matrix))
  message("[DEBUG] combat_correction: batches = ", paste(levels(batch_factor), collapse = ", "))
  message("[DEBUG] combat_correction: table(batch) = ", paste(capture.output(print(table(batch_factor))), collapse = "\n"))
  
  if (any(is.na(data_matrix))) {
    message("[DEBUG] ComBat: Detected missing values, applying KNN imputation temporarily for batch correction.")
    data_matrix <- impute_missing_values(data, method = "knn")
    data_matrix <- as.matrix(data_matrix)
    message("[DEBUG] ComBat: missing values after imputation = ", sum(is.na(data_matrix)))
  }
  row_vars <- apply(data_matrix, 1, var)
  zero_var_rows <- sum(row_vars == 0, na.rm = TRUE)
  if (zero_var_rows > 0) {
    data_matrix <- data_matrix[row_vars > 0, , drop = FALSE]
    message("[DEBUG] ComBat: removed ", zero_var_rows, " zero-variance rows")
  }
  
  message("[DEBUG] ComBat: running ComBat...")
  corrected <- sva::ComBat(dat = data_matrix, batch = batch_factor)
  nan_count <- sum(is.na(corrected))
  inf_count <- sum(!is.finite(as.matrix(corrected)))
  if (nan_count + inf_count > 0) {
    message("[DEBUG] ComBat: produced ", nan_count, " NAs and ", inf_count, " Infs; replacing with 1e-4")
    corrected[is.na(corrected)] <- 1e-4
    corrected[!is.finite(as.matrix(corrected))] <- 1e-4
  }
  message("[DEBUG] combat_correction: completed successfully")
  result <- as.data.frame(corrected)
  rownames(result) <- rownames(data_matrix)
  return(result)
}

apply_missing_filter <- function(data, threshold, mode, sample_info = NULL, sample_names_short = NULL) {
  if (threshold >= 1) return(data)
  fallback_triggered <- FALSE
  unmatched_count <- 0
  message("[DEBUG] apply_missing_filter: mode = ", mode, ", threshold = ", threshold)
  if (mode == "group" && !is.null(sample_info) && "Group" %in% colnames(sample_info) && !is.null(sample_names_short)) {
    si <- sample_info
    si$short <- extract_sample_names(rownames(si))
    group_vec <- si$Group[match(sample_names_short, si$short)]
    na_mask <- is.na(group_vec)
    if (any(na_mask)) {
      unmatched_count <- sum(na_mask)
      message(sprintf("[DEBUG] apply_missing_filter: %d samples not matched to any group, falling back to global mode", unmatched_count))
      missing_frac <- rowMeans(is.na(data))
      keep <- missing_frac <= threshold
      fallback_triggered <- TRUE
    } else {
      groups <- unique(group_vec)
      keep <- rep(FALSE, nrow(data))
      for (g in groups) {
        cols_in_group <- which(group_vec == g)
        if (length(cols_in_group) > 0) {
          missing_frac_group <- rowMeans(is.na(data[, cols_in_group, drop = FALSE]))
          keep <- keep | (missing_frac_group <= threshold)
        }
      }
      message("[DEBUG] apply_missing_filter (group): kept ", sum(keep), " out of ", nrow(data))
    }
  } else {
    missing_frac <- rowMeans(is.na(data))
    keep <- missing_frac <= threshold
    message("[DEBUG] apply_missing_filter (global): kept ", sum(keep), " out of ", nrow(data))
  }
  result <- data[keep, , drop = FALSE]
  attr(result, "fallback_triggered") <- fallback_triggered
  attr(result, "unmatched_count") <- unmatched_count
  return(result)
}

apply_intensity_filter <- function(data, threshold, min_samples) {
  if (threshold <= 0) return(data)
  above_counts <- apply(data, 1, function(x) sum(x > threshold, na.rm = TRUE))
  keep <- above_counts >= min_samples
  message("[DEBUG] apply_intensity_filter: threshold = ", threshold, ", min_samples = ", min_samples, ", kept = ", sum(keep), " out of ", nrow(data))
  data[keep, , drop = FALSE]
}

safe_pca <- function(mat, scale = TRUE) {
  mat <- as.matrix(mat)
  mat[!is.finite(mat)] <- NA
  mat <- mat[complete.cases(mat), , drop = FALSE]
  if (nrow(mat) < 2) return(NULL)
  row_vars <- apply(mat, 1, var, na.rm = TRUE)
  keep <- row_vars > 1e-12
  mat <- mat[keep, , drop = FALSE]
  if (nrow(mat) < 2) return(NULL)
  tryCatch(prcomp(mat, scale. = scale), error = function(e) NULL)
}