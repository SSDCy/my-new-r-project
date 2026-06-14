# global.R - 全局设置、包加载与工具函数

# =====================================================
# 调试信息：开始加载 global.R
# =====================================================
message("[DEBUG] global.R loading started...")

# =====================================================
# 全局选项（确保上传大小限制等设置最先执行）
# =====================================================
options(shiny.maxRequestSize = 200 * 1024^2)   # 允许上传最大 200MB 文件
options(shiny.fullstacktrace = TRUE)
options(warn = -1)
options(jsonlite.keep_vec_names = FALSE)

# =====================================================
# 包加载
# =====================================================
library(shiny)
library(shinythemes)
library(shinyjs)
library(colourpicker)
library(DT)
library(ggplot2)
library(plotly)
library(data.table)
library(dplyr)
library(tidyr)
library(openxlsx)
library(gridExtra)
library(grid)
library(shinycssloaders)
library(bslib)
library(cowplot)
library(VennDiagram)
library(RColorBrewer)
library(ggVennDiagram)
library(UpSetR)
library(pheatmap)
library(ComplexHeatmap)
library(ComplexUpset)          # 新增：用于绘制彩色 UpSet 图
library(scales)
library(limma)
library(sva)
library(matrixStats)
library(jsonlite)
library(reshape2)
library(writexl)
library(readxl)
library(ggseqlogo)

# =====================================================
# 全局字体设置：所有 ggplot2 图表默认使用 Arial
# =====================================================
theme_set(theme_bw(base_family = "Arial"))
message("[DEBUG] ggplot2 global theme set to Arial")

# =====================================================
# 工具函数
# =====================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

volcano_theme <- function() {
  theme_bw(base_size = 14, base_family = "Arial") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      legend.position = "right",
      legend.title = element_blank(),
      legend.background = element_rect(fill = alpha("white", 0.7)),
      panel.grid = element_blank(),
      legend.text = element_text(size = 12)
    )
}

sanitize_name <- function(name, max_len = 31) {
  if (is.null(name) || is.na(name) || name == "") return("Group")
  cleaned <- gsub("[\\\\/:*?\"<>|\\[\\]]", "_", name)
  cleaned <- gsub("[^[:alnum:][:space:]_.-]", "_", cleaned)
  cleaned <- gsub("_+", "_", cleaned)
  cleaned <- trimws(cleaned)
  if (nchar(cleaned) > max_len) cleaned <- substr(cleaned, 1, max_len)
  if (cleaned == "" || cleaned == "_") cleaned <- "Group"
  cleaned
}

make_unique_names <- function(names) {
  clean <- sapply(names, sanitize_name, USE.NAMES = FALSE)
  if (any(duplicated(clean))) {
    counts <- table(clean)
    dup_names <- names(counts[counts > 1])
    for (nm in dup_names) {
      idx <- which(clean == nm)
      for (k in seq_along(idx)) {
        clean[idx[k]] <- paste0(nm, "_", k)
        if (nchar(clean[idx[k]]) > 31)
          clean[idx[k]] <- paste0(substr(nm, 1, 28), "_", k)
      }
    }
  }
  clean <- sapply(clean, function(x) if (nchar(x) > 31) substr(x, 1, 31) else x)
  clean
}

mixedorder <- function(x) {
  if (length(x) < 1) return(integer(0))
  x <- trimws(as.character(x))
  split <- strsplit(x, "(?<=[0-9])(?=[^0-9])|(?<=[^0-9])(?=[0-9])", perl = TRUE)
  keys <- lapply(split, function(parts) {
    nums <- grepl("^[0-9]+$", parts)
    parts[nums] <- sprintf("%010d", as.integer(parts[nums]))
    parts[!nums] <- tolower(parts[!nums])
    paste(parts, collapse = "\r")
  })
  order(unlist(keys))
}

parse_sample_levels <- function(samples, separators = c("-", "_", ".")) {
  if (length(samples) == 0) return(list())
  sep_counts <- sapply(separators, function(sep) sum(grepl(sep, samples, fixed = TRUE)))
  best_sep <- separators[which.max(sep_counts)]
  if (max(sep_counts) == 0) best_sep <- ""
  split_samples <- lapply(samples, function(s) {
    if (is.na(s) || s == "") return(character(0))
    strsplit(s, best_sep, fixed = TRUE)[[1]]
  })
  valid_lengths <- sapply(split_samples, length)
  valid_lengths <- valid_lengths[valid_lengths > 0]
  if (length(valid_lengths) == 0) return(list())
  max_levels <- max(valid_lengths)
  level_options <- list()
  for (i in 1:max_levels) {
    group_names <- sapply(split_samples, function(parts) {
      if (length(parts) < i) return(paste(parts, collapse = best_sep))
      paste(parts[1:i], collapse = best_sep)
    })
    unique_groups <- unique(group_names[group_names != ""])
    if (length(unique_groups) == 0) next
    level_options[[as.character(i)]] <- list(
      level = i,
      separator = best_sep,
      example = paste(head(unique_groups, 3), collapse = ", "),
      groups = unique_groups
    )
  }
  return(level_options)
}

extract_group_prefix <- function(s, level = "default", separator = "-") {
  if (is.na(s) || s == "") return(s)
  if (level == "default") {
    trimmed <- sub("[-_ .]*[0-9]+$", "", s)
    if (is.na(trimmed) || trimmed == "" || nchar(trimmed) < 2) return(s)
    return(trimmed)
  } else {
    level_num <- as.integer(level)
    if (is.na(level_num) || level_num < 1) return(s)
    parts <- strsplit(s, separator, fixed = TRUE)[[1]]
    if (length(parts) == 0 || length(parts) < level_num) return(s)
    return(paste(parts[1:level_num], collapse = separator))
  }
}

handle_missing <- function(mat, filter_threshold = NULL, impute = FALSE) {
  if (!is.null(filter_threshold)) {
    na_frac <- rowMeans(is.na(mat))
    keep <- na_frac <= filter_threshold
    mat <- mat[keep, , drop = FALSE]
  }
  if (impute) mat[is.na(mat)] <- 0.01
  mat
}

default_colors <- function() {
  list(Up = "#FF0000", Down = "#0000FF",
       Increase = "#C00000", Decrease = "#0945A5", NS = "#7f7e83")
}

textInputMax <- function(inputId, label, value = "", maxlength = 500,
                         placeholder = "", width = NULL,
                         allowed_pattern = NULL) {
  input_tag <- tags$input(
    id = inputId, type = "text", class = "form-control",
    value = value, maxlength = maxlength, placeholder = placeholder
  )
  if (!is.null(allowed_pattern)) {
    input_tag <- tagAppendAttributes(
      input_tag,
      oninput = paste0("this.value = this.value.replace(/",
                       allowed_pattern, "/g, '')")
    )
  }
  tags$div(
    class = "form-group shiny-input-container",
    style = if (!is.null(width)) paste0("width: ", validateCssUnit(width), ";"),
    tags$label(label, `for` = inputId),
    input_tag
  )
}

step_indicator <- function(steps, current_step) {
  step_items <- lapply(seq_along(steps), function(i) {
    is_active <- i == current_step
    is_done <- i < current_step
    class <- if (is_active) "step-active" else if (is_done) "step-done" else "step-future"
    tags$div(class = paste("step-item", class),
             tags$span(class = "step-number", i),
             tags$span(class = "step-text", steps[i])
    )
  })
  tags$div(class = "process-steps", do.call(tagList, step_items))
}

# =====================================================
# 样本名处理：点→下划线
# =====================================================
extract_sample_names <- function(cols) {
  short <- sub("^(LFQ intensity |Intensity )", "", cols, ignore.case = TRUE)
  still_full <- which(short == cols)
  if (length(still_full) > 0) {
    short[still_full] <- sub("^(LFQ[._]?intensity[._]?|Intensity[._]?)", "", cols[still_full], ignore.case = TRUE)
  }
  short <- sub("^[^[:alnum:]]+", "", short)
  short <- sub("[^[:alnum:]]+$", "", short)
  short <- gsub("-", ".", short)        # 先转为点
  short <- gsub("\\.", "_", short)      # 再将点替换为下划线
  message("[DEBUG] extract_sample_names: first 3 = ", paste(head(short, 3), collapse = ", "))
  short
}

standardize_sample_name <- function(x) {
  if (is.null(x) || length(x) == 0) return(character(0))
  original <- x
  x <- as.character(x)
  # 将所有分隔符（-、_、空格、点）统一为下划线，但保留单词完整性
  x <- gsub("[-_\\s\\.]+", "_", x)
  x <- gsub("\\.", "_", x)           # 二次保险
  x <- gsub("^_", "", x)
  x <- gsub("_$", "", x)
  message("[DEBUG] standardize_sample_name: first 3: ",
          paste(head(original, 3), collapse = ", "), " -> ",
          paste(head(x, 3), collapse = ", "))
  x
}

# =====================================================
# 缺失值填充函数
# =====================================================
impute_missing_values <- function(data, method = "knn", k = 10, min_value = 1e-4,
                                  quantile_prob = 0.01) {
  message("[DEBUG] impute_missing_values: method = ", method, ", k = ", k, ", min_value = ", min_value, ", quantile_prob = ", quantile_prob)
  
  if (method == "none") {
    attr(data, "actual_method") <- "none"
    return(data)
  }
  
  data_matrix <- as.matrix(data)
  actual_method <- method
  
  if (method == "knn") {
    if (!requireNamespace("impute", quietly = TRUE)) stop("impute package required for KNN")
    suppressMessages({
      impute_result <- impute::impute.knn(data_matrix, k = k)
    })
    data_matrix <- impute_result$data
    actual_method <- "knn"
  } else if (method == "ppca") {
    if (!requireNamespace("pcaMethods", quietly = TRUE)) stop("pcaMethods package required for PPCA")
    orig_rows <- rownames(data)
    orig_cols <- colnames(data)
    storage.mode(data_matrix) <- "numeric"
    
    na_rows <- which(rowSums(is.na(data_matrix)) == ncol(data_matrix))
    na_cols <- which(colSums(is.na(data_matrix)) == nrow(data_matrix))
    constant_cols <- which(apply(data_matrix, 2, var, na.rm = TRUE) == 0)
    remove_cols <- unique(c(na_cols, constant_cols))
    
    clean <- data_matrix
    if (length(na_rows) > 0) clean <- clean[-na_rows, , drop = FALSE]
    if (length(remove_cols) > 0) clean <- clean[, -remove_cols, drop = FALSE]
    
    if (nrow(clean) < 2 || ncol(clean) < 2) {
      result <- impute_missing_values(data, method = "knn", k = k)
      attr(result, "actual_method") <- "knn (fallback from ppca)"
      return(result)
    }
    
    if (anyNA(clean)) {
      clean <- impute_missing_values(as.data.frame(clean), method = "knn", k = min(k, ncol(clean)-1))
      clean <- as.matrix(clean)
      if (anyNA(clean)) {
        result <- impute_missing_values(data, method = "knn", k = k)
        attr(result, "actual_method") <- "knn (fallback from ppca)"
        return(result)
      }
    }
    
    success <- FALSE
    tryCatch({
      pc <- pcaMethods::pca(clean, method = "ppca", nPcs = min(2, ncol(clean)), scale = "uv", center = TRUE)
      imputed_clean <- as.matrix(pcaMethods::completeObs(pc))
      success <- TRUE
    }, error = function(e) {})
    
    if (!success) {
      result <- impute_missing_values(data, method = "knn", k = k)
      attr(result, "actual_method") <- "knn (fallback from ppca)"
      return(result)
    }
    
    full_matrix <- matrix(NA, nrow = nrow(data_matrix), ncol = ncol(data_matrix))
    rownames(full_matrix) <- orig_rows; colnames(full_matrix) <- orig_cols
    row_idx <- setdiff(seq_len(nrow(data_matrix)), na_rows)
    col_idx <- setdiff(seq_len(ncol(data_matrix)), remove_cols)
    full_matrix[row_idx, col_idx] <- imputed_clean
    if (length(na_rows) > 0) full_matrix[na_rows, ] <- min_value
    if (length(na_cols) > 0) full_matrix[, na_cols] <- min_value
    if (length(constant_cols) > 0) for (j in constant_cols) full_matrix[, j] <- data_matrix[, j]
    data_matrix <- full_matrix
    actual_method <- "ppca"
  } else if (method == "minvalue") {
    all_vals <- data_matrix[!is.na(data_matrix) & data_matrix > 0]
    if (length(all_vals) == 0) { global_min <- 1e-4 } else { global_min <- min(all_vals, na.rm = TRUE) }
    factor <- min_value
    actual_fill <- global_min * factor
    message("[DEBUG] minvalue: global_min=", global_min, ", factor=", factor, ", fill=", actual_fill)
    data_matrix[is.na(data_matrix)] <- actual_fill
    actual_method <- "minvalue"
    attr(data_matrix, "global_min_intensity") <- global_min
    attr(data_matrix, "fill_factor") <- factor
    attr(data_matrix, "actual_fill_value") <- actual_fill
  } else if (method == "quantile") {
    qp <- quantile_prob
    if (is.null(qp) || is.na(qp) || qp <= 0 || qp >= 1) qp <- 0.01
    for (j in seq_len(ncol(data_matrix))) {
      col_vals <- data_matrix[, j]
      na_idx <- which(is.na(col_vals))
      if (length(na_idx) == 0) next
      qval <- quantile(col_vals, probs = qp, na.rm = TRUE)
      if (!is.finite(qval) || length(qval) == 0) qval <- min_value
      data_matrix[na_idx, j] <- qval
    }
    actual_method <- "quantile"
  } else stop("Unknown imputation method.")
  
  rownames(data_matrix) <- rownames(data); colnames(data_matrix) <- colnames(data)
  result <- as.data.frame(data_matrix)
  attr(result, "actual_method") <- actual_method
  message("[DEBUG] impute_missing_values: completed with method ", actual_method)
  result
}

message("[DEBUG] global.R loaded successfully (Arial font global)")