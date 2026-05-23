# global.R - 全局设置、包加载与工具函数

# =====================================================
# 全局选项（确保上传大小限制等设置最先执行）
# =====================================================
options(shiny.maxRequestSize = 200 * 1024^2)   # 允许上传最大 200MB 文件
options(shiny.fullstacktrace = TRUE)
options(warn = -1)
options("jsonlite.keep_vec_names" = FALSE)

# =====================================================
# 加载所有可能在 UI/Server 模块中直接或间接调用的包
# 将这些包在全局环境中 attach，确保 source 的每个模块都能找到它们
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
library(ComplexHeatmap)   # 热图备用，确保不报错
library(scales)
library(limma)
library(sva)
library(matrixStats)
library(jsonlite)
library(reshape2)
library(writexl)
library(readxl)

# =====================================================
# 工具函数（与之前完全一致）
# =====================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

volcano_theme <- function() {
  theme_bw(base_size = 14) +
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
  if (impute) {
    mat[is.na(mat)] <- 0.01
  }
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

# 流程步骤指示器
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