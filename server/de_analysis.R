# server/de_analysis.R

# ---------- 输入验证 ----------
validate_inputs <- reactive({
  if (length(rv$comparisons) == 0) { showNotification("Please define at least one comparison", type = "error", duration = 2); return(FALSE) }
  if (is.na(input$fc_up) || input$fc_up < 1) { showNotification("FC up must be >= 1.0", type = "error", duration = 2); return(FALSE) }
  if (is.na(input$fc_down) || input$fc_down <= 0 || input$fc_down >= 1) { showNotification("FC down must be between 0 and 1", type = "error", duration = 2); return(FALSE) }
  if (is.null(input$p_cut) || !(input$p_cut %in% c("0.05","0.1"))) { showNotification("Please select a valid P-value threshold", type = "error", duration = 2); return(FALSE) }
  if (is.na(input$min_treat_valid) || input$min_treat_valid < 1) { showNotification("Min treat valid replicates invalid", type = "error", duration = 2); return(FALSE) }
  if (is.na(input$min_ctrl_valid) || input$min_ctrl_valid < 1) { showNotification("Min ctrl valid replicates invalid", type = "error", duration = 2); return(FALSE) }
  if (is.na(input$min_rep_ttest) || input$min_rep_ttest < 1) { showNotification("Min t-test replicates invalid", type = "error", duration = 2); return(FALSE) }
  if (is.na(input$min_rep_inc) || input$min_rep_inc < 1) { showNotification("Min Increase replicates invalid", type = "error", duration = 2); return(FALSE) }
  if (is.na(input$min_rep_dec) || input$min_rep_dec < 1) { showNotification("Min Decrease replicates invalid", type = "error", duration = 2); return(FALSE) }
  if (is.na(input$min_unique_pep) || input$min_unique_pep < 1 || input$min_unique_pep > 10) { showNotification("Min unique peptides must be between 1 and 10", type = "error", duration = 2); return(FALSE) }
  TRUE
})

# ---------- 唯一肽段过滤 ----------
filtered_data <- reactive({
  req(norm_data_full())
  nd <- norm_data_full()
  
  unique_col <- grep("^Unique peptides$", colnames(nd), value = TRUE)[1]
  
  # 如果没有该列，或者列中全部为 NA，则跳过过滤
  if (is.na(unique_col) || !unique_col %in% colnames(nd)) {
    showNotification("Info: Unique peptides column not found. Skipping peptide filter.", type = "message", duration = 5)
    return(nd)
  }
  
  nd[[unique_col]] <- as.numeric(nd[[unique_col]])
  
  if (all(is.na(nd[[unique_col]]))) {
    showNotification("Info: Unique peptides column contains only missing values. Skipping peptide filter.", type = "message", duration = 5)
    return(nd)
  }
  
  nd_before <- nd
  nd <- filter(nd, .data[[unique_col]] >= input$min_unique_pep)
  
  if (nrow(nd) == 0) {
    showNotification(
      paste0("Warning: Unique peptides filter (>= ", input$min_unique_pep, ") removed all proteins. The filter has been skipped to allow analysis. Consider lowering the threshold in the Parameters tab."),
      type = "warning", duration = 10
    )
    nd <- nd_before
  }
  
  nd
})

# ---------- 安全条件判断辅助函数 ----------
safe_if <- function(cond) {
  if (isTRUE(cond)) return(TRUE)
  if (identical(cond, TRUE)) return(TRUE)
  return(FALSE)
}

# ---------- 差异分析函数 ----------
run_de_analysis <- function(data_subset, treat_cols, ctrl_cols, fc_up, fc_down, p_cut, stat_method) {
  if (nrow(data_subset) == 0) return(data.frame())
  
  treat_mat <- as.matrix(data_subset[, treat_cols, drop = FALSE])
  control_mat <- as.matrix(data_subset[, ctrl_cols, drop = FALSE])
  treat_mat[treat_mat == 0 | is.na(treat_mat)] <- NA
  control_mat[control_mat == 0 | is.na(control_mat)] <- NA
  
  n_treat <- rowSums(!is.na(treat_mat))
  n_control <- rowSums(!is.na(control_mat))
  n_treat[is.na(n_treat)] <- 0
  n_control[is.na(n_control)] <- 0
  
  min_tv <- input$min_treat_valid
  min_cv <- input$min_ctrl_valid
  min_tt <- input$min_rep_ttest
  min_inc <- input$min_rep_inc
  min_dec <- input$min_rep_dec
  
  mean_treat <- rowMeans(treat_mat, na.rm = TRUE)
  mean_control <- rowMeans(control_mat, na.rm = TRUE)
  mean_treat[!is.finite(mean_treat)] <- NA_real_
  mean_control[!is.finite(mean_control)] <- NA_real_
  
  FC <- mean_treat / mean_control
  FC[!is.finite(FC)] <- NA_real_
  
  p <- rep(NA_real_, nrow(data_subset))
  log2FC <- ifelse(is.finite(FC) & FC > 0, log2(FC), NA_real_)
  
  if (stat_method == "t-test") {
    for (i in 1:nrow(data_subset)) {
      if (isTRUE(n_treat[i] >= 2 && n_control[i] >= 2)) {
        p[i] <- tryCatch(t.test(treat_mat[i, !is.na(treat_mat[i, ])], control_mat[i, !is.na(control_mat[i, ])], var.equal = TRUE)$p.value, error = function(e) NA_real_)
      }
    }
  } else if (stat_method == "wilcoxon") {
    for (i in 1:nrow(data_subset)) {
      if (isTRUE(n_treat[i] >= 2 && n_control[i] >= 2)) {
        p[i] <- tryCatch(wilcox.test(treat_mat[i, !is.na(treat_mat[i, ])], control_mat[i, !is.na(control_mat[i, ])], exact = FALSE)$p.value, error = function(e) NA_real_)
      }
    }
  } else if (stat_method == "limma") {
    mat <- cbind(treat_mat, control_mat)
    log2_mat <- log2(mat + 1)
    group <- factor(c(rep("Treat", ncol(treat_mat)), rep("Ctrl", ncol(control_mat))), levels = c("Ctrl", "Treat"))
    design <- model.matrix(~ group)
    fit <- lmFit(log2_mat, design)
    fit <- eBayes(fit)
    res <- topTable(fit, coef = "groupTreat", number = Inf, sort.by = "none")
    p <- res$P.Value
    log2FC <- res$logFC
    FC <- 2^log2FC
  }
  
  n_treat <- as.numeric(n_treat); n_treat[is.na(n_treat)] <- 0
  n_control <- as.numeric(n_control); n_control[is.na(n_control)] <- 0
  FC <- as.numeric(FC); FC[is.na(FC)] <- 0
  p <- as.numeric(p); p[is.na(p)] <- 1
  log2FC <- as.numeric(log2FC); log2FC[is.na(log2FC)] <- 0
  
  reg <- rep(NA_character_, nrow(data_subset))
  reg_note <- rep(NA_character_, nrow(data_subset))
  
  for (i in seq_len(nrow(data_subset))) {
    nt_i <- n_treat[i]; nc_i <- n_control[i]; fc_i <- FC[i]; p_i <- p[i]
    
    cond_decrease <- safe_if(nt_i == 0 && nc_i >= min_dec)
    cond_increase <- safe_if(nc_i == 0 && nt_i >= min_inc)
    cond_ttest    <- safe_if(nt_i >= min_tt && nc_i >= min_tt)
    
    if (cond_decrease) {
      reg[i] <- "Decrease"
      reg_note[i] <- sprintf("仅在对照组检测到：处理组无有效数据，对照组有效重复数 ≥ %d", min_dec)
    } else if (cond_increase) {
      reg[i] <- "Increase"
      reg_note[i] <- sprintf("仅在处理组检测到：对照组无有效数据，处理组有效重复数 ≥ %d", min_inc)
    } else if (cond_ttest) {
      sig_up   <- safe_if(fc_i > fc_up && p_i < p_cut)
      sig_down <- safe_if(fc_i < fc_down && p_i < p_cut)
      if (sig_up) {
        reg[i] <- "Up"
        reg_note[i] <- sprintf("差异显著上调：FC > %.2f 且 P-value < %.2f", fc_up, p_cut)
      } else if (sig_down) {
        reg[i] <- "Down"
        reg_note[i] <- sprintf("差异显著下调：FC < %.2f 且 P-value < %.2f", fc_down, p_cut)
      } else {
        reg[i] <- "NS"
        if (safe_if(fc_i > fc_up)) reg_note[i] <- sprintf("无显著差异：FC > %.2f 但 P-value ≥ %.2f", fc_up, p_cut)
        else if (safe_if(fc_i < fc_down)) reg_note[i] <- sprintf("无显著差异：FC < %.2f 但 P-value ≥ %.2f", fc_down, p_cut)
        else reg_note[i] <- "无显著差异：统计检验失败"
      }
    } else {
      reg[i] <- "NS"
      if (nt_i < min_tt && nc_i < min_tt) reg_note[i] <- sprintf("无显著差异：两组有效重复数均 < %d", min_tt)
      else if (nt_i < min_tt) reg_note[i] <- sprintf("无显著差异：处理组有效重复数 < %d", min_tt)
      else reg_note[i] <- sprintf("无显著差异：对照组有效重复数 < %d", min_tt)
    }
  }
  
  res <- data_subset %>% mutate(
    n_treat = n_treat, n_control = n_control, mean_treat = mean_treat, mean_control = mean_control,
    FC = FC, Pvalue = p, regulation = reg, regulation_note = reg_note,
    log2FC = log2FC,
    log10P = ifelse(!is.na(Pvalue) & Pvalue < 1, -log10(Pvalue + 1e-10), NA_real_)
  )
  
  inc <- which(res$regulation == "Increase"); dec <- which(res$regulation == "Decrease")
  if (length(inc) > 0) { res$log2FC[inc] <- runif(length(inc), 5.5, 7.5); res$log10P[inc] <- runif(length(inc), 2.5, 6) }
  if (length(dec) > 0) { res$log2FC[dec] <- runif(length(dec), -7.5, -5.5); res$log10P[dec] <- runif(length(dec), 2.5, 6) }
  
  attr(res, "counts") <- list(
    Up = sum(res$regulation == "Up", na.rm = TRUE),
    Down = sum(res$regulation == "Down", na.rm = TRUE),
    Increase = length(inc), Decrease = length(dec),
    NS = sum(res$regulation == "NS", na.rm = TRUE)
  )
  return(res)
}

# ---------- 列匹配 ----------
match_norm_columns <- function(samples, nd_colnames, norm_prefix) {
  exact <- paste0(norm_prefix, samples)
  if (all(exact %in% nd_colnames)) return(exact)
  
  matched <- character(length(samples))
  for (i in seq_along(samples)) {
    pat <- paste0("^", norm_prefix, ".*", samples[i], "$")
    hits <- grep(pat, nd_colnames, value = TRUE)
    if (length(hits) >= 1) matched[i] <- hits[1]
    else matched[i] <- NA_character_
  }
  matched <- na.omit(matched)
  return(matched)
}

# ---------- 所有分析结果 ----------
all_analysis_results <- reactive({
  req(input$baseline_sample)
  req(filtered_data(), rv$comparisons)
  req(validate_inputs())
  
  fcu <- input$fc_up; fcd <- input$fc_down; pc <- as.numeric(input$p_cut)
  stat_method <- input$stat_method
  nd <- filtered_data()
  
  if (nrow(nd) == 0) {
    showNotification("No proteins remain after filtering. Adjust unique peptide threshold or preprocessing.", type = "error", duration = 8)
    return(list(raw = rv$raw_data, clean = rv$clean_data, norm = NULL, filtered = nd, unique_col = NULL, results = list()))
  }
  
  unique_col <- grep("^Unique peptides$", colnames(nd), value = TRUE)[1]
  
  norm_prefix <- get_norm_prefix()
  norm_colnames <- colnames(nd)
  
  results <- list()
  missing_comps <- c()
  
  sorted_comp <- sorted_comps()
  for (comp in sorted_comp) {
    treat_group <- comp$treat; ctrl_group <- comp$ctrl; comp_name <- comp$name
    
    if (!treat_group %in% names(rv$groups) || !ctrl_group %in% names(rv$groups)) {
      missing_comps <- c(missing_comps, paste0(comp_name, " (group missing)"))
      next
    }
    if (length(rv$groups[[treat_group]]) == 0 || length(rv$groups[[ctrl_group]]) == 0) {
      missing_comps <- c(missing_comps, paste0(comp_name, " (group empty)"))
      next
    }
    
    treat_samples <- rv$groups[[treat_group]]
    ctrl_samples <- rv$groups[[ctrl_group]]
    treat_cols <- match_norm_columns(treat_samples, norm_colnames, norm_prefix)
    ctrl_cols <- match_norm_columns(ctrl_samples, norm_colnames, norm_prefix)
    
    if (length(treat_cols) == 0 || length(ctrl_cols) == 0) {
      missing_comps <- c(missing_comps, comp_name)
      next
    }
    
    select_cols <- c("Protein IDs", "Majority protein IDs", "Master protein IDs", all_of(unique_col), all_of(treat_cols), all_of(ctrl_cols))
    s <- select(nd, any_of(select_cols))
    
    res <- tryCatch(
      run_de_analysis(s, treat_cols, ctrl_cols, fcu, fcd, pc, stat_method),
      error = function(e) {
        showNotification(paste("Comparison", comp_name, "failed:", e$message), type = "error", duration = 8)
        return(NULL)
      }
    )
    if (!is.null(res)) {
      results[[comp_name]] <- list(data = res, treat = treat_group, ctrl = ctrl_group, name = comp_name)
    } else {
      missing_comps <- c(missing_comps, paste0(comp_name, " (analysis error)"))
    }
  }
  
  if (length(missing_comps) > 0) {
    showNotification(
      paste("Some comparisons missing or failed:", paste(missing_comps, collapse = ", ")),
      type = "warning", duration = 10, id = "missing_cols"
    )
  }
  
  list(raw = rv$raw_data, clean = rv$clean_data, norm = norm_data_full(), filtered = nd, unique_col = unique_col, results = results)
})

# ---------- 更新下拉框 ----------
observe({
  comp_names <- sapply(sorted_comps(), function(c) c$name)
  if (length(comp_names) > 0) {
    updateSelectInput(session, "selected_comparison", choices = comp_names, selected = comp_names[1])
  } else {
    updateSelectInput(session, "selected_comparison", choices = character(0))
  }
  group_names <- names(rv$groups)
  if (length(group_names) > 0) {
    updateSelectInput(session, "heatmap_groups", choices = group_names, selected = NULL)
  } else {
    updateSelectInput(session, "heatmap_groups", choices = character(0))
  }
})