# =============================================================================
# PlastiPrint — Shiny GUI for Microplastic Chemical Fingerprinting Workflow
# =============================================================================
# SPDX-License-Identifier: AGPL-3.0-only
#
# This Shiny application wraps the complete computational workflow from:
#   "Chemical Fingerprints of New vs. Weathered Microplastics:
#    A Machine Learning Approach"
#
# Each workflow step is exposed as a sequential tab with configurable parameters.
# =============================================================================

# ── 0.  Bootstrap: install required packages silently ──────────────────────────
local({
  cran_pkgs <- c(
    # Shiny + UI
    "shiny", "shinydashboard", "shinyjs", "shinyWidgets", "DT", "shinybusy",
    # Data manipulation
    "tidyverse", "readxl", "writexl", "data.table",
    # Visualization
    "crayon", "grid", "gridExtra", "ggplot2", "pheatmap",
    "viridis", "ggpubr", "scales",
    # Statistics
    "stats", "vegan", "cluster", "factoextra", "FactoMineR",
    "ggforce", "MASS", "moments",
    # Machine Learning
    "caret", "randomForest", "ranger", "e1071", "xgboost",
    # Missing data handling
    "missForest", "pROC", "VIM", "softImpute", "RANN",
    # Parallel processing
    "foreach", "doParallel", "parallel",
    # Dimensionality reduction
    "umap"
  )
  new <- cran_pkgs[!cran_pkgs %in% installed.packages()[, "Package"]]
  if (length(new)) install.packages(new, repos = "https://cloud.r-project.org")
})

# ── 1.  Load libraries ────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(shinyjs)
  library(shinyWidgets)
  library(shinybusy)
  library(DT)
  library(tidyverse)
  library(readxl)
  library(writexl)
  library(data.table)
  library(grid)
  library(gridExtra)
  library(stats)
  library(vegan)
  library(cluster)
  library(factoextra)
  library(FactoMineR)
  library(ggforce)
  library(MASS)
  library(caret)
  library(randomForest)
  library(ranger)
  library(e1071)
  library(xgboost)
  library(pheatmap)
  library(viridis)
  library(missForest)
  library(pROC)
  library(VIM)
  library(softImpute)
  library(ggpubr)
  library(RANN)
  library(foreach)
  library(doParallel)
  library(parallel)
  library(ggplot2)
  library(moments)
  library(scales)
  library(umap)
})

# ── 2.  Source helper functions ───────────────────────────────────────────────
# FIX (bug #1): pattern was "Helper.*Functions.*" (plural) but real file is
# "Helper Function ..." (singular). Make pattern singular-and-plural tolerant.
helper_file <- file.path(getwd(), "Helper Function using only RF_Github_13May2026.R")
if (!file.exists(helper_file)) {
  alt <- list.files(path = c(".", ".."),
                    pattern = "Helper.*Function.*\\.R$",
                    full.names = TRUE, recursive = FALSE)
  if (length(alt)) helper_file <- alt[1]
}
if (file.exists(helper_file)) {
  source(helper_file, local = FALSE)
  message("Loaded helper: ", helper_file)
} else {
  warning("Helper functions file not found. Place it in the app directory.")
}

# ── 2a.  Subcategory-aware overrides (unifies handling, answers A5) ──────────
# process_and_filter from the helper drops Subcategory; data_fusion requires it.
# We override both so each preserves Subcategory IF present in the input.
process_and_filter <- function(df_input,
                               class_col = "Plastic_type",
                               threshold = 0.8,
                               filter_mode = "regular80") {

  has_subcat <- "Subcategory" %in% names(df_input)
  pivot_cols <- c("File", "Source", "Polymer", "Feature", "technique",
                  "Values", class_col)
  if (has_subcat) pivot_cols <- c(pivot_cols, "Subcategory")

  group_cols <- c("File", "Source", "Polymer", "technique", "Feature")
  if (has_subcat) group_cols <- c(group_cols, "Subcategory")

  df_wide <- df_input %>%
    dplyr::select(dplyr::all_of(pivot_cols)) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(group_cols, class_col)))) %>%
    dplyr::summarise(Values = mean(Values, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Feature, values_from = Values) %>%
    tibble::column_to_rownames(var = "File") %>%
    dplyr::relocate(dplyr::any_of(c(class_col, "technique", "Source",
                                    "Polymer", "Subcategory")), .before = 1)

  meta_cols <- intersect(c(class_col, "technique", "Source", "Polymer", "Subcategory"),
                         colnames(df_wide))

  feature_data <- df_wide %>% dplyr::select(dplyr::all_of(setdiff(names(df_wide), meta_cols)))
  kept_features <- character(0)

  if (filter_mode == "none") {
    kept_features <- names(feature_data)
  } else if (filter_mode == "regular80") {
    presence_rate <- colMeans(!is.na(feature_data))
    kept_features <- names(presence_rate[presence_rate >= threshold])
  } else if (filter_mode == "modified80") {
    class_split <- split(df_wide, df_wide[[class_col]])
    class_kept_list <- list()
    for (cls in names(class_split)) {
      cls_idx <- which(df_wide[[class_col]] == cls)
      cls_dat <- feature_data[cls_idx, , drop = FALSE]
      cls_rate <- colMeans(!is.na(cls_dat))
      class_kept_list[[cls]] <- names(cls_rate[cls_rate >= threshold])
    }
    kept_features <- Reduce(union, class_kept_list)
  } else if (filter_mode == "inhouse") {
    classes <- unique(df_wide[[class_col]])
    n_feats <- ncol(feature_data)
    class_presence_mat <- matrix(0, nrow = length(classes), ncol = n_feats)
    class_counts_mat   <- matrix(0, nrow = length(classes), ncol = n_feats)
    for (i in seq_along(classes)) {
      cls <- classes[i]
      cls_rows <- which(df_wide[[class_col]] == cls)
      cls_dat  <- feature_data[cls_rows, , drop = FALSE]
      class_presence_mat[i, ] <- colMeans(!is.na(cls_dat))
      class_counts_mat[i, ]   <- colSums(!is.na(cls_dat))
    }
    global_presence    <- colMeans(!is.na(feature_data))
    max_class_presence <- apply(class_presence_mat, 2, max)
    total_data_points  <- colSums(class_counts_mat)
    n_classes_with_data  <- colSums(class_counts_mat > 0)
    n_classes_with_3plus <- colSums(class_counts_mat >= 3)
    keep_r1 <- (global_presence > 0.10) | (max_class_presence >= 0.90)
    keep_r2 <- total_data_points > 1
    is_diagnostic <- max_class_presence >= 0.90
    keep_r3 <- (n_classes_with_data > 1) | is_diagnostic
    keep_r4 <- (n_classes_with_data < 2) | (n_classes_with_3plus >= 2)
    final_mask <- keep_r1 & keep_r2 & keep_r3 & keep_r4
    kept_features <- colnames(feature_data)[final_mask]
  } else {
    stop("Invalid filter_mode. Choose 'regular80', 'modified80', 'inhouse', or 'none'.")
  }

  final_cols <- unique(c(meta_cols, kept_features))
  df_filtered <- df_wide %>% dplyr::select(dplyr::any_of(final_cols))

  list(original_data = df_wide,
       filtered_data = df_filtered,
       kept_features = kept_features)
}

# data_fusion override: tolerate missing Subcategory
data_fusion <- function(df_input, source_feats, class_col = "Plastic_type") {
  has_subcat <- "Subcategory" %in% names(df_input)
  tech_map <- df_input %>% dplyr::distinct(Feature, technique)
  all_techs <- unique(tech_map$technique)

  sel_cols <- c("File", "Source", "Polymer", "Feature", "technique", "Values", class_col)
  if (has_subcat) sel_cols <- c(sel_cols, "Subcategory")
  grp_cols <- c("File", "Source", "Polymer", "technique", "Feature")
  if (has_subcat) grp_cols <- c(grp_cols, "Subcategory")

  df_wide <- df_input %>%
    dplyr::select(dplyr::all_of(sel_cols)) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(grp_cols, class_col)))) %>%
    dplyr::summarise(Values = mean(Values, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Feature, values_from = Values) %>%
    tibble::column_to_rownames(var = "File") %>%
    dplyr::relocate(dplyr::any_of(c(class_col, "technique", "Source",
                                    "Polymer", "Subcategory")), .before = 1)

  meta_cols <- intersect(c(class_col, "technique", "Source", "Polymer", "Subcategory"),
                         colnames(df_wide))

  if (length(all_techs) > 1) {
    feat_cols  <- setdiff(names(df_wide), meta_cols)
    mat        <- as.matrix(df_wide[, feat_cols])
    row_techs  <- df_wide$technique
    for (tech in all_techs) {
      tech_feats  <- tech_map$Feature[tech_map$technique == tech]
      target_cols <- intersect(colnames(mat), tech_feats)
      if (length(target_cols) > 0) {
        target_rows <- which(row_techs != tech)
        if (length(target_rows) > 0) {
          sub_mat <- mat[target_rows, target_cols, drop = FALSE]
          sub_mat[is.na(sub_mat)] <- 0
          mat[target_rows, target_cols] <- sub_mat
        }
      }
    }
    df_wide[, feat_cols] <- mat
  }

  final_feat_cols <- c()
  for (src in source_feats) {
    f <- setdiff(names(src), meta_cols)
    final_feat_cols <- c(final_feat_cols, f)
  }
  requested_cols <- c(meta_cols, unique(final_feat_cols))
  existing_cols  <- intersect(requested_cols, colnames(df_wide))
  df_wide[, existing_cols, drop = FALSE]
}

# Inline helper that is defined in the Rmd (Subcategory-aware now)
print_filtering_stats <- function(results,
                                  meta_cols = c("Plastic_type", "technique",
                                                "Source", "Polymer",
                                                "Subcategory")) {
  df_final <- results$filtered_data
  df_orig  <- results$original_data
  feats_final <- setdiff(names(df_final), meta_cols)
  feats_orig  <- setdiff(names(df_orig),  meta_cols)
  pct_features <- if (length(feats_orig) > 0)
    (length(feats_final) / length(feats_orig)) * 100 else 0
  lines <- character()
  lines <- c(lines, sprintf("Features retained: %d / %d (%.2f%%)",
                            length(feats_final), length(feats_orig), pct_features))
  if (length(feats_final) > 0) {
    vals_final <- df_final[, feats_final, drop = FALSE]
    pct_missing <- (sum(is.na(vals_final)) / prod(dim(vals_final))) * 100
    lines <- c(lines, sprintf("Missing values after reduction: %.2f%%", pct_missing))
  } else {
    lines <- c(lines, "Missing values after reduction: 0% (no features)")
    vals_final <- data.frame(0)
  }
  signal_final <- sum(vals_final, na.rm = TRUE)
  signal_orig  <- sum(df_orig[, feats_orig, drop = FALSE], na.rm = TRUE)
  pct_signal   <- if (signal_orig > 0) (signal_final / signal_orig) * 100 else 0
  lines <- c(lines, sprintf("Signal retained: %.2f%%", pct_signal))
  paste(lines, collapse = "\n")
}

# ── 2b.  Bug-fix override: pairwise_significance_tests (helper line 1019) ────
# Original helper uses undefined `feature_col` inside loop over `feature_name`.
# This silently broke do_pairwise_test = TRUE. Fixed below.
pairwise_significance_tests <- function(input_df, group_col, start_col_index) {
  df_results <- data.frame(Feature = character(),
                           Category_1 = character(),
                           Category_2 = character(),
                           pval = numeric(),
                           stringsAsFactors = FALSE)
  if (length(unique(input_df[[group_col]])) < 2) {
    return(list(uncorrected = df_results, corrected = character(0)))
  }
  group_pairs  <- utils::combn(unique(as.character(input_df[[group_col]])), 2)
  feature_cols <- names(input_df)[sapply(input_df, is.numeric)]

  for (feature_name in feature_cols) {
    for (col in 1:ncol(group_pairs)) {
      p_1 <- group_pairs[1, col]
      p_2 <- group_pairs[2, col]
      vec1 <- as.numeric(input_df[input_df[[group_col]] == p_1, feature_name, drop = TRUE])
      vec2 <- as.numeric(input_df[input_df[[group_col]] == p_2, feature_name, drop = TRUE])

      if (sum(!is.na(vec1)) >= 3 && sum(!is.na(vec2)) >= 3) {
        # Skip if no variance in either group (Shapiro requires variance)
        if (sd(vec1, na.rm = TRUE) == 0 || sd(vec2, na.rm = TRUE) == 0) next
        p_shapiro1 <- tryCatch(shapiro.test(vec1)$p.value, error = function(e) NA)
        p_shapiro2 <- tryCatch(shapiro.test(vec2)$p.value, error = function(e) NA)
        if (is.na(p_shapiro1) || is.na(p_shapiro2)) next

        if (p_shapiro1 < 0.05 || p_shapiro2 < 0.05) {
          pval_test <- tryCatch(suppressWarnings(wilcox.test(vec1, vec2)$p.value),
                                error = function(e) NA)
        } else {
          pval_test <- tryCatch(t.test(vec1, vec2, var.equal = FALSE)$p.value,
                                error = function(e) NA)
        }
        if (!is.na(pval_test)) {
          df_results <- rbind(df_results,
            data.frame(Feature = feature_name, Category_1 = p_1,
                       Category_2 = p_2, pval = as.numeric(pval_test),
                       stringsAsFactors = FALSE))
        }
      }
    }
  }
  if (nrow(df_results) == 0) {
    return(list(uncorrected = df_results, corrected = character(0)))
  }
  df_results$adjusted_pvalue <- p.adjust(df_results$pval, method = "bonferroni")
  significant <- df_results %>%
    dplyr::filter(adjusted_pvalue < 0.05) %>%
    dplyr::arrange(adjusted_pvalue)
  list(uncorrected = df_results,
       corrected   = unique(significant$Feature))
}

# ── 2c.  Bug-fix override: export_rf_results (helper line 2143) ──────────────
# Original hard-coded the 4 top-N (100/50/25/10) and crashes if user changes them.
# Override generates outputs based on whatever top_importance_results are present.
export_rf_results <- function(results, output_dir = ".", fold, prefix) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  fobj <- results[[fold]]
  if (is.null(fobj)) stop("Fold '", fold, "' not found in results.")

  # ---- Feature lists ----
  if (!is.null(fobj$sig_selected_feats))
    write.csv(fobj$sig_selected_feats,
              file.path(output_dir, paste0(prefix, "_pairwise_sig_features.csv")),
              row.names = FALSE)
  if (!is.null(fobj$rfa_selected_feats))
    write.csv(fobj$rfa_selected_feats,
              file.path(output_dir, paste0(prefix, "_rfa_features.csv")),
              row.names = FALSE)

  # Top-N (driven by whatever counts the user supplied)
  top_results <- fobj$top_importance_results %||% list()
  top_names   <- names(top_results)
  top_ns      <- suppressWarnings(as.integer(sub("^top_", "", top_names)))
  if (length(top_results)) {
    for (i in seq_along(top_results)) {
      tn <- top_ns[i]; nm <- top_names[i]
      write.csv(top_results[[nm]]$selected_features,
                file.path(output_dir, paste0(prefix, "_top ", tn, " RF_features.csv")),
                row.names = FALSE)
    }
  }

  # ---- Metrics ----
  metric_lab <- c("MCC (All features)", "MCC (Pairwise)", "MCC (RFA)")
  metric_val <- c(fobj$all_features_mcc %||% NA,
                  fobj$mcc_sig          %||% NA,
                  fobj$mcc_rfa          %||% NA)
  if (length(top_results)) {
    for (i in seq_along(top_results)) {
      metric_lab <- c(metric_lab, sprintf("MCC (Top %d)", top_ns[i]))
      metric_val <- c(metric_val,
        top_results[[top_names[i]]]$eval_metrics_top$MCC_Multiclass %||% NA)
    }
  }
  metrics_df <- data.frame(metric = metric_lab, value = metric_val,
                           stringsAsFactors = FALSE)
  write.csv(metrics_df, file.path(output_dir, paste0(prefix, "_metrics.csv")),
            row.names = FALSE)

  # ---- Confusion-matrix plots ----
  save_plot <- function(p, suffix) {
    if (is.null(p)) return(invisible(NULL))
    ggplot2::ggsave(file.path(output_dir, paste0(prefix, "_", suffix, ".png")),
                    p, width = 15, height = 8, dpi = 300, units = "in")
  }
  save_plot(fobj$conf_mat_all_feats, "all_feats_confusion_matrix")
  save_plot(fobj$conf_mat_sig,       "pairwise_confusion_matrix")
  save_plot(fobj$conf_mat_rfa,       "RFA_confusion_matrix")
  if (length(top_results)) {
    for (i in seq_along(top_results)) {
      save_plot(top_results[[top_names[i]]]$conf_mat_top,
                sprintf("top %d RF_confusion_matrix", top_ns[i]))
    }
  }

  invisible(list(metrics_df = metrics_df, output_dir = output_dir))
}

# Null-coalesce helper used throughout Steps 7/8
`%||%` <- function(a, b) {
  if (is.null(a) || (length(a) == 1 && is.na(a))) b else a
}

# Aggregator: across folds → mean ± SD metrics table.
# Knows which methods were enabled in the run, and which folds yielded features
# vs which produced nothing. Empty cells become "—" with a Notes column.
aggregate_rf_metrics <- function(rf_res_for_dataset, top_counts, methods_enabled) {
  total_folds <- length(rf_res_for_dataset)

  fmt_row <- function(label, vals, enabled) {
    if (!enabled) {
      data.frame(Metric = label, Mean = "—", SD = "—",
                 Folds = sprintf("0 / %d", total_folds),
                 Notes = "Method not run (disabled by user)",
                 stringsAsFactors = FALSE)
    } else {
      n_used <- sum(!is.na(vals))
      if (n_used == 0) {
        data.frame(Metric = label, Mean = "—", SD = "—",
                   Folds = sprintf("0 / %d", total_folds),
                   Notes = "Not computed (no features selected in any fold)",
                   stringsAsFactors = FALSE)
      } else if (n_used < total_folds) {
        data.frame(Metric = label,
                   Mean  = sprintf("%.4f", mean(vals, na.rm = TRUE)),
                   SD    = sprintf("%.4f", sd(vals, na.rm = TRUE)),
                   Folds = sprintf("%d / %d", n_used, total_folds),
                   Notes = sprintf("%d fold(s) had no features selected",
                                   total_folds - n_used),
                   stringsAsFactors = FALSE)
      } else {
        data.frame(Metric = label,
                   Mean  = sprintf("%.4f", mean(vals, na.rm = TRUE)),
                   SD    = sprintf("%.4f", sd(vals, na.rm = TRUE)),
                   Folds = sprintf("%d / %d", total_folds, total_folds),
                   Notes = "", stringsAsFactors = FALSE)
      }
    }
  }

  rows <- list()
  rows$all <- fmt_row("MCC (All features)",
                      sapply(rf_res_for_dataset,
                             function(f) f$all_features_mcc %||% NA),
                      enabled = TRUE)
  rows$sig <- fmt_row("MCC (Pairwise)",
                      sapply(rf_res_for_dataset, function(f) f$mcc_sig %||% NA),
                      enabled = isTRUE(methods_enabled$do_pairwise))
  rows$rfa <- fmt_row("MCC (RFA)",
                      sapply(rf_res_for_dataset, function(f) f$mcc_rfa %||% NA),
                      enabled = isTRUE(methods_enabled$do_rfa))
  for (tn in top_counts) {
    key <- paste0("top_", tn)
    rows[[key]] <- fmt_row(sprintf("MCC (Top %d)", tn),
      sapply(rf_res_for_dataset, function(f)
        f$top_importance_results[[key]]$eval_metrics_top$MCC_Multiclass %||% NA),
      enabled = isTRUE(methods_enabled$do_topn))
  }
  do.call(rbind, rows)
}

# Per-fold feature counts — used for run-time Console Log summary AND
# to detect "method ran but found nothing" at view time.
rf_fold_feature_counts <- function(rf_res_for_dataset, top_counts, methods_enabled) {
  out <- list()
  out$all <- sapply(rf_res_for_dataset,
                    function(f) length(setdiff(names(f$final_imp_norm_train %||% list()),
                                               c("Plastic_type","Source","Polymer",
                                                 "technique","Subcategory"))))
  if (isTRUE(methods_enabled$do_pairwise))
    out$sig <- sapply(rf_res_for_dataset, function(f) length(f$sig_selected_feats %||% character()))
  if (isTRUE(methods_enabled$do_rfa))
    out$rfa <- sapply(rf_res_for_dataset, function(f) length(f$rfa_selected_feats %||% character()))
  if (isTRUE(methods_enabled$do_topn)) {
    for (tn in top_counts) {
      key <- paste0("top_", tn)
      out[[key]] <- sapply(rf_res_for_dataset, function(f)
        length(f$top_importance_results[[key]]$selected_features %||% character()))
    }
  }
  out
}

# Friendly label for a method key (used in console messages)
rf_method_label <- function(method_key) {
  switch(method_key,
    all_feats = "All features",
    sig       = "Pairwise significance test",
    rfa       = "Recursive Feature Addition (RFA)",
    {
      # top_N
      n <- suppressWarnings(as.integer(sub("^top_", "", method_key)))
      if (is.na(n)) method_key else sprintf("Top %d by importance", n)
    })
}

# Small util — parse comma-separated integer/numeric strings safely
.parse_csv_int <- function(txt, default) {
  v <- suppressWarnings(as.integer(trimws(strsplit(txt, ",")[[1]])))
  v <- v[!is.na(v) & v > 0]
  if (length(v) == 0) default else v
}

# =============================================================================
#                               UI
# =============================================================================
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(title = "PlastiPrint v1.0", titleWidth = 260),

  # ── Sidebar ──────────────────────────────────────────────────────────────────
  dashboardSidebar(
    width = 260,
    sidebarMenu(
      id = "main_tabs",
      menuItem("Welcome",           tabName = "welcome",  icon = icon("home")),
      menuItem("Step 1 — Import",   tabName = "step1",    icon = icon("file-import")),
      menuItem("Step 2 — Blank Sub",tabName = "step2",    icon = icon("filter")),
      menuItem("Step 3 — Source",   tabName = "step3",    icon = icon("tags")),
      menuItem("Step 4 — Metals",   tabName = "step4",    icon = icon("flask")),
      menuItem("Step 5 — Labels",   tabName = "step5",    icon = icon("link")),
      menuItem("Step 6 — Features", tabName = "step6",    icon = icon("cut")),
      menuItem("Step 7 — ML",       tabName = "step7",    icon = icon("brain")),
      menuItem("Step 8 — HCA",      tabName = "step8",    icon = icon("project-diagram")),
      menuItem("Step 9 — PCA/UMAP", tabName = "step9",    icon = icon("circle-nodes"))
    ),
    hr(),
    div(style = "padding: 10px; font-size: 11px; color: #aaa;",
        "PlastiPrint wraps the complete",
        br(), "microplastic fingerprinting",
        br(), "computational workflow.")
  ),

  # ── Body ─────────────────────────────────────────────────────────────────────
  dashboardBody(
    useShinyjs(),
    tags$head(tags$style(HTML("
      .content-wrapper { background-color: #f7f9fc; }
      .box { border-top: 3px solid #3c8dbc; }
      .log-box { max-height: 350px; overflow-y: auto;
                 background: #1e1e1e; color: #d4d4d4;
                 padding: 12px; border-radius: 4px;
                 font-family: 'Fira Code', 'Consolas', monospace;
                 font-size: 12px; white-space: pre-wrap; }
      .nav-step-btn { margin-top: 15px; }
      .param-section { margin-bottom: 10px; }
      .small-label .control-label { font-size: 12px; }
      .advanced-collapse { margin-top: 10px; padding: 8px;
                           background: #eef3f8; border-radius: 4px; }
    "))),

    tabItems(
      # ────────────────────── WELCOME ──────────────────────
      tabItem("welcome",
        fluidRow(
          box(width = 12, status = "primary", solidHeader = TRUE,
              title = "Welcome to PlastiPrint — A R Shiny App for Data Processing & Source Tracking of Complex Contaminant Mixtures",
              h4("This GUI was developed from the publication, titled:",
                 tags$b("Chemical Fingerprints of New vs. Weathered Microplastics")),
              p("This GUI walks you through the complete computational workflow for microplastic chemical fingerprinting. Each tab on the left corresponds to a major step."),
              tags$ol(
                tags$li(tags$b("Step 1:"), " Import raw data (ATD-GC-MS / HPLC-TOF-MS) and group compounds"),
                tags$li(tags$b("Step 2:"), " Blank subtraction & negative-peak removal"),
                tags$li(tags$b("Step 3:"), " Assign source (Store-Bought / Environmental)"),
                tags$li(tags$b("Step 4:"), " Import ICP-MS trace-metal data"),
                tags$li(tags$b("Step 5:"), " Merge label information & create dataset combinations"),
                tags$li(tags$b("Step 6:"), " Feature reduction (Regular 80%, Modified 80%, In-House) & data fusion"),
                tags$li(tags$b("Step 7:"), " Random Forest classification with feature-selection comparison"),
                tags$li(tags$b("Step 8:"), " Hierarchical Cluster Analysis (HCA) on RF-processed data")
              ),
              hr(),
              h4("Getting Started"),
              p("1. Set your ", tags$b("project root directory"), " below.",
                " This folder should contain the ", tags$code("Raw data/"),
                " subfolder and the helper-functions R file."),
              fluidRow(
                column(8, textInput("project_dir", "Project Root Directory:",
                                    value = getwd(), width = "100%")),
                column(4, actionButton("set_dir", "Set Directory",
                                       class = "btn-primary",
                                       style = "margin-top:25px;"))
              ),
              verbatimTextOutput("dir_status")
          )
        )
      ),

      # ────────────────────── STEP 1 ──────────────────────
      tabItem("step1",
        fluidRow(
          box(width = 4, status = "primary", solidHeader = TRUE,
              title = "Step 1: Data Import & Compound Grouping",
              h5("1.1 — Import Raw Data"),
              p("Click the button below to select your data files manually."),
              actionButton("open_import_modal", "Import Data",
                           class = "btn-primary", icon = icon("upload")),
              hr(),
              h5("1.2 — Compound Grouping Parameters"),
              div(class = "param-section",
                tags$b("ATD-GC-MS"),
                numericInput("gc_rtthres1", "RT Threshold:", 0.005,
                             min = 0.0001, step = 0.001, width = "60%"),
                numericInput("gc_mzthres", "m/z Threshold:", 0.05,
                             min = 0.0001, step = 0.01, width = "60%")
              ),
              div(class = "param-section",
                tags$b("HPLC-TOF-MS"),
                numericInput("hplc_split_time", "Split Time (min):", 6,
                             min = 1, step = 0.5, width = "60%"),
                numericInput("hplc_rtthres1", "RT Threshold 1 (≤split):", 0.4,
                             min = 0.01, step = 0.05, width = "60%"),
                numericInput("hplc_rtthres2", "RT Threshold 2 (>split):", 0.28,
                             min = 0.01, step = 0.05, width = "60%"),
                numericInput("hplc_mzthres", "m/z Threshold:", 0.0005,
                             min = 0.00001, step = 0.0001, width = "60%")
              ),
              actionButton("run_step1_group", "Group Compounds",
                           class = "btn-success", icon = icon("object-group")),
              hr(),
              h5("1.3 — Benchmark Removal"),
              p("Removes known benchmark/internal-standard compounds."),
              actionButton("run_step1_bench", "Remove Benchmarks",
                           class = "btn-warning", icon = icon("trash")),
              hr(),
              actionButton("go_step2", "Next → Step 2",
                           class = "btn-info nav-step-btn", icon = icon("arrow-right"))
          ),
          box(width = 8, status = "info", solidHeader = TRUE,
              title = "Output / Preview",
              tabsetPanel(
                tabPanel("Console Log",
                         div(class = "log-box", textOutput("step1_log"))),
                tabPanel("GC Preview",  DTOutput("step1_gc_table")),
                tabPanel("HPLC Preview", DTOutput("step1_hplc_table"))
              )
          )
        )
      ),

      # ────────────────────── STEP 2 ──────────────────────
      tabItem("step2",
        fluidRow(
          box(width = 4, status = "primary", solidHeader = TRUE,
              title = "Step 2: Blank Subtraction",
              p("For each feature, the blank mean is subtracted from sample",
                "areas/heights. Values below 3× blank SD are removed."),
              numericInput("blank_sd_mult", "Blank SD Multiplier:", 3,
                           min = 1, max = 10, step = 0.5),
              actionButton("run_step2", "Run Blank Subtraction",
                           class = "btn-primary", icon = icon("filter")),
              hr(),
              actionButton("go_step3", "Next → Step 3",
                           class = "btn-info nav-step-btn", icon = icon("arrow-right"))
          ),
          box(width = 8, status = "info", solidHeader = TRUE,
              title = "Output",
              tabsetPanel(
                tabPanel("Console Log",
                         div(class = "log-box", textOutput("step2_log"))),
                tabPanel("GC After Blank Sub", DTOutput("step2_gc_table")),
                tabPanel("HPLC After Blank Sub", DTOutput("step2_hplc_table"))
              )
          )
        )
      ),

      # ────────────────────── STEP 3 ──────────────────────
      tabItem("step3",
        fluidRow(
          box(width = 4, status = "primary", solidHeader = TRUE,
              title = "Step 3: Assign Source & Rename",
              p("Labels each sample as ", tags$b("Store-Bought"),
                " or ", tags$b("Environmental"),
                " and assigns an instrument technique tag."),
              p("Environmental samples contain 'USE' in their filename."),
              actionButton("run_step3", "Run Step 3",
                           class = "btn-primary", icon = icon("tags")),
              hr(),
              actionButton("go_step4", "Next → Step 4",
                           class = "btn-info nav-step-btn", icon = icon("arrow-right"))
          ),
          box(width = 8, status = "info", solidHeader = TRUE,
              title = "Output",
              tabsetPanel(
                tabPanel("Console Log",
                         div(class = "log-box", textOutput("step3_log"))),
                tabPanel("GC Data", DTOutput("step3_gc_table")),
                tabPanel("HPLC Data", DTOutput("step3_hplc_table"))
              )
          )
        )
      ),

      # ────────────────────── STEP 4 (REFACTORED) ──────────────────────
      tabItem("step4",
        fluidRow(
          box(width = 4, status = "primary", solidHeader = TRUE,
              title = "Step 4: Import Trace-Metal (ICP-MS) Data",
              p("Reads ICP-MS data from Excel files, selects best",
                "analytical mode per element, removes seawater metals."),
              actionButton("open_icp_modal", "Import ICP-MS Files",
                           class = "btn-primary", icon = icon("upload")),
              hr(),
              textInput("metals_remove", "Seawater metals to discard:",
                        "Na, Mg, K, Ca"),
              actionButton("run_step4", "Process Metals",
                           class = "btn-success", icon = icon("flask")),
              hr(),
              actionButton("go_step5", "Next → Step 5",
                           class = "btn-info nav-step-btn", icon = icon("arrow-right"))
          ),
          box(width = 8, status = "info", solidHeader = TRUE,
              title = "Output",
              tabsetPanel(
                tabPanel("Console Log",
                         div(class = "log-box", textOutput("step4_log"))),
                tabPanel("ICP-MS Data", DTOutput("step4_icp_table"))
              )
          )
        )
      ),

      # ────────────────────── STEP 5 ──────────────────────
      tabItem("step5",
        fluidRow(
          box(width = 4, status = "primary", solidHeader = TRUE,
              title = "Step 5: Combine Label Info & Create Datasets",
              p("Merges sample label info (Plastic_type, Subcategory,",
                "Polymer) and creates all single- and multi-instrument",
                "dataset combinations."),
              fileInput("label_file", "Sample Label Excel File (.xlsx):",
                        accept = c(".xlsx", ".xls")),
              p(tags$small("Leave empty to use default file in Raw data/.")),
              actionButton("run_step5", "Run Step 5",
                           class = "btn-primary", icon = icon("link")),
              hr(),
              actionButton("go_step6", "Next → Step 6",
                           class = "btn-info nav-step-btn", icon = icon("arrow-right"))
          ),
          box(width = 8, status = "info", solidHeader = TRUE,
              title = "Output",
              tabsetPanel(
                tabPanel("Console Log",
                         div(class = "log-box", textOutput("step5_log"))),
                tabPanel("Datasets Summary", DTOutput("step5_summary_table")),
                tabPanel("GC Labeled", DTOutput("step5_gc_table"))
              )
          )
        )
      ),

      # ────────────────────── STEP 6 (DROPDOWN RESTRICTED) ──────────────
      tabItem("step6",
        fluidRow(
          box(width = 4, status = "primary", solidHeader = TRUE,
              title = "Step 6: Feature Reduction & Data Fusion",
              h5("6.1–6.3  Feature Reduction"),
              p(tags$small("Single-instrument datasets only. Fused datasets",
                           "are created in Step 6.4 below.")),
              selectInput("feat_dataset", "Dataset to filter:",
                          choices = c("gc", "hplc", "icp"),
                          selected = "gc"),
              radioButtons("feat_mode", "Filtering method:",
                           choices = c("Regular 80% rule" = "regular80",
                                       "Modified 80% rule" = "modified80",
                                       "In-house comprehensive" = "inhouse"),
                           selected = "inhouse"),
              actionButton("run_step6_filter", "Run Feature Reduction",
                           class = "btn-primary", icon = icon("cut")),
              hr(),
              h5("6.4  Data Fusion"),
              p("Applies in-house filtering to GC, HPLC, ICP",
                "individually, then fuses multi-instrument datasets."),
              actionButton("run_step6_fusion", "Run Data Fusion",
                           class = "btn-success", icon = icon("compress")),
              hr(),
              actionButton("go_step7", "Next → Step 7",
                           class = "btn-info nav-step-btn", icon = icon("arrow-right"))
          ),
          box(width = 8, status = "info", solidHeader = TRUE,
              title = "Output",
              tabsetPanel(
                tabPanel("Console Log",
                         div(class = "log-box", textOutput("step6_log"))),
                tabPanel("Filtering Stats",
                         verbatimTextOutput("step6_stats")),
                tabPanel("Filtered Data Preview", DTOutput("step6_table")),
                tabPanel("Fusion Summary", DTOutput("step6_fusion_table"))
              )
          )
        )
      ),

      # ────────────────────── STEP 7 — RANDOM FOREST ──────────────────
      tabItem("step7",
        fluidRow(
          box(width = 4, status = "primary", solidHeader = TRUE,
              title = "Step 7: Random Forest Classification",
              h5("Datasets to analyze"),
              checkboxGroupInput("rf_datasets", NULL,
                                 choices = c("gc_clean" = "gc_clean",
                                             "hplc_clean" = "hplc_clean",
                                             "icp_clean" = "icp_clean",
                                             "gc_hplc_clean" = "gc_hplc_clean",
                                             "gc_icp_clean" = "gc_icp_clean",
                                             "hplc_icp_clean" = "hplc_icp_clean",
                                             "gc_hplc_icp_clean" = "gc_hplc_icp_clean"),
                                 selected = "gc_clean"),
              p(tags$small("Only datasets cleaned in Step 6 will run.")),
              hr(),

              h5("Train/Test Split Mode"),
              radioButtons("rf_split_mode", NULL,
                           choices = c(
                             "Standard CV split" = "cv",
                             "SB-only train, ENV test (Option 1)" = "sb_only",
                             "SB+ENV train, ENV test (Option 2)"  = "sb_env"
                           ),
                           selected = "sb_env"),
              hr(),

              h5("Feature-Selection Methods"),
              checkboxInput("rf_do_pairwise", "Pairwise significance test", TRUE),
              checkboxInput("rf_do_rfa", "Recursive Feature Addition (RFA)", TRUE),
              checkboxInput("rf_do_topn", "Top-N by importance", TRUE),
              checkboxInput("rf_do_imp_screen",
                            "Impute+Normalize screening (needs FactoMineR)", TRUE),
              hr(),

              h5("Main Hyperparameters"),
              numericInput("rf_seed", "Random seed:", 123, min = 1, step = 1, width="60%"),
              numericInput("rf_min_sample", "min_sample_number:", 1,
                           min = 1, step = 1, width="60%"),
              numericInput("rf_K_outer", "K_outer_max (folds):", 5,
                           min = 2, step = 1, width="60%"),
              numericInput("rf_cores", "Number of cores:",
                           value = max(1, parallel::detectCores(logical = FALSE) - 1),
                           min = 1, step = 1, width="60%"),

              # Advanced settings
              tags$details(class = "advanced-collapse",
                tags$summary(tags$b("Advanced settings ▾")),
                textInput("rf_ntree", "ntree_candidates (comma-sep):",
                          "100, 500, 1000, 2500"),
                textInput("rf_top_counts", "top_importance_counts (comma-sep):",
                          "100, 50, 25, 10"),
                selectInput("rf_metric", "Selection metric:",
                            choices = c("Accuracy", "Kappa"),
                            selected = "Accuracy"),
                numericInput("rf_ntree_screen", "ntree_screen:", 500,
                             min = 100, step = 100, width="60%")
              ),
              hr(),

              actionButton("run_step7", "Run Random Forest",
                           class = "btn-primary",
                           icon = icon("play")),
              hr(),

              h5("View results"),
              selectInput("rf_view_dataset", "Dataset:", choices = NULL),
              selectInput("rf_view_fold",    "Fold (for confusion matrix):",
                          choices = NULL),
              selectInput("rf_view_method",  "Confusion-matrix method:",
                          choices = c("All features"  = "all_feats",
                                      "Pairwise"      = "sig",
                                      "RFA"           = "rfa")),
              hr(),

              h5("Export"),
              textInput("rf_export_dir", "Output directory (relative):",
                        "results_rf"),
              textInput("rf_export_prefix", "Prefix:", "rf"),
              actionButton("rf_export_btn", "Export Results",
                           class = "btn-warning",
                           icon = icon("file-export")),
              hr(),

              actionButton("go_step8", "Next → Step 8",
                           class = "btn-info nav-step-btn",
                           icon = icon("arrow-right"))
          ),

          box(width = 8, status = "info", solidHeader = TRUE,
              title = "Output",
              tabsetPanel(id = "step7_tabs",
                tabPanel("Console Log",
                         div(class = "log-box", textOutput("step7_log"))),
                tabPanel("Confusion Matrices",
                         uiOutput("step7_confmat_display")),
                tabPanel("Metrics (mean ± SD)",
                         DTOutput("step7_metrics_df")),
                tabPanel("Selected Features",
                         DTOutput("step7_features_table")),
                tabPanel("Best Imputation+Normalization",
                         DTOutput("step7_impnorm_table"))
              )
          )
        )
      ),

      # ────────────────────── STEP 8 — HCA ────────────────────────────────
      tabItem("step8",
        fluidRow(
          box(width = 4, status = "primary", solidHeader = TRUE,
              title = "Step 8: Hierarchical Cluster Analysis",
              p("Performs HCA on the train+test data from a Step 7 fold,",
                "labeling samples by Plastic_type-Source."),
              hr(),

              h5("Source data from Step 7"),
              selectInput("hca_dataset", "RF dataset:", choices = NULL),
              selectInput("hca_fold",    "RF fold:",    choices = NULL),
              hr(),

              h5("Distance / Linkage"),
              selectInput("hca_distance", "Distance method:",
                          choices = c("Euclidean"        = "euclidean",
                                      "Bray-Curtis"      = "bray",
                                      "Robust Aitchison" = "robust.aitchison",
                                      "Manhattan"        = "manhattan"),
                          selected = "robust.aitchison"),
              selectInput("hca_linkage", "Linkage method:",
                          choices = c("Average"  = "average",
                                      "Complete" = "complete",
                                      "Single"   = "single",
                                      "Ward.D2"  = "ward.D2"),
                          selected = "average"),
              hr(),

              actionButton("run_step8", "Run HCA",
                           class = "btn-primary", icon = icon("project-diagram")),
              hr(),

              downloadButton("hca_download", "Download Dendrogram (PNG)",
                             class = "btn-warning"),
              hr(),

              actionButton("export_all_btn", "Done / Export All",
                           class = "btn-success nav-step-btn",
                           icon = icon("check-circle")),
              hr(),

              actionButton("go_step9", "Next → Step 9",
                           class = "btn-info nav-step-btn",
                           icon = icon("arrow-right"))
          ),
          box(width = 8, status = "info", solidHeader = TRUE,
              title = "Output",
              tabsetPanel(
                tabPanel("Dendrogram",
                         verbatimTextOutput("step8_ccc"),
                         plotOutput("step8_dendro", height = "700px")),
                tabPanel("Console Log",
                         div(class = "log-box", textOutput("step8_log")))
              )
          )
        )
      ),

      # ────────────────────── STEP 9 — PCA / UMAP ─────────────────────────
      tabItem("step9",
        fluidRow(
          box(width = 4, status = "primary", solidHeader = TRUE,
              title = "Step 9: PCA & UMAP",
              p("Projects the train+test data from a Step 7 fold onto PCA and",
                "UMAP, colored by plastic type (viridis) and shaped by source,",
                "with hierarchical-cluster boundaries (from the same feature set)",
                "overlaid on the UMAP as dashed hulls."),
              hr(),

              h5("Source data from Step 7"),
              selectInput("pcaumap_dataset", "RF dataset:", choices = NULL),
              selectInput("pcaumap_fold",    "RF fold:",    choices = NULL),
              selectInput("pcaumap_feature_set", "Feature set:",
                          choices = c("All features" = "all_feats",
                                      "Pairwise"      = "sig",
                                      "RFA"           = "rfa")),
              hr(),

              h5("Hierarchical-cluster overlay"),
              checkboxInput("pcaumap_show_hulls",
                            "Show HCA cluster hulls on UMAP", value = TRUE),
              numericInput("pcaumap_hclust_k", "Number of HCA clusters (k):",
                           value = 6, min = 2, step = 1, width = "60%"),
              hr(),

              actionButton("run_step9", "Run PCA / UMAP",
                           class = "btn-primary", icon = icon("play")),
              hr(),

              downloadButton("pcaumap_download_pca", "Download PCA (PNG)",
                             class = "btn-warning"),
              downloadButton("pcaumap_download_umap", "Download UMAP (PNG)",
                             class = "btn-warning")
          ),
          box(width = 8, status = "info", solidHeader = TRUE,
              title = "Output",
              tabsetPanel(
                tabPanel("PCA",  plotOutput("step9_pca_plot",  height = "650px")),
                tabPanel("UMAP", plotOutput("step9_umap_plot", height = "650px")),
                tabPanel("Console Log",
                         div(class = "log-box", textOutput("step9_log")))
              )
          )
        )
      )

    ) # end tabItems
  ) # end dashboardBody
) # end dashboardPage

# =============================================================================
#                              SERVER
# =============================================================================
server <- function(input, output, session) {

  # ── Reactive storage for intermediate data ──────────────────────────────────
  rv <- reactiveValues(
    # Step 1
    atdgcms_raw = NULL, atdgcms_blank = NULL, atdgcms_step1 = NULL,
    all_hplc = NULL,    atdgcms_grouped = NULL, hplc_grouped = NULL,
    # Step 2
    atdgcms_step2 = NULL, hplc_step2 = NULL,
    # Step 3
    gc = NULL, hplc = NULL,
    # Step 4 uploads
    icp_round1_path = NULL, icp_round2b1_path = NULL, icp_round2b2_path = NULL,
    rec_unc_1_path = NULL, rec_unc_2_1_path = NULL, rec_unc_2_2_path = NULL,
    icp = NULL,
    # Step 5
    sampinfo = NULL,
    gc_labeled = NULL, hplc_labeled = NULL, icp_labeled = NULL,
    gc_hplc = NULL, gc_icp = NULL, hplc_icp = NULL, gc_hplc_icp = NULL,
    # Step 6
    filter_results = NULL,
    gc_clean = NULL, hplc_clean = NULL, icp_clean = NULL,
    gc_hplc_clean = NULL, gc_icp_clean = NULL, hplc_icp_clean = NULL,
    gc_hplc_icp_clean = NULL,
    fusion_summary = NULL,
    # Step 7
    rf_results = list(),           # list keyed by dataset name
    rf_top_counts = NULL,
    rf_methods_enabled = list(do_pairwise = TRUE, do_rfa = TRUE,
                              do_topn = TRUE, do_imp_screen = TRUE),
    rf_export_dir = NULL,
    # Step 8
    hca_object = NULL, hca_dist = NULL, hca_labels = NULL, hca_ccc = NA,
    # Step 9
    step9_pca_plot = NULL, step9_umap_plot = NULL,

    # Logs
    log_step1 = "Awaiting input...",
    log_step2 = "Run Step 1 first.",
    log_step3 = "Run Step 2 first.",
    log_step4 = "Upload ICP-MS files to begin.",
    log_step5 = "Run Steps 1-4 first.",
    log_step6 = "Run Step 5 first.",
    log_step7 = "Run Step 6 first.",
    log_step8 = "Run Step 7 first.",
    log_step9 = "Run Step 7 first.",
    step6_stats_text = ""
  )

  # ── Helper: append to a log ─────────────────────────────────────────────────
  append_log <- function(key, msg) {
    old <- isolate(rv[[key]])
    if (is.null(old) || old == "Awaiting input..." ||
        grepl("^Run Step|^Upload ", old)) {
      rv[[key]] <- msg
    } else {
      rv[[key]] <- paste0(old, "\n", msg)
    }
  }

  # ── Navigation helpers ──────────────────────────────────────────────────────
  observeEvent(input$go_step2, updateTabItems(session, "main_tabs", "step2"))
  observeEvent(input$go_step3, updateTabItems(session, "main_tabs", "step3"))
  observeEvent(input$go_step4, updateTabItems(session, "main_tabs", "step4"))
  observeEvent(input$go_step5, updateTabItems(session, "main_tabs", "step5"))
  observeEvent(input$go_step6, updateTabItems(session, "main_tabs", "step6"))
  observeEvent(input$go_step7, updateTabItems(session, "main_tabs", "step7"))
  observeEvent(input$go_step8, updateTabItems(session, "main_tabs", "step8"))
  observeEvent(input$go_step9, updateTabItems(session, "main_tabs", "step9"))

  # ── Helper: strip the redundant standalone word "plastic" from a type label ─
  clean_plastic_label <- function(x) {
    x <- gsub("(?i)\\bplastic\\b", "", x, perl = TRUE)
    x <- trimws(gsub("\\s+", " ", x))
    ifelse(is.na(x) | x == "", "Unspecified", x)
  }
  source_shapes <- c("Environmental" = 17, "Store-Bought" = 16)

  # ── Set directory ───────────────────────────────────────────────────────────
  observeEvent(input$set_dir, {
    d <- trimws(input$project_dir)
    if (dir.exists(d)) {
      setwd(d)
      output$dir_status <- renderText(paste("✓ Working directory set to:", getwd()))
    } else {
      output$dir_status <- renderText(paste("✗ Directory not found:", d))
    }
  })

  # ===========================================================================
  #  STEP 1  —  Import, Group, Benchmark  (unchanged from original)
  # ===========================================================================
  observeEvent(input$open_import_modal, {
    showModal(modalDialog(
      title = "Import Data Files",
      size = "l",
      easyClose = FALSE,
      fluidRow(
        column(12,
          h4(icon("flask"), " ATD-GC-MS Data (CSV files)"),
          p("Select CSV files for ATD-GC-MS sample data. You can select multiple files."),
          fileInput("gc_sample_files",
                    label = "GC Sample Files (.csv):",
                    multiple = TRUE,
                    accept = c(".csv", "text/csv"),
                    width = "100%"),
          fileInput("gc_blank_files",
                    label = "GC Blank Files (.csv) [Optional]:",
                    multiple = TRUE,
                    accept = c(".csv", "text/csv"),
                    width = "100%"),
          hr(),
          h4(icon("tint"), " HPLC-TOF-MS Data (XLS files)"),
          p("Select XLS files for HPLC-TOF-MS data. You can select multiple files."),
          fileInput("hplc_files",
                    label = "HPLC Files (.xls):",
                    multiple = TRUE,
                    accept = c(".xls", ".xlsx", "application/vnd.ms-excel"),
                    width = "100%"),
          hr(),
          div(style = "background-color: #f5f5f5; padding: 10px; border-radius: 5px;",
            h5(icon("info-circle"), " File Selection Summary:"),
            textOutput("import_summary")
          )
        )
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_import", "Import Selected Files",
                     class = "btn-primary", icon = icon("check"))
      )
    ))
  })

  output$import_summary <- renderText({
    gc_samples <- if (!is.null(input$gc_sample_files)) nrow(input$gc_sample_files) else 0
    gc_blanks  <- if (!is.null(input$gc_blank_files))  nrow(input$gc_blank_files)  else 0
    hplc       <- if (!is.null(input$hplc_files))      nrow(input$hplc_files)      else 0
    paste0("GC Sample files: ", gc_samples, " | ",
           "GC Blank files: ",  gc_blanks,  " | ",
           "HPLC files: ",      hplc)
  })

  observeEvent(input$confirm_import, {
    removeModal()
    rv$log_step1 <- ""
    tryCatch({
      append_log("log_step1", "── Importing ATD-GC-MS data ──")
      read_and_add_filename <- function(file_path, file_name) {
        data.table::fread(file_path, showProgress = FALSE, data.table = FALSE,
                          check.names = TRUE) %>%
          mutate(File = file_name)
      }

      if (!is.null(input$gc_sample_files) && nrow(input$gc_sample_files) > 0) {
        gc_sample_list <- lapply(1:nrow(input$gc_sample_files), function(i) {
          read_and_add_filename(
            input$gc_sample_files$datapath[i],
            gsub(".csv", "", input$gc_sample_files$name[i])
          )
        })
        atdgcms_raw <- data.table::rbindlist(gc_sample_list, use.names = TRUE, fill = TRUE) %>%
          as.data.frame() %>%
          dplyr::select(-any_of(c("Start","End","Width","Base.Peak",
                                  "Cpd","Label","Height","Ions"))) %>%
          dplyr::mutate(File = gsub("_", "-", File)) %>%
          mutate(type = "Sample")
        append_log("log_step1", sprintf("  GC samples: %d rows from %d files",
                                        nrow(atdgcms_raw), nrow(input$gc_sample_files)))
      } else {
        atdgcms_raw <- data.frame()
        append_log("log_step1", "  No GC sample files selected.")
      }

      if (!is.null(input$gc_blank_files) && nrow(input$gc_blank_files) > 0) {
        gc_blank_list <- lapply(1:nrow(input$gc_blank_files), function(i) {
          data.table::fread(input$gc_blank_files$datapath[i],
                           showProgress = FALSE, check.names = TRUE) %>%
            mutate(File = gsub(".csv", "", input$gc_blank_files$name[i]))
        })
        atdgcms_blank <- data.table::rbindlist(gc_blank_list, use.names = TRUE, fill = TRUE) %>%
          as.data.frame() %>%
          dplyr::select(any_of(c("Area","File","m.z","RT"))) %>%
          mutate(type = "Blanks")
        append_log("log_step1", sprintf("  GC blanks : %d rows from %d files",
                                        nrow(atdgcms_blank), nrow(input$gc_blank_files)))
      } else {
        atdgcms_blank <- data.frame()
        append_log("log_step1", "  No GC blank files selected.")
      }

      if (nrow(atdgcms_raw) > 0 || nrow(atdgcms_blank) > 0) {
        atdgcms_step1 <- bind_rows(atdgcms_raw, atdgcms_blank) %>%
          arrange(RT) %>%
          mutate(simplified_file = map_chr(File, ~ {
            parts <- str_split_1(.x, pattern = "-")
            if (length(parts) >= 5) paste(parts[c(4, 5)], collapse = "-") else .x
          }))
        rv$atdgcms_raw   <- atdgcms_raw
        rv$atdgcms_blank <- atdgcms_blank
        rv$atdgcms_step1 <- atdgcms_step1
      }

      # HPLC import
      append_log("log_step1", "\n── Importing HPLC-TOF-MS data ──")
      if (!is.null(input$hplc_files) && nrow(input$hplc_files) > 0) {
        read_hplc_file <- function(file_path, file_name) {
          df <- readxl::read_xls(file_path, skip = 1) %>% as.data.frame()
          if ("m/z" %in% names(df)) df <- df %>% dplyr::rename(m.z = `m/z`)
          df %>%
            dplyr::select(any_of(c("m.z", "RT", "Height", "Area"))) %>%
            mutate(
              File = gsub("_", "-", gsub("\\.(xls|xlsx)$", "", file_name)),
              type = ifelse(str_detect(file_name, regex("blank", ignore_case = TRUE)),
                            "Blanks", "Sample"),
              Day = case_when(
                str_detect(file_name, "Day0") ~ "0",
                str_detect(file_name, "Day2") ~ "2",
                str_detect(file_name, "Day15") ~ "15",
                TRUE ~ "34"
              ),
              Replicate = ifelse(str_detect(file_name, "R1"), "R1", "R2"),
              Batch_number = "manual"
            )
        }
        hplc_list <- lapply(1:nrow(input$hplc_files), function(i) {
          tryCatch(read_hplc_file(input$hplc_files$datapath[i], input$hplc_files$name[i]),
                   error = function(e) {
                     warning(paste("Error reading file:", input$hplc_files$name[i], "-", e$message))
                     data.frame()
                   })
        })
        all_hplc <- data.table::rbindlist(hplc_list, use.names = TRUE, fill = TRUE) %>%
          as.data.frame() %>%
          arrange(RT) %>%
          mutate(simplified_file = map_chr(File, ~ {
            parts <- str_split_1(.x, "-")
            if (length(parts) >= 2) paste(parts[1:min(2, length(parts))], collapse = "-") else .x
          }))
        rv$all_hplc <- all_hplc
        append_log("log_step1", sprintf("  HPLC total: %d rows from %d files",
                                        nrow(all_hplc), nrow(input$hplc_files)))
      } else {
        append_log("log_step1", "  No HPLC files selected.")
      }
      append_log("log_step1", "\n✓ Data import complete.")
    },
    error = function(e) append_log("log_step1", paste("✗ ERROR:", conditionMessage(e))))
  })

  observeEvent(input$run_step1_group, {
    has_gc   <- !is.null(rv$atdgcms_step1) && nrow(rv$atdgcms_step1) > 0
    has_hplc <- !is.null(rv$all_hplc)      && nrow(rv$all_hplc)      > 0
    if (!has_gc && !has_hplc) {
      append_log("log_step1", "\n⚠ No data imported yet. Please import data first.")
      return()
    }
    tryCatch({
      append_log("log_step1", "\n── Grouping compounds ──")
      if (has_gc) {
        rv$atdgcms_grouped <- grouping_comp(
          rv$atdgcms_step1,
          rtthres1 = input$gc_rtthres1, mzthres = input$gc_mzthres, type = "ATDGCMS")
        append_log("log_step1", sprintf("  GC grouped: %d unique features",
                                        n_distinct(rv$atdgcms_grouped$Feature, na.rm = TRUE)))
      } else append_log("log_step1", "  GC: No data to group (skipped)")
      if (has_hplc) {
        rv$hplc_grouped <- grouping_comp(
          rv$all_hplc,
          split_time = input$hplc_split_time,
          rtthres1   = input$hplc_rtthres1, rtthres2 = input$hplc_rtthres2,
          mzthres    = input$hplc_mzthres,  type     = "HPLCTOFMS")
        append_log("log_step1", sprintf("  HPLC grouped: %d unique features",
                                        n_distinct(rv$hplc_grouped$Feature, na.rm = TRUE)))
      } else append_log("log_step1", "  HPLC: No data to group (skipped)")
      append_log("log_step1", "✓ Grouping complete.")
    },
    error = function(e) append_log("log_step1", paste("✗ ERROR:", conditionMessage(e))))
  })

  observeEvent(input$run_step1_bench, {
    has_gc   <- !is.null(rv$atdgcms_grouped) && nrow(rv$atdgcms_grouped) > 0
    has_hplc <- !is.null(rv$hplc_grouped)    && nrow(rv$hplc_grouped)    > 0
    if (!has_gc && !has_hplc) {
      append_log("log_step1", "\n⚠ No grouped data available. Please run 'Group Compounds' first.")
      return()
    }
    tryCatch({
      append_log("log_step1", "\n── Removing benchmarks ──")
      if (has_gc) {
        gc_g <- rv$atdgcms_grouped
        n_before_gc <- nrow(gc_g)
        idx <- which(gc_g$m.z >= 98 & gc_g$m.z <= 98.5 &
                     gc_g$RT > 13.4 & gc_g$RT < 13.5)
        if (length(idx) > 0) gc_g <- gc_g[-idx, ]
        rv$atdgcms_grouped <- gc_g
        append_log("log_step1", sprintf("  GC: removed %d benchmark rows", n_before_gc - nrow(gc_g)))
      } else append_log("log_step1", "  GC: No data to process (skipped)")
      if (has_hplc) {
        hplc_g <- rv$hplc_grouped
        n_before_hplc <- nrow(hplc_g)
        hplc_benchmark_mz <- c(156.0957, 198.1723, 246.2568, 305.1956, 342.1722)
        hplc_idx <- c()
        for (mz in hplc_benchmark_mz) {
          idx <- which(abs(hplc_g$m.z - mz) <= 0.0005)
          hplc_idx <- c(hplc_idx, idx)
        }
        if (length(hplc_idx) > 0) hplc_g <- hplc_g[-hplc_idx, ]
        rv$hplc_grouped <- hplc_g
        append_log("log_step1", sprintf("  HPLC: removed %d benchmark rows",
                                        n_before_hplc - nrow(hplc_g)))
      } else append_log("log_step1", "  HPLC: No data to process (skipped)")
      append_log("log_step1", "✓ Benchmarks removed.")
    },
    error = function(e) append_log("log_step1", paste("✗ ERROR:", conditionMessage(e))))
  })

  output$step1_log <- renderText(rv$log_step1)
  output$step1_gc_table <- renderDT({
    if (is.null(rv$atdgcms_grouped) || nrow(rv$atdgcms_grouped) == 0)
      return(datatable(data.frame(Message = "No GC data available. Import and group data first."),
                       options = list(dom = 't')))
    datatable(head(rv$atdgcms_grouped, 200), options = list(scrollX = TRUE, pageLength = 10))
  })
  output$step1_hplc_table <- renderDT({
    if (is.null(rv$hplc_grouped) || nrow(rv$hplc_grouped) == 0)
      return(datatable(data.frame(Message = "No HPLC data available. Import and group data first."),
                       options = list(dom = 't')))
    datatable(head(rv$hplc_grouped, 200), options = list(scrollX = TRUE, pageLength = 10))
  })

  # ===========================================================================
  #  STEP 2  —  Blank Subtraction  (unchanged)
  # ===========================================================================
  observeEvent(input$run_step2, {
    has_gc   <- !is.null(rv$atdgcms_grouped) && nrow(rv$atdgcms_grouped) > 0
    has_hplc <- !is.null(rv$hplc_grouped)    && nrow(rv$hplc_grouped)    > 0
    if (!has_gc && !has_hplc) {
      rv$log_step2 <- "⚠ No grouped data available. Please complete Step 1 first."
      return()
    }
    rv$log_step2 <- ""
    tryCatch({
      sd_mult <- input$blank_sd_mult
      append_log("log_step2", sprintf("── Blank subtraction (SD multiplier = %.1f) ──", sd_mult))
      if (has_gc) {
        gc_adj <- rv$atdgcms_grouped %>%
          dplyr::group_by(Feature) %>%
          dplyr::mutate(
            has_blank  = any(type == "Blanks"),
            blank_mean = mean(Area[type == "Blanks"], na.rm = TRUE),
            blank_sd   = sd(Area[type == "Blanks"], na.rm = TRUE)) %>%
          ungroup() %>%
          dplyr::filter(type != "Blanks") %>%
          dplyr::mutate(
            Area = if_else(has_blank, Area - blank_mean, Area),
            threshold = if_else(has_blank, sd_mult * blank_sd, 0)) %>%
          dplyr::filter(Area > 0, Area > threshold) %>%
          dplyr::select(-has_blank, -blank_mean, -blank_sd, -threshold)
        rv$atdgcms_step2 <- gc_adj
        append_log("log_step2", sprintf("  GC after blank sub: %d rows", nrow(gc_adj)))
      } else append_log("log_step2", "  GC: No data to process (skipped)")

      if (has_hplc) {
        hplc_adj <- rv$hplc_grouped %>%
          dplyr::group_by(Batch_number, Feature) %>%
          dplyr::mutate(
            has_blank  = any(type == "Blanks"),
            blank_mean = mean(Height[type == "Blanks"], na.rm = TRUE),
            blank_sd   = sd(Height[type == "Blanks"], na.rm = TRUE)) %>%
          ungroup() %>%
          dplyr::filter(type != "Blanks") %>%
          dplyr::mutate(
            Height = if_else(has_blank, Height - blank_mean, Height),
            threshold = if_else(has_blank, sd_mult * blank_sd, 0)) %>%
          dplyr::filter(Height > 0, Height > threshold) %>%
          dplyr::select(-has_blank, -blank_mean, -blank_sd, -threshold)
        rv$hplc_step2 <- hplc_adj
        append_log("log_step2", sprintf("  HPLC after blank sub: %d rows", nrow(hplc_adj)))
      } else append_log("log_step2", "  HPLC: No data to process (skipped)")
      append_log("log_step2", "✓ Blank subtraction complete.")
    },
    error = function(e) append_log("log_step2", paste("✗ ERROR:", conditionMessage(e))))
  })

  output$step2_log <- renderText(rv$log_step2)
  output$step2_gc_table <- renderDT({
    if (is.null(rv$atdgcms_step2) || nrow(rv$atdgcms_step2) == 0)
      return(datatable(data.frame(Message = "No GC data available. Run blank subtraction first."),
                       options = list(dom = 't')))
    datatable(head(rv$atdgcms_step2, 200), options = list(scrollX = TRUE, pageLength = 10))
  })
  output$step2_hplc_table <- renderDT({
    if (is.null(rv$hplc_step2) || nrow(rv$hplc_step2) == 0)
      return(datatable(data.frame(Message = "No HPLC data available. Run blank subtraction first."),
                       options = list(dom = 't')))
    datatable(head(rv$hplc_step2, 200), options = list(scrollX = TRUE, pageLength = 10))
  })

  # ===========================================================================
  #  STEP 3  —  Source Assignment  (unchanged)
  # ===========================================================================
  observeEvent(input$run_step3, {
    has_gc   <- !is.null(rv$atdgcms_step2) && nrow(rv$atdgcms_step2) > 0
    has_hplc <- !is.null(rv$hplc_step2)    && nrow(rv$hplc_step2)    > 0
    if (!has_gc && !has_hplc) {
      rv$log_step3 <- "⚠ No data available from Step 2. Please complete Step 2 first."
      return()
    }
    rv$log_step3 <- ""
    tryCatch({
      append_log("log_step3", "── Assigning Source & Technique ──")
      if (has_gc) {
        gc <- rv$atdgcms_step2 %>%
          dplyr::filter(!is.na(Feature)) %>%
          dplyr::select(File, Feature, m.z, RT, Area) %>%
          mutate(Source = ifelse(str_detect(File, "USE"), "Environmental", "Store-Bought")) %>%
          dplyr::rename(Values = Area) %>%
          mutate(technique = "GC")
        rv$gc <- gc
        append_log("log_step3", sprintf("  GC:   %d rows  |  SB: %d  ENV: %d",
          nrow(gc), sum(gc$Source == "Store-Bought"), sum(gc$Source == "Environmental")))
      } else append_log("log_step3", "  GC: No data to process (skipped)")
      if (has_hplc) {
        hplc <- rv$hplc_step2 %>%
          dplyr::filter(!is.na(Feature)) %>%
          dplyr::select(File, Feature, m.z, RT, Height) %>%
          dplyr::mutate(Source = ifelse(str_detect(File, "USE"), "Environmental", "Store-Bought")) %>%
          dplyr::rename(Values = Height) %>%
          mutate(technique = "HPLC")
        rv$hplc <- hplc
        append_log("log_step3", sprintf("  HPLC: %d rows  |  SB: %d  ENV: %d",
          nrow(hplc), sum(hplc$Source == "Store-Bought"), sum(hplc$Source == "Environmental")))
      } else append_log("log_step3", "  HPLC: No data to process (skipped)")
      append_log("log_step3", "✓ Step 3 complete.")
    },
    error = function(e) append_log("log_step3", paste("✗ ERROR:", conditionMessage(e))))
  })

  output$step3_log <- renderText(rv$log_step3)
  output$step3_gc_table <- renderDT({
    if (is.null(rv$gc) || nrow(rv$gc) == 0)
      return(datatable(data.frame(Message = "No GC data available. Run Step 3 first."),
                       options = list(dom = 't')))
    datatable(head(rv$gc, 200), options = list(scrollX = TRUE, pageLength = 10))
  })
  output$step3_hplc_table <- renderDT({
    if (is.null(rv$hplc) || nrow(rv$hplc) == 0)
      return(datatable(data.frame(Message = "No HPLC data available. Run Step 3 first."),
                       options = list(dom = 't')))
    datatable(head(rv$hplc, 200), options = list(scrollX = TRUE, pageLength = 10))
  })

  # ===========================================================================
  #  STEP 4  —  ICP-MS (REFACTORED to use fileInputs via modal)
  # ===========================================================================
  observeEvent(input$open_icp_modal, {
    showModal(modalDialog(
      title = "Import ICP-MS Files",
      size = "l",
      easyClose = FALSE,
      fluidRow(
        column(12,
          h4(icon("flask"), " Raw ICP-MS Data (3 XLSX files)"),
          p("These are the three raw concentration tables (round 1, round 2 batch 1, round 2 batch 2)."),
          fileInput("icp_round1_file",
                    "Round 1 raw data (.xlsx):",
                    multiple = FALSE,
                    accept = c(".xlsx", ".xls"),
                    width = "100%"),
          fileInput("icp_round2_batch1_file",
                    "Round 2 Batch 1 raw data (.xlsx):",
                    multiple = FALSE,
                    accept = c(".xlsx", ".xls"),
                    width = "100%"),
          fileInput("icp_round2_batch2_file",
                    "Round 2 Batch 2 raw data (.xlsx):",
                    multiple = FALSE,
                    accept = c(".xlsx", ".xls"),
                    width = "100%"),
          hr(),
          h4(icon("clipboard-check"), " Recovery / Uncertainty Files (3 XLSX files)"),
          p("Each contains the Element+Modus / Recovery [%] / U(k=2) columns matching a raw file."),
          fileInput("rec_unc_1_file",
                    "Round 1 recovery+uncertainty (.xlsx):",
                    multiple = FALSE,
                    accept = c(".xlsx", ".xls"),
                    width = "100%"),
          fileInput("rec_unc_2_1_file",
                    "Round 2 Batch 1 recovery+uncertainty (.xlsx):",
                    multiple = FALSE,
                    accept = c(".xlsx", ".xls"),
                    width = "100%"),
          fileInput("rec_unc_2_2_file",
                    "Round 2 Batch 2 recovery+uncertainty (.xlsx):",
                    multiple = FALSE,
                    accept = c(".xlsx", ".xls"),
                    width = "100%"),
          hr(),
          div(style = "background-color: #f5f5f5; padding: 10px; border-radius: 5px;",
            h5(icon("info-circle"), " File Selection Summary:"),
            textOutput("icp_summary")
          )
        )
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_icp_import", "Save Selections",
                     class = "btn-primary", icon = icon("check"))
      )
    ))
  })

  output$icp_summary <- renderText({
    s <- function(x) if (!is.null(x)) basename(x$name) else "(none)"
    paste0(
      "Raw 1: ",         s(input$icp_round1_file),         "\n",
      "Raw 2-B1: ",      s(input$icp_round2_batch1_file),  "\n",
      "Raw 2-B2: ",      s(input$icp_round2_batch2_file),  "\n",
      "Recov 1: ",       s(input$rec_unc_1_file),          "\n",
      "Recov 2-B1: ",    s(input$rec_unc_2_1_file),        "\n",
      "Recov 2-B2: ",    s(input$rec_unc_2_2_file))
  })

  observeEvent(input$confirm_icp_import, {
    rv$icp_round1_path    <- input$icp_round1_file$datapath        %||% NULL
    rv$icp_round2b1_path  <- input$icp_round2_batch1_file$datapath %||% NULL
    rv$icp_round2b2_path  <- input$icp_round2_batch2_file$datapath %||% NULL
    rv$rec_unc_1_path     <- input$rec_unc_1_file$datapath         %||% NULL
    rv$rec_unc_2_1_path   <- input$rec_unc_2_1_file$datapath       %||% NULL
    rv$rec_unc_2_2_path   <- input$rec_unc_2_2_file$datapath       %||% NULL
    removeModal()
    rv$log_step4 <- "✓ ICP-MS file paths registered. Click 'Process Metals' to run."
  })

  observeEvent(input$run_step4, {
    needed <- list(
      "Round 1 raw"    = rv$icp_round1_path,
      "Round 2 B1 raw" = rv$icp_round2b1_path,
      "Round 2 B2 raw" = rv$icp_round2b2_path,
      "Recovery 1"     = rv$rec_unc_1_path,
      "Recovery 2-B1"  = rv$rec_unc_2_1_path,
      "Recovery 2-B2"  = rv$rec_unc_2_2_path)
    missing <- names(needed)[sapply(needed, is.null)]
    if (length(missing)) {
      rv$log_step4 <- paste("⚠ Missing uploads:", paste(missing, collapse = ", "),
                            "\nClick 'Import ICP-MS Files' first.")
      return()
    }
    rv$log_step4 <- ""
    tryCatch({
      append_log("log_step4", "── Importing ICP-MS trace-metal data ──")
      read_trace_metal <- function(path) {
        readxl::read_excel(path) %>% column_to_rownames(var = "File")
      }
      remove_H2_HMI_modus <- function(icpms_df) {
        strings <- colnames(icpms_df)
        pattern <- "^(.*\\[ )([A-Za-z0-9 ]+)( \\])$"
        keep_cols <- strings[!str_detect(str_match(strings, pattern)[,3], "H2 HMI")]
        icpms_df %>% dplyr::select(all_of(keep_cols))
      }
      get_final_metal_names <- function(recovery_uncertainty, target_cols) {
        recovery_uncertainty <- recovery_uncertainty %>% filter(`Element+Modus` %in% target_cols)
        first_numbers <- sub(" .*", "", recovery_uncertainty$`Element+Modus`)
        kept_names <- c()
        for (num in unique(first_numbers)) {
          rows <- recovery_uncertainty[first_numbers == num, ]
          if (nrow(rows) == 1) kept_names <- c(kept_names, rows$`Element+Modus`)
          else {
            best_row <- rows[which.max(rows$`Recovery [%]`), ]
            if (nrow(best_row) == 0) best_row <- rows[which.min(rows$`U(k=2)`), ]
            kept_names <- c(kept_names, best_row$`Element+Modus`)
          }
        }
        kept_names
      }

      icpms_r1   <- read_trace_metal(rv$icp_round1_path)   %>% remove_H2_HMI_modus()
      icpms_r2b1 <- read_trace_metal(rv$icp_round2b1_path) %>% remove_H2_HMI_modus()
      icpms_r2b2 <- read_trace_metal(rv$icp_round2b2_path) %>% remove_H2_HMI_modus()

      rec_unc_1   <- readxl::read_excel(rv$rec_unc_1_path)
      rec_unc_2_1 <- readxl::read_excel(rv$rec_unc_2_1_path)
      rec_unc_2_2 <- readxl::read_excel(rv$rec_unc_2_2_path)

      fn1   <- get_final_metal_names(rec_unc_1,   colnames(icpms_r1))
      fn2_1 <- get_final_metal_names(rec_unc_2_1, colnames(icpms_r2b1))
      fn2_2 <- get_final_metal_names(rec_unc_2_2, colnames(icpms_r2b2))

      icpms_r1   <- icpms_r1   %>% dplyr::select(all_of(fn1))
      icpms_r2b1 <- icpms_r2b1 %>% dplyr::select(all_of(fn2_1))
      icpms_r2b2 <- icpms_r2b2 %>% dplyr::select(all_of(fn2_2))

      common_cols <- Reduce(intersect,
                            list(names(icpms_r1), names(icpms_r2b1), names(icpms_r2b2)))
      combined_icpms <- rbind(icpms_r1[common_cols],
                              icpms_r2b1[common_cols],
                              icpms_r2b2[common_cols])

      metals_remove  <- trimws(unlist(strsplit(input$metals_remove, ",")))
      cols_to_remove <- unique(unlist(sapply(metals_remove,
                              function(x) grep(x, colnames(combined_icpms)))))
      if (length(cols_to_remove) > 0) combined_icpms <- combined_icpms[, -cols_to_remove]

      icp <- combined_icpms %>%
        rownames_to_column(var = "File") %>%
        pivot_longer(cols = 2:ncol(.), names_to = "Feature", values_to = "Values") %>%
        mutate(technique = "ICP",
               Source = ifelse(str_detect(File, "USE"), "Environmental", "Store-Bought"),
               Values = as.numeric(Values))
      rv$icp <- icp
      append_log("log_step4", sprintf("  ICP-MS: %d rows, %d elements",
                                      nrow(icp), n_distinct(icp$Feature)))
      append_log("log_step4", "✓ ICP-MS import complete.")
    },
    error = function(e) append_log("log_step4", paste("✗ ERROR:", conditionMessage(e))))
  })

  output$step4_log <- renderText(rv$log_step4)
  output$step4_icp_table <- renderDT({
    if (is.null(rv$icp) || nrow(rv$icp) == 0)
      return(datatable(data.frame(Message = "No ICP data yet. Upload + process first."),
                       options = list(dom = 't')))
    datatable(head(rv$icp, 200), options = list(scrollX = TRUE, pageLength = 10))
  })

  # ===========================================================================
  #  STEP 5  —  Label merge & Combinations  (unchanged)
  # ===========================================================================
  observeEvent(input$run_step5, {
    has_gc   <- !is.null(rv$gc)   && nrow(rv$gc)   > 0
    has_hplc <- !is.null(rv$hplc) && nrow(rv$hplc) > 0
    has_icp  <- !is.null(rv$icp)  && nrow(rv$icp)  > 0
    if (!has_gc && !has_hplc && !has_icp) {
      rv$log_step5 <- "⚠ No data available. Please complete Steps 1-3 (and optionally Step 4 for ICP)."
      return()
    }
    rv$log_step5 <- ""
    tryCatch({
      append_log("log_step5", "── Merging label information ──")
      label_path <- if (!is.null(input$label_file)) input$label_file$datapath
                    else "./Raw data/Sample Labelling_all data_GC+HPLC+ICP_30Sept2025.xlsx"
      if (!file.exists(label_path)) stop("Label file not found: ", label_path)
      sampinfo <- readxl::read_excel(label_path) %>%
        dplyr::select(c(File, Plastic_type, Subcategory, Polymer))
      rv$sampinfo <- sampinfo

      gc <- hplc <- icp <- NULL
      if (has_gc) {
        gc <- left_join(rv$gc, sampinfo, by = "File")
        rv$gc_labeled <- gc
        append_log("log_step5", sprintf("  GC:   %d rows  (%d plastic types)",
                                        nrow(gc), n_distinct(gc$Plastic_type, na.rm = TRUE)))
      } else append_log("log_step5", "  GC: No data (skipped)")

      if (has_hplc) {
        hplc <- left_join(rv$hplc, sampinfo, by = "File")
        rv$hplc_labeled <- hplc
        append_log("log_step5", sprintf("  HPLC: %d rows  (%d plastic types)",
                                        nrow(hplc), n_distinct(hplc$Plastic_type, na.rm = TRUE)))
      } else append_log("log_step5", "  HPLC: No data (skipped)")

      if (has_icp) {
        icp <- left_join(rv$icp, sampinfo, by = "File")
        rv$icp_labeled <- icp
        append_log("log_step5", sprintf("  ICP:  %d rows  (%d plastic types)",
                                        nrow(icp), n_distinct(icp$Plastic_type, na.rm = TRUE)))
      } else append_log("log_step5", "  ICP: No data (skipped)")

      combinations_created <- c()
      if (has_gc && has_hplc) {
        rv$gc_hplc <- bind_rows(gc, hplc); combinations_created <- c(combinations_created,"gc_hplc")
      }
      if (has_gc && has_icp) {
        rv$gc_icp <- bind_rows(icp, gc %>% dplyr::select(-any_of(c("m.z","RT"))))
        combinations_created <- c(combinations_created,"gc_icp")
      }
      if (has_hplc && has_icp) {
        rv$hplc_icp <- bind_rows(icp, hplc %>% dplyr::select(-any_of(c("m.z","RT"))))
        combinations_created <- c(combinations_created,"hplc_icp")
      }
      if (has_gc && has_hplc && has_icp) {
        rv$gc_hplc_icp <- bind_rows(icp, bind_rows(gc,hplc) %>% dplyr::select(-any_of(c("m.z","RT"))))
        combinations_created <- c(combinations_created,"gc_hplc_icp")
      }

      if (length(combinations_created))
        append_log("log_step5", sprintf("  Combinations created: %s",
                                        paste(combinations_created, collapse = ", ")))
      else
        append_log("log_step5", "  No combinations created (need at least 2 data types)")
      append_log("log_step5", "✓ Step 5 complete.")
    },
    error = function(e) append_log("log_step5", paste("✗ ERROR:", conditionMessage(e))))
  })

  output$step5_log <- renderText(rv$log_step5)
  output$step5_summary_table <- renderDT({
    datasets <- c(); rows <- c(); features <- c()
    addrow <- function(label, x) {
      if (!is.null(x) && nrow(x) > 0) {
        datasets <<- c(datasets, label)
        rows     <<- c(rows, nrow(x))
        features <<- c(features, n_distinct(x$Feature))
      }
    }
    addrow("gc",        rv$gc_labeled)
    addrow("hplc",      rv$hplc_labeled)
    addrow("icp",       rv$icp_labeled)
    addrow("gc_hplc",   rv$gc_hplc)
    addrow("gc_icp",    rv$gc_icp)
    addrow("hplc_icp",  rv$hplc_icp)
    addrow("gc_hplc_icp", rv$gc_hplc_icp)
    if (length(datasets) == 0)
      return(datatable(data.frame(Message = "No data available. Run Step 5 first."),
                       options = list(dom = 't')))
    datatable(data.frame(Dataset = datasets, Rows = rows, Features = features),
              options = list(dom = 't'))
  })
  output$step5_gc_table <- renderDT({
    if (is.null(rv$gc_labeled) || nrow(rv$gc_labeled) == 0)
      return(datatable(data.frame(Message = "No GC data available."), options = list(dom = 't')))
    datatable(head(rv$gc_labeled, 200), options = list(scrollX = TRUE, pageLength = 10))
  })

  # ===========================================================================
  #  STEP 6  —  Feature Reduction & Data Fusion
  # ===========================================================================
  get_labeled_dataset <- function(name) {
    switch(name,
      gc   = rv$gc_labeled,
      hplc = rv$hplc_labeled,
      icp  = rv$icp_labeled,
      NULL)
  }

  observeEvent(input$run_step6_filter, {
    rv$log_step6 <- ""
    tryCatch({
      ds_name <- input$feat_dataset
      ds <- get_labeled_dataset(ds_name)
      if (is.null(ds) || nrow(ds) == 0) {
        rv$log_step6 <- sprintf("⚠ Dataset '%s' is not available. Please complete Step 5 first.", ds_name)
        return()
      }
      mode <- input$feat_mode
      append_log("log_step6", sprintf("── Feature reduction: %s [%s] ──", ds_name, mode))
      results <- shared_feature_filtering(ds, feat_reduc_mode = mode)
      rv$filter_results <- results
      stats_text <- print_filtering_stats(results)
      rv$step6_stats_text <- stats_text
      append_log("log_step6", stats_text)
      append_log("log_step6", "✓ Feature reduction complete.")
    },
    error = function(e) append_log("log_step6", paste("✗ ERROR:", conditionMessage(e))))
  })

  observeEvent(input$run_step6_fusion, {
    has_gc   <- !is.null(rv$gc_labeled)   && nrow(rv$gc_labeled)   > 0
    has_hplc <- !is.null(rv$hplc_labeled) && nrow(rv$hplc_labeled) > 0
    has_icp  <- !is.null(rv$icp_labeled)  && nrow(rv$icp_labeled)  > 0
    if (!has_gc && !has_hplc && !has_icp) {
      append_log("log_step6", "\n⚠ No labeled data available. Please complete Step 5 first.")
      return()
    }
    tryCatch({
      append_log("log_step6", "\n── Running in-house filter + data fusion ──")
      meta <- c("Plastic_type", "technique", "Source", "Polymer", "Subcategory")
      fusion_datasets <- c(); fusion_samples <- c(); fusion_features <- c()

      if (has_gc) {
        rv$gc_clean <- shared_feature_filtering(rv$gc_labeled,   feat_reduc_mode="inhouse")$filtered_data
        fusion_datasets <- c(fusion_datasets, "gc")
        fusion_samples  <- c(fusion_samples, nrow(rv$gc_clean))
        fusion_features <- c(fusion_features, length(setdiff(names(rv$gc_clean), meta)))
        append_log("log_step6", sprintf("  GC: %d samples, %d features",
                                        nrow(rv$gc_clean), length(setdiff(names(rv$gc_clean), meta))))
      }
      if (has_hplc) {
        rv$hplc_clean <- shared_feature_filtering(rv$hplc_labeled, feat_reduc_mode="inhouse")$filtered_data
        fusion_datasets <- c(fusion_datasets, "hplc")
        fusion_samples  <- c(fusion_samples, nrow(rv$hplc_clean))
        fusion_features <- c(fusion_features, length(setdiff(names(rv$hplc_clean), meta)))
        append_log("log_step6", sprintf("  HPLC: %d samples, %d features",
                                        nrow(rv$hplc_clean), length(setdiff(names(rv$hplc_clean), meta))))
      }
      if (has_icp) {
        rv$icp_clean <- shared_feature_filtering(rv$icp_labeled,  feat_reduc_mode="inhouse")$filtered_data
        fusion_datasets <- c(fusion_datasets, "icp")
        fusion_samples  <- c(fusion_samples, nrow(rv$icp_clean))
        fusion_features <- c(fusion_features, length(setdiff(names(rv$icp_clean), meta)))
        append_log("log_step6", sprintf("  ICP: %d samples, %d features",
                                        nrow(rv$icp_clean), length(setdiff(names(rv$icp_clean), meta))))
      }

      if (has_gc && has_hplc && !is.null(rv$gc_hplc)) {
        rv$gc_hplc_clean <- data_fusion(rv$gc_hplc, source_feats=list(rv$gc_clean, rv$hplc_clean))
        fusion_datasets <- c(fusion_datasets, "gc_hplc")
        fusion_samples  <- c(fusion_samples, nrow(rv$gc_hplc_clean))
        fusion_features <- c(fusion_features, length(setdiff(names(rv$gc_hplc_clean), meta)))
      }
      if (has_gc && has_icp && !is.null(rv$gc_icp)) {
        rv$gc_icp_clean <- data_fusion(rv$gc_icp, source_feats=list(rv$gc_clean, rv$icp_clean))
        fusion_datasets <- c(fusion_datasets, "gc_icp")
        fusion_samples  <- c(fusion_samples, nrow(rv$gc_icp_clean))
        fusion_features <- c(fusion_features, length(setdiff(names(rv$gc_icp_clean), meta)))
      }
      if (has_hplc && has_icp && !is.null(rv$hplc_icp)) {
        rv$hplc_icp_clean <- data_fusion(rv$hplc_icp, source_feats=list(rv$hplc_clean, rv$icp_clean))
        fusion_datasets <- c(fusion_datasets, "hplc_icp")
        fusion_samples  <- c(fusion_samples, nrow(rv$hplc_icp_clean))
        fusion_features <- c(fusion_features, length(setdiff(names(rv$hplc_icp_clean), meta)))
      }
      if (has_gc && has_hplc && has_icp && !is.null(rv$gc_hplc_icp)) {
        rv$gc_hplc_icp_clean <- data_fusion(rv$gc_hplc_icp,
                                            source_feats=list(rv$gc_clean, rv$hplc_clean, rv$icp_clean))
        fusion_datasets <- c(fusion_datasets, "gc_hplc_icp")
        fusion_samples  <- c(fusion_samples, nrow(rv$gc_hplc_icp_clean))
        fusion_features <- c(fusion_features, length(setdiff(names(rv$gc_hplc_icp_clean), meta)))
      }

      rv$fusion_summary <- data.frame(
        Dataset = fusion_datasets, Samples = fusion_samples, Features = fusion_features)
      append_log("log_step6", "✓ Data fusion complete. See 'Fusion Summary' tab.")
    },
    error = function(e) append_log("log_step6", paste("✗ ERROR:", conditionMessage(e))))
  })

  output$step6_log   <- renderText(rv$log_step6)
  output$step6_stats <- renderText(rv$step6_stats_text)
  output$step6_table <- renderDT({
    if (is.null(rv$filter_results))
      return(datatable(data.frame(Message="No filter results. Run feature reduction first."),
                       options = list(dom = 't')))
    datatable(head(rv$filter_results$filtered_data, 200),
              options = list(scrollX = TRUE, pageLength = 10))
  })
  output$step6_fusion_table <- renderDT({
    if (is.null(rv$fusion_summary))
      return(datatable(data.frame(Message="No fusion summary. Run data fusion first."),
                       options = list(dom = 't')))
    datatable(rv$fusion_summary, options = list(dom = 't'))
  })

  # ===========================================================================
  #  STEP 7  —  Random Forest
  # ===========================================================================

  # Helper: pull cleaned dataset by name
  get_clean_dataset <- function(name) {
    switch(name,
      gc_clean          = rv$gc_clean,
      hplc_clean        = rv$hplc_clean,
      icp_clean         = rv$icp_clean,
      gc_hplc_clean     = rv$gc_hplc_clean,
      gc_icp_clean      = rv$gc_icp_clean,
      hplc_icp_clean    = rv$hplc_icp_clean,
      gc_hplc_icp_clean = rv$gc_hplc_icp_clean,
      NULL)
  }

  observeEvent(input$run_step7, {
    chosen <- input$rf_datasets
    if (length(chosen) == 0) {
      rv$log_step7 <- "⚠ Please select at least one dataset."
      return()
    }
    # Validate availability
    avail <- chosen[sapply(chosen, function(n) {
      d <- get_clean_dataset(n); !is.null(d) && nrow(d) > 0
    })]
    if (length(avail) == 0) {
      rv$log_step7 <- paste("⚠ None of the selected datasets are available.",
                            "Run Step 6 (Data Fusion) first.")
      return()
    }

    # Resolve split flags
    use_sb_only <- input$rf_split_mode == "sb_only"
    use_sb_env  <- input$rf_split_mode == "sb_env"

    # Parse advanced settings
    ntree_cand <- .parse_csv_int(input$rf_ntree,      c(100, 500, 1000, 2500))
    top_counts <- .parse_csv_int(input$rf_top_counts, c(100, 50, 25, 10))
    rv$rf_top_counts <- top_counts

    rv$log_step7 <- ""
    append_log("log_step7", sprintf("── Running RF for: %s ──", paste(avail, collapse=", ")))
    append_log("log_step7", sprintf("   split=%s | pairwise=%s | rfa=%s | top-n=%s | imp_screen=%s",
        input$rf_split_mode,
        input$rf_do_pairwise, input$rf_do_rfa,
        input$rf_do_topn, input$rf_do_imp_screen))

    # Parallel setup
    n_cores <- max(1, as.integer(input$rf_cores))
    cl <- NULL
    if (n_cores > 1) {
      cl <- tryCatch(parallel::makePSOCKcluster(n_cores), error = function(e) NULL)
      if (!is.null(cl)) doParallel::registerDoParallel(cl)
    }
    on.exit({
      if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE)
      foreach::registerDoSEQ()
    }, add = TRUE)

    rf_results_local <- list()
    n_total <- length(avail)

    # Determine which remove_cols to pass based on data type
    determine_remove_cols <- function(dset) {
      base <- c("Source", "Polymer", "technique")
      if ("Subcategory" %in% names(dset)) base <- c(base, "Subcategory")
      base
    }

    withProgress(message = "Running Random Forest", value = 0, {
      for (i in seq_along(avail)) {
        ds_name <- avail[i]
        dset    <- get_clean_dataset(ds_name)
        setProgress(value = (i - 1) / n_total,
                    detail = sprintf("Dataset %d/%d: %s", i, n_total, ds_name))
        append_log("log_step7", sprintf("\n>> Dataset: %s (%d rows × %d cols)",
                                        ds_name, nrow(dset), ncol(dset)))

        # Coerce Plastic_type to factor up front
        dset$Plastic_type <- as.character(dset$Plastic_type)
        # Drop missing-class rows
        dset <- dset[!is.na(dset$Plastic_type) & dset$Plastic_type != "", ]

        rm_cols <- determine_remove_cols(dset)
        res <- tryCatch({
          run_rf_analysis_manuscript1(
            data                          = dset,
            type_col                      = "Plastic_type",
            remove_cols                   = rm_cols,
            train_proportion              = 0.8,
            use_source_split              = use_sb_env,
            use_store_vs_environmental_split = use_sb_only,
            data_name                     = ds_name,
            group_for_significance        = "Plastic_type",
            do_pairwise_test              = isTRUE(input$rf_do_pairwise),
            do_top_importance_selection   = isTRUE(input$rf_do_topn),
            top_importance_counts         = top_counts,
            do_rfa                        = isTRUE(input$rf_do_rfa),
            do_impute_norm_screen         = isTRUE(input$rf_do_imp_screen),
            seed                          = as.integer(input$rf_seed),
            min_sample_number             = as.integer(input$rf_min_sample),
            K_outer_max                   = as.integer(input$rf_K_outer),
            ntree_candidates              = ntree_cand,
            metric                        = input$rf_metric,
            ntree_screen                  = as.integer(input$rf_ntree_screen))
        }, error = function(e) {
          append_log("log_step7", sprintf("   ✗ ERROR on %s: %s", ds_name, conditionMessage(e)))
          NULL
        })
        if (!is.null(res)) {
          rf_results_local[[ds_name]] <- res
          append_log("log_step7", sprintf("   ✓ Completed %s (%d folds)", ds_name, length(res)))
        }
      }
    })
    rv$rf_results <- rf_results_local
    # Persist which methods were actually enabled in this run (used by the
    # Metrics / Selected-Features / Confusion-Matrix views to decide whether a
    # method was "disabled" or "ran-but-no-features").
    rv$rf_methods_enabled <- list(
      do_pairwise   = isTRUE(input$rf_do_pairwise),
      do_rfa        = isTRUE(input$rf_do_rfa),
      do_topn       = isTRUE(input$rf_do_topn),
      do_imp_screen = isTRUE(input$rf_do_imp_screen)
    )

    # ── Run-time Feature-Selection Summary block in Console Log ──────────────
    append_log("log_step7", "\n── Feature-Selection Summary ──")
    for (ds_name in names(rf_results_local)) {
      res <- rf_results_local[[ds_name]]
      fold_names <- names(res)
      counts <- rf_fold_feature_counts(res, top_counts, rv$rf_methods_enabled)
      append_log("log_step7", sprintf("\n  Dataset: %s  (%d folds: %s)",
                                      ds_name, length(fold_names),
                                      paste(fold_names, collapse = ", ")))
      report <- function(label, vec, enabled) {
        if (!enabled) {
          append_log("log_step7", sprintf("    %-32s : disabled (not run)", label))
          return()
        }
        line <- sprintf("    %-32s : %s", label,
                        paste(sprintf("%s=%d", fold_names, vec), collapse = ", "))
        append_log("log_step7", line)
        zero_folds <- fold_names[vec == 0]
        if (length(zero_folds) > 0) {
          append_log("log_step7",
            sprintf("      ⚠ No features in %s — no MCC / confusion matrix / selected features for these folds.",
                    paste(zero_folds, collapse = ", ")))
        }
      }
      report("All features",                counts$all, TRUE)
      report("Pairwise significance test",  counts$sig %||% rep(0, length(fold_names)),
             rv$rf_methods_enabled$do_pairwise)
      report("Recursive Feature Addition",  counts$rfa %||% rep(0, length(fold_names)),
             rv$rf_methods_enabled$do_rfa)
      if (rv$rf_methods_enabled$do_topn) {
        for (tn in top_counts) {
          report(sprintf("Top %d by importance", tn),
                 counts[[paste0("top_", tn)]] %||% rep(0, length(fold_names)), TRUE)
        }
      } else {
        append_log("log_step7", "    Top-N by importance              : disabled (not run)")
      }
    }
    append_log("log_step7", "\n✓ Step 7 finished.")

    # Update view selectors
    nm <- names(rf_results_local)
    updateSelectInput(session, "rf_view_dataset", choices = nm, selected = nm[1])
    if (length(nm)) {
      folds <- names(rf_results_local[[nm[1]]])
      updateSelectInput(session, "rf_view_fold", choices = folds, selected = folds[1])
    }
    method_choices <- c("All features" = "all_feats")
    if (isTRUE(input$rf_do_pairwise)) method_choices <- c(method_choices, "Pairwise" = "sig")
    if (isTRUE(input$rf_do_rfa))      method_choices <- c(method_choices, "RFA"      = "rfa")
    if (isTRUE(input$rf_do_topn)) {
      for (tn in top_counts) method_choices[sprintf("Top %d", tn)] <- sprintf("top_%d", tn)
    }
    updateSelectInput(session, "rf_view_method", choices = method_choices)
    # Also populate Step 8 selectors
    updateSelectInput(session, "hca_dataset", choices = nm, selected = nm[1])
    if (length(nm)) {
      folds <- names(rf_results_local[[nm[1]]])
      updateSelectInput(session, "hca_fold", choices = folds, selected = folds[1])
    }
    # Also populate Step 9 selectors
    updateSelectInput(session, "pcaumap_dataset", choices = nm, selected = nm[1])
    if (length(nm)) {
      folds <- names(rf_results_local[[nm[1]]])
      updateSelectInput(session, "pcaumap_fold", choices = folds, selected = folds[1])
    }
    updateSelectInput(session, "pcaumap_feature_set", choices = method_choices)
  })

  # When user changes RF dataset selection, update fold choices
  observeEvent(input$rf_view_dataset, {
    nm <- input$rf_view_dataset
    if (is.null(nm) || !nzchar(nm) || is.null(rv$rf_results[[nm]])) return()
    folds <- names(rv$rf_results[[nm]])
    updateSelectInput(session, "rf_view_fold", choices = folds, selected = folds[1])
  })
  observeEvent(input$hca_dataset, {
    nm <- input$hca_dataset
    if (is.null(nm) || !nzchar(nm) || is.null(rv$rf_results[[nm]])) return()
    folds <- names(rv$rf_results[[nm]])
    updateSelectInput(session, "hca_fold", choices = folds, selected = folds[1])
  })
  observeEvent(input$pcaumap_dataset, {
    nm <- input$pcaumap_dataset
    if (is.null(nm) || !nzchar(nm) || is.null(rv$rf_results[[nm]])) return()
    folds <- names(rv$rf_results[[nm]])
    updateSelectInput(session, "pcaumap_fold", choices = folds, selected = folds[1])
  })

  # Render confusion matrix
  # -- Confusion-Matrix display: info-box when empty, else plot -----------------
  # Resolve method status for the currently-selected (dataset, fold, method).
  rf_method_status <- function() {
    ds <- input$rf_view_dataset
    fd <- input$rf_view_fold
    md <- input$rf_view_method
    if (is.null(ds) || is.null(fd) || is.null(md) || !nzchar(md)) return(NULL)
    if (is.null(rv$rf_results[[ds]]) || is.null(rv$rf_results[[ds]][[fd]])) return(NULL)
    fobj <- rv$rf_results[[ds]][[fd]]
    opts <- rv$rf_methods_enabled

    enabled <- switch(md,
      all_feats = TRUE,
      sig       = isTRUE(opts$do_pairwise),
      rfa       = isTRUE(opts$do_rfa),
      isTRUE(opts$do_topn) && grepl("^top_", md))

    feats_n <- switch(md,
      all_feats = length(setdiff(names(fobj$final_imp_norm_train %||% list()),
                                 c("Plastic_type","Source","Polymer",
                                   "technique","Subcategory"))),
      sig       = length(fobj$sig_selected_feats %||% character()),
      rfa       = length(fobj$rfa_selected_feats %||% character()),
      length(fobj$top_importance_results[[md]]$selected_features %||% character()))

    plt <- switch(md,
      all_feats = fobj$conf_mat_all_feats,
      sig       = fobj$conf_mat_sig,
      rfa       = fobj$conf_mat_rfa,
      fobj$top_importance_results[[md]]$conf_mat_top)

    list(ds = ds, fd = fd, md = md,
         label = rf_method_label(md),
         enabled = enabled, feats_n = feats_n, plot = plt)
  }

  # Top-level UI dispatcher for the Confusion Matrices tab
  output$step7_confmat_display <- renderUI({
    st <- rf_method_status()
    if (is.null(st)) return(div(class = "alert alert-info",
                                "Run Step 7 first."))

    if (!st$enabled) {
      return(div(class = "alert alert-secondary", style = "margin: 20px; padding: 20px;",
        h4(icon("info-circle"), sprintf(" %s — method not run", st$label)),
        p("This feature-selection method was disabled in the Run options.",
          " There is no MCC score, no confusion matrix, and no selected-feature list",
          " for it. Re-enable the checkbox in the input panel and re-run Step 7",
          " to populate this view."),
        p(tags$small(sprintf("Dataset: %s  |  Fold: %s", st$ds, st$fd)))
      ))
    }
    if (st$feats_n == 0 || is.null(st$plot)) {
      return(div(class = "alert alert-warning", style = "margin: 20px; padding: 20px;",
        h4(icon("exclamation-triangle"),
           sprintf(" %s — no features selected for this fold", st$label)),
        p(sprintf("The %s ran on %s / %s but did not select any features.",
                  st$label, st$ds, st$fd)),
        p("Because no features were selected, there is ",
          tags$b("no MCC score"), ", ",
          tags$b("no confusion matrix"), ", and ",
          tags$b("no selected-feature list"),
          " for this method on this fold."),
        p(tags$small("Other folds may still have results — switch the 'Fold' selector",
                     " to check. The Metrics tab shows the per-method success rate",
                     " across all folds."))
      ))
    }
    # Otherwise show the plot
    plotOutput("step7_confmat_plot", height = "600px")
  })

  output$step7_confmat_plot <- renderPlot({
    st <- rf_method_status()
    req(st, !is.null(st$plot))
    print(st$plot)
  })

  # Append a Console-Log note when the user picks an empty/disabled method.
  # Debounced via a small ID so we don't spam identical lines.
  observeEvent(list(input$rf_view_dataset, input$rf_view_fold, input$rf_view_method), {
    st <- rf_method_status()
    if (is.null(st)) return()
    key <- paste(st$ds, st$fd, st$md, sep = "/")
    if (!is.null(rv$last_view_log_key) && identical(rv$last_view_log_key, key)) return()

    if (!st$enabled) {
      append_log("log_step7", sprintf(
        "[view] %s / %s — %s is disabled (not run). No MCC / confusion matrix / selected features.",
        st$ds, st$fd, st$label))
      rv$last_view_log_key <- key
    } else if (st$feats_n == 0) {
      append_log("log_step7", sprintf(
        "[view] %s / %s — %s selected 0 features for this fold. No MCC / confusion matrix / selected features.",
        st$ds, st$fd, st$label))
      rv$last_view_log_key <- key
    }
  }, ignoreInit = TRUE)

  # Render aggregated metrics_df (with empty-state Notes column)
  output$step7_metrics_df <- renderDT({
    ds <- input$rf_view_dataset
    if (is.null(ds) || is.null(rv$rf_results[[ds]]))
      return(datatable(data.frame(Message = "Run Step 7 first.")))
    tc <- rv$rf_top_counts %||% c(100, 50, 25, 10)
    m  <- aggregate_rf_metrics(rv$rf_results[[ds]], tc, rv$rf_methods_enabled)
    datatable(m, options = list(dom = 't', scrollX = TRUE),
              rownames = FALSE)
  })

  # Selected features per method — explicit rows for empty/disabled methods
  output$step7_features_table <- renderDT({
    ds <- input$rf_view_dataset; fd <- input$rf_view_fold
    if (is.null(ds) || is.null(fd) || is.null(rv$rf_results[[ds]][[fd]]))
      return(datatable(data.frame(Message = "Run Step 7 first.")))
    fobj <- rv$rf_results[[ds]][[fd]]
    opts <- rv$rf_methods_enabled
    rows <- list()

    add_method_row <- function(method, feats, enabled) {
      rows[[method]] <<- if (!enabled) {
        data.frame(Method = method, N = "—",
                   Features = "Method not run (disabled by user)",
                   stringsAsFactors = FALSE)
      } else if (length(feats) > 0) {
        data.frame(Method = method, N = as.character(length(feats)),
                   Features = paste(feats, collapse = ", "),
                   stringsAsFactors = FALSE)
      } else {
        data.frame(Method = method, N = "0",
                   Features = "(no features selected — not computed for this fold)",
                   stringsAsFactors = FALSE)
      }
    }

    add_method_row("Pairwise", fobj$sig_selected_feats, isTRUE(opts$do_pairwise))
    add_method_row("RFA",      fobj$rfa_selected_feats, isTRUE(opts$do_rfa))
    if (isTRUE(opts$do_topn)) {
      tc <- rv$rf_top_counts %||% c(100, 50, 25, 10)
      for (tn in tc) {
        key <- paste0("top_", tn)
        feats <- fobj$top_importance_results[[key]]$selected_features
        add_method_row(sprintf("Top %d", tn), feats, TRUE)
      }
    } else {
      rows[["Top-N"]] <- data.frame(Method = "Top-N", N = "—",
                                    Features = "Method not run (disabled by user)",
                                    stringsAsFactors = FALSE)
    }
    if (!length(rows))
      return(datatable(data.frame(Message = "No feature lists in this fold.")))
    datatable(do.call(rbind, rows),
              options = list(scrollX = TRUE, pageLength = 10), rownames = FALSE)
  })

  # Best imputation / normalization combo per fold — explicit empty state
  output$step7_impnorm_table <- renderDT({
    ds <- input$rf_view_dataset
    if (is.null(ds) || is.null(rv$rf_results[[ds]]))
      return(datatable(data.frame(Message = "Run Step 7 first.")))
    if (!isTRUE(rv$rf_methods_enabled$do_imp_screen))
      return(datatable(data.frame(
        Note = paste("Impute+Normalize screening was disabled in the Run options.",
                     "No best combination was determined; the default (kNN + zscore",
                     "or whatever the helper sets as fallback) was used silently.")),
        options = list(dom = 't'), rownames = FALSE))
    folds <- rv$rf_results[[ds]]
    rows <- do.call(rbind, lapply(seq_along(folds), function(i) {
      f <- folds[[i]]
      imp  <- f$best_imputation    %||% NA
      norm <- f$best_normalization %||% NA
      data.frame(
        Fold          = names(folds)[i],
        Imputation    = if (is.na(imp))  "—" else as.character(imp),
        Normalization = if (is.na(norm)) "—" else as.character(norm),
        stringsAsFactors = FALSE)
    }))
    datatable(rows, options = list(dom = 't'), rownames = FALSE)
  })

  output$step7_log <- renderText(rv$log_step7)

  # Export Results button
  observeEvent(input$rf_export_btn, {
    ds <- input$rf_view_dataset; fd <- input$rf_view_fold
    if (is.null(ds) || is.null(fd) || is.null(rv$rf_results[[ds]][[fd]])) {
      append_log("log_step7", "\n⚠ No results to export.")
      return()
    }
    outdir <- trimws(input$rf_export_dir)
    if (!nzchar(outdir)) outdir <- "results_rf"
    prefix <- trimws(input$rf_export_prefix)
    if (!nzchar(prefix)) prefix <- "rf"
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
    out <- tryCatch({
      export_rf_results(rv$rf_results[[ds]], output_dir = outdir, fold = fd,
                        prefix = paste0(prefix, "_", ds))
    }, error = function(e) {
      append_log("log_step7", paste("✗ Export ERROR:", conditionMessage(e)))
      NULL
    })
    if (!is.null(out)) {
      append_log("log_step7", sprintf("\n✓ Exported %s / %s -> %s", ds, fd, outdir))
      rv$rf_export_dir <- outdir
    }
  })

  # ===========================================================================
  #  STEP 8  —  Hierarchical Cluster Analysis
  # ===========================================================================
  observeEvent(input$run_step8, {
    ds <- input$hca_dataset
    fd <- input$hca_fold
    if (is.null(ds) || is.null(fd) || is.null(rv$rf_results[[ds]]) ||
        is.null(rv$rf_results[[ds]][[fd]])) {
      rv$log_step8 <- "⚠ No Step 7 result selected. Run Step 7 first."
      return()
    }
    rv$log_step8 <- ""
    tryCatch({
      fold_obj <- rv$rf_results[[ds]][[fd]]
      df <- dplyr::bind_rows(fold_obj$final_imp_norm_train,
                             fold_obj$final_imp_norm_test) %>%
        dplyr::mutate(Plastic_type = gsub("..", " ", Plastic_type, fixed = TRUE)) %>%
        dplyr::mutate(Plastic_type = gsub(".",  " ", Plastic_type, fixed = TRUE)) %>%
        dplyr::mutate(Plastic_type_Source = paste0(Plastic_type, "-", Source))

      hc_df <- df %>%
        dplyr::select(-any_of(c("Subcategory","Source","Plastic_type",
                                "Polymer","technique"))) %>%
        dplyr::select(where(is.numeric)) %>%
        as.data.frame()
      hc_df[!is.finite(as.matrix(hc_df))] <- 0
      keep <- rowSums(hc_df, na.rm = TRUE) > 0
      hc_df_clean <- hc_df[keep, , drop = FALSE]
      labels_clean <- df$Plastic_type_Source[keep]

      dist_method <- input$hca_distance
      link_method <- input$hca_linkage
      append_log("log_step8", sprintf("── HCA: %s / %s (dist=%s, link=%s) ──",
                                      ds, fd, dist_method, link_method))
      append_log("log_step8", sprintf("   %d samples, %d features",
                                      nrow(hc_df_clean), ncol(hc_df_clean)))

      use_vegan <- dist_method %in% c("bray", "robust.aitchison")
      if (use_vegan && any(hc_df_clean < 0)) {
        append_log("log_step8",
          "   Data contains negative values - falling back to Euclidean.")
        d <- dist(hc_df_clean, method = "euclidean")
      } else if (use_vegan) {
        d <- vegan::vegdist(hc_df_clean, method = dist_method)
      } else {
        d <- dist(hc_df_clean, method = dist_method)
      }
      if (any(!is.finite(d))) stop("Distance matrix contains non-finite values.")

      hca_samp <- hclust(d, method = link_method)
      ccc <- suppressWarnings(cor(d, stats::cophenetic(hca_samp)))
      append_log("log_step8", sprintf("   CCC = %.4f", ccc))

      rv$hca_object <- hca_samp
      rv$hca_labels <- labels_clean
      rv$hca_dist   <- d
      rv$hca_ccc    <- ccc
      append_log("log_step8", "✓ HCA complete.")
    },
    error = function(e) append_log("log_step8", paste("✗ ERROR:", conditionMessage(e))))
  })

  output$step8_log <- renderText(rv$log_step8)

  output$step8_ccc <- renderText({
    if (is.na(rv$hca_ccc)) "Run HCA to compute the cophenetic correlation coefficient (CCC)."
    else sprintf("Cophenetic Correlation Coefficient (CCC): %.4f", rv$hca_ccc)
  })

  output$step8_dendro <- renderPlot({
    if (is.null(rv$hca_object)) {
      plot.new(); title(main = "Click 'Run HCA' to generate the dendrogram")
      return()
    }
    par(cex = 1)
    plot(rv$hca_object, labels = rv$hca_labels, hang = -1,
         main = sprintf("HCA Dendrogram  |  CCC = %.3f", rv$hca_ccc))
  })

  # Download dendrogram as PNG
  output$hca_download <- downloadHandler(
    filename = function() {
      paste0("PlastiPrint_HCA_", input$hca_dataset, "_", input$hca_fold,
             "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      if (is.null(rv$hca_object)) {
        png(file, width = 1200, height = 600)
        plot.new(); title(main = "No HCA result")
        dev.off()
        return()
      }
      png(file, width = 1400, height = 700, res = 120)
      par(cex = 1)
      plot(rv$hca_object, labels = rv$hca_labels, hang = -1,
           main = sprintf("HCA Dendrogram  |  CCC = %.3f", rv$hca_ccc))
      dev.off()
    }
  )

  # Final "Done / Export All" button — exports all RF results for all datasets
  observeEvent(input$export_all_btn, {
    if (length(rv$rf_results) == 0) {
      append_log("log_step8", "\n⚠ Nothing to export yet — run Step 7 first.")
      return()
    }
    outdir <- trimws(input$rf_export_dir)
    if (!nzchar(outdir)) outdir <- "results_rf"
    prefix <- trimws(input$rf_export_prefix)
    if (!nzchar(prefix)) prefix <- "rf"
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
    append_log("log_step8", sprintf("\n── Exporting all RF results to %s ──", outdir))
    for (ds in names(rv$rf_results)) {
      folds <- names(rv$rf_results[[ds]])
      for (fd in folds) {
        out <- tryCatch(
          export_rf_results(rv$rf_results[[ds]], output_dir = outdir, fold = fd,
                            prefix = paste0(prefix, "_", ds, "_", fd)),
          error = function(e) {
            append_log("log_step8", sprintf("   ✗ %s/%s: %s", ds, fd, conditionMessage(e))); NULL
          })
        if (!is.null(out))
          append_log("log_step8", sprintf("   ✓ %s/%s exported", ds, fd))
      }
    }
    append_log("log_step8", "✓ Export All complete.")
  })

  # ===========================================================================
  #  STEP 9  —  PCA & UMAP
  # ===========================================================================
  observeEvent(input$run_step9, {
    ds <- input$pcaumap_dataset
    fd <- input$pcaumap_fold
    if (is.null(ds) || is.null(fd) || is.null(rv$rf_results[[ds]]) ||
        is.null(rv$rf_results[[ds]][[fd]])) {
      rv$log_step9 <- "⚠ No Step 7 result selected. Run Step 7 first."
      return()
    }
    rv$log_step9 <- ""
    tryCatch({
      fold_obj <- rv$rf_results[[ds]][[fd]]
      feat_set <- input$pcaumap_feature_set

      df <- dplyr::bind_rows(fold_obj$final_imp_norm_train,
                             fold_obj$final_imp_norm_test) %>%
        dplyr::mutate(Plastic_type = gsub("..", " ", Plastic_type, fixed = TRUE)) %>%
        dplyr::mutate(Plastic_type = gsub(".",  " ", Plastic_type, fixed = TRUE))

      feats <- switch(feat_set,
        all_feats = setdiff(names(df), c("Plastic_type","Source","Polymer",
                                          "technique","Subcategory")),
        sig       = fold_obj$sig_selected_feats,
        rfa       = fold_obj$rfa_selected_feats,
        fold_obj$top_importance_results[[feat_set]]$selected_features)
      feats <- intersect(feats %||% character(), names(df))
      if (length(feats) < 2) stop("Fewer than 2 features available for this feature set.")

      hc_df <- df %>% dplyr::select(all_of(feats)) %>% as.data.frame()
      hc_df[!is.finite(as.matrix(hc_df))] <- 0
      keep <- rowSums(hc_df, na.rm = TRUE) > 0
      hc_df_clean <- hc_df[keep, , drop = FALSE]
      col_var <- apply(hc_df_clean, 2, var, na.rm = TRUE)
      hc_df_clean <- hc_df_clean[, col_var > 0, drop = FALSE]
      if (ncol(hc_df_clean) < 2) stop("Fewer than 2 non-constant features remain.")
      if (nrow(hc_df_clean) < 4) stop("Fewer than 4 samples remain after cleaning.")

      Type <- clean_plastic_label(df$Plastic_type[keep])
      Source <- df$Source[keep]
      top_types <- names(sort(table(Type), decreasing = TRUE))[1:min(9, length(unique(Type)))]
      TypePlot <- ifelse(Type %in% top_types, Type, "Other")

      append_log("log_step9", sprintf("── PCA/UMAP: %s / %s (features=%s) ──", ds, fd, feat_set))
      append_log("log_step9", sprintf("   %d samples, %d features", nrow(hc_df_clean), ncol(hc_df_clean)))

      pca <- prcomp(hc_df_clean, center = TRUE, scale. = TRUE)
      ve <- (pca$sdev^2 / sum(pca$sdev^2))[1:2] * 100
      append_log("log_step9", sprintf("   PCA: PC1+PC2 = %.1f%% + %.1f%% = %.1f%%", ve[1], ve[2], sum(ve)))

      n_neighbors <- min(15, nrow(hc_df_clean) - 1)
      if (n_neighbors < 2) stop("Not enough samples for UMAP (need >= 3).")
      set.seed(123)
      umap_cfg <- umap.defaults
      umap_cfg$n_neighbors <- n_neighbors
      umap_cfg$random_state <- 123
      umap_fit <- umap(scale(hc_df_clean), config = umap_cfg)

      k <- max(2, min(round(input$pcaumap_hclust_k), nrow(hc_df_clean) - 1))
      dist_label <- "Robust Aitchison"
      hca_dist <- tryCatch({
        if (any(hc_df_clean <= 0)) stop("non-positive values")
        vegan::vegdist(hc_df_clean, method = "robust.aitchison")
      }, error = function(e) {
        dist_label <<- "Euclidean (fallback)"
        dist(scale(hc_df_clean), method = "euclidean")
      })
      hclust_obj <- hclust(hca_dist, method = "average")
      hclusters <- cutree(hclust_obj, k = k)
      append_log("log_step9", sprintf("   HCA overlay: %s distance, average linkage, k=%d clusters", dist_label, k))

      plot_meta <- data.frame(Type = Type, TypePlot = TypePlot, Source = Source,
                              HCluster = factor(hclusters))
      pca_df  <- cbind(plot_meta, PC1 = pca$x[, 1], PC2 = pca$x[, 2])
      umap_df <- cbind(plot_meta, UMAP1 = umap_fit$layout[, 1], UMAP2 = umap_fit$layout[, 2])

      p_pca <- ggplot(pca_df, aes(PC1, PC2, color = TypePlot, shape = Source)) +
        geom_point(size = 2.5, alpha = 0.85) +
        scale_color_viridis_d() +
        scale_shape_manual(values = source_shapes) +
        labs(title = sprintf("PCA — %s / %s (%s)", ds, fd, feat_set),
             x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]),
             color = "Type") +
        theme_minimal(base_size = 13)

      p_umap <- ggplot(umap_df, aes(UMAP1, UMAP2))
      if (isTRUE(input$pcaumap_show_hulls)) {
        p_umap <- p_umap +
          ggforce::geom_mark_hull(aes(group = HCluster), fill = "grey50", alpha = 0.08,
                                   color = "grey40", linetype = 2,
                                   expand = unit(3, "mm"), radius = unit(2, "mm"),
                                   concavity = 2.5)
      }
      p_umap <- p_umap +
        geom_point(aes(color = TypePlot, shape = Source), size = 2.5, alpha = 0.9) +
        scale_color_viridis_d() +
        scale_shape_manual(values = source_shapes) +
        labs(title = sprintf("UMAP — %s / %s (%s)%s", ds, fd, feat_set,
                              if (isTRUE(input$pcaumap_show_hulls))
                                sprintf(" | HCA k=%d hulls", k) else ""),
             color = "Type") +
        theme_minimal(base_size = 13)

      rv$step9_pca_plot <- p_pca
      rv$step9_umap_plot <- p_umap
      append_log("log_step9", "✓ PCA/UMAP complete.")
    },
    error = function(e) append_log("log_step9", paste("✗ ERROR:", conditionMessage(e))))
  })

  output$step9_log <- renderText(rv$log_step9)

  output$step9_pca_plot <- renderPlot({
    if (is.null(rv$step9_pca_plot)) {
      plot.new(); title(main = "Click 'Run PCA / UMAP' to generate plots")
      return()
    }
    rv$step9_pca_plot
  })

  output$step9_umap_plot <- renderPlot({
    if (is.null(rv$step9_umap_plot)) {
      plot.new(); title(main = "Click 'Run PCA / UMAP' to generate plots")
      return()
    }
    rv$step9_umap_plot
  })

  output$pcaumap_download_pca <- downloadHandler(
    filename = function() {
      paste0("PlastiPrint_PCA_", input$pcaumap_dataset, "_", input$pcaumap_fold,
             "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      if (is.null(rv$step9_pca_plot)) {
        png(file, width = 1200, height = 700); plot.new(); title(main = "No PCA result"); dev.off()
        return()
      }
      ggsave(file, rv$step9_pca_plot, width = 8.5, height = 5.5, dpi = 150)
    }
  )

  output$pcaumap_download_umap <- downloadHandler(
    filename = function() {
      paste0("PlastiPrint_UMAP_", input$pcaumap_dataset, "_", input$pcaumap_fold,
             "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      if (is.null(rv$step9_umap_plot)) {
        png(file, width = 1200, height = 700); plot.new(); title(main = "No UMAP result"); dev.off()
        return()
      }
      ggsave(file, rv$step9_umap_plot, width = 9.5, height = 5.5, dpi = 150)
    }
  )

} # end server

# =============================================================================
#  Launch
# =============================================================================
shinyApp(ui = ui, server = server)
