# SPDX-License-Identifier: AGPL-3.0-only
#
# Figure: source-tracking performance of the three single-platform fingerprints
# (ATD-GC-MS, HPLC-QToF-MS, ICP-MS/MS) under the two Step 7 train/test modes of
# "Complete workflow_13May2026.Rmd":
#   Option 1  SB-only train -> ENV test   (use_store_vs_environmental_split = TRUE)
#   Option 2  SB+ENV train  -> ENV test   (use_source_split = TRUE)
#
# Run from the project root:  Rscript figure_source_tracking_performance.R
#
# DATA SOURCE
#   The RF result objects (rf_results_SB_train / rf_results_SB_ENV_train) are not
#   saved anywhere in the repository, so MCC values are read from the published
#   Table S8 in the Supporting Spreadsheets (mean and SD of MCC_Multiclass across
#   outer folds, as summarised by the "Table 2" chunk of
#   "Scripts for Tables and Figures (maintext and SI).Rmd").
#   Multi-instrument combinations are deliberately excluded.
#
# WHY THERE ARE NO "correct items / items tested" LABELS
#   1. The result objects are not available, so there are no predictions to count.
#   2. Even if they were, run_rf_analysis_manuscript1() replaces File names with
#      "<Plastic_type>_<technique>_repN" row names, so predictions cannot be
#      traced back to a USE-xx item without changing the helper function.
#   3. Outer folds are caret::createFolds() over replicate/time-point rows, not
#      leave-one-item-out. In Option 2 the other replicates of a held-out
#      environmental item remain in the training set.
#   Instead, the figure shows a structural count taken from the sample labels:
#   how many environmental items belong to a class that is present in the
#   store-bought training set (the only items Option 1 can ever get right).

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})

# ---- Paths ------------------------------------------------------------------
supp_xlsx  <- "Supplementary documents for publication/Computational fingerprinting paper_Supporting Spreadsheets_8March2026_revisedsubmission.xlsx"
raw_dir    <- "Raw Data - for Github only"
label_xlsx <- file.path(raw_dir, "Sample Labelling_all data_GC+HPLC+ICP_26Feb2026.xlsx")
fig_dir    <- "figures"
res_dir    <- "results_rf"
dir.create(fig_dir, showWarnings = FALSE)
dir.create(res_dir, showWarnings = FALSE)

# Feature-selection method reported in the manuscript (best MCC in every row of Table S8)
fs_method <- "Recursive Feature Addition"

tech_map <- c(GC = "ATD-GC-MS", HPLC = "HPLC-QToF-MS", ICP = "ICP-MS/MS")
scen_map <- c(
  "SB only"    = "Universal fingerprint",
  "SB and ENV" = "Suspect–source comparison"
)
result_obj <- c(
  "SB only"    = "rf_results_SB_train (Complete workflow Rmd, Step 7 Option 1)",
  "SB and ENV" = "rf_results_SB_ENV_train (Complete workflow Rmd, Step 7 Option 2)"
)

# ---- 1. Published MCC (Table S8) ---------------------------------------------
s8 <- read_excel(supp_xlsx, sheet = "Table S8", skip = 2) %>%
  dplyr::select(1:5) %>%
  setNames(c("dataset", "scenario_raw", "method", "mcc", "mcc_sd")) %>%
  dplyr::filter(dataset %in% names(tech_map))

s8_all_feats <- s8 %>%
  dplyr::filter(method == "All features") %>%
  dplyr::select(dataset, scenario_raw, mcc_all_features = mcc, mcc_all_features_sd = mcc_sd)

mcc_df <- s8 %>%
  dplyr::filter(method == fs_method) %>%
  dplyr::left_join(s8_all_feats, by = c("dataset", "scenario_raw"))

# ---- 2. Environmental items per platform (from raw data + labels) -----------
# Same exclusions as the ATD-GC-MS import in Step 1.1 of the workflow
gc_excluded <- "USE[-_](01|02|05|06|09|10|11|12|13|15|16|17|18|19|20)([^0-9]|$)|USSB[-_](01|08)([^0-9]|$)"
item_id <- function(x) str_replace(str_extract(x, "US(SB|E)[-_][0-9]+"), "_", "-")

gc_files   <- list.files(file.path(raw_dir, "ATDGCMS"), pattern = "\\.csv$")
gc_files   <- gc_files[!str_detect(gc_files, "^2022") & !str_detect(gc_files, gc_excluded)]
hplc_files <- list.files(file.path(raw_dir, "HPLCTOFMS"), pattern = "\\.xls$", recursive = TRUE)
icp_files  <- list.files(file.path(raw_dir, "ICPMS_Trace metal"), pattern = "rawdata", full.names = TRUE)
icp_ids    <- unlist(lapply(icp_files, function(f) read_excel(f)$File))

raw_items <- bind_rows(
  tibble(dataset = "GC",   item = item_id(gc_files)),
  tibble(dataset = "HPLC", item = item_id(hplc_files)),
  tibble(dataset = "ICP",  item = item_id(icp_ids))
) %>%
  dplyr::filter(!is.na(item)) %>%
  dplyr::distinct()

labels <- read_excel(label_xlsx) %>%
  dplyr::transmute(dataset = technique, item = item_id(File), Plastic_type) %>%
  dplyr::filter(!is.na(item)) %>%
  dplyr::distinct()

items <- raw_items %>%
  dplyr::inner_join(labels, by = c("dataset", "item")) %>%
  dplyr::mutate(source = if_else(str_detect(item, "USE"), "ENV", "SB"))

item_summary <- items %>%
  dplyr::group_by(dataset) %>%
  dplyr::summarise(
    total_items = n_distinct(item[source == "ENV"]),
    items_with_training_class = n_distinct(
      item[source == "ENV" & Plastic_type %in% Plastic_type[source == "SB"]]
    ),
    .groups = "drop"
  )

# ---- 3. Tidy audit table ------------------------------------------------------
perf_plot_df <- mcc_df %>%
  dplyr::left_join(item_summary, by = "dataset") %>%
  dplyr::transmute(
    technique = factor(tech_map[dataset], levels = rev(tech_map)),
    scenario = factor(scen_map[scenario_raw], levels = scen_map),
    training_config = scenario_raw,
    mcc,
    mcc_sd,
    # Not recoverable from the available outputs; see header
    correct_items = NA_integer_,
    total_items,
    items_with_training_class,
    proportion_correct = NA_real_,
    feature_selection_method = method,
    mcc_all_features,
    mcc_all_features_sd,
    result_object = result_obj[scenario_raw],
    value_source = "Supporting Spreadsheets, Table S8 (mean/SD over outer folds)",
    cv_scheme = if_else(
      scenario_raw == "SB only",
      "k-fold over SB rows; every fold tests on the same full ENV set",
      "stratified k-fold over ENV replicate rows (not grouped by item)"
    )
  ) %>%
  dplyr::arrange(technique, scenario)

print(as.data.frame(perf_plot_df %>% dplyr::select(-result_object, -value_source, -cv_scheme)))
readr::write_csv(perf_plot_df, file.path(res_dir, "source_tracking_figure_data.csv"))

# ---- 4. QA --------------------------------------------------------------------
stopifnot(nrow(perf_plot_df) == 6)
stopifnot(all(perf_plot_df$mcc >= -1 & perf_plot_df$mcc <= 1))
stopifnot(all(perf_plot_df$items_with_training_class <= perf_plot_df$total_items))
stopifnot(all(is.na(perf_plot_df$correct_items) |
                perf_plot_df$correct_items <= perf_plot_df$total_items))
stopifnot(!any(str_detect(as.character(perf_plot_df$technique), "-.*-.*-|\\+")))
# Published Table S8 values for the RFA rows (validation only, not used for plotting)
expected <- tibble(
  technique = tech_map[c("GC", "GC", "HPLC", "HPLC", "ICP", "ICP")],
  training_config = rep(c("SB only", "SB and ENV"), 3),
  expected_mcc = c(0.37, 0.51, 0.13, 0.89, 0.08, 0.90)
)
chk <- perf_plot_df %>%
  dplyr::mutate(technique = as.character(technique)) %>%
  dplyr::inner_join(expected, by = c("technique", "training_config"))
stopifnot(nrow(chk) == 6, isTRUE(all.equal(chk$mcc, chk$expected_mcc)))

# ---- 5. Figure ----------------------------------------------------------------
fig_title    <- "Source-tracking performance of three chemical fingerprints"
fig_subtitle <- "Weathered environmental plastics classified by Random Forest (RFA-selected features)"
fig_caption  <- paste0(
  "Points: mean multiclass MCC across outer folds (bars: ± SD), from Supporting Spreadsheets Table S8.\n",
  "Right column: environmental items whose plastic category also occurs among the store-bought training items;\n",
  "items in other categories cannot be classified correctly without environmental training data.\n",
  "Suspect–source folds split replicate and time-point rows, so other replicates of a tested item stay in training."
)

col_univ <- "#6B7280"
col_susp <- "#1F4E79"
shape_vals <- c(21, 19)
names(shape_vals) <- scen_map
col_vals <- c(col_univ, col_susp)
names(col_vals) <- scen_map

seg_df <- perf_plot_df %>%
  dplyr::select(technique, scenario, mcc) %>%
  tidyr::pivot_wider(names_from = scenario, values_from = mcc) %>%
  setNames(c("technique", "x_univ", "x_susp"))

lab_df <- perf_plot_df %>%
  dplyr::mutate(
    label = sprintf("%.2f ± %.2f", mcc, mcc_sd),
    # Universal label below the point, suspect-source label above, so close pairs never collide
    y_nudge = if_else(scenario == scen_map[["SB only"]], -0.3, 0.3)
  )

item_lab_df <- perf_plot_df %>%
  dplyr::distinct(technique, items_with_training_class, total_items) %>%
  dplyr::mutate(label = sprintf("%d / %d", items_with_training_class, total_items))

p_main <- ggplot(perf_plot_df, aes(y = technique)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.7, colour = "grey25") +
  geom_segment(
    data = seg_df,
    aes(x = x_univ, xend = x_susp, y = technique, yend = technique),
    linewidth = 0.8, colour = "grey70"
  ) +
  geom_errorbar(
    aes(xmin = mcc - mcc_sd, xmax = mcc + mcc_sd, colour = scenario),
    orientation = "y", width = 0.12, linewidth = 0.6
  ) +
  geom_point(
    aes(x = mcc, shape = scenario, colour = scenario),
    size = 5.5, stroke = 1.6, fill = "white"
  ) +
  geom_label(
    data = lab_df,
    aes(x = mcc, y = as.numeric(technique) + y_nudge, label = label, colour = scenario),
    size = 4.6, show.legend = FALSE, fill = "white", border.colour = NA,
    label.padding = unit(0.1, "lines")
  ) +
  geom_text(
    data = item_lab_df,
    aes(x = 1.27, y = technique, label = label),
    size = 5.2, colour = "grey20", hjust = 0.5
  ) +
  annotate(
    "text", x = 1.27, y = 3.55, label = "ENV items with a\nstore-bought class",
    size = 3.8, colour = "grey30", lineheight = 0.9, vjust = 0
  ) +
  scale_x_continuous(
    name = "Matthews correlation coefficient (MCC)",
    limits = c(-1, 1.42),
    breaks = c(-1, -0.5, 0, 0.5, 1),
    labels = c("−1", "−0.5", "0", "0.5", "1"),
    expand = expansion(mult = 0)
  ) +
  scale_y_discrete(name = NULL, expand = expansion(add = c(0.5, 0.9))) +
  scale_shape_manual(
    values = shape_vals, name = NULL,
    labels = c(
      "Universal fingerprint\n(store-bought training only)",
      "Suspect–source comparison\n(store-bought + environmental training)"
    )
  ) +
  scale_colour_manual(values = col_vals, guide = "none") +
  guides(shape = guide_legend(override.aes = list(colour = col_vals, size = 5))) +
  coord_cartesian(clip = "off") +
  labs(title = fig_title, subtitle = fig_subtitle) +
  theme_minimal(base_size = 16) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(colour = "grey90"),
    axis.text.y = element_text(size = 17, face = "bold", colour = "grey10"),
    axis.text.x = element_text(size = 13),
    axis.title.x = element_text(size = 14, margin = margin(t = 8)),
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = 12.5, lineheight = 0.95),
    legend.key.width = unit(1.2, "lines"),
    legend.box.spacing = unit(2, "pt"),
    plot.title = element_text(face = "bold", size = 20),
    plot.subtitle = element_text(size = 13.5, colour = "grey30", margin = margin(b = 4)),
    plot.title.position = "plot"
  )

strip_df <- tibble(
  x = 1:3,
  head = c("ATD-GC-MS", "HPLC-QToF-MS", "ICP-MS/MS"),
  body = c("volatile / semivolatile\norganic fingerprint",
           "polar / nonvolatile\norganic fingerprint",
           "elemental\nfingerprint")
)

p_strip <- ggplot(strip_df) +
  annotate("text", x = 2, y = 1.75, label = "Complementary lines of chemical evidence",
           fontface = "bold", size = 4.6, colour = "grey20") +
  geom_text(aes(x = x, y = 1.2, label = head), fontface = "bold", size = 4.2, colour = col_susp) +
  geom_text(aes(x = x, y = 0.55, label = body), size = 3.7, colour = "grey30", lineheight = 0.9) +
  scale_x_continuous(limits = c(0.4, 3.6)) +
  scale_y_continuous(limits = c(0.05, 2.0)) +
  theme_void() +
  theme(plot.background = element_rect(fill = "grey97", colour = NA))

fig <- (p_main / p_strip) +
  plot_layout(heights = c(4.2, 1.3)) +
  plot_annotation(
    caption = fig_caption,
    theme = theme(
      plot.caption = element_text(size = 9.5, colour = "grey35", hjust = 0, lineheight = 1.05),
      plot.caption.position = "plot",
      plot.background = element_rect(fill = "white", colour = NA)
    )
  )

w <- 10; h <- 6.6
ggsave(file.path(fig_dir, "source_tracking_performance.png"), fig, width = w, height = h, dpi = 300, bg = "white")
cat("Saved figure to", fig_dir, "and audit table to", file.path(res_dir, "source_tracking_figure_data.csv"), "\n")
