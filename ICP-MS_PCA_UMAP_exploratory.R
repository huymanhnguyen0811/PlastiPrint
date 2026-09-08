## Exploratory PCA + UMAP on trace-metal (ICP-MS) chemical fingerprints
##
## No pre-tidied sample x feature matrix is checked into this repo (the
## GC/HPLC non-targeted pipeline caches to a gitignored data_cache/ and is
## never committed). The ICP-MS trace-metal data, however, is already a
## rectangular sample x element matrix per raw file, so this script
## reproduces the exact Step 4 "Importing Trace metal data" logic from
## `Complete workflow_13May2026.Rmd` (H2/H2-HMI mode de-duplication, best
## recovery/uncertainty method per element, seawater-metal removal) to
## build a tidy matrix, then runs PCA and UMAP on it.
##
## Run from the repository root.

library(readxl)
library(dplyr)
library(stringr)
library(tibble)
library(tidyr)
library(ggplot2)
library(umap)
library(viridis)
library(randomForest)
library(vegan)
library(dendextend)
library(ggforce)

raw_dir <- "Raw Data - for Github only/ICPMS_Trace metal"
out_dir <- "results_pca_umap"
dir.create(out_dir, showWarnings = FALSE)

## ---- Step 4 (reproduced): import + clean trace metal data ----

read_trace_metal <- function(path) {
  if (!file.exists(path)) stop(paste("File missing:", path))
  readxl::read_excel(path) %>% column_to_rownames(var = "File")
}

icpms_round1 <- read_trace_metal(file.path(raw_dir, "icpms_round1_rawdata_removal USE-01-rep1 and USE-03 (only2observations).xlsx"))
icpms_round2_batch1 <- read_trace_metal(file.path(raw_dir, "icpms_round2_batch1_rawdata_26Feb2026.xlsx"))
icpms_round2_batch2 <- read_trace_metal(file.path(raw_dir, "icpms_round2_batch2_rawdata.xlsx"))

remove_H2_HMI_modus <- function(icpms_df) {
  strings <- colnames(icpms_df)
  pattern <- "^(.*\\[ )([A-Za-z0-9 ]+)( \\])$"
  keep_cols <- strings[!str_detect(str_match(strings, pattern)[, 3], "H2 HMI")]
  icpms_df %>% dplyr::select(all_of(keep_cols))
}

icpms_round1 <- remove_H2_HMI_modus(icpms_round1)
icpms_round2_batch1 <- remove_H2_HMI_modus(icpms_round2_batch1)
icpms_round2_batch2 <- remove_H2_HMI_modus(icpms_round2_batch2)

rec_unc_1 <- readxl::read_excel(file.path(raw_dir, "Trace metal data_Recovery_Uncertainty_round1.xlsx"))
rec_unc_2_1 <- readxl::read_excel(file.path(raw_dir, "Trace metal data_Recovery_Uncertainty_round2_batch1.xlsx"))
rec_unc_2_2 <- readxl::read_excel(file.path(raw_dir, "Trace metal data_Recovery_Uncertainty_round2_batch2.xlsx"))

get_final_metal_names <- function(recovery_uncertainty, target_cols) {
  recovery_uncertainty <- recovery_uncertainty %>% filter(`Element+Modus` %in% target_cols)
  first_numbers <- sub(" .*", "", recovery_uncertainty$`Element+Modus`)
  kept_names <- c()
  for (num in unique(first_numbers)) {
    rows <- recovery_uncertainty[first_numbers == num, ]
    if (nrow(rows) == 1) {
      kept_names <- c(kept_names, rows$`Element+Modus`)
    } else {
      best_row <- rows[which.max(rows$`Recovery [%]`), ]
      if (nrow(best_row) == 0) best_row <- rows[which.min(rows$`U(k=2)`), ]
      kept_names <- c(kept_names, best_row$`Element+Modus`)
    }
  }
  return(kept_names)
}

final_names_1 <- get_final_metal_names(rec_unc_1, colnames(icpms_round1))
final_names_2_1 <- get_final_metal_names(rec_unc_2_1, colnames(icpms_round2_batch1))
final_names_2_2 <- get_final_metal_names(rec_unc_2_2, colnames(icpms_round2_batch2))

final_icpms_round1 <- icpms_round1 %>% dplyr::select(all_of(final_names_1))
final_icpms_round2_batch1 <- icpms_round2_batch1 %>% dplyr::select(all_of(final_names_2_1))
final_icpms_round2_batch2 <- icpms_round2_batch2 %>% dplyr::select(all_of(final_names_2_2))

common_cols <- Reduce(intersect, list(names(final_icpms_round1), names(final_icpms_round2_batch1), names(final_icpms_round2_batch2)))
combined_icpms <- rbind(final_icpms_round1[common_cols], final_icpms_round2_batch1[common_cols], final_icpms_round2_batch2[common_cols])

strings_to_remove <- c("Na", "Mg", "K", "Ca")
cols_to_remove <- unique(unlist(sapply(strings_to_remove, function(x) grep(x, colnames(combined_icpms)))))
if (length(cols_to_remove) > 0) combined_icpms <- combined_icpms[, -cols_to_remove]

combined_icpms[] <- lapply(combined_icpms, function(x) as.numeric(as.character(x)))

cat("Tidy ICP-MS matrix:", nrow(combined_icpms), "samples x", ncol(combined_icpms), "elements\n")

## ---- Attach plastic-type / source labels ----

labels <- read_excel(
  "Raw Data - for Github only/Sample Labelling_all data_GC+HPLC+ICP_26Feb2026.xlsx",
  sheet = "Grouping1"
) %>%
  filter(technique == "ICP") %>%
  distinct(File, .keep_all = TRUE) %>%
  select(File, Plastic_type, Category, Polymer, Source = Type)

meta <- tibble(File = rownames(combined_icpms)) %>%
  left_join(labels, by = "File") %>%
  mutate(
    Source = ifelse(is.na(Source), ifelse(str_detect(File, "USE"), "Environmental", "Store-Bought"), Source),
    Category = ifelse(is.na(Category), "Unlabeled", Category)
  )

## ---- Prepare numeric matrix: drop zero-variance cols, impute, log-transform, scale ----

mat <- as.matrix(combined_icpms)

zero_var <- apply(mat, 2, function(x) var(x, na.rm = TRUE) == 0 || all(is.na(x)))
mat <- mat[, !zero_var, drop = FALSE]

col_medians <- apply(mat, 2, median, na.rm = TRUE)
for (j in seq_len(ncol(mat))) {
  mat[is.na(mat[, j]), j] <- col_medians[j]
}

## ICP-MS intensities span orders of magnitude and include occasional
## non-positive values (blank-subtracted signal) -> shift-and-log10.
shift_and_log <- function(x) {
  min_pos <- min(x[x > 0], na.rm = TRUE)
  x[x <= 0] <- min_pos / 2
  log10(x)
}
mat_log <- apply(mat, 2, shift_and_log)
rownames(mat_log) <- rownames(mat)

stopifnot(identical(rownames(mat_log), meta$File))

## ---- PCA ----

pca <- prcomp(mat_log, center = TRUE, scale. = TRUE)
var_explained <- (pca$sdev^2 / sum(pca$sdev^2))[1:2] * 100

pca_df <- meta %>%
  mutate(PC1 = pca$x[, 1], PC2 = pca$x[, 2])

## ---- UMAP ----

set.seed(123)
umap_cfg <- umap.defaults
umap_cfg$n_neighbors <- min(15, nrow(mat_log) - 1)
umap_cfg$random_state <- 123
umap_fit <- umap(scale(mat_log), config = umap_cfg)

umap_df <- meta %>%
  mutate(UMAP1 = umap_fit$layout[, 1], UMAP2 = umap_fit$layout[, 2])

## ---- Save tidy outputs ----

combined_out <- pca_df %>%
  left_join(umap_df %>% select(File, UMAP1, UMAP2), by = "File")

write.csv(combined_out, file.path(out_dir, "icpms_pca_umap_coordinates.csv"), row.names = FALSE)
write.csv(as.data.frame(mat_log) %>% rownames_to_column("File"),
          file.path(out_dir, "icpms_tidy_log_matrix.csv"), row.names = FALSE)

## ---- Clean "Type" labels: drop the redundant standalone word "plastic" ----

clean_type <- function(x) {
  x <- gsub("(?i)\\bplastic\\b", "", x, perl = TRUE)
  x <- trimws(gsub("\\s+", " ", x))
  ifelse(x == "", "Plastics (misc.)", x)
}

meta <- meta %>% mutate(Type = clean_type(Plastic_type))
pca_df <- pca_df %>% mutate(Type = clean_type(Plastic_type))
umap_df <- umap_df %>% mutate(Type = clean_type(Plastic_type))

top_categories <- names(sort(table(pca_df$Type), decreasing = TRUE))[1:9]
pca_df <- pca_df %>% mutate(TypePlot = ifelse(Type %in% top_categories, Type, "Other"))
umap_df <- umap_df %>% mutate(TypePlot = ifelse(Type %in% top_categories, Type, "Other"))

## ---- Random-Forest feature importance -> top-10 elements ----

set.seed(123)
rf <- randomForest(x = as.data.frame(mat_log), y = as.factor(meta$Plastic_type),
                    ntree = 1000, importance = TRUE)
oob_acc <- mean(rf$predicted == meta$Plastic_type)
cat(sprintf("RF OOB accuracy classifying Type from all %d elements: %.1f%%\n", ncol(mat_log), 100 * oob_acc))

top_feats <- names(sort(importance(rf, type = 2)[, 1], decreasing = TRUE))[1:10]
mat_top <- mat_log[, top_feats, drop = FALSE]

pca_top <- prcomp(mat_top, center = TRUE, scale. = TRUE)
ve_top <- (pca_top$sdev^2 / sum(pca_top$sdev^2))[1:2] * 100
pca_top_df <- meta %>% mutate(PC1 = pca_top$x[, 1], PC2 = pca_top$x[, 2],
                               TypePlot = ifelse(Type %in% top_categories, Type, "Other"))

set.seed(123)
umap_cfg_top <- umap.defaults
umap_cfg_top$n_neighbors <- min(15, nrow(mat_top) - 1)
umap_cfg_top$random_state <- 123
umap_top_fit <- umap(scale(mat_top), config = umap_cfg_top)
umap_top_df <- meta %>% mutate(UMAP1 = umap_top_fit$layout[, 1], UMAP2 = umap_top_fit$layout[, 2],
                                TypePlot = ifelse(Type %in% top_categories, Type, "Other"))

## ---- Plots (viridis) ----

p_pca_source <- ggplot(pca_df, aes(PC1, PC2, color = Source)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_viridis_d(end = 0.8) +
  labs(
    title = "PCA of ICP-MS trace-metal fingerprints",
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2])
  ) +
  theme_minimal(base_size = 13)

source_shapes <- c("Environmental" = 17, "Store-Bought" = 16)

p_pca_type <- ggplot(pca_df, aes(PC1, PC2, color = TypePlot, shape = Source)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_viridis_d() +
  scale_shape_manual(values = source_shapes) +
  labs(
    title = "PCA of ICP-MS trace-metal fingerprints by type",
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2]),
    color = "Type"
  ) +
  theme_minimal(base_size = 13)

p_umap_source <- ggplot(umap_df, aes(UMAP1, UMAP2, color = Source)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_viridis_d(end = 0.8) +
  labs(title = "UMAP of ICP-MS trace-metal fingerprints") +
  theme_minimal(base_size = 13)

p_umap_type <- ggplot(umap_df, aes(UMAP1, UMAP2, color = TypePlot, shape = Source)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_viridis_d() +
  scale_shape_manual(values = source_shapes) +
  labs(title = "UMAP of ICP-MS trace-metal fingerprints by type", color = "Type") +
  theme_minimal(base_size = 13)

p_pca_top <- ggplot(pca_top_df, aes(PC1, PC2, color = TypePlot, shape = Source)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_viridis_d() +
  scale_shape_manual(values = source_shapes) +
  labs(
    title = "PCA (top-10 RF-important elements) by type",
    x = sprintf("PC1 (%.1f%%)", ve_top[1]),
    y = sprintf("PC2 (%.1f%%)", ve_top[2]),
    color = "Type"
  ) +
  theme_minimal(base_size = 13)

p_umap_top <- ggplot(umap_top_df, aes(UMAP1, UMAP2, color = TypePlot, shape = Source)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_color_viridis_d() +
  scale_shape_manual(values = source_shapes) +
  labs(title = "UMAP (top-10 RF-important elements) by type", color = "Type") +
  theme_minimal(base_size = 13)

ggsave(file.path(out_dir, "pca_by_source.png"), p_pca_source, width = 7, height = 5, dpi = 150)
ggsave(file.path(out_dir, "pca_by_type.png"), p_pca_type, width = 8.5, height = 5, dpi = 150)
ggsave(file.path(out_dir, "umap_by_source.png"), p_umap_source, width = 7, height = 5, dpi = 150)
ggsave(file.path(out_dir, "umap_by_type.png"), p_umap_type, width = 8.5, height = 5, dpi = 150)
ggsave(file.path(out_dir, "pca_top10feats_by_type.png"), p_pca_top, width = 8.5, height = 5, dpi = 150)
ggsave(file.path(out_dir, "umap_top10feats_by_type.png"), p_umap_top, width = 8.5, height = 5, dpi = 150)

## ---- Hierarchical Cluster Analysis (HCA) ----
## Reproduces the app's Step 8 logic (PlastiPrint_app-UI, ~line 2220-2240):
## Robust Aitchison distance via vegan::vegdist, falling back to Euclidean
## when negative values are present; average linkage. Run on the same
## top-10 RF-selected elements used for the feature-selected PCA/UMAP above,
## since that mirrors Step 8 operating on Step 7's feature-selected data.

mat_top_raw <- mat[, top_feats, drop = FALSE]  # pre-log, non-negative-checked
hca_labels <- paste0(meta$Type, "-", meta$Source)

if (any(mat_top_raw <= 0)) {
  cat("HCA: non-positive values present -> falling back to Euclidean distance on scaled log data.\n")
  hca_dist <- dist(scale(mat_top), method = "euclidean")
} else {
  cat("HCA: using Robust Aitchison distance (vegan::vegdist).\n")
  hca_dist <- vegan::vegdist(mat_top_raw, method = "robust.aitchison")
}

hca <- hclust(hca_dist, method = "average")
ccc <- suppressWarnings(cor(hca_dist, cophenetic(hca)))
cat(sprintf("HCA cophenetic correlation coefficient (CCC): %.3f\n", ccc))

dend <- as.dendrogram(hca)
labels(dend) <- hca_labels[order.dendrogram(dend)]

source_by_leaf <- meta$Source[order.dendrogram(dend)]
viridis_pal <- viridis(2, end = 0.8)
leaf_cols <- ifelse(source_by_leaf == "Environmental", viridis_pal[1], viridis_pal[2])
dend <- dend %>%
  set("labels_cex", 0.5) %>%
  set("labels_colors", leaf_cols) %>%
  set("branches_lwd", 0.6)

png(file.path(out_dir, "hca_dendrogram.png"), width = 14, height = 30, units = "in", res = 150)
par(mar = c(3, 3, 3, 20))
plot(dend, horiz = TRUE,
     main = sprintf("HCA of ICP-MS trace-metal fingerprints (top-10 elements)\nCCC = %.3f | color = Source", ccc))
legend("topleft", legend = c("Environmental", "Store-Bought"), text.col = viridis_pal, bty = "n", cex = 1.2)
dev.off()

hclusters <- cutree(hca, k = 9)
write.csv(data.frame(File = meta$File, Type = meta$Type, Source = meta$Source,
                      cluster_k9 = hclusters),
          file.path(out_dir, "hca_cluster_assignments.csv"), row.names = FALSE)

## ---- Overlay HCA cluster hulls on the UMAP (same top-10 feature space) ----

umap_top_df <- umap_top_df %>% mutate(HCluster = factor(hclusters))

p_umap_top_hca <- ggplot(umap_top_df, aes(UMAP1, UMAP2)) +
  ggforce::geom_mark_hull(aes(group = HCluster), fill = "grey50", alpha = 0.08,
                           color = "grey40", linetype = 2, expand = unit(3, "mm"),
                           radius = unit(2, "mm"), concavity = 2.5) +
  geom_point(aes(color = TypePlot, shape = Source), size = 2.5, alpha = 0.9) +
  scale_color_viridis_d() +
  scale_shape_manual(values = source_shapes) +
  labs(
    title = "UMAP (top-10 RF-important elements) with HCA clusters (k=9, dashed hulls)",
    color = "Type"
  ) +
  theme_minimal(base_size = 13)

ggsave(file.path(out_dir, "umap_top10feats_with_hclusters.png"), p_umap_top_hca, width = 9.5, height = 5.5, dpi = 150)

cat("\nDone. Outputs written to:", normalizePath(out_dir), "\n")
cat("PC1+PC2 variance explained (all features):", sprintf("%.1f%% + %.1f%% = %.1f%%\n", var_explained[1], var_explained[2], sum(var_explained)))
cat("PC1+PC2 variance explained (top-10 features):", sprintf("%.1f%% + %.1f%% = %.1f%%\n", ve_top[1], ve_top[2], sum(ve_top)))
