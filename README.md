# PlastiPrint — Microplastic Computational Fingerprinting Tool

PlastiPrint is an open-source R Shiny application for chemical fingerprinting of microplastics using machine learning. It wraps the complete computational workflow from the publication (**"Chemical Fingerprints of New vs. Weathered Microplastics: A Machine Learning Approach"**)[https://doi.org/10.1021/acs.est.6c03575] into an interactive GUI with configurable parameters at every step.

## Overview

The app provides a 9-tab guided workflow that takes users from raw instrumental data (ATD-GC-MS, HPLC-QToF-MS, ICP-MS/MS) through data processing, feature selection, ML classification (Random Forest, SVM, XGBoost), and hierarchical cluster analysis (HCA) — all without writing code.

## License

- **Code**: GNU AGPL v3.0 only (AGPL-3.0-only) — see `LICENSE.txt`
- **Data**: CC BY-SA 4.0 — see `Raw Data - for Github only/LICENSE-DATA.txt`

---

## 💻 System Specifications

Developed and tested on:
- **OS**: Windows 11 Pro 64-bit (Build 26200)
- **Processor**: 11th Gen Intel Core i7-11800H @ 2.30 GHz (16 CPUs)
- **Memory**: 16 GB RAM

---

## Prerequisites

- **R ≥ 4.2** — <https://cran.r-project.org>
- **RStudio** (recommended) — <https://posit.co/downloads/>
- **Git** (optional — only needed if you use `git clone`) — <https://git-scm.com/downloads>

---

## 📂 Repository Structure

| File / Folder | Description |
|---|---|
| `PlastiPrint_app-UI_13May2026_version2.R` | The Shiny application (main entry point) |
| `Helper Function using only RF_Github_13May2026.R` | Custom functions for compound grouping, imputation, and feature filtering |
| `Complete workflow_13May2026.Rmd` | The original R Markdown workflow (non-GUI version) |
| `Scripts for Tables and Figures (maintext and SI).Rmd` | Scripts for publication figures and tables |
| `Raw Data - for Github only/` | Raw instrumental data (ATD-GC-MS, HPLC-QToF-MS, ICP-MS/MS) |
| `Supplementary documents for publication/` | Supplementary materials accompanying the publication |
| `R package versions.txt` | Exact package versions for reproducibility |

---

## Installation

### 1. Get the code

You have two options.

#### Option A — Git clone (recommended if you have Git installed)

Open a terminal (PowerShell or Command Prompt on Windows, Terminal on macOS/Linux) and run:

```bash
git clone https://github.com/huymanhnguyen0811/PlastiPrint.git
cd PlastiPrint
```

This downloads the entire repository into a folder called `PlastiPrint/` in your current directory.

#### Option B — Download ZIP (no Git required)

1. Open the repository page in your browser: <https://github.com/huymanhnguyen0811/PlastiPrint>
2. Click the green **`<> Code`** button near the top right.
3. Select **Download ZIP**.
4. Extract the downloaded archive (`PlastiPrint-main.zip`) to a folder you can find easily — for example:
   - Windows: `C:\Users\<your-username>\Documents\PlastiPrint\`
   - macOS / Linux: `~/Documents/PlastiPrint/`
5. Rename the extracted folder from `PlastiPrint-main` to `PlastiPrint` if you'd like a cleaner path.

#### Verify the layout

After either option, the project root should contain (at minimum):

```
PlastiPrint/
├── PlastiPrint_app-UI_13May2026_version2.R
├── Helper Function using only RF_Github_13May2026.R
├── Complete workflow_13May2026.Rmd
├── Raw Data - for Github only/
└── README.md
```

The app file and the helper `.R` file must be at the **same directory level** for the app to find the helper on launch.

---

### 2. Install R packages

The first launch automatically installs any missing CRAN packages (5–10 min on initial run). To pre-install manually:

```r
install.packages(c(
  "shiny", "shinydashboard", "shinyjs", "shinyWidgets", "DT",
  "tidyverse", "readxl", "writexl", "data.table",
  "crayon", "grid", "gridExtra", "ggplot2", "pheatmap",
  "viridis", "ggpubr", "scales",
  "stats", "vegan", "cluster", "factoextra", "FactoMineR",
  "ggforce", "MASS", "moments",
  "caret", "randomForest", "ranger", "e1071", "xgboost",
  "missForest", "pROC", "VIM", "softImpute", "RANN",
  "foreach", "doParallel", "parallel"
))
```

<details>
<summary><strong>Full package list by category</strong></summary>

| Category | Packages |
|---|---|
| **GUI** | shiny, shinydashboard, shinyjs, shinyWidgets, DT |
| **Data** | tidyverse, readxl, writexl, data.table |
| **Visualization** | ggplot2, pheatmap, viridis, ggpubr, scales, gridExtra |
| **Statistics** | vegan, cluster, factoextra, FactoMineR, MASS, moments |
| **ML** | caret, randomForest, ranger, e1071, xgboost |
| **Missing data** | missForest, pROC, VIM, softImpute, RANN |
| **Parallel** | foreach, doParallel, parallel |

For exact versions used during development, see `R package versions.txt`.

</details>

---

### 3. Launch and how to use the PlastiPrint app

## Launching the App

**If you use RStudio:** Open `PlastiPrint_app-UI_13May2026_version2.R` → click **Run App** in the top-right of the editor pane.

**R console:**
```r
setwd("/path/to/PlastiPrint")
shiny::runApp("PlastiPrint_app-UI_13May2026_version2.R")
```

**Terminal:**
```bash
cd /path/to/PlastiPrint
Rscript -e 'shiny::runApp("PlastiPrint_app-UI_13May2026_version2.R", launch.browser = TRUE)'
```

---

## Using PlastiPrint — Step by Step

> 📸 **Screenshots referenced below live in an `images/` folder at the repo root.** If you don't see them rendered on GitHub, create that folder (`mkdir images` from the repo root) and add the PNGs with the filenames given under each step.

Work through the sidebar tabs sequentially:

### Welcome Tab
Set your **Project Root Directory** (the folder containing `Raw Data/`) and click **Set Directory**.

![Welcome tab — set project directory](images/welcome_set_dir.png)

---

### Step 1 — Import & Grouping

Run three sub-steps in order:
1. **Import Data** — reads raw CSV/XLS files from `Raw Data/ATDGCMS` and `Raw Data/HPLCTOFMS`.
2. **Group Compounds** — bins compounds by RT and m/z (adjustable thresholds).
3. **Remove Benchmarks** — strips internal-standard compounds.

![Step 1 — Import file selection modal](images/step1_import_modal.png)

![Step 1 — Grouping configuration on the input panel](images/step1_grouping.png)

![Step 1 — Console log after import + grouping + benchmark removal](images/step1_console.png)

---

### Step 2 — Blank Subtraction

Remove signals below `mean_blank + N × SD_blank` (default SD multiplier: 3×).

![Step 2 — Output console showing row counts after blank subtraction](images/step2_output.png)

---

### Step 3 — Source Assignment

Automatically labels samples as **Store-Bought** or **Environmental** based on filename patterns (files containing `USE` → Environmental, otherwise → Store-Bought).

![Step 3 — Output showing SB / ENV split counts](images/step3_output.png)

---

### Step 4 — Trace Metal Import

Reads three raw ICP-MS Excel files plus three matching recovery/uncertainty Excel files. Click **Import ICP-MS Files** to open the upload modal, then **Process Metals** to run the deduplication and seawater-metal removal. Configurable: which seawater metals to discard (default: Na, Mg, K, Ca).

![Step 4 — ICP-MS file upload modal (3 raw + 3 recovery files)](images/step4_icp_modal.png)

![Step 4 — Processed ICP-MS data preview](images/step4_preview.png)

---

### Step 5 — Label Merge

Merges plastic type, subcategory, and polymer labels from the sample-labelling file. Creates all dataset combinations (gc, hplc, icp, and their fusions).

![Step 5 — Label file selection -> run button -> Datasets summary table showing rows and feature counts per combination](images/step5_summary.png)

---

### Step 6 — Feature Reduction & Data Fusion

- **6.1–6.3 Feature Reduction**: pick a single-instrument dataset (`gc`, `hplc`, or `icp`) and choose a filtering method (regular 80% rule, modified 80% rule, or in-house comprehensive 4-rule filter).
- **6.4 Data Fusion**: applies in-house filtering per instrument, then fuses multi-instrument datasets (`gc_hplc`, `gc_icp`, `hplc_icp`, `gc_hplc_icp`).Selec

![Step 6 — Feature reduction options and filtering statistics and Fusion summary table](images/step6_summary.png)

---

### Step 7 — Random Forest Classification

PlastiPrint's machine-learning step runs a Random Forest classifier on each selected dataset and compares feature-selection strategies side-by-side. All folds run with parallel processing, with a progress spinner indicating the current dataset.

#### Configuration panel (left)

| Section | Options |
|---|---|
| **Datasets to analyze** | Checkbox group: `gc_clean`, `hplc_clean`, `icp_clean`, `gc_hplc_clean`, `gc_icp_clean`, `hplc_icp_clean`, `gc_hplc_icp_clean`. Only datasets cleaned in Step 6 will run; the rest are listed but skipped. |
| **Train/Test split mode** | Radio buttons (one of three): *Standard CV split*, *SB-only train, ENV test (Option 1)*, *SB+ENV train, ENV test (Option 2)*. Default is Option 2. |
| **Feature-Selection Methods** | Four independent toggles: Pairwise significance test, Recursive Feature Addition (RFA), Top-N by importance, and Impute+Normalize screening (requires `FactoMineR`). |
| **Main Hyperparameters** | Random seed (default 123), `min_sample_number` (default 1), `K_outer_max` cross-validation folds (default 5), number of cores (defaults to detected physical cores − 1). |
| **Advanced settings** | Collapsible block: `ntree_candidates` (comma-separated, default `100, 500, 1000, 2500`), `top_importance_counts` (comma-separated, default `100, 50, 25, 10`), selection metric (Accuracy/Kappa), `ntree_screen` (default 500). |

Click **Run Random Forest** to launch the pipeline.

![Step 7 — Full configuration panel](images/step7_config.png)

#### Viewing results

After the run finishes, three dropdowns under **View results** populate automatically:

- **Dataset** — pick which of your selected datasets to view
- **Fold (for confusion matrix)** — pick which cross-validation fold's results to display
- **Confusion-matrix method** — switch between "All features", "Pairwise", "RFA", and each enabled Top-N value

#### Output tabs (right)

| Tab | Content |
|---|---|
| **Console Log** | Run-time log including a **Feature-Selection Summary** block per dataset, showing per-fold feature counts for every method and warnings about folds that produced zero features. |
| **Confusion Matrices** | Renders the confusion matrix for the selected (Dataset, Fold, Method). If a method ran but selected zero features for that fold, an info-box replaces the plot and explains there is no MCC score, no confusion matrix, and no selected-feature list for that combination. |
| **Metrics (mean ± SD)** | Aggregated MCC values across folds with mean, SD, and a "Folds used / total" column. Empty methods show "—" with explanatory Notes. |
| **Selected Features** | The actual features chosen by each method for the currently-viewed fold. Empty methods get an explicit row reading "(no features selected — not computed for this fold)". |
| **Best Imputation+Normalization** | The best (imputation, normalization) pair per fold from the screening step. Shows "—" cleanly if screening was disabled. |

![Step 7 — View Random Forest results](images/step7_viewresults.png)

#### Exporting results

Set the **Output directory** (default `results_rf`) and **Prefix** (default `rf`), then click **Export Results** to write CSV (per-method feature lists, metrics) and PNG (confusion matrices) files for the currently-viewed dataset and fold.

---

### Step 8 — Hierarchical Cluster Analysis (HCA)

HCA operates on the post-processing (imputed + normalized) data from a Step 7 fold. Samples are labeled `Plastic_type-Source` (e.g., `PET-Store-Bought`).

#### Configuration panel (left)

| Section | Options |
|---|---|
| **Source data from Step 7** | Select the **RF dataset** and **RF fold** to analyze. Both selectors auto-populate from your Step 7 results when Step 7 finishes. |
| **Distance method** | Euclidean / Bray-Curtis / Robust Aitchison (default) / Manhattan. Bray-Curtis and Robust Aitchison automatically fall back to Euclidean if the data contain negative values. |
| **Linkage method** | Average (default) / Complete / Single / Ward.D2 |

Click **Run HCA** to build and display the dendrogram. The cophenetic correlation coefficient (CCC) is printed above the plot.

#### Output tabs (right)

| Tab | Content |
|---|---|
| **Dendrogram** | Cluster plot with CCC value printed at the top, samples labeled `Plastic_type-Source`. |
| **Console Log** | Run log including sample/feature counts and notices when a distance method falls back to Euclidean. |

![Step 8 — Dendrogram output with CCC](images/step8_dendro.png)

#### Exporting

- **Download Dendrogram (PNG)** — saves the currently-displayed dendrogram as a timestamped PNG file.
- **Done / Export All** — sweeps through every (dataset, fold) combination from Step 7 and writes the full CSV + PNG output for each to the export directory configured in Step 7.

---

## Troubleshooting

| Issue | Solution |
|---|---|
| "Directory not found" | Check the project root path on the Welcome tab |
| "Helper functions file not found" | Ensure the helper `.R` file is in the same directory as the app |
| Step N fails with "Run Step N-1 first" | Steps must be run sequentially |
| Step 7 is very slow | Reduce imputation/normalization/algorithm combinations, or lower `top_k` / `n_permutations` |
| Out-of-memory error | Close other applications or increase R's memory limit on Windows |
| Package install failure | Run `install.packages("package_name")` manually to see the full error |
