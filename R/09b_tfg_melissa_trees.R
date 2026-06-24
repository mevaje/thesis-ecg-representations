# =============================================================================
# 9b_tfg_melissa_tree_plots.R  --  Tree visualisations
# TFG - Grau en Estadística, UB/UPC
# Author: Melissa Vargas Jerez
# Description: Loads the fitted and pruned classification trees saved in
#              trees_results.rds and produces publication-ready plots for
#              the thesis. Four plots are generated:
#              (1) Tree predicting arrhythmia class on the real dataset
#              (2) Tree predicting arrhythmia class on the augmented dataset
#              (3) Tree predicting arrhythmia class on the synthetic dataset
#              (4) Tree predicting dataset origin on the pooled matrix
#              Additionally, a combined feature importance bar chart is
#              produced for the class prediction trees.
#              Requires trees_results.rds produced by 09_tfg_melissa_trees.R
# =============================================================================

library(here)
library(rpart)
library(rpart.plot)
library(ggplot2)
library(dplyr)
library(tidyr)

# -----------------------------------------------------------------------------
# 0. Colour palette (consistent with EDA and similarity measure scripts)
# -----------------------------------------------------------------------------

DATASET_COLOURS <- c(
  "Real"      = "#2E5EA8",
  "Augmented" = "#E07B39",
  "Synthetic" = "#3BA87A"
)

CLASS_COLOURS <- c(
  "SBRAD" = "#4393C3",
  "SR"    = "#92C5DE",
  "AFIB"  = "#D6604D",
  "STACH" = "#F4A582",
  "SARRH" = "#878787",
  "AFLT"  = "#4DAC26",
  "SVTAC" = "#B8E186"
)

# -----------------------------------------------------------------------------
# 1. Load results
# -----------------------------------------------------------------------------

trees_results <- readRDS(here("trees_results.rds"))

tree_real   <- trees_results$tree_real
tree_aug    <- trees_results$tree_aug
tree_syn    <- trees_results$tree_syn
tree_origin <- trees_results$tree_origin

imp_real     <- trees_results$imp_real
imp_aug      <- trees_results$imp_aug
imp_syn      <- trees_results$imp_syn
imp_combined <- trees_results$imp_combined
imp_origin   <- trees_results$imp_origin

conf_origin      <- trees_results$conf_origin
accuracy_summary <- trees_results$accuracy_summary

# Create figures directory if it does not exist
if (!dir.exists(here("figures"))) dir.create(here("figures"))

# -----------------------------------------------------------------------------
# 2. Helper: save an rpart.plot to PDF
# -----------------------------------------------------------------------------

save_tree_plot <- function(tree_fit, filename, title) {
  pdf(here("figures", filename), width = 12, height = 7)
  rpart.plot(
    tree_fit,
    type        = 2,        # split label above node, class label inside
    extra       = 104,      # show probability and percentage of obs at node
    under       = TRUE,     # place extra info under the node box
    fallen.leaves = TRUE,   # leaves at the bottom for readability
    main        = title,
    cex         = 0.7,      # font size
    box.palette = "auto"    # automatic colour by class
  )
  dev.off()
  message("Saved: ", filename)
}

# -----------------------------------------------------------------------------
# 3. Tree plots: class prediction
# -----------------------------------------------------------------------------

save_tree_plot(
  tree_real,
  "tree_class_real.pdf",
  "Classification tree: arrhythmia class | Real representations"
)

save_tree_plot(
  tree_aug,
  "tree_class_aug.pdf",
  "Classification tree: arrhythmia class | Augmented representations"
)

save_tree_plot(
  tree_syn,
  "tree_class_syn.pdf",
  "Classification tree: arrhythmia class | Synthetic representations"
)

# -----------------------------------------------------------------------------
# 4. Tree plot: origin prediction
# -----------------------------------------------------------------------------

save_tree_plot(
  tree_origin,
  "tree_origin.pdf",
  "Classification tree: dataset origin | Pooled representations (N = 9,000)"
)


