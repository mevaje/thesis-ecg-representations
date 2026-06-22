# =============================================================================
# 09_tfg_melissa_trees.R  --  Classification Trees (Diagnostic Methods)
# TFG - Grau en Estadística, UB/UPC
# Author: Melissa Vargas Jerez
# Description: Fits classification trees as diagnostic tools to investigate
#              the sources of representational differences across the three
#              datasets. Two analyses are performed:
#              (1) Three trees predicting arrhythmia class separately on each
#                  dataset, to compare which activation features each relies on
#                  for class discrimination.
#              (2) One tree predicting dataset origin (real, augmented,
#                  synthetic) on the pooled 9,000 x 256 matrix, to identify
#                  which features best separate the three datasets.
#              Requires feat_real_3k.rds, feat_aug_3k.rds, feat_syn_3k.rds,
#              class_real_3k.rds, class_aug_3k.rds, class_syn_3k.rds
#              produced by 01_tfg_melissa_ear.R
# =============================================================================

library(here)
library(tibble)
library(dplyr)
library(rpart)
library(rpart.plot)

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------

feat_real  <- readRDS(here("feat_real_3k.rds"))
feat_aug   <- readRDS(here("feat_aug_3k.rds"))
feat_syn   <- readRDS(here("feat_syn_3k.rds"))

class_real <- readRDS(here("class_real_3k.rds"))
class_aug  <- readRDS(here("class_aug_3k.rds"))
class_syn  <- readRDS(here("class_syn_3k.rds"))

D <- ncol(feat_aug)
stopifnot(ncol(feat_real) == D, ncol(feat_syn) == D)

# Helper: build a labelled data frame from a feature matrix and a response
build_df <- function(feat_matrix, response) {
  df <- as.data.frame(feat_matrix)
  colnames(df) <- paste0("V", seq_len(ncol(df)))
  df$y <- response
  df
}

# Helper: normalised feature importance from a fitted rpart tree
get_importance <- function(tree_fit, top_n = 20L) {
  imp <- tree_fit$variable.importance
  if (is.null(imp)) return(tibble(Feature = character(), Importance = numeric()))
  imp_norm <- round(100 * imp / sum(imp), 2)
  tibble(
    Feature    = names(imp_norm),
    Importance = as.numeric(imp_norm)
  ) |>
    arrange(desc(Importance)) |>
    slice_head(n = top_n)
}

# Helper: fit an rpart tree with full growth then cross-validated pruning
# following the 1-SE rule (choose the simplest tree within 1 SE of the
# minimum cross-validated error), consistent with the methodology described
# in the Theoretical Background chapter.
fit_pruned_tree <- function(df, seed = 42) {
  set.seed(seed)
  
  # Grow a full tree (cp = 0)
  full_tree <- rpart(y ~ ., data = df, method = "class",
                     control = rpart.control(cp = 0, minsplit = 5))
  
  # Cross-validated error table
  cp_table <- full_tree$cptable
  
  # Minimum xerror and 1-SE threshold
  min_idx   <- which.min(cp_table[, "xerror"])
  threshold <- cp_table[min_idx, "xerror"] + cp_table[min_idx, "xstd"]
  
  # Smallest tree (largest cp) whose xerror is still within the threshold
  best_idx <- min(which(cp_table[, "xerror"] <= threshold))
  best_cp  <- cp_table[best_idx, "CP"]
  
  pruned <- prune(full_tree, cp = best_cp)
  list(tree = pruned, cp_table = as_tibble(as.data.frame(cp_table)))
}

# -----------------------------------------------------------------------------
# 2. Analysis 1: predicting arrhythmia class within each dataset
#    One tree per dataset. Feature importance profiles are compared to
#    identify which activation dimensions each dataset relies on for class
#    discrimination.
# -----------------------------------------------------------------------------

df_real <- build_df(feat_real, class_real)
df_aug  <- build_df(feat_aug,  class_aug)
df_syn  <- build_df(feat_syn,  class_syn)

tree_real <- fit_pruned_tree(df_real)
tree_aug  <- fit_pruned_tree(df_aug)
tree_syn  <- fit_pruned_tree(df_syn)

imp_real <- get_importance(tree_real$tree) |> rename(Importance_Real = Importance)
imp_aug  <- get_importance(tree_aug$tree)  |> rename(Importance_Aug  = Importance)
imp_syn  <- get_importance(tree_syn$tree)  |> rename(Importance_Syn  = Importance)

# Top-20 feature importance per dataset
print(imp_real)
print(imp_aug)
print(imp_syn)

# Joined comparison: all three datasets side by side on the union of top features
imp_combined <- full_join(imp_real, imp_aug,  by = "Feature") |>
  full_join(imp_syn,             by = "Feature") |>
  arrange(desc(coalesce(Importance_Real, 0) +
                 coalesce(Importance_Aug,  0) +
                 coalesce(Importance_Syn,  0)))

print(imp_combined)

# -----------------------------------------------------------------------------
# 3. Analysis 2: predicting dataset origin on the pooled matrix
#    All three N = 3,000 matrices are combined into one 9,000 x 256 input.
#    The response variable is dataset origin (real, augmented, synthetic).
#    High accuracy would indicate that the three representation spaces are
#    linearly separable by a small number of activation features.
# -----------------------------------------------------------------------------

# Build pooled data frame with dataset origin as response
df_pooled <- bind_rows(
  build_df(feat_real, factor(rep("Real",      nrow(feat_real)),
                             levels = c("Real", "Augmented", "Synthetic"))),
  build_df(feat_aug,  factor(rep("Augmented", nrow(feat_aug)),
                             levels = c("Real", "Augmented", "Synthetic"))),
  build_df(feat_syn,  factor(rep("Synthetic", nrow(feat_syn)),
                             levels = c("Real", "Augmented", "Synthetic")))
)

tree_origin <- fit_pruned_tree(df_pooled)

imp_origin <- get_importance(tree_origin$tree, top_n = 20L) |>
  rename(Importance_Origin = Importance)

print(imp_origin)

# In-sample confusion matrix (indicative of separability, not generalisation)
pred_origin <- predict(tree_origin$tree, df_pooled, type = "class")
conf_origin <- table(Predicted = pred_origin, Actual = df_pooled$y)
print(conf_origin)

# Overall accuracy
accuracy_origin <- round(sum(diag(conf_origin)) / sum(conf_origin) * 100, 2)
accuracy_summary <- tibble(
  Model    = "Origin prediction (pooled, in-sample)",
  Accuracy = accuracy_origin
)
print(accuracy_summary)

# -----------------------------------------------------------------------------
# 4. Save results
# -----------------------------------------------------------------------------

trees_results <- list(
  # Analysis 1: class prediction per dataset
  tree_real          = tree_real$tree,
  tree_aug           = tree_aug$tree,
  tree_syn           = tree_syn$tree,
  imp_real           = imp_real,
  imp_aug            = imp_aug,
  imp_syn            = imp_syn,
  imp_combined       = imp_combined,
  
  # Analysis 2: origin prediction on pooled data
  tree_origin        = tree_origin$tree,
  imp_origin         = imp_origin,
  conf_origin        = conf_origin,
  accuracy_summary   = accuracy_summary
)

saveRDS(trees_results, file = here("trees_results.rds"))
message("Results saved to trees_results.rds")