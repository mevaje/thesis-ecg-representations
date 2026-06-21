# =============================================================================
# 07_tfg_melissa_jaccard.R  --  k-NN Jaccard Similarity
# TFG - Grau en Estadística, UB/UPC
# Author: Melissa Vargas Jerez
# Description: Computes k-NN Jaccard similarity along the four-axis
#              comparison framework (real vs. augmented, augmented vs.
#              synthetic; each overall and per arrhythmia class), as
#              described in the Methodology chapter of the thesis.
#              Requires feat_real_3k.rds, feat_aug_3k.rds, feat_syn_3k.rds,
#              class_real_3k.rds, class_aug_3k.rds, class_syn_3k.rds
#              produced by 01_tfg_melissa_ear.R
# =============================================================================

library(here)
library(tibble)
library(dplyr)
library(purrr)

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

# -----------------------------------------------------------------------------
# 2. Preprocessing: row-normalisation to unit length
#    Equivalent to projecting all patients onto the unit hypersphere.
#    Cosine similarity between two unit vectors reduces to their dot product.
#    No mean-centring is applied; cosine similarity is already scale-invariant.
# -----------------------------------------------------------------------------

row_normalise <- function(mat) {
  norms <- sqrt(rowSums(mat^2))
  norms <- pmax(norms, 1e-10)   # guard against zero-norm rows
  mat / norms
}

# -----------------------------------------------------------------------------
# 3. Core functions
# -----------------------------------------------------------------------------

# -- 3.1 k nearest neighbours for all patients in one matrix --
# For unit-normalised rows, cosine similarity = dot product. Returns an
# N x k integer matrix of neighbour indices (self excluded, indices are
# local to the matrix passed in).
compute_knn <- function(mat_norm, k) {
  sim_mat <- tcrossprod(mat_norm)   # N x N cosine similarity matrix
  diag(sim_mat) <- -Inf             # exclude self
  t(apply(sim_mat, 1, function(row) order(row, decreasing = TRUE)[seq_len(k)]))
}

# -- 3.2 Instance-wise Jaccard index --
# v_i^k = |N_R^k(i) cap N_Rp^k(i)| / |N_R^k(i) cup N_Rp^k(i)|
jaccard_instance <- function(set_a, set_b) {
  intersection <- length(intersect(set_a, set_b))
  union        <- length(union(set_a, set_b))
  intersection / union
}

# -- 3.3 k-NN Jaccard score between two (already row-normalised) matrices --
# Requires both matrices to have the same N, since patient i in mat_X and
# patient i in mat_Y must refer to the same comparison slot for the
# instance-wise Jaccard index to be meaningful.
compute_jaccard <- function(mat_X_norm, mat_Y_norm, k) {
  stopifnot(nrow(mat_X_norm) == nrow(mat_Y_norm))
  n <- nrow(mat_X_norm)
  
  knn_X <- compute_knn(mat_X_norm, k)
  knn_Y <- compute_knn(mat_Y_norm, k)
  
  instance_scores <- vapply(
    seq_len(n),
    function(i) jaccard_instance(knn_X[i, ], knn_Y[i, ]),
    numeric(1)
  )
  
  list(
    m_Jac           = mean(instance_scores),
    instance_scores = instance_scores,
    N               = n
  )
}

# -----------------------------------------------------------------------------
# 4. Axis 1 and 3: overall comparisons, for k in {5, 10, 20, 50}
#    k = 10 is the primary value; others show sensitivity to neighbourhood
#    size, following the practice in the literature (Schumacher et al., 2021).
# -----------------------------------------------------------------------------

K_VALUES <- c(5L, 10L, 20L, 50L)

compute_jaccard_overall_by_k <- function(feat_X, feat_Y, k_values = K_VALUES) {
  X_norm <- row_normalise(feat_X)
  Y_norm <- row_normalise(feat_Y)
  
  by_k <- lapply(k_values, function(k) {
    res <- compute_jaccard(X_norm, Y_norm, k)
    list(k = k, m_Jac = res$m_Jac, instance_scores = res$instance_scores,
         N = res$N)
  })
  names(by_k) <- paste0("k", k_values)
  by_k
}

jaccard_axis1_overall <- compute_jaccard_overall_by_k(feat_real, feat_aug)
jaccard_axis3_overall <- compute_jaccard_overall_by_k(feat_aug,  feat_syn)

# -----------------------------------------------------------------------------
# 5. Axis 2 and 4: per-class comparisons, primary k = 10 only
#    Per-class neighbourhoods are computed within each class subset; the
#    neighbour indices returned by compute_knn() are local to that subset.
#    As with the other five measures, each class is matched to the smaller
#    of the two available counts via random subsampling, since the
#    instance-wise comparison requires both matrices to have the same N.
#    Only the primary k = 10 is reported per class, to keep the per-class
#    output (7 classes x 2 axes) a manageable size; the full sensitivity
#    analysis across k is reserved for the overall comparison.
# -----------------------------------------------------------------------------

K_PRIMARY <- 10L

compute_jaccard_per_class <- function(feat_X, class_X, feat_Y, class_Y,
                                      k = K_PRIMARY, min_n = 4L, seed = 42) {
  # min_n must be at least k + 1, since each patient needs k neighbours
  # excluding itself within its own class subset.
  min_n <- max(min_n, k + 1L)
  
  classes <- levels(class_X)
  results <- map_dfr(classes, function(cl) {
    idx_X <- which(class_X == cl)
    idx_Y <- which(class_Y == cl)
    
    if (length(idx_X) < min_n || length(idx_Y) < min_n) {
      warning(sprintf("Class %s has fewer than %d samples in one dataset, skipping.",
                      cl, min_n))
      return(NULL)
    }
    
    n_use <- min(length(idx_X), length(idx_Y))
    set.seed(seed)
    idx_X_use <- sample(idx_X, n_use)
    idx_Y_use <- sample(idx_Y, n_use)
    
    X_norm <- row_normalise(feat_X[idx_X_use, , drop = FALSE])
    Y_norm <- row_normalise(feat_Y[idx_Y_use, , drop = FALSE])
    
    res <- compute_jaccard(X_norm, Y_norm, k)
    
    tibble(Class = cl, N_used = n_use, k = k, m_Jac = res$m_Jac)
  })
  results
}

jaccard_axis2_per_class <- compute_jaccard_per_class(feat_real, class_real,
                                                     feat_aug,  class_aug)

jaccard_axis4_per_class <- compute_jaccard_per_class(feat_aug, class_aug,
                                                     feat_syn, class_syn)

# -----------------------------------------------------------------------------
# 6. Summary tables
# -----------------------------------------------------------------------------

# Overall scores across both axes (1 and 3), for all k values
jaccard_overall_summary <- bind_rows(
  tibble(Axis = "1: Real vs. Augmented",
         k    = K_VALUES,
         m_Jac = round(sapply(jaccard_axis1_overall, `[[`, "m_Jac"), 6)),
  tibble(Axis = "3: Augmented vs. Synthetic",
         k    = K_VALUES,
         m_Jac = round(sapply(jaccard_axis3_overall, `[[`, "m_Jac"), 6))
) |>
  mutate(Primary = ifelse(k == K_PRIMARY, "yes", "no"))

print(jaccard_overall_summary)

# Detailed summary for primary k = 10, both overall axes
jaccard_k10_summary <- tibble(
  Axis  = c("1: Real vs. Augmented", "3: Augmented vs. Synthetic"),
  m_Jac = round(c(jaccard_axis1_overall[["k10"]]$m_Jac,
                  jaccard_axis3_overall[["k10"]]$m_Jac), 6),
  SD    = round(c(sd(jaccard_axis1_overall[["k10"]]$instance_scores),
                  sd(jaccard_axis3_overall[["k10"]]$instance_scores)), 6),
  Median = round(c(median(jaccard_axis1_overall[["k10"]]$instance_scores),
                   median(jaccard_axis3_overall[["k10"]]$instance_scores)), 6)
)

print(jaccard_k10_summary)

# Per-class scores (k = 10 only), both axes combined into one tidy table
jaccard_per_class_summary <- bind_rows(
  jaccard_axis2_per_class |> mutate(Axis = "2: Real vs. Augmented", .before = 1),
  jaccard_axis4_per_class |> mutate(Axis = "4: Augmented vs. Synthetic", .before = 1)
) |>
  mutate(m_Jac = round(m_Jac, 6))

print(jaccard_per_class_summary)

# -----------------------------------------------------------------------------
# 7. Save results
# -----------------------------------------------------------------------------

jaccard_results <- list(
  axis1_overall            = jaccard_axis1_overall,
  axis2_per_class          = jaccard_axis2_per_class,
  axis3_overall            = jaccard_axis3_overall,
  axis4_per_class          = jaccard_axis4_per_class,
  jaccard_overall_summary  = jaccard_overall_summary,
  jaccard_k10_summary      = jaccard_k10_summary,
  jaccard_per_class_summary = jaccard_per_class_summary,
  K_VALUES                 = K_VALUES,
  K_PRIMARY                = K_PRIMARY,
  D                        = D
)

saveRDS(jaccard_results, file = here("jaccard_results.rds"))
message("Results saved to jaccard_results.rds")