# =============================================================================
# 07_tfg_melissa_jaccard.R  --  k-NN Jaccard Similarity
# TFG - Grau en Estadística, UB/UPC
# Author: Melissa Vargas Jerez
# Description: Computes the k-NN Jaccard similarity between the real
#              (augmented) and synthetic representation matrices extracted
#              from the frozen 1D-CNN. Follows the methodology described
#              in Section 4.6.
#              Requires feat_aug_3k.rds and feat_syn_3k.rds produced by
#              01_tfg_melissa_eda.R.
# =============================================================================

library(here)
library(tibble)

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------

feat_aug <- readRDS(here("feat_aug_3k.rds"))   # 3000 x 256, real (augmented)
feat_syn <- readRDS(here("feat_syn_3k.rds"))   # 3000 x 256, synthetic

N <- nrow(feat_aug)
D <- ncol(feat_aug)

stopifnot(nrow(feat_syn) == N, ncol(feat_syn) == D)

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

R_norm  <- row_normalise(feat_aug)
Rp_norm <- row_normalise(feat_syn)

# -----------------------------------------------------------------------------
# 3. Core functions
# -----------------------------------------------------------------------------

# -- 3.1 k nearest neighbours for all patients --
# For unit-normalised rows, cosine similarity = dot product = R_norm %*% t(R_norm).
# Returns an N x k integer matrix of neighbour indices (self excluded).
compute_knn <- function(mat_norm, k) {
  sim_mat <- tcrossprod(mat_norm)   # N x N cosine similarity matrix
  diag(sim_mat) <- -Inf             # exclude self
  # For each row, find the indices of the k largest similarities
  t(apply(sim_mat, 1, function(row) order(row, decreasing = TRUE)[seq_len(k)]))
}

# -- 3.2 Instance-wise Jaccard index --
# v_i^k = |N_R^k(i) ∩ N_Rp^k(i)| / |N_R^k(i) ∪ N_Rp^k(i)|
jaccard_instance <- function(set_a, set_b) {
  intersection <- length(intersect(set_a, set_b))
  union        <- length(union(set_a, set_b))
  intersection / union
}

# -- 3.3 Global k-NN Jaccard score --
# m_Jac^k = (1/N) * sum_i v_i^k
compute_jaccard <- function(R_norm, Rp_norm, k) {
  knn_R  <- compute_knn(R_norm,  k)
  knn_Rp <- compute_knn(Rp_norm, k)
  instance_scores <- vapply(
    seq_len(N),
    function(i) jaccard_instance(knn_R[i, ], knn_Rp[i, ]),
    numeric(1)
  )
  list(
    m_Jac            = mean(instance_scores),
    instance_scores  = instance_scores
  )
}

# -----------------------------------------------------------------------------
# 4. Compute scores for k in {5, 10, 20, 50}
#    k = 10 is the primary value; others show sensitivity to neighbourhood size.
# -----------------------------------------------------------------------------

K_VALUES <- c(5L, 10L, 20L, 50L)

jaccard_by_k <- lapply(K_VALUES, function(k) {
  cat(sprintf("Computing k-NN Jaccard for k = %d ...\n", k))
  res <- compute_jaccard(R_norm, Rp_norm, k)
  list(k = k, m_Jac = res$m_Jac, instance_scores = res$instance_scores)
})
names(jaccard_by_k) <- paste0("k", K_VALUES)

# -----------------------------------------------------------------------------
# 5. Summary tables
# -----------------------------------------------------------------------------

# Global scores across all k values
jaccard_summary <- tibble(
  k          = K_VALUES,
  m_Jac      = round(sapply(jaccard_by_k, `[[`, "m_Jac"), 6),
  Primary    = ifelse(K_VALUES == 10L, "yes", "no")
)

print(jaccard_summary)

# Detailed summary for primary k = 10
k10          <- jaccard_by_k[["k10"]]
k10_scores   <- k10$instance_scores

jaccard_k10_summary <- tibble(
  Metric = c("m_Jac (k = 10, primary)",
             "Standard deviation across patients",
             "Minimum instance score",
             "Median instance score",
             "Maximum instance score"),
  Value  = c(round(k10$m_Jac,       6),
             round(sd(k10_scores),   6),
             round(min(k10_scores),  6),
             round(median(k10_scores), 6),
             round(max(k10_scores),  6))
)

print(jaccard_k10_summary)

# -----------------------------------------------------------------------------
# 6. Save results
# -----------------------------------------------------------------------------

jaccard_results <- list(
  jaccard_by_k        = jaccard_by_k,      # results for all k values
  jaccard_summary     = jaccard_summary,   # global scores table
  jaccard_k10_summary = jaccard_k10_summary,
  K_VALUES            = K_VALUES,
  N                   = N,
  D                   = D
)

saveRDS(jaccard_results, file = here("jaccard_results.rds"))
message("Results saved to jaccard_results.rds")
