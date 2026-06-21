# =============================================================================
# Canonical Correlation Analysis (CCA)
# TFG - Grau en Estadística, UB/UPC
# Author: Melissa Vargas Jerez
# Description: Computes CCA along the four-axis comparison framework
#              (real vs. augmented, augmented vs. synthetic; each overall
#              and per arrhythmia class), as described in the Methodology
#              chapter of the thesis.
# =============================================================================

library(here)
library(tibble)
library(dplyr)

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
# All three matrices are the N = 3,000 (stratified for real) subsamples
# produced by 01_tfg_melissa_ear.R, so that all six similarity measures
# operate on identical observations.

feat_real  <- readRDS(here("feat_real_3k.rds"))    # 3000 x 256, real (reference)
feat_aug   <- readRDS(here("feat_aug_3k.rds"))      # 3000 x 256, real (augmented)
feat_syn   <- readRDS(here("feat_syn_3k.rds"))      # 3000 x 256, synthetic

class_real <- readRDS(here("class_real_3k.rds"))    # factor, length 3000
class_aug  <- readRDS(here("class_aug_3k.rds"))     # factor, length 3000
class_syn  <- readRDS(here("class_syn_3k.rds"))     # factor, length 3000

D <- ncol(feat_aug)
stopifnot(ncol(feat_real) == D, ncol(feat_syn) == D)

# -----------------------------------------------------------------------------
# 2. Core CCA function
#    Computes the full canonical correlation sequence and m_CCA between two
#    representation matrices, following the closed-form SVD solution.
# -----------------------------------------------------------------------------

compute_cca <- function(feat_X, feat_Y) {
  
  N <- nrow(feat_X)
  stopifnot(nrow(feat_Y) == N)
  
  # -- Preprocessing: column-wise mean-centring --
  X <- sweep(feat_X, 2, colMeans(feat_X), "-")
  Y <- sweep(feat_Y, 2, colMeans(feat_Y), "-")
  
  # -- Covariance matrices --
  Sigma_XX <- crossprod(X) / N
  Sigma_YY <- crossprod(Y) / N
  Sigma_XY <- crossprod(X, Y) / N
  
  # -- Symmetric matrix inverse square root via eigendecomposition --
  mat_inv_sqrt <- function(A, ridge = 1e-8) {
    eig <- eigen(A, symmetric = TRUE)
    d   <- pmax(eig$values, 0)
    d_inv_sqrt <- ifelse(d > ridge, 1 / sqrt(d), 0)
    eig$vectors %*% diag(d_inv_sqrt) %*% t(eig$vectors)
  }
  
  Sigma_XX_inv_sqrt <- mat_inv_sqrt(Sigma_XX)
  Sigma_YY_inv_sqrt <- mat_inv_sqrt(Sigma_YY)
  
  # -- Whitened cross-covariance and SVD --
  M     <- Sigma_XX_inv_sqrt %*% Sigma_XY %*% Sigma_YY_inv_sqrt
  svd_M <- svd(M)
  
  canonical_correlations <- pmin(pmax(svd_M$d, 0), 1)
  m_CCA <- mean(canonical_correlations)
  
  list(
    m_CCA                   = m_CCA,
    canonical_correlations  = canonical_correlations,
    N                       = N,
    D                       = ncol(feat_X)
  )
}

# -----------------------------------------------------------------------------
# 3. Axis 1 and 3: overall comparisons
#    Axis 1: Real (reference) vs. Augmented
#    Axis 3: Augmented vs. Synthetic
# -----------------------------------------------------------------------------

cca_axis1_overall <- compute_cca(feat_real, feat_aug)
cca_axis3_overall <- compute_cca(feat_aug,  feat_syn)

# -----------------------------------------------------------------------------
# 4. Axis 2 and 4: per-class comparisons
#    For each arrhythmia class, CCA is computed only on the patients of that
#    class within the relevant pair of matrices.
# -----------------------------------------------------------------------------

compute_cca_per_class <- function(feat_X, class_X, feat_Y, class_Y) {
  classes <- levels(class_X)
  results <- lapply(classes, function(cl) {
    idx_X <- which(class_X == cl)
    idx_Y <- which(class_Y == cl)
    n_use <- min(length(idx_X), length(idx_Y))
    # Use the same number of observations from each class in both matrices,
    # since CCA requires equal N between the two compared matrices.
    res <- compute_cca(feat_X[idx_X[seq_len(n_use)], , drop = FALSE],
                       feat_Y[idx_Y[seq_len(n_use)], , drop = FALSE])
    tibble(Class = cl, N = res$N, m_CCA = res$m_CCA)
  })
  bind_rows(results)
}

cca_axis2_per_class <- compute_cca_per_class(feat_real, class_real,
                                             feat_aug,  class_aug)

cca_axis4_per_class <- compute_cca_per_class(feat_aug, class_aug,
                                             feat_syn, class_syn)

# -----------------------------------------------------------------------------
# 5. Summary tables
# -----------------------------------------------------------------------------

# Overall scores across both axes (1 and 3)
cca_overall_summary <- tibble(
  Axis        = c("1: Real vs. Augmented", "3: Augmented vs. Synthetic"),
  m_CCA       = round(c(cca_axis1_overall$m_CCA, cca_axis3_overall$m_CCA), 6),
  rho_1       = round(c(cca_axis1_overall$canonical_correlations[1],
                        cca_axis3_overall$canonical_correlations[1]), 6),
  rho_D       = round(c(cca_axis1_overall$canonical_correlations[D],
                        cca_axis3_overall$canonical_correlations[D]), 6)
)

print(cca_overall_summary)

# Per-class scores, both axes combined into one tidy table
cca_per_class_summary <- bind_rows(
  cca_axis2_per_class |> mutate(Axis = "2: Real vs. Augmented", .before = 1),
  cca_axis4_per_class |> mutate(Axis = "4: Augmented vs. Synthetic", .before = 1)
) |>
  mutate(m_CCA = round(m_CCA, 6))

print(cca_per_class_summary)

# -----------------------------------------------------------------------------
# 6. Save results
# -----------------------------------------------------------------------------

cca_results <- list(
  axis1_overall        = cca_axis1_overall,
  axis2_per_class      = cca_axis2_per_class,
  axis3_overall        = cca_axis3_overall,
  axis4_per_class      = cca_axis4_per_class,
  cca_overall_summary  = cca_overall_summary,
  cca_per_class_summary = cca_per_class_summary,
  D                    = D
)

saveRDS(cca_results, file = here("cca_results.rds"))