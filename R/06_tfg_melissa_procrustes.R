# =============================================================================
# 06_tfg_melissa_procrustes.R  - Orthogonal Procrustes
# TFG - Grau en Estadística, UB/UPC
# Author: Melissa Vargas Jerez
# Description: Computes the Orthogonal Procrustes distance and similarity
#              along the four-axis comparison framework (real vs. augmented,
#              augmented vs. synthetic; each overall and per arrhythmia
#              class), as described in the Methodology chapter of the
#              thesis.
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
# 2. Core Orthogonal Procrustes function
#    Preprocessing: column-wise mean-centring, then scaling to unit
#    Frobenius norm. Notation: tilde matrices reserved for Procrustes only.
#
#    Cross-product matrix: M = R_tilde' Rp_tilde   (D x D)
#    SVD: M = U Sigma V', optimal rotation Q* = U V'
#    Distance (after unit-norm scaling): m_Proc = sqrt(2 - 2 tr(Sigma))
#    Similarity: sim_Proc = 1 - m_Proc / sqrt(2)
#
#    Like CCA, RV, CKA, and RSA, this measure requires feat_X and feat_Y to
#    have the same number of rows (N), since crossprod(X, Y) is only defined
#    when X and Y share that dimension.
# -----------------------------------------------------------------------------

frob_norm <- function(A) sqrt(sum(A^2))

compute_procrustes <- function(feat_X, feat_Y) {
  
  X_c <- sweep(feat_X, 2, colMeans(feat_X), "-")
  Y_c <- sweep(feat_Y, 2, colMeans(feat_Y), "-")
  
  X_tilde <- X_c / frob_norm(X_c)
  Y_tilde <- Y_c / frob_norm(Y_c)
  
  M     <- crossprod(X_tilde, Y_tilde)   # D x D
  svd_M <- svd(M)
  
  Q_star   <- svd_M$u %*% t(svd_M$v)
  tr_Sigma <- sum(svd_M$d)
  
  m_Proc   <- sqrt(max(0, 2 - 2 * tr_Sigma))
  sim_Proc <- 1 - m_Proc / sqrt(2)
  
  list(
    m_Proc          = m_Proc,
    sim_Proc        = sim_Proc,
    tr_Sigma        = tr_Sigma,
    Q_star          = Q_star,
    singular_values = svd_M$d,
    N_X             = nrow(feat_X),
    N_Y             = nrow(feat_Y)
  )
}

# -----------------------------------------------------------------------------
# 3. Axis 1 and 3: overall comparisons
#    Axis 1: Real (reference) vs. Augmented
#    Axis 3: Augmented vs. Synthetic
# -----------------------------------------------------------------------------

proc_axis1_overall <- compute_procrustes(feat_real, feat_aug)
proc_axis3_overall <- compute_procrustes(feat_aug,  feat_syn)

# -----------------------------------------------------------------------------
# 4. Axis 2 and 4: per-class comparisons
#    For each arrhythmia class, Procrustes is computed only on the patients
#    of that class within the relevant pair of matrices. As with CCA, RV,
#    CKA, and RSA, crossprod(X_tilde, Y_tilde) requires X and Y to have the
#    same number of rows, so each class is matched to the smaller of the two
#    available counts via random subsampling.
# -----------------------------------------------------------------------------

compute_procrustes_per_class <- function(feat_X, class_X, feat_Y, class_Y,
                                         min_n = 4L, seed = 42) {
  classes <- levels(class_X)
  results <- map_dfr(classes, function(cl) {
    idx_X <- which(class_X == cl)
    idx_Y <- which(class_Y == cl)
    
    if (length(idx_X) < min_n || length(idx_Y) < min_n) {
      warning(sprintf("Class %s has fewer than %d samples in one dataset, skipping.",
                      cl, min_n))
      return(NULL)
    }
    
    # crossprod(X_tilde, Y_tilde) requires X and Y to have the same number
    # of rows, so each class is matched to the smaller of the two available
    # counts via random subsampling, consistent with the approach used for
    # CCA, RV, CKA, and RSA.
    n_use <- min(length(idx_X), length(idx_Y))
    set.seed(seed)
    idx_X_use <- sample(idx_X, n_use)
    idx_Y_use <- sample(idx_Y, n_use)
    
    res <- compute_procrustes(feat_X[idx_X_use, , drop = FALSE],
                              feat_Y[idx_Y_use, , drop = FALSE])
    
    tibble(Class = cl, N_used = n_use,
           m_Proc = res$m_Proc, sim_Proc = res$sim_Proc)
  })
  results
}

proc_axis2_per_class <- compute_procrustes_per_class(feat_real, class_real,
                                                     feat_aug,  class_aug)

proc_axis4_per_class <- compute_procrustes_per_class(feat_aug, class_aug,
                                                     feat_syn, class_syn)

# -----------------------------------------------------------------------------
# 5. Summary tables
# -----------------------------------------------------------------------------

# Overall scores across both axes (1 and 3)
proc_overall_summary <- tibble(
  Axis     = c("1: Real vs. Augmented", "3: Augmented vs. Synthetic"),
  m_Proc   = round(c(proc_axis1_overall$m_Proc,   proc_axis3_overall$m_Proc), 6),
  sim_Proc = round(c(proc_axis1_overall$sim_Proc, proc_axis3_overall$sim_Proc), 6)
)

print(proc_overall_summary)

# Per-class scores, both axes combined into one tidy table
proc_per_class_summary <- bind_rows(
  proc_axis2_per_class |> mutate(Axis = "2: Real vs. Augmented", .before = 1),
  proc_axis4_per_class |> mutate(Axis = "4: Augmented vs. Synthetic", .before = 1)
) |>
  mutate(m_Proc = round(m_Proc, 6),
         sim_Proc = round(sim_Proc, 6))

print(proc_per_class_summary)

# -----------------------------------------------------------------------------
# 6. Save results
# -----------------------------------------------------------------------------

procrustes_results <- list(
  axis1_overall          = proc_axis1_overall,
  axis2_per_class        = proc_axis2_per_class,
  axis3_overall           = proc_axis3_overall,
  axis4_per_class        = proc_axis4_per_class,
  proc_overall_summary   = proc_overall_summary,
  proc_per_class_summary = proc_per_class_summary,
  D                      = D
)

saveRDS(procrustes_results, file = here("procrustes_results.rds"))