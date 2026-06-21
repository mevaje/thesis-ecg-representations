# =============================================================================
# RV Coefficient
# TFG - Grau en Estadística, UB/UPC
# Author: Melissa Vargas Jerez
# Description: Computes the RV coefficient along the four-axis comparison
#              framework (real vs. augmented, augmented vs. synthetic; each
#              overall and per arrhythmia class), as described in the
#              Methodology chapter of the thesis.
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
# 2. Core RV coefficient function
#    RV(R, Rp) = tr(KL) / sqrt(tr(K^2) * tr(L^2))
#    where K = R R', L = Rp Rp' are the Gram matrices of the two
#    (column-mean-centred) representation matrices. Unlike CCA, RV does not
#    require the two compared matrices to have equal N.
# -----------------------------------------------------------------------------

compute_rv <- function(feat_X, feat_Y) {
  
  X <- sweep(feat_X, 2, colMeans(feat_X), "-")
  Y <- sweep(feat_Y, 2, colMeans(feat_Y), "-")
  
  K <- tcrossprod(X)   # N_X x N_X
  L <- tcrossprod(Y)   # N_Y x N_Y
  
  tr_KL <- sum(K * L)
  tr_K2 <- sum(K * K)
  tr_L2 <- sum(L * L)
  
  RV <- tr_KL / sqrt(tr_K2 * tr_L2)
  
  list(
    RV    = RV,
    tr_KL = tr_KL,
    tr_K2 = tr_K2,
    tr_L2 = tr_L2,
    N_X   = nrow(feat_X),
    N_Y   = nrow(feat_Y)
  )
}

# -----------------------------------------------------------------------------
# 3. Axis 1 and 3: overall comparisons
#    Axis 1: Real (reference) vs. Augmented
#    Axis 3: Augmented vs. Synthetic
# -----------------------------------------------------------------------------

rv_axis1_overall <- compute_rv(feat_real, feat_aug)
rv_axis3_overall <- compute_rv(feat_aug,  feat_syn)

# -----------------------------------------------------------------------------
# 4. Axis 2 and 4: per-class comparisons
#    For each arrhythmia class, RV is computed only on the patients of that
#    class within the relevant pair of matrices. As with CCA, the RV formula
#    requires tr(KL), which is only defined when K and L are square matrices
#    of equal dimension; the two classes are therefore truncated to the
#    smaller of the two available patient counts before computing K and L.
# -----------------------------------------------------------------------------

compute_rv_per_class <- function(feat_X, class_X, feat_Y, class_Y) {
  classes <- levels(class_X)
  results <- lapply(classes, function(cl) {
    idx_X <- which(class_X == cl)
    idx_Y <- which(class_Y == cl)
    # RV requires K and L to be square matrices of equal dimension (tr(KL)
    # is undefined otherwise), so both classes are truncated to the smaller
    # of the two available counts, consistent with the CCA per-class function.
    n_use <- min(length(idx_X), length(idx_Y))
    res <- compute_rv(feat_X[idx_X[seq_len(n_use)], , drop = FALSE],
                      feat_Y[idx_Y[seq_len(n_use)], , drop = FALSE])
    tibble(Class = cl, N = res$N_X, RV = res$RV)
  })
  bind_rows(results)
}

rv_axis2_per_class <- compute_rv_per_class(feat_real, class_real,
                                           feat_aug,  class_aug)

rv_axis4_per_class <- compute_rv_per_class(feat_aug, class_aug,
                                           feat_syn, class_syn)

# -----------------------------------------------------------------------------
# 5. Summary tables
# -----------------------------------------------------------------------------

# Overall scores across both axes (1 and 3)
rv_overall_summary <- tibble(
  Axis  = c("1: Real vs. Augmented", "3: Augmented vs. Synthetic"),
  RV    = round(c(rv_axis1_overall$RV, rv_axis3_overall$RV), 6)
)

print(rv_overall_summary)

# Per-class scores, both axes combined into one tidy table
rv_per_class_summary <- bind_rows(
  rv_axis2_per_class |> mutate(Axis = "2: Real vs. Augmented", .before = 1),
  rv_axis4_per_class |> mutate(Axis = "4: Augmented vs. Synthetic", .before = 1)
) |>
  mutate(RV = round(RV, 6))

print(rv_per_class_summary)

# -----------------------------------------------------------------------------
# 6. Save results
# -----------------------------------------------------------------------------

rv_results <- list(
  axis1_overall         = rv_axis1_overall,
  axis2_per_class       = rv_axis2_per_class,
  axis3_overall         = rv_axis3_overall,
  axis4_per_class       = rv_axis4_per_class,
  rv_overall_summary    = rv_overall_summary,
  rv_per_class_summary  = rv_per_class_summary,
  D                     = D
)

saveRDS(rv_results, file = here("rv_results.rds"))