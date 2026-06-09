# =============================================================================
# RV Coefficient
# TFG - Grau en Estadística, UB/UPC
# Author: Melissa Vargas Jerez
# Description: Computes the RV coefficient between the real (augmented) and
#              synthetic representation matrices extracted from the frozen
#              1D-CNN.
# =============================================================================

library(here)
library(tibble)

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
# Both matrices are the N = 3,000 subsamples produced by 01_tfg_melissa_eda.R,
# drawn from the same-sized primary matrices (N = 18,027 each) so that all
# six similarity measures operate on identical observations.

feat_aug <- readRDS(here("feat_aug_3k.rds"))   # 3000 x 256, real (augmented)
feat_syn <- readRDS(here("feat_syn_3k.rds"))   # 3000 x 256, synthetic

N <- nrow(feat_aug)
D <- ncol(feat_aug)

stopifnot(nrow(feat_syn) == N, ncol(feat_syn) == D)

# -----------------------------------------------------------------------------
# 2. Preprocessing: column-wise mean-centring
# -----------------------------------------------------------------------------

R  <- sweep(feat_aug, 2, colMeans(feat_aug), "-")
Rp <- sweep(feat_syn, 2, colMeans(feat_syn), "-")

# -----------------------------------------------------------------------------
# 3. Gram matrices
#    K = R R'  (N x N), entry K_ij = dot product of patients i and j in R
#    L = Rp Rp' (N x N), entry L_ij = dot product of patients i and j in Rp
# -----------------------------------------------------------------------------

K <- tcrossprod(R)    # N x N
L <- tcrossprod(Rp)   # N x N

# -----------------------------------------------------------------------------
# 4. RV coefficient
#    RV(R, Rp) = tr(KL) / sqrt(tr(K^2) * tr(L^2))
#
#    Efficient computation:
#      tr(KL)  = sum(K * L)          (elementwise product then sum)
#      tr(K^2) = sum(K * K) = ||K||_F^2
#      tr(L^2) = sum(L * L) = ||L||_F^2
# -----------------------------------------------------------------------------

tr_KL  <- sum(K * L)
tr_K2  <- sum(K * K)
tr_L2  <- sum(L * L)

RV <- tr_KL / sqrt(tr_K2 * tr_L2)

# -----------------------------------------------------------------------------
# 5. Summary table
# -----------------------------------------------------------------------------

rv_summary <- tibble(
  Metric = c("RV coefficient",
             "tr(KL)  — numerator",
             "tr(K^2) — real Gram norm squared",
             "tr(L^2) — synthetic Gram norm squared"),
  Value  = c(round(RV, 6),
             round(tr_KL, 2),
             round(tr_K2, 2),
             round(tr_L2, 2))
)

print(rv_summary)

# -----------------------------------------------------------------------------
# 6. Save results
# -----------------------------------------------------------------------------

rv_results <- list(
  RV         = RV,
  tr_KL      = tr_KL,
  tr_K2      = tr_K2,
  tr_L2      = tr_L2,
  rv_summary = rv_summary,
  N          = N,
  D          = D
)

saveRDS(rv_results, file = here("rv_results.rds"))

