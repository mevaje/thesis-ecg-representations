# =============================================================================
# 06_tfg_melissa_procrustes.R  - Orthogonal Procrustes
# TFG - Grau en Estadística, UB/UPC
# Author: Melissa Vargas Jerez
# Description: Computes the Orthogonal Procrustes distance and similarity
#              between the real (augmented) and synthetic representation
#              matrices extracted from the frozen 1D-CNN.
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
# 2. Preprocessing
#    Step 1: column-wise mean-centring
#    Step 2: scale to unit Frobenius norm
#    Notation: tilde matrices R_tilde and Rp_tilde reserved for Procrustes only
# -----------------------------------------------------------------------------

# Step 1: mean-centre
R_c  <- sweep(feat_aug, 2, colMeans(feat_aug), "-")
Rp_c <- sweep(feat_syn, 2, colMeans(feat_syn), "-")

# Step 2: scale to unit Frobenius norm
# ||A||_F = sqrt(sum(A^2))
frob_norm <- function(A) sqrt(sum(A^2))

R_tilde  <- R_c  / frob_norm(R_c)
Rp_tilde <- Rp_c / frob_norm(Rp_c)

# Verify unit norm (should both equal 1)
stopifnot(abs(frob_norm(R_tilde)  - 1) < 1e-10)
stopifnot(abs(frob_norm(Rp_tilde) - 1) < 1e-10)

# -----------------------------------------------------------------------------
# 3. Orthogonal Procrustes
#
#    Cross-product (covariance) matrix: M = R_tilde' Rp_tilde   (D x D)
#    SVD: M = U Sigma V'
#    Optimal rotation: Q* = U V'
#
#    Distance formula (simplified after unit Frobenius normalisation):
#      m_Proc = sqrt(2 - 2 * tr(Sigma))
#    where tr(Sigma) = sum of singular values of M.
#
#    Similarity (bounded in [0, 1]):
#      sim_Proc = 1 - m_Proc / sqrt(2)
# -----------------------------------------------------------------------------

M       <- crossprod(R_tilde, Rp_tilde)   # D x D = R_tilde' Rp_tilde
svd_M   <- svd(M)

Q_star  <- svd_M$u %*% t(svd_M$v)        # optimal rotation D x D
tr_Sigma <- sum(svd_M$d)                  # sum of singular values

m_Proc  <- sqrt(max(0, 2 - 2 * tr_Sigma))  # clamp for floating-point safety
sim_Proc <- 1 - m_Proc / sqrt(2)

# -----------------------------------------------------------------------------
# 4. Verification: direct Frobenius distance after applying Q*
#    ||R_tilde Q* - Rp_tilde||_F should equal m_Proc
# -----------------------------------------------------------------------------

direct_dist <- frob_norm(R_tilde %*% Q_star - Rp_tilde)
stopifnot(abs(direct_dist - m_Proc) < 1e-6)

# -----------------------------------------------------------------------------
# 5. Summary table
# -----------------------------------------------------------------------------

procrustes_summary <- tibble(
  Metric = c("Procrustes distance  m_Proc  (bounded in [0, sqrt(2)])",
             "Procrustes similarity  sim_Proc  (bounded in [0, 1])",
             "Sum of singular values  tr(Sigma)",
             "Direct Frobenius distance (verification)"),
  Value  = c(round(m_Proc,      6),
             round(sim_Proc,    6),
             round(tr_Sigma,    6),
             round(direct_dist, 6))
)

print(procrustes_summary)

# -----------------------------------------------------------------------------
# 6. Save results
# -----------------------------------------------------------------------------

procrustes_results <- list(
  m_Proc               = m_Proc,
  sim_Proc             = sim_Proc,
  tr_Sigma             = tr_Sigma,
  Q_star               = Q_star,
  singular_values      = svd_M$d,
  procrustes_summary   = procrustes_summary,
  N                    = N,
  D                    = D
)

saveRDS(procrustes_results, file = here("procrustes_results.rds"))
