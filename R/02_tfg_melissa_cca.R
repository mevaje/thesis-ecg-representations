# =============================================================================
# Canonical Correlation Analysis (CCA)
# TFG - Grau en Estadística, UB/UPC
# Author: Melissa Vargas Jerez
# Description: Computes CCA between the real (augmented) and synthetic
#              representation matrices extracted from the frozen 1D-CNN.
#              Follows what described in the Methodology Section of the thesis.
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
# 3. Covariance matrices
# Sigma_RR  = (1/N) R'R
# Sigma_RpRp = (1/N) R''R'
# Sigma_RRp = (1/N) R'R'
# -----------------------------------------------------------------------------

Sigma_RR   <- crossprod(R)  / N          # D x D
Sigma_RpRp <- crossprod(Rp) / N          # D x D
Sigma_RRp  <- crossprod(R, Rp) / N       # D x D

# -----------------------------------------------------------------------------
# 4. Symmetric matrix square root inverse via eigendecomposition
#    For a PSD matrix A = Q diag(lambda) Q',
#    A^{-1/2} = Q diag(1/sqrt(lambda)) Q'
#    A small ridge is added for numerical stability before inversion.
# -----------------------------------------------------------------------------

mat_inv_sqrt <- function(A, ridge = 1e-8) {
  eig <- eigen(A, symmetric = TRUE)
  # Clamp tiny/negative eigenvalues produced by floating-point error
  d   <- pmax(eig$values, 0)
  d_inv_sqrt <- ifelse(d > ridge, 1 / sqrt(d), 0)
  eig$vectors %*% diag(d_inv_sqrt) %*% t(eig$vectors)
}

Sigma_RR_inv_sqrt   <- mat_inv_sqrt(Sigma_RR)
Sigma_RpRp_inv_sqrt <- mat_inv_sqrt(Sigma_RpRp)

# -----------------------------------------------------------------------------
# 5. Whitened cross-covariance matrix and its SVD
#    M = Sigma_RR^{-1/2} Sigma_RRp Sigma_RpRp^{-1/2}
#    Singular values of M are the canonical correlations rho_1 >= ... >= rho_D
#    Canonical weight vectors:
#      a_i = Sigma_RR^{-1/2}   u_i
#      b_i = Sigma_RpRp^{-1/2} v_i
# -----------------------------------------------------------------------------

M   <- Sigma_RR_inv_sqrt %*% Sigma_RRp %*% Sigma_RpRp_inv_sqrt

svd_M <- svd(M)

canonical_correlations <- svd_M$d          # length D, ordered desc
A <- Sigma_RR_inv_sqrt   %*% svd_M$u       # D x D canonical weights for R
B <- Sigma_RpRp_inv_sqrt %*% svd_M$v       # D x D canonical weights for Rp

# Clamp numerical noise, correlations must lie in [0, 1]
canonical_correlations <- pmin(pmax(canonical_correlations, 0), 1)

# -----------------------------------------------------------------------------
# 6. Score
#    m_CCA = (1/D) * sum_{i=1}^{D} rho_i
# -----------------------------------------------------------------------------

m_CCA <- mean(canonical_correlations)

cca_summary <- tibble(
  Metric = c("m_CCA (mean canonical correlation)",
             "rho_1 (maximum)",
             "rho_D (minimum)",
             "Canonical dimensions (D)"),
  Value  = c(round(m_CCA, 6),
             round(canonical_correlations[1], 6),
             round(canonical_correlations[D], 6),
             D)
)

print(cca_summary)
# -----------------------------------------------------------------------------
# 7. Diagnostic: distribution of canonical correlations
# -----------------------------------------------------------------------------


print(round(quantile(canonical_correlations,
                     probs = c(0, 0.25, 0.50, 0.75, 1)), 6))

# -----------------------------------------------------------------------------
# 8. Save results
# -----------------------------------------------------------------------------

cca_results <- list(
  canonical_correlations = canonical_correlations,  # vector of length D
  m_CCA                  = m_CCA,
  A                      = A,                       # canonical weight matrix for R
  B                      = B,                       # canonical weight matrix for Rp
  M                      = M,                       # whitened cross-covariance
  N                      = N,
  D                      = D
)

saveRDS(cca_results, file = here("cca_results.rds"))

