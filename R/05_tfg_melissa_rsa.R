# =============================================================================
# 05_tfg_melissa_rsa.R  -  Representational Similarity Analysis (RSA)
# TFG - Grau en Estadística, UB/UPC
# Author: Melissa Vargas Jerez
# Description: Computes the RSA score between the real (augmented) and
#              synthetic representation matrices extracted from the frozen
#              1D-CNN. Requires feat_aug_3k.rds and feat_syn_3k.rds produced by
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
# 2. Preprocessing: column-wise mean-centring
#    Applied for consistency with the other measures, even though Pearson
#    correlation (used as the inner similarity function) is itself mean-centred.
# -----------------------------------------------------------------------------

R  <- sweep(feat_aug, 2, colMeans(feat_aug), "-")
Rp <- sweep(feat_syn, 2, colMeans(feat_syn), "-")

# -----------------------------------------------------------------------------
# 3. Core RSA functions
# -----------------------------------------------------------------------------

# -- 3.1 Representational similarity matrix (RSM) --
# Computes the N x N matrix of Pearson correlations between every pair of
# row vectors. Entry S[i,j] = Pearson correlation between patients i and j.
# cor() applied to the transpose gives correlations between columns of t(mat),
# which are the rows of mat.
compute_rsm <- function(mat) {
  cor(t(mat), method = "pearson")   # N x N
}

# -- 3.2 RSA score for one subsample --
# Extracts the lower triangular entries of each RSM (excluding the diagonal),
# vectorises them, and returns the Spearman rank correlation between the two
# vectors. Spearman is used as the outer function for robustness to extreme
# values and invariance to monotone transformations of the similarity values.
compute_rsa_subsample <- function(R_sub, Rp_sub) {
  S  <- compute_rsm(R_sub)
  Sp <- compute_rsm(Rp_sub)

  # Lower triangular indices, diagonal excluded
  idx       <- lower.tri(S, diag = FALSE)
  v_S       <- S[idx]
  v_Sp      <- Sp[idx]

  cor(v_S, v_Sp, method = "spearman")
}

# -----------------------------------------------------------------------------
# 4. Repeated subsampling
#    N = 3,000 is feasible but intensive for the full RSM. A subsample of
#    n = 1,000 patients is drawn from each dataset. The procedure is repeated
#    five times with different seeds and results are averaged, following the
#    computational note in Section 4.4.
# -----------------------------------------------------------------------------

N_SUB   <- 1000L
N_REPS  <- 5L
SEEDS   <- c(42L, 7L, 123L, 256L, 999L)

rsa_replicates <- numeric(N_REPS)

for (i in seq_len(N_REPS)) {
  set.seed(SEEDS[i])
  idx_aug <- sample(N, N_SUB)
  idx_syn <- sample(N, N_SUB)

  rsa_replicates[i] <- compute_rsa_subsample(R[idx_aug, ], Rp[idx_syn, ])
}

m_RSA <- mean(rsa_replicates)

# -----------------------------------------------------------------------------
# 5. Summary tables
# -----------------------------------------------------------------------------

# Per-replicate results
rsa_replicates_tbl <- tibble(
  Replicate = seq_len(N_REPS),
  Seed      = SEEDS,
  RSA_score = round(rsa_replicates, 6)
)

print(rsa_replicates_tbl)

# Aggregated summary
rsa_summary <- tibble(
  Metric = c("m_RSA (mean Spearman correlation)",
             "Standard deviation across replicates",
             "Minimum replicate score",
             "Maximum replicate score",
             "Subsamples per replicate (n)",
             "Number of replicates"),
  Value  = c(round(m_RSA, 6),
             round(sd(rsa_replicates), 6),
             round(min(rsa_replicates), 6),
             round(max(rsa_replicates), 6),
             N_SUB,
             N_REPS)
)

print(rsa_summary)

# -----------------------------------------------------------------------------
# 6. Save results
# -----------------------------------------------------------------------------

rsa_results <- list(
  m_RSA              = m_RSA,
  rsa_replicates     = rsa_replicates,
  rsa_replicates_tbl = rsa_replicates_tbl,
  rsa_summary        = rsa_summary,
  N_SUB              = N_SUB,
  N_REPS             = N_REPS,
  SEEDS              = SEEDS,
  N                  = N,
  D                  = D
)

saveRDS(rsa_results, file = here("rsa_results.rds"))
