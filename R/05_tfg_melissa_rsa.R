# =============================================================================
# 05_tfg_melissa_rsa.R  -  Representational Similarity Analysis (RSA)
# TFG - Grau en Estadística, UB/UPC
# Author: Melissa Vargas Jerez
# Description: Computes RSA along the four-axis comparison framework
#              (real vs. augmented, augmented vs. synthetic; each overall
#              and per arrhythmia class), as described in the Methodology
#              chapter of the thesis.
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
# 2. Preprocessing: column-wise mean-centring
#    Applied for consistency with the other measures, even though Pearson
#    correlation (used as the inner similarity function) is itself mean-centred.
# -----------------------------------------------------------------------------

centre_columns <- function(mat) {
  sweep(mat, 2, colMeans(mat), "-")
}

# -----------------------------------------------------------------------------
# 3. Core RSA functions
# -----------------------------------------------------------------------------

# -- 3.1 Representational similarity matrix (RSM) --
# Computes the N x N matrix of Pearson correlations between every pair of
# row vectors. cor() applied to the transpose gives correlations between
# columns of t(mat), which are the rows of mat.
compute_rsm <- function(mat) {
  cor(t(mat), method = "pearson")
}

# -- 3.2 RSA score between two (already centred) matrices --
# Extracts the lower triangular entries of each RSM (excluding the diagonal),
# vectorises them, and returns the Spearman rank correlation between the two
# vectors. Requires both matrices to have the same N, since v(S) and v(S')
# must be the same length to correlate.
compute_rsa <- function(mat_X, mat_Y) {
  stopifnot(nrow(mat_X) == nrow(mat_Y))
  
  S  <- compute_rsm(mat_X)
  Sp <- compute_rsm(mat_Y)
  
  idx  <- lower.tri(S, diag = FALSE)
  v_S  <- S[idx]
  v_Sp <- Sp[idx]
  
  cor(v_S, v_Sp, method = "spearman")
}

# -----------------------------------------------------------------------------
# 4. Axis 1 and 3: overall comparisons, with repeated subsampling
#    At N = 3,000, the full RSM is feasible but intensive. A subsample of
#    n = 1,000 patients is drawn from each dataset, repeated five times with
#    different seeds, and results are averaged.
# -----------------------------------------------------------------------------

N_SUB  <- 1000L
N_REPS <- 5L
SEEDS  <- c(42L, 7L, 123L, 256L, 999L)

compute_rsa_overall <- function(feat_X, feat_Y,
                                n_sub = N_SUB, n_reps = N_REPS, seeds = SEEDS) {
  X <- centre_columns(feat_X)
  Y <- centre_columns(feat_Y)
  
  N_X <- nrow(X)
  N_Y <- nrow(Y)
  
  replicates <- numeric(n_reps)
  for (i in seq_len(n_reps)) {
    set.seed(seeds[i])
    idx_X <- sample(N_X, n_sub)
    idx_Y <- sample(N_Y, n_sub)
    replicates[i] <- compute_rsa(X[idx_X, ], Y[idx_Y, ])
  }
  
  list(
    m_RSA      = mean(replicates),
    replicates = replicates,
    n_sub      = n_sub,
    n_reps     = n_reps,
    seeds      = seeds
  )
}

rsa_axis1_overall <- compute_rsa_overall(feat_real, feat_aug)
rsa_axis3_overall <- compute_rsa_overall(feat_aug,  feat_syn)

# -----------------------------------------------------------------------------
# 5. Axis 2 and 4: per-class comparisons, no subsampling
#    Per-class groups are small enough (well under the n = 1,000 threshold
#    that motivates subsampling at the overall level) that the full class is
#    used directly, with no repeated subsampling. As with the other measures,
#    classes are matched to the smaller of the two available counts via
#    random subsampling, since v(S) and v(S') must be the same length to
#    correlate.
# -----------------------------------------------------------------------------

compute_rsa_per_class <- function(feat_X, class_X, feat_Y, class_Y,
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
    
    n_use <- min(length(idx_X), length(idx_Y))
    set.seed(seed)
    idx_X_use <- sample(idx_X, n_use)
    idx_Y_use <- sample(idx_Y, n_use)
    
    X_c <- centre_columns(feat_X[idx_X_use, , drop = FALSE])
    Y_c <- centre_columns(feat_Y[idx_Y_use, , drop = FALSE])
    
    score <- compute_rsa(X_c, Y_c)
    
    tibble(Class = cl, N_used = n_use, RSA = score)
  })
  results
}

rsa_axis2_per_class <- compute_rsa_per_class(feat_real, class_real,
                                             feat_aug,  class_aug)

rsa_axis4_per_class <- compute_rsa_per_class(feat_aug, class_aug,
                                             feat_syn, class_syn)

# -----------------------------------------------------------------------------
# 6. Summary tables
# -----------------------------------------------------------------------------

# Overall scores across both axes (1 and 3), with replicate spread
rsa_overall_summary <- tibble(
  Axis    = c("1: Real vs. Augmented", "3: Augmented vs. Synthetic"),
  m_RSA   = round(c(rsa_axis1_overall$m_RSA, rsa_axis3_overall$m_RSA), 6),
  SD      = round(c(sd(rsa_axis1_overall$replicates),
                    sd(rsa_axis3_overall$replicates)), 6),
  Min     = round(c(min(rsa_axis1_overall$replicates),
                    min(rsa_axis3_overall$replicates)), 6),
  Max     = round(c(max(rsa_axis1_overall$replicates),
                    max(rsa_axis3_overall$replicates)), 6)
)

print(rsa_overall_summary)

# Per-class scores, both axes combined into one tidy table
rsa_per_class_summary <- bind_rows(
  rsa_axis2_per_class |> mutate(Axis = "2: Real vs. Augmented", .before = 1),
  rsa_axis4_per_class |> mutate(Axis = "4: Augmented vs. Synthetic", .before = 1)
) |>
  mutate(RSA = round(RSA, 6))

print(rsa_per_class_summary)

# -----------------------------------------------------------------------------
# 7. Save results
# -----------------------------------------------------------------------------

rsa_results <- list(
  axis1_overall         = rsa_axis1_overall,
  axis2_per_class       = rsa_axis2_per_class,
  axis3_overall         = rsa_axis3_overall,
  axis4_per_class       = rsa_axis4_per_class,
  rsa_overall_summary   = rsa_overall_summary,
  rsa_per_class_summary = rsa_per_class_summary,
  N_SUB                 = N_SUB,
  N_REPS                = N_REPS,
  SEEDS                 = SEEDS,
  D                     = D
)

saveRDS(rsa_results, file = here("rsa_results.rds"))