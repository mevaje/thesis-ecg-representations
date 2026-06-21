# =============================================================================
# 04_tfg_melissa_cka.R  --  Centered Kernel Alignment
# TFG - Grau en Estadística, UB/UPC
# Author: Melissa Vargas Jerez
# Description: Computes CKA along the four-axis comparison framework
#              (real vs. augmented, augmented vs. synthetic; each overall
#              and per arrhythmia class), as described in the Methodology
#              chapter of the thesis.
#              Requires feat_real_3k.rds, feat_aug_3k.rds, feat_syn_3k.rds,
#              class_real_3k.rds, class_aug_3k.rds, class_syn_3k.rds
#              produced by 01_tfg_melissa_ear.R
# =============================================================================

library(tidyverse)
library(here)

# -----------------------------------------------------------------------------
# 0. Constants (must match 01_tfg_melissa_ear.R)
# -----------------------------------------------------------------------------

N_FEATURES  <- 256L
N_SUBSAMPLE <- 3000L

CLASS_COLOURS <- c(
  "SBRAD" = "#E69F00",
  "SR"    = "#0072B2",
  "AFIB"  = "#D55E00",
  "STACH" = "#009E73",
  "AFLT"  = "#CC79A7",
  "SARRH" = "#F0E442",
  "SVTAC" = "#000000"
)

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------
# All three matrices are the N = 3,000 (stratified for real) subsamples
# produced by 01_tfg_melissa_ear.R, so that all six similarity measures
# operate on identical observations.

feat_real  <- readRDS(here("feat_real_3k.rds"))
feat_aug   <- readRDS(here("feat_aug_3k.rds"))
feat_syn   <- readRDS(here("feat_syn_3k.rds"))

class_real <- readRDS(here("class_real_3k.rds"))
class_aug  <- readRDS(here("class_aug_3k.rds"))
class_syn  <- readRDS(here("class_syn_3k.rds"))

stopifnot(
  ncol(feat_real) == N_FEATURES,
  ncol(feat_aug)  == N_FEATURES,
  ncol(feat_syn)  == N_FEATURES
)

# -----------------------------------------------------------------------------
# 2. Core CKA functions
# -----------------------------------------------------------------------------

# -- 2.1 Column-wise mean-centring --
centre_columns <- function(mat) {
  scale(mat, center = TRUE, scale = FALSE)
}

# -- 2.2 Linear Gram matrix --
# For a centred N x D matrix R, returns the N x N matrix K = R %*% t(R).
gram_linear <- function(mat) {
  tcrossprod(mat)
}

# -- 2.3 Standard HSIC estimator (Kornblith et al., 2019) --
# HSIC(K, L) = 1/(N-1)^2 * tr(K H_N L H_N)
# Efficient computation avoids forming H_N explicitly:
#   tr(K H_N L H_N) = tr(Kc Lc)  where Kc = H_N K H_N
hsic_standard <- function(K, L) {
  n <- nrow(K)
  stopifnot(nrow(L) == n)
  
  centre_gram <- function(G) {
    row_means  <- rowMeans(G)
    col_means  <- colMeans(G)
    grand_mean <- mean(G)
    G - outer(row_means, rep(1, n)) - outer(rep(1, n), col_means) + grand_mean
  }
  
  Kc <- centre_gram(K)
  Lc <- centre_gram(L)
  
  sum(Kc * Lc) / (n - 1)^2
}

# -- 2.4 CKA score --
# Returns a scalar in [0, 1]. Requires feat_X and feat_Y to have equal N,
# since the Gram matrices K and L must be the same size for tr(Kc Lc).
compute_cka <- function(feat_X, feat_Y) {
  stopifnot(nrow(feat_X) == nrow(feat_Y))
  
  X <- centre_columns(feat_X)
  Y <- centre_columns(feat_Y)
  
  K <- gram_linear(X)
  L <- gram_linear(Y)
  
  hsic_KL <- hsic_standard(K, L)
  hsic_KK <- hsic_standard(K, K)
  hsic_LL <- hsic_standard(L, L)
  
  hsic_KL / sqrt(hsic_KK * hsic_LL)
}

# -----------------------------------------------------------------------------
# 3. Axis 1 and 3: overall comparisons
#    Axis 1: Real (reference) vs. Augmented
#    Axis 3: Augmented vs. Synthetic
#    feat_real, feat_aug, and feat_syn all have N = 3,000, so no N-matching
#    is needed for the overall comparison.
# -----------------------------------------------------------------------------

cka_axis1_overall <- compute_cka(feat_real, feat_aug)
cka_axis3_overall <- compute_cka(feat_aug,  feat_syn)

# -----------------------------------------------------------------------------
# 4. Axis 2 and 4: per-class comparisons
#    For each arrhythmia class, CKA is computed only on the patients of that
#    class within the relevant pair of matrices. Class sizes generally
#    differ between datasets (e.g. AFIB is rare in the real reference data
#    but balanced in the augmented data), so each class is matched to the
#    smaller of the two counts via random subsampling, consistent with the
#    approach used for CCA and RV.
# -----------------------------------------------------------------------------

compute_cka_per_class <- function(feat_X, class_X, feat_Y, class_Y,
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
    
    score <- compute_cka(feat_X[idx_X_use, ], feat_Y[idx_Y_use, ])
    
    tibble(Class = cl, N_used = n_use, CKA = score)
  })
  results
}

cka_axis2_per_class <- compute_cka_per_class(feat_real, class_real,
                                             feat_aug,  class_aug)

cka_axis4_per_class <- compute_cka_per_class(feat_aug, class_aug,
                                             feat_syn, class_syn)

# -----------------------------------------------------------------------------
# 5. Summary tables
# -----------------------------------------------------------------------------

# Overall scores across both axes (1 and 3)
cka_overall_summary <- tibble(
  Axis = c("1: Real vs. Augmented", "3: Augmented vs. Synthetic"),
  CKA  = round(c(cka_axis1_overall, cka_axis3_overall), 6)
)

print(cka_overall_summary)

# Per-class scores, both axes combined into one tidy table
cka_per_class_summary <- bind_rows(
  cka_axis2_per_class |> mutate(Axis = "2: Real vs. Augmented", .before = 1),
  cka_axis4_per_class |> mutate(Axis = "4: Augmented vs. Synthetic", .before = 1)
) |>
  mutate(CKA = round(CKA, 6))

print(cka_per_class_summary)

# -----------------------------------------------------------------------------
# 6. Visualisation: per-class CKA, both axes side by side
# -----------------------------------------------------------------------------

plot_cka_per_class <- ggplot(cka_per_class_summary,
                             aes(x = reorder(Class, CKA), y = CKA, fill = Class)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.3) +
  geom_hline(data = cka_overall_summary,
             aes(yintercept = CKA),
             linetype  = "dashed",
             linewidth = 0.5,
             colour    = "grey30") +
  facet_wrap(~Axis, ncol = 1, scales = "free_x") +
  scale_fill_manual(values = CLASS_COLOURS, guide = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  coord_flip() +
  labs(
    title    = "CKA similarity by arrhythmia class",
    subtitle = "Dashed line marks the overall (non class-conditional) score for that axis",
    x        = "Arrhythmia class",
    y        = "CKA score"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    plot.title          = element_text(face = "bold"),
    strip.text          = element_text(face = "bold")
  )

ggsave(here("figures", "cka_per_class.pdf"),
       plot_cka_per_class, width = 8, height = 9)

# -----------------------------------------------------------------------------
# 7. Visualisation: cross-class CKA heatmap (axis 4: augmented vs. synthetic)
# -----------------------------------------------------------------------------
# Computes CKA between every pair of classes (aug class i vs syn class j).
# The diagonal gives the same-class scores already reported in Section 4;
# off-diagonal entries reveal whether the synthetic data mixes classes,
# i.e. whether one synthetic class resembles a *different* augmented class.

build_cross_class_heatmap <- function(feat_X, class_X, feat_Y, class_Y,
                                      min_n = 4L, seed = 42) {
  classes <- levels(class_X)
  class_pairs <- expand.grid(cls_X = classes, cls_Y = classes,
                             stringsAsFactors = FALSE)
  
  scores <- map2_dbl(class_pairs$cls_X, class_pairs$cls_Y, function(cx, cy) {
    idx_X <- which(class_X == cx)
    idx_Y <- which(class_Y == cy)
    if (length(idx_X) < min_n || length(idx_Y) < min_n) return(NA_real_)
    n_use <- min(length(idx_X), length(idx_Y))
    set.seed(seed)
    compute_cka(feat_X[sample(idx_X, n_use), ], feat_Y[sample(idx_Y, n_use), ])
  })
  
  class_pairs |> mutate(CKA = round(scores, 3))
}

cka_heatmap_axis4 <- build_cross_class_heatmap(feat_aug, class_aug,
                                               feat_syn, class_syn)

plot_cka_heatmap <- ggplot(cka_heatmap_axis4,
                           aes(x = cls_Y, y = cls_X, fill = CKA)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", CKA)), size = 3.2, colour = "white") +
  scale_fill_gradient2(
    low      = "#053061",
    mid      = "#4393C3",
    high     = "#FFFFFF",
    midpoint = 0.5,
    limits   = c(0, 1),
    name     = "CKA"
  ) +
  labs(
    title    = "Cross-class CKA matrix (Axis 4: Augmented vs. Synthetic)",
    subtitle = "Rows: augmented classes  |  Columns: synthetic classes",
    x        = "Synthetic class",
    y        = "Augmented class"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title  = element_text(face = "bold"),
    panel.grid  = element_blank()
  )

ggsave(here("figures", "cka_class_heatmap.pdf"),
       plot_cka_heatmap, width = 7, height = 6)

# -----------------------------------------------------------------------------
# 8. Save results for downstream scripts
# -----------------------------------------------------------------------------

cka_results <- list(
  axis1_overall         = cka_axis1_overall,
  axis2_per_class       = cka_axis2_per_class,
  axis3_overall         = cka_axis3_overall,
  axis4_per_class       = cka_axis4_per_class,
  cka_overall_summary   = cka_overall_summary,
  cka_per_class_summary = cka_per_class_summary,
  cka_heatmap_axis4     = cka_heatmap_axis4
)

saveRDS(cka_results, file = here("cka_results.rds"))