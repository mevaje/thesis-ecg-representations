# =============================================================================
# 04_cka.R  --  Centered Kernel Alignment
# TFG - Grau en Estadística, UB/UPC
# Author: Melissa Vargas Jerez
# Description: Computes global and class-conditional CKA between the
#              augmented-real and synthetic representation matrices.
#              Requires feat_aug_3k.rds, feat_syn_3k.rds, class_aug_3k.rds,
#              class_syn_3k.rds produced by 01_tfg_melissa_eda.R
# =============================================================================

library(tidyverse)
library(patchwork)
library(here)

# -----------------------------------------------------------------------------
# 0. Constants (must match 01_tfg_melissa_eda.R)
# -----------------------------------------------------------------------------

N_FEATURES <- 256L
N_SUBSAMPLE <- 3000L          # patients used per dataset

CLASS_LABELS <- c(
  "0" = "SBRAD",
  "1" = "SR",
  "2" = "AFIB",
  "3" = "STACH",
  "4" = "AFLT",
  "5" = "SARRH",
  "6" = "SVTAC"
)

DATASET_COLOURS <- c(
  "Real (augmented)" = "#2E5EA8",
  "Synthetic"        = "#E07B39"
)

CLASS_COLOURS <- c(
  "SBRAD" = "#4393C3",
  "SR"    = "#92C5DE",
  "AFIB"  = "#D6604D",
  "STACH" = "#F4A582",
  "SARRH" = "#878787",
  "AFLT"  = "#4DAC26",
  "SVTAC" = "#B8E186"
)

# -----------------------------------------------------------------------------
# 1. Load objects from EDA script
# -----------------------------------------------------------------------------

feat_aug    <- readRDS(here("feat_aug_3k.rds"))
feat_syn    <- readRDS(here("feat_syn_3k.rds"))
class_aug_f <- readRDS(here("class_aug_3k.rds"))
class_syn_f <- readRDS(here("class_syn_3k.rds"))

# Confirm dimensions
stopifnot(
  nrow(feat_aug) == N_SUBSAMPLE,
  nrow(feat_syn) == N_SUBSAMPLE,
  ncol(feat_aug) == N_FEATURES,
  ncol(feat_syn) == N_FEATURES
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
# Each entry K[i,j] is the dot-product similarity between patients i and j.
gram_linear <- function(mat) {
  tcrossprod(mat)           # equivalent to mat %*% t(mat)
}

# -- 2.3 Standard HSIC estimator (Kornblith et al., 2019) --
# HSIC(K, L) = 1/(N-1)^2 * tr(K H_N L H_N)
# where H_N = I_N - (1/N) 11' is the centering matrix.
# Efficient computation avoids forming H_N explicitly:
#   tr(K H_N L H_N) = tr(Kc Lc)  where Kc = H_N K H_N
# This estimator is always non-negative and bounded in [0, 1] after
# CKA normalisation, consistent with the methodology in Section 4.3.
hsic_standard <- function(K, L) {
  n  <- nrow(K)
  stopifnot(nrow(L) == n)
  
  # Double-centre each Gram matrix: Kc = H_N K H_N
  centre_gram <- function(G) {
    row_means <- rowMeans(G)
    col_means <- colMeans(G)
    grand_mean <- mean(G)
    G - outer(row_means, rep(1, n)) - outer(rep(1, n), col_means) + grand_mean
  }
  
  Kc <- centre_gram(K)
  Lc <- centre_gram(L)
  
  sum(Kc * Lc) / (n - 1)^2
}

# -- 2.4 CKA score --
# Returns a scalar in [0, 1].
compute_cka <- function(R, Rprime) {
  R      <- centre_columns(R)
  Rprime <- centre_columns(Rprime)
  
  K <- gram_linear(R)
  L <- gram_linear(Rprime)
  
  hsic_KL <- hsic_standard(K, L)
  hsic_KK <- hsic_standard(K, K)
  hsic_LL <- hsic_standard(L, L)
  
  hsic_KL / sqrt(hsic_KK * hsic_LL)
}

# -----------------------------------------------------------------------------
# 3. Global CKA score
# -----------------------------------------------------------------------------

cka_global <- compute_cka(feat_aug, feat_syn)

cka_global_summary <- tibble(
  Metric = "CKA (global, augmented vs. synthetic)",
  Value  = round(cka_global, 6)
)

print(cka_global_summary)

# -----------------------------------------------------------------------------
# 4. Class-conditional CKA
# -----------------------------------------------------------------------------
# For each arrhythmia class, extract the submatrix of patients belonging to
# that class from both datasets and compute CKA on those submatrices.
# Note: class factors must share the same levels (defined in EDA script).

classes <- levels(class_aug_f)   # ordered character vector of class names

cka_per_class <- map_dfr(classes, function(cls) {
  
  idx_aug <- which(class_aug_f == cls)
  idx_syn <- which(class_syn_f == cls)
  
  # Skip if either class is empty in one of the datasets
  if (length(idx_aug) < 4 || length(idx_syn) < 4) {
    warning(sprintf("Class %s has too few samples, skipping.", cls))
    return(NULL)
  }
  
  # When class sizes differ, subsample the larger set so both submatrices
  # have the same N. This keeps the Gram matrix square and comparable.
  n_min <- min(length(idx_aug), length(idx_syn))
  set.seed(42)
  idx_aug <- sample(idx_aug, n_min)
  idx_syn <- sample(idx_syn, n_min)
  
  score <- compute_cka(feat_aug[idx_aug, ], feat_syn[idx_syn, ])
  
  tibble(
    class     = cls,
    n_aug     = length(idx_aug),
    n_syn     = length(idx_syn),
    n_used    = n_min,
    cka_score = score
  )
})

print(cka_per_class)

# -----------------------------------------------------------------------------
# 5. Summary table: global + per-class
# -----------------------------------------------------------------------------

cka_summary <- bind_rows(
  tibble(class = "Global", n_aug = nrow(feat_aug), n_syn = nrow(feat_syn),
         n_used = nrow(feat_aug), cka_score = cka_global),
  cka_per_class
) |>
  mutate(cka_score = round(cka_score, 4))

print(cka_summary)

# -----------------------------------------------------------------------------
# 6. Visualisation: bar chart of per-class CKA scores
# -----------------------------------------------------------------------------

# Horizontal reference line at the global score makes it easy to see which
# classes drive the overall similarity up or down.

plot_cka_class <- ggplot(cka_per_class,
                         aes(x = reorder(class, cka_score),
                             y = cka_score,
                             fill = class)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.3) +
  geom_hline(yintercept = cka_global,
             linetype   = "dashed",
             linewidth  = 0.6,
             colour     = "grey30") +
  annotate("text",
           x     = 0.6,
           y     = cka_global + 0.012,
           label = sprintf("Global CKA = %.3f", cka_global),
           hjust = 0,
           size  = 3.2,
           colour = "grey30") +
  scale_fill_manual(values = CLASS_COLOURS, guide = "none") +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, by = 0.1),
                     labels = scales::number_format(accuracy = 0.1)) +
  coord_flip() +
  labs(
    title    = "CKA similarity by arrhythmia class",
    subtitle = "Augmented-real vs. synthetic representations (penultimate layer)",
    x        = "Arrhythmia class",
    y        = "CKA score"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    plot.title         = element_text(face = "bold")
  )

ggsave(here("figures", "cka_per_class.pdf"),
       plot_cka_class, width = 8, height = 5)

# -----------------------------------------------------------------------------
# 7. Visualisation: CKA score matrix across class pairs (optional heatmap)
# -----------------------------------------------------------------------------
# Computes CKA between every pair of classes (aug class i vs syn class j).
# The diagonal gives the same-class scores; off-diagonal shows cross-class
# similarity and can reveal whether the synthetic data mixes classes.

class_pairs <- expand.grid(cls_aug = classes, cls_syn = classes,
                           stringsAsFactors = FALSE)

cka_matrix_vals <- map2_dbl(class_pairs$cls_aug, class_pairs$cls_syn,
                            function(ca, cs) {
                              idx_a <- which(class_aug_f == ca)
                              idx_s <- which(class_syn_f == cs)
                              if (length(idx_a) < 4 || length(idx_s) < 4) return(NA_real_)
                              n_min <- min(length(idx_a), length(idx_s))
                              set.seed(42)
                              compute_cka(feat_aug[sample(idx_a, n_min), ],
                                          feat_syn[sample(idx_s, n_min), ])
                            }
)

cka_matrix_df <- class_pairs |>
  mutate(cka_score = round(cka_matrix_vals, 3))

plot_cka_heatmap <- ggplot(cka_matrix_df,
                           aes(x = cls_syn, y = cls_aug, fill = cka_score)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", cka_score)),
            size = 3.2, colour = "white") +
  scale_fill_gradient2(
    low      = "#053061",
    mid      = "#4393C3",
    high     = "#FFFFFF",
    midpoint = 0.5,
    limits   = c(0, 1),
    name     = "CKA"
  ) +
  labs(
    title    = "Cross-class CKA matrix",
    subtitle = "Rows: augmented-real classes  |  Columns: synthetic classes",
    x        = "Synthetic class",
    y        = "Augmented-real class"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 30, hjust = 1),
    plot.title   = element_text(face = "bold"),
    panel.grid   = element_blank()
  )

ggsave(here("figures", "cka_class_heatmap.pdf"),
       plot_cka_heatmap, width = 7, height = 6)

# -----------------------------------------------------------------------------
# 8. Save results for downstream scripts
# -----------------------------------------------------------------------------

saveRDS(cka_global,         here("cka_global.rds"))
saveRDS(cka_global_summary, here("cka_global_summary.rds"))
saveRDS(cka_per_class,      here("cka_per_class.rds"))
saveRDS(cka_matrix_df,      here("cka_matrix_df.rds"))

message("Results saved. Figures written to figures/")