# =============================================================================
# Exploratory Analysis of Representations
# TFG - Grau en Estadística, UB/UPC
# Author: Melissa Vargas Jerez
# Description: Descriptive analysis of the three representation matrices
#              extracted from the frozen 1D-CNN (real, augmented, synthetic)
# =============================================================================


# libraries
library(tidyverse)
library(patchwork)
library(scales)
library(here)

# Output directories
dir.create(here("figures"), showWarnings = FALSE)
dir.create(here("tables"),  showWarnings = FALSE)

# Shared constants
N_FEATURES    <- 256L
CLASS_COL     <- 257L
N_CLASSES     <- 7L

CLASS_LABELS  <- c(
  "0" = "SBRAD",
  "1" = "SR",
  "2" = "AFIB",
  "3" = "STACH",
  "4" = "AFLT",
  "5" = "SARRH",
  "6" = "SVTAC"
)

# Plots colour palette
DATASET_COLOURS <- c(
  "Real (augmented)"  = "#2E5EA8",
  "Synthetic"         = "#E07B39",
  "Real (reference)"  = "#6DAE6A"
)

CLASS_COLOURS <- c(
  "SBRAD" = "#E69F00",  # orange
  "SR"    = "#0072B2",  # blue
  "AFIB"  = "#D55E00",  # vermillion red
  "STACH" = "#009E73",  # bluish green
  "AFLT"  = "#CC79A7",  # reddish purple / magenta
  "SARRH" = "#F0E442",  # yellow
  "SVTAC" = "#000000"   # black
)

# -----------------------------------------------------------------------------
# 1. Load data
# -----------------------------------------------------------------------------

# Link to download .RData object is in the README section of this GitHub repository
# Assuming RData object and script match current work directory

load(here("data_loaded.RData"))

# -----------------------------------------------------------------------------
# 2. Similarity metrics require a pure numeric matrix, class label excluded
# -----------------------------------------------------------------------------

# Feature matrices: pure numeric N x 256 matrices
feat_aug  <- as.matrix(df_aug[,  seq_len(N_FEATURES)])
feat_syn  <- as.matrix(df_syn[,  seq_len(N_FEATURES)])
feat_real <- as.matrix(df_real[, seq_len(N_FEATURES)])

# Class label vectors (integer 0-6)
# Column V257 is the last column and contains the arrhythmia class (0-6)
class_aug  <- df_aug$V257
class_syn  <- df_syn$V257
class_real <- df_real$V257

# Class label vectors as ordered factor with named levels
class_aug_f  <- factor(CLASS_LABELS[as.character(class_aug)],
                       levels = unname(CLASS_LABELS))
class_syn_f  <- factor(CLASS_LABELS[as.character(class_syn)],
                       levels = unname(CLASS_LABELS))
class_real_f <- factor(CLASS_LABELS[as.character(class_real)],
                       levels = unname(CLASS_LABELS))

# -----------------------------------------------------------------------------
# 3. Dimension check and class distribution
# -----------------------------------------------------------------------------

# Confirms all three matrices have the expected 256-feature structure
dim_summary <- data.frame(
  Dataset    = c("Real (augmented)", "Synthetic", "Real (reference)"),
  Patients   = c(nrow(feat_aug), nrow(feat_syn), nrow(feat_real)),
  Features   = c(ncol(feat_aug), ncol(feat_syn), ncol(feat_real))
)
print(dim_summary, row.names = FALSE)

# Class distribution per dataset to verify balance across arrhythmia categories
class_dist <- rbind(
  table(class_aug_f),
  table(class_syn_f),
  table(class_real_f)
)
rownames(class_dist) <- c("Real (augmented)", "Synthetic", "Real (reference)")
print(class_dist)

# -----------------------------------------------------------------------------
# 4. Distribution of representation norms by class
# -----------------------------------------------------------------------------

compute_norms_df <- function(feat_matrix, class_factor, dataset_name) {
  data.frame(
    norm    = sqrt(rowSums(feat_matrix^2)),
    class   = class_factor,
    dataset = dataset_name
  )
}

norms_aug  <- compute_norms_df(feat_aug,  class_aug_f,  "Real (augmented)")
norms_syn  <- compute_norms_df(feat_syn,  class_syn_f,  "Synthetic")
norms_real <- compute_norms_df(feat_real, class_real_f, "Real (reference)")

norms_all <- bind_rows(norms_aug, norms_syn, norms_real) |>
  mutate(dataset = factor(dataset,
                          levels = c("Real (augmented)",
                                     "Synthetic",
                                     "Real (reference)")))

plot_class_boxplot <- ggplot(norms_all,
                             aes(x = class, y = norm, fill = class)) +
  geom_boxplot(outlier.size = 0.6, outlier.alpha = 0.4, linewidth = 0.3) +
  scale_fill_manual(values = CLASS_COLOURS) +
  facet_wrap(~dataset, ncol = 1, scales = "free_y") +
  labs(
    title = "Distribution of patient representation norms by arrhythmia class",
    x     = "Arrhythmia class",
    y     = "Euclidean norm of representation vector"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title      = element_text(face = "bold"),
    strip.text      = element_text(face = "bold")
  )

ggsave(here("figures", "class_boxplot.pdf"),
       plot_class_boxplot, width = 9, height = 10)

# -----------------------------------------------------------------------------
# 5. Activation sparsity analysis
# -----------------------------------------------------------------------------
# ReLU activations: a value of zero means that neuron was inactive for
# that patient. Sparsity = proportion of zero-valued entries.

compute_sparsity <- function(feat_matrix, class_vector, dataset_name) {
  
  # Overall sparsity per neuron (proportion of patients where neuron = 0)
  sparsity_per_neuron <- colMeans(feat_matrix == 0)
  
  # Sparsity per class
  classes <- sort(unique(class_vector))
  class_sparsity <- sapply(classes, function(cl) {
    mean(feat_matrix[class_vector == cl, ] == 0)
  })
  names(class_sparsity) <- CLASS_LABELS[as.character(classes)]
  
  list(
    dataset           = dataset_name,
    overall_sparsity  = mean(feat_matrix == 0),
    per_neuron        = sparsity_per_neuron,
    per_class         = class_sparsity
  )
}

sparsity_aug  <- compute_sparsity(feat_aug,  class_aug,  "Real (augmented)")
sparsity_syn  <- compute_sparsity(feat_syn,  class_syn,  "Synthetic")
sparsity_real <- compute_sparsity(feat_real, class_real, "Real (reference)")

# Overall sparsity table
sparsity_summary <- data.frame(
  Dataset          = c("Real (augmented)", "Synthetic", "Real (reference)"),
  Overall_Sparsity = paste0(round(c(sparsity_aug$overall_sparsity,
                                    sparsity_syn$overall_sparsity,
                                    sparsity_real$overall_sparsity) * 100, 1), "%")
)
print(sparsity_summary, row.names = FALSE)

# -----------------------------------------------------------------------------
# 6. PCA visualisation of the representation space
# -----------------------------------------------------------------------------

run_pca <- function(feat_matrix, dataset_name) {
  pca_result <- prcomp(feat_matrix, center = TRUE, scale. = FALSE)
  variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  list(
    pca               = pca_result,
    var_explained     = variance_explained,
    cumvar_explained  = cumsum(variance_explained),
    dataset           = dataset_name
  )
}

pca_aug  <- run_pca(feat_aug,  "Real (augmented)")
pca_syn  <- run_pca(feat_syn,  "Synthetic")
pca_real <- run_pca(feat_real, "Real (reference)")


# Scree plot: variance explained by first 30 components
scree_df <- data.frame(
  component = rep(seq_len(30), 3),
  variance  = c(pca_aug$var_explained[1:30],
                pca_syn$var_explained[1:30],
                pca_real$var_explained[1:30]),
  dataset   = rep(c("Real (augmented)", "Synthetic", "Real (reference)"),
                  each = 30)
) |>
  mutate(dataset = factor(dataset,
                          levels = c("Real (augmented)",
                                     "Synthetic",
                                     "Real (reference)")))

plot_scree <- ggplot(scree_df, aes(x = component, y = variance,
                                   colour = dataset)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = DATASET_COLOURS) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(
    title  = "Scree plot: variance explained by the first 30 principal components",
    x      = "Principal component",
    y      = "Proportion of variance explained",
    colour = "Dataset"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title      = element_text(face = "bold"))

ggsave(here("figures", "pca_scree.pdf"),
       plot_scree, width = 9, height = 5)


set.seed(42)
MAX_PLOT_POINTS <- 3000L

# Simple random subsample (suitable for augmented and synthetic)
build_pca_scatter_df <- function(pca_obj, class_factor, dataset_name) {
  n   <- nrow(pca_obj$pca$x)
  idx <- if (n > MAX_PLOT_POINTS) sample(n, MAX_PLOT_POINTS) else seq_len(n)
  data.frame(
    pc1     = pca_obj$pca$x[idx, 1],
    pc2     = pca_obj$pca$x[idx, 2],
    class   = class_factor[idx],
    dataset = dataset_name
  )
}

# Stratified subsample (required for the real reference dataset)
build_pca_scatter_df_stratified <- function(pca_obj, class_factor,
                                            dataset_name,
                                            total_n = MAX_PLOT_POINTS,
                                            min_per_class = 30L) {
  n <- nrow(pca_obj$pca$x)
  if (n <= total_n) {
    idx <- seq_len(n)
  } else {
    classes      <- levels(class_factor)
    class_counts <- table(class_factor)
    # Proportional allocation, then enforce a minimum floor per class
    proportional <- round(total_n * class_counts / sum(class_counts))
    n_per_class  <- pmax(proportional, min_per_class)
    # Cap at the available count for any class smaller than the floor
    n_per_class  <- pmin(n_per_class, class_counts)
    
    idx <- unlist(lapply(classes, function(cl) {
      rows <- which(class_factor == cl)
      sample(rows, n_per_class[cl])
    }))
  }
  data.frame(
    pc1     = pca_obj$pca$x[idx, 1],
    pc2     = pca_obj$pca$x[idx, 2],
    class   = class_factor[idx],
    dataset = dataset_name
  )
}

# PCA sign is arbitrary: align PC1 and PC2 of synthetic and real (reference)
# to real (augmented) by comparing the leading rotation vectors.
# A negative dot product between corresponding eigenvectors means the axis
# points in the opposite direction, so we flip the scores for that component.
# This is purely cosmetic for the scatter plot; no metric is affected.
align_pca_signs <- function(pca_ref, pca_target) {
  for (j in 1:2) {
    if (sum(pca_ref$pca$rotation[, j] * pca_target$pca$rotation[, j]) < 0) {
      pca_target$pca$x[, j]        <- -pca_target$pca$x[, j]
      pca_target$pca$rotation[, j] <- -pca_target$pca$rotation[, j]
    }
  }
  pca_target
}

pca_syn  <- align_pca_signs(pca_aug, pca_syn)
pca_real <- align_pca_signs(pca_aug, pca_real)

scatter_aug  <- build_pca_scatter_df(pca_aug,  class_aug_f,  "Real (augmented)")
scatter_syn  <- build_pca_scatter_df(pca_syn,  class_syn_f,  "Synthetic")
scatter_real <- build_pca_scatter_df_stratified(pca_real, class_real_f,
                                                "Real (reference)")

plot_pca_scatter <- function(scatter_df, pca_obj, dataset_name) {
  pct1 <- round(pca_obj$var_explained[1] * 100, 1)
  pct2 <- round(pca_obj$var_explained[2] * 100, 1)
  ggplot(scatter_df, aes(x = pc1, y = pc2, colour = class)) +
    geom_point(size = 0.6, alpha = 0.5) +
    scale_colour_manual(values = CLASS_COLOURS,
                        name   = "Arrhythmia class") +
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1))) +
    labs(
      title = dataset_name,
      x     = paste0("PC1 (", pct1, "%)"),
      y     = paste0("PC2 (", pct2, "%)")
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title      = element_text(face = "bold"),
          legend.position = "bottom",
          legend.key.size = unit(0.4, "cm"))
}

p_aug  <- plot_pca_scatter(scatter_aug,  pca_aug,  "Real (augmented)")
p_syn  <- plot_pca_scatter(scatter_syn,  pca_syn,  "Synthetic")
p_real <- plot_pca_scatter(scatter_real, pca_real, "Real (reference)")

plot_pca_combined <- p_aug + p_syn + p_real +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(here("figures", "pca_scatter.pdf"),
       plot_pca_combined, width = 15, height = 6)

# -----------------------------------------------------------------------------
# 7. Descriptive statistics of representation geometry
# -----------------------------------------------------------------------------
# Implements magnitude, concentricity, uniformity, tolerance, intrinsic
# dimensionality and sparsity as defined in Klabunde et al. (2023).
# These metrics characterise the geometric structure of each representation
# matrix independently, before any pairwise similarity is computed.

# -- Helper: cosine similarity between two vectors --
cosine_similarity <- function(u, v) {
  sum(u * v) / (sqrt(sum(u^2)) * sqrt(sum(v^2)))
}

# -- 7.1 Magnitude --
# Mean Euclidean length of the mean representation vector.
# Large values indicate strong average activation magnitude.
compute_magnitude <- function(feat_matrix) {
  mean_vec <- colMeans(feat_matrix)
  sqrt(sum(mean_vec^2))
}

# -- 7.2 Concentricity --
# Mean cosine similarity of each patient to the dataset mean representation.
# Values near 1 indicate all patients point in a similar direction.
compute_concentricity <- function(feat_matrix) {
  mean_vec      <- colMeans(feat_matrix)
  instance_sims <- apply(feat_matrix, 1, cosine_similarity, v = mean_vec)
  mean(instance_sims)
}

# -- 7.3 Uniformity --
# How evenly representations are distributed across the unit hypersphere.
# More negative values indicate more uniform use of representational space.
# Computed on a random subsample due to O(N^2) cost.
compute_uniformity <- function(feat_matrix, t = 2, n_sample = 2000L) {
  n   <- nrow(feat_matrix)
  idx <- if (n > n_sample) sample(n, n_sample) else seq_len(n)
  sub <- feat_matrix[idx, ]
  sq_norms    <- rowSums(sub^2)
  sq_dist_mat <- outer(sq_norms, sq_norms, "+") - 2 * tcrossprod(sub)
  sq_dist_mat <- pmax(sq_dist_mat, 0)
  log(mean(exp(-t * sq_dist_mat)))
}

# -- 7.4 Tolerance --
# Mean inner product between unit-normalised representations of the same class.
# Higher values indicate the network groups same-class patients more tightly.
compute_tolerance <- function(feat_matrix, class_vector) {
  row_norms  <- sqrt(rowSums(feat_matrix^2))
  row_norms  <- pmax(row_norms, 1e-10)
  feat_norm  <- feat_matrix / row_norms
  classes    <- sort(unique(class_vector))
  tol_values <- numeric(length(classes))
  for (k in seq_along(classes)) {
    sub_k         <- feat_norm[class_vector == classes[k], , drop = FALSE]
    gram_k        <- tcrossprod(sub_k)
    tol_values[k] <- mean(gram_k)
  }
  class_sizes <- table(class_vector)[as.character(classes)]
  weighted.mean(tol_values, w = as.numeric(class_sizes))
}

# -- 7.5 Intrinsic dimensionality --
# Minimum number of PCA components to explain 95% of variance.
# Indicates how many independent dimensions the network effectively uses.
compute_intrinsic_dim <- function(pca_obj, threshold = 0.95) {
  which(pca_obj$cumvar_explained >= threshold)[1]
}

# -- Compute all statistics for each dataset --
stats_summary <- tibble(
  Dataset              = c("Real (augmented)", "Synthetic", "Real (reference)"),
  N                    = c(nrow(feat_aug), nrow(feat_syn), nrow(feat_real)),
  Magnitude            = c(compute_magnitude(feat_aug),
                           compute_magnitude(feat_syn),
                           compute_magnitude(feat_real)),
  Concentricity        = c(compute_concentricity(feat_aug),
                           compute_concentricity(feat_syn),
                           compute_concentricity(feat_real)),
  Uniformity           = c(compute_uniformity(feat_aug),
                           compute_uniformity(feat_syn),
                           compute_uniformity(feat_real)),
  Tolerance            = c(compute_tolerance(feat_aug,  class_aug),
                           compute_tolerance(feat_syn,  class_syn),
                           compute_tolerance(feat_real, class_real)),
  Intrinsic_Dim_95pct  = c(compute_intrinsic_dim(pca_aug),
                           compute_intrinsic_dim(pca_syn),
                           compute_intrinsic_dim(pca_real)),
  Sparsity             = c(sparsity_aug$overall_sparsity,
                           sparsity_syn$overall_sparsity,
                           sparsity_real$overall_sparsity)
) |>
  mutate(across(where(is.double), \(x) round(x, 4)))

# Clean printed table for console inspection
View(stats_summary)


# =============================================================================
# Save Processed Structures for downstream the rest of the scripts
# =============================================================================

# 1. Save the processed feature matrices 
saveRDS(feat_aug,  file = here("feat_aug.rds"))
saveRDS(feat_syn,  file = here("feat_syn.rds"))
saveRDS(feat_real, file = here("feat_real.rds"))

# 2. Save the engineered class factor vectors
saveRDS(class_aug_f,  file = here("class_aug_f.rds"))
saveRDS(class_syn_f,  file = here("class_syn_f.rds"))
saveRDS(class_real_f, file = here("class_real_f.rds"))

# =============================================================================
# Preprocessing for similarity measures: two-step subsampling
# =============================================================================
# Step 1. Reduce augmented matrix from N = 18,200 to N = 18,027 via simple
#         random subsampling without replacement, so both primary matrices have
#         the same number of rows. The augmented dataset is already balanced
#         across the seven classes so stratification is not required.
# Step 2. Draw a further subsample of N = 3,000 from each primary matrix.
#         This is necessary because measures such as CKA and RSA involve N x N
#         Gram matrices which are infeasible at N = 18,027 on standard hardware.
#         The same subsample is used for all six similarity measures so that
#         results are directly comparable (Kornblith et al., 2019).

set.seed(42)
N_SYN     <- nrow(feat_syn)   # 18,027 — target size for augmented
N_SUBSAMPLE <- 3000L

# -- Step 1: reduce augmented to N_SYN --
idx_aug_full <- sample(nrow(feat_aug), N_SYN)
feat_aug_full  <- feat_aug[idx_aug_full, ]
class_aug_full <- class_aug_f[idx_aug_full]

# -- Step 2: draw N = 3,000 subsample from each --
idx_aug_3k <- sample(N_SYN, N_SUBSAMPLE)
idx_syn_3k <- sample(N_SYN, N_SUBSAMPLE)

feat_aug_3k  <- feat_aug_full[idx_aug_3k, ]
feat_syn_3k  <- feat_syn[idx_syn_3k, ]
class_aug_3k <- class_aug_full[idx_aug_3k]
class_syn_3k <- class_syn_f[idx_syn_3k]

# Confirm dimensions
stopifnot(nrow(feat_aug_3k) == N_SUBSAMPLE,
          nrow(feat_syn_3k) == N_SUBSAMPLE,
          ncol(feat_aug_3k) == N_FEATURES,
          ncol(feat_syn_3k) == N_FEATURES)


print(table(class_aug_3k))
print(table(class_syn_3k))

# 3. Save the 3,000-row subsamples used by all similarity measure scripts
saveRDS(feat_aug_3k,  file = here("feat_aug_3k.rds"))
saveRDS(feat_syn_3k,  file = here("feat_syn_3k.rds"))
saveRDS(class_aug_3k, file = here("class_aug_3k.rds"))
saveRDS(class_syn_3k, file = here("class_syn_3k.rds"))

# -----------------------------------------------------------------------------
# Step 3. Stratified N = 3,000 subsample of the real reference dataset
# -----------------------------------------------------------------------------
# The real reference dataset is heavily imbalanced across classes (see Data
# Description chapter), unlike the augmented and synthetic datasets which are
# already approximately balanced. A simple random subsample would inherit
# this imbalance and under-represent rare classes (e.g. SVTAC). A stratified
# subsample is drawn instead, allocating observations to each class
# proportionally to its share of the real dataset, with a minimum floor per
# class so that even the rarest class retains adequate representation. This
# subsample is required for comparison axes 1 and 2 (real vs. augmented,
# overall and per class) in the four-axis comparison framework.

stratified_subsample_idx <- function(class_vector, total_n,
                                     min_per_class = 30L) {
  classes      <- sort(unique(class_vector))
  class_counts <- table(class_vector)[as.character(classes)]
  
  # Proportional allocation with a minimum floor per class
  proportional <- total_n * class_counts / sum(class_counts)
  n_per_class  <- pmax(floor(proportional), min_per_class)
  n_per_class  <- pmin(n_per_class, class_counts)
  
  # Largest-remainder method: distribute any shortfall/excess from rounding
  # so the total exactly equals total_n, rather than approximately.
  shortfall <- total_n - sum(n_per_class)
  if (shortfall > 0) {
    # Classes with room to grow, ranked by largest fractional remainder
    remainder   <- proportional - floor(proportional)
    can_grow    <- which(n_per_class < class_counts)
    grow_order  <- can_grow[order(remainder[can_grow], decreasing = TRUE)]
    for (i in seq_len(min(shortfall, length(grow_order)))) {
      cl <- grow_order[i]
      n_per_class[cl] <- n_per_class[cl] + 1
    }
  } else if (shortfall < 0) {
    # Should not normally occur given floor() above, but guard anyway
    can_shrink  <- which(n_per_class > min_per_class)
    shrink_order <- can_shrink[order(n_per_class[can_shrink], decreasing = TRUE)]
    for (i in seq_len(min(-shortfall, length(shrink_order)))) {
      cl <- shrink_order[i]
      n_per_class[cl] <- n_per_class[cl] - 1
    }
  }
  
  stopifnot(sum(n_per_class) == total_n)
  
  unlist(lapply(classes, function(cl) {
    rows <- which(class_vector == cl)
    sample(rows, n_per_class[as.character(cl)])
  }))
}

idx_real_3k <- stratified_subsample_idx(class_real, N_SUBSAMPLE)

feat_real_3k  <- feat_real[idx_real_3k, ]
class_real_3k <- class_real_f[idx_real_3k]

stopifnot(ncol(feat_real_3k) == N_FEATURES)

print(table(class_real_3k))

saveRDS(feat_real_3k,  file = here("feat_real_3k.rds"))
saveRDS(class_real_3k, file = here("class_real_3k.rds"))