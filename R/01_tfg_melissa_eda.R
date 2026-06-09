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

# Classes colour palette
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
# 1. Load data
# -----------------------------------------------------------------------------

# Link to download .RData object is in the README section of this repository
# Assuming RData object and script match current work directory

load(here("data", "data_loaded.RData"))

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
# 4. Class distribution plots
# -----------------------------------------------------------------------------

build_class_counts <- function(class_factor, dataset_name) {
  data.frame(
    class      = levels(class_factor),
    count      = as.integer(table(class_factor)),
    proportion = as.numeric(prop.table(table(class_factor))),
    dataset    = dataset_name
  )
}

counts_aug  <- build_class_counts(class_aug_f,  "Real (augmented)")
counts_syn  <- build_class_counts(class_syn_f,  "Synthetic")
counts_real <- build_class_counts(class_real_f, "Real (reference)")

counts_all <- bind_rows(counts_aug, counts_syn, counts_real) |>
  mutate(dataset = factor(dataset,
                          levels = c("Real (augmented)",
                                     "Synthetic",
                                     "Real (reference)")))

plot_class_dist <- ggplot(counts_all,
                          aes(x = class, y = count, fill = dataset)) +
  geom_col(position = "dodge", width = 0.7, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = DATASET_COLOURS) +
  scale_y_continuous(labels = comma_format()) +
  labs(
    title    = "Class distribution across representation matrices",
    subtitle = "Each bar represents the number of patients per arrhythmia class",
    x        = "Arrhythmia class",
    y        = "Number of patients",
    fill     = "Dataset"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position    = "bottom",
    plot.title         = element_text(face = "bold")
  )

ggsave(here("figures", "class_distribution.pdf"),
       plot_class_dist, width = 10, height = 5)

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

# Histogram: distribution of per-neuron sparsity across 256 neurons
sparsity_neuron_df <- data.frame(
  sparsity = c(sparsity_aug$per_neuron,
               sparsity_syn$per_neuron,
               sparsity_real$per_neuron),
  dataset  = rep(c("Real (augmented)", "Synthetic", "Real (reference)"),
                 each = N_FEATURES)
) |>
  mutate(dataset = factor(dataset,
                          levels = c("Real (augmented)",
                                     "Synthetic",
                                     "Real (reference)")))

plot_sparsity_hist <- ggplot(sparsity_neuron_df,
                             aes(x = sparsity, fill = dataset)) +
  geom_histogram(bins = 40, colour = "white", linewidth = 0.2,
                 position = "identity", alpha = 0.75) +
  scale_fill_manual(values = DATASET_COLOURS) +
  scale_x_continuous(labels = percent_format()) +
  facet_wrap(~dataset, ncol = 1) +
  labs(
    title    = "Distribution of per-neuron sparsity across 256 features",
    subtitle = "Proportion of patients for which a given neuron is inactive (value = 0)",
    x        = "Sparsity rate",
    y        = "Number of neurons"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        plot.title      = element_text(face = "bold"))

ggsave(here("figures", "sparsity_histogram.pdf"),
       plot_sparsity_hist, width = 8, height = 8)

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

# Report variance explained by first few components
cat("\n=== PCA: cumulative variance explained ===\n")
for (k in c(2, 5, 10, 20, 50)) {
  cat(sprintf("  First %2d PCs — Aug: %.1f%%  Syn: %.1f%%  Real: %.1f%%\n",
              k,
              pca_aug$cumvar_explained[k]  * 100,
              pca_syn$cumvar_explained[k]  * 100,
              pca_real$cumvar_explained[k] * 100))
}

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

# 2D scatter plots of PC1 vs PC2, coloured by arrhythmia class
# Plot on a random subsample for visual clarity (max 3000 points per dataset)
set.seed(42)
MAX_PLOT_POINTS <- 3000L

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
scatter_real <- build_pca_scatter_df(pca_real, class_real_f, "Real (reference)")

plot_pca_scatter <- function(scatter_df, pca_obj, dataset_name) {
  pct1 <- round(pca_obj$var_explained[1] * 100, 1)
  pct2 <- round(pca_obj$var_explained[2] * 100, 1)
  ggplot(scatter_df, aes(x = pc1, y = pc2, colour = class)) +
    geom_point(size = 0.6, alpha = 0.5) +
    scale_colour_manual(values = CLASS_COLOURS,
                        name   = "Arrhythmia class") +
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


View(stats_summary)


# =============================================================================
# Save Processed Structures for downstream script 02_*
# =============================================================================

# =============================================================================
# Save Processed Structures for downstream scripts (02_*, 03_*)
# =============================================================================

# 1. Save the processed feature matrices 
saveRDS(feat_aug,  file = here("feat_aug.rds"))
saveRDS(feat_syn,  file = here("feat_syn.rds"))
saveRDS(feat_real, file = here("feat_real.rds"))

# 2. Save the engineered class factor vectors
saveRDS(class_aug_f,  file = here("class_aug_f.rds"))
saveRDS(class_syn_f,  file = here("class_syn_f.rds"))
saveRDS(class_real_f, file = here("class_real_f.rds"))
