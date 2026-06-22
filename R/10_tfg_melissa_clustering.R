# =============================================================================
# 10_tfg_melissa_clustering.R  --  k-Medoids Clustering and NMI
# TFG - Grau en Estadística, UB/UPC
# Author: Melissa Vargas Jerez
# Description: Applies k-medoids (PAM) clustering to the full pooled dataset
#              (real + augmented + synthetic, 9,000 x 256) and evaluates the
#              resulting clusters using normalised mutual information (NMI)
#              against arrhythmia class labels and dataset origin labels.
#              Profiling tables cross-tabulate cluster membership against both
#              dimensions to characterise detected differences.
#              Requires feat_real_3k.rds, feat_aug_3k.rds, feat_syn_3k.rds,
#              class_real_3k.rds, class_aug_3k.rds, class_syn_3k.rds
#              produced by 01_tfg_melissa_ear.R
# =============================================================================

library(here)
library(tibble)
library(dplyr)
library(cluster)    # for pam()
library(aricode)    # for NMI()

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
# 2. Build pooled dataset
#    All three N = 3,000 matrices are combined into one 9,000 x 256 matrix.
#    Two label vectors are created alongside it:
#      - class_all:  arrhythmia class (SBRAD, SR, AFIB, ...)
#      - origin_all: dataset of origin (Real, Augmented, Synthetic)
# -----------------------------------------------------------------------------

feat_pooled <- rbind(feat_real, feat_aug, feat_syn)

class_all <- factor(c(as.character(class_real),
                      as.character(class_aug),
                      as.character(class_syn)),
                    levels = levels(class_real))

origin_all <- factor(
  c(rep("Real",      nrow(feat_real)),
    rep("Augmented", nrow(feat_aug)),
    rep("Synthetic",  nrow(feat_syn))),
  levels = c("Real", "Augmented", "Synthetic")
)

stopifnot(nrow(feat_pooled) == length(class_all),
          nrow(feat_pooled) == length(origin_all))

# -----------------------------------------------------------------------------
# 3. k-Medoids clustering (PAM) with k = 7
#    k = 7 matches the number of arrhythmia classes, so that the degree of
#    class recovery can be directly assessed. Euclidean distance is used
#    since the 256-dimensional representation space has a shared coordinate
#    system (all three datasets passed through the same frozen network).
#    PAM is chosen over k-means because it represents each cluster by an
#    actual observation (the medoid), making results directly interpretable
#    as representative patients.
# -----------------------------------------------------------------------------

set.seed(42)
pam_fit <- pam(feat_pooled, k = 7, metric = "euclidean", stand = FALSE)

cluster_labels <- factor(paste0("C", pam_fit$clustering))

# Total within-cluster dissimilarity (objective minimised by PAM)
total_dissimilarity <- tibble(
  Metric = "Total within-cluster dissimilarity (PAM objective)",
  Value  = round(pam_fit$objective["swap"], 4)
)

print(total_dissimilarity)

# -----------------------------------------------------------------------------
# 4. Normalised Mutual Information
#    NMI is computed twice:
#      (a) between cluster assignment and arrhythmia class label
#      (b) between cluster assignment and dataset origin label
#    Values close to 1 indicate strong agreement; values close to 0 indicate
#    independence between the two partitions.
# -----------------------------------------------------------------------------

nmi_class  <- NMI(cluster_labels, class_all)
nmi_origin <- NMI(cluster_labels, origin_all)

nmi_summary <- tibble(
  Comparison                         = c("Clusters vs. Arrhythmia class",
                                         "Clusters vs. Dataset origin"),
  NMI                                = round(c(nmi_class, nmi_origin), 4)
)

print(nmi_summary)

# -----------------------------------------------------------------------------
# 5. Profiling: cross-tabulations for interpretation
# -----------------------------------------------------------------------------

# -- 5.1 Clusters x Arrhythmia class --
# Shows whether each cluster is dominated by one class or mixed.
profile_class <- as.data.frame.matrix(
  table(Cluster = cluster_labels, Class = class_all)
)
print(profile_class)

# Row-proportions: what fraction of each cluster belongs to each class?
profile_class_prop <- round(
  prop.table(table(Cluster = cluster_labels, Class = class_all), margin = 1) * 100,
  1
)
print(profile_class_prop)

# -- 5.2 Clusters x Dataset origin --
# Shows whether each cluster is dominated by one dataset or mixed.
profile_origin <- as.data.frame.matrix(
  table(Cluster = cluster_labels, Origin = origin_all)
)
print(profile_origin)

# Row-proportions: what fraction of each cluster comes from each dataset?
profile_origin_prop <- round(
  prop.table(table(Cluster = cluster_labels, Origin = origin_all), margin = 1) * 100,
  1
)
print(profile_origin_prop)

# -- 5.3 Combined profiling tibble (useful for LaTeX table export) --
profile_combined <- tibble(
  Cluster   = rownames(profile_class),
  SBRAD     = profile_class$SBRAD,
  SR        = profile_class$SR,
  AFIB      = profile_class$AFIB,
  STACH     = profile_class$STACH,
  AFLT      = profile_class$AFLT,
  SARRH     = profile_class$SARRH,
  SVTAC     = profile_class$SVTAC,
  Real      = profile_origin$Real,
  Augmented = profile_origin$Augmented,
  Synthetic = profile_origin$Synthetic,
  Total     = rowSums(profile_class)
)

print(profile_combined)

# -----------------------------------------------------------------------------
# 6. Optimal k investigation (silhouette method)
#    Computes average silhouette width for k = 2 to 10 to check whether
#    k = 7 is supported by the data or whether a different k is preferred.
#    This is reported as a supplementary diagnostic, not the primary result.
# -----------------------------------------------------------------------------

sil_widths <- sapply(2:10, function(k) {
  fit <- pam(feat_pooled, k = k, metric = "euclidean", stand = FALSE)
  fit$silinfo$avg.width
})

sil_summary <- tibble(
  k              = 2:10,
  Avg_Silhouette = round(sil_widths, 4)
)

print(sil_summary)

# -----------------------------------------------------------------------------
# 7. Save results
# -----------------------------------------------------------------------------

clustering_results <- list(
  pam_fit              = pam_fit,
  cluster_labels       = cluster_labels,
  nmi_summary          = nmi_summary,
  nmi_class            = nmi_class,
  nmi_origin           = nmi_origin,
  profile_class        = profile_class,
  profile_class_prop   = profile_class_prop,
  profile_origin       = profile_origin,
  profile_origin_prop  = profile_origin_prop,
  profile_combined     = profile_combined,
  sil_summary          = sil_summary,
  total_dissimilarity  = total_dissimilarity,
  k                    = 7L,
  N_pooled             = nrow(feat_pooled),
  D                    = D
)

saveRDS(clustering_results, file = here("clustering_results.rds"))
message("Results saved to clustering_results.rds")