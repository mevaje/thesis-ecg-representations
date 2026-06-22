# =============================================================================
# 08_tfg_melissa_results_tables.R
# TFG - Grau en Estadística, UB/UPC
# Author: Melissa Vargas Jerez
# Description: Reads all six similarity measure results and generates
#              ready-to-paste LaTeX table code for the Results chapter.
#              Run this after all six measure scripts (02-07) have completed.
# =============================================================================

library(here)
library(dplyr)
library(knitr)

# -----------------------------------------------------------------------------
# 1. Load all results
# -----------------------------------------------------------------------------

cca        <- readRDS(here("cca_results.rds"))
rv         <- readRDS(here("rv_results.rds"))
cka        <- readRDS(here("cka_results.rds"))
rsa        <- readRDS(here("rsa_results.rds"))
procrustes <- readRDS(here("procrustes_results.rds"))
jaccard    <- readRDS(here("jaccard_results.rds"))

# -----------------------------------------------------------------------------
# 2. Table 1: Overall comparison across both axes (Axes 1 and 3)
#    One row per axis, one column per measure. This is the main results table.
# -----------------------------------------------------------------------------

overall_summary <- data.frame(
  Axis = c("Real vs.\\ Augmented", "Augmented vs.\\ Synthetic"),
  CCA  = c(round(cca$cca_overall_summary$m_CCA[1], 3),
           round(cca$cca_overall_summary$m_CCA[2], 3)),
  RV   = c(round(rv$rv_overall_summary$RV[1], 4),
           round(rv$rv_overall_summary$RV[2], 4)),
  CKA  = c(round(cka$cka_overall_summary$CKA[1], 4),
           round(cka$cka_overall_summary$CKA[2], 4)),
  RSA  = c(round(rsa$rsa_overall_summary$m_RSA[1], 4),
           round(rsa$rsa_overall_summary$m_RSA[2], 4)),
  Procrustes = c(round(procrustes$proc_overall_summary$sim_Proc[1], 3),
                 round(procrustes$proc_overall_summary$sim_Proc[2], 3)),
  Jaccard_k10 = c(round(jaccard$jaccard_k10_summary$m_Jac[1], 4),
                  round(jaccard$jaccard_k10_summary$m_Jac[2], 4))
)

cat("\n% ============================================================\n")
cat("% TABLE 1: Overall results across both axes\n")
cat("% Paste into 11_results_discussion.tex\n")
cat("% ============================================================\n\n")

cat(knitr::kable(overall_summary,
                 format  = "latex",
                 booktabs = TRUE,
                 escape  = FALSE,
                 col.names = c("Comparison",
                               "$m_{\\text{CCA}}$",
                               "RV",
                               "CKA",
                               "$m_{\\text{RSA}}$",
                               "Procrustes ($\\text{sim}$)",
                               "Jaccard ($k=10$)"),
                 caption = "Summary of similarity scores across the six measures
                             for the two overall comparison axes. All measures
                             are bounded in $[0,1]$ except RSA which is bounded
                             in $[-1,1]$. Higher values indicate greater
                             similarity in all cases except Procrustes distance
                             $m_{\\text{Proc}}$, which is reported as similarity
                             $1 - m_{\\text{Proc}}/\\sqrt{2}$.",
                 label   = "tab:overall_results"))

# -----------------------------------------------------------------------------
# 3. Table 2: Per-class results for Axis 2 (Real vs. Augmented)
# -----------------------------------------------------------------------------

axis2 <- data.frame(
  Class = cca$cca_per_class_summary$Class[1:7],
  N     = cca$cca_per_class_summary$N[1:7],
  CCA   = round(cca$cca_per_class_summary$m_CCA[1:7], 3),
  RV    = round(rv$rv_per_class_summary$RV[1:7], 4),
  CKA   = round(cka$cka_per_class_summary$CKA[1:7], 4),
  RSA   = round(rsa$rsa_per_class_summary$RSA[1:7], 4),
  Proc  = round(procrustes$proc_per_class_summary$sim_Proc[1:7], 3),
  Jac   = round(jaccard$jaccard_per_class_summary$m_Jac[1:7], 4)
)

cat("\n\n% ============================================================\n")
cat("% TABLE 2: Per-class results, Axis 2 (Real vs. Augmented)\n")
cat("% ============================================================\n\n")

cat(knitr::kable(axis2,
                 format   = "latex",
                 booktabs = TRUE,
                 escape   = FALSE,
                 col.names = c("Class", "$N$",
                               "$m_{\\text{CCA}}$", "RV", "CKA",
                               "$m_{\\text{RSA}}$",
                               "Procrustes", "Jaccard"),
                 caption = "Per-class similarity scores for Axis 2 (real
                             reference versus augmented data). $N$ denotes the
                             number of observations used in the comparison,
                             matched to the smaller of the two class sizes.",
                 label   = "tab:axis2_per_class"))

# -----------------------------------------------------------------------------
# 4. Table 3: Per-class results for Axis 4 (Augmented vs. Synthetic)
# -----------------------------------------------------------------------------

axis4 <- data.frame(
  Class = cca$cca_per_class_summary$Class[8:14],
  N     = cca$cca_per_class_summary$N[8:14],
  CCA   = round(cca$cca_per_class_summary$m_CCA[8:14], 3),
  RV    = round(rv$rv_per_class_summary$RV[8:14], 4),
  CKA   = round(cka$cka_per_class_summary$CKA[8:14], 4),
  RSA   = round(rsa$rsa_per_class_summary$RSA[8:14], 4),
  Proc  = round(procrustes$proc_per_class_summary$sim_Proc[8:14], 3),
  Jac   = round(jaccard$jaccard_per_class_summary$m_Jac[8:14], 4)
)

cat("\n\n% ============================================================\n")
cat("% TABLE 3: Per-class results, Axis 4 (Augmented vs. Synthetic)\n")
cat("% ============================================================\n\n")

cat(knitr::kable(axis4,
                 format   = "latex",
                 booktabs = TRUE,
                 escape   = FALSE,
                 col.names = c("Class", "$N$",
                               "$m_{\\text{CCA}}$", "RV", "CKA",
                               "$m_{\\text{RSA}}$",
                               "Procrustes", "Jaccard"),
                 caption = "Per-class similarity scores for Axis 4 (augmented
                             versus synthetic data). $N$ denotes the number of
                             observations used in the comparison, matched to
                             the smaller of the two class sizes.",
                 label   = "tab:axis4_per_class"))

# -----------------------------------------------------------------------------
# 5. Table 4: k-NN Jaccard sensitivity to k (Axes 1 and 3)
# -----------------------------------------------------------------------------

jaccard_k <- jaccard$jaccard_overall_summary |>
  select(Axis, k, m_Jac) |>
  tidyr::pivot_wider(names_from = k, values_from = m_Jac,
                     names_prefix = "$k=") |>
  rename(Comparison = Axis) |>
  mutate(Comparison = c("Real vs.\\ Augmented",
                        "Augmented vs.\\ Synthetic"))

cat("\n\n% ============================================================\n")
cat("% TABLE 4: k-NN Jaccard sensitivity to k\n")
cat("% ============================================================\n\n")

cat(knitr::kable(jaccard_k,
                 format   = "latex",
                 booktabs = TRUE,
                 escape   = FALSE,
                 col.names = c("Comparison",
                               "$k=5$", "$k=10$", "$k=20$", "$k=50$"),
                 caption = "Global $k$-NN Jaccard similarity for both overall
                             comparison axes across four values of $k$. The
                             primary reported value is $k=10$.",
                 label   = "tab:jaccard_k"))

cat("\n\n% ============================================================\n")
cat("% All tables generated. Copy each block above into Overleaf.\n")
cat("% Add \\usepackage{booktabs} to main.tex if not already present.\n")
cat("% ============================================================\n")