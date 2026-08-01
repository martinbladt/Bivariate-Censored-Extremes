# Loss--ALAE data application
#
# Data: lossalaefull from the CASdatasets package (1,500 observations).
# The loss margin is right-censored; ALAE is fully observed.
# Required packages: CASdatasets and mhazard.
# Run with:
#   Rscript loss_alae.R
#
# Seven diagnostic and estimation figures are written to figures/.

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1L])))
  } else {
    normalizePath(getwd())
  }
}

script_dir <- get_script_dir()
source(file.path(script_dir, "functions.R"))

if (!requireNamespace("CASdatasets", quietly = TRUE)) {
  stop("Package 'CASdatasets' is required.", call. = FALSE)
}

# Form the observed vector (Z1, Z2) and censoring indicators (d1, d2).
data("lossalaefull", package = "CASdatasets")
lossalae_data <- data.frame(
  Z1 = lossalaefull$Loss,
  d1 = as.integer(lossalaefull$Censored == 0),
  Z2 = lossalaefull$ALAE,
  d2 = 1L
)

# Threshold ranges and evaluation grids used in the application.
lossalae_summary <- real_data_functions$generate_real_data_figures(
  dataframe = lossalae_data,
  prefix = "lossalae",
  k_values = seq(50L, 500L, by = 10L),
  curve_k = 450L,
  curve_t_values = 3,
  s_grid = seq(1.1, 5, by = 0.01),
  surface_s_grid = seq(1.05, 5, length.out = 60L),
  surface_t_grid = seq(1.05, 5, length.out = 60L),
  n_B_limits = c(1L, 40L),
  n_A_limits = c(1L, 40L),
  adaptive_by = 1L,
  region_sizes = c(50L, 20L, 10L),
  show_c_curve = FALSE
)

print(lossalae_summary)
