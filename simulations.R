# Simulation studies for the bivariate censored-extremes paper
#
# Statistical software: R with the copula and mhazard packages.
# Run with:
#   Rscript simulations.R
#
# The script generates the four tail-set illustrations, the Monte Carlo
# coverage figures, and four LaTeX summary tables. Default settings use 1,000
# replications in each of 24 scenarios and may require substantial run time.

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1L])))
  } else {
    normalizePath(getwd())
  }
}

script_dir <- get_script_dir()
figure_dir <- file.path(script_dir, "figures")
table_dir <- file.path(script_dir, "tables")

source(file.path(script_dir, "functions.R"))

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)


# Illustrate the tail sets A(1,1), A(s,t), B(s,t), C(s,t), and D(s,t).
generate_set_figures <- function(figure_dir) {
  set.seed(1234)

  n <- 200L
  k <- 130L
  theta <- 1.5
  s <- 1.5
  t <- 2.0
  gamma_X <- c(1, 1)
  censor_props <- c(0.10, 0.10)
  gamma_C <- (1 - censor_props) / censor_props * gamma_X

  df <- simulation_functions$sim_data_pareto(n, theta, gamma_X, gamma_C)
  U1 <- simulation_functions$tail_margin_fit_hill(df$Z1, df$d1, k)$U
  U2 <- simulation_functions$tail_margin_fit_hill(df$Z2, df$d2, k)$U

  # Membership indicators for the observable sets in the estimator.
  A11 <- df$Z1 >= U1 & df$Z2 >= U2
  Ast <- df$Z1 >= s * U1 & df$Z2 >= t * U2
  Bst <- df$Z1 >= s * U1 & df$Z2 >= t * U2 & df$d1 == 1L & df$d2 == 1L
  Cst <- df$Z1 >= s * U1 & df$Z2 >= t * U2 & df$d1 == 1L
  Dst <- df$Z1 >= s * U1 & df$Z2 >= t * U2 & df$d2 == 1L

  panel_defs <- list(
    list(target = Ast, file = "sets_AA.pdf"),
    list(target = Bst, file = "sets_AB.pdf"),
    list(target = Cst, file = "sets_AC.pdf"),
    list(target = Dst, file = "sets_AD.pdf")
  )

  if (interactive()) {
    graphics::par(mar = c(2.2, 2.2, 0.8, 0.8))
    for (panel in panel_defs) {
      simulation_functions$plot_set_panel(
        df, U1, U2, s, t,
        target = panel$target,
        target_col = simulation_functions$tc
      )
    }
  }
  for (panel in panel_defs) {
    grDevices::pdf(file.path(figure_dir, panel$file), width = 5, height = 5)
    graphics::par(mar = c(2.2, 2.2, 0.8, 0.8))
    simulation_functions$plot_set_panel(
      df, U1, U2, s, t,
      target = panel$target,
      target_col = simulation_functions$tc
    )
    grDevices::dev.off()
  }
}


# Monte Carlo size and nominal confidence level.
n_rep <- 1000L
alpha <- 0.05
n_cores <- 10L

# Parameters of the adaptive threshold rule.
region_sizes <- c(50L, 20L, 10L)
n_B_limits <- c(1L, 40L)
n_A_limits <- c(1L, 40L)
adaptive_by <- 2L

# Dependence parameters, sample sizes, evaluation points, and censoring rates.
theta_values <- c(1.5, 3)
n_values <- c(300L, 1000L)
points_list <- list(
  P1 = c(1.5, 1.5),
  P2 = c(1.25, 2.5)
)
censor_list <- list(
  C1 = c(0.05, 0.05),
  C2 = c(0.05, 0.10),
  C3 = c(0.20, 0.20)
)

k_grid_for_n <- function(n) {
  if (n == 300L) {
    seq(40L, 150L, by = 10L)
  } else {
    seq(80L, 300L, by = 20L)
  }
}

theta_tag <- function(theta) {
  if (abs(theta - 1.5) < 1e-12) "1p5" else "3"
}

point_label <- function(point_id) {
  if (point_id == "P1") "$(1.5,1.5)$" else "$(1.25,2.5)$"
}

censor_label <- function(censor_id) {
  switch(
    censor_id,
    C1 = "$(5\\%,5\\%)$",
    C2 = "$(5\\%,10\\%)$",
    C3 = "$(20\\%,20\\%)$"
  )
}

fmt_num <- function(x, digits = 4L) {
  if (!is.finite(x)) {
    "--"
  } else {
    formatC(x, format = "f", digits = digits)
  }
}

write_block_table <- function(block, file) {
  block <- block[order(block$point_id, block$censor_id), ]

  lines <- c(
    "\\begin{tabular}{llrrrrr}",
    "\\hline",
    "$(s,t)$ & censoring & true value & mean estimate & mean $k^\\ast$ & mean bias & mean squared error\\\\",
    "\\hline"
  )

  for (i in seq_len(nrow(block))) {
    lines <- c(
      lines,
      sprintf(
        "%s & %s & %s & %s & %s & %s & %s\\\\",
        point_label(block$point_id[i]),
        censor_label(block$censor_id[i]),
        fmt_num(block$truth[i], 4L),
        fmt_num(block$mean_Fhat[i], 4L),
        fmt_num(block$mean_k_center[i], 1L),
        fmt_num(block$bias[i], 4L),
        fmt_num(block$mse[i], 5L)
      )
    )
  }

  lines <- c(lines, "\\hline", "\\end{tabular}")
  writeLines(lines, con = file)
}

run_one_scenario <- function(theta, n, point_id, censor_id) {
  point <- points_list[[point_id]]
  censor_props <- censor_list[[censor_id]]
  k_values <- k_grid_for_n(n)

  study <- simulation_functions$sensitivity_study_k(
    n_rep = n_rep,
    n = n,
    k_values = k_values,
    s = point[1L],
    t = point[2L],
    theta = theta,
    censor_props = censor_props,
    alpha = alpha,
    n_cores = n_cores,
    progress = TRUE,
    region_sizes = region_sizes,
    n_B_limits = n_B_limits,
    n_A_limits = n_A_limits,
    adaptive_by = adaptive_by
  )

  scenario_name <- sprintf(
    "theta%s_n%d_%s_%s",
    theta_tag(theta), n, tolower(censor_id), tolower(point_id)
  )
  pdf(file.path(figure_dir, paste0(scenario_name, ".pdf")), width = 5, height = 5)
  simulation_functions$plot_sensitivity_k(study, main = "")
  dev.off()

  adaptive <- study$results_table[1L, c(
    "mean_k_center", "usable_reps", "coverage", "mean_m",
    "mean_n_B", "mean_Fhat", "bias", "mse", "truth"
  )]

  data.frame(
    theta = theta,
    n = n,
    point_id = point_id,
    censor_id = censor_id,
    s = point[1L],
    t = point[2L],
    adaptive_plot = file.path("figures", paste0(scenario_name, ".pdf")),
    mean_k_center = adaptive$mean_k_center,
    usable_reps = adaptive$usable_reps,
    coverage = adaptive$coverage,
    mean_m = adaptive$mean_m,
    mean_n_B = adaptive$mean_n_B,
    mean_Fhat = adaptive$mean_Fhat,
    bias = adaptive$bias,
    mse = adaptive$mse,
    truth = adaptive$truth
  )
}

# Generate the set illustrations before running the Monte Carlo study.
generate_set_figures(figure_dir)


# Construct and run the 24 simulation scenarios.
scenario_grid <- expand.grid(
  theta = theta_values,
  n = n_values,
  censor_id = names(censor_list),
  point_id = names(points_list),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
scenario_grid <- scenario_grid[
  order(
    scenario_grid$theta,
    scenario_grid$n,
    scenario_grid$censor_id,
    scenario_grid$point_id
  ),
]

summary_rows <- vector("list", nrow(scenario_grid))

for (i in seq_len(nrow(scenario_grid))) {
  row <- scenario_grid[i, ]
  cat(sprintf(
    "\nScenario %02d/%02d: theta = %s, n = %d, %s, %s\n",
    i, nrow(scenario_grid), row$theta, row$n, row$point_id, row$censor_id
  ))
  summary_rows[[i]] <- run_one_scenario(
    theta = row$theta,
    n = row$n,
    point_id = row$point_id,
    censor_id = row$censor_id
  )
}

summary_table <- do.call(rbind, summary_rows)

block_keys <- unique(summary_table[c("theta", "n")])
for (i in seq_len(nrow(block_keys))) {
  theta <- block_keys$theta[i]
  n <- block_keys$n[i]
  block <- summary_table[summary_table$theta == theta & summary_table$n == n, ]
  write_block_table(
    block,
    file.path(table_dir, sprintf("table_theta%s_n%d.tex", theta_tag(theta), n))
  )
}

print(summary_table, row.names = FALSE, digits = 4)
