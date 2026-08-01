# Estimation and simulation functions for bivariate censored extremes
#
# Statistical convention:
#   Z_j = min(X_j, C_j) is the observed value on margin j;
#   d_j = 1 denotes an uncensored observation and d_j = 0 a censored one;
#   k is the number of upper-order observations used for tail estimation;
#   (s, t) are standardized excess levels above the marginal thresholds.
#
# The Monte Carlo study and the data applications use related calculations
# with slightly different arguments. They are grouped separately to keep the
# notation unambiguous:
#   simulation_functions  calculations used in simulations.R
#   real_data_functions   calculations used in loss_alae.R and uk_data.R


# =============================================================================
# Simulation functions
# =============================================================================

simulation_functions <- local({

require_package <- function(pkg) {
  if (!suppressWarnings(requireNamespace(pkg, quietly = TRUE))) {
    stop(sprintf("Package '%s' is required.", pkg), call. = FALSE)
  }
}


# -------------------------------------------------------------------------
# Data generation and progress display
# -------------------------------------------------------------------------

# Burr model, tail dependence function, and Monte Carlo data generation.
q_burr <- function(u, tau, lambda) {
  ((1 - u)^(-1 / lambda) - 1)^(1 / tau)
}

R_fun <- function(x, y, theta) {
  x + y - (x^theta + y^theta)^(1 / theta)
}

F_true <- function(s, t, gamma1, gamma2, theta) {
  R_fun(s^(-1 / gamma1), t^(-1 / gamma2), theta) / R_fun(1, 1, theta)
}

simulation_parameters <- function(censor_props = c(0.10, 0.10)) {
  gamma_X1 <- 0.2
  gamma_X2 <- 0.4
  tau_X1 <- 10
  tau_X2 <- 5
  lambda_X1 <- 1 / (gamma_X1 * tau_X1)
  lambda_X2 <- 1 / (gamma_X2 * tau_X2)
  lambda_C <- c(0.5, 0.5)
  gamma_C1 <- (1 - censor_props[1L]) / censor_props[1L] * gamma_X1
  gamma_C2 <- (1 - censor_props[2L]) / censor_props[2L] * gamma_X2
  tau_C <- c(1 / (gamma_C1 * lambda_C[1L]), 1 / (gamma_C2 * lambda_C[2L]))

  list(
    gamma_X1 = gamma_X1,
    gamma_X2 = gamma_X2,
    tauX = c(tau_X1, tau_X2),
    lamX = c(lambda_X1, lambda_X2),
    tauC = tau_C,
    lamC = lambda_C
  )
}

sim_data <- function(n, theta, tauX, lamX, tauC, lamC, copula_object = NULL) {
  require_package("copula")

  if (is.null(copula_object)) {
    copula_object <- copula::gumbelCopula(theta, 2)
  }

  U <- copula::rCopula(n, copula_object)
  X1 <- q_burr(U[, 1L], tauX[1L], lamX[1L])
  X2 <- q_burr(U[, 2L], tauX[2L], lamX[2L])
  C1 <- q_burr(stats::runif(n), tauC[1L], lamC[1L])
  C2 <- q_burr(stats::runif(n), tauC[2L], lamC[2L])

  list(
    Z1 = pmin(X1, C1),
    d1 = as.integer(X1 <= C1),
    Z2 = pmin(X2, C2),
    d2 = as.integer(X2 <= C2)
  )
}

format_elapsed <- function(seconds) {
  if (!is.finite(seconds) || seconds < 0) {
    return("?")
  }

  seconds <- as.integer(round(seconds))
  hh <- seconds %/% 3600L
  mm <- (seconds %% 3600L) %/% 60L
  ss <- seconds %% 60L

  if (hh > 0L) {
    sprintf("%02dh:%02dm:%02ds", hh, mm, ss)
  } else {
    sprintf("%02dm:%02ds", mm, ss)
  }
}

run_replications_with_progress <- function(rep_ids, run_one_rep, n_cores, progress, progress_label) {
  n_rep <- length(rep_ids)
  start_time <- proc.time()[3]

  if (!progress || n_rep <= 1L) {
    if (.Platform$OS.type != "windows" && n_cores > 1L) {
      return(parallel::mclapply(rep_ids, run_one_rep, mc.cores = min(n_cores, n_rep)))
    }
    return(lapply(rep_ids, run_one_rep))
  }

  batch_size <- max(min(n_cores, n_rep), ceiling(n_rep / 20L))
  batches <- split(rep_ids, ceiling(seq_along(rep_ids) / batch_size))
  rep_out <- vector("list", n_rep)
  next_pos <- 1L

  cat(sprintf("%s: starting %d replications on %d core(s)\n", progress_label, n_rep, max(1L, min(n_cores, n_rep))))

  for (ids in batches) {
    batch_out <- if (.Platform$OS.type != "windows" && n_cores > 1L) {
      parallel::mclapply(ids, run_one_rep, mc.cores = min(n_cores, length(ids)))
    } else {
      lapply(ids, run_one_rep)
    }

    last_pos <- next_pos + length(ids) - 1L
    rep_out[next_pos:last_pos] <- batch_out
    next_pos <- last_pos + 1L

    done <- last_pos
    elapsed <- proc.time()[3] - start_time
    eta <- elapsed * (n_rep - done) / done

    cat(sprintf(
      "%s: %d/%d (%.1f%%), elapsed %s, eta %s\n",
      progress_label, done, n_rep, 100 * done / n_rep,
      format_elapsed(elapsed), format_elapsed(eta)
    ))
  }

  rep_out
}


# -------------------------------------------------------------------------
# Pointwise confidence interval at fixed k
# -------------------------------------------------------------------------

# Censoring-adjusted marginal tail fits and construction of the tail sample.
tail_margin_fit_c2 <- function(Z, d, k, gamma_method = c("Hill", "EKMI")) {
  gamma_method <- match.arg(gamma_method)
  n <- length(Z)
  ord <- order(Z)
  sZ <- Z[ord]
  sd <- d[ord]
  Z_nk <- sZ[n - k]

  r_all <- rle(sZ)
  ends_all <- cumsum(r_all$lengths)
  event_cum_all <- cumsum(sd)
  events_all <- event_cum_all[ends_all] - c(0, event_cum_all[head(ends_all, -1L)])
  n_risk_all <- n - c(0, cumsum(head(r_all$lengths, -1L)))
  keep_last <- findInterval(Z_nk, r_all$values, rightmost.closed = TRUE)

  S_at <- if (keep_last > 0L) {
    km_factor <- 1 - events_all[seq_len(keep_last)] / n_risk_all[seq_len(keep_last)]
    prod(km_factor[events_all[seq_len(keep_last)] > 0])
  } else {
    1
  }

  block <- (n - k):n
  tail_Z <- sZ[block]
  tail_d <- sd[block]

  if (gamma_method == "Hill") {
    hill_Z <- mean(log(tail_Z[-1L] / tail_Z[1L]))
    p_uncensored <- mean(tail_d[-1L])
    ghat <- hill_Z / p_uncensored
  } else {
    r_tail <- rle(tail_Z)
    ends_tail <- cumsum(r_tail$lengths)
    event_cum_tail <- cumsum(tail_d)
    events_tail <- event_cum_tail[ends_tail] - c(0, event_cum_tail[head(ends_tail, -1L)])
    n_risk_tail <- length(tail_Z) - c(0, cumsum(head(r_tail$lengths, -1L)))
    tail_surv <- cumprod(ifelse(events_tail > 0, 1 - events_tail / n_risk_tail, 1))
    tl <- tail(r_tail$values, k + 1L)
    sl <- tail(tail_surv, k + 1L)
    times <- tl[-1L] / tl[1L]
    wts <- diff(c(0, 1 - sl[-1L] / sl[1L]))
    ghat <- sum(log(times) * wts) / sum(wts)
  }

  list(gamma = ghat, U = Z_nk * (S_at / (k / n))^ghat)
}

make_tail_sample <- function(Z1, d1, Z2, d2, k, gamma_method = c("Hill", "EKMI")) {
  gamma_method <- match.arg(gamma_method)
  U1 <- tail_margin_fit_c2(Z1, d1, k, gamma_method = gamma_method)$U
  U2 <- tail_margin_fit_c2(Z2, d2, k, gamma_method = gamma_method)$U
  Y1_all <- Z1 / U1
  Y2_all <- Z2 / U2
  keep <- (Y1_all >= 1) & (Y2_all >= 1)

  list(
    Y1 = Y1_all[keep],
    Y2 = Y2_all[keep],
    d1 = as.integer(d1[keep]),
    d2 = as.integer(d2[keep]),
    m = sum(keep),
    r_n = sqrt(sum(keep))
  )
}

# Dabrowska point estimator and discrete operator utilities.
estimate_dabrowska_point <- function(tail_sample, s, t) {
  require_package("mhazard")

  T1 <- sort(unique(tail_sample$Y1[tail_sample$d1 == 1L & tail_sample$Y1 <= s]))
  T2 <- sort(unique(tail_sample$Y2[tail_sample$d2 == 1L & tail_sample$Y2 <= t]))
  Fhat <- mhazard:::jointHazLambda(
    tail_sample$Y1, tail_sample$Y2, T1, T2, tail_sample$d1, tail_sample$d2, 1L
  )[, , 1L]

  as.numeric(Fhat[length(T1) + 1L, length(T2) + 1L])
}

thin_nodes <- function(nodes, max_nodes) {
  if (is.null(max_nodes) || length(nodes) <= max_nodes) {
    return(nodes)
  }

  interior <- nodes[-c(1L, length(nodes))]
  keep_n <- max(0L, max_nodes - 2L)

  if (length(interior) == 0L || keep_n == 0L) {
    return(c(nodes[1L], nodes[length(nodes)]))
  }

  idx <- unique(round(seq(1, length(interior), length.out = keep_n)))
  sort(unique(c(nodes[1L], interior[idx], nodes[length(nodes)])))
}

fill_matrix <- function(i, j, weights, n_row, n_col) {
  out <- matrix(0, nrow = n_row, ncol = n_col)
  keep <- i >= 1L & i <= n_row & j >= 1L & j <= n_col

  if (!any(keep)) {
    return(out)
  }

  ii <- i[keep]
  jj <- j[keep]
  ww <- weights[keep]

  for (p in seq_along(ww)) {
    out[ii[p], jj[p]] <- out[ii[p], jj[p]] + ww[p]
  }

  out
}

fill_vector <- function(i, weights, n_cell) {
  out <- numeric(n_cell)
  keep <- i >= 1L & i <= n_cell

  if (!any(keep)) {
    return(out)
  }

  ii <- i[keep]
  ww <- weights[keep]

  for (p in seq_along(ww)) {
    out[ii[p]] <- out[ii[p]] + ww[p]
  }

  out
}

upper_right_cumsum <- function(x) {
  if (nrow(x) > 1L) {
    for (i in (nrow(x) - 1L):1L) x[i, ] <- x[i, ] + x[i + 1L, ]
  }
  if (ncol(x) > 1L) {
    for (j in (ncol(x) - 1L):1L) x[, j] <- x[, j] + x[, j + 1L]
  }
  x
}

upper_left_cumsum <- function(x) {
  if (nrow(x) > 1L) {
    for (i in 2:nrow(x)) x[i, ] <- x[i, ] + x[i - 1L, ]
  }
  if (ncol(x) > 1L) {
    for (j in 2:ncol(x)) x[, j] <- x[, j] + x[, j - 1L]
  }
  x
}

right_cumsum <- function(x) {
  if (ncol(x) > 1L) {
    for (j in (ncol(x) - 1L):1L) x[, j] <- x[, j] + x[, j + 1L]
  }
  x
}

left_cumsum <- function(x) {
  if (ncol(x) > 1L) {
    for (j in 2:ncol(x)) x[, j] <- x[, j] + x[, j - 1L]
  }
  x
}

down_cumsum <- function(x) {
  if (nrow(x) > 1L) {
    for (i in (nrow(x) - 1L):1L) x[i, ] <- x[i, ] + x[i + 1L, ]
  }
  x
}

top_cumsum <- function(x) {
  if (nrow(x) > 1L) {
    for (i in 2:nrow(x)) x[i, ] <- x[i, ] + x[i - 1L, ]
  }
  x
}

make_interval_index <- function(y, nodes) {
  idx <- findInterval(y, nodes, rightmost.closed = FALSE, all.inside = FALSE)
  idx[y == nodes[length(nodes)]] <- length(nodes) - 1L
  idx
}

make_threshold_index <- function(y, nodes) {
  pmax(1L, pmin(findInterval(y, nodes, rightmost.closed = TRUE), length(nodes)))
}

build_tail_operators <- function(tail_sample, s, t, max_u = NULL, max_v = NULL) {
  u <- thin_nodes(sort(unique(c(1, tail_sample$Y1[tail_sample$Y1 <= s], s))), max_u)
  v <- thin_nodes(sort(unique(c(1, tail_sample$Y2[tail_sample$Y2 <= t], t))), max_v)
  n_u <- length(u)
  n_v <- length(v)
  u_index <- make_threshold_index(tail_sample$Y1, u)
  v_index <- make_threshold_index(tail_sample$Y2, v)
  u_cell <- make_interval_index(tail_sample$Y1, u)
  v_cell <- make_interval_index(tail_sample$Y2, v)
  both <- as.numeric(tail_sample$d1 == 1L & tail_sample$d2 == 1L)
  unc1 <- as.numeric(tail_sample$d1 == 1L)
  unc2 <- as.numeric(tail_sample$d2 == 1L)

  H <- upper_right_cumsum(fill_matrix(u_index, v_index, rep(1, tail_sample$m), n_u, n_v)) / tail_sample$m
  K1 <- upper_right_cumsum(fill_matrix(u_index, v_index, both, n_u, n_v)) / tail_sample$m
  K2 <- upper_right_cumsum(fill_matrix(u_index, v_index, unc1, n_u, n_v)) / tail_sample$m
  K3 <- upper_right_cumsum(fill_matrix(u_index, v_index, unc2, n_u, n_v)) / tail_sample$m

  H_inner <- H[-n_u, -n_v, drop = FALSE]
  K1_mass <- fill_matrix(u_cell, v_cell, both, n_u - 1L, n_v - 1L) / tail_sample$m
  K2_u_full <- right_cumsum(fill_matrix(u_cell, v_index, unc1, n_u - 1L, n_v)) / tail_sample$m
  K3_v_full <- down_cumsum(fill_matrix(u_index, v_cell, unc2, n_u, n_v - 1L)) / tail_sample$m
  K2_u_mass <- K2_u_full[, -n_v, drop = FALSE]
  K3_v_mass <- K3_v_full[-n_u, , drop = FALSE]

  L_mass <- matrix(0, nrow = nrow(H_inner), ncol = ncol(H_inner))
  inner_ok <- H_inner > 0
  L_mass[inner_ok] <- (
    K2_u_mass[inner_ok] * K3_v_mass[inner_ok] - H_inner[inner_ok] * K1_mass[inner_ok]
  ) / (H_inner[inner_ok]^2)

  list(
    m = tail_sample$m,
    r_n = tail_sample$r_n,
    u = u,
    v = v,
    u_index = u_index,
    v_index = v_index,
    u_cell = u_cell,
    v_cell = v_cell,
    both = both,
    unc1 = unc1,
    unc2 = unc2,
    H_inner = H_inner,
    H_u_edge = H[-n_u, 1L],
    H_v_edge = H[1L, -n_v],
    K1_mass = K1_mass,
    K2_u_mass = K2_u_mass,
    K3_v_mass = K3_v_mass,
    K2_u_edge = fill_vector(u_cell, unc1, n_u - 1L) / tail_sample$m,
    K3_v_edge = fill_vector(v_cell, unc2, n_v - 1L) / tail_sample$m,
    L_mass = L_mass
  )
}

compute_B_sd <- function(operators) {
  n_u_inner <- nrow(operators$H_inner)
  n_v_inner <- ncol(operators$H_inner)
  u_prefix <- pmin(operators$u_index, n_u_inner)
  v_prefix <- pmin(operators$v_index, n_v_inner)
  loadings <- numeric(operators$m)
  inner_ok <- operators$H_inner > 0
  inv_scale <- 1 / sqrt(operators$m)

  if (any(inner_ok)) {
    weight_13 <- matrix(0, nrow = n_u_inner, ncol = n_v_inner)
    weight_13[inner_ok] <- (
      operators$K1_mass[inner_ok] / operators$H_inner[inner_ok]^2 +
        2 * operators$L_mass[inner_ok] / operators$H_inner[inner_ok]
    )
    prefix_13 <- upper_left_cumsum(weight_13)
    const_13 <- sum(operators$H_inner[inner_ok] * weight_13[inner_ok])
    loadings <- loadings + inv_scale * (prefix_13[cbind(u_prefix, v_prefix)] - const_13)

    weight_2 <- matrix(0, nrow = n_u_inner, ncol = n_v_inner)
    weight_2[inner_ok] <- 1 / operators$H_inner[inner_ok]
    value_2 <- numeric(operators$m)
    both_ok <- operators$both == 1L & operators$u_cell >= 1L & operators$u_cell <= n_u_inner &
      operators$v_cell >= 1L & operators$v_cell <= n_v_inner
    value_2[both_ok] <- weight_2[cbind(operators$u_cell[both_ok], operators$v_cell[both_ok])]
    loadings <- loadings + inv_scale * (value_2 - sum(operators$K1_mass[inner_ok] * weight_2[inner_ok]))

    weight_4 <- matrix(0, nrow = n_u_inner, ncol = n_v_inner)
    weight_4[inner_ok] <- operators$K3_v_mass[inner_ok] / operators$H_inner[inner_ok]^2
    prefix_4 <- left_cumsum(weight_4)
    value_4 <- numeric(operators$m)
    unc1_ok <- operators$unc1 == 1L & operators$u_cell >= 1L & operators$u_cell <= n_u_inner
    value_4[unc1_ok] <- prefix_4[cbind(operators$u_cell[unc1_ok], v_prefix[unc1_ok])]
    loadings <- loadings + inv_scale * (sum(operators$K2_u_mass[inner_ok] * weight_4[inner_ok]) - value_4)

    weight_5 <- matrix(0, nrow = n_u_inner, ncol = n_v_inner)
    weight_5[inner_ok] <- operators$K2_u_mass[inner_ok] / operators$H_inner[inner_ok]^2
    prefix_5 <- top_cumsum(weight_5)
    value_5 <- numeric(operators$m)
    unc2_ok <- operators$unc2 == 1L & operators$v_cell >= 1L & operators$v_cell <= n_v_inner
    value_5[unc2_ok] <- prefix_5[cbind(u_prefix[unc2_ok], operators$v_cell[unc2_ok])]
    loadings <- loadings + inv_scale * (sum(operators$K3_v_mass[inner_ok] * weight_5[inner_ok]) - value_5)
  }

  edge_u_ok <- operators$H_u_edge > 0
  if (any(edge_u_ok)) {
    weight_6 <- numeric(n_u_inner)
    weight_6[edge_u_ok] <- operators$K2_u_edge[edge_u_ok] / operators$H_u_edge[edge_u_ok]^2
    loadings <- loadings + inv_scale * (sum(operators$H_u_edge[edge_u_ok] * weight_6[edge_u_ok]) - cumsum(weight_6)[u_prefix])

    weight_7 <- numeric(n_u_inner)
    weight_7[edge_u_ok] <- 1 / operators$H_u_edge[edge_u_ok]
    value_7 <- numeric(operators$m)
    unc1_ok <- operators$unc1 == 1L & operators$u_cell >= 1L & operators$u_cell <= n_u_inner
    value_7[unc1_ok] <- weight_7[operators$u_cell[unc1_ok]]
    loadings <- loadings + inv_scale * (value_7 - sum(operators$K2_u_edge[edge_u_ok] * weight_7[edge_u_ok]))
  }

  edge_v_ok <- operators$H_v_edge > 0
  if (any(edge_v_ok)) {
    weight_8 <- numeric(n_v_inner)
    weight_8[edge_v_ok] <- operators$K3_v_edge[edge_v_ok] / operators$H_v_edge[edge_v_ok]^2
    loadings <- loadings + inv_scale * (sum(operators$H_v_edge[edge_v_ok] * weight_8[edge_v_ok]) - cumsum(weight_8)[v_prefix])

    weight_9 <- numeric(n_v_inner)
    weight_9[edge_v_ok] <- 1 / operators$H_v_edge[edge_v_ok]
    value_9 <- numeric(operators$m)
    unc2_ok <- operators$unc2 == 1L & operators$v_cell >= 1L & operators$v_cell <= n_v_inner
    value_9[unc2_ok] <- weight_9[operators$v_cell[unc2_ok]]
    loadings <- loadings + inv_scale * (value_9 - sum(operators$K3_v_edge[edge_v_ok] * weight_9[edge_v_ok]))
  }

  sqrt(sum(loadings^2))
}

# Pointwise confidence intervals and cached evaluation over k.
pointwise_ci_dabrowska <- function(
  Z1, d1, Z2, d2, k, s, t,
  alpha = 0.05,
  max_u = NULL,
  max_v = NULL,
  gamma_method = c("Hill", "EKMI")
) {
  gamma_method <- match.arg(gamma_method)
  if (s < 1 || t < 1) {
    stop("The point (s, t) must satisfy s >= 1 and t >= 1.", call. = FALSE)
  }

  tail_sample <- make_tail_sample(Z1, d1, Z2, d2, k, gamma_method = gamma_method)
  if (tail_sample$m < 5L) {
    stop("Too few observations in A_{i,n}(1,1).", call. = FALSE)
  }

  Fhat <- estimate_dabrowska_point(tail_sample, s, t)
  B_sd <- compute_B_sd(build_tail_operators(tail_sample, s, t, max_u = max_u, max_v = max_v))
  n_B <- sum(tail_sample$Y1 >= s & tail_sample$Y2 >= t & tail_sample$d1 == 1L & tail_sample$d2 == 1L)
  z_alpha <- stats::qnorm(1 - alpha / 2)

  list(
    Fhat = Fhat,
    confint = c(
      lower = Fhat / pmax(1 + z_alpha * B_sd / tail_sample$r_n, .Machine$double.eps),
      upper = Fhat / pmax(1 - z_alpha * B_sd / tail_sample$r_n, .Machine$double.eps)
    ),
    tail_sample = tail_sample,
    n_B = n_B
  )
}

make_pointwise_cache <- function(
  Z1, d1, Z2, d2, s, t,
  alpha = 0.05,
  max_u = NULL,
  max_v = NULL,
  gamma_method = c("Hill", "EKMI")
) {
  gamma_method <- match.arg(gamma_method)
  summary_cache <- new.env(parent = emptyenv())
  fit_cache <- new.env(parent = emptyenv())

  summary_at_k <- function(k) {
    key <- as.character(as.integer(k))
    if (!exists(key, envir = summary_cache, inherits = FALSE)) {
      tail_sample <- try(
        make_tail_sample(Z1, d1, Z2, d2, as.integer(k), gamma_method = gamma_method),
        silent = TRUE
      )
      if (inherits(tail_sample, "try-error")) {
        assign(key, list(tail_sample = NULL, ok = FALSE, n_B = NA_real_), envir = summary_cache)
      } else {
        assign(key, list(
          tail_sample = tail_sample,
          ok = tail_sample$m >= 5L,
          n_B = sum(tail_sample$Y1 >= s & tail_sample$Y2 >= t & tail_sample$d1 == 1L & tail_sample$d2 == 1L)
        ), envir = summary_cache)
      }
    }
    get(key, envir = summary_cache, inherits = FALSE)
  }

  fit_at_k <- function(k) {
    key <- as.character(as.integer(k))
    if (!exists(key, envir = fit_cache, inherits = FALSE)) {
      info <- summary_at_k(k)
      if (!info$ok) {
        stop("Too few observations in A_{i,n}(1,1).", call. = FALSE)
      }

      Fhat <- estimate_dabrowska_point(info$tail_sample, s, t)
      B_sd <- compute_B_sd(build_tail_operators(info$tail_sample, s, t, max_u = max_u, max_v = max_v))
      z_alpha <- stats::qnorm(1 - alpha / 2)

      assign(key, list(
        Fhat = Fhat,
        confint = c(
          lower = Fhat / pmax(1 + z_alpha * B_sd / info$tail_sample$r_n, .Machine$double.eps),
          upper = Fhat / pmax(1 - z_alpha * B_sd / info$tail_sample$r_n, .Machine$double.eps)
        ),
        tail_sample = info$tail_sample,
        n_B = info$n_B
      ), envir = fit_cache)
    }
    get(key, envir = fit_cache, inherits = FALSE)
  }

  list(summary_at_k = summary_at_k, fit_at_k = fit_at_k)
}


# -------------------------------------------------------------------------
# Pointwise path in k and adaptive selection
# -------------------------------------------------------------------------

# Threshold paths and adaptive flat-region selection.
pointwise_ci_path_dabrowska <- function(
  Z1, d1, Z2, d2, k_values, s, t,
  alpha = 0.05,
  max_u = NULL,
  max_v = NULL,
  fit_at_k = NULL,
  gamma_method = c("Hill", "EKMI")
) {
  gamma_method <- match.arg(gamma_method)
  k_values <- sort(unique(as.integer(k_values)))
  k_values <- k_values[k_values >= 1L & k_values < length(Z1)]
  if (length(k_values) == 0L) {
    stop("No admissible k values were supplied.", call. = FALSE)
  }
  path <- data.frame(
    k = k_values,
    log_k = log(k_values),
    Fhat = NA_real_,
    log_Fhat = NA_real_,
    lower = NA_real_,
    upper = NA_real_,
    m = NA_real_,
    n_B = NA_real_,
    ok = FALSE
  )

  for (i in seq_along(k_values)) {
    fit <- try(
      if (is.null(fit_at_k)) {
        pointwise_ci_dabrowska(
          Z1 = Z1,
          d1 = d1,
          Z2 = Z2,
          d2 = d2,
          k = k_values[i],
          s = s,
          t = t,
          alpha = alpha,
          max_u = max_u,
          max_v = max_v,
          gamma_method = gamma_method
        )
      } else {
        fit_at_k(k_values[i])
      },
      silent = TRUE
    )

    if (!inherits(fit, "try-error") && is.finite(fit$Fhat) && fit$Fhat > 0) {
      path$Fhat[i] <- fit$Fhat
      path$log_Fhat[i] <- log(fit$Fhat)
      path$lower[i] <- fit$confint["lower"]
      path$upper[i] <- fit$confint["upper"]
      path$m[i] <- fit$tail_sample$m
      path$n_B[i] <- fit$n_B
      path$ok[i] <- TRUE
    }
  }

  path
}

find_first_valid_k <- function(get_summary, k_lo, k_hi) {
  lo <- as.integer(k_lo)
  hi <- as.integer(k_hi)

  while (lo < hi) {
    mid <- (lo + hi) %/% 2L
    if (get_summary(mid)$ok) {
      hi <- mid
    } else {
      lo <- mid + 1L
    }
  }

  if (!get_summary(lo)$ok) {
    stop("No valid k values were available for adaptive selection.", call. = FALSE)
  }

  lo
}

closest_k_to_target <- function(get_summary, field, target, k_lo, k_hi) {
  lo <- as.integer(k_lo)
  hi <- as.integer(k_hi)

  if (get_summary(lo)[[field]] >= target) {
    return(lo)
  }
  if (get_summary(hi)[[field]] <= target) {
    return(hi)
  }

  while (lo + 1L < hi) {
    mid <- (lo + hi) %/% 2L
    if (get_summary(mid)[[field]] < target) {
      lo <- mid
    } else {
      hi <- mid
    }
  }

  candidates <- c(lo, hi)
  values <- vapply(candidates, function(k) get_summary(k)[[field]], numeric(1))
  candidates[which.min(abs(values - target))]
}

range_from_limits <- function(get_summary, field, limits, k_lower, k_upper) {
  out <- c(as.integer(k_lower), as.integer(k_upper))
  if (is.null(limits)) {
    return(out)
  }

  limits <- sort(unique(as.integer(limits[is.finite(limits)])))
  if (length(limits) == 0L) {
    return(out)
  }
  if (length(limits) == 1L) {
    limits <- c(limits, limits)
  }

  out[1L] <- closest_k_to_target(get_summary, field, limits[1L], out[1L], out[2L])
  out[2L] <- closest_k_to_target(get_summary, field, limits[2L], out[1L], out[2L])
  out
}

adaptive_k_grid_from_limits <- function(
  Z1, d1, Z2, d2, s, t,
  n_B_limits = NULL,
  n_A_limits = NULL,
  adaptive_by = 1L,
  summary_at_k = NULL,
  gamma_method = c("Hill", "EKMI")
) {
  gamma_method <- match.arg(gamma_method)
  n <- length(Z1)
  k_upper <- n - 1L
  step <- max(1L, as.integer(adaptive_by))

  get_summary <- function(k) {
    if (!is.null(summary_at_k)) {
      info <- summary_at_k(k)
      return(list(ok = info$ok, m = info$tail_sample$m, n_B = info$n_B))
    }

    tail_sample <- try(
      make_tail_sample(Z1, d1, Z2, d2, as.integer(k), gamma_method = gamma_method),
      silent = TRUE
    )
    if (inherits(tail_sample, "try-error")) {
      return(list(ok = FALSE, m = NA_real_, n_B = NA_real_))
    }
    list(
      ok = tail_sample$m >= 5L,
      m = tail_sample$m,
      n_B = sum(tail_sample$Y1 >= s & tail_sample$Y2 >= t & tail_sample$d1 == 1L & tail_sample$d2 == 1L)
    )
  }

  k_lower <- find_first_valid_k(get_summary, 2L, k_upper)
  base_range <- c(k_lower, k_upper)
  B_range <- range_from_limits(get_summary, "n_B", n_B_limits, base_range[1L], base_range[2L])
  A_range <- range_from_limits(get_summary, "m", n_A_limits, base_range[1L], base_range[2L])
  k_lower <- max(B_range[1L], A_range[1L])
  k_upper <- min(B_range[2L], A_range[2L])

  if (k_lower > k_upper) {
    stop("The n_A_limits and n_B_limits constraints do not yield an overlapping k range.", call. = FALSE)
  }

  k_grid <- seq.int(k_lower, k_upper, by = step)
  if (tail(k_grid, 1L) != k_upper) {
    k_grid <- c(k_grid, k_upper)
  }

  sort(unique(as.integer(k_grid)))
}

best_subregion_of_size <- function(curve, size) {
  size <- min(max(2L, as.integer(size)), nrow(curve))

  if (nrow(curve) <= size) {
    return(curve)
  }

  candidates <- vector("list", nrow(curve) - size + 1L)

  for (i in seq_along(candidates)) {
    j <- i + size - 1L
    x <- curve$k[i:j]
    y <- curve$log_Fhat[i:j]
    x_centered <- x - mean(x)
    y_centered <- y - mean(y)
    slope <- sum(x_centered * y_centered) / sum(x_centered^2)
    rss <- sum((y_centered - slope * x_centered)^2)

    candidates[[i]] <- data.frame(
      start = i,
      end = j,
      abs_slope = abs(slope),
      rss = rss
    )
  }

  candidates <- do.call(rbind, candidates)
  best <- candidates[order(candidates$abs_slope, candidates$rss, -candidates$end)[1L], , drop = FALSE]
  curve[best$start:best$end, , drop = FALSE]
}

select_flat_window_from_path <- function(
  k_path,
  region_sizes = c(50L, 20L, 10L),
  min_m = 0L
) {
  curve <- k_path[k_path$ok & is.finite(k_path$log_Fhat) & k_path$m >= min_m, , drop = FALSE]
  if (nrow(curve) == 0L) {
    stop("No valid k values were available for adaptive selection.", call. = FALSE)
  }
  region_sizes <- unique(as.integer(region_sizes[region_sizes >= 2L]))
  if (length(region_sizes) == 0L) {
    return(curve)
  }

  for (size in region_sizes) {
    curve <- best_subregion_of_size(curve, size)
  }

  curve
}

adaptive_pointwise_from_path <- function(
  k_path,
  region_sizes = c(50L, 20L, 10L),
  min_m = 0L
) {
  window <- select_flat_window_from_path(
    k_path,
    region_sizes = region_sizes,
    min_m = min_m
  )
  center_k <- mean(range(window$k))
  center_idx <- which.min(abs(window$k - center_k))[1L]
  selected <- window[center_idx, , drop = FALSE]

  list(
    window = window,
    k_center = selected$k,
    Fhat = selected$Fhat,
    lower = selected$lower,
    upper = selected$upper,
    m = selected$m,
    n_B = selected$n_B
  )
}


# -------------------------------------------------------------------------
# Sensitivity study for a fixed point (s, t)
# -------------------------------------------------------------------------

# Monte Carlo sensitivity study and coverage summary.
sensitivity_study_k <- function(
  n_rep = 2500,
  n,
  k_values,
  s,
  t,
  theta = 2,
  censor_props = c(0.05, 0.05),
  alpha = 0.05,
  max_u = NULL,
  max_v = NULL,
  seed = 123,
  verbose = FALSE,
  n_cores = NULL,
  progress = TRUE,
  region_sizes = c(50L, 20L, 10L),
  min_m = 0L,
  n_B_limits = NULL,
  n_A_limits = NULL,
  adaptive_by = 1L,
  gamma_method = c("Hill", "EKMI")
) {
  gamma_method <- match.arg(gamma_method)
  pars <- simulation_parameters(censor_props)
  truth <- F_true(s, t, pars$gamma_X1, pars$gamma_X2, theta)
  copula_object <- copula::gumbelCopula(theta, 2)
  k_values <- sort(unique(as.integer(k_values)))
  k_values <- k_values[k_values >= 1L & k_values < n]
  if (length(k_values) == 0L) {
    stop("No admissible k values were supplied.", call. = FALSE)
  }
  n_k <- length(k_values)

  if (is.null(n_cores)) {
    detected_cores <- parallel::detectCores(logical = FALSE)
    n_cores <- if (is.na(detected_cores) || detected_cores < 2L) 1L else detected_cores - 1L
  }

  data_generator <- function(n) {
    sim_data(
      n = n,
      theta = theta,
      tauX = pars$tauX,
      lamX = pars$lamX,
      tauC = pars$tauC,
      lamC = pars$lamC,
      copula_object = copula_object
    )
  }

  run_one_rep <- function(r) {
    set.seed(seed + r)
    dat <- data_generator(n)
    cache <- make_pointwise_cache(
      Z1 = dat$Z1,
      d1 = dat$d1,
      Z2 = dat$Z2,
      d2 = dat$d2,
      s = s,
      t = t,
      alpha = alpha,
      max_u = max_u,
      max_v = max_v,
      gamma_method = gamma_method
    )
    fixed_path <- pointwise_ci_path_dabrowska(
      Z1 = dat$Z1,
      d1 = dat$d1,
      Z2 = dat$Z2,
      d2 = dat$d2,
      k_values = k_values,
      s = s,
      t = t,
      alpha = alpha,
      max_u = max_u,
      max_v = max_v,
      fit_at_k = cache$fit_at_k,
      gamma_method = gamma_method
    )
    fixed_path$covered <- fixed_path$ok & (fixed_path$lower <= truth) & (truth <= fixed_path$upper)

    adaptive_k_values <- try(
      adaptive_k_grid_from_limits(
        Z1 = dat$Z1,
        d1 = dat$d1,
        Z2 = dat$Z2,
        d2 = dat$d2,
        s = s,
        t = t,
        n_B_limits = n_B_limits,
        n_A_limits = n_A_limits,
        adaptive_by = adaptive_by,
        summary_at_k = cache$summary_at_k,
        gamma_method = gamma_method
      ),
      silent = TRUE
    )
    if (inherits(adaptive_k_values, "try-error")) {
      return(list(failed = TRUE))
    }

    adaptive_path <- pointwise_ci_path_dabrowska(
      Z1 = dat$Z1,
      d1 = dat$d1,
      Z2 = dat$Z2,
      d2 = dat$d2,
      k_values = adaptive_k_values,
      s = s,
      t = t,
      alpha = alpha,
      max_u = max_u,
      max_v = max_v,
      fit_at_k = cache$fit_at_k,
      gamma_method = gamma_method
    )
    adaptive <- try(adaptive_pointwise_from_path(adaptive_path, region_sizes = region_sizes, min_m = min_m), silent = TRUE)
    if (inherits(adaptive, "try-error")) {
      return(list(failed = TRUE))
    }

    adaptive_covered <- adaptive$lower <= truth && truth <= adaptive$upper
    if (verbose) {
      cat(sprintf(
        "rep %02d/%02d: adaptive grid [%d, %d], window [%d, %d], Fhat = %.4f, CI = [%.4f, %.4f], covered = %s\n",
        r, n_rep, adaptive_path$k[1L], tail(adaptive_path$k, 1L),
        adaptive$window$k[1L], tail(adaptive$window$k, 1L),
        adaptive$Fhat, adaptive$lower, adaptive$upper,
        ifelse(adaptive_covered, "yes", "no")
      ))
    }

    list(
      failed = FALSE,
      ok = fixed_path$ok,
      covered = fixed_path$covered,
      Fhat = fixed_path$Fhat,
      m = fixed_path$m,
      n_B = fixed_path$n_B,
      adaptive_center = adaptive$k_center,
      adaptive_Fhat = adaptive$Fhat,
      adaptive_m = adaptive$m,
      adaptive_n_B = adaptive$n_B,
      adaptive_covered = adaptive_covered
    )
  }

  rep_out <- run_replications_with_progress(
    rep_ids = seq_len(n_rep),
    run_one_rep = run_one_rep,
    n_cores = n_cores,
    progress = progress,
    progress_label = "k-grid + adaptive"
  )

  failed <- vapply(rep_out, `[[`, logical(1), "failed")
  ok_rep <- do.call(rbind, lapply(rep_out, function(x) if (isTRUE(x$failed)) rep(FALSE, n_k) else x$ok))
  covered_rep <- do.call(rbind, lapply(rep_out, function(x) if (isTRUE(x$failed)) rep(FALSE, n_k) else x$covered))
  Fhat_rep <- do.call(rbind, lapply(rep_out, function(x) if (isTRUE(x$failed)) rep(NA_real_, n_k) else x$Fhat))
  m_rep <- do.call(rbind, lapply(rep_out, function(x) if (isTRUE(x$failed)) rep(NA_real_, n_k) else x$m))
  n_B_rep <- do.call(rbind, lapply(rep_out, function(x) if (isTRUE(x$failed)) rep(NA_real_, n_k) else x$n_B))

  k_table <- data.frame(
    label = as.character(k_values),
    k = k_values,
    mean_k_center = NA_real_,
    usable_reps = colSums(ok_rep),
    coverage = NA_real_,
    mc_halfwidth = NA_real_,
    mean_m = NA_real_,
    mean_n_B = NA_real_,
    mean_Fhat = NA_real_,
    bias = NA_real_,
    mse = NA_real_,
    truth = truth
  )

  for (j in seq_along(k_values)) {
    ok_j <- ok_rep[, j]
    if (any(ok_j)) {
      k_table$coverage[j] <- mean(covered_rep[ok_j, j])
      k_table$mean_m[j] <- mean(m_rep[ok_j, j], na.rm = TRUE)
      k_table$mean_n_B[j] <- mean(n_B_rep[ok_j, j], na.rm = TRUE)
      k_table$mean_Fhat[j] <- mean(Fhat_rep[ok_j, j], na.rm = TRUE)
      k_table$bias[j] <- k_table$mean_Fhat[j] - truth
      k_table$mse[j] <- mean((Fhat_rep[ok_j, j] - truth)^2, na.rm = TRUE)
    }
  }

  adaptive_ok <- !failed
  adaptive_center <- vapply(rep_out, function(x) if (isTRUE(x$failed)) NA_real_ else x$adaptive_center, numeric(1))
  adaptive_covered <- vapply(rep_out, function(x) if (isTRUE(x$failed)) NA else x$adaptive_covered, logical(1))
  adaptive_m <- vapply(rep_out, function(x) if (isTRUE(x$failed)) NA_real_ else x$adaptive_m, numeric(1))
  adaptive_n_B <- vapply(rep_out, function(x) if (isTRUE(x$failed)) NA_real_ else x$adaptive_n_B, numeric(1))
  adaptive_Fhat <- vapply(rep_out, function(x) if (isTRUE(x$failed)) NA_real_ else x$adaptive_Fhat, numeric(1))

  adaptive_row <- data.frame(
    label = "adaptive",
    k = NA_real_,
    mean_k_center = mean(adaptive_center[adaptive_ok], na.rm = TRUE),
    usable_reps = sum(adaptive_ok),
    coverage = mean(adaptive_covered[adaptive_ok], na.rm = TRUE),
    mc_halfwidth = NA_real_,
    mean_m = mean(adaptive_m[adaptive_ok], na.rm = TRUE),
    mean_n_B = mean(adaptive_n_B[adaptive_ok], na.rm = TRUE),
    mean_Fhat = mean(adaptive_Fhat[adaptive_ok], na.rm = TRUE),
    bias = mean(adaptive_Fhat[adaptive_ok], na.rm = TRUE) - truth,
    mse = mean((adaptive_Fhat[adaptive_ok] - truth)^2, na.rm = TRUE),
    truth = truth
  )

  results_table <- rbind(adaptive_row, k_table)
  results_table$mc_halfwidth <- 1.96 * sqrt(
    pmax(results_table$coverage * (1 - results_table$coverage) / results_table$usable_reps, 0)
  )

  list(
    settings = list(
      n_rep = n_rep,
      n = n,
      k_values = k_values,
      s = s,
      t = t,
      theta = theta,
      censor_props = censor_props,
      alpha = alpha,
      region_sizes = region_sizes,
      min_m = min_m,
      n_B_limits = n_B_limits,
      n_A_limits = n_A_limits,
      adaptive_by = adaptive_by
    ),
    truth = truth,
    results_table = results_table,
    adaptive_center = adaptive_center,
    adaptive_ok = adaptive_ok,
    adaptive_Fhat = adaptive_Fhat
  )
}

# Coverage plot used by the simulation study.
plot_sensitivity_k <- function(study, file = NULL, width = 850, height = 500, main = "") {
  draw_plot <- function() {
    x_all <- seq_len(nrow(study$results_table))
    x_fixed <- x_all[-1L]
    ylim <- c(0.7, 1.0)

    graphics::par(mar = c(5, 4, 2, 1))
    graphics::plot(
      x_fixed,
      study$results_table$coverage[-1L],
      type = "b",
      pch = 16,
      xaxt = "n",
      xlim = range(x_all),
      ylim = ylim,
      xlab = "Choice of k",
      ylab = "Coverage probability",
      main = main
    )
    graphics::axis(1, at = x_all, labels = study$results_table$label)
    graphics::abline(h = 0.95, col = "red", lwd = 2, lty = 2)
    graphics::abline(v = 1.5, col = "gray75", lty = 3)
    fixed_ok <- is.finite(study$results_table$mc_halfwidth[-1L]) & study$results_table$mc_halfwidth[-1L] > 0
    if (any(fixed_ok)) {
      graphics::arrows(
        x0 = x_fixed[fixed_ok],
        y0 = study$results_table$coverage[-1L][fixed_ok] - study$results_table$mc_halfwidth[-1L][fixed_ok],
        x1 = x_fixed[fixed_ok],
        y1 = study$results_table$coverage[-1L][fixed_ok] + study$results_table$mc_halfwidth[-1L][fixed_ok],
        angle = 90,
        code = 3,
        length = 0.05
      )
    }
    graphics::points(1, study$results_table$coverage[1L], pch = 16, col = "steelblue4", cex = 1.2)
    if (is.finite(study$results_table$mc_halfwidth[1L]) && study$results_table$mc_halfwidth[1L] > 0) {
      graphics::arrows(
        x0 = 1,
        y0 = study$results_table$coverage[1L] - study$results_table$mc_halfwidth[1L],
        x1 = 1,
        y1 = study$results_table$coverage[1L] + study$results_table$mc_halfwidth[1L],
        angle = 90,
        code = 3,
        length = 0.05,
        col = "steelblue4"
      )
    }
    graphics::legend(
      "bottomright",
      legend = c("fixed k", "adaptive rule", "0.95"),
      col = c("black", "steelblue4", "red"),
      lty = c(1, NA, 2),
      lwd = c(1, NA, 2),
      pch = c(16, 16, NA),
      bty = "n"
    )
  }

  if (is.null(file)) {
    draw_plot()
  } else {
    grDevices::png(file, width = width, height = height)
    on.exit(grDevices::dev.off(), add = TRUE)
    draw_plot()
  }

  invisible(NULL)
}


# -------------------------------------------------------------------------
# Set illustrations
# -------------------------------------------------------------------------

# Pareto model used to illustrate the four observable tail sets.
q_pareto <- function(u, gamma) {
  (1 - u)^(-gamma)
}

tc <- "#27AE60"
bk <- "#E74C3C"

sim_data_pareto <- function(n, theta, gamma_X, gamma_C) {
  require_package("copula")
  U <- copula::rCopula(n, copula::gumbelCopula(theta, 2))
  X1 <- q_pareto(U[, 1L], gamma_X[1L])
  X2 <- q_pareto(U[, 2L], gamma_X[2L])
  C1 <- q_pareto(stats::runif(n), gamma_C[1L])
  C2 <- q_pareto(stats::runif(n), gamma_C[2L])

  data.frame(
    Z1 = pmin(X1, C1),
    d1 = as.integer(X1 <= C1),
    Z2 = pmin(X2, C2),
    d2 = as.integer(X2 <= C2)
  )
}

tail_margin_fit_hill <- function(Z, d, k) {
  n <- length(Z)
  ord <- order(Z)
  sZ <- Z[ord]
  sd <- d[ord]
  Z_nk <- sZ[n - k]

  r_all <- rle(sZ)
  ends_all <- cumsum(r_all$lengths)
  event_cum_all <- cumsum(sd)
  events_all <- event_cum_all[ends_all] - c(0, event_cum_all[head(ends_all, -1L)])
  n_risk_all <- n - c(0, cumsum(head(r_all$lengths, -1L)))
  keep_last <- findInterval(Z_nk, r_all$values, rightmost.closed = TRUE)

  S_at <- if (keep_last > 0L) {
    km_factor <- 1 - events_all[seq_len(keep_last)] / n_risk_all[seq_len(keep_last)]
    prod(km_factor[events_all[seq_len(keep_last)] > 0])
  } else {
    1
  }

  block <- (n - k):n
  tail_Z <- sZ[block]
  tail_d <- sd[block]
  hill_Z <- mean(log(tail_Z[-1L] / tail_Z[1L]))
  p_uncensored <- mean(tail_d[-1L])
  gamma_hat <- hill_Z / p_uncensored

  list(
    gamma = gamma_hat,
    U = Z_nk * (S_at / (k / n))^gamma_hat
  )
}

# Display observations and censoring directions by set membership.
draw_group <- function(x, y, d1, d2, idx, col, dx, dy) {
  complete <- idx & d1 == 1L & d2 == 1L
  censored <- idx & !complete
  right_censored <- idx & d1 == 0L
  up_censored <- idx & d2 == 0L

  if (any(complete)) {
    points(x[complete], y[complete], pch = 16, cex = 0.9, col = col)
  }
  if (any(censored)) {
    points(x[censored], y[censored], pch = 1, cex = 1.0, lwd = 1.3, col = col)
  }
  if (any(right_censored)) {
    arrows(
      x0 = x[right_censored], y0 = y[right_censored],
      x1 = x[right_censored] + dx, y1 = y[right_censored],
      col = col, length = 0.08, lwd = 1.5
    )
  }
  if (any(up_censored)) {
    arrows(
      x0 = x[up_censored], y0 = y[up_censored],
      x1 = x[up_censored], y1 = y[up_censored] + dy,
      col = col, length = 0.08, lwd = 1.5
    )
  }
  if (any(censored)) {
    points(x[censored], y[censored], pch = 1, cex = 1.0, lwd = 1.3, col = col)
  }
}

plot_set_panel <- function(df, U1, U2, s, t, target, target_col) {
  x <- log(df$Z1 / U1)
  y <- log(df$Z2 / U2)
  A11 <- df$Z1 >= U1 & df$Z2 >= U2

  outside_col <- "gray90"
  Aonly_col <- bk

  dx <- 0.05 * diff(range(x))
  dy <- 0.05 * diff(range(y))
  if (!is.finite(dx) || dx == 0) dx <- 0.05
  if (!is.finite(dy) || dy == 0) dy <- 0.05

  outside_idx <- !A11
  Aonly_idx <- A11 & !target

  plot(
    x, y,
    type = "n",
    xlab = "",
    ylab = "",
    main = "",
    xlim = c(min(x), max(x) + 1.5 * dx),
    ylim = c(min(y), max(y) + 1.5 * dy),
    bty = "o",
    xaxt = "n",
    yaxt = "n"
  )

  draw_group(x, y, df$d1, df$d2, outside_idx, outside_col, dx, dy)
  draw_group(x, y, df$d1, df$d2, Aonly_idx, Aonly_col, dx, dy)
  draw_group(x, y, df$d1, df$d2, target, target_col, dx, dy)

  segments(0, 0, max(x) + 1.2 * dx, 0, lty = 2, lwd = 1.2, col = Aonly_col)
  segments(0, 0, 0, max(y) + 1.2 * dy, lty = 2, lwd = 1.2, col = Aonly_col)
  segments(log(s), log(t), max(x) + 1.2 * dx, log(t), lty = 2, lwd = 1.2, col = target_col)
  segments(log(s), log(t), log(s), max(y) + 1.2 * dy, lty = 2, lwd = 1.2, col = target_col)
  box()
}



environment()
})


# =============================================================================
# Real-data functions
# =============================================================================

real_data_functions <- local({

# Store generated figures in a subdirectory of the analysis directory.
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1L])))
} else {
  normalizePath(getwd())
}

# -------------------------------------------------------------------------
# Dabrowska estimator and adaptive threshold selection
# -------------------------------------------------------------------------

require_package <- function(pkg) {
  if (!suppressWarnings(requireNamespace(pkg, quietly = TRUE))) {
    stop(sprintf("Package '%s' is required.", pkg), call. = FALSE)
  }
}

# Censoring-adjusted marginal tail fits and tail-sample construction.
tail_margin_fit_c2 <- function(Z, d, k, gamma_method = c("Hill", "EKMI")) {
  gamma_method <- match.arg(gamma_method)
  n <- length(Z)
  ord <- order(Z)
  sZ <- Z[ord]
  sd <- d[ord]
  Z_nk <- sZ[n - k]

  r_all <- rle(sZ)
  ends_all <- cumsum(r_all$lengths)
  event_cum_all <- cumsum(sd)
  events_all <- event_cum_all[ends_all] - c(0, event_cum_all[head(ends_all, -1L)])
  n_risk_all <- n - c(0, cumsum(head(r_all$lengths, -1L)))
  keep_last <- findInterval(Z_nk, r_all$values, rightmost.closed = TRUE)

  S_at <- if (keep_last > 0L) {
    km_factor <- 1 - events_all[seq_len(keep_last)] / n_risk_all[seq_len(keep_last)]
    prod(km_factor[events_all[seq_len(keep_last)] > 0])
  } else {
    1
  }

  if (gamma_method == "Hill") {
    tail_idx <- (n - k + 1L):n
    hill_z <- mean(log(sZ[tail_idx] / Z_nk))
    p_unc <- mean(sd[tail_idx])
    if (!is.finite(hill_z) || !is.finite(p_unc) || p_unc <= 0) {
      return(list(gamma = NA_real_, U = NA_real_))
    }
    ghat <- hill_z / p_unc
  } else {
    block <- (n - k):n
    tail_Z <- sZ[block]
    tail_d <- sd[block]
    r_tail <- rle(tail_Z)
    ends_tail <- cumsum(r_tail$lengths)
    event_cum_tail <- cumsum(tail_d)
    events_tail <- event_cum_tail[ends_tail] - c(0, event_cum_tail[head(ends_tail, -1L)])
    n_risk_tail <- length(tail_Z) - c(0, cumsum(head(r_tail$lengths, -1L)))
    tail_surv <- cumprod(ifelse(events_tail > 0, 1 - events_tail / n_risk_tail, 1))
    tl <- tail(r_tail$values, k + 1L)
    sl <- tail(tail_surv, k + 1L)
    times <- tl[-1L] / tl[1L]
    wts <- diff(c(0, 1 - sl[-1L] / sl[1L]))
    ghat <- sum(log(times) * wts) / sum(wts)
  }

  list(gamma = ghat, U = Z_nk * (S_at / (k / n))^ghat)
}

make_tail_sample <- function(Z1, d1, Z2, d2, k) {
  U1 <- tail_margin_fit_c2(Z1, d1, k)$U
  U2 <- tail_margin_fit_c2(Z2, d2, k)$U
  Y1_all <- Z1 / U1
  Y2_all <- Z2 / U2
  keep <- Y1_all >= 1 & Y2_all >= 1

  list(
    Y1 = Y1_all[keep],
    Y2 = Y2_all[keep],
    d1 = as.integer(d1[keep]),
    d2 = as.integer(d2[keep]),
    m = sum(keep),
    r_n = sqrt(sum(keep))
  )
}

estimate_dabrowska_point <- function(tail_sample, s, t) {
  require_package("mhazard")

  T1 <- sort(unique(tail_sample$Y1[tail_sample$d1 == 1L & tail_sample$Y1 <= s]))
  T2 <- sort(unique(tail_sample$Y2[tail_sample$d2 == 1L & tail_sample$Y2 <= t]))
  Fhat <- mhazard:::jointHazLambda(
    tail_sample$Y1, tail_sample$Y2, T1, T2, tail_sample$d1, tail_sample$d2, 1L
  )[, , 1L]

  as.numeric(Fhat[length(T1) + 1L, length(T2) + 1L])
}

# Product-form benchmark under tail independence.
km_survival_at_point <- function(Z, d, x) {
  n <- length(Z)
  ord <- order(Z)
  sZ <- Z[ord]
  sd <- d[ord]

  rZ <- rle(sZ)
  ends <- cumsum(rZ$lengths)
  event_cum <- cumsum(sd)
  events <- event_cum[ends] - c(0, event_cum[head(ends, -1L)])
  n_risk <- n - c(0, cumsum(head(rZ$lengths, -1L)))
  keep_last <- findInterval(x, rZ$values, rightmost.closed = TRUE)

  if (keep_last <= 0L) {
    return(1)
  }

  km_factor <- 1 - events[seq_len(keep_last)] / n_risk[seq_len(keep_last)]
  prod(km_factor[events[seq_len(keep_last)] > 0])
}

estimate_independence_point <- function(tail_sample, s, t) {
  km_survival_at_point(tail_sample$Y1, tail_sample$d1, s) *
    km_survival_at_point(tail_sample$Y2, tail_sample$d2, t)
}

# Discrete operator utilities for the plug-in Gaussian variance.
thin_nodes <- function(nodes, max_nodes = NULL) {
  if (is.null(max_nodes) || length(nodes) <= max_nodes) {
    return(nodes)
  }

  interior <- nodes[-c(1L, length(nodes))]
  keep_n <- max(0L, max_nodes - 2L)
  if (length(interior) == 0L || keep_n == 0L) {
    return(c(nodes[1L], nodes[length(nodes)]))
  }

  idx <- unique(round(seq(1, length(interior), length.out = keep_n)))
  sort(unique(c(nodes[1L], interior[idx], nodes[length(nodes)])))
}

fill_matrix <- function(i, j, weights, n_row, n_col) {
  out <- matrix(0, nrow = n_row, ncol = n_col)
  keep <- i >= 1L & i <= n_row & j >= 1L & j <= n_col
  if (!any(keep)) {
    return(out)
  }

  for (p in which(keep)) {
    out[i[p], j[p]] <- out[i[p], j[p]] + weights[p]
  }

  out
}

fill_vector <- function(i, weights, n_cell) {
  out <- numeric(n_cell)
  keep <- i >= 1L & i <= n_cell
  if (!any(keep)) {
    return(out)
  }

  for (p in which(keep)) {
    out[i[p]] <- out[i[p]] + weights[p]
  }

  out
}

upper_right_cumsum <- function(x) {
  if (nrow(x) > 1L) for (i in (nrow(x) - 1L):1L) x[i, ] <- x[i, ] + x[i + 1L, ]
  if (ncol(x) > 1L) for (j in (ncol(x) - 1L):1L) x[, j] <- x[, j] + x[, j + 1L]
  x
}

upper_left_cumsum <- function(x) {
  if (nrow(x) > 1L) for (i in 2:nrow(x)) x[i, ] <- x[i, ] + x[i - 1L, ]
  if (ncol(x) > 1L) for (j in 2:ncol(x)) x[, j] <- x[, j] + x[, j - 1L]
  x
}

right_cumsum <- function(x) {
  if (ncol(x) > 1L) for (j in (ncol(x) - 1L):1L) x[, j] <- x[, j] + x[, j + 1L]
  x
}

left_cumsum <- function(x) {
  if (ncol(x) > 1L) for (j in 2:ncol(x)) x[, j] <- x[, j] + x[, j - 1L]
  x
}

down_cumsum <- function(x) {
  if (nrow(x) > 1L) for (i in (nrow(x) - 1L):1L) x[i, ] <- x[i, ] + x[i + 1L, ]
  x
}

top_cumsum <- function(x) {
  if (nrow(x) > 1L) for (i in 2:nrow(x)) x[i, ] <- x[i, ] + x[i - 1L, ]
  x
}

make_interval_index <- function(y, nodes) {
  idx <- findInterval(y, nodes, rightmost.closed = FALSE, all.inside = FALSE)
  idx[y == nodes[length(nodes)]] <- length(nodes) - 1L
  idx
}

make_threshold_index <- function(y, nodes) {
  pmax(1L, pmin(findInterval(y, nodes, rightmost.closed = TRUE), length(nodes)))
}

build_tail_operators <- function(tail_sample, s, t) {
  u <- thin_nodes(sort(unique(c(1, tail_sample$Y1[tail_sample$Y1 <= s], s))))
  v <- thin_nodes(sort(unique(c(1, tail_sample$Y2[tail_sample$Y2 <= t], t))))
  n_u <- length(u)
  n_v <- length(v)
  u_index <- make_threshold_index(tail_sample$Y1, u)
  v_index <- make_threshold_index(tail_sample$Y2, v)
  u_cell <- make_interval_index(tail_sample$Y1, u)
  v_cell <- make_interval_index(tail_sample$Y2, v)
  both <- as.numeric(tail_sample$d1 == 1L & tail_sample$d2 == 1L)
  unc1 <- as.numeric(tail_sample$d1 == 1L)
  unc2 <- as.numeric(tail_sample$d2 == 1L)

  H <- upper_right_cumsum(fill_matrix(u_index, v_index, rep(1, tail_sample$m), n_u, n_v)) / tail_sample$m
  H_inner <- H[-n_u, -n_v, drop = FALSE]
  K1_mass <- fill_matrix(u_cell, v_cell, both, n_u - 1L, n_v - 1L) / tail_sample$m
  K2_u_full <- right_cumsum(fill_matrix(u_cell, v_index, unc1, n_u - 1L, n_v)) / tail_sample$m
  K3_v_full <- down_cumsum(fill_matrix(u_index, v_cell, unc2, n_u, n_v - 1L)) / tail_sample$m
  K2_u_mass <- K2_u_full[, -n_v, drop = FALSE]
  K3_v_mass <- K3_v_full[-n_u, , drop = FALSE]

  L_mass <- matrix(0, nrow = nrow(H_inner), ncol = ncol(H_inner))
  inner_ok <- H_inner > 0
  L_mass[inner_ok] <- (
    K2_u_mass[inner_ok] * K3_v_mass[inner_ok] - H_inner[inner_ok] * K1_mass[inner_ok]
  ) / (H_inner[inner_ok]^2)

  list(
    m = tail_sample$m,
    r_n = tail_sample$r_n,
    u = u,
    v = v,
    u_index = u_index,
    v_index = v_index,
    u_cell = u_cell,
    v_cell = v_cell,
    both = both,
    unc1 = unc1,
    unc2 = unc2,
    H_inner = H_inner,
    H_u_edge = H[-n_u, 1L],
    H_v_edge = H[1L, -n_v],
    K1_mass = K1_mass,
    K2_u_mass = K2_u_mass,
    K3_v_mass = K3_v_mass,
    K2_u_edge = fill_vector(u_cell, unc1, n_u - 1L) / tail_sample$m,
    K3_v_edge = fill_vector(v_cell, unc2, n_v - 1L) / tail_sample$m,
    L_mass = L_mass
  )
}

compute_B_sd <- function(operators) {
  n_u_inner <- nrow(operators$H_inner)
  n_v_inner <- ncol(operators$H_inner)
  u_prefix <- pmin(operators$u_index, n_u_inner)
  v_prefix <- pmin(operators$v_index, n_v_inner)
  loadings <- numeric(operators$m)
  inner_ok <- operators$H_inner > 0
  inv_scale <- 1 / sqrt(operators$m)

  if (any(inner_ok)) {
    weight_13 <- matrix(0, nrow = n_u_inner, ncol = n_v_inner)
    weight_13[inner_ok] <- operators$K1_mass[inner_ok] / operators$H_inner[inner_ok]^2 +
      2 * operators$L_mass[inner_ok] / operators$H_inner[inner_ok]
    prefix_13 <- upper_left_cumsum(weight_13)
    const_13 <- sum(operators$H_inner[inner_ok] * weight_13[inner_ok])
    loadings <- loadings + inv_scale * (prefix_13[cbind(u_prefix, v_prefix)] - const_13)

    weight_2 <- matrix(0, nrow = n_u_inner, ncol = n_v_inner)
    weight_2[inner_ok] <- 1 / operators$H_inner[inner_ok]
    value_2 <- numeric(operators$m)
    both_ok <- operators$both == 1L & operators$u_cell >= 1L & operators$u_cell <= n_u_inner &
      operators$v_cell >= 1L & operators$v_cell <= n_v_inner
    value_2[both_ok] <- weight_2[cbind(operators$u_cell[both_ok], operators$v_cell[both_ok])]
    loadings <- loadings + inv_scale * (value_2 - sum(operators$K1_mass[inner_ok] * weight_2[inner_ok]))

    weight_4 <- matrix(0, nrow = n_u_inner, ncol = n_v_inner)
    weight_4[inner_ok] <- operators$K3_v_mass[inner_ok] / operators$H_inner[inner_ok]^2
    prefix_4 <- left_cumsum(weight_4)
    value_4 <- numeric(operators$m)
    unc1_ok <- operators$unc1 == 1L & operators$u_cell >= 1L & operators$u_cell <= n_u_inner
    value_4[unc1_ok] <- prefix_4[cbind(operators$u_cell[unc1_ok], v_prefix[unc1_ok])]
    loadings <- loadings + inv_scale * (sum(operators$K2_u_mass[inner_ok] * weight_4[inner_ok]) - value_4)

    weight_5 <- matrix(0, nrow = n_u_inner, ncol = n_v_inner)
    weight_5[inner_ok] <- operators$K2_u_mass[inner_ok] / operators$H_inner[inner_ok]^2
    prefix_5 <- top_cumsum(weight_5)
    value_5 <- numeric(operators$m)
    unc2_ok <- operators$unc2 == 1L & operators$v_cell >= 1L & operators$v_cell <= n_v_inner
    value_5[unc2_ok] <- prefix_5[cbind(u_prefix[unc2_ok], operators$v_cell[unc2_ok])]
    loadings <- loadings + inv_scale * (sum(operators$K3_v_mass[inner_ok] * weight_5[inner_ok]) - value_5)
  }

  edge_u_ok <- operators$H_u_edge > 0
  if (any(edge_u_ok)) {
    weight_6 <- numeric(n_u_inner)
    weight_6[edge_u_ok] <- operators$K2_u_edge[edge_u_ok] / operators$H_u_edge[edge_u_ok]^2
    loadings <- loadings + inv_scale * (sum(operators$H_u_edge[edge_u_ok] * weight_6[edge_u_ok]) - cumsum(weight_6)[u_prefix])

    weight_7 <- numeric(n_u_inner)
    weight_7[edge_u_ok] <- 1 / operators$H_u_edge[edge_u_ok]
    value_7 <- numeric(operators$m)
    unc1_ok <- operators$unc1 == 1L & operators$u_cell >= 1L & operators$u_cell <= n_u_inner
    value_7[unc1_ok] <- weight_7[operators$u_cell[unc1_ok]]
    loadings <- loadings + inv_scale * (value_7 - sum(operators$K2_u_edge[edge_u_ok] * weight_7[edge_u_ok]))
  }

  edge_v_ok <- operators$H_v_edge > 0
  if (any(edge_v_ok)) {
    weight_8 <- numeric(n_v_inner)
    weight_8[edge_v_ok] <- operators$K3_v_edge[edge_v_ok] / operators$H_v_edge[edge_v_ok]^2
    loadings <- loadings + inv_scale * (sum(operators$H_v_edge[edge_v_ok] * weight_8[edge_v_ok]) - cumsum(weight_8)[v_prefix])

    weight_9 <- numeric(n_v_inner)
    weight_9[edge_v_ok] <- 1 / operators$H_v_edge[edge_v_ok]
    value_9 <- numeric(operators$m)
    unc2_ok <- operators$unc2 == 1L & operators$v_cell >= 1L & operators$v_cell <= n_v_inner
    value_9[unc2_ok] <- weight_9[operators$v_cell[unc2_ok]]
    loadings <- loadings + inv_scale * (value_9 - sum(operators$K3_v_edge[edge_v_ok] * weight_9[edge_v_ok]))
  }

  sqrt(sum(loadings^2))
}

# Pointwise paths, confidence intervals, and cached evaluation over k.
pointwise_ci_dabrowska <- function(Z1, d1, Z2, d2, k, s, t, alpha = 0.05) {
  if (s < 1 || t < 1) {
    stop("The point (s, t) must satisfy s >= 1 and t >= 1.", call. = FALSE)
  }

  tail_sample <- make_tail_sample(Z1, d1, Z2, d2, k)
  if (tail_sample$m < 5L) {
    stop("Too few observations in A_{i,n}(1,1).", call. = FALSE)
  }

  Fhat <- estimate_dabrowska_point(tail_sample, s, t)
  B_sd <- compute_B_sd(build_tail_operators(tail_sample, s, t))
  n_B <- sum(tail_sample$Y1 >= s & tail_sample$Y2 >= t & tail_sample$d1 == 1L & tail_sample$d2 == 1L)
  z_alpha <- stats::qnorm(1 - alpha / 2)

  list(
    Fhat = Fhat,
    confint = c(
      lower = Fhat / pmax(1 + z_alpha * B_sd / tail_sample$r_n, .Machine$double.eps),
      upper = Fhat / pmax(1 - z_alpha * B_sd / tail_sample$r_n, .Machine$double.eps)
    ),
    tail_sample = tail_sample,
    n_B = n_B
  )
}

make_pointwise_cache <- function(Z1, d1, Z2, d2, s, t, alpha = 0.05) {
  summary_cache <- new.env(parent = emptyenv())
  fit_cache <- new.env(parent = emptyenv())

  summary_at_k <- function(k) {
    key <- as.character(as.integer(k))
    if (!exists(key, envir = summary_cache, inherits = FALSE)) {
      tail_sample <- make_tail_sample(Z1, d1, Z2, d2, as.integer(k))
      assign(key, list(
        tail_sample = tail_sample,
        ok = tail_sample$m >= 5L,
        n_B = sum(tail_sample$Y1 >= s & tail_sample$Y2 >= t & tail_sample$d1 == 1L & tail_sample$d2 == 1L)
      ), envir = summary_cache)
    }
    get(key, envir = summary_cache, inherits = FALSE)
  }

  fit_at_k <- function(k) {
    key <- as.character(as.integer(k))
    if (!exists(key, envir = fit_cache, inherits = FALSE)) {
      info <- summary_at_k(k)
      if (!info$ok) {
        stop("Too few observations in A_{i,n}(1,1).", call. = FALSE)
      }

      Fhat <- estimate_dabrowska_point(info$tail_sample, s, t)
      B_sd <- compute_B_sd(build_tail_operators(info$tail_sample, s, t))
      z_alpha <- stats::qnorm(1 - alpha / 2)
      assign(key, list(
        Fhat = Fhat,
        confint = c(
          lower = Fhat / pmax(1 + z_alpha * B_sd / info$tail_sample$r_n, .Machine$double.eps),
          upper = Fhat / pmax(1 - z_alpha * B_sd / info$tail_sample$r_n, .Machine$double.eps)
        ),
        tail_sample = info$tail_sample,
        n_B = info$n_B
      ), envir = fit_cache)
    }
    get(key, envir = fit_cache, inherits = FALSE)
  }

  list(summary_at_k = summary_at_k, fit_at_k = fit_at_k)
}

pointwise_ci_path_dabrowska <- function(Z1, d1, Z2, d2, k_values, s, t, fit_at_k) {
  k_values <- sort(unique(as.integer(k_values)))
  k_values <- k_values[k_values >= 1L & k_values < length(Z1)]
  if (length(k_values) == 0L) {
    stop("No admissible k values were supplied.", call. = FALSE)
  }

  path <- data.frame(
    k = k_values,
    log_Fhat = NA_real_,
    Fhat = NA_real_,
    lower = NA_real_,
    upper = NA_real_,
    m = NA_real_,
    n_B = NA_real_,
    ok = FALSE
  )

  for (i in seq_along(k_values)) {
    fit <- try(fit_at_k(k_values[i]), silent = TRUE)
    if (!inherits(fit, "try-error") && is.finite(fit$Fhat) && fit$Fhat > 0) {
      path$log_Fhat[i] <- log(fit$Fhat)
      path$Fhat[i] <- fit$Fhat
      path$lower[i] <- fit$confint["lower"]
      path$upper[i] <- fit$confint["upper"]
      path$m[i] <- fit$tail_sample$m
      path$n_B[i] <- fit$n_B
      path$ok[i] <- TRUE
    }
  }

  path
}

pointwise_path_independence <- function(Z1, d1, Z2, d2, k_values, s, t, summary_at_k) {
  k_values <- sort(unique(as.integer(k_values)))
  k_values <- k_values[k_values >= 1L & k_values < length(Z1)]
  if (length(k_values) == 0L) {
    stop("No admissible k values were supplied.", call. = FALSE)
  }

  path <- data.frame(
    k = k_values,
    log_Fhat = NA_real_,
    Fhat = NA_real_,
    lower = NA_real_,
    upper = NA_real_,
    m = NA_real_,
    n_B = NA_real_,
    ok = FALSE
  )

  for (i in seq_along(k_values)) {
    info <- try(summary_at_k(k_values[i]), silent = TRUE)
    if (inherits(info, "try-error") || !info$ok) {
      next
    }

    Fhat <- estimate_independence_point(info$tail_sample, s, t)
    if (is.finite(Fhat) && Fhat > 0) {
      path$log_Fhat[i] <- log(Fhat)
      path$Fhat[i] <- Fhat
      path$m[i] <- info$tail_sample$m
      path$n_B[i] <- info$n_B
      path$ok[i] <- TRUE
    }
  }

  path
}

# Adaptive threshold-grid and stability-window selection.
find_first_valid_k <- function(get_summary, k_lo, k_hi) {
  lo <- as.integer(k_lo)
  hi <- as.integer(k_hi)
  while (lo < hi) {
    mid <- (lo + hi) %/% 2L
    if (get_summary(mid)$ok) hi <- mid else lo <- mid + 1L
  }
  if (!get_summary(lo)$ok) {
    stop("No valid k values were available for adaptive selection.", call. = FALSE)
  }
  lo
}

closest_k_to_target <- function(get_summary, field, target, k_lo, k_hi) {
  lo <- as.integer(k_lo)
  hi <- as.integer(k_hi)
  if (get_summary(lo)[[field]] >= target) return(lo)
  if (get_summary(hi)[[field]] <= target) return(hi)

  while (lo + 1L < hi) {
    mid <- (lo + hi) %/% 2L
    if (get_summary(mid)[[field]] < target) lo <- mid else hi <- mid
  }

  candidates <- c(lo, hi)
  values <- vapply(candidates, function(k) get_summary(k)[[field]], numeric(1))
  candidates[which.min(abs(values - target))]
}

range_from_limits <- function(get_summary, field, limits, k_lower, k_upper) {
  out <- c(as.integer(k_lower), as.integer(k_upper))
  if (is.null(limits)) {
    return(out)
  }

  limits <- sort(unique(as.integer(limits[is.finite(limits)])))
  if (length(limits) == 0L) {
    return(out)
  }
  if (length(limits) == 1L) {
    limits <- c(limits, limits)
  }

  out[1L] <- closest_k_to_target(get_summary, field, limits[1L], out[1L], out[2L])
  out[2L] <- closest_k_to_target(get_summary, field, limits[2L], out[1L], out[2L])
  out
}

adaptive_k_grid_from_limits <- function(
  Z1, d1, Z2, d2, s, t,
  n_B_limits = NULL,
  n_A_limits = NULL,
  adaptive_by = 1L,
  summary_at_k
) {
  n <- length(Z1)
  k_upper <- n - 1L
  step <- max(1L, as.integer(adaptive_by))

  get_summary <- function(k) {
    info <- summary_at_k(k)
    list(ok = info$ok, m = info$tail_sample$m, n_B = info$n_B)
  }

  k_lower <- find_first_valid_k(get_summary, 2L, k_upper)
  base_range <- c(k_lower, k_upper)
  B_range <- range_from_limits(get_summary, "n_B", n_B_limits, base_range[1L], base_range[2L])
  A_range <- range_from_limits(get_summary, "m", n_A_limits, base_range[1L], base_range[2L])
  k_lower <- max(B_range[1L], A_range[1L])
  k_upper <- min(B_range[2L], A_range[2L])

  if (k_lower > k_upper) {
    stop("The n_A_limits and n_B_limits constraints do not yield an overlapping k range.", call. = FALSE)
  }

  k_grid <- seq.int(k_lower, k_upper, by = step)
  if (tail(k_grid, 1L) != k_upper) {
    k_grid <- c(k_grid, k_upper)
  }
  sort(unique(as.integer(k_grid)))
}

best_subregion_of_size <- function(curve, size) {
  size <- min(max(2L, as.integer(size)), nrow(curve))
  if (nrow(curve) <= size) {
    return(curve)
  }

  candidates <- vector("list", nrow(curve) - size + 1L)
  for (i in seq_along(candidates)) {
    j <- i + size - 1L
    x <- curve$k[i:j]
    y <- curve$log_Fhat[i:j]
    x_centered <- x - mean(x)
    y_centered <- y - mean(y)
    slope <- sum(x_centered * y_centered) / sum(x_centered^2)
    rss <- sum((y_centered - slope * x_centered)^2)
    candidates[[i]] <- data.frame(start = i, end = j, abs_slope = abs(slope), rss = rss)
  }

  candidates <- do.call(rbind, candidates)
  best <- candidates[order(candidates$abs_slope, candidates$rss, -candidates$end)[1L], , drop = FALSE]
  curve[best$start:best$end, , drop = FALSE]
}

adaptive_pointwise_from_path <- function(k_path, region_sizes = c(50L, 20L, 10L)) {
  curve <- k_path[k_path$ok & is.finite(k_path$log_Fhat), , drop = FALSE]
  if (nrow(curve) == 0L) {
    stop("No valid k values were available for adaptive selection.", call. = FALSE)
  }

  region_sizes <- unique(as.integer(region_sizes[region_sizes >= 2L]))
  for (size in region_sizes) {
    curve <- best_subregion_of_size(curve, size)
  }

  center_k <- mean(range(curve$k))
  selected <- curve[which.min(abs(curve$k - center_k))[1L], , drop = FALSE]

  list(
    window = curve,
    k_center = selected$k,
    Fhat = selected$Fhat,
    lower = selected$lower,
    upper = selected$upper,
    m = selected$m,
    n_B = selected$n_B
  )
}

# Fixed and adaptive point estimators used by the applications.
dabrowska_at_point <- function(
  dataframe,
  n_B_limits = NULL,
  n_A_limits = NULL,
  adaptive_by = 1L,
  s,
  t,
  region_sizes = c(50L, 20L, 10L),
  k = NULL
) {
  needed <- c("Z1", "d1", "Z2", "d2")
  missing <- setdiff(needed, names(dataframe))
  if (length(missing) > 0L) {
    stop(sprintf("The dataframe must contain columns %s.", paste(needed, collapse = ", ")), call. = FALSE)
  }

  Z1 <- dataframe$Z1
  d1 <- as.integer(dataframe$d1)
  Z2 <- dataframe$Z2
  d2 <- as.integer(dataframe$d2)

  if (!is.null(k)) {
    fit <- pointwise_ci_dabrowska(Z1, d1, Z2, d2, k = as.integer(k), s = s, t = t)
    return(list(
      adaptive = FALSE,
      k = as.integer(k),
      Fhat = fit$Fhat,
      confint = fit$confint,
      m = fit$tail_sample$m,
      n_B = fit$n_B
    ))
  }

  cache <- make_pointwise_cache(Z1, d1, Z2, d2, s = s, t = t)
  k_grid <- adaptive_k_grid_from_limits(
    Z1, d1, Z2, d2, s = s, t = t,
    n_B_limits = n_B_limits,
    n_A_limits = n_A_limits,
    adaptive_by = adaptive_by,
    summary_at_k = cache$summary_at_k
  )
  k_path <- pointwise_ci_path_dabrowska(Z1, d1, Z2, d2, k_values = k_grid, s = s, t = t, fit_at_k = cache$fit_at_k)
  fit <- adaptive_pointwise_from_path(k_path, region_sizes = region_sizes)

  list(
    adaptive = TRUE,
    k = fit$k_center,
    Fhat = fit$Fhat,
    confint = c(lower = fit$lower, upper = fit$upper),
    m = fit$m,
    n_B = fit$n_B,
    k_grid = k_grid,
    window_k = fit$window$k
  )
}

independence_at_point <- function(
  dataframe,
  n_B_limits = NULL,
  n_A_limits = NULL,
  adaptive_by = 1L,
  s,
  t,
  region_sizes = c(50L, 20L, 10L),
  k = NULL
) {
  needed <- c("Z1", "d1", "Z2", "d2")
  missing <- setdiff(needed, names(dataframe))
  if (length(missing) > 0L) {
    stop(sprintf("The dataframe must contain columns %s.", paste(needed, collapse = ", ")), call. = FALSE)
  }

  Z1 <- dataframe$Z1
  d1 <- as.integer(dataframe$d1)
  Z2 <- dataframe$Z2
  d2 <- as.integer(dataframe$d2)

  if (!is.null(k)) {
    tail_sample <- make_tail_sample(Z1, d1, Z2, d2, k = as.integer(k))
    return(list(
      adaptive = FALSE,
      k = as.integer(k),
      Fhat = estimate_independence_point(tail_sample, s, t),
      m = tail_sample$m,
      n_B = sum(tail_sample$Y1 >= s & tail_sample$Y2 >= t &
                  tail_sample$d1 == 1L & tail_sample$d2 == 1L)
    ))
  }

  cache <- make_pointwise_cache(Z1, d1, Z2, d2, s = s, t = t)
  k_grid <- adaptive_k_grid_from_limits(
    Z1, d1, Z2, d2, s = s, t = t,
    n_B_limits = n_B_limits,
    n_A_limits = n_A_limits,
    adaptive_by = adaptive_by,
    summary_at_k = cache$summary_at_k
  )
  k_path <- pointwise_path_independence(
    Z1, d1, Z2, d2,
    k_values = k_grid,
    s = s,
    t = t,
    summary_at_k = cache$summary_at_k
  )
  fit <- adaptive_pointwise_from_path(k_path, region_sizes = region_sizes)

  list(
    adaptive = TRUE,
    k = fit$k_center,
    Fhat = fit$Fhat,
    m = fit$m,
    n_B = fit$n_B,
    k_grid = k_grid,
    window_k = fit$window$k
  )
}


# -------------------------------------------------------------------------
# Tail-dependence coefficient estimators
# -------------------------------------------------------------------------

# Tail-quantile fits for latent, observed, and censoring components.
tail_quantile_fit_c2 <- function(Z, d, k, gamma_method = c("Hill", "EKMI")) {
  gamma_method <- match.arg(gamma_method)
  n <- length(Z)
  if (k < 2L || k >= n) {
    return(NA_real_)
  }

  ord <- order(Z)
  sZ <- Z[ord]
  sd <- as.integer(d[ord])
  if (sum(sd) == 0L) {
    return(NA_real_)
  }

  Z_nk <- sZ[n - k]

  r_all <- rle(sZ)
  ends_all <- cumsum(r_all$lengths)
  event_cum_all <- cumsum(sd)
  events_all <- event_cum_all[ends_all] - c(0, event_cum_all[head(ends_all, -1L)])
  n_risk_all <- n - c(0, cumsum(head(r_all$lengths, -1L)))
  keep_last <- findInterval(Z_nk, r_all$values, rightmost.closed = TRUE)

  S_at <- if (keep_last > 0L) {
    km_factor <- 1 - events_all[seq_len(keep_last)] / n_risk_all[seq_len(keep_last)]
    prod(km_factor[events_all[seq_len(keep_last)] > 0])
  } else {
    1
  }

  if (gamma_method == "Hill") {
    tail_idx <- (n - k + 1L):n
    hill_z <- mean(log(sZ[tail_idx] / Z_nk))
    p_unc <- mean(sd[tail_idx])
    if (!is.finite(hill_z) || !is.finite(p_unc) || p_unc <= 0) {
      return(NA_real_)
    }
    gamma_hat <- hill_z / p_unc
  } else {
    block <- (n - k):n
    tail_Z <- sZ[block]
    tail_d <- sd[block]
    r_tail <- rle(tail_Z)
    ends_tail <- cumsum(r_tail$lengths)
    event_cum_tail <- cumsum(tail_d)
    events_tail <- event_cum_tail[ends_tail] - c(0, event_cum_tail[head(ends_tail, -1L)])
    n_risk_tail <- length(tail_Z) - c(0, cumsum(head(r_tail$lengths, -1L)))
    tail_surv <- cumprod(ifelse(events_tail > 0, 1 - events_tail / n_risk_tail, 1))
    tl <- tail(r_tail$values, k + 1L)
    sl <- tail(tail_surv, k + 1L)
    times <- tl[-1L] / tl[1L]
    wts <- diff(c(0, 1 - sl[-1L] / sl[1L]))
    if (sum(wts) <= 0 || any(!is.finite(times))) {
      return(NA_real_)
    }

    gamma_hat <- sum(log(times) * wts) / sum(wts)
  }
  if (!is.finite(gamma_hat)) {
    return(NA_real_)
  }

  Z_nk * (S_at / (k / n))^gamma_hat
}

tail_quantile_fit_km <- function(Z, d, k) {
  n <- length(Z)
  if (k < 2L || k >= n) {
    return(NA_real_)
  }

  ord <- order(Z)
  sZ <- Z[ord]
  sd <- as.integer(d[ord])
  if (sum(sd) == 0L) {
    return(NA_real_)
  }

  rZ <- rle(sZ)
  ends <- cumsum(rZ$lengths)
  event_cum <- cumsum(sd)
  events <- event_cum[ends] - c(0, event_cum[head(ends, -1L)])
  n_risk <- n - c(0, cumsum(head(rZ$lengths, -1L)))
  surv <- cumprod(ifelse(events > 0, 1 - events / n_risk, 1))
  idx <- which(surv <= k / n)[1L]
  if (is.na(idx)) {
    return(NA_real_)
  }

  rZ$values[idx]
}

tail_threshold_at_k <- function(Z, d, k, threshold_method = c("evt", "km")) {
  threshold_method <- match.arg(threshold_method)
  if (threshold_method == "km") {
    return(tail_quantile_fit_km(Z, d, k))
  }

  tail_quantile_fit_c2(Z, d, k)
}

dabrowska_joint_survival_at_point <- function(Z1, d1, Z2, d2, u1, u2) {
  if (!is.finite(u1) || !is.finite(u2)) {
    return(NA_real_)
  }

  if (!requireNamespace("mhazard", quietly = TRUE)) {
    stop("Package 'mhazard' is required.", call. = FALSE)
  }

  T1 <- sort(unique(Z1[d1 == 1L & Z1 <= u1]))
  T2 <- sort(unique(Z2[d2 == 1L & Z2 <= u2]))
  S_hat <- mhazard:::jointHazLambda(
    Z1, Z2, T1, T2, as.integer(d1), as.integer(d2), 1L
  )

  if (length(dim(S_hat)) == 3L) {
    return(as.numeric(S_hat[length(T1) + 1L, length(T2) + 1L, 1L]))
  }

  as.numeric(S_hat[length(T1) + 1L, length(T2) + 1L])
}

# Product-limit tail-dependence coefficient.
tail_dependence_coefficient_at_k <- function(dataframe, k, threshold_method = c("evt", "km")) {
  needed <- c("Z1", "d1", "Z2", "d2")
  missing <- setdiff(needed, names(dataframe))
  if (length(missing) > 0L) {
    stop(sprintf("The dataframe must contain columns %s.", paste(needed, collapse = ", ")), call. = FALSE)
  }

  n <- nrow(dataframe)
  k <- as.integer(k)
  threshold_method <- match.arg(threshold_method)
  if (k < 2L || k >= n) {
    stop("The value of k must satisfy 2 <= k < n.", call. = FALSE)
  }

  u1 <- tail_threshold_at_k(dataframe$Z1, dataframe$d1, k, threshold_method)
  u2 <- tail_threshold_at_k(dataframe$Z2, dataframe$d2, k, threshold_method)
  if (!is.finite(u1) || !is.finite(u2)) {
    return(NA_real_)
  }

  S_hat <- dabrowska_joint_survival_at_point(
    dataframe$Z1,
    dataframe$d1,
    dataframe$Z2,
    dataframe$d2,
    u1,
    u2
  )
  if (!is.finite(S_hat)) {
    return(NA_real_)
  }

  n * S_hat / k
}

tail_dependence_coefficient_path <- function(dataframe, k_values, threshold_method = c("evt", "km")) {
  n <- nrow(dataframe)
  k_values <- sort(unique(as.integer(k_values)))
  k_values <- k_values[k_values >= 2L & k_values < n]
  threshold_method <- match.arg(threshold_method)
  if (length(k_values) == 0L) {
    stop("No admissible k values were supplied.", call. = FALSE)
  }

  path <- data.frame(k = k_values, chi_hat = NA_real_, ok = FALSE)

  for (i in seq_along(k_values)) {
    chi_hat <- try(
      tail_dependence_coefficient_at_k(dataframe, k_values[i], threshold_method = threshold_method),
      silent = TRUE
    )
    if (inherits(chi_hat, "try-error") || !is.finite(chi_hat)) {
      next
    }

    path$chi_hat[i] <- chi_hat
    path$ok[i] <- TRUE
  }

  path
}

# Count-based tail-dependence estimates and confidence intervals.
tail_dependence_count_at_k <- function(
  dataframe,
  k,
  threshold_method = c("evt", "km")
) {
  needed <- c("Z1", "d1", "Z2", "d2")
  missing <- setdiff(needed, names(dataframe))
  if (length(missing) > 0L) {
    stop(sprintf("The dataframe must contain columns %s.", paste(needed, collapse = ", ")), call. = FALSE)
  }

  n <- nrow(dataframe)
  k <- as.integer(k)
  threshold_method <- match.arg(threshold_method)
  if (k < 2L || k >= n) {
    stop("The value of k must satisfy 2 <= k < n.", call. = FALSE)
  }

  u1 <- tail_threshold_at_k(dataframe$Z1, dataframe$d1, k, threshold_method)
  u2 <- tail_threshold_at_k(dataframe$Z2, dataframe$d2, k, threshold_method)
  if (!is.finite(u1) || !is.finite(u2)) {
    return(NA_real_)
  }

  n12 <- sum(dataframe$Z1 > u1 & dataframe$Z2 > u2)
  n1 <- sum(dataframe$Z1 > u1)
  n2 <- sum(dataframe$Z2 > u2)

  if (n1 == 0L || n2 == 0L) {
    return(NA_real_)
  }

  k * n12 / (n1 * n2)
}

tail_dependence_count_ci_at_k <- function(
  dataframe,
  k,
  threshold_method = c("evt", "km"),
  alpha = 0.05
) {
  needed <- c("Z1", "d1", "Z2", "d2")
  missing <- setdiff(needed, names(dataframe))
  if (length(missing) > 0L) {
    stop(sprintf("The dataframe must contain columns %s.", paste(needed, collapse = ", ")), call. = FALSE)
  }

  n <- nrow(dataframe)
  k <- as.integer(k)
  threshold_method <- match.arg(threshold_method)
  if (k < 2L || k >= n) {
    stop("The value of k must satisfy 2 <= k < n.", call. = FALSE)
  }

  u1 <- tail_threshold_at_k(dataframe$Z1, dataframe$d1, k, threshold_method)
  u2 <- tail_threshold_at_k(dataframe$Z2, dataframe$d2, k, threshold_method)
  if (!is.finite(u1) || !is.finite(u2)) {
    return(list(chi_hat = NA_real_, confint = c(lower = NA_real_, upper = NA_real_), n12 = NA_real_))
  }

  n12 <- sum(dataframe$Z1 > u1 & dataframe$Z2 > u2)
  n1 <- sum(dataframe$Z1 > u1)
  n2 <- sum(dataframe$Z2 > u2)
  if (n1 == 0L || n2 == 0L) {
    return(list(chi_hat = NA_real_, confint = c(lower = NA_real_, upper = NA_real_), n12 = n12))
  }

  chi_hat <- k * n12 / (n1 * n2)
  r_hat <- sqrt(n12)
  z_alpha <- stats::qnorm(1 - alpha / 2)
  scale_hat <- if (r_hat > 0) z_alpha / r_hat else Inf

  list(
    chi_hat = chi_hat,
    confint = c(
      lower = chi_hat / pmax(1 + scale_hat, .Machine$double.eps),
      upper = chi_hat / pmax(1 - scale_hat, .Machine$double.eps)
    ),
    n12 = n12
  )
}

tail_dependence_count_path <- function(
  dataframe,
  k_values,
  threshold_method = c("evt", "km")
) {
  n <- nrow(dataframe)
  k_values <- sort(unique(as.integer(k_values)))
  k_values <- k_values[k_values >= 2L & k_values < n]
  threshold_method <- match.arg(threshold_method)
  if (length(k_values) == 0L) {
    stop("No admissible k values were supplied.", call. = FALSE)
  }

  path <- data.frame(k = k_values, chi_hat = NA_real_, ok = FALSE)

  for (i in seq_along(k_values)) {
    chi_hat <- try(
      tail_dependence_count_at_k(
        dataframe,
        k_values[i],
        threshold_method = threshold_method
      ),
      silent = TRUE
    )
    if (inherits(chi_hat, "try-error") || !is.finite(chi_hat)) {
      next
    }

    path$chi_hat[i] <- chi_hat
    path$ok[i] <- TRUE
  }

  path
}

tail_dependence_count_ci_path <- function(
  dataframe,
  k_values,
  threshold_method = c("evt", "km"),
  alpha = 0.05
) {
  n <- nrow(dataframe)
  k_values <- sort(unique(as.integer(k_values)))
  k_values <- k_values[k_values >= 2L & k_values < n]
  threshold_method <- match.arg(threshold_method)
  if (length(k_values) == 0L) {
    stop("No admissible k values were supplied.", call. = FALSE)
  }

  path <- data.frame(
    k = k_values,
    chi_hat = NA_real_,
    lower = NA_real_,
    upper = NA_real_,
    n12 = NA_real_,
    ok = FALSE
  )

  for (i in seq_along(k_values)) {
    fit <- try(
      tail_dependence_count_ci_at_k(
        dataframe,
        k_values[i],
        threshold_method = threshold_method,
        alpha = alpha
      ),
      silent = TRUE
    )
    if (inherits(fit, "try-error") || !is.finite(fit$chi_hat)) {
      next
    }

    path$chi_hat[i] <- fit$chi_hat
    path$lower[i] <- fit$confint["lower"]
    path$upper[i] <- fit$confint["upper"]
    path$n12[i] <- fit$n12
    path$ok[i] <- TRUE
  }

  path
}

# Empirical tail dependence of the observed censored vector.
tail_dependence_observed_z_at_k <- function(dataframe, k) {
  needed <- c("Z1", "Z2")
  missing <- setdiff(needed, names(dataframe))
  if (length(missing) > 0L) {
    stop(sprintf("The dataframe must contain columns %s.", paste(needed, collapse = ", ")), call. = FALSE)
  }

  n <- nrow(dataframe)
  k <- as.integer(k)
  if (k < 2L || k >= n) {
    stop("The value of k must satisfy 2 <= k < n.", call. = FALSE)
  }

  u1 <- sort(dataframe$Z1)[n - k]
  u2 <- sort(dataframe$Z2)[n - k]
  sum(dataframe$Z1 > u1 & dataframe$Z2 > u2) / k
}

tail_dependence_observed_z_ci_at_k <- function(dataframe, k, alpha = 0.05) {
  needed <- c("Z1", "Z2")
  missing <- setdiff(needed, names(dataframe))
  if (length(missing) > 0L) {
    stop(sprintf("The dataframe must contain columns %s.", paste(needed, collapse = ", ")), call. = FALSE)
  }

  n <- nrow(dataframe)
  k <- as.integer(k)
  if (k < 2L || k >= n) {
    stop("The value of k must satisfy 2 <= k < n.", call. = FALSE)
  }

  u1 <- sort(dataframe$Z1)[n - k]
  u2 <- sort(dataframe$Z2)[n - k]
  n12 <- sum(dataframe$Z1 > u1 & dataframe$Z2 > u2)
  chi_hat <- n12 / k
  r_hat <- sqrt(n12)
  z_alpha <- stats::qnorm(1 - alpha / 2)
  scale_hat <- if (r_hat > 0) z_alpha / r_hat else Inf

  list(
    chi_hat = chi_hat,
    confint = c(
      lower = chi_hat / pmax(1 + scale_hat, .Machine$double.eps),
      upper = chi_hat / pmax(1 - scale_hat, .Machine$double.eps)
    ),
    n12 = n12
  )
}

tail_dependence_observed_z_path <- function(dataframe, k_values) {
  n <- nrow(dataframe)
  k_values <- sort(unique(as.integer(k_values)))
  k_values <- k_values[k_values >= 2L & k_values < n]
  if (length(k_values) == 0L) {
    stop("No admissible k values were supplied.", call. = FALSE)
  }

  path <- data.frame(k = k_values, chi_hat = NA_real_, ok = FALSE)

  for (i in seq_along(k_values)) {
    chi_hat <- try(tail_dependence_observed_z_at_k(dataframe, k_values[i]), silent = TRUE)
    if (inherits(chi_hat, "try-error") || !is.finite(chi_hat)) {
      next
    }

    path$chi_hat[i] <- chi_hat
    path$ok[i] <- TRUE
  }

  path
}

tail_dependence_observed_z_ci_path <- function(dataframe, k_values, alpha = 0.05) {
  n <- nrow(dataframe)
  k_values <- sort(unique(as.integer(k_values)))
  k_values <- k_values[k_values >= 2L & k_values < n]
  if (length(k_values) == 0L) {
    stop("No admissible k values were supplied.", call. = FALSE)
  }

  path <- data.frame(
    k = k_values,
    chi_hat = NA_real_,
    lower = NA_real_,
    upper = NA_real_,
    n12 = NA_real_,
    ok = FALSE
  )

  for (i in seq_along(k_values)) {
    fit <- try(tail_dependence_observed_z_ci_at_k(dataframe, k_values[i], alpha = alpha), silent = TRUE)
    if (inherits(fit, "try-error") || !is.finite(fit$chi_hat)) {
      next
    }

    path$chi_hat[i] <- fit$chi_hat
    path$lower[i] <- fit$confint["lower"]
    path$upper[i] <- fit$confint["upper"]
    path$n12[i] <- fit$n12
    path$ok[i] <- TRUE
  }

  path
}


# -------------------------------------------------------------------------
# UK portfolio preparation
# -------------------------------------------------------------------------

# Construct one bivariate censored observation for each eligible UK claim.
prepare_confidential_portfolio_dataset <- function(path, cover_pair = c(3L, 35L)) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is required.", call. = FALSE)
  }

  parse_mixed_date <- function(x) {
    x <- trimws(as.character(x))
    out <- as.Date(rep(NA_character_, length(x)))
    is_num <- grepl("^[0-9]+([.][0]+)?$", x)
    out[is_num] <- as.Date(as.numeric(x[is_num]), origin = "1899-12-30")
    is_txt <- !is_num & !is.na(x) & nzchar(x)
    out[is_txt] <- as.Date(toupper(x[is_txt]), format = "%d%b%Y")
    out
  }

  # Read and chronologically order the transactional claim records.
  dat <- as.data.frame(readxl::read_excel(path), stringsAsFactors = FALSE)
  dat$from_dt <- parse_mixed_date(dat[["From Date"]])
  dat$to_dt <- parse_mixed_date(dat[["To Date"]])

  dat <- dat[dat[["Type of cover in the claims"]] %in% cover_pair, ]
  dat <- dat[order(
    dat[["Claim number"]],
    dat[["Type of cover in the claims"]],
    dat$from_dt,
    dat$to_dt,
    na.last = TRUE
  ), ]

  # Retain claim-cover histories containing a positive paid, reserve, or
  # incurred amount, then keep their most recent accounting update.
  grp <- paste(dat[["Claim number"]], dat[["Type of cover in the claims"]], sep = "_")
  active <- ave(
    (dat[["Paid"]] > 0) |
      (dat[["Remaining Case reserve"]] > 0) |
      (dat[["Incurred"]] > 0),
    grp,
    FUN = function(z) rep(any(z, na.rm = TRUE), length(z))
  )
  dat <- dat[active, ]

  grp <- paste(dat[["Claim number"]], dat[["Type of cover in the claims"]], sep = "_")
  last <- dat[!duplicated(grp, fromLast = TRUE), ]
  claims_with_both <- names(which(table(last[["Claim number"]]) == 2L))
  last <- last[last[["Claim number"]] %in% claims_with_both, ]

  # Match the two selected covers within each claim.
  left <- last[last[["Type of cover in the claims"]] == cover_pair[1L], ]
  right <- last[last[["Type of cover in the claims"]] == cover_pair[2L], ]
  left <- left[order(left[["Claim number"]]), ]
  right <- right[order(right[["Claim number"]]), ]

  stopifnot(identical(left[["Claim number"]], right[["Claim number"]]))

  # Paid amounts are the observed values. A zero remaining case reserve marks
  # a settled, and therefore uncensored, claim component.
  out <- data.frame(
    Z1 = left[["Paid"]],
    d1 = as.integer(left[["Remaining Case reserve"]] == 0),
    Z2 = right[["Paid"]],
    d2 = as.integer(right[["Remaining Case reserve"]] == 0)
  )
  out <- out[complete.cases(out), ]
  out <- out[out$Z1 > 0 & out$Z2 > 0, ]

  list(
    dataframe = out,
    cover_pair = as.integer(cover_pair),
    claim_count = nrow(out)
  )
}


# -------------------------------------------------------------------------
# Tail diagnostics and application figures
# -------------------------------------------------------------------------

# Consistent plotting palette used by both data applications.
flat_cols <- list(
  observed = "#2C3E50",
  censored = "#7F8C8D",
  right = "#E67E22",
  up = "#16A085",
  blue = "#2980B9",
  red = "#C0392B",
  purple = "#8E44AD",
  grid = "#ECF0F1",
  band = grDevices::adjustcolor("#2980B9", alpha.f = 0.16),
  band_line = grDevices::adjustcolor("#2980B9", alpha.f = 0.45),
  guide = grDevices::adjustcolor("#424949", alpha.f = 1)
)

figure_dir <- function() {
  out <- file.path(script_dir, "figures")
  if (!dir.exists(out)) {
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
  }
  normalizePath(out, mustWork = TRUE)
}

save_square_pdf <- function(filename, plot_fun, width = 5, height = 5, mar = c(4.1, 5.1, 0.8, 0.8)) {
  grDevices::pdf(
    file.path(figure_dir(), filename),
    width = width,
    height = height,
    useDingbats = FALSE
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = mar, mgp = c(2.6, 0.8, 0), tcl = -0.25)
  plot_fun()
}

# Marginal tail diagnostics.
hill_margin_path <- function(Z, d, k_values) {
  n <- length(Z)
  k_values <- sort(unique(as.integer(k_values)))
  k_values <- k_values[k_values >= 2L & k_values < n]
  if (length(k_values) == 0L) {
    stop("No admissible k values were supplied.", call. = FALSE)
  }

  ord <- order(Z)
  sZ <- Z[ord]
  sd <- as.integer(d[ord])
  logZ <- log(sZ)
  tail_log_sums <- cumsum(rev(logZ))
  tail_unc_sums <- cumsum(rev(sd))

  hill_z <- tail_log_sums[k_values] / k_values - log(sZ[n - k_values])
  p_unc <- tail_unc_sums[k_values] / k_values
  gamma_hat <- hill_z / p_unc
  gamma_hat[!is.finite(gamma_hat) | p_unc <= 0] <- NA_real_

  data.frame(
    k = k_values,
    gamma_hat = gamma_hat,
    p_unc = p_unc,
    ok = is.finite(gamma_hat) & is.finite(p_unc) & p_unc > 0
  )
}

pareto_km_points <- function(Z, d, k_max = NULL) {
  n <- length(Z)
  k_max <- if (is.null(k_max)) n - 1L else min(as.integer(k_max), n - 1L)
  ord <- order(Z)
  sZ <- Z[ord]
  sd <- as.integer(d[ord])

  rZ <- rle(sZ)
  ends <- cumsum(rZ$lengths)
  event_cum <- cumsum(sd)
  events <- event_cum[ends] - c(0, event_cum[head(ends, -1L)])
  n_risk <- n - c(0, cumsum(head(rZ$lengths, -1L)))
  surv <- cumprod(ifelse(events > 0, 1 - events / n_risk, 1))
  S_at <- surv[findInterval(sZ, rZ$values, rightmost.closed = TRUE)]

  idx <- (n - k_max + 1L):n
  out <- data.frame(x = log(sZ[idx]), y = -log(S_at[idx]), d = sd[idx])
  out[is.finite(out$x) & is.finite(out$y), ]
}

plot_pareto_qq_pair <- function(dataframe) {
  p1 <- pareto_km_points(dataframe$Z1, dataframe$d1)
  p2 <- pareto_km_points(dataframe$Z2, dataframe$d2)
  xlim <- range(c(p1$x, p2$x), finite = TRUE)
  x_pad <- 0.04 * diff(xlim)
  if (!is.finite(x_pad) || x_pad == 0) x_pad <- 0.1
  xlim <- xlim + c(-x_pad, x_pad)
  ylim <- c(0, 1.04 * max(c(p1$y, p2$y), na.rm = TRUE))

  plot(
    xlim,
    ylim,
    type = "n",
    xlab = "log(Z)",
    ylab = "-log(KM)",
    bty = "l"
  )
  grid(col = flat_cols$grid, lty = 1)
  points(p1$x, p1$y, pch = 1, cex = 0.72, lwd = 1.1, col = flat_cols$blue)
  points(p2$x, p2$y, pch = 2, cex = 0.72, lwd = 1.1, col = flat_cols$red)
}

# Fixed- and adaptive-threshold estimates used in the figures.
anchor_k_fit <- function(
  dataframe,
  n_B_limits = c(1L, 40L),
  n_A_limits = c(1L, 40L),
  adaptive_by = 1L,
  region_sizes = c(50L, 20L, 10L),
  anchor = c(1.5, 1.5)
) {
  dabrowska_at_point(
    dataframe = dataframe,
    n_B_limits = n_B_limits,
    n_A_limits = n_A_limits,
    adaptive_by = adaptive_by,
    s = anchor[1L],
    t = anchor[2L],
    region_sizes = region_sizes,
    k = NULL
  )
}

fixed_k_tail_sample <- function(dataframe, k) {
  make_tail_sample(
    dataframe$Z1,
    as.integer(dataframe$d1),
    dataframe$Z2,
    as.integer(dataframe$d2),
    k = as.integer(k)
  )
}

fixed_k_curve <- function(dataframe, k, s_grid, t_values, alpha = 0.05) {
  tail_sample <- fixed_k_tail_sample(dataframe, k)
  z_alpha <- stats::qnorm(1 - alpha / 2)
  Fhat <- lower <- upper <- indep <- matrix(
    NA_real_,
    nrow = length(s_grid),
    ncol = length(t_values)
  )

  for (j in seq_along(t_values)) {
    for (i in seq_along(s_grid)) {
      F_ij <- estimate_dabrowska_point(tail_sample, s_grid[i], t_values[j])
      B_sd <- compute_B_sd(build_tail_operators(tail_sample, s_grid[i], t_values[j]))
      scale <- z_alpha * B_sd / tail_sample$r_n
      Fhat[i, j] <- F_ij
      lower[i, j] <- F_ij / pmax(1 + scale, .Machine$double.eps)
      upper[i, j] <- F_ij / pmax(1 - scale, .Machine$double.eps)
      indep[i, j] <- estimate_independence_point(tail_sample, s_grid[i], t_values[j])
    }
  }

  list(
    s_grid = s_grid,
    t_values = t_values,
    Fhat = Fhat,
    lower = lower,
    upper = upper,
    independence = indep
  )
}

fixed_k_surface <- function(dataframe, k, s_grid, t_grid) {
  tail_sample <- fixed_k_tail_sample(dataframe, k)
  zmat <- matrix(NA_real_, nrow = length(s_grid), ncol = length(t_grid))

  for (i in seq_along(s_grid)) {
    for (j in seq_along(t_grid)) {
      zmat[i, j] <- estimate_dabrowska_point(tail_sample, s_grid[i], t_grid[j])
    }
  }

  list(s_grid = s_grid, t_grid = t_grid, z = zmat)
}

# Plotting functions for censored observations and tail diagnostics.
plot_scatter_with_arrows <- function(dataframe) {
  x <- log(dataframe$Z1)
  y <- log(dataframe$Z2)
  full <- dataframe$d1 == 1L & dataframe$d2 == 1L
  any_censored <- !full
  right_censored <- dataframe$d1 == 0L
  up_censored <- dataframe$d2 == 0L
  dx <- 0.05 * diff(range(x))
  dy <- 0.05 * diff(range(y))

  plot(
    x, y,
    type = "n",
    xlab = expression(log(Z[1])),
    ylab = expression(log(Z[2])),
    xlim = c(min(x) - 0.4 * dx, max(x) + 1.8 * dx),
    ylim = c(min(y) - 0.4 * dy, max(y) + 1.8 * dy),
    bty = "l"
  )
  grid(col = flat_cols$grid, lty = 1)
  points(x[full], y[full], pch = 16, cex = 0.72, col = flat_cols$observed)
  points(x[any_censored], y[any_censored], pch = 1, cex = 0.88, lwd = 1.1, col = flat_cols$censored)

  if (any(right_censored)) {
    arrows(
      x0 = x[right_censored], y0 = y[right_censored],
      x1 = x[right_censored] + dx, y1 = y[right_censored],
      col = flat_cols$right, length = 0.08, lwd = 1.5
    )
  }
  if (any(up_censored)) {
    arrows(
      x0 = x[up_censored], y0 = y[up_censored],
      x1 = x[up_censored], y1 = y[up_censored] + dy,
      col = flat_cols$up, length = 0.08, lwd = 1.5
    )
  }
  points(x[any_censored], y[any_censored], pch = 1, cex = 0.88, lwd = 1.1, col = flat_cols$censored)
}

plot_gamma_path <- function(path1, path2, k_star) {
  ok1 <- path1$ok & is.finite(path1$gamma_hat)
  ok2 <- path2$ok & is.finite(path2$gamma_hat)
  y <- c(path1$gamma_hat[ok1], path2$gamma_hat[ok2])

  plot(
    range(c(path1$k[ok1], path2$k[ok2])),
    range(y),
    type = "n",
    xlab = "k",
    ylab = expression(hat(gamma)[j](k)),
    bty = "l"
  )
  grid(col = flat_cols$grid, lty = 1)
  abline(v = k_star, col = flat_cols$guide, lty = 3, lwd = 2)
  lines(path1$k, path1$gamma_hat, col = flat_cols$blue, lwd = 2)
  lines(path2$k, path2$gamma_hat, col = flat_cols$red, lwd = 2, lty = 2)
}

plot_noncensoring_path <- function(path1, path2, k_star) {
  plot(
    range(c(path1$k, path2$k)),
    c(0, 1),
    type = "n",
    xlab = "k",
    ylab = expression(hat(p)[j](k)),
    bty = "l"
  )
  grid(col = flat_cols$grid, lty = 1)
  lines(path1$k, path1$p_unc, col = flat_cols$blue, lwd = 2)
  lines(path2$k, path2$p_unc, col = flat_cols$red, lwd = 2, lty = 2)
}

plot_surface <- function(surface) {
  persp(
    x = surface$s_grid,
    y = surface$t_grid,
    z = surface$z,
    theta = 35,
    phi = 28,
    expand = 0.65,
    shade = 0.35,
    ticktype = "detailed",
    col = "#AED6F1",
    border = "#5DADE2",
    xlab = "s",
    ylab = "t",
    zlab = "F_D"
  )
}

plot_fixed_t_curves <- function(curves) {
  y <- c(curves$Fhat, curves$lower, curves$upper, curves$independence)
  y <- y[is.finite(y)]
  pad <- 0.04 * diff(range(y))
  if (!is.finite(pad) || pad == 0) {
    pad <- 0.02
  }

  plot(
    range(curves$s_grid),
    range(y) + c(-pad, pad),
    type = "n",
    xlab = "s",
    ylab = expression(hat(F)[D]),
    bty = "l"
  )
  grid(col = flat_cols$grid, lty = 1)

  for (j in seq_along(curves$t_values)) {
    ok <- is.finite(curves$lower[, j]) & is.finite(curves$upper[, j])
    if (any(ok)) {
      polygon(
        c(curves$s_grid[ok], rev(curves$s_grid[ok])),
        c(curves$upper[ok, j], rev(curves$lower[ok, j])),
        col = flat_cols$band,
        border = NA
      )
      lines(curves$s_grid, curves$lower[, j], col = flat_cols$band_line, lwd = 0.9, lty = j)
      lines(curves$s_grid, curves$upper[, j], col = flat_cols$band_line, lwd = 0.9, lty = j)
    }
    lines(curves$s_grid, curves$Fhat[, j], col = flat_cols$blue, lwd = 2, lty = j)
    lines(curves$s_grid, curves$independence[, j], col = flat_cols$red, lwd = 2, lty = j)
  }
}

plot_chi_paths <- function(chi_z, chi_x, chi_c = NULL, k_star = NULL) {
  plot(
    range(chi_x$k),
    c(0, 1),
    type = "n",
    xlab = "k",
    ylab = expression(hat(chi)[X](k)),
    bty = "l"
  )
  grid(col = flat_cols$grid, lty = 1)
  lines(chi_x$k, chi_x$chi_hat, col = flat_cols$blue, lwd = 2)
}

# Compute the threshold diagnostics, tail estimates, and application figures.
generate_real_data_figures <- function(
  dataframe,
  prefix,
  k_values,
  curve_k = NULL,
  curve_t_values,
  s_grid,
  surface_s_grid,
  surface_t_grid,
  n_B_limits = c(1L, 40L),
  n_A_limits = c(1L, 40L),
  adaptive_by = 1L,
  region_sizes = c(50L, 20L, 10L),
  show_c_curve = TRUE
) {
  # Select k at the anchor point (1.5, 1.5).
  anchor_fit <- anchor_k_fit(
    dataframe = dataframe,
    n_B_limits = n_B_limits,
    n_A_limits = n_A_limits,
    adaptive_by = adaptive_by,
    region_sizes = region_sizes
  )
  k_star <- as.integer(anchor_fit$k)
  curve_k <- if (is.null(curve_k)) k_star else as.integer(curve_k)

  # Marginal diagnostics, sections of F_D, its surface, and tail dependence.
  gamma1 <- hill_margin_path(dataframe$Z1, dataframe$d1, k_values)
  gamma2 <- hill_margin_path(dataframe$Z2, dataframe$d2, k_values)
  curves <- fixed_k_curve(dataframe, curve_k, s_grid = s_grid, t_values = curve_t_values)
  surface <- fixed_k_surface(dataframe, k_star, s_grid = surface_s_grid, t_grid = surface_t_grid)
  chi_z <- tail_dependence_observed_z_path(dataframe, k_values)
  chi_x <- tail_dependence_count_path(dataframe, k_values, threshold_method = "evt")
  chi_c <- NULL

  if (show_c_curve) {
    censoring_data <- data.frame(
      Z1 = dataframe$Z1,
      d1 = 1L - as.integer(dataframe$d1),
      Z2 = dataframe$Z2,
      d2 = 1L - as.integer(dataframe$d2)
    )
    chi_c <- tail_dependence_count_path(censoring_data, k_values, threshold_method = "km")
  }

  # Write the seven figures used for each application.
  save_square_pdf(
    paste0(prefix, "_scatter.pdf"),
    function() plot_scatter_with_arrows(dataframe),
    mar = c(4.2, 5.3, 1.0, 1.0)
  )
  save_square_pdf(
    paste0(prefix, "_gamma_hill.pdf"),
    function() plot_gamma_path(gamma1, gamma2, k_star)
  )
  save_square_pdf(
    paste0(prefix, "_pareto_qq.pdf"),
    function() plot_pareto_qq_pair(dataframe),
    mar = c(4.2, 5.3, 1.0, 1.0)
  )
  save_square_pdf(
    paste0(prefix, "_noncensoring.pdf"),
    function() plot_noncensoring_path(gamma1, gamma2, k_star)
  )
  save_square_pdf(
    paste0(prefix, "_surface.pdf"),
    function() plot_surface(surface)
  )
  save_square_pdf(
    paste0(prefix, "_curves.pdf"),
    function() plot_fixed_t_curves(curves)
  )
  save_square_pdf(
    paste0(prefix, "_chi_paths.pdf"),
    function() plot_chi_paths(chi_z, chi_x, chi_c, k_star)
  )

  list(
    k_star = k_star,
    curve_k = curve_k,
    anchor_window = anchor_fit$window_k,
    n = nrow(dataframe),
    curve_t_values = curve_t_values
  )
}


environment()
})
