############################
## Quantile-Only Regression: Inference and Partial Identification under Quantile Disclosure
## Political Control and Children in Need Assessment Duration
############################
# Primary design:
#   - all matched council-by-year observations;
#   - running variable = Labour seat share in percentage points minus 50;
#   - fixed windows of +/-3, +/-4, and +/-5 percentage points;
#   - +/-4 is the main window; +/-3 and +/-5 are sensitivity analyses;
#   - X = (intercept, Labour-at-least-50).
# A supplemental sensitivity specification adds reporting-year fixed effects while
# retaining the same fixed windows, assessment-count weights, and quantile
# reconstruction.
#
# The script reports (i) qlm interpolation estimates and (ii) the simultaneous
# confidence region constructed by qlm_region under explicit support-upper-
# bound sensitivity assumptions.


project_root <- normalizePath("..", mustWork = TRUE)
input_path <- file.path(
  project_root,
  "data", "derived", "children_in_need_labour_rd_all_years.csv"
)
output_dir <- file.path(project_root, "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_path)) {
  stop("Analysis data are missing. Run R/01_build_analysis_data.R first.")
}
if (!requireNamespace("qonlyreg", quietly = TRUE)) {
  stop("qonlyreg is not installed. Install it from the fixed commit specified in README.md.")
}
if (!requireNamespace("lpSolve", quietly = TRUE)) {
  stop("lpSolve, which is required to project the simultaneous confidence region onto coefficient directions, is not installed.")
}

# No finite administrative maximum is documented for assessment duration.
# These are sensitivity assumptions, not facts about the data-generating
# process.  Every value is larger than the largest disclosed Q3 (117 days).
upper_bounds <- c(150, 175, 200)
simulation_draws <- 10000L

data <- read.csv(input_path, stringsAsFactors = FALSE)
logical_columns <- c(
  "analysis_eligible", "labour_at_least_50", "labour_strict_majority",
  "labour_exactly_50"
)
for (nm in logical_columns) data[[nm]] <- as.logical(data[[nm]])

tau_grid <- c(0.25, 0.50, 0.75)
bandwidths_pp <- c(3.0, 4.0, 5.0)
bandwidth_role <- function(bw) {
  if (bw == 4.0) return("main")
  if (bw < 4.0) return("narrower_window_sensitivity")
  "wider_window_sensitivity"
}
analysis_specifications <- data.frame(
  treatment_definition = c(
    "at_least_50_main",
    "strictly_above_50_sensitivity",
    "at_least_50_main"
  ),
  treatment_column = c(
    "labour_at_least_50",
    "labour_strict_majority",
    "labour_at_least_50"
  ),
  model_specification = c(
    "no_year_fixed_effects",
    "no_year_fixed_effects",
    "reporting_year_fixed_effects_sensitivity"
  ),
  reporting_year_fixed_effects = c(FALSE, FALSE, TRUE),
  stringsAsFactors = FALSE
)

selected_data <- function(bandwidth_pp) {
  keep <- data$analysis_eligible %in% TRUE &
    !is.na(data$running_variable_pp) &
    abs(data$running_variable_pp) <= bandwidth_pp + 1e-12
  d <- data[keep, , drop = FALSE]
  d[order(d$reporting_year, d$new_la_code), , drop = FALSE]
}

design_matrix <- function(d, treatment_column, reporting_year_fixed_effects) {
  treatment <- as.integer(d[[treatment_column]])
  X <- cbind(`(Intercept)` = 1, labour_majority = treatment)
  if (reporting_year_fixed_effects) {
    reporting_year <- factor(d$reporting_year)
    year_matrix <- stats::model.matrix(~reporting_year)[, -1, drop = FALSE]
    colnames(year_matrix) <- paste0(
      "reporting_year_", levels(reporting_year)[-1]
    )
    X <- cbind(X, year_matrix)
  }
  if (qr(X)$rank < ncol(X)) {
    stop(
      "A fixed-window design matrix is rank deficient. Check threshold-side ",
      "and reporting-year support."
    )
  }
  X
}

# qonlyreg point estimator with the known deterministic lower
# support endpoint included in the interpolation grid.  
# The upper support endpoint is unknown, so the quantile function remains constant above the largest
# disclosed quantile.  Known endpoints have zero sampling variance.
qlm_known_lower <- function(
    qhat, tau, X, n_g, y_lower = 0, density_floor = 1e-8) {
  qhat <- as.matrix(qhat)
  X <- as.matrix(X)
  tau <- as.numeric(tau)
  n_g <- as.numeric(n_g)
  G <- nrow(qhat)
  K <- ncol(qhat)
  
  stopifnot(
    length(tau) == K,
    length(n_g) == G,
    nrow(X) == G,
    K >= 2,
    all(diff(tau) > 0),
    all(tau > 0 & tau < 1),
    all(apply(qhat, 1, function(z) all(diff(z) >= 0))),
    all(qhat[, 1] >= y_lower)
  )
  
  # Trapezoidal weights on the augmented lower grid, followed by 
  # constant-tail rule above the last disclosed quantile.
  alpha <- numeric(K)
  alpha[1] <- tau[1] / 2 + (tau[2] - tau[1]) / 2
  if (K >= 3) {
    for (k in 2:(K - 1)) {
      alpha[k] <- (tau[k] - tau[k - 1]) / 2 +
        (tau[k + 1] - tau[k]) / 2
    }
  }
  alpha[K] <- (tau[K] - tau[K - 1]) / 2 + (1 - tau[K])
  lower_endpoint_weight <- tau[1] / 2
  if (abs(lower_endpoint_weight + sum(alpha) - 1) > 1e-12) {
    stop("Known-lower-support interpolation weights do not sum to one.")
  }
  
  m_hat <- as.numeric(y_lower * lower_endpoint_weight + qhat %*% alpha)
  cov_tau <- outer(tau, tau, function(a, b) pmin(a, b) - a * b)
  var_m_hat <- numeric(G)
  fhat <- matrix(NA_real_, nrow = G, ncol = K)
  
  for (g in seq_len(G)) {
    qprime <- numeric(K)
    # Use the known lower endpoint when estimating the slope at Q1.
    qprime[1] <- (qhat[g, 2] - y_lower) / (tau[2] - 0)
    qprime[K] <- (qhat[g, K] - qhat[g, K - 1]) /
      (tau[K] - tau[K - 1])
    if (K >= 3) {
      for (k in 2:(K - 1)) {
        qprime[k] <- (qhat[g, k + 1] - qhat[g, k - 1]) /
          (tau[k + 1] - tau[k - 1])
      }
    }
    qprime <- pmax(qprime, density_floor)
    fhat[g, ] <- 1 / qprime
    quantile_cov <- cov_tau / outer(fhat[g, ], fhat[g, ], "*") / n_g[g]
    var_m_hat[g] <- as.numeric(t(alpha) %*% quantile_cov %*% alpha)
  }
  
  W <- diag(n_g, nrow = G)
  A <- crossprod(X, W %*% X)
  beta_hat <- solve(A, crossprod(X, W %*% m_hat))
  meat <- matrix(0, nrow = ncol(X), ncol = ncol(X))
  for (g in seq_len(G)) {
    xg <- matrix(X[g, ], ncol = 1)
    meat <- meat + n_g[g]^2 * var_m_hat[g] * (xg %*% t(xg))
  }
  A_inv <- solve(A)
  vcov_beta <- A_inv %*% meat %*% A_inv
  se_beta <- sqrt(diag(vcov_beta))
  if (!is.null(colnames(X))) {
    names(beta_hat) <- colnames(X)
    names(se_beta) <- colnames(X)
    dimnames(vcov_beta) <- list(colnames(X), colnames(X))
  }
  
  list(
    coefficients = beta_hat,
    se = se_beta,
    vcov = vcov_beta,
    m_hat = m_hat,
    var_m_hat = var_m_hat,
    quantile_weights = alpha,
    lower_endpoint_weight = lower_endpoint_weight,
    fhat = fhat
  )
}

# Project a polyhedron A beta (>= or <=) b onto each coefficient.  lpSolve
# requires nonnegative variables, so each unrestricted coefficient is written
# as beta_plus - beta_minus.
project_linear_region <- function(const_mat, const_dir, const_rhs) {
  const_mat <- as.matrix(const_mat)
  p <- ncol(const_mat)
  free_mat <- cbind(const_mat, -const_mat)
  zero_objective <- rep(0, 2 * p)
  feasible <- lpSolve::lp(
    direction = "min",
    objective.in = zero_objective,
    const.mat = free_mat,
    const.dir = const_dir,
    const.rhs = const_rhs
  )
  if (feasible$status != 0) {
    return(list(status = feasible$status, bounds = matrix(
      NA_real_, nrow = p, ncol = 2,
      dimnames = list(colnames(const_mat), c("lower", "upper"))
    )))
  }
  
  out <- matrix(
    NA_real_, nrow = p, ncol = 2,
    dimnames = list(colnames(const_mat), c("lower", "upper"))
  )
  for (j in seq_len(p)) {
    objective <- rep(0, 2 * p)
    objective[j] <- 1
    objective[p + j] <- -1
    lower <- lpSolve::lp(
      direction = "min", objective.in = objective,
      const.mat = free_mat, const.dir = const_dir, const.rhs = const_rhs
    )
    upper <- lpSolve::lp(
      direction = "max", objective.in = objective,
      const.mat = free_mat, const.dir = const_dir, const.rhs = const_rhs
    )
    if (lower$status == 0 && upper$status == 0) {
      out[j, ] <- c(lower$objval, upper$objval)
    }
  }
  list(status = 0L, bounds = out)
}

fit_point <- function(
    bandwidth_pp, treatment_name, treatment_column, model_specification,
    reporting_year_fixed_effects) {
  d <- selected_data(bandwidth_pp)
  X <- design_matrix(d, treatment_column, reporting_year_fixed_effects)
  qhat <- as.matrix(d[c("q25_days", "median_days", "q75_days")])
  fit <- qlm_known_lower(
    qhat = qhat,
    tau = tau_grid,
    X = X,
    n_g = d$assessments_completed,
    y_lower = 0
  )
  if (max(abs(fit$m_hat - d$qonlyreg_mean_days)) > 1e-10) {
    stop("The qlm reconstructed mean values do not match the audit column.")
  }
  
  treatment <- X[, "labour_majority"]
  data.frame(
    bandwidth_pp = bandwidth_pp,
    bandwidth_role = bandwidth_role(bandwidth_pp),
    treatment_definition = treatment_name,
    model_specification = model_specification,
    reporting_year_fixed_effects = reporting_year_fixed_effects,
    specification = if (reporting_year_fixed_effects) {
      "fixed_window_intercept_plus_treatment_plus_reporting_year_fixed_effects_no_kernel"
    } else {
      "fixed_window_intercept_plus_treatment_no_kernel"
    },
    observations = nrow(d),
    unique_councils = length(unique(d$new_la_code)),
    treated = sum(treatment == 1),
    control = sum(treatment == 0),
    exactly_50 = sum(d$labour_exactly_50),
    estimate_days = unname(fit$coefficients["labour_majority"]),
    qlm_standard_error = unname(fit$se["labour_majority"]),
    known_support_lower_days = 0,
    lower_endpoint_weight = fit$lower_endpoint_weight,
    q25_weight = fit$quantile_weights[1],
    median_weight = fit$quantile_weights[2],
    q75_weight = fit$quantile_weights[3],
    point_estimate_role = "qonlyreg point estimator with known lower support endpoint",
    assessment_count_weighted = TRUE,
    qonlyreg_version = as.character(utils::packageVersion("qonlyreg")),
    stringsAsFactors = FALSE
  )
}

point_results <- do.call(rbind, lapply(bandwidths_pp, function(bw) {
  do.call(rbind, lapply(seq_len(nrow(analysis_specifications)), function(j) {
    spec <- analysis_specifications[j, ]
    fit_point(
      bandwidth_pp = bw,
      treatment_name = spec$treatment_definition,
      treatment_column = spec$treatment_column,
      model_specification = spec$model_specification,
      reporting_year_fixed_effects = spec$reporting_year_fixed_effects
    )
  }))
}))
write.csv(
  point_results,
  file.path(output_dir, "qonlyreg_labour_fixed_window_point_estimates.csv"),
  row.names = FALSE,
  na = ""
)

# Each assumed upper support bound must be at least as large as the largest disclosed third quartile. 
# This is only a compatibility check; the bounds are sensitivity assumptions, not documented maxima or
# exact support endpoints.
eligible <- data[data$analysis_eligible %in% TRUE, , drop = FALSE]
qhat_all <- as.matrix(eligible[c("q25_days", "median_days", "q75_days")])
max_disclosed_q3 <- max(qhat_all[, 3])
if (any(upper_bounds < max_disclosed_q3)) {
  stop(
    "Every support upper bound must be at least the largest disclosed Q3 (",
    max_disclosed_q3, ")."
  )
}

simultaneous_region_one <- function(
    bandwidth_pp,
    treatment_name,
    treatment_column,
    model_specification,
    reporting_year_fixed_effects,
    y_upper_bound,
    gaussian_seed) {
  
  d <- selected_data(bandwidth_pp)
  X <- design_matrix(d,treatment_column, reporting_year_fixed_effects)
  qhat <- as.matrix(d[c("q25_days", "median_days", "q75_days")])
  
  # The lower value 0 is an exact support endpoint.
  # The upper value is an assumed support bound and is not used as a finite-difference boundary node.
  set.seed(gaussian_seed)
  region <- qonlyreg::qlm_region(
    qhat = qhat,
    tau = tau_grid,
    X = X,
    n_g = d$assessments_completed,
    y_lower = 0,
    y_upper_bound = y_upper_bound,
    level = 0.95,
    B = simulation_draws
  )
  
  confidence_projection <- project_linear_region(
    const_mat = region$const_mat,
    const_dir = region$const_dir,
    const_rhs = region$const_rhs
  )
  
  treatment <- X[, "labour_majority"]
  
  summary_row <- data.frame(
    bandwidth_pp = bandwidth_pp,
    bandwidth_role = bandwidth_role(bandwidth_pp),
    treatment_definition = treatment_name,
    model_specification = model_specification,
    reporting_year_fixed_effects = reporting_year_fixed_effects,
    observations = nrow(d),
    unique_councils = length(unique(d$new_la_code)),
    treated = sum(treatment == 1),
    control = sum(treatment == 0),
    exactly_50 = sum(d$labour_exactly_50),
    support_lower_days = 0,
    support_upper_days = y_upper_bound,
    support_status = "lower_exact_endpoint_upper_assumed_bound",
    simultaneous_region_level = region$level,
    simultaneous_region_status =
      if (confidence_projection$status == 0) {
        "nonempty"
      } else {
        "empty"
      },
    
    intercept_projection_lower_95 = confidence_projection$bounds["(Intercept)", "lower"],
    intercept_projection_upper_95 = confidence_projection$bounds["(Intercept)", "upper"],
    treatment_projection_lower_95 = confidence_projection$bounds["labour_majority", "lower"],
    treatment_projection_upper_95 = confidence_projection$bounds["labour_majority", "upper"],
    simultaneous_critical_value = region$critical_value,
    simulation_draws = simulation_draws,
    exact_mean_model =
      if (reporting_year_fixed_effects) {
        paste0("mean_ct = intercept + treatment_effect * LabourMajority_ct + reporting_year_fixed_effects")
      } else {
        paste0("mean_ct = intercept + treatment_effect * LabourMajority_ct")
      },
    
    group_sample_sizes_used_in_covariance = TRUE,
    qonlyreg_version = as.character(utils::packageVersion("qonlyreg")),
    stringsAsFactors = FALSE
  )
  
  constraints <- data.frame(
    bandwidth_pp = bandwidth_pp,
    bandwidth_role = bandwidth_role(bandwidth_pp),
    treatment_definition = treatment_name,
    model_specification = model_specification,
    reporting_year_fixed_effects =
      reporting_year_fixed_effects,
    
    support_lower_days = 0,
    support_upper_days = y_upper_bound,
    support_lower_type = "exact_endpoint",
    support_upper_type = "assumed_bound",
    
    simultaneous_region_level = region$level,
    simultaneous_critical_value = region$critical_value,
    constraint_group = seq_len(nrow(d)),
    new_la_code = d$new_la_code,
    la_name = d$la_name,
    reporting_year = d$reporting_year,
    council_composition_year =
      d$council_composition_year,
    running_variable_pp = d$running_variable_pp,
    treatment = treatment,
    intercept_coefficient = X[, "(Intercept)"],
    treatment_coefficient = X[, "labour_majority"],
    estimated_mean_lower = region$bounds$m_lower,
    estimated_mean_upper = region$bounds$m_upper,
    lower_adjusted_95 = region$lower_adjusted,
    upper_adjusted_95 = region$upper_adjusted,
    stringsAsFactors = FALSE
  )
  
  all_year_columns <- paste0(
    "reporting_year_",
    sort(unique(data$reporting_year))[-1]
  )
  
  for (nm in all_year_columns) {
    constraints[[paste0(nm, "_coefficient")]] <-
      if (nm %in% colnames(X)) {
        X[, nm]
      } else {
        0
      }
  }
  
  list(
    summary = summary_row,
    constraints = constraints
  )
}

region_results <- list()
result_counter <- 0L
design_counter <- 0L
for (bw in bandwidths_pp) {
  for (j in seq_len(nrow(analysis_specifications))) {
    design_counter <- design_counter + 1L
    spec <- analysis_specifications[j, ]
    # The same Gaussian draws are used for all upper bounds
    # within a given bandwidth/model specification.
    gaussian_seed <- 20260806L + design_counter
    for (upper in upper_bounds) {
      result_counter <- result_counter + 1L
      region_results[[result_counter]] <-
        simultaneous_region_one(
          bandwidth_pp = bw,
          treatment_name = spec$treatment_definition,
          treatment_column = spec$treatment_column,
          model_specification = spec$model_specification,
          reporting_year_fixed_effects =
            spec$reporting_year_fixed_effects,
          y_upper_bound = upper,
          gaussian_seed = gaussian_seed
        )
    }
  }
}

region_summary <- do.call(rbind, lapply(region_results, `[[`, "summary"))
region_constraints <- do.call(rbind, lapply(region_results, `[[`, "constraints"))
write.csv(
  region_summary,
  file.path(output_dir, "qonlyreg_labour_fixed_window_simultaneous_region_summary.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  region_constraints,
  file.path(output_dir, "qonlyreg_labour_fixed_window_simultaneous_region_constraints.csv"),
  row.names = FALSE,
  na = ""
)

main_point <- point_results[
  point_results$bandwidth_pp == 4 &
    point_results$treatment_definition == "at_least_50_main" &
    point_results$model_specification == "no_year_fixed_effects",
  ,
  drop = FALSE
]
message(
  "Main point estimate with known lower support (+/-4pp, Labour >=50%): ",
  sprintf("%.3f days, N=%d", main_point$estimate_days, main_point$observations)
)
message("Point estimates and simultaneous-confidence-region outputs written to output/.")
