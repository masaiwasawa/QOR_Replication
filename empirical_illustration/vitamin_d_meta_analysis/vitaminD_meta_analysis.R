############################
## Quantile-Only Regression: Inference and Partial Identification under Quantile Disclosure
## Vitamin D meta-analysis using qonlyreg
############################

library(qonlyreg)
library(dplyr)
library(readr)
library(knitr)
library(kableExtra)

df <- read_csv("vitaminD_meta_analysis_table1.csv", show_col_types = FALSE)

U <- 150
alpha_level <- 0.05
z <- qnorm(1 - alpha_level / 2)
B <- 10000

## SD approximation from range for mean_range study
sd_from_range <- function(min_y, max_y, n) {
  (max_y - min_y) / (2 * qnorm((n - 0.375) / (n + 0.25)))
}

# quantile-based mean
qmean_hat <- function(q0, q50, q100) {
  0.25 * q0 + 0.50 * q50 + 0.25 * q100
}

## Var(alpha' Qhat) = alpha' Sigma_Q alpha
## Sigma_Q[i,j] = {min(tau_i,tau_j)-tau_i*tau_j} Q'(tau_i) Q'(tau_j) / n
qmean_var_hat <- function(q0, q50, q100, n) {
  tau <- c(0, 0.5, 1)
  alpha <- c(0.25, 0.50, 0.25)
  q <- c(q0, q50, q100)

  ## finite-difference estimate of Q'(tau)
  qprime <- c(
    (q[2] - q[1]) / (tau[2] - tau[1]),
    (q[3] - q[1]) / (tau[3] - tau[1]),
    (q[3] - q[2]) / (tau[3] - tau[2])
  )

  Sigma <- outer(tau, tau, function(a, b) pmin(a, b) - a * b)
  Sigma <- Sigma * outer(qprime, qprime) / n

  as.numeric(t(alpha) %*% Sigma %*% alpha)
}

dat <- df %>%
  filter(stat_type != "odds_ratio") %>%
  rowwise() %>%
  mutate(
    ## -----------------------------
    ## Point-ID mean estimates
    ## -----------------------------
    case_mhat = case_when(
      stat_type == "median_range" ~ qmean_hat(case_min, case_median, case_max),
      stat_type %in% c("mean_sd", "mean_range") ~ case_mean
    ),
    control_mhat = case_when(
      stat_type == "median_range" ~ qmean_hat(control_min, control_median, control_max),
      stat_type %in% c("mean_sd", "mean_range") ~ control_mean
    ),
    theta_hat = case_mhat - control_mhat,

    ## -----------------------------
    ## Point-ID variance estimates
    ## -----------------------------
    case_var_mhat = case_when(
      stat_type == "median_range" ~ qmean_var_hat(case_min, case_median, case_max, case_n),
      stat_type == "mean_sd" ~ case_sd^2 / case_n,
      stat_type == "mean_range" ~ sd_from_range(case_min, case_max, case_n)^2 / case_n
    ),
    control_var_mhat = case_when(
      stat_type == "median_range" ~ qmean_var_hat(control_min, control_median, control_max, control_n),
      stat_type == "mean_sd" ~ control_sd^2 / control_n,
      stat_type == "mean_range" ~ sd_from_range(control_min, control_max, control_n)^2 / control_n
    ),

    theta_var = case_var_mhat + control_var_mhat,
    theta_se = sqrt(theta_var),
    n_g = case_n + control_n
  ) %>%
  ungroup()

############################
## Point estimation using qlm()
############################

theta_disclosure <- dat %>%
  transmute(
    study = study_id,
    type = "mean",
    value = theta_hat,
    se = theta_se,
    n = n_g
  )

X_theta <- matrix(1, nrow = nrow(theta_disclosure), ncol = 1)
colnames(X_theta) <- "intercept"

fit <- qlm(
  disclosure = theta_disclosure,
  X = X_theta
)

beta_hat <- unname(coef(fit)[1])
beta_se <- unname(sqrt(vcov(fit)[1, 1]))
beta_ci <- beta_hat + c(-1, 1) * z * beta_se

############################
## Arm-level disclosure for partial identification
############################

case_arm <- dat %>%
  transmute(
    study = study_id,
    arm = "case",
    stat_type = stat_type,
    n = case_n,
    median = case_median,
    mean = case_mean,
    sd = case_when(
      stat_type == "mean_sd" ~ case_sd,
      stat_type == "mean_range" ~ sd_from_range(case_min, case_max, case_n),
      TRUE ~ NA_real_
    )
  )

control_arm <- dat %>%
  transmute(
    study = study_id,
    arm = "control",
    stat_type = stat_type,
    n = control_n,
    median = control_median,
    mean = control_mean,
    sd = case_when(
      stat_type == "mean_sd" ~ control_sd,
      stat_type == "mean_range" ~ sd_from_range(control_min, control_max, control_n),
      TRUE ~ NA_real_
    )
  )

arm_df <- bind_rows(case_arm, control_arm)

median_rows <- arm_df %>%
  filter(stat_type == "median_range") %>%
  transmute(
    study = study,
    arm = arm,
    type = "quantile",
    tau = 0.5,
    value = median,
    n = n
  )

mean_rows <- arm_df %>%
  filter(stat_type %in% c("mean_sd", "mean_range")) %>%
  transmute(
    study = study,
    arm = arm,
    type = "mean",
    tau = NA_real_,
    value = mean,
    n = n,
    sd = sd
  )

arm_disclosure <- bind_rows(median_rows, mean_rows)

############################
## Plug-in identified region using qlm_contrast_bounds()
############################

contrast_bounds <- qlm_contrast_bounds(
  disclosure = arm_disclosure,
  treatment = "case",
  control = "control",
  y_lower = 0,
  y_upper_bound = U
)

id_region <- c(
  lower = max(contrast_bounds$theta_lower),
  upper = min(contrast_bounds$theta_upper)
)

############################
## Feasible confidence region using qlm_contrast_region()
############################

X_contrast <- matrix(1, nrow = nrow(dat), ncol = 1)
colnames(X_contrast) <- "intercept"

set.seed(123)

contrast_region <- qlm_contrast_region(
  disclosure = arm_disclosure,
  X = X_contrast,
  treatment = "case",
  control = "control",
  y_lower = 0,
  y_upper_bound = U,
  level = 1 - alpha_level,
  B = B
)

ci_region <- confint(contrast_region)
simul_conf_region <- ci_region$confint["intercept", ]

############################
## Output
############################

cat("\nPoint-ID beta_hat =", beta_hat, "\n")
cat("Point-ID SE       =", beta_se, "\n")
cat("Point-ID CI       = [", beta_ci[1], ", ", beta_ci[2], "]\n")

if (id_region[1] <= id_region[2]) {
  cat("ID region = [", id_region[1], ", ", id_region[2], "]\n")
} else {
  cat("ID region is empty: lower =", id_region[1],
      ", upper =", id_region[2], "\n")
}

cat("Simultaneous critical value q_hat =",
    contrast_region$critical_value, "\n")
cat("Simultaneous confidence region = [",
    simul_conf_region[1], ", ", simul_conf_region[2], "]\n")

############################
## Table
############################

results_tbl <- tibble(
  Method = c(
    "Point estimator",
    "Plug-in identified region",
    "Feasible confidence region"
  ),
  Estimate = c(
    sprintf("%.2f", beta_hat),
    "",
    ""
  ),
  Interval = c(
    sprintf("[%.2f, %.2f]", beta_ci[1], beta_ci[2]),
    if (id_region[1] <= id_region[2]) {
      sprintf("[%.2f, %.2f]", id_region[1], id_region[2])
    } else {
      "Empty set"
    },
    sprintf("[%.2f, %.2f]", simul_conf_region[1], simul_conf_region[2])
  )
)

tex_table <-
  kable(
    results_tbl,
    format = "latex",
    booktabs = TRUE,
    align = c("l", "c", "c"),
    caption = "Vitamin D levels and tuberculosis: an application to meta-analysis",
    label = "vitd_main",
    escape = FALSE
  ) %>%
  kable_styling(
    latex_options = c("hold_position")
  )

writeLines(
  as.character(tex_table),
  "table_vitaminD_meta.tex"
)
