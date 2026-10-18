########################
# Quantile-Only Regression: Inference and Partial Identification under Quantile Disclosure
# Empirical Illustration: Healthcare-Associated Infections and Medical Expenditure
#########################
# Outcome:
#   Total medical expenditure per inpatient, HAI minus non-HAI, in CNY.
#
# Model:
#   m_{sa} = alpha_s + beta * HAI_{sa},
#
# where s indexes the seven study(-year) comparisons and a indexes the
# HAI/control arm. After profiling out the study fixed effects, beta is a
# sample-size-weighted average of the seven reconstructed mean differences.
#
# Li et al. (2013-2015): use the original paired-difference SE reported in Li et al. (2017).
#
# Other studies:
#   when paired covariance information is unavailable, use the working
#   assumption Cov(HAI summary estimator, control summary estimator) = 0.
#
# Inputs:
#   study_data.csv
#
# Required package:
#   qonlyreg, including qlm() and qlm_contrast_region()
#
# Outputs:
#   table_liu_main.csv
#   table_liu_support_sensitivity.csv
#


library("qonlyreg")


# -------------------------------------------------------------------
# 1. settings
# -------------------------------------------------------------------

d <- read.csv("study_data.csv", na.strings = c("", "NA"), stringsAsFactors = FALSE)

tau <- c(0.25, 0.50, 0.75)
level <- 0.95 # confidence level for both Wald CI and simultaneous CR
B <- 100000L # number of simulation draws for the calculation of simultaneous critical value in qlm_contrast_region()
seed <- 20260929L

# Main assumed upper bound on the outcome support.
# Alternative upper bounds are examined in the sensitivity analysis.
U_main <- 200000

# For the sensitivity analysis of the upper support bound
U_grid <- c(70000, 80000, 90000, 100000, 200000, 500000)


# -------------------------------------------------------------------
# 2. Construct qonlyreg disclosure data
# -------------------------------------------------------------------

make_arm_rows <- function(i, arm, group_id) {

    n <- d[[paste0("n_", arm)]][i]

    if (d$format[i] == "quartiles") {
        q <- as.numeric( d[i, paste0(c("q1_", "median_", "q3_"), arm)])
        stopifnot(
            all(is.finite(q)),
            all(diff(q) > 0)
        )

        data.frame(
            group = group_id,
            study_id = d$study_id[i],
            arm = arm,
            type = "quantile",
            tau = tau,
            value = q,
            n = rep(n, length(tau)),
            sd = NA_real_,
            stringsAsFactors = FALSE
        )

    } else {
        m <- d[[paste0("mean_", arm)]][i]
        s <- d[[paste0("sd_", arm)]][i]

        stopifnot(
            is.finite(m),
            is.finite(s),
            s >= 0
        )

        data.frame(
            group = group_id,
            study_id = d$study_id[i],
            arm = arm,
            type = "mean",
            tau = NA_real_,
            value = m,
            n = n,
            sd = s,
            stringsAsFactors = FALSE
        )
    }
}


# qlm() requires each arm to be a separate group.
point_disclosure <- do.call(
    rbind,
    lapply(seq_len(nrow(d)), function(i) {
        rbind(
            make_arm_rows(
                i = i,
                arm = "control",
                group_id = paste0(d$study_id[i], "__control")
            ),
            make_arm_rows(
                i = i,
                arm = "hai",
                group_id = paste0(d$study_id[i], "__hai")
            )
        )
    })
)

# qlm_contrast_region() matches HAI and control arms by study_id.
contrast_disclosure <- point_disclosure
contrast_disclosure$group <- contrast_disclosure$study_id


# -------------------------------------------------------------------
# 3. Point estimator: study fixed effects + common HAI coefficient
# -------------------------------------------------------------------

arm_map <- unique(
  point_disclosure[c("group", "study_id", "arm")]
)

# Seven study fixed effects, no overall intercept, plus HAI indicator.
X_study <- sapply(
    d$study_id,
    function(id) as.numeric(arm_map$study_id == id)
)

X_study <- as.matrix(X_study)
colnames(X_study) <- paste0("study_", d$study_id)

X_point <- cbind(
    X_study,
    HAI = as.numeric(arm_map$arm == "hai")
)

fit_point <- qlm(
    disclosure = point_disclosure,
    X = X_point
)


### For inference: for Li et al. (2013-2015), use the original paired-difference SE reported in Li et al. (2017).
arm_ids <- as.character(arm_map$group)
v_arm <- fit_point$var_m_hat
names(v_arm) <- arm_ids

study_var <- numeric(nrow(d))
w_raw <- numeric(nrow(d))

for (i in seq_len(nrow(d))) {
    id_h <- paste0(d$study_id[i], "__hai")
    id_c <- paste0(d$study_id[i], "__control")

    # Li: use the paired-difference SE from the primary paper.
    if (!is.na(d$paired_se_usd[i])) {
        study_var[i] <- (d$paired_se_usd[i] * d$usd_to_cny[i])^2
    } else {
        # Working zero cross-arm covariance assumption.
        study_var[i] <- v_arm[id_h] + v_arm[id_c]
    }

    # Profiling the study fixed effect from
    # n_H (m_H-alpha-beta)^2 + n_C (m_C-alpha)^2
    # gives a study-difference weight proportional to:
    w_raw[i] <-
        d$n_hai[i] * d$n_control[i] /
        (d$n_hai[i] + d$n_control[i])
}

w <- w_raw / sum(w_raw)

beta_hat <- as.numeric(coef(fit_point)[ncol(X_point), 1])
beta_se <- sqrt(sum(w^2 * study_var))

zcrit <- qnorm(1 - (1 - level) / 2)

beta_ci <- c(
    beta_hat - zcrit * beta_se,
    beta_hat + zcrit * beta_se
)

# qlm() estimates exactly the same study-FE coefficient from the arm means.
beta_qlm <- as.numeric(coef(fit_point)[ncol(X_point), 1])

stopifnot(
    abs(beta_hat - beta_qlm) < 1e-6
)


# -------------------------------------------------------------------
# 4. Fixed-K contrast region
# -------------------------------------------------------------------
psd_root <- function(S) {
    S <- (S + t(S)) / 2
    ee <- eigen(S, symmetric = TRUE)
    tol <- 1e-8 * max(1, max(abs(ee$values)))
    if (min(ee$values) < -tol) {
        stop("Covariance matrix is not positive semidefinite.")
    }
    ee$vectors %*%
        diag(sqrt(pmax(ee$values, 0)), nrow(S)) %*%
        t(ee$vectors)
}


make_region <- function(U) {

    stopifnot(
        length(U) == 1,
        is.finite(U),
        U > 0
    )

    li_ids <- d$study_id[!is.na(d$paired_se_usd)]
    nonli_ids <- d$study_id[is.na(d$paired_se_usd)]

    # qonlyreg handles all non-Li studies under the zero cross-arm
    # covariance working assumption built into qlm_contrast_region().
    nonli_disclosure <- contrast_disclosure[
        contrast_disclosure$study_id %in% nonli_ids,
        ,
        drop = FALSE
    ]

    X_nonli <- matrix(
        1,
        nrow = length(nonli_ids),
        ncol = 1
    )
    colnames(X_nonli) <- "beta0"

    # We need qonlyreg's estimated contrast bounds and covariance matrix.
    # Its simulated critical value is not used, because Li's covariance
    # blocks are replaced below by the reported paired-difference variance.
    set.seed(seed)

    pkg_region <- qlm_contrast_region(
      disclosure = nonli_disclosure,
      X = X_nonli,
      treatment = "hai",
      control = "control",
      y_lower_bound = 0,
      y_upper_bound = U,
      level = level,
      B = 100
    )

    # qlm_contrast_region stores Omega for sqrt(N)*(estimated bounds - bounds).
    # Divide by the N used in that package call to recover the actual
    # sampling covariance of the estimated bound endpoints.
    N_nonli <- sum(pkg_region$n_g)
    V_nonli <- pkg_region$Omega / N_nonli

    G <- nrow(d)
    theta_lower <- numeric(G)
    theta_upper <- numeric(G)

    # Actual sampling covariance of
    # (lower_1, upper_1, ..., lower_G, upper_G).
    V <- matrix(0, nrow = 2 * G, ncol = 2 * G)

    for (i in seq_len(G)) {

        id <- d$study_id[i]
        idx <- c(2 * i - 1, 2 * i)

        if (id %in% li_ids) {
            # Exact mean disclosure: singleton mean-difference constraint.
            # Use Liu Table IV arm means for the point constraint.
            theta_i <- d$mean_hai[i] - d$mean_control[i]

            theta_lower[i] <- theta_i
            theta_upper[i] <- theta_i

            # Primary-paper paired-difference variance.
            var_i <- (d$paired_se_usd[i] * d$usd_to_cny[i])^2
            V[idx, idx] <- var_i * matrix(1, 2, 2)

        } else {
            j <- match(id, pkg_region$id)
            
            theta_lower[i] <- pkg_region$bounds$theta_lower[j]
            theta_upper[i] <- pkg_region$bounds$theta_upper[j]

            jdx <- c(2 * j - 1, 2 * j)

            V[idx, idx] <- V_nonli[jdx, jdx, drop = FALSE]
        }
    }

    stopifnot(
        all(theta_lower <= theta_upper)
    )

    s <- sqrt(diag(V))

    if (any(!is.finite(s)) || any(s <= 0)) {
        stop("Non-finite or non-positive standard error in contrast region.")
    }

    # Simultaneous Gaussian critical value using the corrected covariance, including Li's paired-difference blocks.
    R <- V / outer(s, s)
    root <- psd_root(R)

    set.seed(seed)

    Z <- matrix(
        rnorm(B * nrow(R)),
        nrow = B,
        ncol = nrow(R)
    ) %*% t(root)

    Tsim <- apply(abs(Z), 1, max)

    critical_value <- as.numeric(
        quantile(
            Tsim,
            probs = level,
            type = 8
        )
    )

    lower_adj <- theta_lower - critical_value * s[seq(1, 2 * G, by = 2)]
    upper_adj <- theta_upper + critical_value * s[seq(2, 2 * G, by = 2)]

    plugin_lower <- max(theta_lower)
    plugin_upper <- min(theta_upper)

    region_lower <- max(lower_adj)
    region_upper <- min(upper_adj)

    list(
        critical_value = critical_value,
        theta_lower = theta_lower,
        theta_upper = theta_upper,
        plugin_lower = plugin_lower,
        plugin_upper = plugin_upper,
        plugin_empty = plugin_lower > plugin_upper,
        region_lower = region_lower,
        region_upper = region_upper,
        region_empty = region_lower > region_upper
    )
}



# -------------------------------------------------------------------
# 5. Manuscript tables
# -------------------------------------------------------------------

fmt_interval <- function(lower, upper, empty = FALSE) {

    if (empty) {
        return("Empty set")
    }

    sprintf("[%.2f, %.2f]", lower, upper)
}


main_region <- make_region(U_main)

table_main <- data.frame(
    Method = c(
        "Point estimator",
        "Estimated identified set",
        "Feasible confidence region"
    ),
    Estimate = c(
        round(beta_hat, 2),
        NA,
        NA
    ),
    Interval = c(
        fmt_interval(beta_ci[1], beta_ci[2]),
        fmt_interval(
            main_region$plugin_lower,
            main_region$plugin_upper,
            main_region$plugin_empty
        ),
        fmt_interval(
            main_region$region_lower,
            main_region$region_upper,
            main_region$region_empty
        )
    ),
    stringsAsFactors = FALSE
)

write.csv(
    table_main,
    "table_liu_main.csv",
    row.names = FALSE,
    na = ""
)


support_list <- lapply(
    U_grid,
    function(U) {

        r <- make_region(U)

        data.frame(
            Upper_support_CNY = as.integer(U),
            Plugin_identified_region = fmt_interval(
                r$plugin_lower,
                r$plugin_upper,
                r$plugin_empty
            ),
            Feasible_95pct_region = fmt_interval(
                r$region_lower,
                r$region_upper,
                r$region_empty
            ),
            stringsAsFactors = FALSE
        )
    }
)

table_support <- do.call(rbind, support_list)

write.csv(
    table_support,
    "table_liu_support_sensitivity.csv",
    row.names = FALSE
)


# -------------------------------------------------------------------
# 6. Console output
# -------------------------------------------------------------------

cat("\nEmpirical Illustration: Healthcare-Associated Infections and Medical Expenditure\n")
cat("------------------------------------------------\n")
cat(sprintf(
    "Point estimate: %.2f CNY\n",
    beta_hat
))
cat(sprintf(
    "Working SE: %.2f\n",
    beta_se
))
cat(sprintf(
    "Working 95%% confidence interval: [%.2f, %.2f]\n",
    beta_ci[1],
    beta_ci[2]
))
cat(sprintf(
    "Main support assumption: U = %.0f CNY\n",
    U_main
))
cat(
    "Estimated identified set: ",
    table_main$Interval[2],
    "\n",
    sep = ""
)
cat(
    "Feasible 95% confidence region: ",
    table_main$Interval[3],
    "\n\n",
    sep = ""
)
cat("Saved manuscript tables:\n")
cat("  table_liu_main.csv\n")
cat("  table_liu_support_sensitivity.csv\n")

