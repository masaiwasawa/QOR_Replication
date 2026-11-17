########################
# Quantile-Only Regression: Inference and Partial Identification under Quantile Disclosure
# Supplement D4: Partial identification inference under truncated-normal errors
############################

rm(list = ls())
library(qonlyreg)
library(dplyr)
library(tidyr)
library(ggplot2)
library(doParallel)
library(foreach)
library(doRNG)

set.seed(123)

############################
## Helper functions
############################

make_grid <- function(K, type = "even") {
  p <- seq(1, K) / (K + 1)

  if (type == "even") {
    tau <- p
  } else if (type == "right_tail_coarse_mild") {
    tau <- p^1.5
  } else if (type == "right_tail_coarse_severe") {
    tau <- p^3
  } else {
    stop("Unknown grid type")
  }

  sort(tau)
}

trunc_lower <- -2
trunc_upper <- 3
generate_data <- function(N, G = 10, beta0 = 1, beta1 = 2,
                          a = trunc_lower,
                          b = trunc_upper) {

  X_grid <- seq(-1, 1, length.out = G)
  g <- rep(1:G, length.out = N)
  Xg <- X_grid[g]

  prob_a <- pnorm(a)
  prob_b <- pnorm(b)
  normalizing_constant <- prob_b - prob_a
  trunc_mean <- (dnorm(a) - dnorm(b)) / normalizing_constant
  trunc_var <- 1 +
    (a * dnorm(a) - b * dnorm(b)) / normalizing_constant -
    trunc_mean^2
  trunc_sd <- sqrt(trunc_var)

  z_trunc <- qnorm(prob_a + runif(N) * normalizing_constant)
  sigma_g <- exp(0.5 * Xg)
  u <- sigma_g * (z_trunc - trunc_mean) / trunc_sd
  Y <- beta0 + beta1 * Xg + u

  y_lower <- beta0 + beta1 * X_grid +
    exp(0.5 * X_grid) * (a - trunc_mean) / trunc_sd
  y_upper <- beta0 + beta1 * X_grid +
    exp(0.5 * X_grid) * (b - trunc_mean) / trunc_sd

  data.frame(
    Y = Y,
    g = g,
    Xg = Xg,
    sigma_g = sigma_g,
    y_lower_g = y_lower[g],
    y_upper_g = y_upper[g]
  )
}

make_disclosed_quantiles <- function(df, tau) {

  group_tbl <- df %>%
    group_by(g, Xg) %>%
    summarise(
      qhat = list(as.numeric( stats::quantile(Y, probs = tau, type = 1))),
      n_g = dplyr::n(),
      y_lower = first(y_lower_g),
      y_upper = first(y_upper_g),
      .groups = "drop"
    )

  list(
    qhat = do.call(rbind, group_tbl$qhat),
    X = stats::model.matrix(~ Xg, data = group_tbl),
    n_g = group_tbl$n_g,
    y_lower = group_tbl$y_lower,
    y_upper = group_tbl$y_upper,
    group_tbl = group_tbl
  )
}




############################
## Helpers for true identified set coverage
############################
# true identified interval for the mean
true_bounds_truncnorm <- function(X_grid, tau, beta0 = 1, beta1 = 2,
                                 a = -2, b = 3) {
  sigma_g <- exp(0.5 * X_grid)

  prob_a <- pnorm(a)
  prob_b <- pnorm(b)
  normalizing_constant <- prob_b - prob_a
  trunc_mean <- (dnorm(a) - dnorm(b)) / normalizing_constant
  trunc_var <- 1 + (a * dnorm(a) - b * dnorm(b)) / normalizing_constant - trunc_mean^2
  trunc_sd <- sqrt(trunc_var)

  # True quantiles under the standardized truncated-normal errors.
  Q_true <- sapply(tau, function(t) {
    beta0 + beta1 * X_grid + sigma_g *
      (qnorm(prob_a + t * normalizing_constant) - trunc_mean) / trunc_sd
  })

  tau_full <- c(0, tau, 1)
  delta_tau <- diff(tau_full)

  # Exact support endpoints for the standardized truncated-normal errors.
  y_lower <- beta0 + beta1 * X_grid + sigma_g * (a - trunc_mean) / trunc_sd
  y_upper <- beta0 + beta1 * X_grid + sigma_g * (b - trunc_mean) / trunc_sd

  weights_lower <- delta_tau[-1]
  weights_upper <- delta_tau[-length(delta_tau)]

  m_lower <- delta_tau[1] * y_lower + as.numeric(Q_true %*% weights_lower)
  m_upper <- as.numeric(Q_true %*% weights_upper) + delta_tau[length(delta_tau)] * y_upper

  list(
    m_lower = m_lower,
    m_upper = m_upper,
    y_lower = y_lower,
    y_upper = y_upper
  )
}

# vertices of identification region of two-dimensional parameters
get_vertices_2d <- function(X, lower, upper, tol = 1e-8) {
  A <- rbind(X, -X)
  b <- c(upper, -lower)

  m <- nrow(A)
  vertices <- list()

  for (i in seq_len(m - 1)) {
    for (j in (i + 1):m) {
      Aij <- rbind(A[i, ], A[j, ])
      bij <- c(b[i], b[j])

      if (abs(det(Aij)) < tol) next

      v <- as.numeric(solve(Aij, bij))

      if (all(A %*% v <= b + tol)) {
        vertices[[length(vertices) + 1]] <- v
      }
    }
  }

  if (length(vertices) == 0) {
    return(matrix(numeric(0), ncol = 2))
  }

  V <- unique(round(do.call(rbind, vertices), 10))

  center <- colMeans(V)
  angle <- atan2(V[, 2] - center[2], V[, 1] - center[1])
  V[order(angle), , drop = FALSE]
}


region_contains_vertices <- function(region, vertices, tol = 1e-6) {
  if (nrow(vertices) == 0) return(FALSE)

  lhs <- region$const_mat %*% t(vertices)

  ok <- logical(ncol(lhs))

  for (v in seq_len(ncol(lhs))) {
    ok[v] <- all(
      (region$const_dir == ">=" & lhs[, v] >= region$const_rhs - tol) |
        (region$const_dir == "<=" & lhs[, v] <= region$const_rhs + tol)
    )
  }

  all(ok)
}



############################
## Simulation design
############################

N_list <- c(1000, 5000, 10000)
K_list <- c(3, 5, 7, 10, 15, 20, 50)
grid_types <- c("even", "right_tail_coarse_mild", "right_tail_coarse_severe")

G <- 10
M <- 10000
level <- 0.95
B_gauss <- 10000

design_data <- expand.grid(
  N = N_list,
  rep = 1:M
)

block_size <- 100
blocks <- split(seq_len(nrow(design_data)),
                ceiling(seq_len(nrow(design_data)) / block_size))


############################
## Run simulation in parallel
############################

start_time <- Sys.time()

n_cores <- max(1, parallel::detectCores() - 1)

cl <- makeCluster(n_cores)
registerDoParallel(cl)
registerDoRNG(123)

results <- foreach(
  b = blocks,
  .combine = "rbind",
  .packages = c("qonlyreg", "dplyr"),
  .export = c(
    "generate_data",
    "trunc_lower",
    "trunc_upper",
    "make_grid",
    "make_disclosed_quantiles",
    "true_bounds_truncnorm",
    "get_vertices_2d",
    "region_contains_vertices",
    "K_list",
    "grid_types",
    "level",
    "B_gauss",
    "G"
  )
) %dorng% {

  block_out <- lapply(b, function(i) {

    N_i <- design_data$N[i]
    rep_i <- design_data$rep[i]

    ## Generate common data once for this (N, rep)
    df_i <- generate_data(
      N = N_i,
      G = G,
      beta0 = 1,
      beta1 = 2
    )

    do.call(
      rbind,
      lapply(K_list, function(K_i) {

        do.call(
          rbind,
          lapply(grid_types, function(grid_i) {

            tau <- make_grid(K_i, grid_i)

            summaries <- make_disclosed_quantiles(
              df = df_i,
              tau = tau
            )

            region <- qlm_region(
              qhat = summaries$qhat,
              tau = tau,
              X = summaries$X,
              n_g = summaries$n_g,
              y_lower = summaries$y_lower,
              y_upper = summaries$y_upper,
              level = level,
              B = B_gauss
            )

            X_grid <- seq(-1, 1, length.out = G)

            true_bounds <- true_bounds_truncnorm(
              X_grid = X_grid,
              tau = tau,
              beta0 = 1,
              beta1 = 2
            )

            X_true <- stats::model.matrix(~ X_grid)

            vertices_true <- get_vertices_2d(
              X = X_true,
              lower = true_bounds$m_lower,
              upper = true_bounds$m_upper
            )

            cover_region <- as.numeric(
              region_contains_vertices(region, vertices_true)
            )

            ci <- confint(region)$confint

            beta0_lower <- ci[1, "lower"]
            beta0_upper <- ci[1, "upper"]

            beta1_lower <- ci[2, "lower"]
            beta1_upper <- ci[2, "upper"]

            out <- c(
              beta0_lower = beta0_lower,
              beta0_upper = beta0_upper,
              beta1_lower = beta1_lower,
              beta1_upper = beta1_upper,

              cover_region = cover_region,
              n_vertices_true = nrow(vertices_true),

              cover_beta0 = as.numeric(
                beta0_lower <= 1 & 1 <= beta0_upper
              ),

              cover_beta1 = as.numeric(
                beta1_lower <= 2 & 2 <= beta1_upper
              ),

              length_beta0 = beta0_upper - beta0_lower,
              length_beta1 = beta1_upper - beta1_lower,

              critical_value = region$critical_value,

              Delta = max(
                c(
                  tau[1],
                  diff(tau),
                  1 - tau[length(tau)]
                )
              )
            )

            data.frame(
              N = N_i,
              K = K_i,
              grid_type = grid_i,
              rep = rep_i,
              t(out)
            )
          })
        )
      })
    )
  })

  do.call(rbind, block_out)
}

stopCluster(cl)
registerDoSEQ()

end_time <- Sys.time()
elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))


dir.create("results_compact_truncnorm_qlm_partial", showWarnings = FALSE)
saveRDS(results, "results_compact_truncnorm_qlm_partial/results.rds")


############################
## Summaries
# results <- readRDS("results_compact_truncnorm_qlm_partial/results.rds")
############################
summary_tbl <- results %>%
  group_by(N, K, grid_type) %>%
  summarise(
    Delta = mean(Delta),
    coverage_region = mean(cover_region),
    coverage_beta0 = mean(cover_beta0),
    coverage_beta1 = mean(cover_beta1),
    mean_length_beta0 = mean(length_beta0),
    mean_length_beta1 = mean(length_beta1),
    mean_critical_value = mean(critical_value),
    mean_n_vertices_true = mean(n_vertices_true),
    .groups = "drop"
  )
print(summary_tbl)

write.csv(
  summary_tbl,
  "results_compact_truncnorm_qlm_partial/summary_table.csv",
  row.names = FALSE
)


saveRDS(summary_tbl, "results_compact_truncnorm_qlm_partial/summary_tbl.rds")

############################
## Table for parameter wise coverage
############################

library(knitr)
library(kableExtra)

tex_table <- summary_tbl %>%
  mutate(across(where(is.numeric), ~ round(., 4))) %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    caption = "Partial identification inference under truncated-normal errors",
    label = "tab:partial_identification_simulation_truncnorm",
    linesep = ""
  ) %>%
  kable_styling(latex_options = c("hold_position", "scale_down"))

cat(
  tex_table,
  file = "results_compact_truncnorm_qlm_partial/table_parameterCR.tex"
)


############################
## Table for simultaneous confidence region
############################

for (grid_i in unique(summary_tbl$grid_type)) {
  grid_label <- dplyr::recode(
    grid_i,
    "even" = "even",
    "right_tail_coarse_mild" = "mildly right-tail-coarse",
    "right_tail_coarse_severe" = "severely right-tail-coarse"
  )

  coverage_region_tbl_grid <- summary_tbl %>%
    filter(grid_type == grid_i) %>%
    select(
      N, K,
      coverage_region,
      mean_length_beta0,
      mean_length_beta1
    ) %>%
    mutate(
      K = paste0("$K = ", K, "$"),
      coverage_region = round(coverage_region, 4),
      mean_length_beta0 = round(mean_length_beta0, 4),
      mean_length_beta1 = round(mean_length_beta1, 4)
    ) %>%
    pivot_wider(
      names_from = K,
      values_from = c(
        coverage_region,
        mean_length_beta0,
        mean_length_beta1
      ),
      names_glue = "{.value}_{K}"
    ) %>%
    arrange(N)

  coverage_region_tex_grid <- coverage_region_tbl_grid %>%
    kable(
      format = "latex",
      booktabs = TRUE,
      align = "c",
      escape = FALSE,
      row.names = FALSE,
      col.names = c(
        "$N$",
        rep(c("$K = 3$", "$K = 5$", "$K = 7$","$K = 10$", "$K = 15$", "$K = 20$", "$K = 50$"), 3)
      ),
      caption = paste0(
        "\\label{tab:partial_identification_region_coverage_truncnorm_",
        grid_i,
        "} ",
        "Coverage of the simultaneous confidence region under truncated-normal errors and ",
        grid_label,
        " grids"
      ),
      linesep = ""
    ) %>%
    add_header_above(
      c(
        " " = 1,
        "Coverage" = 3,
        "Length of $\\\\beta_0$ interval" = 3,
        "Length of $\\\\beta_1$ interval" = 3
      ),
      escape = FALSE
    ) %>%
    kable_styling(
      latex_options = "hold_position",
      full_width = FALSE
    )
  coverage_region_tex_grid <- gsub("\\\\addlinespace\n", "", coverage_region_tex_grid)

  coverage_region_tex_grid <- sub(
    "\\\\begin\\{table\\}",
    "\\\\begin{table}\n\\\\centering\n\\\\begin{threeparttable}",
    coverage_region_tex_grid
  )

  coverage_region_tex_grid <- sub(
    "\\\\end\\{tabular\\}",
    paste0(
      "\\\\end{tabular}\n",
      "\\\\begin{tablenotes}[flushleft]\n",
      "\\\\footnotesize\n",
      "\\\\item Notes: ",
      "Coverage reports the empirical coverage probability of the simultaneous ",
      "confidence region for the true identified set. ",
      "The length columns report the average lengths of the projection intervals.\n",
      "\\\\end{tablenotes}"
    ),
    coverage_region_tex_grid
  )

  coverage_region_tex_grid <- sub(
    "\\\\end\\{table\\}",
    "\\\\end{threeparttable}\n\\\\end{table}",
    coverage_region_tex_grid
  )

  cat(
    coverage_region_tex_grid,
    file = paste0(
      "results_compact_truncnorm_qlm_partial/table_region_coverage_",
      grid_i,
      ".tex"
    )
  )
}





############################
## Combined coverage table for all grid designs
############################

out_dir <- "results_compact_truncnorm_qlm_partial"
M_table <- 10000L

required_cols <- c("N", "K", "grid_type", "coverage_region")
if (!all(required_cols %in% names(summary_tbl))) {
  stop("summary_tbl must contain N, K, grid_type, and coverage_region.")
}

table_data <- as.data.frame(summary_tbl[, required_cols])

if (anyDuplicated(table_data[c("grid_type", "N", "K")])) {
  stop("Duplicate grid_type, N, K combinations were found.")
}

if (anyNA(table_data) ||
    any(!is.finite(table_data$coverage_region)) ||
    any(table_data$coverage_region < 0 |
        table_data$coverage_region > 1)) {
  stop("The table contains missing values or invalid coverage probabilities.")
}

# Panel order and labels.
grid_labels <- c(
  even = "Evenly spaced grid",
  right_tail_coarse_mild = "Mildly right-tail-coarse grid",
  right_tail_coarse_severe = "Severely right-tail-coarse grid"
)

if (!setequal(as.character(unique(table_data$grid_type)),
              names(grid_labels))) {
  stop("The results must contain the three expected grid designs.")
}

# Obtain N and K from the results.
K_values <- sort(unique(table_data$K))
N_values <- sort(unique(table_data$N))

n_columns <- length(K_values) + 1L
column_spec <- paste0("l", strrep("r", length(K_values)))

# For example, 10000 becomes 10{,}000 in LaTeX.
format_integer <- function(x) {
  format(x, big.mark = "{,}", scientific = FALSE, trim = TRUE)
}

# Table header.
tex_lines <- c(
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\begin{threeparttable}",
  "\\caption{Coverage of the simultaneous confidence region under truncated-normal errors.}",
  "\\label{tab:coverage_region_compact}",
  "\\setlength{\\tabcolsep}{5pt}",
  paste0("\\begin{tabular}{", column_spec, "}"),
  "\\toprule",
  paste0(
    paste(
      c("\\(N\\)", paste0("\\(K=", K_values, "\\)")),
      collapse = " & "
    ),
    " \\\\"
  ),
  "\\midrule"
)

# One panel for each grid design.
for (grid_i in names(grid_labels)) {
  
  if (grid_i != names(grid_labels)[1L]) {
    tex_lines <- c(tex_lines, "\\addlinespace")
  }
  
  tex_lines <- c(
    tex_lines,
    paste0(
      "\\multicolumn{", n_columns, "}{l}{",
      grid_labels[[grid_i]], "} \\\\"
    )
  )
  
  for (N_i in N_values) {
    
    d <- table_data[
      table_data$grid_type == grid_i & table_data$N == N_i,
      ,
      drop = FALSE
    ]
    
    # Match the values to the common K column order.
    idx <- match(K_values, d$K)
    
    if (anyNA(idx)) {
      stop(
        "Missing K values for grid ", grid_i,
        " and N = ", N_i, "."
      )
    }
    
    coverage_values <- formatC(
      d$coverage_region[idx],
      format = "f",
      digits = 4,
      decimal.mark = "."
    )
    
    tex_lines <- c(
      tex_lines,
      paste0(
        "\\(", format_integer(N_i), "\\) & ",
        paste(coverage_values, collapse = " & "),
        " \\\\"
      )
    )
  }
}

# Table footer and notes.
tex_lines <- c(
  tex_lines,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{tablenotes}[flushleft]",
  "\\footnotesize",
  paste0(
    "\\item Notes: Entries report the empirical probability that the ",
    "nominal 95\\% simultaneous confidence region contains the entire ",
    "true identified set, based on \\(", format_integer(M_table),
    "\\) Monte Carlo replications."
  ),
  "\\end{tablenotes}",
  "\\end{threeparttable}",
  "\\end{table}"
)

# Save the combined table.
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(
  out_dir,
  "table_region_coverage_all_grids.tex"
)

writeLines(tex_lines, con = out_file, useBytes = TRUE)

message("Saved: ", out_file)



############################
## Plots for parameter wise CI
############################

summary_tbl_plot <- summary_tbl %>%
  arrange(N) %>%
  mutate(
    n_g = N / G,
    facet_label = paste0("N = ", N, ", n_g = ", n_g),
    facet_label = factor(
      facet_label,
      levels = unique(facet_label)
    ),
    Delta2 = Delta^2
  )

grid_scales <- list(
  scale_color_manual(values = c(
    "even" = "black",
    "right_tail_coarse_mild" = "blue",
    "right_tail_coarse_severe" = "red"
  )),
  scale_linetype_manual(values = c(
    "even" = "solid",
    "right_tail_coarse_mild" = "dashed",
    "right_tail_coarse_severe" = "dotted"
  )),
  scale_shape_manual(values = c(
    "even" = 16,
    "right_tail_coarse_mild" = 17,
    "right_tail_coarse_severe" = 15
  ))
)

facet_theme <- theme(
  strip.background = element_blank(),
  strip.text = element_text(face = "plain")
)

legend_theme <- theme(
  legend.position = "bottom"
)

p_cov_beta1 <- summary_tbl_plot %>%
  arrange(N, grid_type, K) %>%
  ggplot(aes(
    x = K,
    y = coverage_beta1,
    color = grid_type,
    linetype = grid_type,
    shape = grid_type
  )) +
  geom_line(aes(group = grid_type)) +
  geom_point(size = 2) +
  geom_hline(yintercept = level, linetype = "longdash") +
  facet_wrap(~ facet_label) +
  labs(
    x = "Number of disclosed quantiles K",
    y = "Coverage probability for beta1",
    color = "Grid type",
    linetype = "Grid type",
    shape = "Grid type"
  ) +
  grid_scales +
  theme_bw() +
  facet_theme +
  legend_theme

ggsave(
  "results_compact_truncnorm_qlm_partial/plot_coverage_beta1.png",
  p_cov_beta1,
  width = 8,
  height = 6,
  dpi = 300
)

p_len_beta1 <- summary_tbl_plot %>%
  arrange(N, grid_type, K) %>%
  ggplot(aes(
    x = K,
    y = mean_length_beta1,
    color = grid_type,
    linetype = grid_type,
    shape = grid_type
  )) +
  geom_line(aes(group = grid_type)) +
  geom_point(size = 2) +
  facet_wrap(~ facet_label, nrow = 1, scales = "fixed") +
  labs(
    x = "Number of disclosed quantiles K",
    y = "Average confidence interval length for beta1",
    color = "Grid type",
    linetype = "Grid type",
    shape = "Grid type"
  ) +
  grid_scales +
  theme_bw() +
  facet_theme +
  legend_theme

ggsave(
  "results_compact_truncnorm_qlm_partial/plot_length_beta1.png",
  p_len_beta1,
  width = 8,
  height = 6,
  dpi = 300
)




############################
## Plots for simultaneous confidence regions
############################

polygon_df <- function(vertices, set_name, N, K) {
  data.frame(
    beta0 = vertices[, 1],
    beta1 = vertices[, 2],
    set = set_name,
    N = N,
    K = K
  )
}

region_vertices_2d <- function(region, tol = 1e-8) {
  A <- region$const_mat
  dir <- region$const_dir
  rhs <- region$const_rhs

  A_leq <- A
  b_leq <- rhs

  ge_idx <- dir == ">="
  A_leq[ge_idx, ] <- -A_leq[ge_idx, ]
  b_leq[ge_idx] <- -b_leq[ge_idx]

  m <- nrow(A_leq)
  vertices <- list()

  for (i in seq_len(m - 1)) {
    for (j in (i + 1):m) {
      Aij <- rbind(A_leq[i, ], A_leq[j, ])
      bij <- c(b_leq[i], b_leq[j])

      if (abs(det(Aij)) < tol) next

      v <- as.numeric(solve(Aij, bij))

      if (all(A_leq %*% v <= b_leq + tol)) {
        vertices[[length(vertices) + 1]] <- v
      }
    }
  }

  if (length(vertices) == 0) {
    return(matrix(numeric(0), ncol = 2))
  }

  V <- unique(round(do.call(rbind, vertices), 10))

  center <- colMeans(V)
  angle <- atan2(V[, 2] - center[2], V[, 1] - center[1])

  V[order(angle), , drop = FALSE]
}


make_region_polygon_data <- function(N_plot,
                                     K_list_plot,
                                     grid_plot = "even",
                                     G = 10,
                                     beta0 = 1,
                                     beta1 = 2,
                                     level = 0.95,
                                     B_gauss = 10000,
                                     seed = 123) {

  set.seed(seed)

  df_plot <- generate_data(
    N = N_plot,
    G = G,
    beta0 = beta0,
    beta1 = beta1
  )

  X_grid <- seq(-1, 1, length.out = G)
  X_true <- stats::model.matrix(~ X_grid)

  poly_list <- lapply(K_list_plot, function(K_plot) {

    tau_plot <- make_grid(K_plot, grid_plot)

    summaries_plot <- make_disclosed_quantiles(
      df = df_plot,
      tau = tau_plot
    )

    region_plot <- qlm_region(
      qhat = summaries_plot$qhat,
      tau = tau_plot,
      X = summaries_plot$X,
      n_g = summaries_plot$n_g,
      y_lower = summaries_plot$y_lower,
      y_upper = summaries_plot$y_upper,
      level = level,
      B = B_gauss
    )

    true_bounds_plot <- true_bounds_truncnorm(
      X_grid = X_grid,
      tau = tau_plot,
      beta0 = beta0,
      beta1 = beta1
    )

    vertices_true_plot <- get_vertices_2d(
      X = X_true,
      lower = true_bounds_plot$m_lower,
      upper = true_bounds_plot$m_upper
    )

    vertices_conf_plot <- region_vertices_2d(region_plot)

    rbind(
      polygon_df(
        vertices_conf_plot,
        set_name = "Confidence region",
        N = N_plot,
        K = K_plot
      ),
      polygon_df(
        vertices_true_plot,
        set_name = "True identified set",
        N = N_plot,
        K = K_plot
      )
    )
  })

  do.call(rbind, poly_list)
}




make_region_polygon_data_NK <- function(N_list_plot,
                                        K_list_plot = c(3, 5, 7, 10, 15, 20, 50),
                                        grid_plot = "even",
                                        G = 10,
                                        beta0 = 1,
                                        beta1 = 2,
                                        level = 0.95,
                                        B_gauss = 10000,
                                        seed = 123) {

  do.call(
    rbind,
    lapply(N_list_plot, function(N_plot) {
      make_region_polygon_data(
        N_plot = N_plot,
        K_list_plot = K_list_plot,
        grid_plot = grid_plot,
        G = G,
        beta0 = beta0,
        beta1 = beta1,
        level = level,
        B_gauss = B_gauss,
        seed = seed
      )
    })
  )
}


plot_region_N_by_K <- function(N_list_plot,
                               K_list_plot = c(3, 5, 7, 10, 15, 20, 50),
                               grid_plot = "even",
                               G = 10,
                               beta0 = 1,
                               beta1 = 2,
                               level = 0.95,
                               B_gauss = 10000,
                               seed = 123,
                               xlim = c(0.5, 1.5),
                               ylim = c(1.5, 2.5)) {

  df_poly <- make_region_polygon_data_NK(
    N_list_plot = N_list_plot,
    K_list_plot = K_list_plot,
    grid_plot = grid_plot,
    G = G,
    beta0 = beta0,
    beta1 = beta1,
    level = level,
    B_gauss = B_gauss,
    seed = seed
  )

  df_poly$N_label <- factor(
    paste0("N = ", df_poly$N),
    levels = paste0("N = ", N_list_plot)
  )

  df_poly$K_label <- factor(
    paste0("K = ", df_poly$K),
    levels = paste0("K = ", K_list_plot)
  )

  true_point <- expand.grid(
    N = N_list_plot,
    K = K_list_plot
  ) %>%
    mutate(
      beta0 = beta0,
      beta1 = beta1,
      N_label = factor(
        paste0("N = ", N),
        levels = paste0("N = ", N_list_plot)
      ),
      K_label = factor(
        paste0("K = ", K),
        levels = paste0("K = ", K_list_plot)
      )
    )

  ggplot(
    df_poly,
    aes(x = beta0, y = beta1, fill = set, color = set)
  ) +
    geom_polygon(alpha = 0.22, linewidth = 0.6) +
    geom_point(
      data = true_point,
      aes(x = beta0, y = beta1),
      inherit.aes = FALSE,
      size = 1.8
    ) +
    facet_grid(N_label ~ K_label) +
    coord_cartesian(
      xlim = xlim,
      ylim = ylim
    ) +
    scale_fill_manual(values = c(
      "Confidence region" = "#F4A6A6",
      "True identified set" = "#8FD3D3"
    )) +
    scale_color_manual(values = c(
      "Confidence region" = "#E56B6F",
      "True identified set" = "#00A6A6"
    )) +
    labs(
      x = expression(beta[0]),
      y = expression(beta[1]),
      fill = "",
      color = "",
      title = "True identified set and simultaneous confidence region",
      subtitle = paste0("Grid type = ", grid_plot)
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      strip.background = element_blank(),
      strip.text = element_text(face = "plain"),
      plot.title = element_text(face = "plain")
    )
}


############################
## Save polygon figures by N by K
############################
for(grid_type in grid_types){
  p_region_NK <- plot_region_N_by_K(
    N_list_plot = N_list,
    K_list_plot = K_list,
    grid_plot = grid_type,
    G = G,
    beta0 = 1,
    beta1 = 2,
    level = level,
    B_gauss = B_gauss,
    seed = 101,
    xlim = c(0, 3),
    ylim = c(0, 4)
  )

  ggsave(
    filename = paste0("results_compact_truncnorm_qlm_partial/plot_region_polygon_N_by_K_",grid_type,".png"),
    plot = p_region_NK,
    width = 10,
    height = 10,
    dpi = 300
  )
}




############################
## Meta information
############################
actual_M <- max(design_data$rep)
meta_txt <- sprintf(
  "Partial identification simulation summary (truncated-normal errors, a = -2, b = 3)\n------------------------------------------\nReplications: %d\nGaussian draws: %d\nCores: %d\nElapsed time (seconds): %.2f\nStart time: %s\nEnd time: %s\n",
  actual_M,
  B_gauss,
  n_cores,
  elapsed_time,
  format(start_time, "%Y-%m-%d %H:%M:%S"),
  format(end_time, "%Y-%m-%d %H:%M:%S")
)

writeLines(meta_txt, "results_compact_truncnorm_qlm_partial/meta_info.txt")
