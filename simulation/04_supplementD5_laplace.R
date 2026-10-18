########################
# Quantile-Only Regression: Inference and Partial Identification under Quantile Disclosure
# Simulation in Supplement D5
############################

rm(list = ls())
library(qonlyreg)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(doParallel)
library(foreach)
library(doRNG)
library(patchwork)

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


generate_data <- function(N, G = 10, beta0 = 1, beta1 = 2) {
  X_grid <- seq(-1, 1, length.out = G)
  g <- rep(1:G, length.out = N)
  Xg <- X_grid[g]

  sigma_g <- exp(0.5 * Xg)
  ## Laplace error with Var(u | Xg) = sigma_g^2
  e <- rexp(N, rate = 1) - rexp(N, rate = 1)  # Laplace(0, 1), variance = 2
  u <- sigma_g * e / sqrt(2)
  Y <- beta0 + beta1 * Xg + u

  data.frame(Y = Y, g = g, Xg = Xg)
}

make_disclosed_quantiles <- function(df, tau) {

  group_tbl <- df %>%
    dplyr::group_by(g, Xg) %>%
    dplyr::summarise(
      qhat = list(as.numeric(stats::quantile(Y, probs = tau, type = 1))),
      n_g = dplyr::n(),
      .groups = "drop"
    )

  list(
    qhat = do.call(rbind, group_tbl$qhat),
    X = stats::model.matrix(~ Xg, data = group_tbl),
    n_g = group_tbl$n_g,
    group_tbl = group_tbl
  )
}

############################
## Simulation design
############################

N_list <- c(1000, 2000, 5000, 10000, 20000)
K_list <- c(3, 5, 7, 10, 15, 20, 50)
G <- 10
grid_types <- c("even", "right_tail_coarse_mild", "right_tail_coarse_severe")
M <- 10000
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
  .packages = c("dplyr", "qonlyreg"),
  .export = c(
    "generate_data",
    "make_grid",
    "make_disclosed_quantiles",
    "K_list",
    "grid_types"
  )
) %dorng% {

  block_out <- lapply(b, function(i) {

    N_i <- design_data$N[i]
    rep_i <- design_data$rep[i]

    df_i <- generate_data(N = N_i)

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

            fit <- qlm(
              qhat = summaries$qhat,
              tau = tau,
              X = summaries$X,
              n_g = summaries$n_g
            )

            beta_hat <- unname(coef(fit))
            se_beta <- unname(sqrt(diag(vcov(fit))))

            beta_hat_0 <- beta_hat[1]
            beta_hat_1 <- beta_hat[2]
            se_beta1 <- se_beta[2]

            ci_lower_beta1 <- unname(beta_hat_1 - 1.96 * se_beta1)
            ci_upper_beta1 <- unname(beta_hat_1 + 1.96 * se_beta1)

            out <- c(
              beta0_hat = beta_hat_0,
              beta1_hat = beta_hat_1,
              se_beta1 = se_beta1,
              ci_lower_beta1 = ci_lower_beta1,
              ci_upper_beta1 = ci_upper_beta1,
              cover_beta1 = as.numeric(
                ci_lower_beta1 <= 2 & 2 <= ci_upper_beta1
              ),
              error_beta1 = beta_hat_1 - 2,
              abs_error_beta1 = abs(beta_hat_1 - 2),
              sq_error_beta1 = (beta_hat_1 - 2)^2,
              Delta = max(c(tau[1], diff(tau), 1 - tau[length(tau)]))
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


dir.create("results_nonsmooth", showWarnings = FALSE)
############################
## Summary table
# results <- readRDS("results_nonsmooth/results.rds")
############################

summary_tbl <- results %>%
  group_by(N, K, grid_type) %>%
  summarise(
    Delta = mean(Delta),
    bias_beta1 = mean(error_beta1),
    abs_bias_beta1 = abs(mean(error_beta1)),
    empirical_sd_beta1 = sd(beta1_hat),
    mean_plugin_se_beta1 = mean(se_beta1, na.rm = TRUE),
    mse_beta1 = mean(sq_error_beta1),
    mae_beta1 = mean(abs_error_beta1),
    coverage_beta1 = mean(cover_beta1, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_tbl)
write.csv(summary_tbl, "results_nonsmooth/summary_table.csv", row.names = FALSE)

library(knitr)
library(kableExtra)

tex_table <- summary_tbl %>%
  mutate(across(where(is.numeric), ~ round(., 4))) %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    caption = "Simulation results",
    label = "tab:simulation_results",
    linesep = ""
  ) %>%
  kable_styling(latex_options = c("hold_position", "scale_down"))

cat(tex_table, file = "results_nonsmooth/summary_table.tex")

############################
## Data for plots
############################

facet_levels <- summary_tbl %>%
  distinct(N) %>%
  arrange(N) %>%
  mutate(
    n_g = N / G,
    facet_label = paste0("N = ", N, ", n_g = ", n_g)
  ) %>%
  pull(facet_label)

summary_tbl_plot <- summary_tbl %>%
  mutate(
    Delta2 = Delta^2,
    n_g = N / G,
    facet_label_raw = paste0("N = ", N, ", n_g = ", n_g),
    facet_label = factor(facet_label_raw, levels = facet_levels)
  )

facet_theme <- theme(
  strip.background = element_blank(),
  strip.text = element_text(face = "plain")
)


############################
## Common scales
############################

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

legend_theme <- theme(
  legend.position = "bottom"
)

############################
## Plot 1: Bias vs K
############################

p1 <- summary_tbl_plot %>%
  arrange(N, grid_type, K) %>%
  ggplot(aes(
    x = K,
    y = abs_bias_beta1,
    color = grid_type,
    linetype = grid_type,
    shape = grid_type
  )) +
  geom_line(aes(group = grid_type)) +
  geom_point(size = 2) +
  facet_wrap(~ facet_label, scales = "free_y") +
  labs(
    x = "Number of disclosed quantiles K",
    y = "Absolute bias of beta1",
    color = "Grid type",
    linetype = "Grid type",
    shape = "Grid type"
  ) +
  grid_scales +
  theme_bw() +
  facet_theme +
  legend_theme

ggsave("results_nonsmooth/plot1_bias.png", p1, width = 8, height = 6, dpi = 300)

p1_fixed <- summary_tbl_plot %>%
  arrange(N, grid_type, K) %>%
  ggplot(aes(
    x = K,
    y = abs_bias_beta1,
    color = grid_type,
    linetype = grid_type,
    shape = grid_type
  )) +
  geom_line(aes(group = grid_type)) +
  geom_point(size = 2) +
  facet_wrap(~ facet_label, scales = "fixed") +
  labs(
    x = "Number of disclosed quantiles K",
    y = "Absolute bias of beta1",
    color = "Grid type",
    linetype = "Grid type",
    shape = "Grid type"
  ) +
  grid_scales +
  theme_bw() +
  facet_theme +
  legend_theme

ggsave("results_nonsmooth/plot1_bias_fixed_y.png", p1_fixed, width = 8, height = 6, dpi = 300)

############################
## Plot 2: MSE vs K
############################

p2 <- summary_tbl_plot %>%
  arrange(N, grid_type, K) %>%
  ggplot(aes(
    x = K,
    y = mse_beta1,
    color = grid_type,
    linetype = grid_type,
    shape = grid_type
  )) +
  geom_line(aes(group = grid_type)) +
  geom_point(size = 2) +
  facet_wrap(~ facet_label, scales = "free_y") +
  labs(
    x = "Number of disclosed quantiles K",
    y = "MSE of beta1",
    color = "Grid type",
    linetype = "Grid type",
    shape = "Grid type"
  ) +
  grid_scales +
  theme_bw() +
  facet_theme +
  legend_theme

ggsave("results_nonsmooth/plot2_mse.png", p2, width = 8, height = 6, dpi = 300)

p2_fixed <- summary_tbl_plot %>%
  arrange(N, grid_type, K) %>%
  ggplot(aes(
    x = K,
    y = mse_beta1,
    color = grid_type,
    linetype = grid_type,
    shape = grid_type
  )) +
  geom_line(aes(group = grid_type)) +
  geom_point(size = 2) +
  facet_wrap(~ facet_label, scales = "fixed") +
  labs(
    x = "Number of disclosed quantiles K",
    y = "MSE of beta1",
    color = "Grid type",
    linetype = "Grid type",
    shape = "Grid type"
  ) +
  grid_scales +
  theme_bw() +
  facet_theme +
  legend_theme

ggsave("results_nonsmooth/plot2_mse_fixed_y.png", p2_fixed, width = 8, height = 6, dpi = 300)

############################
## Plot 3: Bias vs Delta^2
############################

p3 <- summary_tbl_plot %>%
  arrange(N, grid_type, Delta2) %>%
  ggplot(aes(
    x = Delta2,
    y = abs_bias_beta1,
    color = grid_type,
    linetype = grid_type,
    shape = grid_type
  )) +
  geom_line(aes(group = grid_type)) +
  geom_point(size = 2) +
  scale_x_log10() +
  facet_wrap(~ facet_label, scales = "free_y") +
  labs(
    x = expression(log(Delta(T[K])^2)),
    y = "Absolute bias of beta1",
    color = "Grid type",
    linetype = "Grid type",
    shape = "Grid type"
  ) +
  grid_scales +
  theme_bw() +
  facet_theme +
  legend_theme

ggsave("results_nonsmooth/plot3_delta.png", p3, width = 8, height = 6, dpi = 300)

p3_fixed <- summary_tbl_plot %>%
  arrange(N, grid_type, Delta2) %>%
  ggplot(aes(
    x = Delta2,
    y = abs_bias_beta1,
    color = grid_type,
    linetype = grid_type,
    shape = grid_type
  )) +
  geom_line(aes(group = grid_type)) +
  geom_point(size = 2) +
  scale_x_log10() +
  facet_wrap(~ facet_label, scales = "fixed") +
  labs(
    x = expression(log(Delta(T[K])^2)),
    y = "Absolute bias of beta1",
    color = "Grid type",
    linetype = "Grid type",
    shape = "Grid type"
  ) +
  grid_scales +
  theme_bw() +
  facet_theme +
  legend_theme

ggsave("results_nonsmooth/plot3_delta_fixed_y.png", p3_fixed, width = 8, height = 6, dpi = 300)

############################
## Plot 4: Coverage vs K
############################

p4 <- summary_tbl_plot %>%
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
  geom_hline(yintercept = 0.95, linetype = "longdash", color = "black") +
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

ggsave("results_nonsmooth/plot4_coverage.png", p4, width = 8, height = 6, dpi = 300)

############################
## Plot 5: Empirical SD vs plug-in SE
############################

p5 <- summary_tbl_plot %>%
  ggplot(aes(
    x = empirical_sd_beta1,
    y = mean_plugin_se_beta1,
    color = grid_type,
    shape = grid_type
  )) +
  geom_point(size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "longdash", color = "black") +
  facet_wrap(~ facet_label) +
  labs(
    x = "Empirical SD of beta1_hat",
    y = "Mean plug-in SE",
    color = "Grid type",
    shape = "Grid type"
  ) +
  grid_scales[c(1, 3)] +
  theme_bw() +
  facet_theme +
  legend_theme

ggsave("results_nonsmooth/plot5_sd_vs_se.png", p5, width = 8, height = 6, dpi = 300)



#####################
## Figure for Supplemental Material D.5
#####################
N_main <- c(2000, 20000)

summary_tbl_main <- summary_tbl_plot %>%
  filter(N %in% N_main)


p1_main <- summary_tbl_main %>%
  arrange(N, grid_type, K) %>%
  ggplot(aes(x = K,
             y = abs_bias_beta1,
             color = grid_type,
             linetype = grid_type,
             shape = grid_type
  )) +
  geom_line(aes(group = grid_type)) +
  geom_point(size = 2) +
  facet_wrap(~ facet_label, scales = "fixed") +
  labs(
    x = "Number of disclosed quantiles K",
    y = "Absolute bias of beta1",
    color = "Grid type",
    linetype = "Grid type",
    shape = "Grid type"
  ) +
  grid_scales +
  theme_bw() +
  facet_theme +
  theme(legend.position = "bottom")

ggsave("results_nonsmooth/plot_main_bias.png", p1_main, width = 7, height = 4, dpi = 300)

p2_main <- summary_tbl_main %>%
  arrange(N, grid_type, K) %>%
  ggplot(aes(
    x = K,
    y = mse_beta1,
    color = grid_type,
    linetype = grid_type,
    shape = grid_type
  )) +
  geom_line(aes(group = grid_type)) +
  geom_point(size = 2) +
  facet_wrap(~ facet_label, scales = "fixed") +
  labs(
    x = "Number of disclosed quantiles K",
    y = "MSE of beta1",
    color = "Grid type",
    linetype = "Grid type",
    shape = "Grid type"
  ) +
  grid_scales +
  theme_bw() +
  facet_theme +
  theme(legend.position = "bottom")

ggsave("results_nonsmooth/plot_main_mse.png", p2_main, width = 7, height = 4, dpi = 300)


p4_main <- summary_tbl_main %>%
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
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
  facet_wrap(~ facet_label) +
  labs(
    x = "Number of disclosed quantiles K",
    y = "Coverage probability",
    color = "Grid type",
    linetype = "Grid type",
    shape = "Grid type"
  ) +
  grid_scales +
  theme_bw() +
  facet_theme +
  theme(legend.position = "bottom")

ggsave("results_nonsmooth/plot_main_coverage.png", p4_main, width = 7, height = 4, dpi = 300)


p3_main <- summary_tbl_main %>%
  arrange(N, grid_type, Delta2) %>%
  ggplot(aes(
    x = Delta2,
    y = abs_bias_beta1,
    color = grid_type,
    linetype = grid_type,
    shape = grid_type
  )) +
  geom_point(size = 2) +
  geom_line(aes(group = grid_type)) +
  scale_x_log10() +
  facet_wrap(~ facet_label, scales = "fixed") +
  labs(
    x = expression(log(Delta(T[K])^2)),
    y = "Absolute bias of beta1",
    color = "Grid type",
    linetype = "Grid type",
    shape = "Grid type"
  ) +
  grid_scales +
  theme_bw() +
  facet_theme +
  theme(legend.position = "bottom")

ggsave("results_nonsmooth/plot_main_delta.png", p3_main, width = 7, height = 4, dpi = 300)



library(patchwork)
p_combined <- (
  p1_main + p2_main
) / p4_main +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("results_nonsmooth/plot_main_combined.png", p_combined, width = 8, height = 7, dpi = 300)

# coverage probability
coverage_even_tbl <- summary_tbl %>%
  filter(grid_type == "even", K %in% K_list, N %in% N_main) %>%
  select(N, K, coverage_beta1) %>%
  mutate(
    coverage_beta1 = sprintf("%.3f", coverage_beta1)
  ) %>%
  pivot_wider(
    names_from = K,
    values_from = coverage_beta1,
    names_prefix = "K = "
  ) %>%
  arrange(N)

coverage_even_tex <- coverage_even_tbl %>%
  kable(
    format = "latex",
    booktabs = TRUE,
    escape = FALSE,
    caption = "Coverage probability under evenly spaced grids",
    label = "coverage_even_nonsmooth",
    align = "c",
    linesep = ""
  ) %>%
  kable_styling(
    latex_options = c("hold_position"),
    font_size = 9
  )

cat(coverage_even_tex, file = "results_nonsmooth/table_coverage_even.tex")



############################
## Meta information and RDS files
############################

meta_txt <- sprintf(
  "Simulation summary\n------------------\nReplications: %d\nCores: %d\nElapsed time (seconds): %.2f\nStart time: %s\nEnd time: %s\n",
  M,
  n_cores,
  elapsed_time,
  format(start_time, "%Y-%m-%d %H:%M:%S"),
  format(end_time, "%Y-%m-%d %H:%M:%S")
)

writeLines(meta_txt, "results_nonsmooth/meta_info.txt")

saveRDS(results, "results_nonsmooth/results.rds")
saveRDS(summary_tbl, "results_nonsmooth/summary_tbl.rds")

######################
## Grid points
######################

K_grid_plot <- 50

grid_df <- do.call(rbind, lapply(grid_types, function(type) {
  tau <- make_grid(K_grid_plot, type)
  data.frame(
    tau = tau,
    index = 1:K_grid_plot,
    grid_type = type
  )
}))

p_grid <- ggplot(grid_df, aes(x = tau, y = grid_type, color = grid_type)) +
  geom_point(size = 2) +
  scale_y_discrete(limits = rev(grid_types)) +
  labs(
    x = expression(tau),
    y = "",
    title = paste("Quantile grid (K =", K_grid_plot, ")")
  ) +
  theme_bw() +
  theme(legend.position = "none")

ggsave("results_nonsmooth/plot_grid_points_K50.png", p_grid, width = 8, height = 3, dpi = 300)

