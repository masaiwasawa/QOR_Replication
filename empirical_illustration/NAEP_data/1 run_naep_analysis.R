############################
## Quantile-Only Regression: Inference and Partial Identification under Quantile Disclosure
## NAEP empirical analysis
############################

rm(list = ls())

library(qonlyreg)
library(dplyr)
library(tidyr)
library(ggplot2)

############################
## 1. Load data
############################

df0 <- read.csv(
  "naep_analysis_data.csv",
  check.names = FALSE
)

df0 <- df0 %>%
  mutate(
    state = as.character(state),
    eds_status = as.character(eds_status),
    noneds = if_else(
      eds_status == "Not economically disadvantaged",
      1,
      0
    )
  ) %>%
  filter(
    eds_status %in% c(
      "Economically disadvantaged",
      "Not economically disadvantaged"
    ),
    !is.na(n_g)
  )

############################
## 2. Benchmark regression
## using published means
############################

benchmark_fit <- lm(
  mean_score ~ 0 + factor(state) + noneds,
  data = df0,
  weights = n_g
)

benchmark_beta <- coef(benchmark_fit)["noneds"]

summary(benchmark_fit)

############################
## 3. Prepare qonlyreg input
############################

make_qonly_input <- function(
    data,
    grid = c("median", "quartiles", "rich")
) {

  grid <- match.arg(grid)


  ############################
  ## disclosed quantiles
  ############################
  if (grid == "median") {
    tau <- c(0.50)
    qcols <- c("q50")
  }

  if (grid == "quartiles") {
    tau <- c(0.25, 0.50, 0.75)
    qcols <- c("q25", "q50", "q75")
  }

  if (grid == "rich") {
    tau <- c(0.10, 0.25, 0.50, 0.75, 0.90)
    qcols <- c("q10", "q25", "q50", "q75", "q90")
  }

  data2 <- data %>%
    arrange(state, eds_status)

  ############################
  ## quantile matrix
  ############################
  qhat <- as.matrix(data2[, qcols])
  storage.mode(qhat) <- "numeric"

  ############################
  ## regression design
  ############################
  X <- model.matrix(
    ~ 0 + factor(state) + noneds,
    data = data2
  )

  list(
    grid = grid,
    data = data2,
    qhat = qhat,
    tau = tau,
    X = X,
    n_g = data2$n_g,
    y_lower_bound = rep(0, nrow(data2)),
    y_upper_bound = rep(500, nrow(data2))
  )
}

############################
## 4. Main analysis function
############################

run_analysis <- function(grid) {
  
  obj <- make_qonly_input(df0, grid)
  
  beta_name <- "noneds"
  beta_idx <- which(colnames(obj$X) == beta_name)
  
  ############################
  ## Reconstruction-based point estimate
  ############################
  
  fit <- qlm(
    qhat = obj$qhat,
    tau = obj$tau,
    X = obj$X,
    n_g = obj$n_g
  )
  
  beta_hat <- unname(
    coef(fit)[beta_idx]
  )
  
  ############################
  ## Median-only case
  ##
  ## The point estimate is available, but neither
  ## the conventional CI nor the partial-ID confidence
  ## region is reported under support bounds alone.
  ############################
  
  if (length(obj$tau) == 1L) {
    
    return(
      tibble(
        grid = grid,
        beta_hat = beta_hat,
        ci_point_lower = NA_real_,
        ci_point_upper = NA_real_,
        proj_confset_lower = NA_real_,
        proj_confset_upper = NA_real_,
        benchmark_beta = benchmark_beta,
        benchmark_in_point_ci = NA,
        benchmark_in_proj_confset = NA
      )
    )
  }
  
  ############################
  ## Conventional confidence interval
  ############################
  
  ci_point <- confint(fit)
  
  ci_point_lower <- ci_point[beta_idx, 1]
  ci_point_upper <- ci_point[beta_idx, 2]
  
  ############################
  ## Partial-ID confidence region
  ############################
  
  region <- qlm_region(
    qhat = obj$qhat,
    tau = obj$tau,
    X = obj$X,
    n_g = obj$n_g,
    y_lower_bound = obj$y_lower_bound,
    y_upper_bound = obj$y_upper_bound,
    level = 0.95,
    B = 10000
  )
  
  proj_ci <- confint(region)$confint
  
  proj_confset_lower <- proj_ci[beta_idx, "lower"]
  proj_confset_upper <- proj_ci[beta_idx, "upper"]
  
  ############################
  ## Return results
  ############################
  
  tibble(
    grid = grid,
    beta_hat = beta_hat,
    ci_point_lower = ci_point_lower,
    ci_point_upper = ci_point_upper,
    proj_confset_lower = proj_confset_lower,
    proj_confset_upper = proj_confset_upper,
    benchmark_beta = benchmark_beta,
    benchmark_in_point_ci =
      benchmark_beta >= ci_point_lower &
      benchmark_beta <= ci_point_upper,
    benchmark_in_proj_confset =
      benchmark_beta >= proj_confset_lower &
      benchmark_beta <= proj_confset_upper
  )
}

############################
## 5. Run analyses
############################

res_median <- run_analysis("median")
res_quartiles <- run_analysis("quartiles")
res_rich <- run_analysis("rich")

results_all <- bind_rows(
  res_median,
  res_quartiles,
  res_rich
)

print(results_all)

############################
## 6. Save results
############################

write.csv(
  results_all,
  "naep_qonlyreg_results.csv",
  row.names = FALSE
)

saveRDS(
  results_all,
  "naep_qonlyreg_results.rds"
)

############################
## 7. Plot
############################
p_confset <- results_all %>%
  filter(!is.na(proj_confset_lower)) %>%
  ggplot(aes(x = grid)) +
  geom_errorbar(
    aes(
      ymin = proj_confset_lower,
      ymax = proj_confset_upper
    ),
    width = 0.20,
    linewidth = 1.1
  ) +
  geom_hline(yintercept = benchmark_beta, linetype = "dashed") +
  labs(
    x = "Quantile disclosure",
    y = "Coefficient on non-EDS",
    title = "Projection of simultaneous confidence region"
  ) +
  theme_bw()

p_point <- results_all %>%
  filter(!is.na(beta_hat)) %>%
  ggplot(aes(x = grid, y = beta_hat)) +
  geom_point(size = 3) +
  geom_errorbar(
    data = function(d) {
      d[!is.na(d$ci_point_lower), , drop = FALSE]
    },
    aes(
      ymin = ci_point_lower,
      ymax = ci_point_upper
    ),
    width = 0.15,
    linewidth = 0.9
  ) +
  geom_hline(yintercept = benchmark_beta, linetype = "dashed") +
  labs(
    x = "Quantile disclosure",
    y = "Coefficient on non-EDS",
    title = "Point estimates and confidence intervals"
  ) +
  theme_bw()

p_confset
p_point


library(patchwork)
p_confset / p_point


############################
## Table for paper
############################

library(dplyr)
library(knitr)
library(kableExtra)

############################
## 1. Prepare table
############################

table_results <- results_all %>%
  mutate(

    Disclosure = recode(
      grid,
      median = "Median",
      quartiles = "Quartiles",
      rich = "Five percentiles"
    ),

    point_estimate = if_else(
      is.na(beta_hat),
      "--",
      sprintf("%.1f", beta_hat)
    ),

    point_ci = if_else(
      is.na(ci_point_lower),
      "--",
      paste0(
        "[",
        sprintf("%.1f", ci_point_lower),
        ", ",
        sprintf("%.1f", ci_point_upper),
        "]"
      )
    ),

    point_covered = if_else(
      is.na(ci_point_lower),
      "--",
      if_else(
        benchmark_beta >= ci_point_lower &
          benchmark_beta <= ci_point_upper,
        "Yes",
        "No"
      )
    ),

    projection_ci = if_else(
      is.na(proj_confset_lower),
      "--",
      paste0(
        "[",
        sprintf("%.1f", proj_confset_lower),
        ", ",
        sprintf("%.1f", proj_confset_upper),
        "]"
      )
    ),
    
    projection_width = if_else(
      is.na(proj_confset_lower),
      "--",
      sprintf(
        "%.1f",
        proj_confset_upper - proj_confset_lower
      )
    ),
    
    projection_covered = if_else(
      is.na(proj_confset_lower),
      "--",
      if_else(
        benchmark_beta >= proj_confset_lower &
          benchmark_beta <= proj_confset_upper,
        "Yes",
        "No"
        )
      ) 
    )%>%
  select(
    Disclosure,
    Estimate = point_estimate,
    CI = point_ci,
    `Point covered` = point_covered,
    `Projection CI` = projection_ci,
    Width = projection_width,
    `Partial covered` = projection_covered
  )

############################
## 2. LaTeX table
############################

note_text <- paste0(
  "Benchmark estimate is ",
  sprintf("%.1f", benchmark_beta),
  ". The benchmark estimate is obtained by regressing the published group mean scores ",
  "on the non-EDS indicator and state fixed effects. “Covered” indicates whether the benchmark estimate
  based on the published mean scores is contained
  in the corresponding interval."
)

colnames(table_results) <- c(
  "Disclosure",
  "Estimate",
  "CI",
  "Covered",
  "Projection CI",
  "Width",
  "Covered"
)



tex_table <- kable(
  table_results,
  format = "latex",
  booktabs = TRUE,
  escape = FALSE,
  align = c("l", "c", "c", "c", "c", "c", "c"),
  caption = "Quantile-only regression results using NAEP data",
  label = "naep_qonlyreg"
) %>%
  add_header_above(
    c(" " = 1, "Point" = 3, "Partial" = 3)
  ) %>%
  kable_styling(
    latex_options = c("hold_position"),
    full_width = FALSE
  )

tex_table <- sub(
  "\\\\centering\n",
  "\\\\centering\n\\\\begin{threeparttable}\n",
  tex_table
)

tex_table <- sub(
  "\\\\end\\{tabular\\}",
  paste0(
    "\\\\end{tabular}\n",
    "\\\\begin{tablenotes}[flushleft]\n",
    "\\\\footnotesize\n",
    "\\\\item Notes: ",
    note_text,
    "\n\\\\end{tablenotes}"
  ),
  tex_table
)
tex_table <- sub(
  "\\\\end\\{table\\}",
  "\\\\end{threeparttable}\n\\\\end{table}",
  tex_table
)


############################
## 3. Save and display
############################

writeLines(
  tex_table,
  "table_naep_qonlyreg.tex"
)

tex_table
