# Create the manuscript and supplemental tables for the empirical illustration.
# All reported calculations use the saved analysis data and estimates produced
# by R/01_build_analysis_data.R and R/02_run_qonlyreg_rd.R.

project_root <- normalizePath("..", mustWork = TRUE)

analysis_path <- file.path(
  project_root, "data", "derived", "children_in_need_labour_rd_all_years.csv"
)
point_path <- file.path(
  project_root, "output", "qonlyreg_labour_fixed_window_point_estimates.csv"
)
region_path <- file.path(
  project_root,
  "output", "qonlyreg_labour_fixed_window_simultaneous_region_summary.csv"
)
output_dir <- file.path(project_root, "output", "paper_tables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(analysis_path, point_path, region_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop(
    "Required inputs are missing. Run R/01_build_analysis_data.R and ",
    "R/02_run_qonlyreg_rd.R first. Missing: ",
    paste(missing_files, collapse = ", ")
  )
}

analysis <- read.csv(analysis_path, stringsAsFactors = FALSE)
point <- read.csv(point_path, stringsAsFactors = FALSE)
region <- read.csv(region_path, stringsAsFactors = FALSE)

logical_columns <- c(
  "analysis_eligible", "labour_at_least_50", "labour_strict_majority",
  "labour_exactly_50"
)
for (nm in logical_columns) analysis[[nm]] <- as.logical(analysis[[nm]])

window_order <- c(3, 4, 5)
window_roles <- c(
  `3` = "Narrower-window sensitivity",
  `4` = "Main specification",
  `5` = "Wider-window sensitivity"
)
upper_bounds <- c(150, 175, 200)

format_number <- function(x, digits = 2) {
  ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
}

format_integer <- function(x) {
  ifelse(is.na(x), "", formatC(x, format = "d", big.mark = ","))
}

format_interval <- function(lower, upper, digits = 2) {
  paste0("[", format_number(lower, digits), ", ", format_number(upper, digits), "]")
}

# Escape table cell text before inserting it into a LaTeX tabular environment.
latex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "\\\\&", x, fixed = TRUE)
  x <- gsub("%", "\\\\%", x, fixed = TRUE)
  x <- gsub("#", "\\\\#", x, fixed = TRUE)
  x <- gsub("_", "\\\\_", x, fixed = TRUE)
  x
}

write_latex_table <- function(
    path, caption, label, headers, rows, alignment, note = NULL,
    font_command = "\\small") {
  stopifnot(length(headers) == ncol(rows), nchar(alignment) == ncol(rows))
  lines <- c(
    "\\begin{table}[!h]",
    "\\centering",
    "\\begin{threeparttable}",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    font_command,
    paste0("\\begin{tabular}{", alignment, "}"),
    "\\toprule",
    paste0(paste(headers, collapse = " & "), " \\\\"),
    "\\midrule"
  )
  body <- apply(rows, 1, function(z) {
    paste0(paste(latex_escape(z), collapse = " & "), " \\\\")
  })
  lines <- c(lines, body, "\\bottomrule", "\\end{tabular}")
  if (!is.null(note) && nzchar(note)) {
    lines <- c(
      lines,
      "\\begin{tablenotes}[flushleft]",
      "\\footnotesize",
      paste0("\\item Notes: ", note),
      "\\end{tablenotes}"
    )
  }
  lines <- c(lines, "\\end{threeparttable}", "\\end{table}")
  writeLines(lines, path, useBytes = TRUE)
}

assemble_results_table <- function(treatment_definition, model_specification) {
  p <- point[
    point$treatment_definition == treatment_definition &
      point$model_specification == model_specification &
      point$bandwidth_pp %in% window_order,
    ,
    drop = FALSE
  ]
  p <- p[match(window_order, p$bandwidth_pp), , drop = FALSE]
  if (nrow(p) != length(window_order) || anyNA(p$estimate_days)) {
    stop("Point-estimate results are incomplete for ", treatment_definition, ".")
  }

  out <- data.frame(
    bandwidth_pp = p$bandwidth_pp,
    role = unname(window_roles[as.character(p$bandwidth_pp)]),
    council_years = p$observations,
    unique_councils = p$unique_councils,
    treated = p$treated,
    control = p$control,
    exactly_50 = p$exactly_50,
    estimate_days = p$estimate_days,
    qlm_standard_error = p$qlm_standard_error,
    stringsAsFactors = FALSE
  )

  for (support_upper in upper_bounds) {
    r <- region[
      region$treatment_definition == treatment_definition &
        region$model_specification == model_specification &
        region$bandwidth_pp %in% window_order &
        region$support_upper_days == support_upper,
      ,
      drop = FALSE
    ]
    r <- r[match(window_order, r$bandwidth_pp), , drop = FALSE]
    if (nrow(r) != length(window_order)) {
      stop("Simultaneous-region results are incomplete for ", treatment_definition, ".")
    }
    out[[paste0("simultaneous_95_lower_U", support_upper)]] <-
      r$treatment_projection_lower_95
    out[[paste0("simultaneous_95_upper_U", support_upper)]] <-
      r$treatment_projection_upper_95
  }
  out
}

format_results_for_latex <- function(x) {
  data.frame(
    Window = paste0("$\\pm ", format_number(x$bandwidth_pp, 0), "$"),
    `Council-years` = format_integer(x$council_years),
    Councils = format_integer(x$unique_councils),
    `Treated / control` = paste0(x$treated, " / ", x$control),
    `Estimate (SE)` = paste0(
      format_number(x$estimate_days, 3), " (",
      format_number(x$qlm_standard_error, 3), ")"
    ),
    `U=150` = format_interval(
      x$simultaneous_95_lower_U150, x$simultaneous_95_upper_U150
    ),
    `U=175` = format_interval(
      x$simultaneous_95_lower_U175, x$simultaneous_95_upper_U175
    ),
    `U=200` = format_interval(
      x$simultaneous_95_lower_U200, x$simultaneous_95_upper_U200
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

# Main table: fixed-window estimates and simultaneous confidence regions.
main_results <- assemble_results_table(
  "at_least_50_main", "no_year_fixed_effects"
)
write.csv(
  main_results,
  file.path(output_dir, "table_main_results.csv"),
  row.names = FALSE,
  na = ""
)
main_latex <- format_results_for_latex(main_results)
write_latex_table(
  file.path(output_dir, "table_main_results.tex"),
  caption = paste0(
    "Labour majority and assessment duration: fixed-window qonlyreg results"
  ),
  label = "tab:qonlyreg_main",
  headers = c(
    "$h$", "$|\\mathcal G_h|$", "Co.", "T / C",
    "Est. (SE)", "$U=150$", "$U=175$", "$U=200$"
  ),
  rows = main_latex,
  alignment = "cccccccc",
  note = paste0(
    "$h$ is the window half-width in percentage points; ",
    "$|\\mathcal G_h|$ and Co.\\ are the numbers of council-year groups ",
    "and councils; T/C gives treated/control group counts. ",
    "Est.\\ (SE) reports the Labour-majority coefficient (standard error). ",
    "Brackets give projections of the 95\\% simultaneous confidence ",
    "region under $Q_g(0)=0$ and upper support bound $U$ (days). ",
    "The main specification uses $h=4$."
  )
)

# Supplemental sensitivity table: treatment begins strictly above 50 percent.
strict_results <- assemble_results_table(
  "strictly_above_50_sensitivity", "no_year_fixed_effects"
)
write.csv(
  strict_results,
  file.path(output_dir, "table_strict_majority_sensitivity.csv"),
  row.names = FALSE,
  na = ""
)
strict_latex <- format_results_for_latex(strict_results)
write_latex_table(
  file.path(output_dir, "table_strict_majority_sensitivity.tex"),
  caption = "Sensitivity to defining treatment as strictly above 50 percent",
  label = "tab:qonlyreg_strict_majority",
  headers = c(
    "Window", "$|\\mathcal G_h|$", "Councils", "T / C",
    "Estimate (SE)", "$U=150$", "$U=175$", "$U=200$"
  ),
  rows = strict_latex,
  alignment = "cccccccc",
  note = paste0(
    "The treatment is a Labour seat share strictly greater than 50 percent, so ",
    "observations exactly at 50 percent are assigned to control. All other ",
    "definitions match Table~\\ref{tab:qonlyreg_main}."
  )
)

# Appendix sensitivity table: add reporting-year fixed effects.
year_fe_results <- assemble_results_table(
  "at_least_50_main", "reporting_year_fixed_effects_sensitivity"
)
write.csv(
  year_fe_results,
  file.path(output_dir, "table_year_fixed_effects_sensitivity.csv"),
  row.names = FALSE,
  na = ""
)
year_fe_latex <- format_results_for_latex(year_fe_results)
write_latex_table(
  file.path(output_dir, "table_year_fixed_effects_sensitivity.tex"),
  caption = "Sensitivity to reporting-year fixed effects",
  label = "tab:qonlyreg_year_fixed_effects",
  headers = c(
    "Window", "$|\\mathcal G_h|$", "Councils", "T / C",
    "Estimate (SE)", "$U=150$", "$U=175$", "$U=200$"
  ),
  rows = year_fe_latex,
  alignment = "cccccccc",
  note = paste0(
    "The treatment is a Labour seat share of at least 50 percent. Each model ",
    "adds reporting-year indicators to the intercept and treatment indicator; ",
    "2017 is the omitted reporting year. The quantile reconstruction, fixed ",
    "windows, assessment-count weights, and support assumptions match ",
    "Table~\\ref{tab:qonlyreg_main}. Parentheses contain qlm standard errors, ",
    "and brackets report projections of the 95 percent simultaneous confidence ",
    "region."
  )
)

# Decompose the main-window reconstructed-mean difference by disclosed quantile.
main_data <- analysis[
  analysis$analysis_eligible %in% TRUE &
    !is.na(analysis$running_variable_pp) &
    abs(analysis$running_variable_pp) <= 4 + 1e-12,
  ,
  drop = FALSE
]
main_data <- main_data[order(main_data$reporting_year, main_data$new_la_code), ]
if (!nrow(main_data) || !all(c(FALSE, TRUE) %in% main_data$labour_at_least_50)) {
  stop("The main window does not contain observations on both sides of 50 percent.")
}

weighted_side_average <- function(values, treated) {
  stats::weighted.mean(
    values[main_data$labour_at_least_50 == treated],
    main_data$assessments_completed[main_data$labour_at_least_50 == treated]
  )
}

component_names <- c(
  "Known lower endpoint Q(0)", "Q1 (25th percentile)",
  "Median (50th percentile)", "Q3 (75th percentile)"
)
component_variables <- c(NA, "q25_days", "median_days", "q75_days")
reconstruction_weights <- c(0.125, 0.250, 0.250, 0.375)
control_average <- numeric(length(component_names))
treated_average <- numeric(length(component_names))
for (j in seq_along(component_names)) {
  if (is.na(component_variables[j])) {
    control_average[j] <- 0
    treated_average[j] <- 0
  } else {
    values <- main_data[[component_variables[j]]]
    control_average[j] <- weighted_side_average(values, FALSE)
    treated_average[j] <- weighted_side_average(values, TRUE)
  }
}
difference <- treated_average - control_average
contribution <- reconstruction_weights * difference

control_mean <- weighted_side_average(main_data$qonlyreg_mean_days, FALSE)
treated_mean <- weighted_side_average(main_data$qonlyreg_mean_days, TRUE)
mean_difference <- treated_mean - control_mean
if (abs(sum(contribution) - mean_difference) > 1e-10) {
  stop("Quantile contributions do not reproduce the reconstructed-mean difference.")
}

quantile_decomposition <- data.frame(
  component = c(component_names, "Reconstructed mean"),
  reconstruction_weight = c(reconstruction_weights, 1),
  control_weighted_average_days = c(control_average, control_mean),
  treated_weighted_average_days = c(treated_average, treated_mean),
  treated_minus_control_days = c(difference, mean_difference),
  contribution_to_mean_difference_days = c(contribution, mean_difference),
  stringsAsFactors = FALSE
)
write.csv(
  quantile_decomposition,
  file.path(output_dir, "table_quantile_decomposition.csv"),
  row.names = FALSE,
  na = ""
)
quantile_latex <- data.frame(
  Component = quantile_decomposition$component,
  Weight = format_number(quantile_decomposition$reconstruction_weight, 3),
  Control = format_number(quantile_decomposition$control_weighted_average_days, 2),
  Treated = format_number(quantile_decomposition$treated_weighted_average_days, 2),
  Difference = format_number(quantile_decomposition$treated_minus_control_days, 2),
  Contribution = format_number(
    quantile_decomposition$contribution_to_mean_difference_days, 2
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write_latex_table(
  file.path(output_dir, "table_quantile_decomposition.tex"),
  caption = "Decomposition of the reconstructed-mean difference in the main window",
  label = "tab:qonlyreg_quantile_decomposition",
  headers = c(
    "Component", "Weight", "Control", "Treated", "Difference", "Contribution"
  ),
  rows = quantile_latex,
  alignment = "cccccc",
  note = paste0(
    "Control and treated columns are assessment-count-weighted side averages ",
    "within the $\\pm4$ percentage-point window. Difference is treated minus ",
    "control. Contribution equals the reconstruction weight times the difference. ",
    "The known lower endpoint contributes zero by construction."
  )
)

# Re-estimate the main-window weighted difference after excluding each council.
estimate_weighted_difference <- function(d) {
  treated <- d$labour_at_least_50
  if (!all(c(FALSE, TRUE) %in% treated)) return(NA_real_)
  stats::weighted.mean(
    d$qonlyreg_mean_days[treated], d$assessments_completed[treated]
  ) - stats::weighted.mean(
    d$qonlyreg_mean_days[!treated], d$assessments_completed[!treated]
  )
}

baseline_estimate <- estimate_weighted_difference(main_data)
council_codes <- unique(main_data$new_la_code)
loo_rows <- lapply(council_codes, function(code) {
  removed <- main_data[main_data$new_la_code == code, , drop = FALSE]
  retained <- main_data[main_data$new_la_code != code, , drop = FALSE]
  estimate <- estimate_weighted_difference(retained)
  side <- unique(ifelse(
    removed$labour_at_least_50, "treated", "control"
  ))
  data.frame(
    excluded_council_code = code,
    excluded_council = removed$la_name[1],
    threshold_side = if (length(side) == 1) side else "both across years",
    council_years_removed = nrow(removed),
    assessments_removed = sum(removed$assessments_completed),
    leave_one_council_out_estimate_days = estimate,
    change_from_baseline_days = estimate - baseline_estimate,
    absolute_change_days = abs(estimate - baseline_estimate),
    stringsAsFactors = FALSE
  )
})
loo_full <- do.call(rbind, loo_rows)
loo_full <- loo_full[order(-loo_full$absolute_change_days, loo_full$excluded_council), ]
rownames(loo_full) <- NULL
negative_loo_count <- sum(
  loo_full$leave_one_council_out_estimate_days < 0,
  na.rm = TRUE
)
write.csv(
  loo_full,
  file.path(output_dir, "table_leave_one_council_out_full.csv"),
  row.names = FALSE,
  na = ""
)

top_count <- min(5L, nrow(loo_full))
loo_compact <- rbind(
  data.frame(
    excluded_council = "None (baseline)",
    threshold_side = "--",
    council_years_removed = 0,
    assessments_removed = 0,
    leave_one_council_out_estimate_days = baseline_estimate,
    change_from_baseline_days = 0,
    stringsAsFactors = FALSE
  ),
  loo_full[seq_len(top_count), c(
    "excluded_council", "threshold_side", "council_years_removed",
    "assessments_removed", "leave_one_council_out_estimate_days",
    "change_from_baseline_days"
  )]
)
loo_compact$excluded_council[
  loo_compact$excluded_council == "Bristol, City of"
] <- "Bristol"

write.csv(
  loo_compact,
  file.path(output_dir, "table_leave_one_council_out.csv"),
  row.names = FALSE,
  na = ""
)
loo_latex <- data.frame(
  `Excluded council` = loo_compact$excluded_council,
  Side = loo_compact$threshold_side,
  `Council-years removed` = format_integer(loo_compact$council_years_removed),
  `Assessments removed` = format_integer(loo_compact$assessments_removed),
  Estimate = format_number(
    loo_compact$leave_one_council_out_estimate_days, 3
  ),
  Change = format_number(loo_compact$change_from_baseline_days, 3),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write_latex_table(
  file.path(output_dir, "table_leave_one_council_out.tex"),
  caption = "Leave-one-council-out sensitivity in the main window",
  label = "tab:qonlyreg_leave_one_council_out",
  headers = c(
    "Excluded council", "Side", "Council-years", "Assessments", "Estimate", "Change"
  ),
  rows = loo_latex,
  alignment = "cccccc",
  note = paste0(
    "The first row is the baseline $\\pm4$ estimate. The remaining rows show ",
    "the five councils producing the largest absolute change when all of that ",
    "council's observations are removed. Change is relative to the baseline. ",
    "Across all ", nrow(loo_full), " exclusions, ", negative_loo_count,
    " estimates remain negative. ",
    "The complete council-level results are supplied in ",
    "\\texttt{table\\_leave\\_one\\_council\\_out\\_full.csv}."
  )
)

# Save a small audit file describing the manuscript placement of each output.
manifest <- data.frame(
  manuscript_location = c(
    "Main text", "Main text", "Supplement", "Supplement", "Supplement", "Supplement",
    "Reproducibility supplement"
  ),
  artifact = c(
    "table_main_results.csv and .tex",
    "rd_quantiles_and_reconstructed_mean_main_window_4pp.png and .pdf",
    "rd_quantiles_main_window_4pp.png and .pdf",
    "table_leave_one_council_out.csv and .tex",
    "table_year_fixed_effects_sensitivity.csv and .tex",
    "table_strict_majority_sensitivity.csv and .tex",
    "table_quantile_decomposition.csv and .tex; table_leave_one_council_out_full.csv"
  ),
  purpose = c(
    "Primary point estimates and simultaneous confidence regions",
    "All disclosed quantiles and the reconstructed mean",
    "Quantile-specific diagnostic panels",
    "Council-level influence sensitivity",
    "Reporting-year fixed-effects sensitivity",
    "Treatment-definition sensitivity",
    "Numerical audit and complete influence results"
  ),
  stringsAsFactors = FALSE
)
write.csv(
  manifest,
  file.path(output_dir, "paper_output_manifest.csv"),
  row.names = FALSE
)

message("Paper tables written to: ", output_dir)
message(
  "Main-window estimate: ", format_number(baseline_estimate, 3),
  " days; council-level leave-one-out range: [",
  format_number(min(loo_full$leave_one_council_out_estimate_days, na.rm = TRUE), 3),
  ", ",
  format_number(max(loo_full$leave_one_council_out_estimate_days, na.rm = TRUE), 3),
  "]."
)
