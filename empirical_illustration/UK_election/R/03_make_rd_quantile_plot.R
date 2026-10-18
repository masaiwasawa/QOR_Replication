# RD visualization of the disclosed assessment-duration quartiles and the
# qonlyreg reconstructed mean.

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
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2 is required to create the figures.")
}

data <- read.csv(input_path, stringsAsFactors = FALSE)
data$analysis_eligible <- as.logical(data$analysis_eligible)
data$labour_at_least_50 <- as.logical(data$labour_at_least_50)

main_bandwidth_pp <- 4.0
plot_data <- data[
  data$analysis_eligible %in% TRUE &
    !is.na(data$running_variable_pp) &
    abs(data$running_variable_pp) <= main_bandwidth_pp + 1e-12,
  ,
  drop = FALSE
]
plot_data <- plot_data[order(plot_data$reporting_year, plot_data$new_la_code), ]

if (!nrow(plot_data)) stop("There are no observations to plot within the +/-4 percentage-point window.")

# Use the same small horizontal displacement in all three panels.  It only
# reveals overlapping council-years and does not alter the analysis variable.
set.seed(20260807L)
plot_data$running_variable_plot <- plot_data$running_variable_pp +
  stats::runif(nrow(plot_data), min = -0.035, max = 0.035)
below <- plot_data$running_variable_pp < 0
above <- plot_data$running_variable_pp > 0
exact <- abs(plot_data$running_variable_pp) <= 1e-12
plot_data$running_variable_plot[below] <- pmin(
  plot_data$running_variable_plot[below], -0.01
)
plot_data$running_variable_plot[above] <- pmax(
  plot_data$running_variable_plot[above], 0.01
)
plot_data$running_variable_plot[exact] <- stats::runif(sum(exact), 0.01, 0.06)

make_panel <- function(variable, label) {
  data.frame(
    new_la_code = plot_data$new_la_code,
    la_name = plot_data$la_name,
    reporting_year = plot_data$reporting_year,
    running_variable_pp = plot_data$running_variable_pp,
    running_variable_plot = plot_data$running_variable_plot,
    labour_at_least_50 = plot_data$labour_at_least_50,
    assessments_completed = plot_data$assessments_completed,
    quantile = label,
    duration_days = plot_data[[variable]],
    stringsAsFactors = FALSE
  )
}

long <- rbind(
  make_panel("q25_days", "25th percentile"),
  make_panel("median_days", "Median"),
  make_panel("q75_days", "75th percentile")
)
long$quantile <- factor(
  long$quantile,
  levels = c(
    "25th percentile",
    "Median",
    "75th percentile"
  )
)
long$treatment_group <- factor(
  ifelse(
    long$labour_at_least_50,
    "Labour seat share >= 50%",
    "Labour seat share < 50%"
  ),
  levels = c("Labour seat share < 50%", "Labour seat share >= 50%")
)

quantile_segments <- do.call(rbind, lapply(levels(long$quantile), function(q) {
  do.call(rbind, lapply(levels(long$treatment_group), function(group) {
    z <- long[long$quantile == q & long$treatment_group == group, ]
    is_control <- group == "Labour seat share < 50%"
    weighted_average <- stats::weighted.mean(
      z$duration_days,
      z$assessments_completed
    )
    data.frame(
      quantile = q,
      treatment_group = group,
      x = if (is_control) -main_bandwidth_pp else 0.06,
      xend = if (is_control) -0.06 else main_bandwidth_pp,
      y = weighted_average,
      yend = weighted_average,
      label_x = if (is_control) -3.86 else 0.14,
      label = sprintf("%.2f", weighted_average),
      stringsAsFactors = FALSE
    )
  }))
}))
quantile_segments$quantile <- factor(
  quantile_segments$quantile,
  levels = levels(long$quantile)
)
quantile_segments$treatment_group <- factor(
  quantile_segments$treatment_group,
  levels = levels(long$treatment_group)
)

# Supplemental three-panel display.
p_three_panel <- ggplot2::ggplot(
  long,
  ggplot2::aes(
    x = running_variable_plot,
    y = duration_days,
    size = assessments_completed,
    colour = treatment_group,
    shape = treatment_group
  )
) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.45, colour = "#333333") +
  ggplot2::geom_point(alpha = 0.66) +
  ggplot2::geom_segment(
    data = quantile_segments,
    ggplot2::aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend,
      colour = treatment_group
    ),
    inherit.aes = FALSE,
    show.legend = FALSE,
    linewidth = 1.05
  ) +
  ggplot2::geom_text(
    data = quantile_segments,
    ggplot2::aes(
      x = label_x,
      y = y,
      label = label,
      colour = treatment_group
    ),
    inherit.aes = FALSE,
    show.legend = FALSE,
    hjust = 0,
    vjust = -0.55,
    size = 2.9
  ) +
  ggplot2::facet_wrap(~quantile, nrow = 1) +
  ggplot2::scale_x_continuous(
    breaks = seq(-4, 4, by = 2),
    limits = c(-4.15, 4.15),
    expand = ggplot2::expansion(mult = 0)
  ) +
  ggplot2::scale_colour_manual(
    name = "Threshold side",
    breaks = c(
      "Labour seat share < 50%",
      "Labour seat share >= 50%"
    ),
    values = c(
      "Labour seat share < 50%" = "#F8766D",
      "Labour seat share >= 50%" = "#00BFC4"
    )
  ) +
  
  ggplot2::scale_shape_manual(
    name = "Threshold side",
    breaks = c(
      "Labour seat share < 50%",
      "Labour seat share >= 50%"
    ),
    values = c(
      "Labour seat share < 50%" = 16,
      "Labour seat share >= 50%" = 17
    )
  ) +  
  ggplot2::scale_size_area(
    max_size = 4.6,
    breaks = c(1000, 3000, 6000),
    name = "Completed assessments"
  ) +
  ggplot2::labs(
    title = "Assessment-duration quantiles around the Labour seat-share threshold",
    x = "Labour seat share minus 50 (percentage points)",
    y = "Assessment completion time (days)",
    caption = paste0(
      "Council-year observations within the main +/-4 percentage-point window. ",
      "Horizontal jitter is for display only.\n",
      "Point area is proportional to completed assessments; horizontal segments are assessment-count-weighted side averages.\n",
      "Observations exactly at 50% are assigned to the treated group. No fitted curves are imposed."
    )
  ) +
  ggplot2::theme_bw(base_size = 10.5) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "#F2F2F2", colour = "#777777"),
    strip.text = ggplot2::element_text(face = "bold"),
    plot.title = ggplot2::element_text(face = "bold", hjust = 0),
    plot.caption = ggplot2::element_text(hjust = 0, size = 8),
    legend.position = "bottom",
    legend.box = "horizontal"
  ) +
  ggplot2::guides(
    colour = ggplot2::guide_legend(
      order = 1,
      override.aes = list(size = 3, alpha = 1)
    ),
    shape = ggplot2::guide_legend(
      order = 1,
      override.aes = list(size = 3, alpha = 1)
    ),
    size = ggplot2::guide_legend(
      order = 2,
      override.aes = list(alpha = 0.66)
    )
  )

mean_points <- data.frame(
  running_variable_plot = plot_data$running_variable_plot,
  assessments_completed = plot_data$assessments_completed,
  duration_days = plot_data$qonlyreg_mean_days,
  stringsAsFactors = FALSE
)

control <- !plot_data$labour_at_least_50
treated <- plot_data$labour_at_least_50
weighted_mean_control <- stats::weighted.mean(
  plot_data$qonlyreg_mean_days[control],
  plot_data$assessments_completed[control]
)
weighted_mean_treated <- stats::weighted.mean(
  plot_data$qonlyreg_mean_days[treated],
  plot_data$assessments_completed[treated]
)
mean_segments <- data.frame(
  x = c(-main_bandwidth_pp, 0.06),
  xend = c(-0.06, main_bandwidth_pp),
  y = c(weighted_mean_control, weighted_mean_treated),
  yend = c(weighted_mean_control, weighted_mean_treated),
  label_x = c(-3.86, 0.14),
  label = sprintf("%.2f", c(weighted_mean_control, weighted_mean_treated))
)

control_assessments <- sum(plot_data$assessments_completed[control])
treated_assessments <- sum(plot_data$assessments_completed[treated])

# Figure 2: all disclosed quantiles in one panel, with the reconstructed mean
# and only its assessment-count-weighted side averages added in black.
p_combined <- ggplot2::ggplot() +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.45, colour = "#333333") +
  ggplot2::geom_point(
    data = long,
    ggplot2::aes(
      x = running_variable_plot,
      y = duration_days,
      size = assessments_completed,
      colour = quantile,
      shape = quantile
    ),
    alpha = 0.52
  ) +
  ggplot2::geom_point(
    data = mean_points,
    ggplot2::aes(
      x = running_variable_plot,
      y = duration_days,
      size = assessments_completed,
      colour = "Reconstructed mean",
      shape = "Reconstructed mean"
    ),
    alpha = 0.82
  ) +
  ggplot2::geom_segment(
    data = mean_segments,
    ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE,
    linewidth = 1.15,
    colour = "#111111"
  ) +
  ggplot2::geom_text(
    data = mean_segments,
    ggplot2::aes(x = label_x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = -0.55,
    size = 3.1,
    colour = "#111111"
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq(-4, 4, by = 2),
    limits = c(-4.15, 4.15),
    expand = ggplot2::expansion(mult = 0)
  ) +
  ggplot2::scale_colour_manual(
    name = "Statistic",
    breaks = c(
      "25th percentile",
      "Median",
      "75th percentile",
      "Reconstructed mean"
    ),
    values = c(
      "25th percentile" = "#D55E00",
      "Median" = "#009E73",
      "75th percentile" = "#0072B2",
      "Reconstructed mean" = "#111111"
    )
  ) +
  
  ggplot2::scale_shape_manual(
    name = "Statistic",
    breaks = c(
      "25th percentile",
      "Median",
      "75th percentile",
      "Reconstructed mean"
    ),
    values = c(
      "25th percentile" = 16,   # circle
      "Median" = 17,            # triangle
      "75th percentile" = 15,   # square
      "Reconstructed mean" = 18 # diamond
    )
  ) +
  ggplot2::scale_size_area(
    max_size = 4.6,
    breaks = c(1000, 3000, 6000),
    labels = function(x) format(x, big.mark = ",", scientific = FALSE),
    name = "Assessments"
  ) +
  ggplot2::labs(
    title = "Disclosed quantiles and reconstructed mean around the Labour threshold",
    x = "Labour seat share minus 50 (percentage points)",
    y = "Assessment completion time (days)",
    caption = paste0(
      "Main +/-4 percentage-point window. Left: ", sum(control),
      " council-years and ", format(control_assessments, big.mark = ","),
      " assessments; right: ", sum(treated), " council-years and ",
      format(treated_assessments, big.mark = ","), " assessments.\n",
      "Orange circles, green triangles, and blue squares denote the disclosed ",
      "25th percentile, median, and 75th percentile, respectively. ",
      "Black diamonds denote qonlyreg reconstructed means; black segments are ",
      "their assessment-count-weighted side averages. ",
      "Horizontal jitter is for display only.\n",
      "The disclosed quantiles are council-year summaries, not individual-level ",
      "quantile discontinuities. Observations exactly at 50% are treated."
    )
  ) +
  ggplot2::theme_bw(base_size = 10.5) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold", hjust = 0),
    plot.caption = ggplot2::element_text(hjust = 0, size = 8),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.spacing.x = grid::unit(0.12, "cm")
  ) +
  ggplot2::guides(
    colour = ggplot2::guide_legend(
      order = 1,
      nrow = 1,
      title.position = "left",
      override.aes = list(
        shape = c(16, 17),
        size = 3,
        alpha = 1
      )
    ),
    shape = "none",
    size = ggplot2::guide_legend(
      order = 2,
      nrow = 1,
      title.position = "left",
      override.aes = list(
        shape = 16,
        colour = "#555555",
        alpha = 0.66
      )
    )
  )

p_combined <- p_combined +
  ggplot2::labs(title = NULL, caption = NULL)

three_panel_png <- file.path(output_dir, "rd_quantiles_main_window_4pp.png")
three_panel_pdf <- file.path(output_dir, "rd_quantiles_main_window_4pp.pdf")
combined_png <- file.path(
  output_dir,
  "rd_quantiles_and_reconstructed_mean_main_window_4pp.png"
)
combined_pdf <- file.path(
  output_dir,
  "rd_quantiles_and_reconstructed_mean_main_window_4pp.pdf"
)
ggplot2::ggsave(
  three_panel_png, p_three_panel,
  width = 10, height = 5.2, units = "in", dpi = 300
)
ggplot2::ggsave(
  three_panel_pdf, p_three_panel,
  width = 10, height = 5.2, units = "in"
)
ggplot2::ggsave(
  combined_png, p_combined,
  width = 9.2, height = 5.4, units = "in", dpi = 300
)
ggplot2::ggsave(
  combined_pdf, p_combined,
  width = 9.2, height = 5.4, units = "in"
)

message(
  "Two RD quantile figures written to output/ (N=", nrow(plot_data),
  ", councils=", length(unique(plot_data$new_la_code)), ")."
)
