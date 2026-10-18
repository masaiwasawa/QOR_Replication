# Build the all-council-by-year Labour-seat-share RD file from the untouched
# public CSV files in data/raw.  The design follows the sample construction of
# Alonso and Andrews (2020): annual council composition, not election years
# only, and Labour seat share as the running variable.

project_root <- normalizePath("..", mustWork = TRUE)
raw_dir <- file.path(project_root, "data", "raw")
derived_dir <- file.path(project_root, "data", "derived")
dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- file.path(
  raw_dir,
  c(
    "children_in_need_assessment_duration_2014_2025.csv",
    "council_composition_2016_2025.csv"
  )
)
if (any(!file.exists(required_files))) {
  stop("Raw data files are missing. Run R/00_download_raw_data.R first.")
}

dfe <- read.csv(required_files[1], check.names = FALSE, stringsAsFactors = FALSE)
composition <- read.csv(required_files[2], check.names = FALSE, stringsAsFactors = FALSE)

required_dfe <- c(
  "time_period", "geographic_level", "new_la_code", "la_name",
  "Assessments_completed", "25th_Percentile_days", "Median_days",
  "75th_Percentile_days"
)
required_composition <- c(
  "year", "authority", "total", "con", "lab", "ld", "green", "ukip",
  "ref", "pc", "snp", "other", "majority"
)
stopifnot(
  all(required_dfe %in% names(dfe)),
  all(required_composition %in% names(composition))
)

authority_key <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("[’']", "", x)
  x <- gsub("&", " and ", x, fixed = TRUE)
  x <- gsub("[^a-z0-9]+", " ", x)
  x <- trimws(gsub(" +", " ", x))
  
  aliases <- c(
    "bristol city of" = "bristol",
    "kingston upon hull city of" = "kingston upon hull",
    "herefordshire county of" = "herefordshire",
    "county durham" = "durham",
    "buckinghamshire cc" = "buckinghamshire",
    "city of london corporation" = "city of london",
    "east riding" = "east riding of yorkshire"
  )
  replacement <- unname(aliases[x])
  x[!is.na(replacement)] <- replacement[!is.na(replacement)]
  x
}

numeric_or_na <- function(x) suppressWarnings(as.numeric(x))
within_band <- function(x, bandwidth_pp) abs(x) <= bandwidth_pp + 1e-12

# DfE time_period is the calendar year in which the reporting year ends.
# Example: time_period 2025 covers April 2024--March 2025.  The post-May-2024
# council composition snapshot is therefore joined to reporting year 2025.
dfe <- dfe[dfe$geographic_level == "Local authority", , drop = FALSE]
dfe$authority_key <- authority_key(dfe$la_name)
dfe$reporting_year <- as.integer(dfe$time_period)
dfe$council_composition_year <- dfe$reporting_year - 1L
dfe$fiscal_year <- paste0(
  dfe$council_composition_year,
  "/",
  substr(dfe$reporting_year, 3, 4)
)

dfe$assessments_completed_raw <- dfe$Assessments_completed
dfe$q25_days_raw <- dfe$`25th_Percentile_days`
dfe$median_days_raw <- dfe$Median_days
dfe$q75_days_raw <- dfe$`75th_Percentile_days`
dfe$assessments_completed <- numeric_or_na(dfe$Assessments_completed)
dfe$q25_days <- numeric_or_na(dfe$`25th_Percentile_days`)
dfe$median_days <- numeric_or_na(dfe$Median_days)
dfe$q75_days <- numeric_or_na(dfe$`75th_Percentile_days`)
dfe$quartiles_available <- complete.cases(dfe[c("q25_days", "median_days", "q75_days")])
dfe$assessment_count_available <- !is.na(dfe$assessments_completed) & dfe$assessments_completed > 0
dfe$quartiles_monotone <- with(
  dfe,
  is.na(q25_days) | is.na(median_days) | is.na(q75_days) |
    (q25_days <= median_days & median_days <= q75_days)
)
if (!all(dfe$quartiles_monotone)) {
  stop("The DfE data contain numeric rows that do not satisfy Q1 <= median <= Q3.")
}

composition$authority_key <- authority_key(composition$authority)
composition$year <- as.integer(composition$year)
composition$reporting_year <- composition$year + 1L
names(composition)[names(composition) == "total"] <- "council_total_seats"
names(composition)[names(composition) == "majority"] <- "source_control_label"

party_columns <- c("con", "lab", "ld", "green", "ukip", "ref", "pc", "snp", "other")
composition$party_seat_sum <- rowSums(composition[party_columns], na.rm = FALSE)
composition$seat_sum_matches_total <- composition$party_seat_sum == composition$council_total_seats
if (!all(composition$seat_sum_matches_total)) {
  stop("Some rows have party seat totals that do not match the reported total number of seats.")
}

# Percentage-point scale: h=5 means Labour seat share is in [45%, 55%].
composition$labour_seat_share <- with(composition, lab / council_total_seats)
composition$labour_seat_share_pct <- 100 * composition$labour_seat_share
composition$running_variable <- composition$labour_seat_share - 0.5
composition$running_variable_pp <- composition$labour_seat_share_pct - 50

# Alonso and Andrews (2020) code 50% as treated.  Because exactly 50% is not
# a strict majority in the ordinary institutional sense, both indicators are
# retained; labour_at_least_50 is the prespecified main definition.
composition$labour_at_least_50 <- composition$labour_seat_share >= 0.5 - 1e-12
composition$labour_strict_majority <- with(
  composition,
  lab > council_total_seats / 2
)
composition$labour_exactly_50 <- abs(composition$running_variable_pp) <= 1e-10
composition$labour_seats_relative_to_half <- with(
  composition,
  lab - council_total_seats / 2
)

# Only composition years that can be matched to an observed DfE reporting
# year are relevant.  The inner join automatically retains English upper-tier
# and single-tier authorities because those are the units in the DfE file.
composition_use <- composition[
  composition$reporting_year >= min(dfe$reporting_year) &
    composition$reporting_year <= max(dfe$reporting_year),
  ,
  drop = FALSE
]

analysis <- merge(
  composition_use,
  dfe,
  by = c("reporting_year", "authority_key"),
  all = FALSE,
  sort = FALSE
)

# The assessment duration has known lower support endpoint Q(0)=0.  It is
# included as a deterministic interpolation node.  Above Q3, where the upper
# support endpoint is unknown, the Section 3 constant-tail rule is retained.
analysis$known_support_lower_days <- 0
analysis$qonlyreg_mean_days <- with(
  analysis,
  0.125 * known_support_lower_days +
    0.250 * q25_days + 0.250 * median_days + 0.375 * q75_days
)
analysis$analysis_eligible <- with(
  analysis,
  quartiles_available & assessment_count_available
)
analysis$within_3pp <- within_band(analysis$running_variable_pp, 3.0)
analysis$within_4pp <- within_band(analysis$running_variable_pp, 4.0)
analysis$within_5pp <- within_band(analysis$running_variable_pp, 5.0)

keep <- c(
  "new_la_code", "la_name", "authority_key", "authority",
  "year", "reporting_year", "fiscal_year", "council_total_seats",
  party_columns, "party_seat_sum", "seat_sum_matches_total",
  "source_control_label", "labour_seat_share", "labour_seat_share_pct",
  "running_variable", "running_variable_pp", "labour_at_least_50",
  "labour_strict_majority", "labour_exactly_50",
  "labour_seats_relative_to_half", "assessments_completed_raw",
  "q25_days_raw", "median_days_raw", "q75_days_raw", "assessments_completed",
  "q25_days", "median_days", "q75_days", "quartiles_available",
  "assessment_count_available", "qonlyreg_mean_days", "analysis_eligible",
  "known_support_lower_days",
  "within_3pp", "within_4pp", "within_5pp"
)
analysis <- analysis[keep]
names(analysis)[names(analysis) == "year"] <- "council_composition_year"
analysis <- analysis[order(analysis$reporting_year, analysis$la_name), , drop = FALSE]

if (anyDuplicated(analysis[c("reporting_year", "new_la_code")])) {
  stop("The analysis data contain duplicate local authority-reporting year observations.")
}

analysis_path <- file.path(derived_dir, "children_in_need_labour_rd_all_years.csv")
write.csv(analysis, analysis_path, row.names = FALSE, na = "")

# Record DfE authority-years in the overlap period that could not be matched to
# an annual composition snapshot.  This makes boundary/reorganisation losses
# visible rather than silently dropping them.
dfe_overlap <- dfe[
  dfe$council_composition_year >= min(composition_use$year) &
    dfe$council_composition_year <= max(composition_use$year),
  ,
  drop = FALSE
]
matched_id <- paste(analysis$reporting_year, analysis$authority_key, sep = "::")
dfe_id <- paste(dfe_overlap$reporting_year, dfe_overlap$authority_key, sep = "::")
unmatched_dfe <- dfe_overlap[
  !(dfe_id %in% matched_id),
  c("new_la_code", "la_name", "authority_key", "reporting_year",
    "council_composition_year", "fiscal_year")
]
write.csv(
  unmatched_dfe,
  file.path(derived_dir, "unmatched_dfe_to_composition.csv"),
  row.names = FALSE,
  na = ""
)

bandwidths_pp <- c(3.0, 4.0, 5.0)
bandwidth_role <- function(bw) {
  if (bw == 4.0) return("main")
  if (bw < 4.0) return("narrower_window_sensitivity")
  "wider_window_sensitivity"
}
summary_rows <- lapply(bandwidths_pp, function(bw) {
  z <- analysis[
    analysis$analysis_eligible & within_band(analysis$running_variable_pp, bw),
    ,
    drop = FALSE
  ]
  data.frame(
    bandwidth_pp = bw,
    bandwidth_role = bandwidth_role(bw),
    observations = nrow(z),
    unique_councils = length(unique(z$new_la_code)),
    labour_at_least_50 = sum(z$labour_at_least_50),
    below_50 = sum(!z$labour_at_least_50),
    exactly_50 = sum(z$labour_exactly_50),
    composition_year_min = if (nrow(z)) min(z$council_composition_year) else NA,
    composition_year_max = if (nrow(z)) max(z$council_composition_year) else NA,
    reporting_year_min = if (nrow(z)) min(z$reporting_year) else NA,
    reporting_year_max = if (nrow(z)) max(z$reporting_year) else NA
  )
})
build_summary <- do.call(rbind, summary_rows)
write.csv(
  build_summary,
  file.path(derived_dir, "labour_rd_build_summary.csv"),
  row.names = FALSE,
  na = ""
)

message("All-year matched authority-years: ", nrow(analysis))
message("Eligible authority-years: ", sum(analysis$analysis_eligible))
message("Analysis data written to data/derived/children_in_need_labour_rd_all_years.csv")
