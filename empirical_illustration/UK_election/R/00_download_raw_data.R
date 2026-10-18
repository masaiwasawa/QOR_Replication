# Download the two public source files used by the analysis and verify that
# their contents match the recorded SHA-256 digests.

project_root <- normalizePath("..", mustWork = TRUE)
raw_dir <- file.path(project_root, "data", "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

force_download <- FALSE

sources <- data.frame(
  file = c(
    "children_in_need_assessment_duration_2014_2025.csv",
    "council_composition_2016_2025.csv"
  ),
  source = c(
    "Department for Education: Characteristics of children in need",
    "electionresults.uk: annual council composition snapshot"
  ),
  version = c(
    "DfE data-set 3fb5a93d-1eb3-4f3a-bdd0-5caea4556ce6",
    "Git commit 29f587a2e761357c34cb828b0204861d8ec1c015"
  ),
  url = c(
    paste0(
      "https://explore-education-statistics.service.gov.uk/data-catalogue/",
      "data-set/3fb5a93d-1eb3-4f3a-bdd0-5caea4556ce6/csv"
    ),
    paste0(
      "https://raw.githubusercontent.com/fsargent/electionresults.uk/",
      "29f587a2e761357c34cb828b0204861d8ec1c015/",
      "source-data/council/history2016-2025.csv"
    )
  ),
  expected_sha256 = c(
    "275abab2383a7ca6da1223e49ca2facd589f7f45dd4b428f170c7fe3d85ca453",
    "18a99629b8b762cd7c71a2a3fff217f20a5490eed1e4fe0b2d5df9befeb88e51"
  ),
  stringsAsFactors = FALSE
)

sha256_file <- function(path) {
  if (requireNamespace("openssl", quietly = TRUE)) {
    con <- file(path, open = "rb")
    on.exit(close(con), add = TRUE)
    return(as.character(openssl::sha256(con)))
  }
  
  sha256sum <- Sys.which("sha256sum")
  if (nzchar(sha256sum)) {
    out <- system2(sha256sum, path, stdout = TRUE)
    return(strsplit(out[1], "[[:space:]]+")[[1]][1])
  }
  
  shasum <- Sys.which("shasum")
  if (nzchar(shasum)) {
    out <- system2(shasum, c("-a", "256", path), stdout = TRUE)
    return(strsplit(out[1], "[[:space:]]+")[[1]][1])
  }
  
  stop("Unable to compute SHA-256. R package 'openssl', sha256sum, or shasum is required.")
}

for (i in seq_len(nrow(sources))) {
  destination <- file.path(raw_dir, sources$file[i])
  
  if (file.exists(destination) && !force_download) {
    actual <- unname(tolower(sha256_file(destination)))
    if (!isTRUE(actual == sources$expected_sha256[i])) {
      stop(
        "The SHA-256 hash of the existing file does not match the recorded value: ", destination,
        "\nSet force_download <- TRUE if you intend to download the file again."
      )
    }
    message("Verified existing file: ", sources$file[i])
    next
  }
  
  temporary <- tempfile(tmpdir = raw_dir, fileext = ".download")
  tryCatch({
    message("Downloading: ", sources$url[i])
    download.file(sources$url[i], temporary, mode = "wb", method = "libcurl", quiet = TRUE)
    
    actual <- unname(tolower(sha256_file(temporary)))
    if (!isTRUE(actual == sources$expected_sha256[i])) {
      stop(
        "The SHA-256 hash of the downloaded file does not match the recorded value: ", sources$file[i],
        "\nexpected: ", sources$expected_sha256[i],
        "\nactual:   ", actual,
        "\nThe file at the public source may have been updated."
      )
    }
    
    if (file.exists(destination)) {
      unlink(destination)
    }
    if (!file.rename(temporary, destination)) {
      stop("Unable to move the downloaded file to its destination: ", destination)
    }
  }, finally = {
    if (file.exists(temporary)) unlink(temporary)
  })
  message("Downloaded and verified: ", sources$file[i])
}

sources$bytes <- file.info(file.path(raw_dir, sources$file))$size
sources$verified_sha256 <- unname(vapply(
  file.path(raw_dir, sources$file), function(path) unname(tolower(sha256_file(path))), character(1)
))
sources$verified_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)

write.csv(
  sources,
  file.path(raw_dir, "download_manifest.csv"),
  row.names = FALSE,
  na = ""
)

message("Raw-data manifest written to data/raw/download_manifest.csv")
