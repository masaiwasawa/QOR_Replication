############################
## Quantile-Only Regression: Inference and Partial Identification under Quantile Disclosure
## NAEP empirical analysis
## Get effective df from NAEP NDE regression API
############################

library(httr)
library(jsonlite)
library(dplyr)
library(purrr)
library(tibble)

############################
## 1. State list
############################

state_codes <- c(
  "AL","AK","AZ","AR","CA","CO","CT","DE","DC",
  "FL","GA","HI","ID","IL","IN","IA","KS","KY",
  "LA","ME","MD","MA","MI","MN","MS","MO","MT",
  "NE","NV","NH","NJ","NM","NY","NC","ND","OH",
  "OK","OR","PA","PR","RI","SC","SD","TN","TX",
  "UT","VT","VA","WA","WV","WI","WY"
)

state_lookup <- c(
  AL = "Alabama",
  AK = "Alaska",
  AZ = "Arizona",
  AR = "Arkansas",
  CA = "California",
  CO = "Colorado",
  CT = "Connecticut",
  DE = "Delaware",
  DC = "District of Columbia",
  FL = "Florida",
  GA = "Georgia",
  HI = "Hawaii",
  ID = "Idaho",
  IL = "Illinois",
  IN = "Indiana",
  IA = "Iowa",
  KS = "Kansas",
  KY = "Kentucky",
  LA = "Louisiana",
  ME = "Maine",
  MD = "Maryland",
  MA = "Massachusetts",
  MI = "Michigan",
  MN = "Minnesota",
  MS = "Mississippi",
  MO = "Missouri",
  MT = "Montana",
  NE = "Nebraska",
  NV = "Nevada",
  NH = "New Hampshire",
  NJ = "New Jersey",
  NM = "New Mexico",
  NY = "New York",
  NC = "North Carolina",
  ND = "North Dakota",
  OH = "Ohio",
  OK = "Oklahoma",
  OR = "Oregon",
  PA = "Pennsylvania",
  PR = "Puerto Rico",
  RI = "Rhode Island",
  SC = "South Carolina",
  SD = "South Dakota",
  TN = "Tennessee",
  TX = "Texas",
  UT = "Utah",
  VT = "Vermont",
  VA = "Virginia",
  WA = "Washington",
  WV = "West Virginia",
  WI = "Wisconsin",
  WY = "Wyoming"
)
############################
## 2. Base payload
############################
payload_txt_base <- paste(
  readLines("payload_full.json", warn = FALSE),
  collapse = ""
)

make_payload_txt <- function(state_code) {
  state_label <- unname(state_lookup[state_code])
  
  gsub(
    '"selectedJurisdictions":\\[\\{"id":"AK","label":"Alaska","code":"AK"\\}\\]',
    paste0(
      '"selectedJurisdictions":[{"id":"',
      state_code,
      '","label":"',
      state_label,
      '","code":"',
      state_code,
      '"}]'
    ),
    payload_txt_base
  )
}


############################
## 3. API function
############################

get_state_df <- function(state_code) {
  
  Sys.sleep(0.5)
  
  payload_txt <- make_payload_txt(state_code)
  
  res <- httr::POST(
    url = "https://www.nationsreportcard.gov/ndecore/api/regAnalysis",
    httr::add_headers(
      `Content-Type` = "application/json;charset=utf-8",
      Accept = "application/json, text/plain, */*",
      Origin = "https://www.nationsreportcard.gov",
      Referer = "https://www.nationsreportcard.gov/ndecore/xplore/NDE"
    ),
    body = payload_txt,
    encode = "raw"
  )
  
  txt <- httr::content(res, as = "text", encoding = "UTF-8")
  
  if (httr::status_code(res) != 200) {
    stop(paste0("HTTP ", httr::status_code(res), ": ", txt))
  }
  
  json <- jsonlite::fromJSON(txt, simplifyVector = FALSE)
  
  error_df <- json$tabularDataList$`6`$tabDataRow$`2`$cellNumericData
  beta_eds <- json$tabularDataList$`11`$tabDataRow$`4`$cellNumericData
  se_eds   <- json$tabularDataList$`11`$tabDataRow$`5`$cellNumericData
  
  tibble::tibble(
    state = unname(state_lookup[state_code]),
    state_code = state_code,
    error_df = error_df,
    beta_eds = beta_eds,
    se_eds = se_eds,
    n_eff_total = error_df + 3  
    )
}

############################
## 4. Loop over states
############################
get_state_df_safe <- safely(get_state_df)

tmp <- vector("list", length(state_codes))

for (i in seq_along(state_codes)) {
  
  cat(
    "Running",
    i,
    "of",
    length(state_codes),
    ":",
    state_codes[i],
    "\n"
  )
  
  tmp[[i]] <- get_state_df_safe(state_codes[i])
}

############################
## 5. Successful results
############################

successful <- tmp %>%
  map("result") %>%
  compact() %>%
  bind_rows()

print(successful)


write.csv(
  successful,
  "naep_effective_sample_sizes.csv",
  row.names = FALSE
)
############################
## 6. Errors (optional)
############################

errors <- tibble(
  state = state_codes,
  error = map_chr(
    tmp,
    ~ {
      if (is.null(.x$error)) {
        ""
      } else {
        .x$error$message
      }
    }
  )
)
print(errors)

############################
## 7. make data
############################
############################
## 7A. Load NAEP percentile data
############################

library(tidyr)

naep_df <- read.csv(
  "NAEP_2024_math_grade8_eds.csv",
  check.names = FALSE,
  fileEncoding = "CP932",
  na.strings = c("#", "‡", "Ѓц")
)

naep_df <- naep_df %>%
  rename(
    state = Jurisdiction,
    eds_status = `Economically disadvantaged status`,
    mean_score = `Average scale score`,
    q10 = `10th percentile`,
    q25 = `25th percentile`,
    q50 = `50th percentile`,
    q75 = `75th percentile`,
    q90 = `90th percentile`
  ) %>%
  mutate(
    state = na_if(state, "")
  ) %>%
  fill(state, .direction = "down") %>%
  filter(
    state != "National",
    state != "DoDEA",
    eds_status %in% c(
      "Economically disadvantaged",
      "Not economically disadvantaged"
    )
  ) %>%
  mutate(
    eds = if_else(eds_status == "Economically disadvantaged", 1, 0)
  )

############################
## 7B. Load NAEP sample size share data
############################
share_df <- read.csv(
  "NAEP_2024_math_grade8_eds_share.csv",
  check.names = FALSE,
  na.strings = c("#", "‡")
)


share_df <- share_df %>%
  rename(
    state = Jurisdiction,
    share_eds = `Economically disadvantaged`,
    share_noneds = `Not economically disadvantaged`
  ) %>%
  mutate(
    share_eds = as.numeric(share_eds),
    share_noneds = as.numeric(share_noneds)
  ) %>%
  mutate(
    total_share = share_eds + share_noneds,
    pi_eds = share_eds / total_share,
    pi_noneds = share_noneds / total_share
  )


############################
## 7C. Merge effective n_g
############################

naep_df <- naep_df %>%
  left_join(successful, by = "state") %>%
  left_join(share_df, by = "state") %>%
  mutate(
    n_g = if_else(
      eds == 1,
      pi_eds * n_eff_total,
      pi_noneds * n_eff_total
    )
  )

############################
## 8. Save
############################
write.csv(
  naep_df,
  "naep_analysis_data.csv",
  row.names = FALSE
)
