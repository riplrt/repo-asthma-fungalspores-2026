############################################################
## Time-Stratified Case-Crossover Analysis
## Focus:  ED and Inpatient Asthma Events (2018 & 2020)
## Tests short-term trigger effects of cumulative spores (MA3/MA7)
## Controls for season/trend by construction
############################################################

# 0) Packages -------------------------------------------------------------

library(tidyverse)
library(readxl)
library(janitor)
library(lubridate)
library(zoo)
library(survival)    # For conditional logistic regression
library(broom)

# 1) SETTINGS -------------------------------------------------------------

# Data paths
claims_path_cases  <- "data/2022_asma_From_date_ver 3.xlsx"
claims_sheet_cases <- 1

spore_files_sj_cases <- c("data/Datos esporas y polen 24-08-19 a 31-12-22 Corregidas 12-09-25.xlsx")

# Years to analyze (especially 2018 and 2020)
years_to_run_cases <- c(2018, 2020)

# Weather data (pre-downloaded)
weather_path_cases <- "data/sju_weather_2018_2021.csv"

# Settings for ED/Inpatient definition
clinical_settings_cases <- c("ED", "Inpatient", "Outpatient", "Urgent")
acute_settings_cases    <- c("ED", "Inpatient")  # Define cases

# Moving average windows for cumulative spore exposure
ma_windows_cases <- c(3, 7)  # MA3 and MA7

# Case-crossover parameters
# Referent days:  same month, same day-of-week, bidirectional
max_referent_weeks_cases <- 3  # Look ±3 weeks for referents

# 2) HELPER FUNCTIONS -----------------------------------------------------

parse_date_any_cases <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30"))
  x <- as.character(x)
  d <- suppressWarnings(ymd(x))
  d <- coalesce(d, suppressWarnings(mdy(x)))
  d <- coalesce(d, suppressWarnings(dmy(x)))
  coalesce(d, suppressWarnings(as.Date(x)))
}

care_setting_fun_cases <- function(hcg_setting, hcg_mr_line_desc) {
  s <- tolower(as.character(hcg_setting))
  l <- tolower(as.character(hcg_mr_line_desc))
  
  if (grepl("ed|emergency", s) || grepl("\\bed\\b|er", l)) return("ED")
  if (grepl("inpatient|hospital", s) || grepl("inpatient|hospital", l)) return("Inpatient")
  if (grepl("urgent", s) || grepl("urgent", l)) return("Urgent")
  if (grepl("office|clinic|outpatient", s) || grepl("office|clinic", l)) return("Outpatient")
  if (grepl("dme|supplies", l)) return("DME")
  "Other"
}

ensure_care_setting_cases <- function(dat_cases) {
  if ("care_setting_cases" %in% names(dat_cases)) return(dat_cases)
  
  if ("care_setting" %in% names(dat_cases)) {
    return(dat_cases %>% mutate(care_setting_cases = care_setting))
  }
  
  needed_cases <- c("hcg_setting", "hcg_mr_line_desc")
  if (! all(needed_cases %in% names(dat_cases))) {
    stop("No care_setting and missing hcg_setting/hcg_mr_line_desc to create it.")
  }
  
  dat_cases %>%
    mutate(care_setting_cases = mapply(care_setting_fun_cases, hcg_setting, hcg_mr_line_desc))
}

ensure_asthma_primary_cases <- function(dat_cases) {
  if ("asthma_primary_cases" %in% names(dat_cases)) return(dat_cases)
  
  if ("asthma_primary" %in% names(dat_cases)) {
    return(dat_cases %>% mutate(asthma_primary_cases = asthma_primary))
  }
  
  if (!"icd_diag_01_primary" %in% names(dat_cases)) {
    stop("Missing asthma_primary and icd_diag_01_primary; create asthma_primary before modeling.")
  }
  
  dat_cases %>%
    mutate(
      icd_diag_01_primary = str_to_upper(str_trim(as.character(icd_diag_01_primary))),
      asthma_primary_cases = if_else(str_detect(icd_diag_01_primary, "^J45|^J46"), "asthma", "non-asthma")
    )
}

read_claims_cases <- function(path_cases, sheet_cases = 1) {
  dat_cases <- read_excel(path_cases, sheet = sheet_cases) %>%
    clean_names()
  
  if (!"from_date" %in% names(dat_cases)) stop("Missing from_date (after clean_names()).")
  
  dat_cases %>%
    mutate(
      from_date = parse_date_any_cases(from_date),
      year_cases = lubridate::year(from_date)
    )
}

read_spores_sj_cases <- function(file_cases) {
  x_cases <- read_excel(file_cases) %>% clean_names()
  
  date_col_cases   <- names(x_cases)[1]
  spore_col_cases  <- names(x_cases)[str_detect(names(x_cases), "esporas|spores")]
  pollen_col_cases <- names(x_cases)[str_detect(names(x_cases), "polen|pollen")]
  
  if (length(spore_col_cases) == 0) stop("Could not find spores column in:  ", file_cases)
  
  x_cases %>%
    transmute(
      date = parse_date_any_cases(.data[[date_col_cases]]),
      spores_total_cases = suppressWarnings(as.numeric(.data[[spore_col_cases[1]]])),
      pollen_total_cases = if (length(pollen_col_cases) > 0) suppressWarnings(as.numeric(.data[[pollen_col_cases[1]]])) else NA_real_
    ) %>%
    filter(!is.na(date))
}

# 3) LOAD & PREPARE DATA --------------------------------------------------

cat("\n=== Loading Claims Data ===\n")
claims_raw_cases <- read_claims_cases(claims_path_cases, claims_sheet_cases) %>%
  ensure_care_setting_cases() %>%
  ensure_asthma_primary_cases()

# Filter to San Juan region, asthma primary diagnosis
claims_df_cases <- claims_raw_cases %>%
  mutate(
    from_date = as.Date(from_date),
    year_cases = year(from_date),
    health_region_cases = if_else(municipality_desc == "Caguas", "Caguas", "San Juan"),
    health_region_cases = factor(health_region_cases, levels = c("San Juan", "Caguas"))
  )

# Define ED/Inpatient cases
claims_cases <- claims_df_cases %>%
  filter(
    health_region_cases == "San Juan",
    asthma_primary_cases == "asthma",
    care_setting_cases %in% acute_settings_cases,
    year_cases %in% years_to_run_cases
  ) %>%
  mutate(
    case_id_cases = row_number(),  # Unique case identifier
    month_cases = month(from_date),
    dow_cases = wday(from_date),    # Day of week (1=Sunday)
    ym_cases = floor_date(from_date, "month")
  )

cat(sprintf("Total ED/Inpatient cases identified: %d\n", nrow(claims_cases)))
cat("\nCases by year:\n")
print(table(claims_cases$year_cases))
cat("\nCases by care setting:\n")
print(table(claims_cases$care_setting_cases))

# 4) LOAD ENVIRONMENTAL DATA ----------------------------------------------

cat("\n=== Loading Spores Data ===\n")
spores_raw_cases <- map_dfr(spore_files_sj_cases, read_spores_sj_cases) %>%
  distinct(date, .keep_all = TRUE) %>%
  arrange(date) %>%
  filter(year(date) %in% years_to_run_cases)

cat(sprintf("Spore measurements available:  %d days\n", sum(! is.na(spores_raw_cases$spores_total_cases))))

# Calculate moving averages (MA3, MA7)
spores_ma_cases <- spores_raw_cases %>%
  arrange(date) %>%
  mutate(
    spores_ma3_cases = rollmean(spores_total_cases, k = 3, fill = NA, align = "right"),
    spores_ma7_cases = rollmean(spores_total_cases, k = 7, fill = NA, align = "right")
  )

cat("\n=== Loading Weather Data ===\n")
weather_cases <- read_csv(weather_path_cases, show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(date = as.Date(date)) %>%
  filter(year(date) %in% years_to_run_cases)

# Calculate weather moving averages
weather_ma_cases <- weather_cases %>%
  arrange(date) %>%
  mutate(
    tmax_ma3_cases = rollmean(tmax_c, k = 3, fill = NA, align = "right"),
    tmax_ma7_cases = rollmean(tmax_c, k = 7, fill = NA, align = "right"),
    rh_ma3_cases = rollmean(rh_mean, k = 3, fill = NA, align = "right"),
    rh_ma7_cases = rollmean(rh_mean, k = 7, fill = NA, align = "right"),
    prcp_ma3_cases = rollmean(prcp_mm, k = 3, fill = NA, align = "right"),
    prcp_ma7_cases = rollmean(prcp_mm, k = 7, fill = NA, align = "right")
  )

# 5) CREATE CASE-CROSSOVER DATASET ----------------------------------------

cat("\n=== Creating Case-Crossover Dataset ===\n")

# Function to find referent days for each case
find_referents_cases <- function(case_date, case_month, case_dow, max_weeks = 3) {
  # Bidirectional:  look before and after
  possible_dates <- seq(case_date - weeks(max_weeks), 
                        case_date + weeks(max_weeks), 
                        by = "1 day")
  
  # Filter to same month and same day of week, exclude the case day itself
  referents <- possible_dates[
    month(possible_dates) == case_month &
    wday(possible_dates) == case_dow &
    possible_dates != case_date
  ]
  
  return(referents)
}

# Create case-crossover dataset
crossover_list_cases <- list()

for (i in 1:nrow(claims_cases)) {
  case_row <- claims_cases[i, ]
  
  # Find referent days
  referent_dates <- find_referents_cases(
    case_row$from_date,
    case_row$month_cases,
    case_row$dow_cases,
    max_weeks = max_referent_weeks_cases
  )
  
  # Case observation
  case_obs <- tibble(
    case_id_cases = case_row$case_id_cases,
    date = case_row$from_date,
    is_case_cases = 1,
    year_cases = case_row$year_cases,
    stratum_id_cases = i  # Each case forms its own stratum
  )
  
  # Referent observations
  if (length(referent_dates) > 0) {
    ref_obs <- tibble(
      case_id_cases = case_row$case_id_cases,
      date = referent_dates,
      is_case_cases = 0,
      year_cases = year(referent_dates),
      stratum_id_cases = i
    )
    
    crossover_list_cases[[i]] <- bind_rows(case_obs, ref_obs)
  } else {
    crossover_list_cases[[i]] <- case_obs
  }
}

crossover_data_cases <- bind_rows(crossover_list_cases)

cat(sprintf("Total case-crossover observations: %d\n", nrow(crossover_data_cases)))
cat(sprintf("Unique strata (cases): %d\n", n_distinct(crossover_data_cases$stratum_id_cases)))

# 6) MERGE EXPOSURES ------------------------------------------------------

cat("\n=== Merging Exposure Data ===\n")

# Merge spores and weather
crossover_full_cases <- crossover_data_cases %>%
  left_join(
    spores_ma_cases %>% 
      select(date, spores_total_cases, spores_ma3_cases, spores_ma7_cases),
    by = "date"
  ) %>%
  left_join(
    weather_ma_cases %>% 
      select(date, tmax_c, rh_mean, prcp_mm, wind_ms, dewpoint_c,
             tmax_ma3_cases, tmax_ma7_cases, 
             rh_ma3_cases, rh_ma7_cases,
             prcp_ma3_cases, prcp_ma7_cases),
    by = "date"
  )

# Check completeness
complete_cases_ma3 <- crossover_full_cases %>%
  filter(! is.na(spores_ma3_cases) & !is.na(tmax_ma3_cases) & !is.na(rh_ma3_cases))

complete_cases_ma7 <- crossover_full_cases %>%
  filter(!is.na(spores_ma7_cases) & !is.na(tmax_ma7_cases) & !is.na(rh_ma7_cases))

cat(sprintf("Complete observations for MA3 analysis: %d\n", nrow(complete_cases_ma3)))
cat(sprintf("Complete observations for MA7 analysis: %d\n", nrow(complete_cases_ma7)))

# 7) CONDITIONAL LOGISTIC REGRESSION --------------------------------------

cat("\n=== Running Conditional Logistic Regression Models ===\n")

# Helper function to run clogit and extract results
run_clogit_cases <- function(data_cases, exposure_var, weather_vars, lag_label) {
  
  # Build formula
  formula_str <- paste0(
    "is_case_cases ~ ", exposure_var, " + ",
    paste(weather_vars, collapse = " + "),
    " + strata(stratum_id_cases)"
  )
  
  # Fit model
  model_cases <- tryCatch(
    clogit(as.formula(formula_str), data = data_cases),
    error = function(e) {
      message(sprintf("Error fitting %s model: %s", lag_label, e$message))
      return(NULL)
    }
  )
  
  if (is.null(model_cases)) {
    return(tibble(
      lag = lag_label,
      term = exposure_var,
      or_cases = NA_real_,
      ci_low_cases = NA_real_,
      ci_high_cases = NA_real_,
      p_value_cases = NA_real_,
      n_strata = NA_integer_,
      n_obs = NA_integer_
    ))
  }
  
  # Extract results
  results_tidy <- tidy(model_cases, exponentiate = TRUE, conf.int = TRUE)
  
  exposure_result <- results_tidy %>%
    filter(term == exposure_var) %>%
    transmute(
      lag = lag_label,
      term = term,
      or_cases = estimate,
      ci_low_cases = conf.low,
      ci_high_cases = conf.high,
      p_value_cases = p.value,
      n_strata = model_cases$n,
      n_obs = nrow(data_cases)
    )
  
  return(exposure_result)
}

# Initialize results list
results_list_cases <- list()

# For each year
for (yr in years_to_run_cases) {
  
  cat(sprintf("\n--- Year %d ---\n", yr))
  
  # MA3 models
  data_ma3_yr <- complete_cases_ma3 %>% filter(year_cases == yr)
  
  cat(sprintf("MA3 - Strata: %d, Observations:  %d\n", 
              n_distinct(data_ma3_yr$stratum_id_cases), 
              nrow(data_ma3_yr)))
  
  if (nrow(data_ma3_yr) > 0) {
    result_ma3 <- run_clogit_cases(
      data_ma3_yr,
      exposure_var = "spores_ma3_cases",
      weather_vars = c("tmax_ma3_cases", "rh_ma3_cases", "prcp_ma3_cases"),
      lag_label = sprintf("%d_MA3", yr)
    )
    results_list_cases[[length(results_list_cases) + 1]] <- result_ma3
  }
  
  # MA7 models
  data_ma7_yr <- complete_cases_ma7 %>% filter(year_cases == yr)
  
  cat(sprintf("MA7 - Strata: %d, Observations: %d\n", 
              n_distinct(data_ma7_yr$stratum_id_cases), 
              nrow(data_ma7_yr)))
  
  if (nrow(data_ma7_yr) > 0) {
    result_ma7 <- run_clogit_cases(
      data_ma7_yr,
      exposure_var = "spores_ma7_cases",
      weather_vars = c("tmax_ma7_cases", "rh_ma7_cases", "prcp_ma7_cases"),
      lag_label = sprintf("%d_MA7", yr)
    )
    results_list_cases[[length(results_list_cases) + 1]] <- result_ma7
  }
}

# Combine results
results_table_cases <- bind_rows(results_list_cases)

# 8) DISPLAY & SAVE RESULTS -----------------------------------------------

cat("\n=== RESULTS ===\n\n")
print(results_table_cases, n = Inf)

# Create output directory
dir.create("statistical_results", showWarnings = FALSE, recursive = TRUE)

# Save results
output_file_cases <- "statistical_results/case_crossover_ed_inpatient_results_cases.csv"
write_csv(results_table_cases, output_file_cases)
cat(sprintf("\nResults saved to:  %s\n", output_file_cases))

# Additional summary by exposure window
cat("\n=== Summary by Exposure Window ===\n")
results_summary_cases <- results_table_cases %>%
  separate(lag, into = c("year", "window"), sep = "_") %>%
  group_by(window) %>%
  summarise(
    n_analyses = n(),
    mean_or_cases = mean(or_cases, na.rm = TRUE),
    median_or_cases = median(or_cases, na.rm = TRUE),
    n_significant = sum(p_value_cases < 0.05, na.rm = TRUE),    .groups = "drop"
  )

print(results_summary_cases)

# 9) DIAGNOSTIC PLOTS -----------------------------------------------------

cat("\n=== Generating Diagnostic Plots ===\n")

# Plot:  ORs with confidence intervals
p_forest_cases <- results_table_cases %>%
  filter(! is.na(or_cases)) %>%
  ggplot(aes(x = lag, y = or_cases)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_low_cases, ymax = ci_high_cases), width = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(
    title = "Case-Crossover Analysis:  Spore Exposure and ED/Inpatient Asthma Events",
    subtitle = "Time-stratified design with same month, same day-of-week matching",
    x = "Year and Exposure Window",
    y = "Odds Ratio (95% CI)",
    caption = "Adjusted for temperature, humidity, and precipitation moving averages"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

ggsave("exploratory_results/case_crossover_forest_plot_cases.png", 
       p_forest_cases, width = 10, height = 6, dpi = 300)

cat("\nForest plot saved to: exploratory_results/case_crossover_forest_plot_cases.png\n")

# Summary statistics
cat("\n=== Data Summary ===\n")
cat(sprintf("Total ED/Inpatient cases analyzed: %d\n", 
            n_distinct(crossover_full_cases$case_id_cases)))
cat(sprintf("Date range: %s to %s\n", 
            min(crossover_full_cases$date), 
            max(crossover_full_cases$date)))

crossover_full_cases %>%
  group_by(year_cases, is_case_cases) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = is_case_cases, values_from = n, 
              names_prefix = "is_case_") %>%
  rename(year = year_cases, cases = is_case_1, referents = is_case_0) %>%
  print()

cat("\n=== Analysis Complete ===\n")