#!/usr/bin/env Rscript
############################################################
## Healthcare Utilization and Age-Denominator Context
## Focus:  Encounters by period, care setting and age stratum
## Descriptive tables for interpreting the ED/inpatient fraction (no models):
##   - Per-day encounters by period and setting; telehealth ascertainment
##   - Age-stratum descriptives and fraction change vs pre-2020
##   - Adult 40+ encounter share by period
##   - 2020 Census age context for the six municipalities (enrollment
##     denominators unavailable)
##
## Outputs (CSV) written to: R2_utilization_interpretation/tables_drafts/
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(lubridate)
  library(stringr)
  library(tidyr)
  library(tidycensus)
})

out_dir <- "Journal_Aerobiology_Submission_response/R2_utilization_interpretation/tables_drafts"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

muni6 <- c("SAN JUAN", "GUAYNABO", "CAROLINA", "TRUJILLO ALTO", "CANOVANAS", "LOIZA")
clinical_settings <- c("ED", "Inpatient", "Outpatient", "Urgent")

claims <- read_csv("data/2022_asma_From_date_ver 3.csv", show_col_types = FALSE) |>
  mutate(
    from_date = suppressWarnings(mdy(FROM_DATE)),
    year = year(from_date),
    muni = toupper(MUNICIPALITY_DESC),
    icd1 = toupper(trimws(ICD_DIAG_01_PRIMARY)),
    hcg_setting = as.character(HCG_SETTING),
    hcg_mr_line_desc = as.character(HCG_MR_LINE_DESC),
    line_group_desc = as.character(LINE_GROUP_DESC),
    form_type = as.character(FORM_TYPE),
    care_setting = case_when(
      grepl("ed|emergency", tolower(hcg_setting)) | grepl("\\bed\\b|\\ber\\b", tolower(hcg_mr_line_desc)) ~ "ED",
      grepl("inpatient|hospital", tolower(hcg_setting)) | grepl("inpatient|hospital", tolower(hcg_mr_line_desc)) ~ "Inpatient",
      grepl("urgent", tolower(hcg_setting)) | grepl("urgent", tolower(hcg_mr_line_desc)) ~ "Urgent",
      grepl("office|clinic|outpatient", tolower(hcg_setting)) | grepl("office|clinic", tolower(hcg_mr_line_desc)) ~ "Outpatient",
      grepl("dme|supplies", tolower(hcg_mr_line_desc)) ~ "DME",
      TRUE ~ "Other"
    ),
    asthma_primary = grepl("^J45|^J46", icd1),
    telehealth_flag = grepl("telehealth|telemedicine|virtual|audio|video|remote|phone",
                            paste(hcg_setting, hcg_mr_line_desc, line_group_desc, form_type),
                            ignore.case = TRUE),
    period = case_when(
      year < 2020 ~ "pre2020",
      year == 2020 ~ "y2020",
      year >= 2021 ~ "y2021plus"
    ),
    age_num = suppressWarnings(as.numeric(gsub("\\+", "", as.character(Age)))),
    age_group = cut(
      age_num,
      breaks = c(-Inf, 4, 11, 17, 39, 64, 74, 84, Inf),
      labels = c("0-4", "5-11", "12-17", "18-39", "40-64", "65-74", "75-84", "85+"),
      right = TRUE
    )
  ) |>
  filter(!is.na(from_date), muni %in% muni6)

asthma <- claims |>
  filter(asthma_primary, care_setting %in% clinical_settings, !is.na(age_group)) |>
  mutate(
    acute_flag = care_setting %in% c("ED", "Inpatient"),
    other_flag = care_setting %in% c("Outpatient", "Urgent")
  )

period_days <- asthma |>
  group_by(period) |>
  summarise(
    period_start = min(from_date, na.rm = TRUE),
    period_end = max(from_date, na.rm = TRUE),
    n_days = as.numeric(period_end - period_start) + 1,
    .groups = "drop"
  )

tab1 <- asthma |>
  group_by(period) |>
  summarise(
    ed_inpatient_n = sum(acute_flag),
    other_n = sum(other_flag),
    total_n = n(),
    telehealth_n = sum(telehealth_flag, na.rm = TRUE),
    telehealth_other_n = sum(telehealth_flag & other_flag, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(period_days, by = "period") |>
  mutate(
    ed_inpatient_per_day = ed_inpatient_n / n_days,
    other_per_day = other_n / n_days,
    total_per_day = total_n / n_days,
    ed_inpatient_encounter_fraction = ed_inpatient_n / (ed_inpatient_n + other_n),
    telehealth_share_total = telehealth_n / total_n,
    telehealth_share_other = telehealth_other_n / pmax(other_n, 1)
  )

write_csv(tab1, file.path(out_dir, "Table_R2_period_setting_per_day_draft.csv"))

age_desc <- asthma |>
  group_by(period, age_group) |>
  summarise(
    encounters = n(),
    ed_inpatient = sum(acute_flag),
    .groups = "drop"
  ) |>
  left_join(period_days, by = "period") |>
  group_by(period) |>
  mutate(
    encounters_per_day = encounters / n_days,
    ed_inpatient_per_day = ed_inpatient / n_days,
    ed_inpatient_fraction = ed_inpatient / encounters,
    share_of_encounters = encounters / sum(encounters)
  ) |>
  ungroup()

write_csv(age_desc, file.path(out_dir, "Table_R2_age_stratum_period_descriptives_draft.csv"))

age_change <- age_desc |>
  dplyr::select(period, age_group, ed_inpatient_fraction, share_of_encounters) |>
  pivot_wider(
    names_from = period,
    values_from = c(ed_inpatient_fraction, share_of_encounters)
  ) |>
  mutate(
    delta_fraction_y2020_vs_pre2020 = ed_inpatient_fraction_y2020 - ed_inpatient_fraction_pre2020,
    delta_fraction_y2021plus_vs_pre2020 = ed_inpatient_fraction_y2021plus - ed_inpatient_fraction_pre2020,
    delta_share_y2020_vs_pre2020 = share_of_encounters_y2020 - share_of_encounters_pre2020,
    delta_share_y2021plus_vs_pre2020 = share_of_encounters_y2021plus - share_of_encounters_pre2020
  )

write_csv(age_change, file.path(out_dir, "Table_R2_age_stratum_fraction_change_vs_pre2020_draft.csv"))

tele_detail <- asthma |>
  mutate(
    tele_source_hcg_setting = grepl("telehealth|telemedicine|virtual|audio|video|remote|phone", hcg_setting, ignore.case = TRUE),
    tele_source_line_desc = grepl("telehealth|telemedicine|virtual|audio|video|remote|phone", hcg_mr_line_desc, ignore.case = TRUE),
    tele_source_line_group = grepl("telehealth|telemedicine|virtual|audio|video|remote|phone", line_group_desc, ignore.case = TRUE),
    tele_source_form_type = grepl("telehealth|telemedicine|virtual|audio|video|remote|phone", form_type, ignore.case = TRUE)
  ) |>
  group_by(period) |>
  summarise(
    n_total = n(),
    tele_any = sum(telehealth_flag),
    tele_from_hcg_setting = sum(tele_source_hcg_setting),
    tele_from_line_desc = sum(tele_source_line_desc),
    tele_from_line_group = sum(tele_source_line_group),
    tele_from_form_type = sum(tele_source_form_type),
    .groups = "drop"
  )

write_csv(tele_detail, file.path(out_dir, "Table_R2_telehealth_ascertainment_draft.csv"))

adult40 <- age_desc |>
  mutate(adult40plus = age_group %in% c("40-64", "65-74", "75-84", "85+")) |>
  group_by(period) |>
  summarise(
    encounters_total = sum(encounters),
    encounters_age40plus = sum(encounters[adult40plus]),
    share_age40plus = encounters_age40plus / encounters_total,
    .groups = "drop"
  )

write_csv(adult40, file.path(out_dir, "Table_R2_adult40plus_encounter_share_by_period_draft.csv"))

api_key <- Sys.getenv("CENSUS_API_KEY")
if (api_key != "") {
  p12_vars <- c(sprintf("P12_%03dN", 1:49))

  c_age <- tryCatch(
    get_decennial(
      geography = "county",
      state = "72",
      variables = p12_vars,
      year = 2020,
      sumfile = "dhc",
      key = api_key
    ),
    error = function(e) e
  )

  if (!inherits(c_age, "error")) {
    c_age6 <- c_age |>
      mutate(municipio = str_remove(NAME, " Municipio, Puerto Rico")) |>
      filter(municipio %in% c("San Juan", "Guaynabo", "Carolina", "Trujillo Alto", "Canóvanas", "Loíza"))

    wide <- c_age6 |>
      dplyr::select(municipio, variable, value) |>
      pivot_wider(names_from = variable, values_from = value)

    census_broad <- wide |>
      summarise(
        total = sum(P12_001N),
        `0-17` = sum(P12_003N + P12_004N + P12_005N + P12_006N + P12_027N + P12_028N + P12_029N + P12_030N),
        `18-39` = sum(P12_007N + P12_008N + P12_009N + P12_010N + P12_011N + P12_012N + P12_013N +
                        P12_031N + P12_032N + P12_033N + P12_034N + P12_035N + P12_036N + P12_037N),
        `40-64` = sum(P12_014N + P12_015N + P12_016N + P12_017N + P12_018N + P12_019N +
                        P12_038N + P12_039N + P12_040N + P12_041N + P12_042N + P12_043N),
        `65+` = sum(P12_020N + P12_021N + P12_022N + P12_023N + P12_024N + P12_025N +
                      P12_044N + P12_045N + P12_046N + P12_047N + P12_048N + P12_049N)
      ) |>
      pivot_longer(cols = -total, names_to = "age_group_broad", values_to = "census_n") |>
      mutate(census_share = census_n / total)

    write_csv(census_broad, file.path(out_dir, "Table_R2_2020_census_age_context_six_municipalities_broad_draft.csv"))

    encounter_broad_2020 <- age_desc |>
      filter(period == "y2020") |>
      mutate(age_group_broad = case_when(
        age_group %in% c("0-4", "5-11", "12-17") ~ "0-17",
        age_group == "18-39" ~ "18-39",
        age_group == "40-64" ~ "40-64",
        age_group %in% c("65-74", "75-84", "85+") ~ "65+"
      )) |>
      group_by(age_group_broad) |>
      summarise(encounters = sum(encounters), .groups = "drop") |>
      mutate(encounter_share_2020 = encounters / sum(encounters))

    comp_broad <- census_broad |>
      left_join(encounter_broad_2020, by = "age_group_broad") |>
      mutate(diff_pct_points = (encounter_share_2020 - census_share) * 100)

    write_csv(comp_broad, file.path(out_dir, "Table_R2_2020_encounter_vs_census_age_context_broad_draft.csv"))
  }
}

message("Utilization interpretation drafts written to: ", out_dir)
