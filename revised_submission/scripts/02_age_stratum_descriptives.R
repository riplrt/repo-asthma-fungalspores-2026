#!/usr/bin/env Rscript
############################################################
## Age-Stratum Descriptives by Pandemic Phase
## Focus:  Encounters per day, ED/inpatient fraction, share of all encounters
## Supports the age-composition discussion requested in the Figure S1/S2 comment
##
## Outputs (CSV) written to: R2_table6_regenerated/tables_drafts/
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(lubridate)
})

out_dir <- "Journal_Aerobiology_Submission_response/R2_table6_regenerated/tables_drafts"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

san_juan_munis <- c("SAN JUAN", "GUAYNABO", "CAROLINA", "TRUJILLO ALTO", "CANOVANAS", "LOIZA")
clinical_settings <- c("ED", "Inpatient", "Outpatient", "Urgent")
age_breaks <- c(0, 4, 11, 17, 39, 64, 74, 84, Inf)
age_labels <- c("0-4", "5-11", "12-17", "18-39", "40-64", "65-74", "75-84", "85+")

care_setting_fun <- function(hcg_setting, hcg_mr_line_desc) {
  s <- tolower(as.character(hcg_setting))
  l <- tolower(as.character(hcg_mr_line_desc))
  if (grepl("ed|emergency", s) || grepl("\\bed\\b|\\ber\\b", l)) return("ED")
  if (grepl("inpatient|hospital", s) || grepl("inpatient|hospital", l)) return("Inpatient")
  if (grepl("urgent", s) || grepl("urgent", l)) return("Urgent")
  if (grepl("office|clinic|outpatient", s) || grepl("office|clinic", l)) return("Outpatient")
  "Other"
}

claims <- read_csv("data/2022_asma_From_date_ver 3.csv", show_col_types = FALSE) |>
  mutate(
    date = suppressWarnings(mdy(FROM_DATE)),
    municipality = toupper(trimws(MUNICIPALITY_DESC)),
    icd1 = toupper(trimws(ICD_DIAG_01_PRIMARY)),
    care_setting = mapply(care_setting_fun, HCG_SETTING, HCG_MR_LINE_DESC),
    age_num = suppressWarnings(as.numeric(gsub("\\+", "", as.character(Age)))),
    age_group = cut(age_num, breaks = age_breaks, labels = age_labels, right = TRUE, include.lowest = TRUE)
  ) |>
  filter(!is.na(date), year(date) %in% 2018:2021, municipality %in% san_juan_munis,
         grepl("^J45|^J46", icd1), care_setting %in% clinical_settings, !is.na(age_group)) |>
  mutate(phase = factor(case_when(year(date) < 2020 ~ "pre2020",
                                  year(date) == 2020 ~ "y2020",
                                  TRUE ~ "y2021"),
                        levels = c("pre2020", "y2020", "y2021")))

phase_days <- tibble(phase = factor(c("pre2020", "y2020", "y2021"), levels = c("pre2020", "y2020", "y2021")),
                     n_days = c(730, 366, 365))

desc <- claims |>
  group_by(phase, age_group) |>
  summarise(encounters = n(), ed_inpatient = sum(care_setting %in% c("ED", "Inpatient")), .groups = "drop") |>
  left_join(phase_days, by = "phase") |>
  group_by(phase) |>
  mutate(
    encounters_per_day = encounters / n_days,
    ed_inpatient_per_day = ed_inpatient / n_days,
    ed_inpatient_fraction = ed_inpatient / encounters,
    share_of_encounters = encounters / sum(encounters),
    share_of_ed_inpatient = ed_inpatient / sum(ed_inpatient)
  ) |>
  ungroup()

write_csv(desc, file.path(out_dir, "age_stratum_descriptives_by_phase.csv"))
