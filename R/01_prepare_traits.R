# R/01_prepare_traits.R
# -------------------------------------------------------
# Prepare and clean GROWELL trait data
# Outputs → results/01_traits/
# -------------------------------------------------------
source("R/00_load_packages.R")
cfg <- yaml::read_yaml("config/config.yml")
out <- cfg$output$s01
dir.create(out,               showWarnings = FALSE, recursive = TRUE)
dir.create(cfg$output$processed, showWarnings = FALSE, recursive = TRUE)

message("Loading trait data...")
traits_raw <- read_csv(cfg$data$traits, show_col_types = FALSE)
traits_raw$User_ID <- gsub("GROW-00+", "", traits_raw$User_ID)

# ── Sleep composite scores ────────────────────────────────
base_sleep <- c("covid_sleep_base","covid_sleep_less_base","covid_sleep_more_base",
                "covid_wakeup_base","covid_15_base","covid_7_base","covid_rested_base")
wk36_sleep <- c("covid_sleep_wk36","covid_sleep_less_wk36","covid_sleep_more_wk36",
                "covid_wakeup_wk36","covid_15_wk36","covid_7_wk36","covid_rested_wk36")
mth6_sleep <- c("covid_sleep_mth6","covid_sleep_less_mth6","covid_sleep_more_mth6",
                "covid_wakeup_mth6","covid_15_mth6","covid_7_mth6","covid_rested_mth6")
traits_raw[, base_sleep[1:5]]  <- !traits_raw[, base_sleep[1:5]]
traits_raw[, wk36_sleep[1:5]]  <- !traits_raw[, wk36_sleep[1:5]]
traits_raw[, mth6_sleep[1:5]]  <- !traits_raw[, mth6_sleep[1:5]]
traits_raw$sleep_base <- rowSums(traits_raw[, base_sleep])
traits_raw$sleep_wk36 <- rowSums(traits_raw[, wk36_sleep])
traits_raw$sleep_mth6 <- rowSums(traits_raw[, mth6_sleep])

traits_raw <- as.data.frame(traits_raw)
rownames(traits_raw) <- traits_raw$User_ID

# ── Logical → 0/1 factor ─────────────────────────────────
lv <- sapply(traits_raw, is.logical)
traits_raw[, lv] <- lapply(traits_raw[, lv], function(x)
  factor(x, levels = c(FALSE, TRUE), labels = c("0","1")))

# ── Character → ordered factor ────────────────────────────
traits_raw$group <- factor(traits_raw$group,
  levels = c("Control","Intervention"), labels = c(0,1))
traits_raw$race_eth_new <- factor(traits_raw$race_eth_new,
  levels = c("Multi-racial/Other","Hispanic","Asian","Black","White"),
  labels = c(0,1,2,3,4))
traits_raw$pregtype_conception_base <- factor(traits_raw$pregtype_conception_base,
  levels = c("Medical Help","Naturally"), labels = c(0,1))
traits_raw$pregtype_planned_base <- factor(traits_raw$pregtype_planned_base,
  levels = c("No","Yes"), labels = c(0,1))
traits_raw$moth_marit_stat_base <- factor(traits_raw$moth_marit_stat_base,
  levels = c("Divorced","Married","Separated","Single"), labels = c(0,1,2,3))
traits_raw$mother_s_education_base <- factor(traits_raw$mother_s_education_base,
  levels = c("High School","Post-Baccalaureate","Some College"), labels = c(0,1,2))
traits_raw$father_s_race_base <- factor(traits_raw$father_s_race_base,
  levels = c("Asian/South Pacific","Black or African American","Other","White"),
  labels = c(0,1,2,3))
traits_raw$father_s_hispanic_base <- factor(traits_raw$father_s_hispanic_base,
  levels = c("Hispanic or Latino","Not Hispanic or Latino"), labels = c(0,1))
traits_raw$father_s_education_base <- factor(traits_raw$father_s_education_base,
  levels = c("High School","Post-Baccalaureate","Some College"), labels = c(0,1,2))
traits_raw$second_hand_rules_base <- factor(traits_raw$second_hand_rules_base,
  levels = c("Allowed anywhere","Allowed in some places or at sometimes","Not allowed anywhere"),
  labels = c(0,1,2))
traits_raw$how_willing_are_you_to_mak_base <- factor(traits_raw$how_willing_are_you_to_mak_base,
  levels = c("Not at all willing","Somewhat willing","Uncertain","Very willing"),
  labels = c(0,1,2,3))
traits_raw$covid_test_yes_base <- factor(traits_raw$covid_test_yes_base,
  levels = c("No","Yes"), labels = c(0,1))
traits_raw$how_willing_are_you_to_mak_wk26 <- factor(traits_raw$how_willing_are_you_to_mak_wk26,
  levels = c("Somewhat willing","Uncertain","Very willing"), labels = c(0,1,2))
traits_raw$how_willing_are_you_to_mak_wk36 <- factor(traits_raw$how_willing_are_you_to_mak_wk36,
  levels = c("Not at all willing","Somewhat unwilling","Somewhat willing","Uncertain","Very willing"),
  labels = c(0,1,2,3,4))
traits_raw$infant_gender_mth3 <- factor(traits_raw$infant_gender_mth3,
  levels = c("Female","Male"), labels = c(0,1))
traits_raw$how_willing_are_you_to_mak_mth3 <- factor(traits_raw$how_willing_are_you_to_mak_mth3,
  levels = c("Somewhat unwilling","Somewhat willing","Uncertain","Very willing"),
  labels = c(0,1,2,3))
traits_raw$how_willing_are_you_to_mak_mth6 <- factor(traits_raw$how_willing_are_you_to_mak_mth6,
  levels = c("Somewhat unwilling","Somewhat willing","Uncertain","Very willing"),
  labels = c(0,1,2,3))
traits_raw$covid_test_yes_mth6 <- factor(traits_raw$covid_test_yes_mth6,
  levels = c("I don't know","No","Yes"), labels = c(NA,0,1))
traits_raw$covid_test_yes_wk36 <- factor(traits_raw$covid_test_yes_wk36,
  levels = c("I don't know","No","Yes"), labels = c(NA,0,1))

# ── Remove uninformative variables ───────────────────────
remove_vars <- c(
  "apo_diagnosis","enroll_race_v2_demo","enroll_race_v2_base",
  "enroll_hispanic_v2_base","enroll_hispanic_v2_demo","enroll_hispanic_sub_v2_demo",
  "enroll_gender_base","adh_base","number_of_therapeutic_abor_92de80_base",
  "postpart_visit_remote","number_of_ectopic_pregnanc_4b2a0e_base",
  "pregtype_assistance_base","num_spontaneous_abortions_base",
  "feeding_plan_base","moth_US_birth_base","mother_s_occu_base",
  "enroll_hispanic_sub_v2_base","father_s_occupation_base",
  "postpart_visit_mth3","feeding_other_type__1_mth3",
  "feeding_other_type__2_mth3","feeding_other_type__3_mth3",
  "postpart_visit_mth6","infant_gender_mth6",
  "feeding_other_type__1_mth6","feeding_other_type__2_mth6",
  "feeding_other_type__3_mth6","postpart_visit_bbmr",
  "infant_gender_bbmr","postpart_weight_bbmr",
  "infant_discharge_bbmr","infant_disposition_bbmr",
  "urbanicity","metro_area","father_s_origin_base",
  "adh_wk26","adh_wk36","adh_mth3","adh_mth6",
  "postpart_bp_mth6","postpart_bp_bbmr",
  "breastfeeding_amount","breastfeeding_continue",
  "breastfeeding_partial_age_bcheck",
  "do_any_of_your_coworkers_s_base",
  base_sleep, wk36_sleep, mth6_sleep
)
traits_clean <- traits_raw[, !(names(traits_raw) %in% remove_vars)]

# ── Derived variables ─────────────────────────────────────
traits_clean <- traits_clean %>%
  mutate(
    White_vs_nonwhite = case_when(race_eth_new == 4 ~ 1, TRUE ~ 0),
    apo_total = case_when(apo == 1 ~ 1, is.na(apo) ~ 0, TRUE ~ 0)
  )
traits_clean <- traits_clean[, !colnames(traits_clean) %in% c("apo","User_ID")]
colnames(traits_clean) <- gsub("^apo_total$", "apo", colnames(traits_clean))

# ── Convert all factors to numeric ───────────────────────
traits_clean[] <- lapply(traits_clean, function(x)
  if (is.factor(x)) as.numeric(as.character(x)) else x)
traits_clean <- traits_clean[order(rownames(traits_clean)), ]

# ── Demographic data with BMI by timepoint ───────────────
#group, age,White_vs_nonwhite, gest_age_at_birth,  parity_base,
demographic_data <- traits_clean %>%
  dplyr::select(weight_base, weight_wk26, weight_wk36, weight_mth3, weight_mth6,height_base,
                bmi, bmi_base, gwg, egwg, ppwr, ppwr_e, 
                apo, preterm, apo_hdp, apo_gdm, apo_other) %>%
  mutate(
    bmi_wk26 = weight_wk26 / (height_base^2),
    bmi_wk36 = weight_wk36 / (height_base^2),
    bmi_mth3 = weight_mth3 / (height_base^2),
    bmi_mth6 = weight_mth6 / (height_base^2)
  )

demographic_data <- demographic_data %>%
  dplyr::select(weight_base, weight_wk26, weight_wk36, weight_mth3, weight_mth6,height_base,
                bmi, bmi_base,bmi_wk26,bmi_wk36,bmi_mth3,bmi_mth6, gwg, egwg, ppwr, ppwr_e, 
                apo, preterm, apo_hdp, apo_gdm, apo_other)

# ── Save outputs ──────────────────────────────────────────
saveRDS(traits_clean,     file.path(cfg$output$processed, "traits_cleaned.rds"))
saveRDS(demographic_data, file.path(cfg$output$processed, "demographic_data_bmi.rds"))
write.xlsx(demographic_data,
           file.path(out, "demographic_data_bmi.xlsx"), rowNames = TRUE)
write.csv(traits_clean,
          file.path(out, "traits_cleaned.csv"), row.names = TRUE)

message("Script 01 complete → ", out)