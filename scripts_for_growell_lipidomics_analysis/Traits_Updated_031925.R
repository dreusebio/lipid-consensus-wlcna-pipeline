#This is the reanalysis of the lipidomics. Trait file used is new because Sebastian had changed the variables. 
#So all the analysis need to be repeated..
#Analysis was also repeated because we added normalizations by median, base 10 log transformation and pareto scaling
#load data
#setwd to GROWELL_folder
library(caret)
library(randomForest)
library(xgboost)
library(ggplot2)
library(readr)
library(dplyr)
library(WGCNA)
library(writexl)
library(reshape2)
library(pheatmap)
library(viridis)
library(ComplexHeatmap)
library(circlize)

#=====================================================================================
#
#  Code chunk 1a- Load your trait file. Convert variables to numeric for WGCNA prep
#
#=====================================================================================
# Load your data
ID_Updated_Combined_Traits_apo <- read_csv("Identified_lipids_used_in_metaboanalyst_to_normalize/traits_taken_box_03_19_25/growell_curated_deid_george.csv")
dim(ID_Updated_Combined_Traits_apo)
#decided to remove prefix for GROW so that I can see numbers better
ID_Updated_Combined_Traits_apo$User_ID <- gsub("GROW-00+", "", ID_Updated_Combined_Traits_apo$User_ID)


colnames(ID_Updated_Combined_Traits_apo)
table(ID_Updated_Combined_Traits_apo$apo)
str(ID_Updated_Combined_Traits_apo$apo)

base_sleep <- c("covid_sleep_base", "covid_sleep_less_base", "covid_sleep_more_base","covid_wakeup_base",                   
                "covid_15_base",  "covid_7_base",   "covid_rested_base")
wk_36_sleep <- c("covid_sleep_wk36", "covid_sleep_less_wk36", "covid_sleep_more_wk36", 
                 "covid_wakeup_wk36", "covid_15_wk36", "covid_7_wk36", "covid_rested_wk36" )
mth6_sleep <- c("covid_sleep_mth6", "covid_sleep_less_mth6", "covid_sleep_more_mth6",                 
                "covid_wakeup_mth6", "covid_15_mth6" ,"covid_7_mth6", "covid_rested_mth6")

# Reverse scoring
ID_Updated_Combined_Traits_apo[, base_sleep[1:5]] <- !ID_Updated_Combined_Traits_apo[, base_sleep[1:5]]
ID_Updated_Combined_Traits_apo[, wk_36_sleep[1:5]] <- !ID_Updated_Combined_Traits_apo[, wk_36_sleep[1:5]]
ID_Updated_Combined_Traits_apo[, mth6_sleep[1:5]] <- !ID_Updated_Combined_Traits_apo[, mth6_sleep[1:5]]
# Compute total score
ID_Updated_Combined_Traits_apo$sleep_base <- rowSums(ID_Updated_Combined_Traits_apo[, base_sleep])
ID_Updated_Combined_Traits_apo$sleep_wk36 <- rowSums(ID_Updated_Combined_Traits_apo[, wk_36_sleep])
ID_Updated_Combined_Traits_apo$sleep_mth6 <- rowSums(ID_Updated_Combined_Traits_apo[, mth6_sleep])

# Set rownames
ID_Updated_Combined_Traits_apo <- as.data.frame(ID_Updated_Combined_Traits_apo)
rownames(ID_Updated_Combined_Traits_apo) <- ID_Updated_Combined_Traits_apo$User_ID

# Convert logical variables to factors with levels 0 for FALSE and 1 for TRUE, and labels "FALSE" and "TRUE"
logical_vars <- sapply(ID_Updated_Combined_Traits_apo, is.logical)
ID_Updated_Combined_Traits_apo[, logical_vars] <- lapply(ID_Updated_Combined_Traits_apo[, logical_vars], function(x) factor(x, levels = c(FALSE, TRUE), labels = c("0", "1")))


# Convert specific character variables to factors and set levels
ID_Updated_Combined_Traits_apo$group <- factor(ID_Updated_Combined_Traits_apo$group, levels = c("Control", "Intervention"), labels = c(0, 1))
ID_Updated_Combined_Traits_apo$race_eth_new <- factor(ID_Updated_Combined_Traits_apo$race_eth_new, levels = c("Multi-racial/Other", "Hispanic", "Asian",  "Black", "White"), labels = c(0, 1, 2, 3, 4))
#ID_Updated_Combined_Traits_apo$enroll_hispanic_v2_base <- factor(ID_Updated_Combined_Traits_apo$enroll_hispanic_v2_base, levels = c("Hispanic or Latino", "Not Hispanic or Latino"), labels = c(0, 1))
ID_Updated_Combined_Traits_apo$pregtype_conception_base <- factor(ID_Updated_Combined_Traits_apo$pregtype_conception_base, levels = c("Medical Help", "Naturally"), labels = c(0, 1))
ID_Updated_Combined_Traits_apo$pregtype_planned_base <- factor(ID_Updated_Combined_Traits_apo$pregtype_planned_base, levels = c("No", "Yes"), labels = c(0, 1))
ID_Updated_Combined_Traits_apo$moth_marit_stat_base <- factor(ID_Updated_Combined_Traits_apo$moth_marit_stat_base, levels = c("Divorced", "Married", "Separated", "Single"), labels = c(0, 1, 2, 3))
ID_Updated_Combined_Traits_apo$mother_s_education_base <- factor(ID_Updated_Combined_Traits_apo$mother_s_education_base, levels = c("High School", "Post-Baccalaureate", "Some College"), labels = c(0, 1, 2))
ID_Updated_Combined_Traits_apo$father_s_race_base <- factor(ID_Updated_Combined_Traits_apo$father_s_race_base, levels = c("Asian/South Pacific", "Black or African American", "Other", "White"), labels = c(0, 1, 2, 3))
ID_Updated_Combined_Traits_apo$father_s_hispanic_base <- factor(ID_Updated_Combined_Traits_apo$father_s_hispanic_base, levels = c("Hispanic or Latino", "Not Hispanic or Latino"), labels = c(0, 1))
ID_Updated_Combined_Traits_apo$father_s_education_base <- factor(ID_Updated_Combined_Traits_apo$father_s_education_base, levels = c("High School", "Post-Baccalaureate", "Some College"), labels = c(0, 1, 2))
ID_Updated_Combined_Traits_apo$second_hand_rules_base <- factor(ID_Updated_Combined_Traits_apo$second_hand_rules_base, levels = c("Allowed anywhere", "Allowed in some places or at sometimes", "Not allowed anywhere"), labels = c(0, 1, 2))
ID_Updated_Combined_Traits_apo$how_willing_are_you_to_mak_base <- factor(ID_Updated_Combined_Traits_apo$how_willing_are_you_to_mak_base, levels = c("Not at all willing", "Somewhat willing", "Uncertain", "Very willing"), labels = c(0, 1, 2, 3))
ID_Updated_Combined_Traits_apo$covid_test_yes_base <- factor(ID_Updated_Combined_Traits_apo$covid_test_yes_base, levels = c("No", "Yes"), labels = c(0, 1))
ID_Updated_Combined_Traits_apo$how_willing_are_you_to_mak_wk26 <- factor(ID_Updated_Combined_Traits_apo$how_willing_are_you_to_mak_wk26, levels = c("Somewhat willing", "Uncertain", "Very willing"), labels = c(0, 1, 2))
ID_Updated_Combined_Traits_apo$how_willing_are_you_to_mak_wk36 <- factor(ID_Updated_Combined_Traits_apo$how_willing_are_you_to_mak_wk36, levels = c("Not at all willing", "Somewhat unwilling", "Somewhat willing", "Uncertain", "Very willing"), labels = c(0, 1, 2, 3, 4))
ID_Updated_Combined_Traits_apo$infant_gender_mth3 <- factor(ID_Updated_Combined_Traits_apo$infant_gender_mth3, levels = c("Female", "Male"), labels = c(0, 1))
ID_Updated_Combined_Traits_apo$how_willing_are_you_to_mak_mth3 <- factor(ID_Updated_Combined_Traits_apo$how_willing_are_you_to_mak_mth3, levels = c("Somewhat unwilling", "Somewhat willing", "Uncertain", "Very willing"), labels = c(0, 1, 2, 3))
ID_Updated_Combined_Traits_apo$how_willing_are_you_to_mak_mth6 <- factor(ID_Updated_Combined_Traits_apo$how_willing_are_you_to_mak_mth6, levels = c("Somewhat unwilling", "Somewhat willing", "Uncertain", "Very willing"), labels = c(0, 1, 2, 3))
ID_Updated_Combined_Traits_apo$covid_test_yes_mth6 <- factor(ID_Updated_Combined_Traits_apo$covid_test_yes_mth6, levels = c("I don't know", "No", "Yes"), labels = c(NA, 0, 1))
ID_Updated_Combined_Traits_apo$covid_test_yes_wk36 <- factor(ID_Updated_Combined_Traits_apo$covid_test_yes_wk36, levels = c("I don't know", "No", "Yes"), labels = c(NA, 0, 1))

colnames(ID_Updated_Combined_Traits_apo)
# Remove specified variables 
remove_vars <- c( "apo_diagnosis", "enroll_race_v2_demo", "enroll_race_v2_base", "enroll_hispanic_v2_base", "enroll_hispanic_v2_demo", "enroll_hispanic_sub_v2_demo", 
                 "enroll_gender_base", "adh_base", "number_of_therapeutic_abor_92de80_base", "postpart_visit_remote", 
                 "number_of_ectopic_pregnanc_4b2a0e_base", "pregtype_assistance_base", "num_spontaneous_abortions_base", 
                 "feeding_plan_base", "moth_US_birth_base", "mother_s_occu_base", "enroll_hispanic_sub_v2_base", 
                 "father_s_occupation_base", "postpart_visit_mth3", "feeding_other_type__1_mth3", 
                 "feeding_other_type__2_mth3", "feeding_other_type__3_mth3", "postpart_visit_mth6", 
                 "infant_gender_mth6", "feeding_other_type__1_mth6", "feeding_other_type__2_mth6", 
                 "feeding_other_type__3_mth6", "postpart_visit_bbmr", "infant_gender_bbmr", 
                 "postpart_weight_bbmr", "infant_discharge_bbmr", "infant_disposition_bbmr",
                 "urbanicity", "metro_area", "father_s_origin_base", "adh_wk26", "adh_wk36", "adh_mth3", "adh_mth6",
                 "postpart_bp_mth6",     "postpart_bp_bbmr",     "breastfeeding_amount","breastfeeding_continue",
                 "breastfeeding_partial_age_bcheck", "covid_sleep_base", "covid_sleep_less_base", "covid_sleep_more_base","covid_wakeup_base",                   
                 "covid_15_base",  "covid_7_base",   "covid_rested_base", "covid_sleep_wk36", "covid_sleep_less_wk36", "covid_sleep_more_wk36", 
                 "covid_wakeup_wk36", "covid_15_wk36", "covid_7_wk36", "covid_rested_wk36" ,"covid_sleep_mth6", "covid_sleep_less_mth6", "covid_sleep_more_mth6",                 
                 "covid_wakeup_mth6", "covid_15_mth6" ,"covid_7_mth6", "covid_rested_mth6", "do_any_of_your_coworkers_s_base")

ID_Updated_Combined_Traits_apo_ana <- ID_Updated_Combined_Traits_apo[, !(names(ID_Updated_Combined_Traits_apo) %in% remove_vars)]

colnames(ID_Updated_Combined_Traits_apo_ana)

# View the updated dataframe
head(ID_Updated_Combined_Traits_apo_ana)

# Create a new variable for White_vs_nonwhite using case_when
ID_Updated_Combined_Traits_apo_ana <- ID_Updated_Combined_Traits_apo_ana %>%
  mutate(
    White_vs_nonwhite = case_when(
      race_eth_new == 4 ~ 1,        # If Race is "White", White_vs_nonwhite is 1
      TRUE ~ 0                    # If Race is anything else, White_vs_nonwhite is 0
    )
  )


# Create a new variable adverse_pregnancy_outcome using case_when
ID_Updated_Combined_Traits_apo_ana <- ID_Updated_Combined_Traits_apo_ana %>%
  mutate(apo_total = case_when(
    apo == 1 ~ 1,                 # If apo is 1, then adverse_pregnancy_outcome is 1
    is.na(apo) ~ 0,               # If apo is NA, then adverse_pregnancy_outcome is 0
    TRUE ~ 0                      # In all other cases, adverse_pregnancy_outcome is 0
  ))
#I have realized I dont these other variables, what Sebastian coded for is enough and the WGCNA plot looks ok.

#ID_Updated_Combined_Traits_apo_ana <- ID_Updated_Combined_Traits_apo_ana %>%
#  mutate(
 #   adverse_pregnancy_outcome_hdp = case_when(
  #    apo_hdp == 1 ~ 1,               # If apo_hdp is 1, then adverse_pregnancy_outcome_hdp is 1
   #   is.na(apo_hdp) ~ 0,             # If apo_hdp is NA, then adverse_pregnancy_outcome_hdp is 0
    #  TRUE ~ 0                        # In all other cases, adverse_pregnancy_outcome_hdp is 0
#    ),
 #   adverse_pregnancy_outcome_gdm = case_when(
  #    apo_gdm == 1 ~ 1,               # If apo_gdm is 1, then adverse_pregnancy_outcome_gdm is 1
   #   is.na(apo_gdm) ~ 0,             # If apo_gdm is NA, then adverse_pregnancy_outcome_gdm is 0
    #  TRUE ~ 0                        # In all other cases, adverse_pregnancy_outcome_gdm is 0
#    ),
 #   adverse_pregnancy_outcome_other = case_when(
  #    apo_other == 1 ~ 1,             # If apo_other is 1, then adverse_pregnancy_outcome_other is 1
   #   is.na(apo_other) ~ 0,           # If apo_other is NA, then adverse_pregnancy_outcome_other is 0
    #  TRUE ~ 0                        # In all other cases, adverse_pregnancy_outcome_other is 0
#    )
  #)

colnames(ID_Updated_Combined_Traits_apo_ana)

# View the updated dataframe
head(ID_Updated_Combined_Traits_apo_ana)


# Function to print the levels of factor variables along with their labels
print_factor_levels <- function(data) {
  factor_vars <- names(data)[sapply(data, is.factor)]
  for (var in factor_vars) {
    cat("Levels of", var, ":\n")
    print(levels(data[[var]]))
    cat("Labels of", var, ":\n")
    print(names(levels(data[[var]])))
    cat("\n")
  }
}

# Verify the changes by printing factor levels
print_factor_levels(ID_Updated_Combined_Traits_apo_ana)


# Convert factors to numeric 
ID_Updated_Combined_Traits_apo_ana[] <- lapply(ID_Updated_Combined_Traits_apo_ana, function(x) {
  if (is.factor(x)) as.numeric(as.character(x)) else x
})

# Check for any remaining character variables
character_vars <- names(ID_Updated_Combined_Traits_apo_ana)[sapply(ID_Updated_Combined_Traits_apo_ana, is.character)]
print(character_vars)  # Should return character(0) if no character variables are left

rownames(ID_Updated_Combined_Traits_apo_ana)
colnames(ID_Updated_Combined_Traits_apo_ana)

write.csv(ID_Updated_Combined_Traits_apo_ana, "ID_Updated_Combined_Traits.csv", row.names = TRUE)

#remove USER_ID
ID_Updated_Combined_Traits_apo_analy <- 
  ID_Updated_Combined_Traits_apo_analy[, !colnames(ID_Updated_Combined_Traits_apo_analy) %in% "User_ID"]

ID_Updated_Combined_Traits_apo_analy <- 
  ID_Updated_Combined_Traits_apo_analy[, !colnames(ID_Updated_Combined_Traits_apo_analy) %in% "apo"]


#rename apo_total to apo 

colnames(ID_Updated_Combined_Traits_apo_analy) <- gsub("^apo_total$", "apo", colnames(ID_Updated_Combined_Traits_apo_analy))


colnames(ID_Updated_Combined_Traits_apo_analy)

# Create a new dataframe with reordered columns
ID_Updated_Combined_Traits_apo_analysis <- ID_Updated_Combined_Traits_apo_analy %>%
  select(
    # Group 1: Demographic and Basic Information
    group, age, gest_age, gest_age_at_birth, preterm, height_base, weight_base, weight_wk26, weight_wk36, weight_mth3, weight_mth6, bmi, bmi_base, parity_base,gwg, 
    egwg,ppwr,ppwr_e, apo , apo_hdp, apo_gdm, apo_other,         
    
    
    # Group 2: Parental and Socioeconomic Information
    race_eth_new, White_vs_nonwhite, total_number_of_pregnancie_30f336_base, total_number_abortions_base, pregtype_conception_base, pregtype_planned_base,
    other_children_base, moth_marit_stat_base, mother_s_education_base, father_s_age_at_enrollment_base, fathers_country_of_birth_us_base, father_s_race_base,
    father_s_hispanic_base, father_s_education_base, father_currently_employed_base,low_income, low_access,
    
    # Group 3: Smoking and Environmental Exposures
    second_hand_base, second_hand_days_base, second_hand_rules_base, PM2.5_7Days_base, PM2.5_30Days_base, PM2.5_90Days_base,
    PM2.5_7Days_wk26, PM2.5_30Days_wk26, PM2.5_90Days_wk26, PM2.5_7Days_wk36, PM2.5_30Days_wk36, PM2.5_90Days_wk36, PM2.5_7Days_mth3,
    PM2.5_30Days_mth3, PM2.5_90Days_mth3, PM2.5_7Days_mth6, PM2.5_30Days_mth6, PM2.5_90Days_mth6, do_any_close_friends_or_re_base,
    
    # Group 4: Dietary Information
    KCAL_base, HEI2015C9_FATTYACID_base, HEI2015C12_SFAT_base, HEI2015_TOTAL_SCORE_base, KCAL_wk36, HEI2015C9_FATTYACID_wk36, HEI2015C12_SFAT_wk36,
    HEI2015_TOTAL_SCORE_wk36, KCAL_mth6, HEI2015C9_FATTYACID_mth6, HEI2015C12_SFAT_mth6, HEI2015_TOTAL_SCORE_mth6, reap_base, reap_wk26, reap_wk36, reap_mth3, reap_mth6,
    
    # Group 5: Physical Activity
    ppaq_total_base, ppaq_light_base, ppaq_mod_base, ppaq_house_base, ppaq_total_mth3, ppaq_light_mth3, ppaq_mod_mth3, ppaq_house_mth3,
    ppaq_total_mth6, ppaq_light_mth6, ppaq_mod_mth6, ppaq_house_mth6,
    
    # Group 6: Psychological Measures
    epds_base, epds_3a_base, epds_mth3, epds_3a_mth3, epds_mth6, epds_3a_mth6,
    sleep_base, sleep_wk36, sleep_mth6, rhodes_base,rhodes_wk26, rhodes_wk36, 
    
    # Group 7: COVID-related Information
    covid_essential_base, covid_test_yes_base, covid_schedule_base, covid_mask_base, covid_visit_base, covid_wipe_base,
    covid_essential_wk36, covid_test_yes_wk36, covid_schedule_wk36, covid_mask_wk36, covid_visit_wk36, covid_wipe_wk36,
    covid_essential_mth6, covid_test_yes_mth6, covid_schedule_mth6,
    
    # Group 8: Infant Information
    infant_gender_mth3, feeding_mth3, ex_feeding_mth3,feeding_mth6, ex_feeding_mth6, infant_weight_bbmr, breastfeeding,
    
    # Group 9: Tasks and Goals
    n_weights_pre, n_weights_post, n_weights_total, tasks_comp_pre, tasks_fail_pre, total_tasks_pre, prop_tasks_pre, tasks_comp_post, tasks_fail_post, total_tasks_post, prop_tasks_post,
    n_goals_pre, n_achieved_pre, prop_achieved_pre, n_goals_post, n_achieved_post, prop_achieved_post,
    do_you_prepare_your_own_fo_base, ever_have_trouble_being_ab_base, follow_a_special_diet_eat_base, how_willing_are_you_to_mak_base,
    do_you_prepare_your_own_fo_wk26, ever_have_trouble_being_ab_wk26, follow_a_special_diet_eat_wk26, how_willing_are_you_to_mak_wk26,
    do_you_prepare_your_own_fo_wk36, ever_have_trouble_being_ab_wk36, follow_a_special_diet_eat_wk36, how_willing_are_you_to_mak_wk36,
    do_you_prepare_your_own_fo_mth3, ever_have_trouble_being_ab_mth3, follow_a_special_diet_eat_mth3, how_willing_are_you_to_mak_mth3,
    do_you_prepare_your_own_fo_mth6, ever_have_trouble_being_ab_mth6, follow_a_special_diet_eat_mth6, how_willing_are_you_to_mak_mth6
    
  )

# View the new dataframe
View(ID_Updated_Combined_Traits_apo_analysis)

#check to see if you have all the variables. 
# Find variables in ID_Updated_Combined_Traits that are not in ID_Updated_Combined_Traits_apo_analysis
vars_only_in_original <- setdiff(colnames(ID_Updated_Combined_Traits_apo_analy), colnames(ID_Updated_Combined_Traits_apo_analysis))

# Find variables in ID_Updated_Combined_Traits_apo_analysis that are not in ID_Updated_Combined_Traits
vars_only_in_reordered <- setdiff(colnames(ID_Updated_Combined_Traits_apo_analysis), colnames(ID_Updated_Combined_Traits_apo_analy))

# Display the results
cat("Variables only in ID_Updated_Combined_Traits:\n")
print(vars_only_in_original)

cat("\nVariables only in ID_Updated_Combined_Traits_apo_analysis:\n")
print(vars_only_in_reordered)

#save your file
# Sort dataset by row names in ascending order
ID_Updated_Combined_Traits_apo_analysis <- ID_Updated_Combined_Traits_apo_analysis[order(rownames(ID_Updated_Combined_Traits_apo_analysis)), ]


write.csv(ID_Updated_Combined_Traits_apo_analysis, "ID_Updated_Combined_Traits_apo_analysis.csv", row.names = TRUE)

#=====================================================================================
#
#  Code chunk 1b- Separate Variables for analysis
#
#=====================================================================================
# Load your data

rownames(ID_Updated_Combined_Traits_apo_analysis)

# Sort dataset by row names in ascending order
ID_Updated_Combined_Traits_apo_analysis <- ID_Updated_Combined_Traits_apo_analysis[order(rownames(ID_Updated_Combined_Traits_apo_analysis)), ]

# View first few rows
head(ID_Updated_Combined_Traits_apo_analysis)

library(dplyr)

#  ID_Updated_Combined_Traits_apo_analysis is your dataframe
WGCNA_trial_variables <- ID_Updated_Combined_Traits_apo_analysis %>%
  select(group, age, gwg, egwg, ppwr, ppwr_e, father_s_education_base,
         do_any_close_friends_or_re_base, covid_schedule_base, 
         follow_a_special_diet_eat_mth6, weight_wk36, weight_mth3, weight_mth6,
         sleep_base, total_tasks_pre, ever_have_trouble_being_ab_wk26, covid_test_yes_mth6, 
         father_s_hispanic_base, gest_age, gest_age_at_birth, preterm, apo,
         moth_marit_stat_base, low_income, HEI2015C12_SFAT_mth6, tasks_fail_pre, prop_tasks_pre,
         how_willing_are_you_to_mak_wk36, apo_hdp, bmi_base, apo_hdp, ppaq_total_base, ppaq_mod_base, 
         ppaq_total_mth3, ppaq_light_mth3, ppaq_mod_mth3, ppaq_house_mth3, covid_mask_wk36,
         infant_gender_mth3, apo_gdm, rhodes_wk26, covid_schedule_mth6, weight_base, apo_other, fathers_country_of_birth_us_base, breastfeeding,
         n_weights_pre, n_weights_total, do_you_prepare_your_own_fo_wk26, do_you_prepare_your_own_fo_mth3)

# View the selected data
head(WGCNA_trial_variables)


# Group 1: Demographic and Basic Information #gest_age, for now removed
#gest_age, this variable is gest_age when signing up, I dont think it is
#informative
demographic_data <- ID_Updated_Combined_Traits_apo_analysis %>%
  select(group, age, gest_age_at_birth, preterm, height_base, weight_base, 
         weight_wk26, weight_wk36, weight_mth3, weight_mth6, bmi, bmi_base, parity_base, 
         gwg, egwg, ppwr, ppwr_e, apo , apo_hdp, apo_gdm, apo_other)

write.csv(demographic_data, "demographic_data.csv", row.names = TRUE)



#group 1b 
#add bmi for the different timepoints
demographic_data_bmi <- demographic_data %>%
  mutate(
    bmi_wk26 = weight_wk26 / (height_base^2),
    bmi_wk36 = weight_wk36 / (height_base^2),
    bmi_mth3 = weight_mth3 / (height_base^2),
    bmi_mth6 = weight_mth6 / (height_base^2)
  )

library(dplyr)

demographic_data_bmi <- demographic_data_bmi %>%
  select(
    weight_base, weight_wk26, weight_wk36, weight_mth3, weight_mth6,height_base, 
    # Move all bmi-related columns here
    bmi, bmi_base, bmi_wk26, bmi_wk36, bmi_mth3, bmi_mth6, 
    # Then keep the rest
    gwg, egwg, ppwr, ppwr_e,
    apo, preterm, apo_hdp, apo_gdm, apo_other
  )


write.xlsx(demographic_data_bmi, file = "demographic_data_bmi.xlsx", rowNames = TRUE)

# Group 2: Parental and Socioeconomic Information
parental_socioeconomic_data <- ID_Updated_Combined_Traits_apo_analysis %>%
  select(race_eth_new, White_vs_nonwhite, total_number_of_pregnancie_30f336_base, 
         total_number_abortions_base, pregtype_conception_base, pregtype_planned_base, 
         other_children_base, moth_marit_stat_base, mother_s_education_base, 
         father_s_age_at_enrollment_base, fathers_country_of_birth_us_base, 
         father_s_race_base, father_s_hispanic_base, father_s_education_base, 
         father_currently_employed_base, low_income, low_access)

# Group 3: Smoking and Environmental Exposures
smoking_environmental_data <- ID_Updated_Combined_Traits_apo_analysis %>%
  select(second_hand_base, second_hand_days_base, second_hand_rules_base, 
         PM2.5_7Days_base, PM2.5_30Days_base, PM2.5_90Days_base, PM2.5_7Days_wk26, 
         PM2.5_30Days_wk26, PM2.5_90Days_wk26, PM2.5_7Days_wk36, PM2.5_30Days_wk36, 
         PM2.5_90Days_wk36, PM2.5_7Days_mth3, PM2.5_30Days_mth3, PM2.5_90Days_mth3, 
         PM2.5_7Days_mth6, PM2.5_30Days_mth6, PM2.5_90Days_mth6, do_any_close_friends_or_re_base)

# Group 4: Dietary Information
dietary_data <- ID_Updated_Combined_Traits_apo_analysis %>%
  select(KCAL_base, HEI2015C9_FATTYACID_base, HEI2015C12_SFAT_base, 
         HEI2015_TOTAL_SCORE_base, KCAL_wk36, HEI2015C9_FATTYACID_wk36, 
         HEI2015C12_SFAT_wk36, HEI2015_TOTAL_SCORE_wk36, KCAL_mth6, 
         HEI2015C9_FATTYACID_mth6, HEI2015C12_SFAT_mth6, HEI2015_TOTAL_SCORE_mth6, 
         reap_base, reap_wk26, reap_wk36, reap_mth3, reap_mth6)

# Group 5: Physical Activity
physical_activity_data <- ID_Updated_Combined_Traits_apo_analysis %>%
  select(ppaq_total_base, ppaq_light_base, ppaq_mod_base, ppaq_house_base, 
         ppaq_total_mth3, ppaq_light_mth3, ppaq_mod_mth3, ppaq_house_mth3, 
         ppaq_total_mth6, ppaq_light_mth6, ppaq_mod_mth6, ppaq_house_mth6)

# Group 6: Psychological Measures
psychological_measures_data <- ID_Updated_Combined_Traits_apo_analysis %>%
  select(epds_base, epds_3a_base, epds_mth3, epds_3a_mth3, epds_mth6, epds_3a_mth6, 
         sleep_base, sleep_wk36, sleep_mth6, rhodes_base, rhodes_wk26, rhodes_wk36)

# Group 7: COVID-related Information
covid_data <- ID_Updated_Combined_Traits_apo_analysis %>%
  select(covid_essential_base, covid_test_yes_base, covid_schedule_base, 
         covid_mask_base, covid_visit_base, covid_wipe_base, covid_essential_wk36, 
         covid_test_yes_wk36, covid_schedule_wk36, covid_mask_wk36, covid_visit_wk36, 
         covid_wipe_wk36, covid_essential_mth6, covid_test_yes_mth6, covid_schedule_mth6)

# Group 8: Infant Information
infant_data <- ID_Updated_Combined_Traits_apo_analysis %>%
  select(infant_gender_mth3, feeding_mth3, ex_feeding_mth3, feeding_mth6, 
         ex_feeding_mth6, infant_weight_bbmr, breastfeeding)

# Group 9: Tasks and Goals
tasks_goals_data <- ID_Updated_Combined_Traits_apo_analysis %>%
  select(n_weights_pre, n_weights_post, n_weights_total, tasks_comp_pre, 
         tasks_fail_pre, total_tasks_pre, prop_tasks_pre, tasks_comp_post, 
         tasks_fail_post, total_tasks_post, prop_tasks_post, n_goals_pre, 
         n_achieved_pre, prop_achieved_pre, n_goals_post, n_achieved_post, 
         prop_achieved_post)

write.csv(tasks_goals_data, "tasks_goals_data.csv", row.names = TRUE)


# Group 10: Readiness to change
readines_data <- ID_Updated_Combined_Traits_apo_analysis   %>%
  select (do_you_prepare_your_own_fo_base, 
         ever_have_trouble_being_ab_base, follow_a_special_diet_eat_base, 
         how_willing_are_you_to_mak_base, do_you_prepare_your_own_fo_wk26, 
         ever_have_trouble_being_ab_wk26, follow_a_special_diet_eat_wk26, 
         how_willing_are_you_to_mak_wk26, do_you_prepare_your_own_fo_wk36, 
         ever_have_trouble_being_ab_wk36, follow_a_special_diet_eat_wk36, 
         how_willing_are_you_to_mak_wk36, do_you_prepare_your_own_fo_mth3, 
         ever_have_trouble_being_ab_mth3, follow_a_special_diet_eat_mth3, 
         how_willing_are_you_to_mak_mth3, do_you_prepare_your_own_fo_mth6, 
         ever_have_trouble_being_ab_mth6, follow_a_special_diet_eat_mth6, 
         how_willing_are_you_to_mak_mth6)

write.csv(readines_data, "readines_data.csv", row.names = TRUE)

# Check each dataset
head(demographic_data)
head(demographic_data_bmi)
head(parental_socioeconomic_data)
head(smoking_environmental_data)
head(dietary_data)
head(physical_activity_data)
head(psychological_measures_data)
head(covid_data)
head(infant_data)
head(tasks_goals_data)
head(readines_data)

#This trait file changes depending with the question
pheno_ID_pre <- list(B = list(data = demographic_data_bmi),
                     TSTE = list(data = demographic_data_bmi),
                     PP = list(data = demographic_data_bmi))
#=====================================================================================
#
#  Code chunk 1c- Make plots to show APO vs other traits of interest 
#
#=====================================================================================
# Group 1: Demographic and Basic Information
demographic_data <- ID_Updated_Combined_Traits_apo_analysis %>%
  select(group, age, gest_age, gest_age_at_birth, preterm, height_base, weight_base, 
         weight_wk26, weight_wk36, weight_mth3, weight_mth6, bmi, bmi_base, parity_base, 
         gwg, egwg, ppwr, ppwr_e, apo , apo_hdp, apo_gdm, apo_other)

#add bmi for the different timepoints
demographic_data <- demographic_data %>%
mutate(
  bmi_wk26 = weight_wk26 / (height_base^2),
  bmi_wk36 = weight_wk36 / (height_base^2),
  bmi_mth3 = weight_mth3 / (height_base^2),
  bmi_mth6 = weight_mth6 / (height_base^2)
)

# change apo to factor
demographic_data <- demographic_data %>%
  mutate(
    apo = factor(apo, levels = c(0, 1), labels = c("No APO", "APO"))
  )
# plot with no stats
library(ggplot2)
library(dplyr)
# plot 1. This worked by it have some grid lines and had to add the jitter dots
# they were also black. 
p <- ggplot(demographic_data, aes(x = apo, y = gwg, fill = apo)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "white", outlier.color = "black") +
  scale_fill_manual(values = c("No APO" = "#3366CC", "APO" = "#FF3366")) +
  labs(
    title = "APO Effect on Magenta Module - GWG",
    x = "APO Status",
    y = "Gestational Weight Gain"
  ) +
  theme_minimal()

print(p)

p <- p + 
  geom_jitter(width = 0.2, shape = 16, alpha = 0.6, color = "black")

# Optionally save
# ggsave("APO_vs_Magenta_GWG_Boxplot.pdf", plot = p, height = 7, width = 6)

library(ggplot2)
library(dplyr)
# plot 2. This worked by it didnt put a square panel
# Ensure 'apo' is properly formatted
demographic_data <- demographic_data %>%
  mutate(apo = factor(apo, levels = c(0, 1), labels = c("No APO", "APO")))

# Custom color mapping
trait_colors <- c("No APO" = "#3366CC", "APO" = "#FF3366")

# Plot
p <- ggplot(demographic_data, aes(x = apo, y = gwg, fill = apo)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, color = "black") +
  geom_jitter(aes(color = apo), width = 0.15, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = trait_colors) +
  scale_color_manual(values = trait_colors) +
  labs(
    title = "GWG by APO Status",
    x = "APO Status",
    y = "Gestational Weight Gain (GWG)"
  ) +
  theme_classic(base_size = 14) +  # white background
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(p)


library(ggplot2)
library(dplyr)
#plot 3. This worked. Depending with variable you might need to adjust 
#coord_fixed(ratio = 0.1,. increase or decrease and see what works
# Ensure 'apo' is formatted as a factor
demographic_data <- demographic_data %>%
  mutate(apo = factor(apo, levels = c(0, 1), labels = c("No APO", "APO")))

# Define colors
trait_colors <- c("No APO" = "#3366CC", "APO" = "red")

# Calculate GWG y-axis range
ymin <- 0
ymax <- ceiling(max(demographic_data$gwg, na.rm = TRUE) + 1)

# Plot GWG with square axis area
p <- ggplot(demographic_data, aes(x = apo, y = gwg, fill = apo)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.6, color = "black") +
  geom_jitter(aes(color = apo), width = 0.15, size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = trait_colors, name = "APO Status") +
  scale_color_manual(values = trait_colors, name = "APO Status") +
  labs(
    title = "GWG by APO Status",
    x = "APO Status",
    y = "Gestational Weight Gain (GWG)"
  ) +
  #coord_fixed(ratio = 0.1, ylim = c(ymin, ymax), clip = "on") +  # make axis box square
  coord_fixed(ratio = 0.1, ylim = c(-1, 20), clip = "on") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# Save or print
print(p)
ggsave("GWG_APO_Square_Axis_Boxplot.pdf", plot = p, width = 8, height = 8)

#Now I want to add some stats to it. I chose Wilcox because most likely it is 
# not normally distributed

continuous_traits <- c(
  "age", "gest_age", "gest_age_at_birth", "height_base", "weight_base", 
  "weight_wk26", "weight_wk36", "weight_mth3", "weight_mth6", 
  "bmi", "bmi_base", "parity_base", "gwg", "ppwr", 
)
#I added a calculate of bmi to the begining and then add the other traits
continuous_traits <- c(
  "age", "gest_age", "gest_age_at_birth", "height_base", "weight_base", 
  "weight_wk26", "weight_wk36", "weight_mth3", "weight_mth6", 
  "bmi", "bmi_base", "parity_base", "gwg", "ppwr", 
  "bmi_wk26", "bmi_wk36", "bmi_mth3", "bmi_mth6"
)


t_test_results <- lapply(continuous_traits, function(trait) {
  formula <- as.formula(paste(trait, "~ apo"))
  result <- t.test(formula, data = demographic_data)
  tibble(
    Trait = trait,
    t_statistic = result$statistic,
    p_value = result$p.value,
    mean_No_APO = result$estimate[["mean in group No APO"]],
    mean_APO = result$estimate[["mean in group APO"]]
  )
})

# Combine into a single table
t_test_df <- bind_rows(t_test_results)

print(t_test_df)

#Now I want to try add some stats to this
#plot 1
#if you did factor, trait color and ymin no need to repeat
demographic_data <- demographic_data %>%
  mutate(apo = factor(apo, levels = c(0, 1), labels = c("No APO", "APO")))

# Define colors
trait_colors <- c("No APO" = "#3366CC", "APO" = "red")

# Calculate GWG y-axis range
ymin <- 0
ymax <- ceiling(max(demographic_data$gwg, na.rm = TRUE) + 1)

# Extract stats for GWG
gwg_test <- t_test_df %>% filter(Trait == "gwg")
p_val <- gwg_test$p_value
stars <- ifelse(p_val < 0.01, "**", ifelse(p_val < 0.05, "*", ""))
title_text <- paste0("GWG by APO Status (p = ", signif(p_val, 3), ")", stars)

# Final plot
p <- ggplot(demographic_data, aes(x = apo, y = gwg, fill = apo)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.6, color = "black") +
  geom_jitter(aes(color = apo), width = 0.15, size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = trait_colors, name = "APO Status") +
  scale_color_manual(values = trait_colors, name = "APO Status") +
  labs(
    title = title_text,
    x = "APO Status",
    y = "Gestational Weight Gain (GWG)"
  ) +
  coord_fixed(ratio = 0.1, ylim = c(-1, 20), clip = "on") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(p)
ggsave("GWG_APO_Square_Axis__stats_Boxplot.pdf", plot = p, width = 8, height = 8)

library(ggplot2)
library(dplyr)
#keep in mind that plot 2 and 3 both work, they just have different utilities

#Now I want to try add a function that goes through all my samples
#plot 2 The problem with the function is it generate all plots for my traits but
# i cant just go and manipulate one if there was an error of some sort.
# Colors
#if you did factor, trait color and ymin no need to repeat
demographic_data <- demographic_data %>%
  mutate(apo = factor(apo, levels = c(0, 1), labels = c("No APO", "APO")))

# Define colors
trait_colors <- c("No APO" = "#3366CC", "APO" = "red")

# Loop over each continuous trait and plot
for (trait in continuous_traits) {
  
  # Safely skip if trait is missing
  if (!(trait %in% names(demographic_data))) next
  
  # Extract t-test p-value
  test_row <- t_test_df %>% filter(Trait == trait)
  if (nrow(test_row) == 0) next
  
  p_val <- test_row$p_value
  stars <- ifelse(p_val < 0.01, "**", ifelse(p_val < 0.05, "*", ""))
  title_text <- paste0(trait, " by APO Status (p = ", signif(p_val, 3), ")", stars)
  
  # Set reasonable y-axis limits
  trait_vals <- demographic_data[[trait]]
  ylim_min <- floor(min(trait_vals, na.rm = TRUE)) - 1
  ylim_max <- ceiling(max(trait_vals, na.rm = TRUE)) + 1
  
  # Plot
  p <- ggplot(demographic_data, aes_string(x = "apo", y = trait, fill = "apo")) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.6, color = "black") +
    geom_jitter(aes_string(color = "apo"), width = 0.15, size = 1.5, alpha = 0.7) +
    scale_fill_manual(values = trait_colors, name = "APO Status") +
    scale_color_manual(values = trait_colors, name = "APO Status") +
    labs(
      title = title_text,
      x = "APO Status",
      y = trait
    ) +
    coord_fixed(ratio = 0.1, ylim = c(ylim_min, ylim_max), clip = "on") +
    theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  print(p)
  
  # Save
  ggsave(paste0("Boxplot_", trait, "_APO.pdf"), plot = p, width = 12, height = 12)
}
#plot3
#make function  that allows you to manipulate plots if 
#you want to adjust something
plot_trait_boxplot <- function(data, trait, t_test_df,
                               colors = c("No APO" = "#3366CC", "APO" = "red"),
                               y_limits = NULL,
                               ratio = 0.1,
                               title_prefix = "",
                               width = 6,
                               height = 6,
                               save_path = NULL) {
  
  # Check if trait exists
  if (!(trait %in% names(data))) {
    warning(paste("Trait", trait, "not found in dataset."))
    return(NULL)
  }
  
  # Extract t-test info
  test_row <- t_test_df %>% filter(Trait == trait)
  if (nrow(test_row) == 0) {
    warning(paste("No t-test result for trait:", trait))
    return(NULL)
  }
  
  p_val <- test_row$p_value
  stars <- ifelse(p_val < 0.01, "**", ifelse(p_val < 0.05, "*", ""))
  title_text <- paste0(title_prefix, trait, " by APO Status (p = ", signif(p_val, 3), ")", stars)
  
  # Y-axis limits
  if (is.null(y_limits)) {
    values <- data[[trait]]
    y_min <- floor(min(values, na.rm = TRUE)) - 1
    y_max <- ceiling(max(values, na.rm = TRUE)) + 1
    y_limits <- c(y_min, y_max)
  }
  
  # Build plot
  p <- ggplot(data, aes_string(x = "apo", y = trait, fill = "apo")) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.6, color = "black") +
    geom_jitter(aes_string(color = "apo"), width = 0.15, size = 1.5, alpha = 0.7) +
    scale_fill_manual(values = colors, name = "APO Status") +
    scale_color_manual(values = colors, name = "APO Status") +
    labs(
      title = title_text,
      x = "APO Status",
      y = trait
    ) +
    coord_fixed(ratio = ratio, ylim = y_limits, clip = "on") +
    theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  # Save or print
  if (!is.null(save_path)) {
    ggsave(save_path, plot = p, width = width, height = height)
  } else {
    print(p)
  }
  
  return(p)
}

#default
plot_trait_boxplot(
  data = demographic_data,
  trait = "gwg",
  t_test_df = t_test_df
)

#adjusting limits custom version
plot_trait_boxplot(
  data = demographic_data,
  trait = "gwg",
  t_test_df = t_test_df,
  colors = c("No APO" = "turquoise", "APO" = "red"),
  y_limits = c(-5, 20),
  ratio = 0.1,
  title_prefix = "Adjusted: ",
  save_path = "Adjusted_Boxplot_GWG_APO.pdf"
)

#plot 4
#now the goal is to be able to put traits into panels/facets for publications
#
library(ggplot2)
library(dplyr)
library(patchwork)

facet_trait_boxplots <- function(data, traits, t_test_df,
                                 colors = c("No APO" = "#3366CC", "APO" = "red"),
                                 panel_labels = LETTERS,
                                 ratio = 0.1,
                                 y_limits_list = NULL,
                                 title_prefix = "",
                                 width = 10,
                                 height = 10,
                                 ncol = 2,
                                 save_path = NULL) {
  
  plot_list <- list()
  
  for (i in seq_along(traits)) {
    trait <- traits[i]
    
    if (!(trait %in% names(data))) {
      warning(paste("Skipping missing trait:", trait))
      next
    }
    
    # Get t-test row
    test_row <- t_test_df %>% filter(Trait == trait)
    if (nrow(test_row) == 0) {
      warning(paste("No t-test result for:", trait))
      next
    }
    
    p_val <- test_row$p_value
    stars <- ifelse(p_val < 0.01, "**", ifelse(p_val < 0.05, "*", ""))
    title_text <- paste0(panel_labels[i], ". ", title_prefix, trait, " (p = ", signif(p_val, 3), ")", stars)
    
    # Define y-limits
    if (!is.null(y_limits_list) && !is.null(y_limits_list[[trait]])) {
      ylims <- y_limits_list[[trait]]
    } else {
      vals <- data[[trait]]
      ylims <- c(floor(min(vals, na.rm = TRUE)) - 1, ceiling(max(vals, na.rm = TRUE)) + 1)
    }
    
    # Plot
    p <- ggplot(data, aes_string(x = "apo", y = trait, fill = "apo")) +
      geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.6, color = "black") +
      geom_jitter(aes_string(color = "apo"), width = 0.15, size = 1.5, alpha = 0.7) +
      scale_fill_manual(values = colors, name = "APO Status") +
      scale_color_manual(values = colors, name = "APO Status") +
      labs(
        title = title_text,
        x = "APO Status",
        y = trait
      ) +
      coord_cartesian(ylim = ylims, clip = "on") + #this was added to ensure sameheight despite different ranges
     # coord_fixed(ratio = ratio, ylim = ylims, clip = "on") + #you can remove coord_cartesian if you dont want it to have same height with different ranges and use this part
      theme_bw(base_size = 12) +
      theme(
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5)
      )
    
    plot_list[[i]] <- p
  }
  
  # Combine panels
  combined_plot <- wrap_plots(plot_list, ncol = ncol)  
  
  # Save or print
  if (!is.null(save_path)) {
    ggsave(save_path, plot = combined_plot, width = width, height = height)
  } else {
    print(combined_plot)
  }
  
  return(combined_plot)
}

#examples lets do bmi
bmi_traits <- c("bmi","bmi_base", "bmi_wk26", "bmi_wk36", "bmi_mth3", "bmi_mth6")

facet_trait_boxplots(
  data = demographic_data,
  traits = bmi_traits,
  t_test_df = t_test_df,
  save_path = "Panel_BMI_APO.pdf"
)


#this can help adjust width and height and y-limits
y_limits_list <- setNames(
  rep(list(c(25, 50)), length(weight_traits)),
  weight_traits
)

facet_trait_boxplots(
  data = demographic_data,
  traits = bmi_traits,
  t_test_df = t_test_df,
  y_limits_list = y_limits_list,
  save_path = "Panel_BMI_APO.pdf",
  ncol = 3,
  width = 13,
  height = 6.5
)

#weight
weight_traits <- c( "weight_base",  "weight_wk26",  "weight_wk36", "weight_mth3",
  "weight_mth6"
)


#weight
# I want to adjust  ylimits. i checked min and max
# Define y-limits for each weight trait
y_limits_list <- setNames(
  rep(list(c(60, 135)), length(weight_traits)),
  weight_traits
)

facet_trait_boxplots(
  data = demographic_data,
  traits = weight_traits,
  t_test_df = t_test_df,
  y_limits_list = y_limits_list,
  ratio = 0.02,
  ncol = 3,
  save_path = "Panel_weight_APO.pdf",
  width = 13,
  height = 6.5
)

#pregnancy traits
y_limits_list <- list(
  "parity_base" = c(0, 8),              # taller because wider y-range
  "gest_age_at_birth" = c(30, 43),      # also tall
  "gwg" = c(-5, 20),                     # medium height
  "ppwr" = c(-5, 20)                    # shorter panel
)

pregnancy_traits <- c("parity_base",      
  "gest_age_at_birth", "gwg", "ppwr"
)

facet_trait_boxplots(
  data = demographic_data,
  traits = pregnancy_traits,
  t_test_df = t_test_df,
  #y_limits_list = y_limits_list,
  ratio = 0.07,
  save_path = "Panel_Preg_traits_APO.pdf",
  width = 13,
  height = 6.5
)


#=====================================================================================
#
#  Code chunk 1c- Separate lipidomics variables for  analysis
#
#=====================================================================================
#load data
identified_lipids_baseline_data_normalized <- read_csv("Identified_lipids_used_in_metaboanalyst_to_normalize/Identified_baseline/identified_lipids_baseline_data_normalized.csv")
identified_lipids_TP36_data_normalized <- read_csv("Identified_lipids_used_in_metaboanalyst_to_normalize/identified_TP36-38/identified_lipids_TP36_data_normalized.csv")
identified_postpartum_data_normalized <- read_csv("Identified_lipids_used_in_metaboanalyst_to_normalize/Identified_Postpartum/identified_postpartum_data_normalized.csv")

dim(identified_lipids_baseline_data_normalized)
# Convert to data frame and transpose
identified_lipids_baseline_data_normalized_transposed <- as.data.frame(t(identified_lipids_baseline_data_normalized))
identified_lipids_TP36_data_normalized_transposed <- as.data.frame(t(identified_lipids_TP36_data_normalized))
identified_postpartum_data_normalized_transposed <- as.data.frame(t(identified_postpartum_data_normalized))

# Optionally, rename columns based on original row names
colnames(identified_lipids_baseline_data_normalized_transposed) <- identified_lipids_baseline_data_normalized_transposed[1, ]
identified_lipids_baseline_data_normalized_transposed <- identified_lipids_baseline_data_normalized_transposed[-1, ]

colnames(identified_lipids_TP36_data_normalized_transposed) <- identified_lipids_TP36_data_normalized_transposed[1, ]
identified_lipids_TP36_data_normalized_transposed <- identified_lipids_TP36_data_normalized_transposed[-1, ]

colnames(identified_postpartum_data_normalized_transposed) <- identified_postpartum_data_normalized_transposed[1, ]
identified_postpartum_data_normalized_transposed <- identified_postpartum_data_normalized_transposed[-1, ]

#sort the datasets by rownames

identified_lipids_baseline_data_normalized_transposed <- 
  identified_lipids_baseline_data_normalized_transposed[order(rownames(identified_lipids_baseline_data_normalized_transposed)), ]

head(identified_lipids_baseline_data_normalized_transposed)

identified_lipids_TP36_data_normalized_transposed <- 
  identified_lipids_TP36_data_normalized_transposed[order(rownames(identified_lipids_TP36_data_normalized_transposed)), ]

head(identified_lipids_TP36_data_normalized_transposed)

identified_postpartum_data_normalized_transposed <- 
  identified_postpartum_data_normalized_transposed[order(rownames(identified_postpartum_data_normalized_transposed)), ]

head(identified_postpartum_data_normalized_transposed)


# Save the transposed data (optional)
write.csv(identified_lipids_baseline_data_normalized_transposed, "identified_lipids_baseline_data_normalized_transposed.csv", row.names = TRUE)
write.csv(identified_lipids_TP36_data_normalized_transposed, "identified_lipids_TP36_data_normalized_transposed.csv", row.names = TRUE)
write.csv(identified_postpartum_data_normalized_transposed, "identified_postpartum_data_normalized_transposed.csv", row.names = TRUE)

#=====================================================================================
#
#  Code chunk 1d- Separate lipidomics variables for  analysis using metaboanalyst on outcomes 
#
#=====================================================================================
colnames(identified_lipids_baseline_data_normalized_transposed)
colnames(identified_lipids_baseline_data_normalized_transposed)
colnames(identified_lipids_baseline_data_normalized_transposed)
rownames(identified_lipids_baseline_data_normalized_transposed)

head(demographic_data)
demographic_data$egwg_status
# Recode 1 to "EGWG_YES", everything else to "EGWG_NO"
demographic_data$egwg_status <- ifelse(demographic_data$egwg == 1, "EGWG_YES", "EGWG_NO")
demographic_data$ppwr_status <- ifelse(demographic_data$ppwr_e == 1, "PPWR_YES", "PPWR_NO")
demographic_data$apo_status <- ifelse(demographic_data$apo == 1, "APO_YES","APO_NO")
demographic_data$apo_hdp_status <- ifelse(demographic_data$apo_hdp == 1, "APO_HDP_YES", "APO_HDP_NO")
demographic_data$apo_gdm_status <- ifelse(demographic_data$apo_gdm == 1, "APO_GDM_YES", "APO_GDM_NO")
demographic_data$apo_other_status <- ifelse(demographic_data$apo_other == 1, "APO_OTHER_YES","APO_OTHER_NO" )

# Create a new dataset with only the *_status variables
status_vars <- demographic_data %>%
  select(egwg_status, ppwr_status, apo_status, apo_hdp_status, apo_gdm_status, apo_other_status)

# Remove those variables from the original dataset
demographic_data <- demographic_data %>%
  select(-egwg_status, -ppwr_status, -apo_status, -apo_hdp_status, -apo_gdm_status, -apo_other_status)

# egwg_status first
identified_lipids_baseline_with_egwg <- cbind(egwg_status = status_vars$egwg_status,
                                              identified_lipids_baseline_data_normalized_transposed_ana)

identified_lipids_TP36_with_egwg <- cbind(egwg_status = status_vars$egwg_status,
                                          identified_lipids_TP36_data_normalized_transposed_ana)

identified_postpartum_with_egwg <- cbind(egwg_status = status_vars$egwg_status,
                                         identified_postpartum_data_normalized_transposed_ana)

# ppwr_status first
identified_lipids_baseline_with_ppwr <- cbind(ppwr_status = status_vars$ppwr_status,
                                              identified_lipids_baseline_data_normalized_transposed_ana)

identified_lipids_TP36_with_ppwr <- cbind(ppwr_status = status_vars$ppwr_status,
                                          identified_lipids_TP36_data_normalized_transposed_ana)

identified_postpartum_with_ppwr <- cbind(ppwr_status = status_vars$ppwr_status,
                                         identified_postpartum_data_normalized_transposed_ana)

# apo_status first
identified_lipids_baseline_with_apo <- cbind(apo_status = status_vars$apo_status,
                                             identified_lipids_baseline_data_normalized_transposed_ana)

identified_lipids_TP36_with_apo <- cbind(apo_status = status_vars$apo_status,
                                         identified_lipids_TP36_data_normalized_transposed_ana)

identified_postpartum_with_apo <- cbind(apo_status = status_vars$apo_status,
                                        identified_postpartum_data_normalized_transposed_ana)

# apo_hdp_status first
identified_lipids_baseline_with_apo_hdp <- cbind(apo_hdp_status = status_vars$apo_hdp_status,
                                                 identified_lipids_baseline_data_normalized_transposed_ana)

identified_lipids_TP36_with_apo_hdp <- cbind(apo_hdp_status = status_vars$apo_hdp_status,
                                             identified_lipids_TP36_data_normalized_transposed_ana)

identified_postpartum_with_apo_hdp <- cbind(apo_hdp_status = status_vars$apo_hdp_status,
                                            identified_postpartum_data_normalized_transposed_ana)

# apo_gdm_status first
identified_lipids_baseline_with_apo_gdm <- cbind(apo_gdm_status = status_vars$apo_gdm_status,
                                                 identified_lipids_baseline_data_normalized_transposed_ana)

identified_lipids_TP36_with_apo_gdm <- cbind(apo_gdm_status = status_vars$apo_gdm_status,
                                             identified_lipids_TP36_data_normalized_transposed_ana)

identified_postpartum_with_apo_gdm <- cbind(apo_gdm_status = status_vars$apo_gdm_status,
                                            identified_postpartum_data_normalized_transposed_ana)

# apo_other_status first
identified_lipids_baseline_with_apo_other <- cbind(apo_other_status = status_vars$apo_other_status,
                                                   identified_lipids_baseline_data_normalized_transposed_ana)

identified_lipids_TP36_with_apo_other <- cbind(apo_other_status = status_vars$apo_other_status,
                                               identified_lipids_TP36_data_normalized_transposed_ana)

identified_postpartum_with_apo_other <- cbind(apo_other_status = status_vars$apo_other_status,
                                              identified_postpartum_data_normalized_transposed_ana)


# Write CSV files for egwg_status
write.csv(identified_lipids_baseline_with_egwg, "identified_lipids_baseline_with_egwg.csv", row.names = TRUE)
write.csv(identified_lipids_TP36_with_egwg, "identified_lipids_TP36_with_egwg.csv", row.names = TRUE)
write.csv(identified_postpartum_with_egwg, "identified_postpartum_with_egwg.csv", row.names = TRUE)

# Write CSV files for ppwr_status
write.csv(identified_lipids_baseline_with_ppwr, "identified_lipids_baseline_with_ppwr.csv", row.names = TRUE)
write.csv(identified_lipids_TP36_with_ppwr, "identified_lipids_TP36_with_ppwr.csv", row.names = TRUE)
write.csv(identified_postpartum_with_ppwr, "identified_postpartum_with_ppwr.csv", row.names = TRUE)

# Write CSV files for apo_status
write.csv(identified_lipids_baseline_with_apo, "identified_lipids_baseline_with_apo.csv", row.names = TRUE)
write.csv(identified_lipids_TP36_with_apo, "identified_lipids_TP36_with_apo.csv", row.names = TRUE)
write.csv(identified_postpartum_with_apo, "identified_postpartum_with_apo.csv", row.names = TRUE)

# Write CSV files for apo_hdp_status
write.csv(identified_lipids_baseline_with_apo_hdp, "identified_lipids_baseline_with_apo_hdp.csv", row.names = TRUE)
write.csv(identified_lipids_TP36_with_apo_hdp, "identified_lipids_TP36_with_apo_hdp.csv", row.names = TRUE)
write.csv(identified_postpartum_with_apo_hdp, "identified_postpartum_with_apo_hdp.csv", row.names = TRUE)

# Write CSV files for apo_gdm_status
write.csv(identified_lipids_baseline_with_apo_gdm, "identified_lipids_baseline_with_apo_gdm.csv", row.names = TRUE)
write.csv(identified_lipids_TP36_with_apo_gdm, "identified_lipids_TP36_with_apo_gdm.csv", row.names = TRUE)
write.csv(identified_postpartum_with_apo_gdm, "identified_postpartum_with_apo_gdm.csv", row.names = TRUE)

# Write CSV files for apo_other_status
write.csv(identified_lipids_baseline_with_apo_other, "identified_lipids_baseline_with_apo_other.csv", row.names = TRUE)
write.csv(identified_lipids_TP36_with_apo_other, "identified_lipids_TP36_with_apo_other.csv", row.names = TRUE)
write.csv(identified_postpartum_with_apo_other, "identified_postpartum_with_apo_other.csv", row.names = TRUE)



#load data

#for analysis remove the label/group column. This is necessary when using metaboanalysis
#but we will not need it for now.It is the first column
# Optionally, rename columns based on original row names
colnames(identified_lipids_baseline_data_normalized_transposed)
colnames(identified_lipids_baseline_data_normalized_transposed)
colnames(identified_lipids_baseline_data_normalized_transposed)
rownames(identified_lipids_baseline_data_normalized_transposed)


identified_lipids_baseline_data_normalized_transposed_ana <- 
  identified_lipids_baseline_data_normalized_transposed[, !colnames(identified_lipids_baseline_data_normalized_transposed) %in% "Label"]

identified_lipids_TP36_data_normalized_transposed_ana <- 
  identified_lipids_TP36_data_normalized_transposed[, !colnames(identified_lipids_TP36_data_normalized_transposed) %in% "Label"]

identified_postpartum_data_normalized_transposed_ana <- 
  identified_postpartum_data_normalized_transposed[, !colnames(identified_postpartum_data_normalized_transposed) %in% "Label"]

#convert data to numeric which means add the X to some lipid names eg16.16.dimethylprostaglandin.A1 
#is X16.16.dimethylprostaglandin.A1 because I got the error below
#Error in blockwiseIndividualTOMs(multiExpr = multiExpr, checkMissingData = checkMissingData,  : 
 #                                  REAL() can only be applied to a 'numeric', not a 'character'
str(colnames(identified_lipids_baseline_data_normalized_transposed_ana))
str(colnames(identified_lipids_TP36_data_normalized_transposed_ana))
str(colnames(identified_postpartum_data_normalized_transposed_ana))

#
identified_lipids_baseline_data_normalized_transposed_ana <- data.frame(
  lapply(identified_lipids_baseline_data_normalized_transposed_ana, function(x) as.numeric(as.character(x))),
  row.names = rownames(identified_lipids_baseline_data_normalized_transposed_ana)
)

identified_lipids_TP36_data_normalized_transposed_ana <- data.frame(
  lapply(identified_lipids_TP36_data_normalized_transposed_ana, function(x) as.numeric(as.character(x))),
  row.names = rownames(identified_lipids_TP36_data_normalized_transposed_ana)
)

identified_postpartum_data_normalized_transposed_ana <- data.frame(
  lapply(identified_postpartum_data_normalized_transposed_ana, function(x) as.numeric(as.character(x))),
  row.names = rownames(identified_postpartum_data_normalized_transposed_ana)
)
#check it worked
class(identified_lipids_baseline_data_normalized_transposed_ana[[1]])  # Should be "data.frame"
class(identified_lipids_TP36_data_normalized_transposed_ana[[1]])  # Should be "data.frame"
class(identified_postpartum_data_normalized_transposed_ana[[1]])  # Should be "data.frame"

#check again to see if it worked
rownames(identified_lipids_baseline_data_normalized_transposed_ana)
colnames(identified_lipids_baseline_data_normalized_transposed_ana)
colnames(identified_lipids_TP36_data_normalized_transposed_ana)
colnames(identified_postpartum_data_normalized_transposed_ana)

multiExpr_ID <- list( Identified_Baseline = list(data = as.data.frame(identified_lipids_baseline_data_normalized_transposed_ana)),
                      Identified_TP36_38weeks = list(data = as.data.frame(identified_lipids_TP36_data_normalized_transposed_ana)),
                      Identified_Postpartum = list(data = as.data.frame(identified_postpartum_data_normalized_transposed_ana)))

exprSize_ID = checkSets(multiExpr_ID)
All_nSets_ID = exprSize_ID$nSets

#=====================================================================================
#
#  Code chunk 2- Run WGCNA and do the scale free network topology to choose the power
#
#=====================================================================================

# Load the WGCNA package
library(WGCNA)

#Start running WGCNA
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
nSets_ID = checkSets(multiExpr_ID)$nSets

setLabels_ID = c("10-16 weeks", "36-38weeks", "Postpartum") 
powers_ID = c(seq(4,10,by=1), seq(12,40, by=2))
powerTables_ID = vector(mode = "list", length = nSets_ID)
for (set in 1:nSets_ID) {
  powerTables_ID[[set]] = list(data = pickSoftThreshold(multiExpr_ID[[set]]$data, blockSize = 1000,
                                                        powerVector = powers_ID,
                                                        corFnc = "bicor",
                                                        networkType = "signed",
                                                        verbose = 2)[[2]])
}
collectGarbage()

colors_ID = c("black", "red", "green") 
plotCols_ID = c(2,5,6,7)
colNames_ID = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
                "Max connectivity")

ylim_ID = matrix(NA, nrow = 2, ncol = 4)
for (set in 1:nSets_ID) {
  for (col in 1:length(plotCols_ID)) {
    ylim_ID[1, col] = min(ylim_ID[1, col], powerTables_ID[[set]]$data[, plotCols_ID[col]], na.rm = TRUE)
    ylim_ID[2, col] = max(ylim_ID[2, col], powerTables_ID[[set]]$data[, plotCols_ID[col]], na.rm = TRUE)
  }
}

sizeGrWindow(8, 6)
pdf(file = "scaleFreeAnalysis_signed_ID.pdf", wi = 8, he = 6)
par(mfcol = c(2,2))
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7
for (col in 1:length(plotCols_ID)) {
  for (set in 1:nSets_ID) {
    if (set == 1) {
      plot(powerTables_ID[[set]]$data[,1], -sign(powerTables_ID[[set]]$data[,3])*powerTables_ID[[set]]$data[,2],
           xlab = "Soft Threshold (power)", ylab = colNames_ID[col], type = "n", ylim = ylim_ID[, col],
           main = colNames_ID[col])
      abline(h = c(0.80, 0.90), col = "red")
      addGrid()
    }
    if (col == 1) {
      text(powerTables_ID[[set]]$data[,1], -sign(powerTables_ID[[set]]$data[,3])*powerTables_ID[[set]]$data[,2],
           labels = powers_ID, cex = cex1, col = colors_ID[set])
    } else {
      text(powerTables_ID[[set]]$data[,1], powerTables_ID[[set]]$data[,plotCols_ID[col]],
           labels = powers_ID, cex = cex1, col = colors_ID[set])
    }
    if (col == 1) {
      legend("bottomright", legend = setLabels_ID, col = colors_ID, pch = 20)
    } else {
      legend("topright", legend = setLabels_ID, col = colors_ID, pch = 20)
    }
  }
}
dev.off()

### soft power to be used on all dataset has to be the same so I chose 24. This based on cutoff being 0.8 on the scale free topology model fit

### e.g exp <- multiExpr$Baseline$data exp <- multiExpr$TP36_38weeks$data exp <- multiExpr$Postpartum$data

#consensusMods_ID <- blockwiseConsensusModules(multiExpr_ID, checkMissingData = FALSE, maxBlockSize = 3000, corType = "bicor",
#                                             maxPOutliers = 0.1, power = 24, networkType = "signed", 
#                                            checkPower = FALSE, TOMType = "signed", 
#                                           networkCalibration = "full quantile", saveConsensusTOMs = TRUE,
#                                          deepSplit = 4, mergeCutHeight = 0.25, verbose = 5)
#table(consensusMods_ID$colors) %>% sort(decreasing = TRUE)
#module.dist_ID <- as.data.frame(table(consensusMods_ID$colors) %>% sort(decreasing = TRUE))
#colnames(module.dist_ID) <- c("Module", "Lipids")
#write.csv(module.dist_ID,"no_minimummodule.distribution_consensus_signed_ID.csv")

#Consensus Module analysis. I added minModule size to get more modules
consensusMods_ID <- blockwiseConsensusModules(multiExpr_ID, checkMissingData = FALSE, maxBlockSize = 3000, corType = "bicor",
                                              maxPOutliers = 0.1, power = 24, networkType = "signed", 
                                              checkPower = FALSE, TOMType = "signed", 
                                              networkCalibration = "full quantile", saveConsensusTOMs = TRUE,
                                              deepSplit = 4, minModuleSize = 5, mergeCutHeight = 0.1, verbose = 5)

table(consensusMods_ID$colors) %>% sort(decreasing = TRUE)
module.dist_ID <- as.data.frame(table(consensusMods_ID$colors) %>% sort(decreasing = TRUE))
colnames(module.dist_ID) <- c("Module", "Lipids")
write.csv(module.dist_ID,"module.distribution_consensus_signed_ID.csv")


# Plot Merged Gene Dendrogram with Modules
pdf("cluster.dendrogram.wgcnatrial_ID.pdf", width = 10, height = 5)
sizeGrWindow(10, 5)
plotDendroAndColors(dendro = consensusMods_ID$dendrograms[[1]], colors = consensusMods_ID$colors, 
                    groupLabels = "Modules", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, 
                    guideHang = 0.05, marAll = c(1, 5, 1, 0), main = "", cex.colorLabels = 1.3)
dev.off()

# Baseline Cluster Modules by Eigenlipids and Plot Dendrogram
METree_ID <- (1 - bicor(consensusMods_ID$multiMEs$Identified_Baseline$data, maxPOutliers = 0.1)) %>% as.dist %>% 
  hclust(method = "average")
pdf("Baseline_Module_Eigennode_Dendrogram_ID.pdf", height = 5, width = 10)
sizeGrWindow(height = 5, width = 10)
par(mar = c(0, 5, 1, 1))
plot(METree_ID, main = "", xlab = "", sub = "", ylim = c(0, 1), cex = 0.6)
abline(h = 0.1, col = "red")
dev.off()

# 36-38 weeks Cluster Modules by EigenLipids and Plot Dendrogram
METree_ID <- (1 - bicor(consensusMods_ID$multiMEs$Identified_TP36_38weeks$data, maxPOutliers = 0.1)) %>% as.dist %>% 
  hclust(method = "average")
pdf("36-38_weeks_Module_Eigennode_Dendrogram_ID.pdf", height = 5, width = 10)
sizeGrWindow(height = 5, width = 10)
par(mar = c(0, 5, 1, 1))
plot(METree_ID, main = "", xlab = "", sub = "", ylim = c(0, 1), cex = 0.6)
abline(h = 0.1, col = "red")
dev.off()

# Postpartum Cluster Modules by EigenLipids and Plot Dendrogram
METree_ID <- (1 - bicor(consensusMods_ID$multiMEs$Identified_Postpartum$data, maxPOutliers = 0.1)) %>% as.dist %>% 
  hclust(method = "average")
pdf("Postpartum_Module_Eigennode_Dendrogram_ID.pdf", height = 5, width = 10)
sizeGrWindow(height = 5, width = 10)
par(mar = c(0, 5, 1, 1))
plot(METree_ID, main = "", xlab = "", sub = "", ylim = c(0, 1), cex = 0.6)
abline(h = 0.1, col = "red")
dev.off()

consensusMEs_ID <- consensusOrderMEs(consensusMods_ID$multiMEs)
pdf(file = "Consensus_WGCNA_Eigennode_Networks_ID.pdf", width = 8, height = 7)
sizeGrWindow(width = 8, height = 7)
par(cex = 0.8)
plotEigengeneNetworks(consensusMEs_ID, setLabels = c("10-16 weeks", "36-38 weeks", "Postpartum"), 
                      plotDendrograms = FALSE, marHeatmap = c(3, 3, 2, 1), zlimPreservation = c(0.5, 1), 
                      xLabelsAngle = 90)
dev.off()

#=====================================================================================
#
#  Code chunk 3- Do Module Membership
#
#=====================================================================================

# Different section for module membership
#Baseline
moduleMembership_ID <- mtd.mapply(bicorAndPvalue, multiExpr_ID, consensusMEs_ID, 
                                  MoreArgs = list(alternative = "two.sided", use = "pairwise.complete.obs", 
                                                  maxPOutliers = 0.1))
Consensus_Modules_Baseline_ID_Probe_Module_Membership <- as.data.frame(moduleMembership_ID$Identified_Baseline$data$bicor)
colnames(Consensus_Modules_Baseline_ID_Probe_Module_Membership) <- gsub(pattern = "ME", replacement = "", x = colnames(Consensus_Modules_Baseline_ID_Probe_Module_Membership), fixed = TRUE)
Consensus_Modules_Baseline_ID_Probe_Module_Membership$Probe <- rownames(Consensus_Modules_Baseline_ID_Probe_Module_Membership)
Consensus_Modules_Baseline_ID_Probe_Module_Membership$Module <- consensusMods_ID$colors
write.table(Consensus_Modules_Baseline_ID_Probe_Module_Membership, "Consensus_Modules_Baseline_ID_Probe_Module_Membership.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write_xlsx(Consensus_Modules_Baseline_ID_Probe_Module_Membership, "Consensus_Modules_Baseline_ID_Probe_Module_Membership.xlsx")

#Week 36
Consensus_Modules_36_38_ID_Probe_Module_Membership <- as.data.frame(moduleMembership_ID$Identified_TP36_38weeks$data$bicor)
colnames(Consensus_Modules_36_38_ID_Probe_Module_Membership) <- gsub(pattern = "ME", replacement = "", x = colnames(Consensus_Modules_36_38_ID_Probe_Module_Membership), fixed = TRUE)
Consensus_Modules_36_38_ID_Probe_Module_Membership$Probe <- rownames(Consensus_Modules_36_38_ID_Probe_Module_Membership)
Consensus_Modules_36_38_ID_Probe_Module_Membership$Module <- consensusMods_ID$colors
write.table(Consensus_Modules_36_38_ID_Probe_Module_Membership, "Consensus_Modules_36_38_ID_Probe_Module_Membership.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write_xlsx(Consensus_Modules_36_38_ID_Probe_Module_Membership, "Consensus_Modules_36_38_ID_Probe_Module_Membership.xlsx")

#Postpartum
Consensus_Modules_Postpartum_ID_Probe_Module_Membership <- as.data.frame(moduleMembership_ID$Identified_Postpartum$data$bicor)
colnames(Consensus_Modules_Postpartum_ID_Probe_Module_Membership) <- gsub(pattern = "ME", replacement = "", x = colnames(Consensus_Modules_Postpartum_ID_Probe_Module_Membership), fixed = TRUE)
Consensus_Modules_Postpartum_ID_Probe_Module_Membership$Probe <- rownames(Consensus_Modules_Postpartum_ID_Probe_Module_Membership)
Consensus_Modules_Postpartum_ID_Probe_Module_Membership$Module <- consensusMods_ID$colors
write.table(Consensus_Modules_Postpartum_ID_Probe_Module_Membership, "Consensus_Modules_Postpartum_ID_Probe_Module_Membership.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write_xlsx(Consensus_Modules_Postpartum_ID_Probe_Module_Membership, "Consensus_Modules_Postpartum_ID_Probe_Module_Membership.xlsx")


#=====================================================================================
#
#  Code chunk 4- Identify which lipids are in specific modules
#
#=====================================================================================
#I want to first extract the inCHIKey from the original file data then I want to combine
# it with the PubChem_CID. I used the Pubchem ID exchange to get the PubChem CID.

Raw_data <- read_excel("Raw_lipidomics.xlsx")
colnames(Raw_data)
rownames(Raw_data)

Lipid_PubChem_CID_InCHIKey_Sheet1_ <- read_csv("Lipid_PubChem_CID_InCHIKey(Sheet1).csv")

# Remove the first 8 rows
Raw_data_sh <- Raw_data[-c(1:7), ]

# Set the first row as column names
colnames(Raw_data_sh) <- Raw_data_sh[1, ]

# Remove the first row
Raw_data_sho <- Raw_data_sh[-1, ]

# View the cleaned dataset
head(Raw_data_sho)

# Select the specific columns to create a new dataset
Raw_data_selected <- Raw_data_sho %>%
  select(annotation, InChiKey)

# View the first few rows of the new dataset
head(Raw_data_selected)


# Sort by 'annotation' (or any other column you prefer)
Raw_data_sorted <- Raw_data_selected %>%
  arrange(annotation)

# View sorted dataset
head(Raw_data_sorted)

# Select rows 1 to 451
Raw_data_identified_lipids <- Raw_data_sorted %>%
  slice(1:451)

# View the first few rows
head(Raw_data_identified_lipids)


# Perform full join on the 'InChiKey' column
Full_Joined_Data <- full_join(Lipid_PubChem_CID_InCHIKey_Sheet1_, Raw_data_identified_lipids, by = "InChiKey")

# View the first few rows of the merged dataset
head(Full_Joined_Data)

#how many are distinct/unique
Full_Joined_Data %>%
  distinct(annotation) %>%
  count()

# Create new dataset with distinct annotation names
Distinct_Annotations <- Full_Joined_Data %>%
  distinct(annotation, .keep_all = TRUE)

# View the first few rows
head(Distinct_Annotations)

# Sort the dataset alphabetically by annotation
Distinct_Annotations_sorted <- Distinct_Annotations %>%
  arrange(annotation)

# View the first few rows
head(Distinct_Annotations_sorted)
rownames(Distinct_Annotations_sorted)
Distinct_Annotations_sorted$rownames <- 1:451

#now lets combine with the module membership
head(Consensus_Modules_Baseline_ID_Probe_Module_Membership)
View(Consensus_Modules_Baseline_ID_Probe_Module_Membership)
colnames(Consensus_Modules_Baseline_ID_Probe_Module_Membership)

#first create new dataset with Proble and Module
Probe_Module_Only <- Consensus_Modules_Baseline_ID_Probe_Module_Membership %>%
  select(Probe, Module)

# Add new column with numbers 1 to 451
Probe_Module_Only$rownames <- 1:451

#merge the two datasets
# Make sure both datasets have a column named 'rownames' and are numeric
Distinct_Annotations_sorted$rownames <- as.numeric(Distinct_Annotations_sorted$rownames)
Probe_Module_Only$rownames <- as.numeric(Probe_Module_Only$rownames)

# Perform full join using the 'rownames' column
Lipid_annotations_module <- full_join(Distinct_Annotations_sorted, Probe_Module_Only, by = "rownames")

write.csv(Lipid_annotations_module, "Lipid_annotations_module_before_correction.csv", row.names = TRUE)

#Now load the corrected version plus load the Ref_Met results. This is website that 
#took the annotation name and helped match the lipids. I tried to use super-class based on
#Lipid_Maps

#load data
Ref_Met_Results <- read_csv("Ref_Met_Results(Sheet1).csv")
Lipid_annotations_module_corrected <- read_csv("Lipid_annotations_module_corrected.csv")

colnames(Ref_Met_Results)
#merge the datasets using annotations columns
Ref_Met_Lipid_annotations_Module <- full_join(Ref_Met_Results, Lipid_annotations_module_corrected, by =  "annotation")

colnames(Ref_Met_Lipid_annotations_Module) <- gsub("[ _]+", "_", colnames(Ref_Met_Lipid_annotations_Module))

colnames(Ref_Met_Lipid_annotations_Module)

Ref_Met_Lipid_annotations_Module$Main_class


#save
write.csv(Ref_Met_Lipid_annotations_Module, "Ref_Met_Lipid_annotations_Module.csv", row.names = TRUE)

lipid_ref_met_class <- Ref_Met_Lipid_annotations_Module %>%
  select(Super_class, Main_class, Sub_class, Probe)

row.names(lipid_ref_met_class) <- lipid_ref_met_class$Probe
colnames(lipid_ref_met_class)
lipid_ref_met_class$Probe

# Code to assign the lipid class and module to see how different classes of lipids are distributed in different modules
# Load data
# First we do it for Baseline

# Merge with module membership
Consensus_Baseline_lipid_class_ID <- merge(Consensus_Modules_Baseline_ID_Probe_Module_Membership, lipid_ref_met_class, by = "Probe")
head(Consensus_Baseline_lipid_class_ID) # Check the merged data
colnames(Consensus_Baseline_lipid_class_ID)
write_xlsx(Consensus_Baseline_lipid_class_ID, "Lipid_Class_Summary_all_classes_and_modules_ID.xlsx")

#summary
# Create a named list: one element per main class
lipid_list_by_class <- split(Consensus_Baseline_lipid_class_ID$Probe,
                             Consensus_Baseline_lipid_class_ID$Main_class)

# Example: View lipids in the 'Triradylglycerols' class
lipid_list_by_class[["Triradylglycerols"]]


View(Consensus_Modules_Baseline_ID_Probe_Module_Membership)
View(lipid_ref_met_class)
View(Consensus_Baseline_lipid_class_ID)


# Summarize the number of lipids per class within each module for super class
lipid_summary_Super_class_ID <- Consensus_Baseline_lipid_class_ID %>%
  group_by(Module, Super_class) %>%
  summarise(Count = n()) %>%
  ungroup()

# Convert to a data frame
lipid_summary_Super_class_ID <- as.data.frame(lipid_summary_Super_class_ID)
head(lipid_summary_Super_class_ID) # Check the summary data
write_xlsx(lipid_summary_Super_class_ID, "Lipid_Module_Super_Class_Summary_ID.xlsx")

library(reshape2)
library(pheatmap)
library(viridis)
library(ComplexHeatmap)
library(circlize)
install.packages("circlize")
BiocManager::install("circlize")
BiocManager::install("ComplexHeatmap")

# Create matrix
super_matrix <- acast(lipid_summary_Super_class_ID, Super_class ~ Module, value.var = "Count", fill = 0)

# Create color scale
super_colors <- colorRamp2(c(0, 1, max(super_matrix)), c("white", "#44AA99", "#117733"))

# Save PDF
pdf("Heatmap_Super_Class_Distribution.pdf", width = 12, height = 9)

# Plot heatmap
Heatmap(super_matrix,
        name = "Count",
        col = super_colors,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        column_title = "Super Class Distribution Across Modules",
        cell_fun = function(j, i, x, y, width, height, fill) {
          val <- super_matrix[i, j]
          grid.text(
            label = sprintf("%.0f", val),
            x = x, y = y,
            gp = gpar(
              col = ifelse(val == 0, "black", "black"), 
              fontsize = 10
            )
          )
        })

dev.off()


# Summarize the number of lipids per class within each module for main class
lipid_summary_Main_class_ID <- Consensus_Baseline_lipid_class_ID %>%
  group_by(Module, Main_class) %>%
  summarise(Count = n()) %>%
  ungroup()

# Convert to a data frame
lipid_summary_Main_class_ID <- as.data.frame(lipid_summary_Main_class_ID)
head(lipid_summary_Main_class_ID) # Check the summary data
write_xlsx(lipid_summary_Main_class_ID, "Lipid_Module_Main_Class_Summary_ID.xlsx")



# Create a color palette for values 1–89 (white for 0)
heatmap_colors <- colorRamp2(c(0, 1, 89), c("white", "#44AA99", "#117733"))  # 
  
# Open PDF
pdf("Heatmap_Main_Class_Distribution.pdf", width = 12, height = 9)

# Plot heatmap
Heatmap(main_matrix,
        name = "Count",
        col = heatmap_colors,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        column_title = "Main Class Distribution Across Modules",
        cell_fun = function(j, i, x, y, width, height, fill) {
          val <- main_matrix[i, j]
          grid.text(
            label = sprintf("%.0f", val),
            x = x, y = y,
            gp = gpar(
              col = ifelse(val == 0, "black", "black"),  # conditional text color
              fontsize = 10
            )
          )
        })

# Close PDF
dev.off()

# Summarize the number of lipids per class within each module for subclass
lipid_summary_Sub_class_ID <- Consensus_Baseline_lipid_class_ID %>%
  group_by(Module, Sub_class) %>%
  summarise(Count = n()) %>%
  ungroup()

# Convert to a data frame
lipid_summary_Sub_class_ID <- as.data.frame(lipid_summary_Sub_class_ID)
head(lipid_summary_Sub_class_ID) # Check the summary data
write_xlsx(lipid_summary_Sub_class_ID, "Lipid_Module_Sub_Class_Summary_ID.xlsx")


# Create matrix
sub_matrix <- acast(lipid_summary_Sub_class_ID, Sub_class ~ Module, value.var = "Count", fill = 0)

# Create color scale
sub_colors <- colorRamp2(c(0, 1, max(sub_matrix)), c("white", "#44AA99", "#117733"))

# Save PDF
pdf("Heatmap_Sub_Class_Distribution.pdf", width = 12, height = 9)

# Plot heatmap
Heatmap(sub_matrix,
        name = "Count",
        col = sub_colors,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        column_title = "Sub Class Distribution Across Modules",
        cell_fun = function(j, i, x, y, width, height, fill) {
          val <- sub_matrix[i, j]
          grid.text(
            label = sprintf("%.0f", val),
            x = x, y = y,
            gp = gpar(
              col = ifelse(val == 0, "black", "black"), 
              fontsize = 10
            )
          )
        })

dev.off()

# Open a PDF device
pdf("Heatmap_Sub_Class_Distribution.pdf", width = 12, height = 9)


# Define a custom color palette
custom_colors <- colorRampPalette(c("white", "yellow", "orange"))(47)

# Open PDF
pdf("Heatmap_Sub_Class_Distribution.pdf", width = 12, height = 9)

# Create matrix
sub_matrix <- acast(lipid_summary_Sub_class_ID, Sub_class ~ Module, value.var = "Count", fill = 0)

# Plot
pheatmap(sub_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = custom_colors,
         display_numbers = TRUE,
         number_color = "black",  # works well with custom colors
         main = "Heatmap: Sub Class Distribution Across Modules",
         legend_breaks = c(0, 24, 47),
         legend_labels = c("0", "24", "47"),
         fontsize_number = 10)

# Close PDF device
dev.off()

# Create a bar plot for each module
# This plot was ugly, see heatmap below
ggplot(lipid_summary_ID, aes(x = Module, y = Count, fill = Lipid_class)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(title = "Lipid Class Distribution in Each Module", 
       x = "Module", 
       y = "Count of Lipids", 
       fill = "Lipid Class") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(reshape2)
library(pheatmap)
library(viridis)
# Create a matrix for the heatmap
lipid_matrix_ID <- acast(lipid_summary_ID, Lipid_class ~ Module, value.var = "Count", fill = 0)

# Create a custom color palette
# Define a custom color palette to highlight non-zero values distinctly
custom_colors <- c("white", viridis(100))

# Generate the heatmap
pheatmap(lipid_matrix_ID, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         color = custom_colors,
         display_numbers = TRUE,
         number_color = "red",
         main = "Heatmap of Lipid Class Distribution Across Modules",
         legend_breaks = c(0, 1, 10, 50, 100, 200, 300, 400, 500, 600, 700),
         legend_labels = c("0", "1", "10", "50", "100", "200", "300", "400", "500", "600", "700"),
         fontsize_number = 10)


#=====================================================================================
#
#  Code chunk 5- Find your Hub_lipids
#
#=====================================================================================

# Find the hub lipid start with Baseline
rownames(Consensus_Modules_Baseline_ID_Probe_Module_Membership) <- Consensus_Modules_Baseline_ID_Probe_Module_Membership$Probe 
colnames(Consensus_Modules_Baseline_ID_Probe_Module_Membership)

# Note that Probe means the lipid name and module means the color. 
# Baseline for just 1
hubLipids_Base_ID <- sapply(colnames(Consensus_Modules_Baseline_ID_Probe_Module_Membership)[!colnames(Consensus_Modules_Baseline_ID_Probe_Module_Membership) %in%  c("Probe", "Module")], function(x){
  temp_Base_ID <- Consensus_Modules_Baseline_ID_Probe_Module_Membership[Consensus_Modules_Baseline_ID_Probe_Module_Membership$Module == x,]
  temp_Base_ID$Probe[temp_Base_ID[, x] == max(temp_Base_ID[, x])] %>% as.character %>% unique %>% sort})
hubLipids_Base_ID

write.csv(hubLipids_Base_ID, "ConsensushubLipids_Baseline_ID.csv")

# Set rownames to lipid names
rownames(Consensus_Modules_Baseline_ID_Probe_Module_Membership) <- Consensus_Modules_Baseline_ID_Probe_Module_Membership$Probe 

# Extract top 10 hub lipids per module with rank
hubLipids_Baseline_long <- lapply(
  colnames(Consensus_Modules_Baseline_ID_Probe_Module_Membership)[
    !colnames(Consensus_Modules_Baseline_ID_Probe_Module_Membership) %in% c("Probe", "Module")
  ],
  function(mod) {
    temp <- Consensus_Modules_Baseline_ID_Probe_Module_Membership[
      Consensus_Modules_Baseline_ID_Probe_Module_Membership$Module == mod, 
    ]
    
    top_indices <- order(temp[[mod]], decreasing = TRUE)[1:min(10, nrow(temp))]
    
    data.frame(
      Module = mod,
      HubLipid = temp$Probe[top_indices],
      Rank = seq_along(top_indices),
      stringsAsFactors = FALSE
    )
  }
) %>% bind_rows()

# Save to CSV
write.csv(hubLipids_Baseline_long, "ConsensushubLipids_Baseline_Top10.csv", row.names = FALSE)

# Find the hub lipid then week 36-38
rownames(Consensus_Modules_36_38_ID_Probe_Module_Membership) <- Consensus_Modules_36_38_ID_Probe_Module_Membership$Probe 
rownames(Consensus_Modules_36_38_ID_Probe_Module_Membership)
colnames(Consensus_Modules_36_38_ID_Probe_Module_Membership)

# 36-38 weeks just 1 lipid
hubLipids_Week36_ID <- sapply(colnames(Consensus_Modules_36_38_ID_Probe_Module_Membership)[!colnames(Consensus_Modules_36_38_ID_Probe_Module_Membership) %in%  c("Probe", "Module")], function(x){
  temp_Week36_ID <- Consensus_Modules_36_38_ID_Probe_Module_Membership[Consensus_Modules_36_38_ID_Probe_Module_Membership$Module == x,]
  temp_Week36_ID$Probe[temp_Week36_ID[, x] == max(temp_Week36_ID[, x])] %>% as.character %>% unique %>% sort})
hubLipids_Week36_ID

write.csv(hubLipids_Week36_ID, "ConsensushubLipids_Week36_ID.csv")

# Ensure probe names are rownames
rownames(Consensus_Modules_36_38_ID_Probe_Module_Membership) <- Consensus_Modules_36_38_ID_Probe_Module_Membership$Probe 

# Extract top 10 hub lipids per module
hubLipids_36_38_long <- lapply(
  colnames(Consensus_Modules_36_38_ID_Probe_Module_Membership)[
    !colnames(Consensus_Modules_36_38_ID_Probe_Module_Membership) %in% c("Probe", "Module")
  ],
  function(mod) {
    temp <- Consensus_Modules_36_38_ID_Probe_Module_Membership[
      Consensus_Modules_36_38_ID_Probe_Module_Membership$Module == mod, 
    ]
    
    top_indices <- order(temp[[mod]], decreasing = TRUE)[1:min(10, nrow(temp))]
    
    data.frame(
      Module = mod,
      HubLipid = temp$Probe[top_indices],
      Rank = seq_along(top_indices),
      stringsAsFactors = FALSE
    )
  }
) %>% bind_rows()

# Save to CSV
write.csv(hubLipids_36_38_long, "ConsensushubLipids_Week36_38_Top10.csv", row.names = FALSE)




# Find the hub lipid Postpartum
rownames(Consensus_Modules_Postpartum_ID_Probe_Module_Membership) <- Consensus_Modules_Postpartum_ID_Probe_Module_Membership$Probe 
colnames(Consensus_Modules_Postpartum_ID_Probe_Module_Membership)

# Postpartum
hubLipids_Month3_ID <- sapply(colnames(Consensus_Modules_Postpartum_ID_Probe_Module_Membership)[!colnames(Consensus_Modules_Postpartum_ID_Probe_Module_Membership) %in%  c("Probe", "Module")], function(x){
  temp_Month3_ID <- Consensus_Modules_Postpartum_ID_Probe_Module_Membership[Consensus_Modules_Postpartum_ID_Probe_Module_Membership$Module == x,]
  temp_Month3_ID$Probe[temp_Month3_ID[, x] == max(temp_Month3_ID[, x])] %>% as.character %>% unique %>% sort})
hubLipids_Month3_ID

write.csv(hubLipids_Month3_ID, "ConsensushubLipids_Month3_ID.csv")

# Set rownames to lipid names
rownames(Consensus_Modules_Postpartum_ID_Probe_Module_Membership) <- Consensus_Modules_Postpartum_ID_Probe_Module_Membership$Probe 

# Extract top 10 hub lipids per module with rank
hubLipids_Month3_long <- lapply(
  colnames(Consensus_Modules_Postpartum_ID_Probe_Module_Membership)[
    !colnames(Consensus_Modules_Postpartum_ID_Probe_Module_Membership) %in% c("Probe", "Module")
  ],
  function(mod) {
    temp <- Consensus_Modules_Postpartum_ID_Probe_Module_Membership[
      Consensus_Modules_Postpartum_ID_Probe_Module_Membership$Module == mod, 
    ]
    
    top_indices <- order(temp[[mod]], decreasing = TRUE)[1:min(10, nrow(temp))]
    
    data.frame(
      Module = mod,
      HubLipid = temp$Probe[top_indices],
      Rank = seq_along(top_indices),
      stringsAsFactors = FALSE
    )
  }
) %>% bind_rows()

# Save to CSV
write.csv(hubLipids_Month3_long, "ConsensushubLipids_Month3_Top10.csv", row.names = FALSE)




#=====================================================================================
#
#  Code chunk 5b- Find your Hub_lipids Plots
#
#=====================================================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
#Baseline


# Get top 10 hub lipids per module at Baseline
hubLipids_Base_ID <- lapply(
  colnames(Consensus_Modules_Baseline_ID_Probe_Module_Membership)[
    !colnames(Consensus_Modules_Baseline_ID_Probe_Module_Membership) %in% c("Probe", "Module")
  ],
  function(mod) {
    if (mod == "grey") return(NULL)
    temp <- Consensus_Modules_Baseline_ID_Probe_Module_Membership[
      Consensus_Modules_Baseline_ID_Probe_Module_Membership$Module == mod, 
    ]
    top_indices <- order(temp[[mod]], decreasing = TRUE)[1:min(10, nrow(temp))]
    data.frame(Module = mod, HubLipid = temp$Probe[top_indices], Rank = seq_along(top_indices))
  }
) %>% bind_rows()

hubLipids_Base_ID
# Format
hubLipids_Base_ID$Module <- fct_inorder(hubLipids_Base_ID$Module)
hubLipids_Base_ID$HubLipid <- substr(hubLipids_Base_ID$HubLipid, 1, 15)
module_labels <- distinct(hubLipids_Base_ID, Module)
fill_colors <- setNames(levels(hubLipids_Base_ID$Module), levels(hubLipids_Base_ID$Module))
text_colors <- fill_colors
text_colors[c("lightcyan","green","cyan", "lightyellow", "yellow","grey60", "pink", "tan", "lightgreen", "salmon", "greenyellow")] <- "black"

# Plot
ggplot(hubLipids_Base_ID, aes(x = Rank, y = Module)) +
  geom_tile(data = module_labels, aes(x = 0.5, y = Module, fill = Module),
            width = 0.4, height = 0.8, inherit.aes = FALSE) +
  geom_text(aes(label = HubLipid, color = Module), size = 3.5, fontface = "bold") +
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = text_colors) +
  scale_x_continuous(breaks = 1:10, name = "Rank", expand = expansion(add = c(0.5, 0.5))) +
  scale_y_discrete(name = "Module") +
  labs(title = "Top 10 Hub Lipids per Module (10-16 weeks)") +
  theme_minimal() +
  theme(axis.text.y = element_text(face = "bold"), panel.grid = element_blank(), legend.position = "none")

ggsave("Top10_Baseline_HubLipids_PerModule.pdf", width = 19.19, height = 11.26, units = "in")


#36-38 weeks
hubLipids_Week36_ID <- lapply(
  colnames(Consensus_Modules_36_38_ID_Probe_Module_Membership)[
    !colnames(Consensus_Modules_36_38_ID_Probe_Module_Membership) %in% c("Probe", "Module")
  ],
  function(mod) {
    if (mod == "grey") return(NULL)
    temp <- Consensus_Modules_36_38_ID_Probe_Module_Membership[
      Consensus_Modules_36_38_ID_Probe_Module_Membership$Module == mod, 
    ]
    top_indices <- order(temp[[mod]], decreasing = TRUE)[1:min(10, nrow(temp))]
    data.frame(Module = mod, HubLipid = temp$Probe[top_indices], Rank = seq_along(top_indices))
  }
) %>% bind_rows()


hubLipids_Week36_ID$Module <- fct_inorder(hubLipids_Week36_ID$Module)
hubLipids_Week36_ID$HubLipid <- substr(hubLipids_Week36_ID$HubLipid, 1, 15)
module_labels <- distinct(hubLipids_Week36_ID, Module)
fill_colors <- setNames(levels(hubLipids_Week36_ID$Module), levels(hubLipids_Week36_ID$Module))
text_colors <- fill_colors
text_colors[c("lightcyan", "lightyellow","green","cyan", "yellow", "pink","grey60", "tan", "lightgreen", "salmon", "greenyellow")] <- "black"

ggplot(hubLipids_Week36_ID, aes(x = Rank, y = Module)) +
  geom_tile(data = module_labels, aes(x = 0.5, y = Module, fill = Module),
            width = 0.4, height = 0.8, inherit.aes = FALSE) +
  geom_text(aes(label = HubLipid, color = Module), size = 3.5, fontface = "bold") +
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = text_colors) +
  scale_x_continuous(breaks = 1:10, name = "Rank", expand = expansion(add = c(0.5, 0.5))) +
  scale_y_discrete(name = "Module") +
  labs(title = "Top 10 Hub Lipids per Module (36–38 Weeks)") +
  theme_minimal() +
  theme(axis.text.y = element_text(face = "bold"), panel.grid = element_blank(), legend.position = "none")

ggsave("Top10_Week36_HubLipids_PerModule.pdf", width = 19.19, height = 11.26, units = "in")






#postpartum
hubLipids_Month3_ID <- lapply(
  colnames(Consensus_Modules_Postpartum_ID_Probe_Module_Membership)[
    !colnames(Consensus_Modules_Postpartum_ID_Probe_Module_Membership) %in% c("Probe", "Module")
  ],
  function(mod) {
    if (mod == "grey") return(NULL)
    temp <- Consensus_Modules_Postpartum_ID_Probe_Module_Membership[
      Consensus_Modules_Postpartum_ID_Probe_Module_Membership$Module == mod, 
    ]
    top_indices <- order(temp[[mod]], decreasing = TRUE)[1:min(10, nrow(temp))]
    data.frame(Module = mod, HubLipid = temp$Probe[top_indices], Rank = seq_along(top_indices))
  }
) %>% bind_rows()


hubLipids_Month3_ID$Module <- fct_inorder(hubLipids_Month3_ID$Module)
hubLipids_Month3_ID$HubLipid <- substr(hubLipids_Month3_ID$HubLipid, 1, 15)
module_labels <- distinct(hubLipids_Month3_ID, Module)
fill_colors <- setNames(levels(hubLipids_Month3_ID$Module), levels(hubLipids_Month3_ID$Module))
text_colors <- fill_colors
text_colors[c("lightcyan", "lightyellow", "yellow","green","cyan", "pink", "tan", "grey60", "lightgreen", "salmon", "greenyellow")] <- "black"

ggplot(hubLipids_Month3_ID, aes(x = Rank, y = Module)) +
  geom_tile(data = module_labels, aes(x = 0.5, y = Module, fill = Module),
            width = 0.4, height = 0.8, inherit.aes = FALSE) +
  geom_text(aes(label = HubLipid, color = Module), size = 3.5, fontface = "bold") +
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = text_colors) +
  scale_x_continuous(breaks = 1:10, name = "Rank", expand = expansion(add = c(0.5, 0.5))) +
  scale_y_discrete(name = "Module") +
  labs(title = "Top 10 Hub Lipids per Module (3 Months Postpartum)") +
  theme_minimal() +
  theme(axis.text.y = element_text(face = "bold"), panel.grid = element_blank(), legend.position = "none")

ggsave("Top10_Month3_HubLipids_PerModule.pdf", width = 19.19, height = 11.26, units = "in")


#=====================================================================================
#
#  Code chunk 6- Setup module-trait correlations
#
#=====================================================================================

moduleTraitCor_ID = list()
moduleTraitPvalue_ID = list()

#Error in .useNThreads() : could not find function ".useNThreads" so i decided to specifically ask it to use WGCNA
#for (set in 1:All_nSets_ID) {
 # moduleTraitCor_ID[[set]] = cor(consensusMEs_ID[[set]]$data, pheno_ID_pre[[set]]$data, method = "pearson", use = "p")
  #moduleTraitPvalue_ID[[set]] = corPvalueFisher(moduleTraitCor_ID[[set]], exprSize_ID$nSamples[set], twoSided = TRUE)
#}

for (set in 1:All_nSets_ID) {
  moduleTraitCor_ID[[set]] = WGCNA::cor(
    consensusMEs_ID[[set]]$data, 
    pheno_ID_pre[[set]]$data, 
    method = "pearson", 
    use = "p"
  )
  
  moduleTraitPvalue_ID[[set]] = WGCNA::corPvalueFisher(
    moduleTraitCor_ID[[set]], 
    exprSize_ID$nSamples[set], 
    twoSided = TRUE
  )
}


#=====================================================================================
#
#  Code chunk 7- Perform Module-trait correlations
#
#=====================================================================================

sizeGrWindow(12, 7)
pdf(file = "Module_Combined_Trait_Preterm_Correlation_Identified_Baseline.pdf", wi = 12, he = 7)
set = 1
textMatrix_ID =  apply(moduleTraitPvalue_ID[[set]], 2, function(x) {sapply(x, function(y) {ifelse(y < 0.05, "*", "")})})
dim(textMatrix_ID) = dim(moduleTraitCor_ID[[set]])
par(mar = c(14, 8.8, 3, 2.2))
labeledHeatmap(Matrix = moduleTraitCor_ID[[set]], xLabels = colnames(moduleTraitCor_ID[[set]]), yLabels = rownames(moduleTraitCor_ID[[set]]), 
               ySymbols = gsub("ME", "", rownames(moduleTraitCor_ID[[set]])), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = textMatrix_ID, setStdMargins = TRUE, cex.text = 2, textAdj = c(0.5, 0.8), 
               zlim = c(-1, 1), main = "Identified 10-16 weeks Consensus Module Trait Correlations Spearman Correlation", cex.lab.y = 0.9,  cex.lab.x = 0.9, plotLegend = TRUE, legendLabel = "Spearman Correlation Coefficient")
dev.off()


# Extract correlation and p-value matrices
library(openxlsx)

set = 1
cor_matrix <- moduleTraitCor_ID[[set]]
pval_matrix <- moduleTraitPvalue_ID[[set]]

# Optionally round values for readability
cor_matrix_rounded <- round(cor_matrix, 3)
pval_matrix_rounded <- signif(pval_matrix, 3)

# Create a workbook and add both matrices as separate sheets
wb <- createWorkbook()
addWorksheet(wb, "Correlations")
addWorksheet(wb, "P-values")

writeData(wb, sheet = "Correlations", cor_matrix_rounded, rowNames = TRUE)
writeData(wb, sheet = "P-values", pval_matrix_rounded, rowNames = TRUE)

# Save to file
saveWorkbook(wb, "ModuleTraitCor_Pval_Baseline.xlsx", overwrite = TRUE)

set = 2
textMatrix_ID = apply(moduleTraitPvalue_ID[[set]], 2, function(x) {sapply(x, function(y) {ifelse(y < 0.05, "*", "")})})
dim(textMatrix_ID) = dim(moduleTraitCor_ID[[set]])
sizeGrWindow(12, 7)
pdf(file = "Module_Combined_Trait_Preterm_Correlation_Identified_36-38.pdf", wi = 12, he = 7)
par(mar = c(14, 8.8, 3, 2.2))
labeledHeatmap(Matrix = moduleTraitCor_ID[[set]], xLabels = colnames(moduleTraitCor_ID[[set]]), yLabels = rownames(moduleTraitCor_ID[[set]]), 
               ySymbols = gsub("ME", "", rownames(moduleTraitCor_ID[[set]])), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = textMatrix_ID, setStdMargins = TRUE, cex.text = 2, textAdj = c(0.5, 0.8), 
               zlim = c(-1, 1), main = "Identified 36-38 weeks Consensus Module Trait Correlations Spearman Correlation", cex.lab.y = 0.9,cex.lab.x = 0.9, plotLegend = TRUE, legendLabel = "Spearman Correlation Coefficient")
dev.off()

set = 2
cor_matrix <- moduleTraitCor_ID[[set]]
pval_matrix <- moduleTraitPvalue_ID[[set]]

# Optionally round values for readability
cor_matrix_rounded <- round(cor_matrix, 3)
pval_matrix_rounded <- signif(pval_matrix, 3)

# Create a workbook and add both matrices as separate sheets
wb <- createWorkbook()
addWorksheet(wb, "Correlations")
addWorksheet(wb, "P-values")

writeData(wb, sheet = "Correlations", cor_matrix_rounded, rowNames = TRUE)
writeData(wb, sheet = "P-values", pval_matrix_rounded, rowNames = TRUE)

# Save to file
saveWorkbook(wb, "ModuleTraitCor_Pval_36-38.xlsx", overwrite = TRUE)

set = 3
textMatrix_ID = apply(moduleTraitPvalue_ID[[set]], 2, function(x) {sapply(x, function(y) {ifelse(y < 0.05, "*", "")})})
dim(textMatrix_ID) = dim(moduleTraitCor_ID[[set]])
sizeGrWindow(12, 7)
pdf(file = "Module_Combined_Trait__PretermCorrelation_Identified_Postpartum.pdf", wi = 12, he = 7)
par(mar = c(14, 8.8, 3, 2.2))
labeledHeatmap(Matrix = moduleTraitCor_ID[[set]], xLabels = colnames(moduleTraitCor_ID[[set]]), yLabels = rownames(moduleTraitCor_ID[[set]]), 
               ySymbols = gsub("ME", "", rownames(moduleTraitCor_ID[[set]])), colorLabels = FALSE, colors = blueWhiteRed(50), 
               textMatrix = textMatrix_ID, setStdMargins = TRUE, cex.text = 2, textAdj = c(0.5, 0.8), 
               zlim = c(-1, 1), main = "Identified Postpartum Consensus Module Trait Correlations Spearman Correlation", cex.lab.y = 0.9,cex.lab.x = 0.9, plotLegend = TRUE, legendLabel = "Spearman Correlation Coefficient")
dev.off()

set = 3
cor_matrix <- moduleTraitCor_ID[[set]]
pval_matrix <- moduleTraitPvalue_ID[[set]]

# Optionally round values for readability
cor_matrix_rounded <- round(cor_matrix, 3)
pval_matrix_rounded <- signif(pval_matrix, 3)

# Create a workbook and add both matrices as separate sheets
wb <- createWorkbook()
addWorksheet(wb, "Correlations")
addWorksheet(wb, "P-values")

writeData(wb, sheet = "Correlations", cor_matrix_rounded, rowNames = TRUE)
writeData(wb, sheet = "P-values", pval_matrix_rounded, rowNames = TRUE)

# Save to file
saveWorkbook(wb, "ModuleTraitCor_Pval_Postpartum.xlsx", overwrite = TRUE)

#=====================================================================================
#
#  Code chunk 8- Prepare data to create dataset that has all of the ConsensusMEs and Traits together
#
#=====================================================================================



# Prepare data to create dataset that has all of the ConsensusMEs and Traits together
ConsensusModuleEigenValues_ID <- as.data.frame(consensusMEs_ID)
ID_Combined_Traits_with_preterm 
demographic_data

library(tibble)  # Load tibble for rownames_to_column function
library(dplyr)       # For data manipulation and joining
library(writexl)     # For writing to Excel

# Convert row names to a column in each dataframe
ConsensusModuleEigenValues_ID_m <- tibble(RowNames = rownames(ConsensusModuleEigenValues_ID), ConsensusModuleEigenValues_ID)
demographic_data <- tibble(RowNames = rownames(demographic_data), demographic_data)

# Sequentially merge data frames
merged_data_ID_adverse <- ConsensusModuleEigenValues_ID_m %>%
  left_join(demographic_data, by = "RowNames") 

# sort your data by RowNames
# Sort the merged dataset by User_ID
merged_data_ID_adverse <- merged_data_ID_adverse[order(merged_data_ID_adverse$RowNames), ]
head(merged_data_ID_adverse)

# Write the merged data to an Excel file
write_xlsx(merged_data_ID_adverse, "MergedData_ID_adverse.xlsx")

colnames(merged_data_ID_adverse)
#=====================================================================================
#
#  Code chunk 9- Prepare data to create dataset that has all of the Lipids and Traits together
#
#=====================================================================================

#remove unidentified lipids
Identified_Baseline <- as.data.frame(Baseline[,-c(1,2,454:2653)])
Identified_TP36_38weeks<- as.data.frame(TP36_38weeks[,-c(1,2,454:2653)])
Identified_Postpartum <- as.data.frame(Postpartum[,-c(1,2,454:2653)])

#rename your colnames
colnames(Identified_Baseline)
colnames(Identified_Baseline) <- paste("baseline", colnames(Identified_Baseline), sep = "_")

colnames(Identified_TP36_38weeks)
colnames(Identified_TP36_38weeks) <- paste("TP36_38weeks", colnames(Identified_TP36_38weeks), sep = "_")

colnames(Identified_Postpartum)
colnames(Identified_Postpartum) <- paste("Postpartum", colnames(Identified_Postpartum), sep = "_")

rownames(Identified_Baseline)
rownames(Identified_TP36_38weeks)
rownames(Identified_Postpartum)
rownames(demographic_data)
# Convert row names to a column in each dataframe

Identified_Baseline <- tibble(RowNames = rownames(Identified_Baseline), Identified_Baseline)
Identified_TP36_38weeks <- tibble(RowNames = rownames(Identified_TP36_38weeks), Identified_TP36_38weeks)
Identified_Postpartum <- tibble(RowNames = rownames(Identified_Postpartum), Identified_Postpartum)
demographic_data <- tibble(RowNames = rownames(demographic_data), demographic_data)


# Sequentially merge data frames
merged_data_ID_lipids_traits <- Identified_Baseline %>%
  left_join(Identified_TP36_38weeks, by = "RowNames") %>%
  left_join(Identified_Postpartum, by = "RowNames") %>%
  left_join(demographic_data, by = "RowNames")

# sort your data by RowNames
# Sort the merged dataset by User_ID
merged_data_ID_lipids_traits_adverse <- merged_data_ID_lipids_traits[order(merged_data_ID_lipids_traits$RowNames), ]

# Write the merged data to an Excel file
write_xlsx(merged_data_ID_lipids_traits_adverse, "MergedData_ID_lipids_traits_adverse.xlsx")




#=====================================================================================
#
#  Code chunk - Modified module-trait correlations
#
#=====================================================================================
#for make correlation with WGCNA
#make the trait correlation first
plotMEtraitCorWGCNA <- function(cor_matrix, p_matrix = NULL, 
                                moduleOrder = NULL, traitOrder = NULL, 
                                file = NULL, width = 12, height = 7,
                                main_title = "Module-Trait Correlation Heatmap",
                                label.size = 2, star_cutoff = 0.05,
                                margins = c(14, 8.8, 3, 2.2)) {
  
  # Reorder matrix if needed
  if (!is.null(moduleOrder)) cor_matrix <- cor_matrix[moduleOrder, , drop = FALSE]
  if (!is.null(traitOrder)) cor_matrix <- cor_matrix[, traitOrder, drop = FALSE]
  
  # Prepare significance stars
  if (!is.null(p_matrix)) {
    if (!is.null(moduleOrder)) p_matrix <- p_matrix[moduleOrder, , drop = FALSE]
    if (!is.null(traitOrder)) p_matrix <- p_matrix[, traitOrder, drop = FALSE]
    
    textMatrix <- apply(p_matrix, 2, function(x) {
      sapply(x, function(p) ifelse(p < star_cutoff, "*", ""))
    })
    dim(textMatrix) <- dim(cor_matrix)
  } else {
    textMatrix <- matrix("", nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))
  }
  
  # Output to PDF if file is specified
  if (!is.null(file)) pdf(file = file, width = width, height = height)
  
  # Draw the heatmap
  par(mar = margins)
  labeledHeatmap(Matrix = cor_matrix,
                 xLabels = colnames(cor_matrix),
                 yLabels = rownames(cor_matrix),
                 ySymbols = gsub("^ME", "", rownames(cor_matrix)),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = TRUE,
                 cex.text = label.size,
                 textAdj = c(0.5, 0.8),
                 zlim = c(-1, 1),
                 main = main_title,
                 cex.lab.x = 0.9, cex.lab.y = 0.9,
                 plotLegend = TRUE,
                 legendLabel = "Pearson Correlation")
  
  if (!is.null(file)) dev.off()
}


# Timepoint 1: Baseline
plotMEtraitCorWGCNA(
  cor_matrix = moduleTraitCor_ID[[1]],
  p_matrix = moduleTraitPvalue_ID[[1]],
  file = "Module_Combined_Trait_Preterm_Correlation_Identified_Baseline.pdf",
  main_title = "Identified Baseline Consensus Module Trait Correlations",
  margins = c(14, 8.8, 3, 2.2)
)

# Timepoint 2: 36-38 weeks
plotMEtraitCorWGCNA(
  cor_matrix = moduleTraitCor_ID[[2]],
  p_matrix = moduleTraitPvalue_ID[[2]],
  file = "Module_Combined_Trait_Preterm_Correlation_Identified_36-38.pdf",
  main_title = "Identified 36-38 Weeks Consensus Module Trait Correlations",
  margins = c(14, 8.8, 3, 2.2)
)

# Timepoint 3: Postpartum
plotMEtraitCorWGCNA(
  cor_matrix = moduleTraitCor_ID[[3]],
  p_matrix = moduleTraitPvalue_ID[[3]],
  file = "Module_Combined_Trait_Preterm_Correlation_Identified_Postpartum.pdf",
  main_title = "Identified Postpartum Consensus Module Trait Correlations",
  margins = c(14, 8.8, 3, 2.2)
)

#=====================================================================================
#
#  Code chunk - Dot plots and Scatter plots for module-trait correlations
#
#=====================================================================================

#New code for plot that will add significance to the results. 
# This is only one timepoint at a time.

plotMEtraitDot <- function(moduleEigennode, 
                           trait, 
                           traitLabels = NULL, 
                           colors = NULL,
                           corr_val = NULL,
                           p_val = NULL,
                           ylim = NULL,
                           xlab = "Trait",
                           ylab = "Module Eigennode",
                           title = NULL,
                           file = NULL,
                           width = 6,
                           height = 4) {
  df <- data.frame(
    Trait = as.factor(trait),
    ME = moduleEigennode
  )
  
  if (!is.null(traitLabels)) {
    levels(df$Trait) <- names(traitLabels)
  }
  
  if (!is.null(colors)) {
    fill_colors <- colors
  } else {
    fill_colors <- scales::hue_pal()(length(levels(df$Trait)))
    names(fill_colors) <- levels(df$Trait)
  }
  
  # Append p-value stars to title
  stars <- ifelse(!is.null(p_val), 
                  ifelse(p_val < 0.01, "**", ifelse(p_val < 0.05, "*", "")), "")
  
  final_title <- paste0(title, ifelse(!is.null(p_val), paste0(" (p=", signif(p_val, 3), ")", stars), ""))
  
  p <- ggplot(df, aes(x = Trait, y = ME, fill = Trait)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
    scale_fill_manual(values = fill_colors) +
    labs(title = final_title, x = xlab, y = ylab) +
    coord_cartesian(ylim = ylim) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    )
  
  if (!is.null(file)) {
    ggsave(file, plot = p, width = width, height = height)
  } else {
    print(p)
  }
}

plotMEtraitScatter <- function(moduleEigennode, 
                               trait, 
                               corr_val = NULL,
                               p_val = NULL,
                               xlab = "Trait", 
                               ylab = "Module Eigennode", 
                               title = NULL,
                               ylim = NULL,
                               file = NULL,
                               width = 6,
                               height = 4) {
  
  df <- data.frame(
    Trait = trait,
    ME = moduleEigennode
  )
  
  # Title with correlation and p-value
  stars <- ifelse(!is.null(p_val), 
                  ifelse(p_val < 0.01, "**", ifelse(p_val < 0.05, "*", "")), "")
  
  final_title <- if (!is.null(corr_val) & !is.null(p_val)) {
    paste0(title, " (r = ", signif(corr_val, 3), 
           ", p = ", signif(p_val, 3), ")", stars)
  } else {
    title
  }
  
  p <- ggplot(df, aes(x = Trait, y = ME)) +
    geom_point(color = "#2c7fb8", size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "#d73027", linetype = "dashed") +
    labs(title = final_title, x = xlab, y = ylab) +
    coord_cartesian(ylim = ylim) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  
  if (!is.null(file)) {
    ggsave(file, plot = p, width = width, height = height)
  } else {
    print(p)
  }
}

# Get stat info
corr_val <- moduleTraitCor_ID[[1]]["MEmagenta", "apo"]
p_val <- moduleTraitPvalue_ID[[1]]["MEmagenta", "apo"]

# Dot plot
plotMEtraitDot(
  moduleEigennode = ConsensusModuleEigenValues_ID$Identified_Baseline.data.MEmagenta, 
  trait = demographic_data$apo,
  traitLabels = c("No APO" = 0, "APO" = 1),
  colors = c("No APO" = "#3366CC", "APO" = "#FF3366"),
  ylim = c(-0.4, 0.4),
  ylab = "Magenta Module Eigennode",
  title = "Magenta vs APO (Baseline)",
  corr_val = corr_val,
  p_val = p_val,
  file = "Baseline_Magenta_APO_Dotplot.pdf"
)

corr_val <- moduleTraitCor_ID[[2]]["MEpurple", "bmi"]
p_val <- moduleTraitPvalue_ID[[2]]["MEpurple", "bmi"]

plotMEtraitScatter(
  moduleEigennode = ConsensusModuleEigenValues_ID$Identified_TP36_38weeks.data.MEpurple,
  trait = demographic_data$bmi,
  xlab = "BMI",
  ylab = "Purple Module Eigennode",
  title = "Purple vs BMI (36-38 weeks)",
  corr_val = corr_val,
  p_val = p_val,
  file = "Purple_BMI_Scatter_36weeks.pdf"
)

#=====================================================================================
#
#  Code chunk - Dot plots and Scatter plots for module-trait correlations
#
#=====================================================================================


#New Code takes into account the time point and the trait of interest
# also adds the correlation coefficient and p-value

library(ggplot2)
library(patchwork)

plotMEtraitTimepointsWithStats <- function(moduleName,
                                           timepointEigennodes,
                                           traitData,
                                           corList, pvalList,
                                           traitLabels = NULL,
                                           colors = NULL,
                                           isCategorical = FALSE,
                                           ylim = NULL,
                                           xlab = "Trait",
                                           ylab = "Module Eigennode",
                                           title_prefix = "",
                                           file = NULL,
                                           width = 12,
                                           height = 5) {
  
  plot_list <- list()
  
  for (tp in names(timepointEigennodes)) {
    Eigennode <- timepointEigennodes[[tp]]
    trait <- traitData[[tp]]
    corr_val <- corList[[tp]]
    p_val <- pvalList[[tp]]
    
    stars <- ifelse(p_val < 0.01, "**", ifelse(p_val < 0.05, "*", ""))
    
    title_text <- if (isCategorical) {
      paste0(title_prefix, tp, " (p = ", signif(p_val, 3), ")", stars)
    } else {
      paste0(title_prefix, tp, " (r = ", signif(corr_val, 3), 
             ", p = ", signif(p_val, 3), ")", stars)
    }
    
    df <- data.frame(Trait = trait, ME = Eigennode)
    
    if (isCategorical) {
      df$Trait <- as.factor(df$Trait)
      if (!is.null(traitLabels)) {
        levels(df$Trait) <- names(traitLabels)
      }
      fill_colors <- if (!is.null(colors)) colors else scales::hue_pal()(length(levels(df$Trait)))
      p <- ggplot(df, aes(x = Trait, y = ME, fill = Trait)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.6) +
        geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
        scale_fill_manual(values = fill_colors)
    } else {
      p <- ggplot(df, aes(x = Trait, y = ME)) +
        geom_point(color = "#2c7fb8", size = 2, alpha = 0.7) +
        geom_smooth(method = "lm", se = TRUE, color = "#d73027", linetype = "dashed")
    }
    
    p <- p +
      labs(title = title_text, x = xlab, y = ylab) +
      coord_cartesian(ylim = ylim) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 11),
        panel.grid = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0.7, 0.5, 0.7, 0.5), "cm")
      )
    
    plot_list[[tp]] <- p
  }
  
  combined <- wrap_plots(plot_list, ncol = 3)
  
  if (!is.null(file)) {
    ggsave(file, plot = combined, width = width, height = height)
  } else {
    print(combined)
  }
}

#example usage dot plot
plotMEtraitTimepointsWithStats(
  moduleName = "magenta",
  timepointEigennodes = list(
    Baseline   = ConsensusModuleEigenValues_ID$Identified_Baseline.data.MEmagenta,
    TP36_38    = ConsensusModuleEigenValues_ID$Identified_TP36_38weeks.data.MEmagenta,
    Postpartum = ConsensusModuleEigenValues_ID$Identified_Postpartum.data.MEmagenta
  ),
  traitData = list(
    Baseline   = demographic_data$apo,
    TP36_38    = demographic_data$apo,
    Postpartum = demographic_data$apo
  ),
  corList = list(
    Baseline   = moduleTraitCor_ID[[1]]["MEmagenta", "apo"],
    TP36_38    = moduleTraitCor_ID[[2]]["MEmagenta", "apo"],
    Postpartum = moduleTraitCor_ID[[3]]["MEmagenta", "apo"]
  ),
  pvalList = list(
    Baseline   = moduleTraitPvalue_ID[[1]]["MEmagenta", "apo"],
    TP36_38    = moduleTraitPvalue_ID[[2]]["MEmagenta", "apo"],
    Postpartum = moduleTraitPvalue_ID[[3]]["MEmagenta", "apo"]
  ),
  traitLabels = c("No APO" = 0, "APO" = 1),
  colors = c("No APO" = "#3366CC", "APO" = "#FF3366"),
  isCategorical = TRUE,
  ylim = c(-0.4, 0.4),
  ylab = "Magenta Module Eigennode",
  title_prefix = "Magenta vs APO – ",
  file = "Magenta_vs_APO_3TP_Annotated.pdf"
)

#example usage scatter plot
plotMEtraitTimepointsWithStats(
  moduleName = "purple",
  timepointEigennodes = list(
    Baseline   = ConsensusModuleEigenValues_ID$Identified_Baseline.data.MEpurple,
    TP36_38    = ConsensusModuleEigenValues_ID$Identified_TP36_38weeks.data.MEpurple,
    Postpartum = ConsensusModuleEigenValues_ID$Identified_Postpartum.data.MEpurple
  ),
  traitData = list(
    Baseline   = demographic_data$bmi,
    TP36_38    = demographic_data$bmi,
    Postpartum = demographic_data$bmi
  ),
  corList = list(
    Baseline   = moduleTraitCor_ID[[1]]["MEpurple", "bmi"],
    TP36_38    = moduleTraitCor_ID[[2]]["MEpurple", "bmi"],
    Postpartum = moduleTraitCor_ID[[3]]["MEpurple", "bmi"]
  ),
  pvalList = list(
    Baseline   = moduleTraitPvalue_ID[[1]]["MEpurple", "bmi"],
    TP36_38    = moduleTraitPvalue_ID[[2]]["MEpurple", "bmi"],
    Postpartum = moduleTraitPvalue_ID[[3]]["MEpurple", "bmi"]
  ),
  isCategorical = FALSE,
  xlab = "BMI",
  ylab = "Purple Module Eigennode",
  ylim = c(-0.4, 0.4),
  title_prefix = "Purple vs BMI – ",
  file = "Purple_vs_BMI_3TP_Scatterplot.pdf"
)

#=====================================================================================
#
#  Code chunk - Dot plots and Scatter plots for module-trait correlations
#
#=====================================================================================

#This will automate everything
#You just need to select your module traits of interest
autoPlotMEtraitTimepoints <- function(modules, traits, 
                                      EigennodeDF, traitData,
                                      corMatrices, pvalMatrices,
                                      traitInfo,
                                      ylims = NULL,
                                      outdir = "ModuleTraitPlots") {
  dir.create(outdir, showWarnings = FALSE)
  
  for (module in modules) {
    for (trait in traits) {
      message("Processing ", module, " vs ", trait)
      
      isCategorical <- traitInfo[[trait]]$type == "categorical"
      traitLabels <- traitInfo[[trait]]$labels
      colors <- traitInfo[[trait]]$colors
      ylim <- ylims[[module]]
      
      # Extract ME data by column name
      MEdata <- list(
        Baseline   = EigennodeDF[[paste0("Identified_Baseline.data.ME", module)]],
        TP36_38    = EigennodeDF[[paste0("Identified_TP36_38weeks.data.ME", module)]],
        Postpartum = EigennodeDF[[paste0("Identified_Postpartum.data.ME", module)]]
      )
      
      # Trait values
      traitVals <- list(
        Baseline   = traitData$Baseline[[trait]],
        TP36_38    = traitData$TP36_38[[trait]],
        Postpartum = traitData$Postpartum[[trait]]
      )
      
      # Correlation/pval
      corList <- list(
        Baseline   = corMatrices[[1]][paste0("ME", module), trait],
        TP36_38    = corMatrices[[2]][paste0("ME", module), trait],
        Postpartum = corMatrices[[3]][paste0("ME", module), trait]
      )
      pvalList <- list(
        Baseline   = pvalMatrices[[1]][paste0("ME", module), trait],
        TP36_38    = pvalMatrices[[2]][paste0("ME", module), trait],
        Postpartum = pvalMatrices[[3]][paste0("ME", module), trait]
      )
      
      # Output file name
      file_name <- paste0(outdir, "/", module, "_vs_", trait, "_3TP.pdf")
      
      # Plot using previously defined function
      plotMEtraitTimepointsWithStats(
        moduleName = module,
        timepointEigennodes = MEdata,
        traitData = traitVals,
        corList = corList,
        pvalList = pvalList,
        traitLabels = traitLabels,
        colors = colors,
        isCategorical = isCategorical,
        ylim = ylim,
        xlab = trait,
        ylab = paste(module, "Module Eigennode"),
        title_prefix = paste0(module, " vs ", trait, " - "),
        file = file_name
      )
    }
  }
}




autoPlotMEtraitTimepoints(
  modules = c("magenta", "purple", "brown", "turquoise","greenyellow","salmon","brown","yellow"),
  traits = c("apo"#, "bmi","group","bmi_base","gwg","egwg","ppwr","ppwr_e"
             ),
  EigennodeDF = ConsensusModuleEigenValues_ID,
  traitData = traitData,
  corMatrices = moduleTraitCor_ID,
  pvalMatrices = moduleTraitPvalue_ID,
  traitInfo = traitInfo,
  ylims = ylims,
  outdir = "ModuleTrait_Timepoint_Plots"
)

#=====================================================================================
#
#  Code chunk -This is the current Code to use for your figures Dot plots and Scatter plots for module-trait correlations
#
#=====================================================================================

#if the code above fails. You can use the code below. 
#trying to add square boarder so updated the code above
#at the moment commented out letter those can be added later on 
library(ggplot2)
library(gridExtra)
library(grid)

plotMEtraitTimepointsWithStats <- function(moduleName, timepointEigennodes, traitData,
                                           corList, pvalList,
                                           traitLabels = NULL, colors = NULL,
                                           isCategorical = TRUE,
                                           ylim = NULL, xlab = "Trait", ylab = NULL,
                                           title_prefix = NULL,
                                           file = NULL, width = 11, height = 4) {
  plots <- list()
  tp_labels <- c("Baseline", "TP36_38", "Postpartum")
 # panel_labels <- c("A", "B", "C")
  
  for (i in seq_along(tp_labels)) {
    tp <- tp_labels[i]
    Eigennode <- timepointEigennodes[[tp]]
    trait <- traitData[[tp]]
    
    if (is.null(Eigennode) || is.null(trait) || length(Eigennode) != length(trait)) {
      warning(paste("Skipping", tp, "due to invalid data"))
      next
    }
    
    #df <- data.frame(Trait = trait, ME = Eigennode)
    df <- data.frame(Trait = trait, ME = Eigennode)
    df <- na.omit(df)
    
    corr_val <- corList[[tp]]
    p_val <- pvalList[[tp]]
    stars <- ifelse(p_val < 0.01, "**", ifelse(p_val < 0.05, "*", ""))
    title <- if (isCategorical) {
      paste0(title_prefix, tp, " (p = ", signif(p_val, 3), ")", stars)
    } else {
      paste0(title_prefix, tp, " (r = ", signif(corr_val, 3), ", p = ", signif(p_val, 3), ")", stars)
    }
    
    # Base plot
    p <- if (isCategorical) {
      df$Trait <- factor(df$Trait)
      if (!is.null(traitLabels)) levels(df$Trait) <- names(traitLabels)
      ggplot(df, aes(x = Trait, y = ME, fill = Trait)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.6) +
        geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
        scale_fill_manual(values = colors) +
        labs(title = title, x = xlab, y = ylab) +
        coord_cartesian(ylim = ylim) +
        theme_minimal(base_size = 11) +
        theme(
          panel.grid = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.title.position = "plot",
          plot.margin = ggplot2::margin(10, 5, 5, 5)  # top, right, bottom, left
        )
      
      
    } else {
      ggplot(df, aes(x = Trait, y = ME)) +
        geom_point(color = "#2c7fb8", size = 2, alpha = 0.7) +
        geom_smooth(method = "lm", se = TRUE, color = "#d73027", linetype = "dashed") +
        labs(title = title, x = xlab, y = ylab) +
        coord_cartesian(ylim = ylim) +
        theme_minimal(base_size = 11) +
        theme(
          panel.grid = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.title.position = "plot",
          plot.margin = ggplot2::margin(10, 5, 5, 5)  # top, right, bottom, left
        )
      
      
    }
    
    # Add panel letter # i removed the letters for now
 #   p <- p + annotation_custom(
  #    grob = textGrob(panel_labels[i], x = unit(0, "npc"), y = unit(1, "npc"),
   #                   just = c("left", "top"), gp = gpar(fontsize = 14, fontface = "bold")),
    #  xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    #)
    
    plots[[tp]] <- p
  }
  
  # Save the arranged plot
  if (length(plots) > 0 && !is.null(file)) {
    ggsave(file, plot = grid.arrange(grobs = plots, ncol = 3), width = width, height = height)
  }
}
#this fuction will automate everything but it is dependent on the helper function
#before called plotMEtraitTimepointsWithStats
autoPlotMEtraitTimepoints <- function(modules, traits, 
                                      EigennodeDF, traitData,
                                      corMatrices, pvalMatrices,
                                      traitInfo, ylims = NULL,
                                      outdir = "ModuleTraitPlots") {
  dir.create(outdir, showWarnings = FALSE)
  
  for (module in modules) {
    for (trait in traits) {
      isCategorical <- traitInfo[[trait]]$type == "categorical"
      traitLabels <- traitInfo[[trait]]$labels
      colors <- traitInfo[[trait]]$colors
      ylim <- if (!is.null(ylims)) ylims[[module]] else NULL
      
      # Extract Eigennodes using correct names
      MEdata <- list(
        Baseline   = EigennodeDF[[paste0("Identified_Baseline.data.ME", module)]],
        TP36_38    = EigennodeDF[[paste0("Identified_TP36_38weeks.data.ME", module)]],
        Postpartum = EigennodeDF[[paste0("Identified_Postpartum.data.ME", module)]]
      )
      
      # Trait values
      traitVals <- list(
        Baseline   = traitData$Baseline[[trait]],
        TP36_38    = traitData$TP36_38[[trait]],
        Postpartum = traitData$Postpartum[[trait]]
      )
      
      # Correlations & p-values
      corList <- list(
        Baseline   = corMatrices[[1]][paste0("ME", module), trait],
        TP36_38    = corMatrices[[2]][paste0("ME", module), trait],
        Postpartum = corMatrices[[3]][paste0("ME", module), trait]
      )
      pvalList <- list(
        Baseline   = pvalMatrices[[1]][paste0("ME", module), trait],
        TP36_38    = pvalMatrices[[2]][paste0("ME", module), trait],
        Postpartum = pvalMatrices[[3]][paste0("ME", module), trait]
      )
      
      file_name <- paste0(outdir, "/", module, "_vs_", trait, "_3TP.pdf")
      
      message("Plotting ", module, " vs ", trait)
      plotMEtraitTimepointsWithStats(
        moduleName = module,
        timepointEigennodes = MEdata,
        traitData = traitVals,
        corList = corList,
        pvalList = pvalList,
        traitLabels = traitLabels,
        colors = colors,
        isCategorical = isCategorical,
        ylim = ylim,
        xlab = trait,
        ylab = paste(module, "Module Eigennode"),
        title_prefix = paste0(module, " vs ", trait, " – "),
        file = file_name
      )
    }
  }
}


#usage examples i.e simple and straight forward
autoPlotMEtraitTimepoints(
  modules = c("magenta", "purple", "brown","tan","grey60","red","blue","lightcyan","lightyellow","cyan","black","midnightblue"),
  traits = c("apo", "bmi","ppwr_e","weight_base","weight_wk36","weight_mth3"), #"ppwr",
  EigennodeDF = ConsensusModuleEigenValues_ID,
  traitData = traitData,
  corMatrices = moduleTraitCor_ID,
  pvalMatrices = moduleTraitPvalue_ID,
  traitInfo = traitInfo,
  ylims = ylims,
  outdir = "ModuleTrait_Timepoint_Plots"
)

#usage examples in case you want to adjust some thing
modules <- c("magenta", "purple", "brown", "turquoise","greenyellow","salmon","brown","yellow")
#modules <- c("magenta", "purple", "brown")
traits <- c("apo")#, "bmi")

traitInfo <- list(
  apo = list(type = "categorical", 
             labels = c("No APO" = 0, "APO" = 1), 
             colors = c("No APO" = "#3366CC", "APO" = "#FF3366")),
  bmi = list(type = "continuous")
)

traitData <- list(
  Baseline = demographic_data,
  TP36_38 = demographic_data,
  Postpartum = demographic_data
)

ylims <- list(
  magenta = c(-0.4, 0.4),
  purple = c(-0.4, 0.4),
  brown = c(-0.4, 0.4)
)

autoPlotMEtraitTimepoints(
  modules = modules,
  traits = traits,
  EigennodeDF = ConsensusModuleEigenValues_ID,
  traitData = traitData,
  corMatrices = moduleTraitCor_ID,
  pvalMatrices = moduleTraitPvalue_ID,
  traitInfo = traitInfo,
  ylims = ylims,
  outdir = "ModuleTrait_Timepoint_Plots"
)

#tan and group
modules <- c("tan")
traits <- c("group")

traitInfo <- list(
  group = list(type = "categorical", 
             labels = c("Control" = 0, "Intervention" = 1), 
             colors = c("Control" = "#3366CC", "Intervention" = "#FF3366"))
)

traitData <- list(
  Baseline = demographic_data,
  TP36_38 = demographic_data,
  Postpartum = demographic_data
)

ylims <- list(
  magenta = c(-0.4, 0.6),
  purple = c(-0.4, 0.6),
  brown = c(-0.4, 0.6)
)

autoPlotMEtraitTimepoints(
  modules = modules,
  traits = traits,
  EigennodeDF = ConsensusModuleEigenValues_ID,
  traitData = traitData,
  corMatrices = moduleTraitCor_ID,
  pvalMatrices = moduleTraitPvalue_ID,
  traitInfo = traitInfo,
  ylims = ylims,
  outdir = "ModuleTrait_Timepoint_Plots"
)

#=====================================================================================
#
#  Code chunk -This is the current Code to use for PPWR_E or other variables with NAs for
#   your figures Dot plots and Scatter plots
#   for module-trait correlations. This has no letter labels for panels
#
#=====================================================================================
plotMEtraitTimepointsWithStats <- function(moduleName, timepointEigennodes, traitData,
                                           corList, pvalList,
                                           traitLabels = NULL, colors = NULL,
                                           isCategorical = TRUE,
                                           ylim = NULL, xlab = "Trait", ylab = NULL,
                                           title_prefix = NULL,
                                           file = NULL, width = 11, height = 4) {
  plots <- list()
  tp_labels <- c("Baseline", "TP36_38", "Postpartum")
  
  for (i in seq_along(tp_labels)) {
    tp <- tp_labels[i]
    Eigennode <- timepointEigennodes[[tp]]
    trait <- traitData[[tp]]
    
    if (is.null(Eigennode) || is.null(trait) || length(Eigennode) != length(trait)) {
      warning(paste("Skipping", tp, "due to invalid data"))
      next
    }
    
    # Combine and filter NAs
    df <- data.frame(Trait = trait, ME = Eigennode)
    df <- na.omit(df)
    
    if (nrow(df) == 0) {
      warning(paste("Skipping", tp, "due to all NA data"))
      next
    }
    
    corr_val <- corList[[tp]]
    p_val <- pvalList[[tp]]
    stars <- ifelse(p_val < 0.01, "**", ifelse(p_val < 0.05, "*", ""))
    title <- if (isCategorical) {
      paste0(title_prefix, tp, " (p = ", signif(p_val, 3), ")", stars)
    } else {
      paste0(title_prefix, tp, " (r = ", signif(corr_val, 3), ", p = ", signif(p_val, 3), ")", stars)
    }
    
    # Generate plot
    p <- if (isCategorical) {
      df$Trait <- factor(df$Trait)
      if (!is.null(traitLabels)) levels(df$Trait) <- names(traitLabels)
      ggplot(df, aes(x = Trait, y = ME, fill = Trait)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.6) +
        geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
        scale_fill_manual(values = colors) +
        labs(title = title, x = xlab, y = ylab) +
        coord_cartesian(ylim = ylim) +
        theme_minimal(base_size = 11) +
        theme(
          panel.grid = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.title.position = "plot",
          plot.margin = ggplot2::margin(10, 5, 5, 5)
        )
    } else {
      ggplot(df, aes(x = Trait, y = ME)) +
        geom_point(color = "#2c7fb8", size = 2, alpha = 0.7) +
        geom_smooth(method = "lm", se = TRUE, color = "#d73027", linetype = "dashed") +
        labs(title = title, x = xlab, y = ylab) +
        coord_cartesian(ylim = ylim) +
        theme_minimal(base_size = 11) +
        theme(
          panel.grid = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.title.position = "plot",
          plot.margin = ggplot2::margin(10, 5, 5, 5)
        )
    }
    
    plots[[tp]] <- p
  }
  
  # Save arranged plots
  if (length(plots) > 0 && !is.null(file)) {
    ggsave(file, plot = gridExtra::grid.arrange(grobs = plots, ncol = 3), width = width, height = height)
  }
}


#this fuction will automate everything but it is dependent on the helper function
#before called plotMEtraitTimepointsWithStats
autoPlotMEtraitTimepoints <- function(modules, traits, 
                                      EigennodeDF, traitData,
                                      corMatrices, pvalMatrices,
                                      traitInfo, ylims = NULL,
                                      outdir = "ModuleTraitPlots") {
  dir.create(outdir, showWarnings = FALSE)
  
  for (module in modules) {
    for (trait in traits) {
      isCategorical <- traitInfo[[trait]]$type == "categorical"
      traitLabels <- traitInfo[[trait]]$labels
      colors <- traitInfo[[trait]]$colors
      ylim <- if (!is.null(ylims)) ylims[[module]] else NULL
      
      # Raw eigennodes for each timepoint
      MEdata <- list(
        Baseline   = EigennodeDF[[paste0("Identified_Baseline.data.ME", module)]],
        TP36_38    = EigennodeDF[[paste0("Identified_TP36_38weeks.data.ME", module)]],
        Postpartum = EigennodeDF[[paste0("Identified_Postpartum.data.ME", module)]]
      )
      
      # Raw trait values for each timepoint
      traitVals_raw <- list(
        Baseline   = traitData$Baseline[[trait]],
        TP36_38    = traitData$TP36_38[[trait]],
        Postpartum = traitData$Postpartum[[trait]]
      )
      
      # Filter NA values and align Eigennodes with trait data
      traitVals <- list()
      MEdata_filtered <- list()
      
      for (tp in c("Baseline", "TP36_38", "Postpartum")) {
        trait_tp <- traitVals_raw[[tp]]
        me_tp <- MEdata[[tp]]
        
        if (!is.null(trait_tp) && !is.null(me_tp) && length(trait_tp) == length(me_tp)) {
          valid_idx <- complete.cases(trait_tp, me_tp)
          traitVals[[tp]] <- trait_tp[valid_idx]
          MEdata_filtered[[tp]] <- me_tp[valid_idx]
        } else {
          traitVals[[tp]] <- NULL
          MEdata_filtered[[tp]] <- NULL
        }
      }
      
      # Correlations & p-values
      corList <- list(
        Baseline   = corMatrices[[1]][paste0("ME", module), trait],
        TP36_38    = corMatrices[[2]][paste0("ME", module), trait],
        Postpartum = corMatrices[[3]][paste0("ME", module), trait]
      )
      pvalList <- list(
        Baseline   = pvalMatrices[[1]][paste0("ME", module), trait],
        TP36_38    = pvalMatrices[[2]][paste0("ME", module), trait],
        Postpartum = pvalMatrices[[3]][paste0("ME", module), trait]
      )
      
      file_name <- paste0(outdir, "/", module, "_vs_", trait, "_3TP.pdf")
      
      message("Plotting ", module, " vs ", trait)
      plotMEtraitTimepointsWithStats(
        moduleName = module,
        timepointEigennodes = MEdata_filtered,
        traitData = traitVals,
        corList = corList,
        pvalList = pvalList,
        traitLabels = traitLabels,
        colors = colors,
        isCategorical = isCategorical,
        ylim = ylim,
        xlab = trait,
        ylab = paste(module, "Module Eigennode"),
        title_prefix = paste0(module, " vs ", trait, " – "),
        file = file_name
      )
    }
  }
}



#ppwr_e and modules for the manuscript
modules = c("magenta", "purple", "brown","tan","grey60","red","blue","lightcyan","lightyellow","cyan","black","midnightblue")
traits <- c("ppwr_e")
traits <- c("ppwr","ppwr_e","weight_base","weight_wk36","weight_mth3","apo_hdp","apo_gdm")

traitInfo <- list(
  ppwr_e = list(type = "categorical", 
                labels = c("No PPWR" = 0, "PPWR" = 1), 
                colors = c("No PPWR" = "#3366CC", "PPWR" = "#FF3366")),
  apo_hdp = list(type = "categorical", 
                labels = c("No HDP" = 0, "HDP" = 1), 
                colors = c("No HDP" = "#3366CC", "HDP" = "#FF3366")),
  apo_gdm = list(type = "categorical", 
                 labels = c("No GDM" = 0, "GDM" = 1), 
                 colors = c("No GDM" = "#3366CC", "GDM" = "#FF3366")),
  
  ppwr = list(type = "continuous"),
  weight_base = list(type = "continuous"),
  weight_wk36 = list(type = "continuous"),
  weight_mth3 = list(type = "continuous")
)

traitData <- list(
  Baseline = demographic_data,
  TP36_38 = demographic_data,
  Postpartum = demographic_data
)

ylims <- list(
  magenta = c(-0.4, 0.4),
  purple = c(-0.4, 0.4),
  tan = c(-0.4, 0.4),
  grey60 = c(-0.4, 0.4),
  red = c(-0.4, 0.4),
  blue = c(-0.4, 0.4),
  lightcyan = c(-0.4, 0.4),
  lightyellow = c(-0.4, 0.4),
  cyan = c(-0.4, 0.4),
  black = c(-0.4, 0.4),
  midnightblue = c(-0.4, 0.4),
  brown = c(-0.4, 0.4)
)

autoPlotMEtraitTimepoints(
  modules = modules,
  traits = traits,
  EigennodeDF = ConsensusModuleEigenValues_ID,
  traitData = traitData,
  corMatrices = moduleTraitCor_ID,
  pvalMatrices = moduleTraitPvalue_ID,
  traitInfo = traitInfo,
  ylims = ylims,
  outdir = "ModuleTrait_Timepoint_Plots"
)





















#=====================================================================================
#
#  Code chunk - Temporary code_subject to removal 
#             - Dot plots and Scatter plots for module-trait correlations 
#
#=====================================================================================

#This code does take into account having multiple timepoints.
library(ggplot2)
library(patchwork)

plotMEtraitTimepoints <- function(moduleName,
                                  timepointEigennodes,# Named list: e.g. list(Baseline = ..., TP36_38 = ..., Postpartum = ...)
                                  traitData,          # Named list: same structure
                                  traitLabels = NULL, # Only for dot plots
                                  colors = NULL,
                                  isCategorical = FALSE,
                                  ylim = NULL,
                                  xlab = "Trait",
                                  ylab = "Module Eigennode",
                                  title_prefix = "",
                                  file = NULL,
                                  width = 12,
                                  height = 5) {
  
  plot_list <- list()
  
  for (tp in names(timepointEigennodes)) {
    df <- data.frame(
      Trait = traitData[[tp]],
      ME = timepointEigennodes[[tp]]
    )
    
    title_text <- paste0(title_prefix, tp)
    
    if (isCategorical) {
      df$Trait <- as.factor(df$Trait)
      if (!is.null(traitLabels)) {
        levels(df$Trait) <- names(traitLabels)
      }
      p <- ggplot(df, aes(x = Trait, y = ME, fill = Trait)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.6) +
        geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
        scale_fill_manual(values = colors) +
        labs(title = title_text, x = xlab, y = ylab)
    } else {
      p <- ggplot(df, aes(x = Trait, y = ME)) +
        geom_point(color = "#2c7fb8", size = 2, alpha = 0.7) +
        geom_smooth(method = "lm", se = TRUE, color = "#d73027", linetype = "dashed") +
        labs(title = title_text, x = xlab, y = ylab)
    }
    
    p <- p +
      coord_cartesian(ylim = ylim) +
      theme_bw(base_size = 11) +
      theme(
        plot.title = element_text(size = 11, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0.7, 0.5, 0.7, 0.5), "cm")
      )
    
    plot_list[[tp]] <- p
  }
  
  combined <- wrap_plots(plot_list, ncol = 3) +
    plot_layout(guides = "collect") &
    theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
  
  if (!is.null(file)) {
    ggsave(file, plot = combined, width = width, height = height)
  } else {
    print(combined)
  }
}

plotMEtraitTimepoints(
  moduleName = "magenta",
  timepointEigennodes = list(
    Baseline = ConsensusModuleEigenValues_ID$Identified_Baseline.data.MEmagenta,
    TP36_38  = ConsensusModuleEigenValues_ID$Identified_TP36_38weeks.data.MEmagenta,
    Postpartum = ConsensusModuleEigenValues_ID$Identified_Postpartum.data.MEmagenta
  ),
  traitData = list(
    Baseline = demographic_data$apo,
    TP36_38  = demographic_data$apo,
    Postpartum = demographic_data$apo
  ),
  traitLabels = c("No APO" = 0, "APO" = 1),
  colors = c("No APO" = "#3366CC", "APO" = "#FF3366"),
  isCategorical = TRUE,
  ylim = c(-0.4, 0.4),
  ylab = "Magenta Module Eigennode",
  title_prefix = "APO Effect on Magenta Module - ",
  file = "APO_vs_Magenta_3TP.pdf",
  height = 7 # slightly taller for room
)


#examples of usage
plotMEtraitTimepoints(
  moduleName = "magenta",
  timepointEigennodes = list(
    Baseline = ConsensusModuleEigenValues_ID$Identified_Baseline.data.MEmagenta,
    TP36_38  = ConsensusModuleEigenValues_ID$Identified_TP36_38weeks.data.MEmagenta,
    Postpartum = ConsensusModuleEigenValues_ID$Identified_Postpartum.data.MEmagenta
  ),
  traitData = list(
    Baseline = demographic_data$apo,
    TP36_38  = demographic_data$apo,
    Postpartum = demographic_data$apo
  ),
  traitLabels = c("No APO" = 0, "APO" = 1),
  colors = c("No APO" = "#3366CC", "APO" = "#FF3366"),
  isCategorical = TRUE,
  ylim = c(-0.4, 0.4),
  ylab = "Magenta Module Eigennode",
  title_prefix = "Magenta Module vs APO - ",
  file = "Magenta_vs_APO_3TP_Dotplot.pdf"
)

#scatter plot usage
plotMEtraitTimepoints(
  moduleName = "blue",
  timepointEigennodes = list(
    Baseline = ConsensusModuleEigenValues_ID$Identified_Baseline.data.MEblue,
    TP36_38  = ConsensusModuleEigenValues_ID$Identified_TP36_38weeks.data.MEblue,
    Postpartum = ConsensusModuleEigenValues_ID$Identified_Postpartum.data.MEblue
  ),
  traitData = list(
    Baseline = demographic_data$bmi,
    TP36_38  = demographic_data$bmi,
    Postpartum = demographic_data$bmi
  ),
  isCategorical = FALSE,
  ylim = c(-0.8, 0.8),
  xlab = "BMI",
  ylab = "Blue Module Eigennode",
  title_prefix = "Blue Module vs BMI - ",
  file = "Blue_vs_BMI_3TP_Scatterplot.pdf"
)



#New
#This code doesnt specify whethere there are multiple timepoints. So it will just
#take the time point you have and the ME by trait. See Usage below
plotMEtraitDot <- function(moduleEigennode, 
                           trait, 
                           traitLabels = NULL, 
                           colors = NULL,
                           ylim = NULL,
                           xlab = "Trait",
                           ylab = "Module Eigennode",
                           title = NULL,
                           file = NULL,
                           width = 6,
                           height = 4) {
  df <- data.frame(
    Trait = as.factor(trait),
    ME = moduleEigennode
  )
  
  # Assign levels and colors if provided
  if (!is.null(traitLabels)) {
    levels(df$Trait) <- names(traitLabels)
  }
  if (!is.null(colors)) {
    fill_colors <- colors
  } else {
    fill_colors <- scales::hue_pal()(length(levels(df$Trait)))
    names(fill_colors) <- levels(df$Trait)
  }
  
  p <- ggplot(df, aes(x = Trait, y = ME, fill = Trait)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
    scale_fill_manual(values = fill_colors) +
    labs(title = title, x = xlab, y = ylab) +
    coord_cartesian(ylim = ylim) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")  # top, right, bottom, left
    )
  
  if (!is.null(file)) {
    ggsave(file, plot = p, width = width, height = height)
  } else {
    print(p)
  }
}


plotMEtraitScatter <- function(moduleEigennode, 
                               trait, 
                               xlab = "Trait", 
                               ylab = "Module Eigennode", 
                               title = NULL,
                               ylim = NULL,
                               file = NULL,
                               width = 6,
                               height = 4) {
  df <- data.frame(
    Trait = trait,
    ME = moduleEigennode
  )
  
  p <- ggplot(df, aes(x = Trait, y = ME)) +
    geom_point(color = "#2c7fb8", size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "#d73027", linetype = "dashed") +
    labs(title = title, x = xlab, y = ylab) +
    coord_cartesian(ylim = ylim) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    )
  
  if (!is.null(file)) {
    ggsave(file, plot = p, width = width, height = height)
  } else {
    print(p)
  }
}

# Dotplot: Magenta module vs APO
plotMEtraitDot(
  moduleEigennode = ConsensusModuleEigenValues_ID$Identified_Baseline.data.MEmagenta, 
  trait = demographic_data$apo,
  traitLabels = c("No APO" = 0, "APO" = 1),
  colors = c("No APO" = "#3366CC", "APO" = "#FF3366"),
  ylim = c(-0.4, 0.4),
  ylab = "Magenta Module Eigennode",
  title = "Baseline Magenta Module vs APO",
  file = "Baseline_Magenta_No_vs_APO_Dotplot.pdf"
)

# Scatterplot: Magenta module vs BMI
plotMEtraitScatter(
  moduleEigennode = ConsensusModuleEigenValues_ID$Identified_TP36_38weeks.data.MEpurple,
  trait = demographic_data$bmi,
  xlab = "BMI",
  ylab = "Purple Module Eigennode",
  title = "Purple Module vs BMI (36-38 weeks)",
  file = "Purple_Module_vs_BMI_Scatterplot.pdf"
)



#this is semi-automated
#This has no panel to it but it is side by side
library(ggplot2)
library(gridExtra)
library(grid)

plotMEtraitTimepointsWithStats <- function(moduleName, timepointEigennodes, traitData,
                                           corList, pvalList,
                                           traitLabels = NULL, colors = NULL,
                                           isCategorical = TRUE,
                                           ylim = NULL, xlab = "Trait", ylab = NULL,
                                           title_prefix = NULL,
                                           file = NULL, width = 11, height = 4) {
  plots <- list()
  tp_labels <- c("Baseline", "TP36_38", "Postpartum")
  panel_labels <- c("A", "B", "C")
  
  for (i in seq_along(tp_labels)) {
    tp <- tp_labels[i]
    Eigennode <- timepointEigennodes[[tp]]
    trait <- traitData[[tp]]
    
    if (is.null(Eigennode) || is.null(trait) || length(Eigennode) != length(trait)) {
      warning(paste("Skipping", tp, "due to invalid data"))
      next
    }
    
    df <- data.frame(Trait = trait, ME = Eigennode)
    corr_val <- corList[[tp]]
    p_val <- pvalList[[tp]]
    stars <- ifelse(p_val < 0.01, "**", ifelse(p_val < 0.05, "*", ""))
    title <- if (isCategorical) {
      paste0(title_prefix, tp, " (p = ", signif(p_val, 3), ")", stars)
    } else {
      paste0(title_prefix, tp, " (r = ", signif(corr_val, 3), ", p = ", signif(p_val, 3), ")", stars)
    }
    
    p <- if (isCategorical) {
      df$Trait <- factor(df$Trait)
      if (!is.null(traitLabels)) levels(df$Trait) <- names(traitLabels)
      ggplot(df, aes(x = Trait, y = ME, fill = Trait)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.6) +
        geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
        scale_fill_manual(values = colors) +
        labs(title = title, x = xlab, y = ylab) +
        coord_cartesian(ylim = ylim) +
        theme_minimal(base_size = 11) +
        theme(panel.grid = element_blank(),
              legend.position = "none",
              plot.title = element_text(size = 10, face = "bold"))
    } else {
      ggplot(df, aes(x = Trait, y = ME)) +
        geom_point(color = "#2c7fb8", size = 2, alpha = 0.7) +
        geom_smooth(method = "lm", se = TRUE, color = "#d73027", linetype = "dashed") +
        labs(title = title, x = xlab, y = ylab) +
        coord_cartesian(ylim = ylim) +
        theme_minimal(base_size = 11) +
        theme(panel.grid = element_blank(),
              plot.title = element_text(size = 10, face = "bold"))
    }
    
    # Add panel label (A, B, C)
    p <- p + annotation_custom(
      grob = textGrob(panel_labels[i], x = unit(0, "npc"), y = unit(1, "npc"),
                      just = c("left", "top"), gp = gpar(fontsize = 14, fontface = "bold")),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    )
    
    plots[[tp]] <- p
  }
  
  if (length(plots) > 0 && !is.null(file)) {
    ggsave(file, plot = grid.arrange(grobs = plots, ncol = 3), width = width, height = height)
  }
}

#now trying side by side with automation #this is one on top of the other 

autoPlotMEtraitTimepoints <- function(modules, traits, 
                                      EigennodeDF, traitData,
                                      corMatrices, pvalMatrices,
                                      traitInfo, # named list: trait_name → list(type, labels/colors if categorical)
                                      ylims = NULL,
                                      outdir = "ModuleTraitPlots") {
  dir.create(outdir, showWarnings = FALSE)
  
  for (module in modules) {
    for (trait in traits) {
      message("Processing ", module, " vs ", trait)
      isCategorical <- traitInfo[[trait]]$type == "categorical"
      traitLabels <- traitInfo[[trait]]$labels
      colors <- traitInfo[[trait]]$colors
      ylim <- ylims[[module]]
      
      # Extract Eigennodes using unified column names
      MEdata <- list(
        Baseline   = EigennodeDF[[paste0("Identified_Baseline.data.ME", module)]],
        TP36_38    = EigennodeDF[[paste0("Identified_TP36_38weeks.data.ME", module)]],
        Postpartum = EigennodeDF[[paste0("Identified_Postpartum.data.ME", module)]]
      )
      
      # Trait values
      traitVals <- list(
        Baseline   = traitData$Baseline[[trait]],
        TP36_38    = traitData$TP36_38[[trait]],
        Postpartum = traitData$Postpartum[[trait]]
      )
      
      # Correlations & p-values
      corList <- list(
        Baseline   = corMatrices[[1]][paste0("ME", module), trait],
        TP36_38    = corMatrices[[2]][paste0("ME", module), trait],
        Postpartum = corMatrices[[3]][paste0("ME", module), trait]
      )
      pvalList <- list(
        Baseline   = pvalMatrices[[1]][paste0("ME", module), trait],
        TP36_38    = pvalMatrices[[2]][paste0("ME", module), trait],
        Postpartum = pvalMatrices[[3]][paste0("ME", module), trait]
      )
      
      # Combine plots side-by-side
      plots <- list()
      timepoints <- names(MEdata)
      for (tp in timepoints) {
        Eigennode <- MEdata[[tp]]
        trait_vec <- traitVals[[tp]]
        
        if (is.null(Eigennode) || is.null(trait_vec)) next
        if (length(Eigennode) != length(trait_vec)) next
        
        corr_val <- corList[[tp]]
        p_val <- pvalList[[tp]]
        stars <- ifelse(p_val < 0.01, "**", ifelse(p_val < 0.05, "*", ""))
        
        title_text <- if (isCategorical) {
          paste0(tp, " (p = ", signif(p_val, 3), ")", stars)
        } else {
          paste0(tp, " (r = ", signif(corr_val, 3), ", p = ", signif(p_val, 3), ")", stars)
        }
        
        df <- data.frame(Trait = trait_vec, ME = Eigennode)
        
        if (isCategorical) {
          df$Trait <- as.factor(df$Trait)
          if (!is.null(traitLabels)) levels(df$Trait) <- names(traitLabels)
          
          p <- ggplot(df, aes(x = Trait, y = ME, fill = Trait)) +
            geom_boxplot(outlier.shape = NA, alpha = 0.6) +
            geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
            scale_fill_manual(values = colors) +
            labs(title = title_text, x = trait, y = paste(module, "Module Eigennode")) +
            coord_cartesian(ylim = ylim) +
            theme_minimal(base_size = 11) +
            theme(
              legend.position = "none",
              panel.grid = element_blank(),
              plot.background = element_rect(fill = "white", color = NA)
            )
        } else {
          p <- ggplot(df, aes(x = Trait, y = ME)) +
            geom_point(color = "#2c7fb8", size = 2, alpha = 0.7) +
            geom_smooth(method = "lm", se = TRUE, color = "#d73027", linetype = "dashed") +
            labs(title = title_text, x = trait, y = paste(module, "Module Eigennode")) +
            coord_cartesian(ylim = ylim) +
            theme_minimal(base_size = 11) +
            theme(
              panel.grid = element_blank(),
              plot.background = element_rect(fill = "white", color = NA)
            )
        }
        
        plots[[tp]] <- p
      }
      
      if (length(plots) > 0) {
        file_name <- paste0(outdir, "/", module, "_vs_", trait, "_3TP.pdf")
        ggsave(file_name, plot = patchwork::wrap_plots(plots, nrow = 1), width = 12, height = 4)
      }
    }
  }
}

autoPlotMEtraitTimepoints(
  modules = c("magenta", "purple", "brown"),
  traits = c("apo", "bmi"),
  EigennodeDF = ConsensusModuleEigenValues_ID,
  traitData = traitData,
  corMatrices = moduleTraitCor_ID,
  pvalMatrices = moduleTraitPvalue_ID,
  traitInfo = traitInfo,
  ylims = ylims,
  outdir = "ModuleTrait_Timepoint_Plots"
)

#=====================================================================================
#
#  Code chunk 10- Prepare data to create dataset that has all of the hubLipids and Traits together
#
#=====================================================================================

#convert hublipids to dataframe so that you can see the colnames
hubLipids_Base_ID <- as.data.frame(hubLipids_Base_ID)
hubLipids_Week36_ID<- as.data.frame(hubLipids_Week36_ID)
hubLipids_Month3_ID <- as.data.frame(hubLipids_Month3_ID)

# Convert row names to a column in each dataframe
hubLipids_Base_ID <- tibble(RowNames = rownames(hubLipids_Base_ID), hubLipids_Base_ID)
hubLipids_Week36_ID <- tibble(RowNames = rownames(hubLipids_Week36_ID), hubLipids_Week36_ID)
hubLipids_Month3_ID <- tibble(RowNames = rownames(hubLipids_Month3_ID), hubLipids_Month3_ID)

colnames(hubLipids_Base_ID)
colnames(hubLipids_Week36_ID)
colnames(hubLipids_Month3_ID)

#make new dataframe for merging
hubLipids_Base_ID_m <- hubLipids_Base_ID
hubLipids_Week36_ID_m <- hubLipids_Week36_ID
hubLipids_Month3_ID_m <- hubLipids_Month3_ID
 
head(hubLipids_Base_ID_m) 
#add timepoints to colnames
colnames(hubLipids_Base_ID_m)[2:4] <- paste0(colnames(hubLipids_Base_ID)[2:4], "_baseline")
colnames(hubLipids_Base_ID_m)
colnames(hubLipids_Week36_ID_m)[2:4] <- paste0(colnames(hubLipids_Week36_ID)[2:4], "_TP336_38weeks")
colnames(hubLipids_Week36_ID_m)
colnames(hubLipids_Month3_ID_m)[2:4] <- paste0(colnames(hubLipids_Month3_ID)[2:4], "_Postpartum")
colnames(hubLipids_Month3_ID_m)


# Sequentially merge data frames
merged_data_ID_hub_lipids <- hubLipids_Base_ID_m %>%
  left_join(hubLipids_Week36_ID_m, by = "RowNames") %>%
  left_join(hubLipids_Month3_ID_m, by = "RowNames") 

# sort your data by RowNames if data has User_ID rownames otherwise leave it like that
#merged_data_ID_hub_lipids <- merged_data_ID_hub_lipids[order(merged_data_ID_hub_lipids$RowNames), ]
head(merged_data_ID_hub_lipids)

# Write the merged data to an Excel file
write_xlsx(merged_data_ID_hub_lipids, "MergedData_ID_hublipids_names.xlsx")


# Combine all three vectors into one
combined_lipids <- c(merged_data_ID_hub_lipids$HubLipid_baseline,
                     merged_data_ID_hub_lipids$HubLipid_TP336_38weeks,
                     merged_data_ID_hub_lipids$HubLipid_Postpartum)

# Remove duplicates to get unique lipid names
unique_lipids <- unique(combined_lipids)


# Subset the baseline dataset based on these unique lipid names
#check it worked
class(identified_lipids_baseline_data_normalized_transposed_ana[[1]])  # Should be "data.frame"
class(identified_lipids_TP36_data_normalized_transposed_ana[[1]])  # Should be "data.frame"
class(identified_postpartum_data_normalized_transposed_ana[[1]])  # Should be "data.frame"
str(ID_Updated_Combined_Traits_apo_analysis)

# get ready to subset
subset_baseline_hub_ID <- identified_lipids_baseline_data_normalized_transposed_ana[, colnames(identified_lipids_baseline_data_normalized_transposed_ana) %in% unique_lipids]
subset_TP36_38weeks_hub_ID <- identified_lipids_TP36_data_normalized_transposed_ana[, colnames(identified_lipids_TP36_data_normalized_transposed_ana) %in% unique_lipids]
subset_month3_hub_ID <- identified_postpartum_data_normalized_transposed_ana[, colnames(identified_postpartum_data_normalized_transposed_ana) %in% unique_lipids]

#add cholesterol and others 
subset_baseline_hub_ID$Cholesterol <- identified_lipids_baseline_data_normalized_transposed_ana$Cholesterol
subset_TP36_38weeks_hub_ID$Cholesterol <- identified_lipids_TP36_data_normalized_transposed_ana$Cholesterol
subset_month3_hub_ID$Cholesterol <- identified_postpartum_data_normalized_transposed_ana$Cholesterol

#FAFHA
subset_baseline_hub_ID$FAHFA.18.0.20.2 <- identified_lipids_baseline_data_normalized_transposed_ana$FAHFA.18.0.20.2
subset_TP36_38weeks_hub_ID$FAHFA.18.0.20.2 <- identified_lipids_TP36_data_normalized_transposed_ana$FAHFA.18.0.20.2
subset_month3_hub_ID$FAHFA.18.0.20.2 <- identified_postpartum_data_normalized_transposed_ana$FAHFA.18.0.20.2

subset_baseline_hub_ID$FAHFA.38.4 <- identified_lipids_baseline_data_normalized_transposed_ana$FAHFA.38.4
subset_TP36_38weeks_hub_ID$FAHFA.38.4 <- identified_lipids_TP36_data_normalized_transposed_ana$FAHFA.38.4
subset_month3_hub_ID$FAHFA.38.4 <- identified_postpartum_data_normalized_transposed_ana$FAHFA.38.4

#methyl_oxopentoic
subset_baseline_hub_ID$X4.methyl.2.oxopentanoic.acid <- identified_lipids_baseline_data_normalized_transposed_ana$X4.methyl.2.oxopentanoic.acid
subset_TP36_38weeks_hub_ID$X4.methyl.2.oxopentanoic.acid <- identified_lipids_TP36_data_normalized_transposed_ana$X4.methyl.2.oxopentanoic.acid
subset_month3_hub_ID$X4.methyl.2.oxopentanoic.acid <- identified_postpartum_data_normalized_transposed_ana$X4.methyl.2.oxopentanoic.acid

#Oleic_acid
subset_baseline_hub_ID$FA.18.1..oleic.acid. <- identified_lipids_baseline_data_normalized_transposed_ana$FA.18.1..oleic.acid.
subset_TP36_38weeks_hub_ID$FA.18.1..oleic.acid. <- identified_lipids_TP36_data_normalized_transposed_ana$FA.18.1..oleic.acid.
subset_month3_hub_ID$FA.18.1..oleic.acid. <- identified_postpartum_data_normalized_transposed_ana$FA.18.1..oleic.acid.

#linoleic_acid
subset_baseline_hub_ID$FA.18.2..linoleic.acid. <- identified_lipids_baseline_data_normalized_transposed_ana$FA.18.2..linoleic.acid.
subset_TP36_38weeks_hub_ID$FA.18.2..linoleic.acid. <- identified_lipids_TP36_data_normalized_transposed_ana$FA.18.2..linoleic.acid.
subset_month3_hub_ID$FA.18.2..linoleic.acid. <- identified_postpartum_data_normalized_transposed_ana$FA.18.2..linoleic.acid.

# View the subsetted baseline dataset
print(subset_baseline_hub_ID)

#rename your colnames
colnames(subset_baseline_hub_ID)
colnames(subset_baseline_hub_ID) <- paste("baseline", colnames(subset_baseline_hub_ID), sep = "_")
colnames(subset_baseline_hub_ID)

colnames(subset_TP36_38weeks_hub_ID)
colnames(subset_TP36_38weeks_hub_ID) <- paste("TP36_38weeks", colnames(subset_TP36_38weeks_hub_ID), sep = "_")
colnames(subset_TP36_38weeks_hub_ID)

colnames(subset_month3_hub_ID)
colnames(subset_month3_hub_ID) <- paste("Postpartum", colnames(subset_month3_hub_ID), sep = "_")
colnames(subset_month3_hub_ID)


# Convert row names to a column in each dataframe
#make sure that rownames you see before conversion is growell_ID
rownames(subset_baseline_hub_ID)
#check if it wasnt done already
#subset_baseline_hub_ID$baseline_RowNames
rownames(ID_Updated_Combined_Traits_apo_analysis)
subset_baseline_hub_ID <- tibble(RowNames = rownames(subset_baseline_hub_ID), subset_baseline_hub_ID)
subset_TP36_38weeks_hub_ID <- tibble(RowNames = rownames(subset_TP36_38weeks_hub_ID), subset_TP36_38weeks_hub_ID)
subset_month3_hub_ID <- tibble(RowNames = rownames(subset_month3_hub_ID), subset_month3_hub_ID)
ID_Combined_Traits <- tibble(RowNames = rownames(ID_Updated_Combined_Traits_apo_analysis), ID_Updated_Combined_Traits_apo_analysis)

# Sequentially merge data frames
merged_data_ID_hublipids_traits <- subset_baseline_hub_ID %>%
  left_join(subset_TP36_38weeks_hub_ID, by = "RowNames") %>%
  left_join(subset_month3_hub_ID, by = "RowNames") %>%
  left_join(ID_Combined_Traits, by = "RowNames") 

head(merged_data_ID_hublipids_traits)

# Sequentially merge data frames. This is for a trial analysis. I want to see if it will work. 
#merged_data_ID_hublipids_traits_trial_analysis <- ID_Combined_Traits  %>%
 # left_join(subset_baseline_hub_ID, by = "RowNames") %>%
  #left_join(subset_TP36_38weeks_hub_ID, by = "RowNames") %>%
  #left_join(subset_month3_hub_ID, by = "RowNames") 


# sort your data by RowNames
# Sort the merged dataset by User_ID
merged_data_ID_hublipids_traits <- merged_data_ID_hublipids_traits[order(merged_data_ID_hublipids_traits$RowNames), ]
head(merged_data_ID_hublipids_traits)

# Write the merged data to an Excel file
write_xlsx(merged_data_ID_hublipids_traits, "Pareto_Norm_MergedData_ID_hublipids_traits.xlsx")


#=====================================================================================
#
#  Code chunk - Box plots and ANOVA for the hubLipids and Traits together
#
#=====================================================================================
#Time to do ANOVA #this for doing it for one sample long way
rownames(merged_data_ID_hublipids_traits) <- merged_data_ID_hublipids_traits$RowNames
merged_data_ID_hublipids_traits$baseline_TG.50.1
merged_data_ID_hublipids_traits$TP36_38weeks_TG.50.1
merged_data_ID_hublipids_traits$Postpartum_TG.50.1
merged_data_ID_hublipids_traits$apo
library(tidyr)

# Reshape just the relevant columns
long_data <- merged_data_ID_hublipids_traits[, c("RowNames", "apo", 
                                                 "baseline_TG.50.1", 
                                                 "TP36_38weeks_TG.50.1", 
                                                 "Postpartum_TG.50.1")]

long_data <- pivot_longer(
  data = long_data,
  cols = c(baseline_TG.50.1, TP36_38weeks_TG.50.1, Postpartum_TG.50.1),
  names_to = "Timepoint",
  values_to = "TG_50_1"
)

head(long_data)


long_data$Timepoint <- factor(long_data$Timepoint, 
                              levels = c("baseline_TG.50.1", "TP36_38weeks_TG.50.1", "Postpartum_TG.50.1"))
long_data$apo <- factor(long_data$apo, 
                         levels = c(0, 1), 
                         labels = c("NO_APO", "APO"))

head(long_data)
str(long_data)

# Ensure you have a subject identifier
long_data$RowNames <- factor(rep(1:(nrow(merged_data_ID_hublipids_traits)), each = 3))  # or use your real ID column

# Run repeated measures ANOVA
model <- aov(TG_50_1 ~ apo * Timepoint + Error(RowNames/Timepoint), data = long_data)

summary(model)

#since result was significant next is posthoc test
library(emmeans)

# Run pairwise comparisons for apo
emmeans_model <- emmeans(model, ~ apo)
pairs(emmeans_model)

#run pairwise difference at each Timepoint


# Run pairwise comparisons at each timepoint
emmeans_model_by_time <- emmeans(model, ~ apo | Timepoint)

# Now show pairwise comparisons
pairs(emmeans_model_by_time)


pairs(emmeans_model_by_time, adjust = "tukey")
#results are the same because you only have two groups per
#timepoint instead of 3 or more.



#Function to do ANOVA for one sample or multiple
#this funcntnio works but its using t-tests not the emeans
#so im using t-test. you can change the part about adjusting within the fucntion
#this is for one lipid at a time

run_repeated_measures_anova <- function(data,
                                        lipid_name,
                                        id_col = "RowNames",
                                        apo_col = "apo",
                                        time_levels = c("baseline", "TP36_38weeks", "Postpartum")) {
  library(tidyr)
  library(dplyr)
  
  # Prepare column names
  expected_cols <- c(id_col, apo_col, paste0(time_levels, "_", lipid_name))
  missing_cols <- setdiff(expected_cols, colnames(data))
  
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Reshape to long format
  df <- as.data.frame(data[, expected_cols])
  colnames(df) <- c("ID", "apo", time_levels)
  
  df <- df %>% filter(rowSums(is.na(select(., all_of(time_levels)))) < length(time_levels))
  
  long_df <- pivot_longer(df,
                          cols = all_of(time_levels),
                          names_to = "Timepoint",
                          values_to = "value")
  
  # Set factors
  long_df$Timepoint <- factor(long_df$Timepoint, levels = time_levels)
  long_df$apo <- factor(long_df$apo, levels = c(0, 1), labels = c("NO_APO", "APO"))
  long_df$ID <- factor(long_df$ID)
  
  # Run repeated measures ANOVA
  model <- aov(value ~ apo * Timepoint + Error(ID/Timepoint), data = long_df)
  
  cat("\n======================\n")
  cat("ANOVA for:", lipid_name, "\n")
  cat("======================\n")
  print(summary(model))
  
  # Pairwise comparisons per timepoint
  cat("\nPairwise comparisons (NO_APO vs APO) at each timepoint:\n")
  for (tp in time_levels) {
    tp_data <- long_df %>% filter(Timepoint == tp)
    res <- pairwise.t.test(tp_data$value, tp_data$apo, p.adjust.method = "holm")
    cat("\n🧪 Timepoint:", tp, "\n")
    print(res)
  }
}

run_repeated_measures_anova(
  data = merged_data_ID_hublipids_traits,
  lipid_name = "TG.50.1"
)


#Run repeated measures ANOVA across multiple lipids
#Extract main effect p-values from the ANOVA summary
#Extract per-timepoint unadjusted pairwise p-values
#Output everything in a single tidy results table
run_repeated_measures_anova_multi <- function(data,
                                              lipid_names,
                                              id_col = "RowNames",
                                              apo_col = "apo",
                                              time_levels = c("baseline", "TP36_38weeks", "Postpartum")) {
  library(tidyr)
  library(dplyr)
 #library(broom)
  
  # Ensure data is clean
  data <- as.data.frame(data)
  attr(data, "ts") <- NULL
  class(data) <- "data.frame"
  
  # Prepare results list
  all_results <- list()
  
  for (lipid_name in lipid_names) {
    expected_cols <- c(id_col, apo_col, paste0(time_levels, "_", lipid_name))
    missing_cols <- setdiff(expected_cols, colnames(data))
    if (length(missing_cols) > 0) {
      warning(paste("Skipping", lipid_name, "- missing columns:", paste(missing_cols, collapse = ", ")))
      next
    }
    
    df <- data[, expected_cols]
    colnames(df) <- c("ID", "apo", time_levels)
    df <- df %>% filter(rowSums(is.na(select(., all_of(time_levels)))) < length(time_levels))
    
    long_df <- pivot_longer(df,
                            cols = all_of(time_levels),
                            names_to = "Timepoint",
                            values_to = "value")
    
    long_df$Timepoint <- factor(long_df$Timepoint, levels = time_levels)
    long_df$apo <- factor(long_df$apo, levels = c(0, 1), labels = c("NO_APO", "APO"))
    long_df$ID <- factor(long_df$ID)
    
    # Run ANOVA
    model <- aov(value ~ apo * Timepoint + Error(ID/Timepoint), data = long_df)
    smry <- summary(model)
    
    # Extract p-values from ANOVA
    main_effects <- tibble(
      Lipid = lipid_name,
      apo_p = tryCatch(smry[[1]][[1]]["apo", "Pr(>F)"], error = function(e) NA),
      timepoint_p = tryCatch(smry[[2]][[1]]["Timepoint", "Pr(>F)"], error = function(e) NA),
      interaction_p = tryCatch(smry[[2]][[1]]["apo:Timepoint", "Pr(>F)"], error = function(e) NA)
    )
    
    # Per-timepoint pairwise comparisons
    pairwise_df <- lapply(time_levels, function(tp) {
      tp_data <- long_df %>% filter(Timepoint == tp)
      pw <- tryCatch(
        pairwise.t.test(tp_data$value, tp_data$apo, p.adjust.method = "none")$p.value[1, 1],
        error = function(e) NA
      )
      tibble(Lipid = lipid_name, Timepoint = tp, Pairwise_p = pw)
    }) %>% bind_rows()
    
    # Merge ANOVA and pairwise results
    merged <- pairwise_df %>%
      left_join(main_effects, by = "Lipid")
    
    all_results[[lipid_name]] <- merged
  }
  
  # Combine all results
  result_df <- bind_rows(all_results)
  
  # Return results
  return(result_df)
}

lipid_list <- c("TG.52.1", "TG.48.1", "TG.50.0", "TG.48.0","TG.46.0","TG.54.3",
                "TG.54.2", "TG.52.2", "TG.50.2", "TG.50.1", "TG.48.2","FAHFA.18.0.20.2",
                "FAHFA.38.4", "X4.methyl.2.oxopentanoic.acid", "FA.18.1..oleic.acid.",
                "FA.18.2..linoleic.acid.")
results <- run_repeated_measures_anova_multi(
  data = merged_data_ID_hublipids_traits,
  lipid_names = lipid_list
)

# View or save
print(results)
write.csv(results, "repeated_measures_anova_results.csv", row.names = FALSE)


#plot1
library(ggplot2)

ggplot(long_data, aes(x = Timepoint, y = TG_50_1, color = apo, group = RowNames)) +
  geom_line(alpha = 0.3) +
  stat_summary(aes(group = apo), fun = mean, geom = "line", size = 1.5) +
  labs(title = "TG.50.1 over Time by APO Status",
       y = "Lipid Intensity", x = "Timepoint") +
  theme_minimal()

#plot2
library(ggplot2)
library(ggpubr)  # for stat_pvalue_manual

# Rebuild the dataset if needed
# long_data should already have: RowNames, Timepoint, TG_50_1, apo

# Basic boxplot
p <- ggplot(long_data, aes(x = Timepoint, y = TG_50_1, fill = apo)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.7) +
  geom_jitter(aes(color = apo), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 1, alpha = 0.7) +
  labs(
    title = "Repeated Measures: TG.50.1 across Timepoints",
    y = "TG.50.1 Intensity",
    x = "Timepoint"
  ) +
  theme_minimal() +
  theme(legend.position = "top")
p

#plot 3 # this worked but it had some grid lines, no sqaures, no main effect, the colors were a bit blunt

library(ggplot2)
library(grid)
library(gridExtra)

plotAPOtimepoints <- function(data, pvalues_list, title_prefix = "TG.50.1: ", file = NULL, width = 12, height = 4) {
  plots <- list()
  tp_labels <- c("baseline_TG.50.1", "TP36_38weeks_TG.50.1", "Postpartum_TG.50.1")
  pretty_labels <- c("Baseline", "TP36_38", "Postpartum")
  panel_labels <- c("A", "B", "C")
  
  for (i in seq_along(tp_labels)) {
    tp <- tp_labels[i]
    df_tp <- subset(data, Timepoint == tp)
    
    if (nrow(df_tp) == 0) {
      warning(paste("Skipping", tp, "due to missing data"))
      next
    }
    
    p_val <- pvalues_list[[tp]]
    stars <- ifelse(p_val < 0.001, "***", ifelse(p_val < 0.01, "**", ifelse(p_val < 0.05, "*", "ns")))
    subtitle <- paste0(title_prefix, pretty_labels[i], " (p = ", signif(p_val, 3), ") ", stars)
    
    p <- ggplot(df_tp, aes(x = apo, y = TG_50_1, fill = apo)) +
      geom_boxplot(position = position_dodge(0.8), width = 0.7, outlier.shape = NA, alpha = 0.7) +
      geom_jitter(aes(color = apo), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
                  size = 1.5, alpha = 0.7) +
      scale_fill_manual(values = c("NO_APO" = "skyblue", "APO" = "salmon")) +
      scale_color_manual(values = c("NO_APO" = "skyblue4", "APO" = "firebrick")) +
      labs(
        title = panel_labels[i],
        subtitle = subtitle,
        x = "APO Status",
        y = "TG.50.1 Intensity"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11, face = "italic"),
        axis.text.x = element_text(angle = 20, vjust = 0.8, hjust = 0.8),
        legend.position = "none",
        strip.text = element_text(size = 13, face = "bold")
      )
    
    plots[[i]] <- p
  }
  
  if (!is.null(file)) {
    ggsave(file, plot = grid.arrange(grobs = plots, ncol = 3), width = width, height = height)
  } else {
    grid.arrange(grobs = plots, ncol = 3)
  }
}

pvalues_list <- list(
  "baseline_TG.50.1" = 0.001821792,
  "TP36_38weeks_TG.50.1" = 0.861976444,
  "Postpartum_TG.50.1" = 0.008845569
)

plotAPOtimepoints(long_data, pvalues_list)

#plot5 this worked. 
library(ggplot2)
library(gridExtra)
library(grid)

plotAPOtimepoints <- function(data, pvalues_list, main_effect_p = NULL, title_prefix = "TG.50.1: ", file = NULL, width = 12, height = 4) {
  plots <- list()
  tp_labels <- c("baseline_TG.50.1", "TP36_38weeks_TG.50.1", "Postpartum_TG.50.1")
  pretty_labels <- c("TP10_16", "TP36_38", "TP_Postpartum")
#  pretty_labels <- c("Baseline", "TP36_38", "Postpartum")
  panel_labels <- c("A", "B", "C")
  
  for (i in seq_along(tp_labels)) {
    tp <- tp_labels[i]
    df_tp <- subset(data, Timepoint == tp)
    
    if (nrow(df_tp) == 0) {
      warning(paste("Skipping", tp, "due to missing data"))
      next
    }
    
    p_val <- pvalues_list[[tp]]
    stars <- ifelse(p_val < 0.001, "***", ifelse(p_val < 0.01, "**", ifelse(p_val < 0.05, "*", "ns")))
    subtitle <- paste0(title_prefix, pretty_labels[i], " (p = ", signif(p_val, 3), ") ", stars)
    
    p <- ggplot(df_tp, aes(x = apo, y = TG_50_1, fill = apo)) +
      geom_boxplot(position = position_dodge(0.8), width = 0.7, outlier.shape = NA, alpha = 0.7, color = "black") +
      geom_jitter(aes(color = apo), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
                  size = 1.5, alpha = 0.7) +
      scale_fill_manual(values = c("NO_APO" = "turquoise", "APO" = "red")) +  # nice colors
      scale_color_manual(values = c("NO_APO" = "turquoise", "APO" = "red")) +  # darker matching line colors
      labs(
        title = panel_labels[i],
        subtitle = subtitle,
        x = "APO Status",
        y = "TG.50.1 Intensity"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        panel.background = element_rect(fill = "white", color = "black", size = 1),  # white background and black frame
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11, face = "italic"),
        axis.text.x = element_text(angle = 20, vjust = 0.8, hjust = 0.8),
        legend.position = "none",
        strip.text = element_text(size = 13, face = "bold")
      )
    
    plots[[i]] <- p
  }
  
  # Arrange plots
  combined_plot <- arrangeGrob(grobs = plots, ncol = 3)
  
  # If main model p-value is provided, add it on top
  if (!is.null(main_effect_p)) {
    main_effect_label <- textGrob(
      paste0("Main effect of APO: p = ", signif(main_effect_p, 3)),
      gp = gpar(fontsize = 13, fontface = "italic")
    )
    
    final_plot <- grid.arrange(main_effect_label, combined_plot, nrow = 2, heights = c(0.1, 1))
  } else {
    final_plot <- combined_plot
  }
  
  # Save if file path is given
  if (!is.null(file)) {
    ggsave(file, plot = final_plot, width = width, height = height)
  } else {
    grid.draw(final_plot)
  }
}

pvalues_list <- list(
  "baseline_TG.50.1" = 0.001821792,
  "TP36_38weeks_TG.50.1" = 0.861976444,
  "Postpartum_TG.50.1" = 0.008845569
)

main_effect_p <- 0.00279  # from your ANOVA model

plotAPOtimepoints(long_data, pvalues_list, main_effect_p = main_effect_p)

library(ggplot2)
library(gridExtra)
library(grid)

#plot 5 this also worked but different colors and saving options enhanced. 
#current go to option
plotAPOtimepoints <- function(data, pvalues_list, main_effect_p = NULL, title_prefix = "TG.50.1: ", file = NULL, width = 12, height = 4) {
  plots <- list()
  tp_labels <- c("baseline_TG.50.1", "TP36_38weeks_TG.50.1", "Postpartum_TG.50.1")
  pretty_labels <- c("TP10_16", "TP36_38", "TP_Postpartum")
  #pretty_labels <- c("Baseline", "TP36_38", "Postpartum")
  panel_labels <- c("A", "B", "C")
  
  for (i in seq_along(tp_labels)) {
    tp <- tp_labels[i]
    df_tp <- subset(data, Timepoint == tp)
    
    if (nrow(df_tp) == 0) {
      warning(paste("Skipping", tp, "due to missing data"))
      next
    }
    
    p_val <- pvalues_list[[tp]]
    stars <- ifelse(p_val < 0.001, "***", ifelse(p_val < 0.01, "**", ifelse(p_val < 0.05, "*", "ns")))
    subtitle <- paste0(title_prefix, pretty_labels[i], " (p = ", signif(p_val, 3), ") ", stars)
    
    p <- ggplot(df_tp, aes(x = apo, y = TG_50_1, fill = apo)) +
      geom_boxplot(position = position_dodge(0.8), width = 0.7, outlier.shape = NA, alpha = 0.7, color = "black") +
      geom_jitter(aes(color = apo), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
                  size = 1.5, alpha = 0.7) +
      scale_fill_manual(values = c("NO_APO" = "#1E90FF", "APO" = "red")) +
      scale_color_manual(values = c("NO_APO" = "#1E90FF", "APO" = "red")) +
      labs(
        title = panel_labels[i],
        subtitle = subtitle,
        x = "APO Status",
        y = "TG.50.1 Intensity"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        panel.background = element_rect(fill = "white", color = "black", size = 1),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11, face = "italic"),
        axis.text.x = element_text(angle = 20, vjust = 0.8, hjust = 0.8),
        legend.position = "none",
        strip.text = element_text(size = 13, face = "bold")
      )
    
    plots[[i]] <- p
  }
  
  # Arrange plots
  combined_plot <- arrangeGrob(grobs = plots, ncol = 3)
  
  # Add main effect p-value title
  if (!is.null(main_effect_p)) {
    main_effect_label <- textGrob(
      paste0("Main effect of APO: p = ", signif(main_effect_p, 3)),
      gp = gpar(fontsize = 13, fontface = "italic")
    )
    
    final_plot <- grid.arrange(main_effect_label, combined_plot, nrow = 2, heights = c(0.1, 1))
  } else {
    final_plot <- combined_plot
  }
  
  # Save if file path is given
  if (!is.null(file)) {
    message("Saving PDF to: ", file)
    ggsave(filename = file, plot = final_plot, device = "pdf", width = width, height = height)
  } else {
    grid.draw(final_plot)
  }
}

#no saving
plotAPOtimepoints(long_data, pvalues_list, main_effect_p = main_effect_p)

#save pdf
plotAPOtimepoints(long_data, pvalues_list, main_effect_p = main_effect_p, file = "TG50_1_APO_Timepoints.pdf")

#adjust height and width
plotAPOtimepoints(long_data, pvalues_list, main_effect_p = main_effect_p, file = "TG50_1_APO_Timepoints_big.pdf", width = 14, height = 5)


#=====================================================================================
#
#  Code chunk - Box plots and ANOVA for the hubLipids and Traits together
#.  Plot 6 - Current Code for making figures for APOs but make sure all stats was calculated
#   the lipid of interest
#=====================================================================================


#plot 6, I want to try create a function that does this for a number of lipids together
#This worked. This will be go to to make multiple images
#adjusted functions to make sure that y-axis is the same

library(ggplot2)
library(gridExtra)
library(grid)

plotAPOtimepoints <- function(data, pvalues_list, main_effect_p = NULL, lipid_name = "TG.50.1",
                              title_prefix = "TG.50.1: ", file = NULL, width = 12, height = 4) {
  plots <- list()
  tp_labels <- paste0(c("baseline", "TP36_38weeks", "Postpartum"), "_", lipid_name)
  pretty_labels <- c("TP10_16", "TP36_38", "TP_Postpartum")
  #panel_labels <- c("A", "B", "C")
  
  # Determine consistent y-axis limits
  global_ymin <- min(data[[lipid_name]], na.rm = TRUE)
  global_ymax <- max(data[[lipid_name]], na.rm = TRUE)
  
  for (i in seq_along(tp_labels)) {
    tp <- tp_labels[i]
    df_tp <- subset(data, Timepoint == tp)
    
    if (nrow(df_tp) == 0) {
      warning(paste("Skipping", tp, "due to missing data"))
      next
    }
    
    p_val <- pvalues_list[[tp]]
    stars <- ifelse(p_val < 0.001, "***", ifelse(p_val < 0.01, "**", ifelse(p_val < 0.05, "*", "ns")))
    subtitle <- paste0(title_prefix, pretty_labels[i], " (p = ", signif(p_val, 3), ") ", stars)
    
    p <- ggplot(df_tp, aes_string(x = "apo", y = lipid_name, fill = "apo")) +
      geom_boxplot(position = position_dodge(0.8), width = 0.7, outlier.shape = NA, alpha = 0.7, color = "black") +
      geom_jitter(aes_string(color = "apo"), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
                  size = 1.5, alpha = 0.7) +
      scale_y_continuous(limits = c(global_ymin, global_ymax)) +   # <--- ADDED LINE
      scale_fill_manual(values = c("NO_APO" = "#1E90FF", "APO" = "red")) +
      scale_color_manual(values = c("NO_APO" = "#1E90FF", "APO" = "red")) +
      labs(
        #   title = panel_labels[i],
            subtitle = subtitle,
        x = "APO Status",
        y = paste(lipid_name, "Intensity")
      ) +
      theme_minimal(base_size = 12) +
      theme(
        panel.background = element_rect(fill = "white", color = "black", size = 1),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11, face = "italic"),
        axis.text.x = element_text(angle = 20, vjust = 0.8, hjust = 0.8),
        legend.position = "none",
        strip.text = element_text(size = 13, face = "bold")
      )
    
    plots[[i]] <- p
  }
  
  combined_plot <- arrangeGrob(grobs = plots, ncol = 3)
  
  if (!is.null(main_effect_p)) {
    main_effect_label <- textGrob(
      paste0("Main effect of APO: p = ", signif(main_effect_p, 3)),
      gp = gpar(fontsize = 13, fontface = "italic")
    )
    final_plot <- grid.arrange(main_effect_label, combined_plot, nrow = 2, heights = c(0.1, 1))
  } else {
    final_plot <- combined_plot
  }
  
  if (!is.null(file)) {
    message("Saving PDF to: ", file)
    ggsave(filename = file, plot = final_plot, device = "pdf", width = width, height = height)
  } else {
    grid.draw(final_plot)
  }
}


plot_lipid_boxplot_with_pvalues <- function(lipid_name,
                                            data,
                                            results_tibble,
                                            output_file = NULL,
                                            width = 12,
                                            height = 5) {
  library(dplyr)
  library(tidyr)
  library(tibble)
  
  timepoints <- c("baseline", "TP36_38weeks", "Postpartum")
  full_names <- paste0(timepoints, "_", lipid_name)
  
  plot_data <- data %>%
    select(RowNames, apo, all_of(full_names)) %>%
    pivot_longer(cols = starts_with(timepoints),
                 names_to = "Timepoint", values_to = "value") %>%
    mutate(
      apo = factor(apo, levels = c(0, 1), labels = c("NO_APO", "APO")),
      Timepoint = factor(Timepoint, levels = full_names)
    )
  colnames(plot_data)[colnames(plot_data) == "value"] <- lipid_name
  
  # Extract p-values
  lipid_res <- results_tibble %>% filter(Lipid == lipid_name)
  
  pvalues_list <- lipid_res %>%
    mutate(Timepoint = paste0(Timepoint, "_", lipid_name)) %>%
    { setNames(.$Pairwise_p, .$Timepoint) }
  
  main_effect_p <- unique(lipid_res$apo_p)[1]
  
  # Plot
  plotAPOtimepoints(
    data = plot_data,
    pvalues_list = pvalues_list,
    main_effect_p = main_effect_p,
    lipid_name = lipid_name,
    title_prefix = paste0(lipid_name, ": "),
    file = output_file,
    width = width,
    height = height
  )
}

lipids <- c("TG.50.1", "TG.48.1", "TG.54.3")

lipid_list <- c("TG.52.1", "TG.48.1", "TG.50.0", "TG.48.0","TG.46.0","TG.54.3",
                "TG.54.2", "TG.52.2", "TG.50.2", "TG.50.1", "TG.48.2","FAHFA.18.0.20.2",
                "FAHFA.38.4", "X4.methyl.2.oxopentanoic.acid", "FA.18.1..oleic.acid.",
                "FA.18.2..linoleic.acid.")

for (lip in lipid_list) {
  plot_lipid_boxplot_with_pvalues(
    lipid_name = lip,
    data = merged_data_ID_hublipids_traits,
    results_tibble = results,
    output_file = paste0("Panel_", lip, "_APO_Boxplot.pdf")
  )
}



#one lipids
plot_lipid_boxplot_with_pvalues(
  lipid_name = "X4.methyl.2.oxopentanoic.acid",
  data = merged_data_ID_hublipids_traits,
  results_tibble = results,
  height = 7,# slightly taller for room,
  width = 15,
  output_file = "X4.methyl.2.oxopentanoic.acid_Boxplot_Panel.pdf"
)

#=====================================================================================
#
#  Code chunk - Box plots and ANOVA for the hubLipids and PPWR traits together
#
#=====================================================================================
library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(lme4)
library(lmerTest)
library(emmeans)

#convert the variable to a factor
merged_data_ID_hublipids_traits$ppwr_e <- factor(merged_data_ID_hublipids_traits$ppwr_e,
                                                 levels = c(0, 1), labels = c("NO_PPWR", "PPWR"))
#this is the function for the repeated measures with ppwr_e variables
run_repeated_measures_anova_multi <- function(data,
                                              lipid_names,
                                              id_col = "RowNames",
                                              group_col = "ppwr_e",
                                              time_levels = c("baseline", "TP36_38weeks", "Postpartum")) {
  all_results <- list()
  
  for (lipid_name in lipid_names) {
    expected_cols <- c(id_col, group_col, paste0(time_levels, "_", lipid_name))
    if (!all(expected_cols %in% colnames(data))) next
    
    #df <- data[, expected_cols]
    df <- data[, expected_cols] %>%
      filter(!is.na(.data[[group_col]]))  # exclude NAs *only* here
    colnames(df) <- c("ID", "group", time_levels)
    df <- df %>% filter(rowSums(is.na(select(., all_of(time_levels)))) < length(time_levels))
    
    long_df <- pivot_longer(df, cols = all_of(time_levels),
                            names_to = "Timepoint", values_to = "value") %>%
      mutate(Timepoint = factor(Timepoint, levels = time_levels),
             group = factor(group),
             ID = factor(ID))
    
    model <- aov(value ~ group * Timepoint + Error(ID/Timepoint), data = long_df)
    smry <- summary(model)
    
    main_effects <- tibble(
      Lipid = lipid_name,
      group_p = tryCatch(smry[[1]][[1]]["group", "Pr(>F)"], error = function(e) NA),
      timepoint_p = tryCatch(smry[[2]][[1]]["Timepoint", "Pr(>F)"], error = function(e) NA),
      interaction_p = tryCatch(smry[[2]][[1]]["group:Timepoint", "Pr(>F)"], error = function(e) NA)
    )
    
    pairwise_df <- lapply(time_levels, function(tp) {
      tp_data <- long_df %>% filter(Timepoint == tp)
      pw <- tryCatch(
        pairwise.t.test(tp_data$value, tp_data$group, p.adjust.method = "none")$p.value[1, 1],
        error = function(e) NA
      )
      tibble(Lipid = lipid_name, Timepoint = tp, Pairwise_p = pw)
    }) %>% bind_rows()
    
    merged <- left_join(pairwise_df, main_effects, by = "Lipid")
    all_results[[lipid_name]] <- merged
  }
  
  bind_rows(all_results)
}

#linear mixed model function
run_repeated_measures_lmer <- function(data,
                                       lipid_name,
                                       id_col = "RowNames",
                                       cont_var = "ppwr",
                                       time_levels = c("baseline", "TP36_38weeks", "Postpartum")) {
  #df <- data[, c(id_col, cont_var, paste0(time_levels, "_", lipid_name))]
  df <- data[, c(id_col, cont_var, paste0(time_levels, "_", lipid_name))] %>%
    filter(!is.na(.data[[cont_var]]))  # Exclude only here
  
  colnames(df) <- c("ID", "ppwr", time_levels)
  
  df <- df %>% filter(rowSums(is.na(select(., all_of(time_levels)))) < length(time_levels))
  
  long_df <- pivot_longer(df, cols = all_of(time_levels),
                          names_to = "Timepoint", values_to = "value") %>%
    mutate(Timepoint = factor(Timepoint, levels = time_levels),
           ID = factor(ID))
  
  model <- lmer(value ~ ppwr * Timepoint + (1|ID), data = long_df)
  
  cat("\n======================\n")
  cat("LME for:", lipid_name, "\n")
  cat("======================\n")
  print(anova(model))
  
  return(model)
}


#plotAPOtimepoints <- function(data, pvalues_list, main_effect_p = NULL, lipid_name = "TG.50.1",
                            title_prefix = "TG.50.1: ", file = NULL, width = 12, height = 4) {
  plots <- list()
  tp_labels <- paste0(c("baseline", "TP36_38weeks", "Postpartum"), "_", lipid_name)
  pretty_labels <- c("TP10_16", "TP36_38", "TP_Postpartum")
  #panel_labels <- c("A", "B", "C")
  
  # Determine consistent y-axis limits
  global_ymin <- min(data[[lipid_name]], na.rm = TRUE)
  global_ymax <- max(data[[lipid_name]], na.rm = TRUE)
  
  for (i in seq_along(tp_labels)) {
    tp <- tp_labels[i]
    df_tp <- subset(data, Timepoint == tp)
    
    if (nrow(df_tp) == 0) {
      warning(paste("Skipping", tp, "due to missing data"))
      next
    }
    
    p_val <- pvalues_list[[tp]]
    stars <- ifelse(p_val < 0.001, "***", ifelse(p_val < 0.01, "**", ifelse(p_val < 0.05, "*", "ns")))
    subtitle <- paste0(title_prefix, pretty_labels[i], " (p = ", signif(p_val, 3), ") ", stars)
    
    p <- ggplot(df_tp, aes_string(x = "group", y = lipid_name, fill = "group")) +
      geom_boxplot(position = position_dodge(0.8), width = 0.7, outlier.shape = NA, alpha = 0.7, color = "black") +
      geom_jitter(aes_string(color = "group"), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
                  size = 1.5, alpha = 0.7)  +
      scale_y_continuous(limits = c(global_ymin, global_ymax)) +   # <--- ADDED LINE
      scale_fill_manual(values = c("NO_APO" = "#1E90FF", "APO" = "red")) +
      scale_color_manual(values = c("NO_APO" = "#1E90FF", "APO" = "red")) +
      labs(
        #   title = panel_labels[i],
        subtitle = subtitle,
        x = "APO Status",
        y = paste(lipid_name, "Intensity")
      ) +
      theme_minimal(base_size = 12) +
      theme(
        panel.background = element_rect(fill = "white", color = "black", size = 1),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11, face = "italic"),
        axis.text.x = element_text(angle = 20, vjust = 0.8, hjust = 0.8),
        legend.position = "none",
        strip.text = element_text(size = 13, face = "bold")
      )
    
    plots[[i]] <- p
  }
  
  combined_plot <- arrangeGrob(grobs = plots, ncol = 3)
  
  if (!is.null(main_effect_p)) {
    main_effect_label <- textGrob(
      paste0("Main effect of APO: p = ", signif(main_effect_p, 3)),
      gp = gpar(fontsize = 13, fontface = "italic")
    )
    final_plot <- grid.arrange(main_effect_label, combined_plot, nrow = 2, heights = c(0.1, 1))
  } else {
    final_plot <- combined_plot
  }
  
  if (!is.null(file)) {
    message("Saving PDF to: ", file)
    ggsave(filename = file, plot = final_plot, device = "pdf", width = width, height = height)
  } else {
    grid.draw(final_plot)
  }
}

plotAPOtimepoints <- function(data, pvalues_list, main_effect_p = NULL, lipid_name = "TG.50.1",
                              title_prefix = "TG.50.1: ", file = NULL, width = 12, height = 4) {
  plots <- list()
  tp_labels <- paste0(c("baseline", "TP36_38weeks", "Postpartum"), "_", lipid_name)
  pretty_labels <- c("TP10_16", "TP36_38", "TP_Postpartum")
  
  global_ymin <- min(data[[lipid_name]], na.rm = TRUE)
  global_ymax <- max(data[[lipid_name]], na.rm = TRUE)
  
  for (i in seq_along(tp_labels)) {
    tp <- tp_labels[i]
    df_tp <- subset(data, Timepoint == tp)
    
    if (nrow(df_tp) == 0) {
      warning(paste("Skipping", tp, "due to missing data"))
      next
    }
    
    p_val <- pvalues_list[[tp]]
    stars <- ifelse(p_val < 0.001, "***", ifelse(p_val < 0.01, "**", ifelse(p_val < 0.05, "*", "ns")))
    subtitle <- paste0(title_prefix, pretty_labels[i], " (p = ", signif(p_val, 3), ") ", stars)
    
    group_colors <- c("NO_APO" = "#1E90FF", "APO" = "red",
                      "NO_PPWR" = "#1E90FF", "PPWR" = "red")
    
    # Limit to only the groups present in the current subset
    present_groups <- unique(df_tp$group)
    group_colors <- group_colors[names(group_colors) %in% present_groups]
    
    
    p <- ggplot(df_tp, aes(x = group, y = .data[[lipid_name]], fill = group)) +
      geom_boxplot(position = position_dodge(0.8), width = 0.7, outlier.shape = NA, alpha = 0.7, color = "black") +
      geom_jitter(aes(color = group), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
                  size = 1.5, alpha = 0.7) +
      scale_y_continuous(limits = c(global_ymin, global_ymax)) +
      scale_fill_manual(values = group_colors) +
      scale_color_manual(values = group_colors) +
      labs(
        subtitle = subtitle,
        x = "Group",
        y = paste(lipid_name, "Intensity")
      ) +
      theme_minimal(base_size = 12) +
      theme(
        panel.background = element_rect(fill = "white", color = "black", size = 1),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11, face = "italic"),
        axis.text.x = element_text(angle = 20, vjust = 0.8, hjust = 0.8),
        legend.position = "none",
        strip.text = element_text(size = 13, face = "bold")
      )
    
    plots[[i]] <- p
  }
  
  combined_plot <- arrangeGrob(grobs = plots, ncol = 3)
  
  if (!is.null(main_effect_p)) {
    main_effect_label <- textGrob(
      paste0("Main effect of group: p = ", signif(main_effect_p, 3)),
      gp = gpar(fontsize = 13, fontface = "italic")
    )
    final_plot <- grid.arrange(main_effect_label, combined_plot, nrow = 2, heights = c(0.1, 1))
  } else {
    final_plot <- combined_plot
  }
  
  if (!is.null(file)) {
    message("Saving PDF to: ", file)
    ggsave(filename = file, plot = final_plot, device = "pdf", width = width, height = height)
  } else {
    grid.draw(final_plot)
  }
}

#function to plot the categorical variables. Make sure you run the helperfunction plotAPOtimepoints()
plot_lipid_boxplot_with_pvalues <- function(lipid_name,
                                            data,
                                            results_tibble,
                                            group_col = "ppwr_e",
                                            time_levels = c("baseline", "TP36_38weeks", "Postpartum"),
                                            output_file = NULL,
                                            width = 12,
                                            height = 5) {
  full_names <- paste0(time_levels, "_", lipid_name)
  
  plot_data <- data %>%
    filter(!is.na(.data[[group_col]])) %>%
    select(RowNames, all_of(group_col), all_of(full_names)) %>%
    pivot_longer(cols = all_of(full_names),
                 names_to = "Timepoint", values_to = "value") %>%
    mutate(
      group = factor(.data[[group_col]]),
      Timepoint = factor(Timepoint, levels = paste0(time_levels, "_", lipid_name))
    )
  
  colnames(plot_data)[colnames(plot_data) == "value"] <- lipid_name
  
  pvalues_list <- results_tibble %>%
    filter(Lipid == lipid_name) %>%
    mutate(Timepoint = paste0(Timepoint, "_", lipid_name)) %>%
    { setNames(.$Pairwise_p, .$Timepoint) }
  
  main_effect_p <- unique(results_tibble$group_p[results_tibble$Lipid == lipid_name])[1]
  
  plotAPOtimepoints(
    data = plot_data,
    pvalues_list = pvalues_list,
    main_effect_p = main_effect_p,
    lipid_name = lipid_name,
    title_prefix = paste0(lipid_name, ": "),
    file = output_file,
    width = width,
    height = height
  )
}

#plot continous straight#
plot_continuous_trait_lipid <- function(lipid_name, data,
                                        trait = "ppwr", time_levels = c("baseline", "TP36_38weeks", "Postpartum"),
                                        width = 14, height = 5, output_file = NULL) {
  long_df <- data %>%
    filter(!is.na(.data[[trait]])) %>%
    select(RowNames, all_of(trait), paste0(time_levels, "_", lipid_name)) %>%
    pivot_longer(cols = paste0(time_levels, "_", lipid_name),
                 names_to = "Timepoint", values_to = "value") %>%
    mutate(Timepoint = factor(Timepoint, levels = paste0(time_levels, "_", lipid_name)))
  
  p <- ggplot(long_df, aes_string(x = trait, y = "value")) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    facet_wrap(~Timepoint, scales = "free_y") +
    labs(
      title = paste("Association of", lipid_name, "with", trait),
      x = trait,
      y = paste(lipid_name, "Intensity")
    ) +
    theme_minimal()
  
  if (!is.null(output_file)) {
    ggsave(output_file, plot = p, width = width, height = height)
  } else {
    print(p)
  }
}

#example of usage for categorical
lipid_list <- c("TG.42.1","TG.60.2", "TG.64.2",
                "TG.46.0" ,"TG.46.1","TG.58.3","TG.58.8" ,"TG.58.9")

results_ppwr_e <- run_repeated_measures_anova_multi(
  data = merged_data_ID_hublipids_traits,
  lipid_names = lipid_list,
  group_col = "ppwr_e"
)

for (lip in lipid_list) {
  plot_lipid_boxplot_with_pvalues(
    lipid_name = lip,
    data = merged_data_ID_hublipids_traits,
    results_tibble = results_ppwr_e,
    group_col = "ppwr_e",
    output_file = paste0("Panel_", lip, "_PPWR_Boxplot.pdf")
  )
}
#export and save the results

# Save to CSV
write.csv(results_ppwr_e, "Repeated_Measures_ANOVA_ppwr_e_results.csv", row.names = FALSE)

#for continous variables
for (lip in lipid_list) {
  run_repeated_measures_lmer(
    data = merged_data_ID_hublipids_traits,
    lipid_name = lip
  )
  
  plot_continuous_trait_lipid(
    lipid_name = lip,
    data = merged_data_ID_hublipids_traits,
    trait = "ppwr",
    output_file = paste0("Panel_", lip, "_PPWR_Continuous_Regression.pdf")
  )
}

#save results
#this didnt work. Problem for another day
extract_lmer_summary <- function(model, lipid_name) {
  smry <- anova(model)
  smry_df <- as.data.frame(smry)
  smry_df$Effect <- rownames(smry_df)
  smry_df$Lipid <- lipid_name
  rownames(smry_df) <- NULL
  smry_df <- smry_df[, c("Lipid", "Effect", "Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")]
  return(smry_df)
}

all_lmer_stats <- list()

for (lip in lipid_list) {
  model <- run_repeated_measures_lmer(
    data = merged_data_ID_hublipids_traits,
    lipid_name = lip
  )
  
  lmer_stats <- extract_lmer_summary(model, lip)
  all_lmer_stats[[lip]] <- lmer_stats
}

lmer_results_df <- bind_rows(all_lmer_stats)

# Save to CSV
write.csv(lmer_results_df, "Repeated_Measures_LMER_ppwr_results.csv", row.names = FALSE)


#=====================================================================================
#
#  Code chunk 11- MEs and traits plotted 
#
#=====================================================================================

#===========================MEpurple vs adverse outcome==========================================================

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# Convert group to a factor
merged_data_ID_adverse$ExperimentalGroup.x <- factor(merged_data_ID_adverse$ExperimentalGroup.x, levels = c(0, 1), labels = c("Control", "Intervention"))

# Convert new_adverse_outcome to a factor with levels "No" and "Yes"
merged_data_ID_adverse$new_adverse_outcome <- factor(merged_data_ID_adverse$new_adverse_outcome, levels = c(0, 1), labels = c("No", "Yes"))

# Convert data to long format and rename to data_long_MEpurple
data_long_MEpurple <- merged_data_ID_adverse %>%
  select(RowNames, ExperimentalGroup.x, new_adverse_outcome, 
         Identified_Baseline.data.MEpurple, Identified_TP36_38weeks.data.MEpurple, Identified_Postpartum.data.MEpurple) %>%
  gather(key = "time_point", value = "MEpurple", 
         Identified_Baseline.data.MEpurple, Identified_TP36_38weeks.data.MEpurple, Identified_Postpartum.data.MEpurple)

# Rename time points for better readability and order them
data_long_MEpurple$time_point <- factor(data_long_MEpurple$time_point, 
                                        levels = c("Identified_Baseline.data.MEpurple", "Identified_TP36_38weeks.data.MEpurple", "Identified_Postpartum.data.MEpurple"),
                                        labels = c("Baseline", "36-38 Weeks", "Postpartum"))

# Function to perform Wilcoxon test for adverse outcome and get p-values
get_p_value_adverse <- function(data, time_point) {
  data_subset <- data %>% filter(time_point == !!time_point)
  if (length(unique(data_subset$new_adverse_outcome)) > 1) {
    test <- wilcox.test(MEpurple ~ new_adverse_outcome, data = data_subset, paired = FALSE)
    return(test$p.value)
  } else {
    return(NA)
  }
}

# Function to perform Wilcoxon test for intervention/control group and get p-values for each level of adverse outcome
get_p_value_group <- function(data, time_point) {
  data_subset <- data %>% filter(time_point == !!time_point)
  p_values <- data_subset %>%
    group_by(new_adverse_outcome) %>%
    summarise(p_value = ifelse(length(unique(ExperimentalGroup.x)) > 1, wilcox.test(MEpurple ~ ExperimentalGroup.x)$p.value, NA))
  return(p_values)
}

# Calculate p-values for each time point and each level of adverse outcome
p_values_group <- do.call(rbind, lapply(unique(data_long_MEpurple$time_point), function(tp) {
  get_p_value_group(data_long_MEpurple, tp) %>%
    mutate(time_point = tp)
}))

# Calculate p-values for adverse outcome comparisons
p_values_adverse <- sapply(unique(data_long_MEpurple$time_point), function(tp) get_p_value_adverse(data_long_MEpurple, tp))

# Add significance annotations based on p-values
signif_labels_adverse <- ifelse(p_values_adverse < 0.05, "*", "")
signif_labels_group <- ifelse(p_values_group$p_value < 0.05, "*", "")

# Create a data frame for annotation positions and labels for adverse outcome comparisons
annotation_df_adverse <- data.frame(
  time_point = unique(data_long_MEpurple$time_point),
  x = 1.5,
  y = 0.9,
  label = signif_labels_adverse
)

# Create a data frame for annotation positions and labels for group comparisons within adverse outcome levels
annotation_df_group <- data.frame(
  time_point = p_values_group$time_point,
  new_adverse_outcome = p_values_group$new_adverse_outcome,
  x = 1.5,
  y = 0.9,
  label = signif_labels_group
)

# Create the first set of box plots: MEpurple by Adverse Outcome
plot1 <- ggplot(data_long_MEpurple, aes(x = time_point, y = MEpurple, fill = new_adverse_outcome)) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(title = "MEpurple Values by Adverse Outcome and Time Point",
       x = "Time Point",
       y = "MEpurple Value",
       fill = "Adverse Outcome") +
  theme_minimal() +
  ylim(-0.6, 1) +
  geom_text(data = annotation_df_adverse, aes(x = time_point, y = y, label = label), size = 5, inherit.aes = FALSE, position = position_dodge(0.8))

# Create the second set of box plots: MEpurple by Adverse Outcome and Experimental Group
plot2 <- ggplot(data_long_MEpurple, aes(x = interaction(new_adverse_outcome, time_point), y = MEpurple, fill = ExperimentalGroup.x)) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(title = "MEpurple Values by Experimental Group within Adverse Outcome and Time Point",
       x = "Adverse Outcome and Time Point",
       y = "MEpurple Value",
       fill = "Experimental Group") +
  theme_minimal() +
  ylim(-0.6, 1) +
  scale_fill_manual(values = c("Control" = "blue", "Intervention" = "red")) +
  geom_text(data = annotation_df_group, aes(x = interaction(new_adverse_outcome, time_point), y = y, label = label), size = 5, inherit.aes = FALSE, position = position_dodge(0.8))

# Display both plots
grid.arrange(plot1, plot2, ncol = 1)

#===========================MEred vs adverse outcome==========================================================

# Analysis for MEred

# Convert data to long format and rename to data_long_MEred
data_long_MEred <- merged_data_ID_adverse %>%
  select(RowNames, ExperimentalGroup.x, new_adverse_outcome, 
         Identified_Baseline.data.MEred, Identified_TP36_38weeks.data.MEred, Identified_Postpartum.data.MEred) %>%
  gather(key = "time_point", value = "MEred", 
         Identified_Baseline.data.MEred, Identified_TP36_38weeks.data.MEred, Identified_Postpartum.data.MEred)

# Rename time points for better readability and order them
data_long_MEred$time_point <- factor(data_long_MEred$time_point, 
                                     levels = c("Identified_Baseline.data.MEred", "Identified_TP36_38weeks.data.MEred", "Identified_Postpartum.data.MEred"),
                                     labels = c("Baseline", "36-38 Weeks", "Postpartum"))

# Function to perform Wilcoxon test for adverse outcome and get p-values
get_p_value_adverse_MEred <- function(data, time_point) {
  data_subset <- data %>% filter(time_point == !!time_point)
  if (length(unique(data_subset$new_adverse_outcome)) > 1) {
    test <- wilcox.test(MEred ~ new_adverse_outcome, data = data_subset, paired = FALSE)
    return(test$p.value)
  } else {
    return(NA)
  }
}

# Function to perform Wilcoxon test for intervention/control group and get p-values for each level of adverse outcome
get_p_value_group_MEred <- function(data, time_point) {
  data_subset <- data %>% filter(time_point == !!time_point)
  p_values <- data_subset %>%
    group_by(new_adverse_outcome) %>%
    summarise(p_value = ifelse(length(unique(ExperimentalGroup.x)) > 1, wilcox.test(MEred ~ ExperimentalGroup.x)$p.value, NA))
  return(p_values)
}

# Calculate p-values for each time point and each level of adverse outcome
p_values_group_MEred <- do.call(rbind, lapply(unique(data_long_MEred$time_point), function(tp) {
  get_p_value_group_MEred(data_long_MEred, tp) %>%
    mutate(time_point = tp)
}))

# Calculate p-values for adverse outcome comparisons
p_values_adverse_MEred <- sapply(unique(data_long_MEred$time_point), function(tp) get_p_value_adverse_MEred(data_long_MEred, tp))

# Add significance annotations based on p-values
signif_labels_adverse_MEred <- ifelse(p_values_adverse_MEred < 0.05, "*", "")
signif_labels_group_MEred <- ifelse(p_values_group_MEred$p_value < 0.05, "*", "")

# Create a data frame for annotation positions and labels for adverse outcome comparisons
annotation_df_adverse_MEred <- data.frame(
  time_point = unique(data_long_MEred$time_point),
  x = 1.5,
  y = 0.9,
  label = signif_labels_adverse_MEred
)

# Create a data frame for annotation positions and labels for group comparisons within adverse outcome levels
annotation_df_group_MEred <- data.frame(
  time_point = p_values_group_MEred$time_point,
  new_adverse_outcome = p_values_group_MEred$new_adverse_outcome,
  x = 1.5,
  y = 0.9,
  label = signif_labels_group_MEred
)

# Create the first set of box plots: MEred by Adverse Outcome
plot_MEred_1 <- ggplot(data_long_MEred, aes(x = time_point, y = MEred, fill = new_adverse_outcome)) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(title = "MEred Values by Adverse Outcome and Time Point",
       x = "Time Point",
       y = "MEred Value",
       fill = "Adverse Outcome") +
  theme_minimal() +
  ylim(-0.5, 1) +
  geom_text(data = annotation_df_adverse_MEred, aes(x = time_point, y = y, label = label), size = 5, inherit.aes = FALSE, position = position_dodge(0.8))

# Create the second set of box plots: MEred by Adverse Outcome and Experimental Group
plot_MEred_2 <- ggplot(data_long_MEred, aes(x = interaction(new_adverse_outcome, time_point), y = MEred, fill = ExperimentalGroup.x)) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(title = "MEred Values by Experimental Group within Adverse Outcome and Time Point",
       x = "Adverse Outcome and Time Point",
       y = "MEred Value",
       fill = "Experimental Group") +
  theme_minimal() +
  ylim(-0.5, 1) +
  scale_fill_manual(values = c("Control" = "blue", "Intervention" = "red")) +
  geom_text(data = annotation_df_group_MEred, aes(x = interaction(new_adverse_outcome, time_point), y = y, label = label), size = 5, inherit.aes = FALSE, position = position_dodge(0.8))

# Display both plots
grid.arrange(plot_MEred_1, plot_MEred_2, ncol = 1)

#===========================MEturquoise vs adverse outcome==========================================================
# Analysis for MEturquoise

# Convert data to long format and rename to data_long_MEturquoise
data_long_MEturquoise <- merged_data_ID_adverse %>%
  select(RowNames, ExperimentalGroup.x, new_adverse_outcome, 
         Identified_Baseline.data.MEturquoise, Identified_TP36_38weeks.data.MEturquoise, Identified_Postpartum.data.MEturquoise) %>%
  gather(key = "time_point", value = "MEturquoise", 
         Identified_Baseline.data.MEturquoise, Identified_TP36_38weeks.data.MEturquoise, Identified_Postpartum.data.MEturquoise)

# Rename time points for better readability and order them
data_long_MEturquoise$time_point <- factor(data_long_MEturquoise$time_point, 
                                           levels = c("Identified_Baseline.data.MEturquoise", "Identified_TP36_38weeks.data.MEturquoise", "Identified_Postpartum.data.MEturquoise"),
                                           labels = c("Baseline", "36-38 Weeks", "Postpartum"))

# Function to perform Wilcoxon test for adverse outcome and get p-values
get_p_value_adverse_MEturquoise <- function(data, time_point) {
  data_subset <- data %>% filter(time_point == !!time_point)
  if (length(unique(data_subset$new_adverse_outcome)) > 1) {
    test <- wilcox.test(MEturquoise ~ new_adverse_outcome, data = data_subset, paired = FALSE)
    return(test$p.value)
  } else {
    return(NA)
  }
}

# Function to perform Wilcoxon test for intervention/control group and get p-values for each level of adverse outcome
get_p_value_group_MEturquoise <- function(data, time_point) {
  data_subset <- data %>% filter(time_point == !!time_point)
  p_values <- data_subset %>%
    group_by(new_adverse_outcome) %>%
    summarise(p_value = ifelse(length(unique(ExperimentalGroup.x)) > 1, wilcox.test(MEturquoise ~ ExperimentalGroup.x)$p.value, NA))
  return(p_values)
}

# Calculate p-values for each time point and each level of adverse outcome
p_values_group_MEturquoise <- do.call(rbind, lapply(unique(data_long_MEturquoise$time_point), function(tp) {
  get_p_value_group_MEturquoise(data_long_MEturquoise, tp) %>%
    mutate(time_point = tp)
}))

# Calculate p-values for adverse outcome comparisons
p_values_adverse_MEturquoise <- sapply(unique(data_long_MEturquoise$time_point), function(tp) get_p_value_adverse_MEturquoise(data_long_MEturquoise, tp))

# Add significance annotations based on p-values
signif_labels_adverse_MEturquoise <- ifelse(p_values_adverse_MEturquoise < 0.05, "*", "")
signif_labels_group_MEturquoise <- ifelse(p_values_group_MEturquoise$p_value < 0.05, "*", "")

# Create a data frame for annotation positions and labels for adverse outcome comparisons
annotation_df_adverse_MEturquoise <- data.frame(
  time_point = unique(data_long_MEturquoise$time_point),
  x = 1.5,
  y = 0.9,
  label = signif_labels_adverse_MEturquoise
)

# Create a data frame for annotation positions and labels for group comparisons within adverse outcome levels
annotation_df_group_MEturquoise <- data.frame(
  time_point = p_values_group_MEturquoise$time_point,
  new_adverse_outcome = p_values_group_MEturquoise$new_adverse_outcome,
  x = 1.5,
  y = 0.9,
  label = signif_labels_group_MEturquoise
)

# Create the first set of box plots: MEturquoise by Adverse Outcome
plot_MEturquoise_1 <- ggplot(data_long_MEturquoise, aes(x = time_point, y = MEturquoise, fill = new_adverse_outcome)) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(title = "MEturquoise Values by Adverse Outcome and Time Point",
       x = "Time Point",
       y = "MEturquoise Value",
       fill = "Adverse Outcome") +
  theme_minimal() +
  ylim(-0.65, 1) +
  geom_text(data = annotation_df_adverse_MEturquoise, aes(x = time_point, y = y, label = label), size = 5, inherit.aes = FALSE, position = position_dodge(0.8))

# Create the second set of box plots: MEturquoise by Adverse Outcome and Experimental Group
plot_MEturquoise_2 <- ggplot(data_long_MEturquoise, aes(x = interaction(new_adverse_outcome, time_point), y = MEturquoise, fill = ExperimentalGroup.x)) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(title = "MEturquoise Values by Experimental Group within Adverse Outcome and Time Point",
       x = "Adverse Outcome and Time Point",
       y = "MEturquoise Value",
       fill = "Experimental Group") +
  theme_minimal() +
  ylim(-0.65, 1) +
  scale_fill_manual(values = c("Control" = "blue", "Intervention" = "red")) +
  geom_text(data = annotation_df_group_MEturquoise, aes(x = interaction(new_adverse_outcome, time_point), y = y, label = label), size = 5, inherit.aes = FALSE, position = position_dodge(0.8))

# Display both plots
grid.arrange(plot_MEturquoise_1, plot_MEturquoise_2, ncol = 1)

#===========================MElightyellow vs adverse outcome==========================================================

# Analysis for MElightyellow

# Convert data to long format and rename to data_long_MElightyellow
data_long_MElightyellow <- merged_data_ID_adverse %>%
  select(RowNames, ExperimentalGroup.x, new_adverse_outcome, 
         Identified_Baseline.data.MElightyellow, Identified_TP36_38weeks.data.MElightyellow, Identified_Postpartum.data.MElightyellow) %>%
  gather(key = "time_point", value = "MElightyellow", 
         Identified_Baseline.data.MElightyellow, Identified_TP36_38weeks.data.MElightyellow, Identified_Postpartum.data.MElightyellow)

# Rename time points for better readability and order them
data_long_MElightyellow$time_point <- factor(data_long_MElightyellow$time_point, 
                                             levels = c("Identified_Baseline.data.MElightyellow", "Identified_TP36_38weeks.data.MElightyellow", "Identified_Postpartum.data.MElightyellow"),
                                             labels = c("Baseline", "36-38 Weeks", "Postpartum"))

# Function to perform Wilcoxon test for adverse outcome and get p-values
get_p_value_adverse_MElightyellow <- function(data, time_point) {
  data_subset <- data %>% filter(time_point == !!time_point)
  if (length(unique(data_subset$new_adverse_outcome)) > 1) {
    test <- wilcox.test(MElightyellow ~ new_adverse_outcome, data = data_subset, paired = FALSE)
    return(test$p.value)
  } else {
    return(NA)
  }
}

# Function to perform Wilcoxon test for intervention/control group and get p-values for each level of adverse outcome
get_p_value_group_MElightyellow <- function(data, time_point) {
  data_subset <- data %>% filter(time_point == !!time_point)
  p_values <- data_subset %>%
    group_by(new_adverse_outcome) %>%
    summarise(p_value = ifelse(length(unique(ExperimentalGroup.x)) > 1, wilcox.test(MElightyellow ~ ExperimentalGroup.x)$p.value, NA))
  return(p_values)
}

# Calculate p-values for each time point and each level of adverse outcome
p_values_group_MElightyellow <- do.call(rbind, lapply(unique(data_long_MElightyellow$time_point), function(tp) {
  get_p_value_group_MElightyellow(data_long_MElightyellow, tp) %>%
    mutate(time_point = tp)
}))

# Calculate p-values for adverse outcome comparisons
p_values_adverse_MElightyellow <- sapply(unique(data_long_MElightyellow$time_point), function(tp) get_p_value_adverse_MElightyellow(data_long_MElightyellow, tp))

# Add significance annotations based on p-values
signif_labels_adverse_MElightyellow <- ifelse(p_values_adverse_MElightyellow < 0.05, "*", "")
signif_labels_group_MElightyellow <- ifelse(p_values_group_MElightyellow$p_value < 0.05, "*", "")

# Create a data frame for annotation positions and labels for adverse outcome comparisons
annotation_df_adverse_MElightyellow <- data.frame(
  time_point = unique(data_long_MElightyellow$time_point),
  x = 1.5,
  y = 0.9,
  label = signif_labels_adverse_MElightyellow
)

# Create a data frame for annotation positions and labels for group comparisons within adverse outcome levels
annotation_df_group_MElightyellow <- data.frame(
  time_point = p_values_group_MElightyellow$time_point,
  new_adverse_outcome = p_values_group_MElightyellow$new_adverse_outcome,
  x = 1.5,
  y = 0.9,
  label = signif_labels_group_MElightyellow
)

# Create the first set of box plots: MElightyellow by Adverse Outcome
plot_MElightyellow_1 <- ggplot(data_long_MElightyellow, aes(x = time_point, y = MElightyellow, fill = new_adverse_outcome)) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(title = "MElightyellow Values by Adverse Outcome and Time Point",
       x = "Time Point",
       y = "MElightyellow Value",
       fill = "Adverse Outcome") +
  theme_minimal() +
  ylim(-0.5, 1) +
  geom_text(data = annotation_df_adverse_MElightyellow, aes(x = time_point, y = y, label = label), size = 5, inherit.aes = FALSE, position = position_dodge(0.8))

# Create the second set of box plots: MElightyellow by Adverse Outcome and Experimental Group
plot_MElightyellow_2 <- ggplot(data_long_MElightyellow, aes(x = interaction(new_adverse_outcome, time_point), y = MElightyellow, fill = ExperimentalGroup.x)) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(title = "MElightyellow Values by Experimental Group within Adverse Outcome and Time Point",
       x = "Adverse Outcome and Time Point",
       y = "MElightyellow Value",
       fill = "Experimental Group") +
  theme_minimal() +
  ylim(-0.5, 1) +
  scale_fill_manual(values = c("Control" = "blue", "Intervention" = "red")) +
  geom_text(data = annotation_df_group_MElightyellow, aes(x = interaction(new_adverse_outcome, time_point), y = y, label = label), size = 5, inherit.aes = FALSE, position = position_dodge(0.8))

# Display both plots
grid.arrange(plot_MElightyellow_1, plot_MElightyellow_2, ncol = 1)

#===========================MEgrey60 vs adverse outcome==========================================================

# Analysis for MEgrey60

# Convert data to long format and rename to data_long_MEgrey60
data_long_MEgrey60 <- merged_data_ID_adverse %>%
  select(RowNames, ExperimentalGroup.x, new_adverse_outcome, 
         Identified_Baseline.data.MEgrey60, Identified_TP36_38weeks.data.MEgrey60, Identified_Postpartum.data.MEgrey60) %>%
  gather(key = "time_point", value = "MEgrey60", 
         Identified_Baseline.data.MEgrey60, Identified_TP36_38weeks.data.MEgrey60, Identified_Postpartum.data.MEgrey60)

# Rename time points for better readability and order them
data_long_MEgrey60$time_point <- factor(data_long_MEgrey60$time_point, 
                                        levels = c("Identified_Baseline.data.MEgrey60", "Identified_TP36_38weeks.data.MEgrey60", "Identified_Postpartum.data.MEgrey60"),
                                        labels = c("Baseline", "36-38 Weeks", "Postpartum"))

# Function to perform Wilcoxon test for adverse outcome and get p-values
get_p_value_adverse_MEgrey60 <- function(data, time_point) {
  data_subset <- data %>% filter(time_point == !!time_point)
  if (length(unique(data_subset$new_adverse_outcome)) > 1) {
    test <- wilcox.test(MEgrey60 ~ new_adverse_outcome, data = data_subset, paired = FALSE)
    return(test$p.value)
  } else {
    return(NA)
  }
}

# Function to perform Wilcoxon test for intervention/control group and get p-values for each level of adverse outcome
get_p_value_group_MEgrey60 <- function(data, time_point) {
  data_subset <- data %>% filter(time_point == !!time_point)
  p_values <- data_subset %>%
    group_by(new_adverse_outcome) %>%
    summarise(p_value = ifelse(length(unique(ExperimentalGroup.x)) > 1, wilcox.test(MEgrey60 ~ ExperimentalGroup.x)$p.value, NA))
  return(p_values)
}

# Calculate p-values for each time point and each level of adverse outcome
p_values_group_MEgrey60 <- do.call(rbind, lapply(unique(data_long_MEgrey60$time_point), function(tp) {
  get_p_value_group_MEgrey60(data_long_MEgrey60, tp) %>%
    mutate(time_point = tp)
}))

# Calculate p-values for adverse outcome comparisons
p_values_adverse_MEgrey60 <- sapply(unique(data_long_MEgrey60$time_point), function(tp) get_p_value_adverse_MEgrey60(data_long_MEgrey60, tp))

# Add significance annotations based on p-values
signif_labels_adverse_MEgrey60 <- ifelse(p_values_adverse_MEgrey60 < 0.05, "*", "")
signif_labels_group_MEgrey60 <- ifelse(p_values_group_MEgrey60$p_value < 0.05, "*", "")

# Create a data frame for annotation positions and labels for adverse outcome comparisons
annotation_df_adverse_MEgrey60 <- data.frame(
  time_point = unique(data_long_MEgrey60$time_point),
  x = 1.5,
  y = 0.9,
  label = signif_labels_adverse_MEgrey60
)

# Create a data frame for annotation positions and labels for group comparisons within adverse outcome levels
annotation_df_group_MEgrey60 <- data.frame(
  time_point = p_values_group_MEgrey60$time_point,
  new_adverse_outcome = p_values_group_MEgrey60$new_adverse_outcome,
  x = 1.5,
  y = 0.9,
  label = signif_labels_group_MEgrey60
)

# Create the first set of box plots: MEgrey60 by Adverse Outcome
plot_MEgrey60_1 <- ggplot(data_long_MEgrey60, aes(x = time_point, y = MEgrey60, fill = new_adverse_outcome)) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(title = "MEgrey60 Values by Adverse Outcome and Time Point",
       x = "Time Point",
       y = "MEgrey60 Value",
       fill = "Adverse Outcome") +
  theme_minimal() +
  ylim(-0.5, 1) +
  geom_text(data = annotation_df_adverse_MEgrey60, aes(x = time_point, y = y, label = label), size = 5, inherit.aes = FALSE, position = position_dodge(0.8))

# Create the second set of box plots: MEgrey60 by Adverse Outcome and Experimental Group
plot_MEgrey60_2 <- ggplot(data_long_MEgrey60, aes(x = interaction(new_adverse_outcome, time_point), y = MEgrey60, fill = ExperimentalGroup.x)) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(title = "MEgrey60 Values by Experimental Group within Adverse Outcome and Time Point",
       x = "Adverse Outcome and Time Point",
       y = "MEgrey60 Value",
       fill = "Experimental Group") +
  theme_minimal() +
  ylim(-0.5, 1) +
  scale_fill_manual(values = c("Control" = "blue", "Intervention" = "red")) +
  geom_text(data = annotation_df_group_MEgrey60, aes(x = interaction(new_adverse_outcome, time_point), y = y, label = label), size = 5, inherit.aes = FALSE, position = position_dodge(0.8))

# Display both plots
grid.arrange(plot_MEgrey60_1, plot_MEgrey60_2, ncol = 1)


#===========================MEdarkgrey vs adverse outcome==========================================================

# Analysis for MEdarkgrey

# Convert data to long format and rename to data_long_MEdarkgrey
data_long_MEdarkgrey <- merged_data_ID_adverse %>%
  select(RowNames, ExperimentalGroup.x, new_adverse_outcome, 
         Identified_Baseline.data.MEdarkgrey, Identified_TP36_38weeks.data.MEdarkgrey, Identified_Postpartum.data.MEdarkgrey) %>%
  gather(key = "time_point", value = "MEdarkgrey", 
         Identified_Baseline.data.MEdarkgrey, Identified_TP36_38weeks.data.MEdarkgrey, Identified_Postpartum.data.MEdarkgrey)

# Rename time points for better readability and order them
data_long_MEdarkgrey$time_point <- factor(data_long_MEdarkgrey$time_point, 
                                          levels = c("Identified_Baseline.data.MEdarkgrey", "Identified_TP36_38weeks.data.MEdarkgrey", "Identified_Postpartum.data.MEdarkgrey"),
                                          labels = c("Baseline", "36-38 Weeks", "Postpartum"))

# Function to perform Wilcoxon test for adverse outcome and get p-values
get_p_value_adverse_MEdarkgrey <- function(data, time_point) {
  data_subset <- data %>% filter(time_point == !!time_point)
  if (length(unique(data_subset$new_adverse_outcome)) > 1) {
    test <- wilcox.test(MEdarkgrey ~ new_adverse_outcome, data = data_subset, paired = FALSE)
    return(test$p.value)
  } else {
    return(NA)
  }
}

# Function to perform Wilcoxon test for intervention/control group and get p-values for each level of adverse outcome
get_p_value_group_MEdarkgrey <- function(data, time_point) {
  data_subset <- data %>% filter(time_point == !!time_point)
  p_values <- data_subset %>%
    group_by(new_adverse_outcome) %>%
    summarise(p_value = ifelse(length(unique(ExperimentalGroup.x)) > 1, wilcox.test(MEdarkgrey ~ ExperimentalGroup.x)$p.value, NA))
  return(p_values)
}

# Calculate p-values for each time point and each level of adverse outcome
p_values_group_MEdarkgrey <- do.call(rbind, lapply(unique(data_long_MEdarkgrey$time_point), function(tp) {
  get_p_value_group_MEdarkgrey(data_long_MEdarkgrey, tp) %>%
    mutate(time_point = tp)
}))

# Calculate p-values for adverse outcome comparisons
p_values_adverse_MEdarkgrey <- sapply(unique(data_long_MEdarkgrey$time_point), function(tp) get_p_value_adverse_MEdarkgrey(data_long_MEdarkgrey, tp))

# Add significance annotations based on p-values
signif_labels_adverse_MEdarkgrey <- ifelse(p_values_adverse_MEdarkgrey < 0.05, "*", "")
signif_labels_group_MEdarkgrey <- ifelse(p_values_group_MEdarkgrey$p_value < 0.05, "*", "")

# Create a data frame for annotation positions and labels for adverse outcome comparisons
annotation_df_adverse_MEdarkgrey <- data.frame(
  time_point = unique(data_long_MEdarkgrey$time_point),
  x = 1.5,
  y = 0.9,
  label = signif_labels_adverse_MEdarkgrey
)

# Create a data frame for annotation positions and labels for group comparisons within adverse outcome levels
annotation_df_group_MEdarkgrey <- data.frame(
  time_point = p_values_group_MEdarkgrey$time_point,
  new_adverse_outcome = p_values_group_MEdarkgrey$new_adverse_outcome,
  x = 1.5,
  y = 0.9,
  label = signif_labels_group_MEdarkgrey
)

# Create the first set of box plots: MEdarkgrey by Adverse Outcome
plot_MEdarkgrey_1 <- ggplot(data_long_MEdarkgrey, aes(x = time_point, y = MEdarkgrey, fill = new_adverse_outcome)) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(title = "MEdarkgrey Values by Adverse Outcome and Time Point",
       x = "Time Point",
       y = "MEdarkgrey Value",
       fill = "Adverse Outcome") +
  theme_minimal() +
  ylim(-0.9, 1) +
  geom_text(data = annotation_df_adverse_MEdarkgrey, aes(x = time_point, y = y, label = label), size = 5, inherit.aes = FALSE, position = position_dodge(0.8))

# Create the second set of box plots: MEdarkgrey by Adverse Outcome and Experimental Group
plot_MEdarkgrey_2 <- ggplot(data_long_MEdarkgrey, aes(x = interaction(new_adverse_outcome, time_point), y = MEdarkgrey, fill = ExperimentalGroup.x)) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(title = "MEdarkgrey Values by Experimental Group within Adverse Outcome and Time Point",
       x = "Adverse Outcome and Time Point",
       y = "MEdarkgrey Value",
       fill = "Experimental Group") +
  theme_minimal() +
  ylim(-0.9, 1) +
  scale_fill_manual(values = c("Control" = "blue", "Intervention" = "red")) +
  geom_text(data = annotation_df_group_MEdarkgrey, aes(x = interaction(new_adverse_outcome, time_point), y = y, label = label), size = 5, inherit.aes = FALSE, position = position_dodge(0.8))

# Display both plots
grid.arrange(plot_MEdarkgrey_1, plot_MEdarkgrey_2, ncol = 1)

#===========================MEdarkgrey vs adverse outcome==========================================================

# Analysis for MEdarkturquoise

# Convert data to long format and rename to data_long_MEdarkturquoise
data_long_MEdarkturquoise <- merged_data_ID_adverse %>%
  select(RowNames, ExperimentalGroup.x, new_adverse_outcome, 
         Identified_Baseline.data.MEdarkturquoise, Identified_TP36_38weeks.data.MEdarkturquoise, Identified_Postpartum.data.MEdarkturquoise) %>%
  gather(key = "time_point", value = "MEdarkturquoise", 
         Identified_Baseline.data.MEdarkturquoise, Identified_TP36_38weeks.data.MEdarkturquoise, Identified_Postpartum.data.MEdarkturquoise)

# Rename time points for better readability and order them
data_long_MEdarkturquoise$time_point <- factor(data_long_MEdarkturquoise$time_point, 
                                               levels = c("Identified_Baseline.data.MEdarkturquoise", "Identified_TP36_38weeks.data.MEdarkturquoise", "Identified_Postpartum.data.MEdarkturquoise"),
                                               labels = c("Baseline", "36-38 Weeks", "Postpartum"))

# Function to perform Wilcoxon test for adverse outcome and get p-values
get_p_value_adverse_MEdarkturquoise <- function(data, time_point) {
  data_subset <- data %>% filter(time_point == !!time_point)
  if (length(unique(data_subset$new_adverse_outcome)) > 1) {
    test <- wilcox.test(MEdarkturquoise ~ new_adverse_outcome, data = data_subset, paired = FALSE)
    return(test$p.value)
  } else {
    return(NA)
  }
}

# Function to perform Wilcoxon test for intervention/control group and get p-values for each level of adverse outcome
get_p_value_group_MEdarkturquoise <- function(data, time_point) {
  data_subset <- data %>% filter(time_point == !!time_point)
  p_values <- data_subset %>%
    group_by(new_adverse_outcome) %>%
    summarise(p_value = ifelse(length(unique(ExperimentalGroup.x)) > 1, wilcox.test(MEdarkturquoise ~ ExperimentalGroup.x)$p.value, NA))
  return(p_values)
}

# Calculate p-values for each time point and each level of adverse outcome
p_values_group_MEdarkturquoise <- do.call(rbind, lapply(unique(data_long_MEdarkturquoise$time_point), function(tp) {
  get_p_value_group_MEdarkturquoise(data_long_MEdarkturquoise, tp) %>%
    mutate(time_point = tp)
}))

# Calculate p-values for adverse outcome comparisons
p_values_adverse_MEdarkturquoise <- sapply(unique(data_long_MEdarkturquoise$time_point), function(tp) get_p_value_adverse_MEdarkturquoise(data_long_MEdarkturquoise, tp))

# Add significance annotations based on p-values
signif_labels_adverse_MEdarkturquoise <- ifelse(p_values_adverse_MEdarkturquoise < 0.05, "*", "")
signif_labels_group_MEdarkturquoise <- ifelse(p_values_group_MEdarkturquoise$p_value < 0.05, "*", "")

# Create a data frame for annotation positions and labels for adverse outcome comparisons
annotation_df_adverse_MEdarkturquoise <- data.frame(
  time_point = unique(data_long_MEdarkturquoise$time_point),
  x = 1.5,
  y = 0.9,
  label = signif_labels_adverse_MEdarkturquoise
)

# Create a data frame for annotation positions and labels for group comparisons within adverse outcome levels
annotation_df_group_MEdarkturquoise <- data.frame(
  time_point = p_values_group_MEdarkturquoise$time_point,
  new_adverse_outcome = p_values_group_MEdarkturquoise$new_adverse_outcome,
  x = 1.5,
  y = 0.9,
  label = signif_labels_group_MEdarkturquoise
)

# Create the first set of box plots: MEdarkturquoise by Adverse Outcome
plot_MEdarkturquoise_1 <- ggplot(data_long_MEdarkturquoise, aes(x = time_point, y = MEdarkturquoise, fill = new_adverse_outcome)) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(title = "MEdarkturquoise Values by Adverse Outcome and Time Point",
       x = "Time Point",
       y = "MEdarkturquoise Value",
       fill = "Adverse Outcome") +
  theme_minimal() +
  ylim(-0.9, 1) +
  geom_text(data = annotation_df_adverse_MEdarkturquoise, aes(x = time_point, y = y, label = label), size = 5, inherit.aes = FALSE, position = position_dodge(0.8))

# Create the second set of box plots: MEdarkturquoise by Adverse Outcome and Experimental Group
plot_MEdarkturquoise_2 <- ggplot(data_long_MEdarkturquoise, aes(x = interaction(new_adverse_outcome, time_point), y = MEdarkturquoise, fill = ExperimentalGroup.x)) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(title = "MEdarkturquoise Values by Experimental Group within Adverse Outcome and Time Point",
       x = "Adverse Outcome and Time Point",
       y = "MEdarkturquoise Value",
       fill = "Experimental Group") +
  theme_minimal() +
  ylim(-0.9, 1) +
  scale_fill_manual(values = c("Control" = "blue", "Intervention" = "red")) +
  geom_text(data = annotation_df_group_MEdarkturquoise, aes(x = interaction(new_adverse_outcome, time_point), y = y, label = label), size = 5, inherit.aes = FALSE, position = position_dodge(0.8))

# Display both plots
grid.arrange(plot_MEdarkturquoise_1, plot_MEdarkturquoise_2, ncol = 1)

#===========================MEgreenyellow vs adverse outcome==========================================================

# Analysis for MEgreenyellow

# Convert data to long format and rename to data_long_MEgreenyellow
data_long_MEgreenyellow <- merged_data_ID_adverse %>%
  select(RowNames, ExperimentalGroup.x, new_adverse_outcome, 
         Identified_Baseline.data.MEgreenyellow, Identified_TP36_38weeks.data.MEgreenyellow, Identified_Postpartum.data.MEgreenyellow) %>%
  gather(key = "time_point", value = "MEgreenyellow", 
         Identified_Baseline.data.MEgreenyellow, Identified_TP36_38weeks.data.MEgreenyellow, Identified_Postpartum.data.MEgreenyellow)

# Rename time points for better readability and order them
data_long_MEgreenyellow$time_point <- factor(data_long_MEgreenyellow$time_point, 
                                             levels = c("Identified_Baseline.data.MEgreenyellow", "Identified_TP36_38weeks.data.MEgreenyellow", "Identified_Postpartum.data.MEgreenyellow"),
                                             labels = c("Baseline", "36-38 Weeks", "Postpartum"))

# Function to perform Wilcoxon test for adverse outcome and get p-values
get_p_value_adverse_MEgreenyellow <- function(data, time_point) {
  data_subset <- data %>% filter(time_point == !!time_point)
  if (length(unique(data_subset$new_adverse_outcome)) > 1) {
    test <- wilcox.test(MEgreenyellow ~ new_adverse_outcome, data = data_subset, paired = FALSE)
    return(test$p.value)
  } else {
    return(NA)
  }
}

# Function to perform Wilcoxon test for intervention/control group and get p-values for each level of adverse outcome
get_p_value_group_MEgreenyellow <- function(data, time_point) {
  data_subset <- data %>% filter(time_point == !!time_point)
  p_values <- data_subset %>%
    group_by(new_adverse_outcome) %>%
    summarise(p_value = ifelse(length(unique(ExperimentalGroup.x)) > 1, wilcox.test(MEgreenyellow ~ ExperimentalGroup.x)$p.value, NA))
  return(p_values)
}

# Calculate p-values for each time point and each level of adverse outcome
p_values_group_MEgreenyellow <- do.call(rbind, lapply(unique(data_long_MEgreenyellow$time_point), function(tp) {
  get_p_value_group_MEgreenyellow(data_long_MEgreenyellow, tp) %>%
    mutate(time_point = tp)
}))

# Calculate p-values for adverse outcome comparisons
p_values_adverse_MEgreenyellow <- sapply(unique(data_long_MEgreenyellow$time_point), function(tp) get_p_value_adverse_MEgreenyellow(data_long_MEgreenyellow, tp))

# Add significance annotations based on p-values
signif_labels_adverse_MEgreenyellow <- ifelse(p_values_adverse_MEgreenyellow < 0.05, "*", "")
signif_labels_group_MEgreenyellow <- ifelse(p_values_group_MEgreenyellow$p_value < 0.05, "*", "")

# Create a data frame for annotation positions and labels for adverse outcome comparisons
annotation_df_adverse_MEgreenyellow <- data.frame(
  time_point = unique(data_long_MEgreenyellow$time_point),
  x = 1.5,
  y = 0.9,
  label = signif_labels_adverse_MEgreenyellow
)

# Create a data frame for annotation positions and labels for group comparisons within adverse outcome levels
annotation_df_group_MEgreenyellow <- data.frame(
  time_point = p_values_group_MEgreenyellow$time_point,
  new_adverse_outcome = p_values_group_MEgreenyellow$new_adverse_outcome,
  x = 1.5,
  y = 0.9,
  label = signif_labels_group_MEgreenyellow
)

# Create the first set of box plots: MEgreenyellow by Adverse Outcome
plot_MEgreenyellow_1 <- ggplot(data_long_MEgreenyellow, aes(x = time_point, y = MEgreenyellow, fill = new_adverse_outcome)) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(title = "MEgreenyellow Values by Adverse Outcome and Time Point",
       x = "Time Point",
       y = "MEgreenyellow Value",
       fill = "Adverse Outcome") +
  theme_minimal() +
  ylim(-0.9, 1) +
  geom_text(data = annotation_df_adverse_MEgreenyellow, aes(x = time_point, y = y, label = label), size = 5, inherit.aes = FALSE, position = position_dodge(0.8))

# Create the second set of box plots: MEgreenyellow by Adverse Outcome and Experimental Group
plot_MEgreenyellow_2 <- ggplot(data_long_MEgreenyellow, aes(x = interaction(new_adverse_outcome, time_point), y = MEgreenyellow, fill = ExperimentalGroup.x)) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(title = "MEgreenyellow Values by Experimental Group within Adverse Outcome and Time Point",
       x = "Adverse Outcome and Time Point",
       y = "MEgreenyellow Value",
       fill = "Experimental Group") +
  theme_minimal() +
  ylim(-0.9, 1) +
  scale_fill_manual(values = c("Control" = "blue", "Intervention" = "red")) +
  geom_text(data = annotation_df_group_MEgreenyellow, aes(x = interaction(new_adverse_outcome, time_point), y = y, label = label), size = 5, inherit.aes = FALSE, position = position_dodge(0.8))

# Display both plots
grid.arrange(plot_MEgreenyellow_1, plot_MEgreenyellow_2, ncol = 1)

#===========================MEgreen vs adverse outcome==========================================================

# Analysis for MEgreen

# Convert data to long format and rename to data_long_MEgreen
data_long_MEgreen <- merged_data_ID_adverse %>%
  select(RowNames, ExperimentalGroup.x, new_adverse_outcome, 
         Identified_Baseline.data.MEgreen, Identified_TP36_38weeks.data.MEgreen, Identified_Postpartum.data.MEgreen) %>%
  gather(key = "time_point", value = "MEgreen", 
         Identified_Baseline.data.MEgreen, Identified_TP36_38weeks.data.MEgreen, Identified_Postpartum.data.MEgreen)

# Rename time points for better readability and order them
data_long_MEgreen$time_point <- factor(data_long_MEgreen$time_point, 
                                       levels = c("Identified_Baseline.data.MEgreen", "Identified_TP36_38weeks.data.MEgreen", "Identified_Postpartum.data.MEgreen"),
                                       labels = c("Baseline", "36-38 Weeks", "Postpartum"))

# Function to perform Wilcoxon test for adverse outcome and get p-values
get_p_value_adverse_MEgreen <- function(data, time_point) {
  data_subset <- data %>% filter(time_point == !!time_point)
  if (length(unique(data_subset$new_adverse_outcome)) > 1) {
    test <- wilcox.test(MEgreen ~ new_adverse_outcome, data = data_subset, paired = FALSE)
    return(test$p.value)
  } else {
    return(NA)
  }
}

# Function to perform Wilcoxon test for intervention/control group and get p-values for each level of adverse outcome
get_p_value_group_MEgreen <- function(data, time_point) {
  data_subset <- data %>% filter(time_point == !!time_point)
  p_values <- data_subset %>%
    group_by(new_adverse_outcome) %>%
    summarise(p_value = ifelse(length(unique(ExperimentalGroup.x)) > 1, wilcox.test(MEgreen ~ ExperimentalGroup.x)$p.value, NA))
  return(p_values)
}

# Calculate p-values for each time point and each level of adverse outcome
p_values_group_MEgreen <- do.call(rbind, lapply(unique(data_long_MEgreen$time_point), function(tp) {
  get_p_value_group_MEgreen(data_long_MEgreen, tp) %>%
    mutate(time_point = tp)
}))

# Calculate p-values for adverse outcome comparisons
p_values_adverse_MEgreen <- sapply(unique(data_long_MEgreen$time_point), function(tp) get_p_value_adverse_MEgreen(data_long_MEgreen, tp))

# Add significance annotations based on p-values
signif_labels_adverse_MEgreen <- ifelse(p_values_adverse_MEgreen < 0.05, "*", "")
signif_labels_group_MEgreen <- ifelse(p_values_group_MEgreen$p_value < 0.05, "*", "")

# Create a data frame for annotation positions and labels for adverse outcome comparisons
annotation_df_adverse_MEgreen <- data.frame(
  time_point = unique(data_long_MEgreen$time_point),
  x = 1.5,
  y = 0.9,
  label = signif_labels_adverse_MEgreen
)

# Create a data frame for annotation positions and labels for group comparisons within adverse outcome levels
annotation_df_group_MEgreen <- data.frame(
  time_point = p_values_group_MEgreen$time_point,
  new_adverse_outcome = p_values_group_MEgreen$new_adverse_outcome,
  x = 1.5,
  y = 0.9,
  label = signif_labels_group_MEgreen
)

# Create the first set of box plots: MEgreen by Adverse Outcome
plot_MEgreen_1 <- ggplot(data_long_MEgreen, aes(x = time_point, y = MEgreen, fill = new_adverse_outcome)) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(title = "MEgreen Values by Adverse Outcome and Time Point",
       x = "Time Point",
       y = "MEgreen Value",
       fill = "Adverse Outcome") +
  theme_minimal() +
  ylim(-0.9, 1) +
  geom_text(data = annotation_df_adverse_MEgreen, aes(x = time_point, y = y, label = label), size = 5, inherit.aes = FALSE, position = position_dodge(0.8))

# Create the second set of box plots: MEgreen by Adverse Outcome and Experimental Group
plot_MEgreen_2 <- ggplot(data_long_MEgreen, aes(x = interaction(new_adverse_outcome, time_point), y = MEgreen, fill = ExperimentalGroup.x)) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(title = "MEgreen Values by Experimental Group within Adverse Outcome and Time Point",
       x = "Adverse Outcome and Time Point",
       y = "MEgreen Value",
       fill = "Experimental Group") +
  theme_minimal() +
  ylim(-0.9, 1) +
  scale_fill_manual(values = c("Control" = "blue", "Intervention" = "red")) +
  geom_text(data = annotation_df_group_MEgreen, aes(x = interaction(new_adverse_outcome, time_point), y = y, label = label), size = 5, inherit.aes = FALSE, position = position_dodge(0.8))

# Display both plots
grid.arrange(plot_MEgreen_1, plot_MEgreen_2, ncol = 1)

#===========================MEmagenta vs adverse outcome==========================================================


# Analysis for MEmagenta

# Convert data to long format and rename to data_long_MEmagenta
data_long_MEmagenta <- merged_data_ID_adverse %>%
  select(RowNames, ExperimentalGroup.x, new_adverse_outcome, 
         Identified_Baseline.data.MEmagenta, Identified_TP36_38weeks.data.MEmagenta, Identified_Postpartum.data.MEmagenta) %>%
  gather(key = "time_point", value = "MEmagenta", 
         Identified_Baseline.data.MEmagenta, Identified_TP36_38weeks.data.MEmagenta, Identified_Postpartum.data.MEmagenta)

# Rename time points for better readability and order them
data_long_MEmagenta$time_point <- factor(data_long_MEmagenta$time_point, 
                                         levels = c("Identified_Baseline.data.MEmagenta", "Identified_TP36_38weeks.data.MEmagenta", "Identified_Postpartum.data.MEmagenta"),
                                         labels = c("Baseline", "36-38 Weeks", "Postpartum"))

# Function to perform Wilcoxon test for adverse outcome and get p-values
get_p_value_adverse_MEmagenta <- function(data, time_point) {
  data_subset <- data %>% filter(time_point == !!time_point)
  if (length(unique(data_subset$new_adverse_outcome)) > 1) {
    test <- wilcox.test(MEmagenta ~ new_adverse_outcome, data = data_subset, paired = FALSE)
    return(test$p.value)
  } else {
    return(NA)
  }
}

# Function to perform Wilcoxon test for intervention/control group and get p-values for each level of adverse outcome
get_p_value_group_MEmagenta <- function(data, time_point) {
  data_subset <- data %>% filter(time_point == !!time_point)
  p_values <- data_subset %>%
    group_by(new_adverse_outcome) %>%
    summarise(p_value = ifelse(length(unique(ExperimentalGroup.x)) > 1, wilcox.test(MEmagenta ~ ExperimentalGroup.x)$p.value, NA))
  return(p_values)
}

# Calculate p-values for each time point and each level of adverse outcome
p_values_group_MEmagenta <- do.call(rbind, lapply(unique(data_long_MEmagenta$time_point), function(tp) {
  get_p_value_group_MEmagenta(data_long_MEmagenta, tp) %>%
    mutate(time_point = tp)
}))

# Calculate p-values for adverse outcome comparisons
p_values_adverse_MEmagenta <- sapply(unique(data_long_MEmagenta$time_point), function(tp) get_p_value_adverse_MEmagenta(data_long_MEmagenta, tp))

# Add significance annotations based on p-values
signif_labels_adverse_MEmagenta <- ifelse(p_values_adverse_MEmagenta < 0.05, "*", "")
signif_labels_group_MEmagenta <- ifelse(p_values_group_MEmagenta$p_value < 0.05, "*", "")

# Create a data frame for annotation positions and labels for adverse outcome comparisons
annotation_df_adverse_MEmagenta <- data.frame(
  time_point = unique(data_long_MEmagenta$time_point),
  x = 1.5,
  y = 0.9,
  label = signif_labels_adverse_MEmagenta
)

# Create a data frame for annotation positions and labels for group comparisons within adverse outcome levels
annotation_df_group_MEmagenta <- data.frame(
  time_point = p_values_group_MEmagenta$time_point,
  new_adverse_outcome = p_values_group_MEmagenta$new_adverse_outcome,
  x = 1.5,
  y = 0.9,
  label = signif_labels_group_MEmagenta
)

# Create the first set of box plots: MEmagenta by Adverse Outcome
plot_MEmagenta_1 <- ggplot(data_long_MEmagenta, aes(x = time_point, y = MEmagenta, fill = new_adverse_outcome)) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(title = "MEmagenta Values by Adverse Outcome and Time Point",
       x = "Time Point",
       y = "MEmagenta Value",
       fill = "Adverse Outcome") +
  theme_minimal() +
  ylim(-0.9, 1) +
  geom_text(data = annotation_df_adverse_MEmagenta, aes(x = time_point, y = y, label = label), size = 5, inherit.aes = FALSE, position = position_dodge(0.8))

# Create the second set of box plots: MEmagenta by Adverse Outcome and Experimental Group
plot_MEmagenta_2 <- ggplot(data_long_MEmagenta, aes(x = interaction(new_adverse_outcome, time_point), y = MEmagenta, fill = ExperimentalGroup.x)) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(title = "MEmagenta Values by Experimental Group within Adverse Outcome and Time Point",
       x = "Adverse Outcome and Time Point",
       y = "MEmagenta Value",
       fill = "Experimental Group") +
  theme_minimal() +
  ylim(-0.9, 1) +
  scale_fill_manual(values = c("Control" = "blue", "Intervention" = "red")) +
  geom_text(data = annotation_df_group_MEmagenta, aes(x = interaction(new_adverse_outcome, time_point), y = y, label = label), size = 5, inherit.aes = FALSE, position = position_dodge(0.8))

# Display both plots
grid.arrange(plot_MEmagenta_1, plot_MEmagenta_2, ncol = 1)

#===========================MElightyellow vs preterm adverse outcome==========================================================


# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# Convert group to a factor
merged_data_ID_adverse$ExperimentalGroup.x <- factor(merged_data_ID_adverse$ExperimentalGroup.x, levels = c(0, 1), labels = c("Control", "Intervention"))

# Convert preterm_adverse_outcome to a factor with levels "No" and "Yes"
merged_data_ID_adverse$preterm_adverse_outcome <- factor(merged_data_ID_adverse$preterm_adverse_outcome, levels = c(0, 1), labels = c("No", "Yes"))

# Convert data to long format and rename to data_long_MElightyellow
data_long_MElightyellow <- merged_data_ID_adverse %>%
  select(RowNames, ExperimentalGroup.x, preterm_adverse_outcome, 
         Identified_Baseline.data.MElightyellow, Identified_TP36_38weeks.data.MElightyellow, Identified_Postpartum.data.MElightyellow) %>%
  gather(key = "time_point", value = "MElightyellow", 
         Identified_Baseline.data.MElightyellow, Identified_TP36_38weeks.data.MElightyellow, Identified_Postpartum.data.MElightyellow)

# Rename time points for better readability and order them
data_long_MElightyellow$time_point <- factor(data_long_MElightyellow$time_point, 
                                             levels = c("Identified_Baseline.data.MElightyellow", "Identified_TP36_38weeks.data.MElightyellow", "Identified_Postpartum.data.MElightyellow"),
                                             labels = c("Baseline", "36-38 Weeks", "Postpartum"))

# Function to perform Wilcoxon test for preterm_adverse_outcome and get p-values
get_p_value_adverse_MElightyellow <- function(data, time_point) {
  data_subset <- data %>% filter(time_point == !!time_point)
  if (length(unique(data_subset$preterm_adverse_outcome)) > 1) {
    test <- wilcox.test(MElightyellow ~ preterm_adverse_outcome, data = data_subset, paired = FALSE)
    return(test$p.value)
  } else {
    return(NA)
  }
}

# Function to perform Wilcoxon test for intervention/control group and get p-values for each level of preterm_adverse_outcome
get_p_value_group_MElightyellow <- function(data, time_point) {
  data_subset <- data %>% filter(time_point == !!time_point)
  p_values <- data_subset %>%
    group_by(preterm_adverse_outcome) %>%
    summarise(p_value = ifelse(length(unique(ExperimentalGroup.x)) > 1, wilcox.test(MElightyellow ~ ExperimentalGroup.x)$p.value, NA))
  return(p_values)
}

# Calculate p-values for each time point and each level of preterm_adverse_outcome
p_values_group_MElightyellow <- do.call(rbind, lapply(unique(data_long_MElightyellow$time_point), function(tp) {
  get_p_value_group_MElightyellow(data_long_MElightyellow, tp) %>%
    mutate(time_point = tp)
}))

# Calculate p-values for preterm_adverse_outcome comparisons
p_values_adverse_MElightyellow <- sapply(unique(data_long_MElightyellow$time_point), function(tp) get_p_value_adverse_MElightyellow(data_long_MElightyellow, tp))

# Add significance annotations based on p-values
signif_labels_adverse_MElightyellow <- ifelse(p_values_adverse_MElightyellow < 0.05, "*", "")
signif_labels_group_MElightyellow <- ifelse(p_values_group_MElightyellow$p_value < 0.05, "*", "")

# Create a data frame for annotation positions and labels for preterm_adverse_outcome comparisons
annotation_df_adverse_MElightyellow <- data.frame(
  time_point = unique(data_long_MElightyellow$time_point),
  x = 1.5,
  y = 0.9,
  label = signif_labels_adverse_MElightyellow
)

# Create a data frame for annotation positions and labels for group comparisons within preterm_adverse_outcome levels
annotation_df_group_MElightyellow <- data.frame(
  time_point = p_values_group_MElightyellow$time_point,
  preterm_adverse_outcome = p_values_group_MElightyellow$preterm_adverse_outcome,
  x = 1.5,
  y = 0.9,
  label = signif_labels_group_MElightyellow
)

# Create the first set of box plots: MElightyellow by Preterm Adverse Outcome
plot_MElightyellow_1 <- ggplot(data_long_MElightyellow, aes(x = time_point, y = MElightyellow, fill = preterm_adverse_outcome)) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(title = "MElightyellow Values by Preterm Adverse Outcome and Time Point",
       x = "Time Point",
       y = "MElightyellow Value",
       fill = "Preterm Adverse Outcome") +
  theme_minimal() +
  ylim(-0.5, 1) +
  geom_text(data = annotation_df_adverse_MElightyellow, aes(x = time_point, y = y, label = label), size = 5, inherit.aes = FALSE, position = position_dodge(0.8))

# Create the second set of box plots: MElightyellow by Preterm Adverse Outcome and Experimental Group
plot_MElightyellow_2 <- ggplot(data_long_MElightyellow, aes(x = interaction(preterm_adverse_outcome, time_point), y = MElightyellow, fill = ExperimentalGroup.x)) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(title = "MElightyellow Values by Experimental Group within Preterm Adverse Outcome and Time Point",
       x = "Preterm Adverse Outcome and Time Point",
       y = "MElightyellow Value",
       fill = "Experimental Group") +
  theme_minimal() +
  ylim(-0.5, 1) +
  scale_fill_manual(values = c("Control" = "blue", "Intervention" = "red")) +
  geom_text(data = annotation_df_group_MElightyellow, aes(x = interaction(preterm_adverse_outcome, time_point), y = y, label = label), size = 5, inherit.aes = FALSE, position = position_dodge(0.8))

# Display both plots
grid.arrange(plot_MElightyellow_1, plot_MElightyellow_2, ncol = 1)

#===========================MEgreenyellow vs preterm adverse outcome==========================================================

# Analysis for MEgreenyellow

# Convert data to long format and rename to data_long_MEgreenyellow
data_long_MEgreenyellow <- merged_data_ID_adverse %>%
  select(RowNames, ExperimentalGroup.x, preterm_adverse_outcome, 
         Identified_Baseline.data.MEgreenyellow, Identified_TP36_38weeks.data.MEgreenyellow, Identified_Postpartum.data.MEgreenyellow) %>%
  gather(key = "time_point", value = "MEgreenyellow", 
         Identified_Baseline.data.MEgreenyellow, Identified_TP36_38weeks.data.MEgreenyellow, Identified_Postpartum.data.MEgreenyellow)

# Rename time points for better readability and order them
data_long_MEgreenyellow$time_point <- factor(data_long_MEgreenyellow$time_point, 
                                             levels = c("Identified_Baseline.data.MEgreenyellow", "Identified_TP36_38weeks.data.MEgreenyellow", "Identified_Postpartum.data.MEgreenyellow"),
                                             labels = c("Baseline", "36-38 Weeks", "Postpartum"))

# Function to perform Wilcoxon test for preterm_adverse_outcome and get p-values
get_p_value_adverse_MEgreenyellow <- function(data, time_point) {
  data_subset <- data %>% filter(time_point == !!time_point)
  if (length(unique(data_subset$preterm_adverse_outcome)) > 1) {
    test <- wilcox.test(MEgreenyellow ~ preterm_adverse_outcome, data = data_subset, paired = FALSE)
    return(test$p.value)
  } else {
    return(NA)
  }
}

# Function to perform Wilcoxon test for intervention/control group and get p-values for each level of preterm_adverse_outcome
get_p_value_group_MEgreenyellow <- function(data, time_point) {
  data_subset <- data %>% filter(time_point == !!time_point)
  p_values <- data_subset %>%
    group_by(preterm_adverse_outcome) %>%
    summarise(p_value = ifelse(length(unique(ExperimentalGroup.x)) > 1, wilcox.test(MEgreenyellow ~ ExperimentalGroup.x)$p.value, NA))
  return(p_values)
}

# Calculate p-values for each time point and each level of preterm_adverse_outcome
p_values_group_MEgreenyellow <- do.call(rbind, lapply(unique(data_long_MEgreenyellow$time_point), function(tp) {
  get_p_value_group_MEgreenyellow(data_long_MEgreenyellow, tp) %>%
    mutate(time_point = tp)
}))

# Calculate p-values for preterm_adverse_outcome comparisons
p_values_adverse_MEgreenyellow <- sapply(unique(data_long_MEgreenyellow$time_point), function(tp) get_p_value_adverse_MEgreenyellow(data_long_MEgreenyellow, tp))

# Add significance annotations based on p-values
signif_labels_adverse_MEgreenyellow <- ifelse(p_values_adverse_MEgreenyellow < 0.05, "*", "")
signif_labels_group_MEgreenyellow <- ifelse(p_values_group_MEgreenyellow$p_value < 0.05, "*", "")

# Create a data frame for annotation positions and labels for preterm_adverse_outcome comparisons
annotation_df_adverse_MEgreenyellow <- data.frame(
  time_point = unique(data_long_MEgreenyellow$time_point),
  x = 1.5,
  y = 0.9,
  label = signif_labels_adverse_MEgreenyellow
)

# Create a data frame for annotation positions and labels for group comparisons within preterm_adverse_outcome levels
annotation_df_group_MEgreenyellow <- data.frame(
  time_point = p_values_group_MEgreenyellow$time_point,
  preterm_adverse_outcome = p_values_group_MEgreenyellow$preterm_adverse_outcome,
  x = 1.5,
  y = 0.9,
  label = signif_labels_group_MEgreenyellow
)

# Create the first set of box plots: MEgreenyellow by Preterm Adverse Outcome
plot_MEgreenyellow_1 <- ggplot(data_long_MEgreenyellow, aes(x = time_point, y = MEgreenyellow, fill = preterm_adverse_outcome)) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(title = "MEgreenyellow Values by Preterm Adverse Outcome and Time Point",
       x = "Time Point",
       y = "MEgreenyellow Value",
       fill = "Preterm Adverse Outcome") +
  theme_minimal() +
  ylim(-0.8, 1) +
  geom_text(data = annotation_df_adverse_MEgreenyellow, aes(x = time_point, y = y, label = label), size = 5, inherit.aes = FALSE, position = position_dodge(0.8))

# Create the second set of box plots: MEgreenyellow by Preterm Adverse Outcome and Experimental Group
plot_MEgreenyellow_2 <- ggplot(data_long_MEgreenyellow, aes(x = interaction(preterm_adverse_outcome, time_point), y = MEgreenyellow, fill = ExperimentalGroup.x)) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(title = "MEgreenyellow Values by Experimental Group within Preterm Adverse Outcome and Time Point",
       x = "Preterm Adverse Outcome and Time Point",
       y = "MEgreenyellow Value",
       fill = "Experimental Group") +
  theme_minimal() +
  ylim(-0.8, 1) +
  scale_fill_manual(values = c("Control" = "blue", "Intervention" = "red")) +
  geom_text(data = annotation_df_group_MEgreenyellow, aes(x = interaction(preterm_adverse_outcome, time_point), y = y, label = label), size = 5, inherit.aes = FALSE, position = position_dodge(0.8))

# Display both plots
grid.arrange(plot_MEgreenyellow_1, plot_MEgreenyellow_2, ncol = 1)

#===========================MEmidnightblue vs preterm adverse outcome==========================================================

# Analysis for MEmidnightblue

# Convert data to long format and rename to data_long_MEmidnightblue
data_long_MEmidnightblue <- merged_data_ID_adverse %>%
  select(RowNames, ExperimentalGroup.x, preterm_adverse_outcome, 
         Identified_Baseline.data.MEmidnightblue, Identified_TP36_38weeks.data.MEmidnightblue, Identified_Postpartum.data.MEmidnightblue) %>%
  gather(key = "time_point", value = "MEmidnightblue", 
         Identified_Baseline.data.MEmidnightblue, Identified_TP36_38weeks.data.MEmidnightblue, Identified_Postpartum.data.MEmidnightblue)

# Rename time points for better readability and order them
data_long_MEmidnightblue$time_point <- factor(data_long_MEmidnightblue$time_point, 
                                              levels = c("Identified_Baseline.data.MEmidnightblue", "Identified_TP36_38weeks.data.MEmidnightblue", "Identified_Postpartum.data.MEmidnightblue"),
                                              labels = c("Baseline", "36-38 Weeks", "Postpartum"))

# Function to perform Wilcoxon test for preterm_adverse_outcome and get p-values
get_p_value_adverse_MEmidnightblue <- function(data, time_point) {
  data_subset <- data %>% filter(time_point == !!time_point)
  if (length(unique(data_subset$preterm_adverse_outcome)) > 1) {
    test <- wilcox.test(MEmidnightblue ~ preterm_adverse_outcome, data = data_subset, paired = FALSE)
    return(test$p.value)
  } else {
    return(NA)
  }
}

# Function to perform Wilcoxon test for intervention/control group and get p-values for each level of preterm_adverse_outcome
get_p_value_group_MEmidnightblue <- function(data, time_point) {
  data_subset <- data %>% filter(time_point == !!time_point)
  p_values <- data_subset %>%
    group_by(preterm_adverse_outcome) %>%
    summarise(p_value = ifelse(length(unique(ExperimentalGroup.x)) > 1, wilcox.test(MEmidnightblue ~ ExperimentalGroup.x)$p.value, NA))
  return(p_values)
}

# Calculate p-values for each time point and each level of preterm_adverse_outcome
p_values_group_MEmidnightblue <- do.call(rbind, lapply(unique(data_long_MEmidnightblue$time_point), function(tp) {
  get_p_value_group_MEmidnightblue(data_long_MEmidnightblue, tp) %>%
    mutate(time_point = tp)
}))

# Calculate p-values for preterm_adverse_outcome comparisons
p_values_adverse_MEmidnightblue <- sapply(unique(data_long_MEmidnightblue$time_point), function(tp) get_p_value_adverse_MEmidnightblue(data_long_MEmidnightblue, tp))

# Add significance annotations based on p-values
signif_labels_adverse_MEmidnightblue <- ifelse(p_values_adverse_MEmidnightblue < 0.05, "*", "")
signif_labels_group_MEmidnightblue <- ifelse(p_values_group_MEmidnightblue$p_value < 0.05, "*", "")

# Create a data frame for annotation positions and labels for preterm_adverse_outcome comparisons
annotation_df_adverse_MEmidnightblue <- data.frame(
  time_point = unique(data_long_MEmidnightblue$time_point),
  x = 1.5,
  y = 0.9,
  label = signif_labels_adverse_MEmidnightblue
)

# Create a data frame for annotation positions and labels for group comparisons within preterm_adverse_outcome levels
annotation_df_group_MEmidnightblue <- data.frame(
  time_point = p_values_group_MEmidnightblue$time_point,
  preterm_adverse_outcome = p_values_group_MEmidnightblue$preterm_adverse_outcome,
  x = 1.5,
  y = 0.9,
  label = signif_labels_group_MEmidnightblue
)

# Create the first set of box plots: MEmidnightblue by Preterm Adverse Outcome
plot_MEmidnightblue_1 <- ggplot(data_long_MEmidnightblue, aes(x = time_point, y = MEmidnightblue, fill = preterm_adverse_outcome)) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(title = "MEmidnightblue Values by Preterm Adverse Outcome and Time Point",
       x = "Time Point",
       y = "MEmidnightblue Value",
       fill = "Preterm Adverse Outcome") +
  theme_minimal() +
  ylim(-0.8, 1) +
  geom_text(data = annotation_df_adverse_MEmidnightblue, aes(x = time_point, y = y, label = label), size = 5, inherit.aes = FALSE, position = position_dodge(0.8))

# Create the second set of box plots: MEmidnightblue by Preterm Adverse Outcome and Experimental Group
plot_MEmidnightblue_2 <- ggplot(data_long_MEmidnightblue, aes(x = interaction(preterm_adverse_outcome, time_point), y = MEmidnightblue, fill = ExperimentalGroup.x)) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(title = "MEmidnightblue Values by Experimental Group within Preterm Adverse Outcome and Time Point",
       x = "Preterm Adverse Outcome and Time Point",
       y = "MEmidnightblue Value",
       fill = "Experimental Group") +
  theme_minimal() +
  ylim(-0.8, 1) +
  scale_fill_manual(values = c("Control" = "blue", "Intervention" = "red")) +
  geom_text(data = annotation_df_group_MEmidnightblue, aes(x = interaction(preterm_adverse_outcome, time_point), y = y, label = label), size = 5, inherit.aes = FALSE, position = position_dodge(0.8))

# Display both plots
grid.arrange(plot_MEmidnightblue_1, plot_MEmidnightblue_2, ncol = 1)










#===========================Old chunk of code==========================================================

