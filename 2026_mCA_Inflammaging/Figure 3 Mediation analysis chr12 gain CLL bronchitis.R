# Load required libraries
library(dplyr)        # Load dplyr first to avoid conflicts
library(mediation)
library(survival)
library(survminer)
library(ggplot2)
library(readxl)
library(data.table)
library(lubridate)
library(tidyr)
library(stringr)

# Ensure dplyr functions are used (in case of conflicts)
select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate

# ==========================================
# LOAD AND PREPARE YOUR UKBB DATA FOR CHR12 → CLL → CHRONIC BRONCHITIS ANALYSIS
# WITH EXACT SAME FILTERS AS COMPETING RISK ANALYSIS
# ==========================================

cat("Loading UKBB data for chr12 gain → CLL → chronic bronchitis analysis...\n")

# Set censor date
censor_date <- as.Date("2023-07-14")

# ===========================
# NEW: Load analyzed samples (SAME AS COMPETING RISK CODE)
# ===========================
analyzed <- read_excel("~/Desktop/UK BB/ukb4777.samples_analyzed.xlsx") %>%
  mutate(ID = as.character(ID))

cat("Total analyzed samples:", nrow(analyzed), "\n")

# ===========================
# Load mCA/cohort data with STRICTER DATE FILTERING
# ===========================
mca_info_raw <- read_excel("~/Desktop/UK BB/Age sex dob collection date.xlsx") %>%
  mutate(
    ID = as.character(ID),
    mCA_date = mdy(`Date of mCA assessment`)
  ) %>%
  filter(ID %in% analyzed$ID)  # FILTER TO ANALYZED SAMPLES ONLY

# Check how many invalid dates before filtering
cat("\n=== Checking mCA assessment dates ===\n")
cat("Total rows before filtering:", nrow(mca_info_raw), "\n")
cat("Missing mCA dates:", sum(is.na(mca_info_raw$mCA_date)), "\n")
cat("Dates = 1900-01-01:", sum(mca_info_raw$mCA_date == as.Date("1900-01-01"), na.rm = TRUE), "\n")
cat("Dates <= 1901-01-01:", sum(mca_info_raw$mCA_date <= as.Date("1901-01-01"), na.rm = TRUE), "\n")

# Filter out invalid dates (STRICTER FILTER: > 1901-01-01)
mca_info <- mca_info_raw %>%
  filter(!is.na(mCA_date),                           # Exclude missing dates
         mCA_date != as.Date("1900-01-01"),         # Exclude 1900-01-01 explicitly
         mCA_date > as.Date("1901-01-01")) %>%      # STRICTER: Exclude dates <= 1901-01-01
  select(ID, mCA_date, Sex, `Age at blood collection`) %>%
  rename(age = `Age at blood collection`) %>%
  filter(!is.na(age), !is.na(Sex)) %>%
  mutate(
    sex = case_when(
      Sex == 0 ~ "Female",
      Sex == 1 ~ "Male",
      TRUE ~ NA_character_
    ),
    sex = factor(sex, levels = c("Female","Male"))
  ) %>%
  filter(!is.na(sex)) %>%
  select(ID, mCA_date, age, sex) %>%
  distinct(ID, .keep_all = TRUE)  # Ensure unique IDs

cat("\nAfter filtering:\n")
cat("Valid participants:", nrow(mca_info), "\n")
cat("Date range:", as.character(min(mca_info$mCA_date)), "to", as.character(max(mca_info$mCA_date)), "\n")

# ===========================
# Load CNV data and determine mCA status (SAME LOGIC AS COMPETING RISK CODE)
# ===========================
cnv_all <- read_excel("~/Desktop/UK BB/CNV_data_19808.xlsx") %>%
  mutate(ID = as.character(ID)) %>%
  filter(ID %in% analyzed$ID)  # Filter to analyzed samples

# Summarize CNV status for each patient
cnv_summary <- cnv_all %>%
  group_by(ID) %>%
  summarise(
    has_known_cnv = any(COPY_CHANGE %in% c("gain", "loss", "neutral")),
    has_only_unknown = all(COPY_CHANGE == "unknown"),
    .groups = "drop"
  )

# Determine mCA status (SAME LOGIC AS COMPETING RISK CODE)
cnv_status <- analyzed %>%
  left_join(cnv_summary, by = "ID") %>%
  mutate(
    has_known_cnv = ifelse(is.na(has_known_cnv), FALSE, has_known_cnv),
    has_only_unknown = ifelse(is.na(has_only_unknown), FALSE, has_only_unknown),
    mCA_status = case_when(
      has_known_cnv ~ "mCA",
      !has_known_cnv & !has_only_unknown ~ "no_mCA",
      has_only_unknown ~ "unknown_only"
    )
  ) %>%
  filter(mCA_status %in% c("mCA", "no_mCA")) %>%  # EXCLUDE patients with only "unknown"
  select(ID, mCA_status)

cat("\n=== mCA Status Distribution ===\n")
print(table(cnv_status$mCA_status))
cat("Total patients with definitive mCA results:", nrow(cnv_status), "\n")

# Merge mCA status with mca_info
mca_info <- mca_info %>%
  inner_join(cnv_status, by = "ID")

cat("\nFinal cohort after excluding 'unknown' CNV patients:", nrow(mca_info), "\n")
cat("mCA status distribution:\n")
print(table(mca_info$mCA_status))

# Create chr12_gain variable
chr12_gain_ids <- cnv_all %>%
  filter(CHR == 12, COPY_CHANGE == "gain") %>%
  distinct(ID) %>%
  pull(ID)

# Add chr12_gain indicator
mca_info <- mca_info %>%
  mutate(
    chr12_gain = ifelse(ID %in% chr12_gain_ids, 1, 0),
    No_mCA = ifelse(mCA_status == "no_mCA", 1, 0)  # 1 = no mCA, 0 = has some mCA
  )

# ===========================
# Load diagnosis data
# ===========================
diag_data <- fread("~/Desktop/UK BB/Diagnosis code date.tsv", sep = "\t") %>%
  mutate(
    patient_id = as.character(patient_id),
    admission_date = as.Date(admission_date),
    diagnosis = str_replace_all(toupper(str_trim(as.character(diagnosis))), "\\.", "")
  ) %>%
  filter(!is.na(admission_date), !is.na(diagnosis), !is.na(patient_id))

# ===========================
# Calculate baseline disease burden (SAME AS COMPETING RISK CODE)
# ===========================
cat("\n=== Calculating baseline disease burden ===\n")

# Convert to data.table for speed
setDT(diag_data)
mca_info_dt <- as.data.table(mca_info)
setkey(diag_data, patient_id)
setkey(mca_info_dt, ID)

# Filter diagnosis data to only include patients in mca_info FIRST
cat("Filtering diagnosis data to study participants...\n")
diag_data_filtered <- diag_data[patient_id %in% mca_info_dt$ID]
cat("Filtered diagnosis records:", nrow(diag_data_filtered), "\n")

# Join with mCA dates
cat("Joining with mCA dates...\n")
diag_with_mca <- diag_data_filtered[mca_info_dt[, .(ID, mCA_date)], 
                                    on = .(patient_id = ID), 
                                    nomatch = 0]

# Calculate baseline burden
cat("Calculating baseline disease counts...\n")
baseline_disease_burden <- diag_with_mca[admission_date <= mCA_date, 
                                         .(n_baseline_diseases = uniqueN(diagnosis),
                                           first_disease_date = min(admission_date)),
                                         by = patient_id
]

cat("Baseline calculations complete for", nrow(baseline_disease_burden), "patients\n")

# Convert back to data.frame for compatibility
baseline_disease_burden <- as.data.frame(baseline_disease_burden)

# Add to mca_info
mca_info <- mca_info %>%
  left_join(baseline_disease_burden, by = c("ID" = "patient_id")) %>%
  mutate(
    n_baseline_diseases = replace_na(n_baseline_diseases, 0),
    time_from_first_disease = if_else(
      !is.na(first_disease_date),
      as.numeric(difftime(mCA_date, first_disease_date, units = "days")) / 365.25,
      NA_real_
    )
  )

cat("Baseline disease burden calculated\n")
cat("Mean baseline diseases - chr12_gain:", 
    round(mean(mca_info$n_baseline_diseases[mca_info$chr12_gain == 1], na.rm = TRUE), 2), "\n")
cat("Mean baseline diseases - No_mCA:", 
    round(mean(mca_info$n_baseline_diseases[mca_info$No_mCA == 1], na.rm = TRUE), 2), "\n\n")

# Compare baseline disease burden between groups
baseline_comparison <- mca_info %>%
  filter(chr12_gain == 1 | No_mCA == 1) %>%
  mutate(group = ifelse(chr12_gain == 1, "chr12_gain", "No_mCA")) %>%
  group_by(group) %>%
  summarise(
    n = n(),
    mean_baseline = mean(n_baseline_diseases),
    sd_baseline = sd(n_baseline_diseases),
    median_baseline = median(n_baseline_diseases),
    mean_time_from_first_dx = mean(time_from_first_disease, na.rm = TRUE),
    sd_time_from_first_dx = sd(time_from_first_disease, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== Baseline Disease Burden Comparison ===\n")
print(baseline_comparison)

# ===========================
# CLL - INCIDENT CASES ONLY (EXCLUDE PRE-EXISTING CODES)
# ===========================
cat("\n=== Identifying INCIDENT CLL cases ===\n")

# Pre-existing CLL codes (C911 and all subcodes)
pre_existing_cll <- diag_data %>%
  filter(grepl("^C911", diagnosis)) %>%
  inner_join(mca_info %>% select(ID, mCA_date), 
             by = c("patient_id" = "ID"),
             relationship = "many-to-one") %>%
  filter(admission_date <= mCA_date) %>%
  distinct(patient_id, diagnosis)

cat("Patients with pre-existing CLL codes:", n_distinct(pre_existing_cll$patient_id), "\n")

# Post-mCA CLL diagnoses
post_mca_cll <- diag_data %>%
  filter(grepl("^C911", diagnosis)) %>%
  inner_join(mca_info %>% select(ID, mCA_date), 
             by = c("patient_id" = "ID"),
             relationship = "many-to-one") %>%
  filter(admission_date > mCA_date)

cat("Post-mCA CLL diagnoses:", nrow(post_mca_cll), "\n")

# Exclude pre-existing CLL codes (anti_join removes matching codes)
incident_cll <- post_mca_cll %>%
  anti_join(pre_existing_cll, by = c("patient_id", "diagnosis")) %>%
  group_by(patient_id) %>%
  summarise(first_cll_date = min(admission_date), .groups = "drop")

cat("INCIDENT CLL cases (new codes only):", nrow(incident_cll), "patients\n")

# ===========================
# Chronic Bronchitis - INCIDENT CASES ONLY (EXCLUDE PRE-EXISTING CODES)
# ===========================
cat("\n=== Identifying INCIDENT Chronic Bronchitis cases ===\n")

# Pre-existing chronic bronchitis codes (J44, J431, J439)
pre_existing_cb <- diag_data %>%
  filter(grepl("^J44|^J431$|^J439$", diagnosis)) %>%
  inner_join(mca_info %>% select(ID, mCA_date), 
             by = c("patient_id" = "ID"),
             relationship = "many-to-one") %>%
  filter(admission_date <= mCA_date) %>%
  distinct(patient_id, diagnosis)

cat("Patients with pre-existing chronic bronchitis codes:", n_distinct(pre_existing_cb$patient_id), "\n")

# Post-mCA chronic bronchitis diagnoses
post_mca_cb <- diag_data %>%
  filter(grepl("^J44|^J431$|^J439$", diagnosis)) %>%
  inner_join(mca_info %>% select(ID, mCA_date), 
             by = c("patient_id" = "ID"),
             relationship = "many-to-one") %>%
  filter(admission_date > mCA_date)

cat("Post-mCA chronic bronchitis diagnoses:", nrow(post_mca_cb), "\n")

# Exclude pre-existing chronic bronchitis codes
incident_cb <- post_mca_cb %>%
  anti_join(pre_existing_cb, by = c("patient_id", "diagnosis")) %>%
  group_by(patient_id) %>%
  summarise(first_cb_date = min(admission_date), .groups = "drop")

cat("INCIDENT chronic bronchitis cases (new codes only):", nrow(incident_cb), "patients\n")

# ===========================
# Combine all data
# ===========================
cohort <- mca_info %>%
  left_join(incident_cll, by = c("ID" = "patient_id")) %>%
  left_join(incident_cb, by = c("ID" = "patient_id")) %>%
  mutate(
    # CLL status - INCIDENT cases only (excludes pre-existing codes)
    cll_status = ifelse(!is.na(first_cll_date), 1, 0),
    # Chronic bronchitis status - INCIDENT cases only (excludes pre-existing codes)
    cb_status = ifelse(!is.na(first_cb_date), 1, 0)
  )

cat("\n=== Cohort Summary ===\n")
cat("Total cohort:", nrow(cohort), "participants\n")
cat("chr12_gain cases:", sum(cohort$chr12_gain), "\n")
cat("No_mCA cases:", sum(cohort$No_mCA), "\n")
cat("INCIDENT CLL cases:", sum(cohort$cll_status), "\n") 
cat("INCIDENT chronic bronchitis cases:", sum(cohort$cb_status), "\n\n")

# ==========================================
# TIME-TO-EVENT MEDIATION ANALYSIS
# ==========================================

# Create survival times for analysis
analysis_data <- cohort %>%
  filter(chr12_gain == 1 | No_mCA == 1) %>%  # Only chr12_gain carriers + No_mCA patients
  filter(complete.cases(chr12_gain, cll_status, cb_status, age, sex, n_baseline_diseases)) %>%
  mutate(
    # Time to CLL (in years from mCA date)
    time_to_cll = ifelse(cll_status == 1, 
                         as.numeric(first_cll_date - mCA_date)/365.25,
                         as.numeric(censor_date - mCA_date)/365.25),
    cll_event = cll_status,
    
    # Time to chronic bronchitis (in years from mCA date)
    time_to_cb = ifelse(cb_status == 1,
                        as.numeric(first_cb_date - mCA_date)/365.25,
                        as.numeric(censor_date - mCA_date)/365.25),
    cb_event = cb_status,
    
    # Ensure positive follow-up times
    time_to_cll = pmax(time_to_cll, 0.001),
    time_to_cb = pmax(time_to_cb, 0.001)
  )

cat("\n=== Analysis Dataset ===\n")
cat("Total participants:", nrow(analysis_data), "\n")
cat("chr12_gain carriers:", sum(analysis_data$chr12_gain), "\n")
cat("No_mCA controls:", sum(analysis_data$No_mCA), "\n")
cat("INCIDENT CLL events:", sum(analysis_data$cll_event), "\n")
cat("INCIDENT chronic bronchitis events:", sum(analysis_data$cb_event), "\n")
cat("Median follow-up time:", round(median(analysis_data$time_to_cb), 2), "years\n\n")

# ==========================================
# METHOD 1: SEQUENTIAL COX MODELS
# ==========================================

cat("=== METHOD 1: SEQUENTIAL COX MODELS ===\n")

# Step 1: chr12_gain -> time to CLL
cat("Step 1: chr12_gain -> CLL (Cox model, adjusted for age, sex, baseline disease burden)\n")
cll_cox <- coxph(Surv(time_to_cll, cll_event) ~ chr12_gain + age + sex + n_baseline_diseases, 
                 data = analysis_data)
print(summary(cll_cox))

# Step 2: chr12_gain + CLL status -> time to chronic bronchitis
cat("\nStep 2: chr12_gain + CLL -> chronic bronchitis (Cox model, adjusted for age, sex, baseline disease burden)\n")
cb_cox <- coxph(Surv(time_to_cb, cb_event) ~ chr12_gain + cll_status + age + sex + n_baseline_diseases, 
                data = analysis_data)
print(summary(cb_cox))

# Extract hazard ratios and confidence intervals
cll_hr <- exp(coef(cll_cox)["chr12_gain"])
cll_ci <- exp(confint(cll_cox)["chr12_gain", ])
cll_p <- summary(cll_cox)$coefficients["chr12_gain", "Pr(>|z|)"]

cb_direct_hr <- exp(coef(cb_cox)["chr12_gain"])
cb_direct_ci <- exp(confint(cb_cox)["chr12_gain", ])
cb_direct_p <- summary(cb_cox)$coefficients["chr12_gain", "Pr(>|z|)"]

cll_to_cb_hr <- exp(coef(cb_cox)["cll_status"])
cll_to_cb_ci <- exp(confint(cb_cox)["cll_status", ])
cll_to_cb_p <- summary(cb_cox)$coefficients["cll_status", "Pr(>|z|)"]

# Total effect model (chr12_gain -> chronic bronchitis without CLL adjustment)
cat("\nTotal Effect: chr12_gain -> chronic bronchitis (adjusted for age, sex, baseline disease burden, without CLL adjustment)\n")
cb_total_cox <- coxph(Surv(time_to_cb, cb_event) ~ chr12_gain + age + sex + n_baseline_diseases, 
                      data = analysis_data)
print(summary(cb_total_cox))

total_hr <- exp(coef(cb_total_cox)["chr12_gain"])
total_ci <- exp(confint(cb_total_cox)["chr12_gain", ])
total_p <- summary(cb_total_cox)$coefficients["chr12_gain", "Pr(>|z|)"]

# ==========================================
# SUMMARY RESULTS TABLE
# ==========================================

cat("\n=== TIME-TO-EVENT MEDIATION RESULTS SUMMARY ===\n")

# Create results table
survival_results <- data.frame(
  Effect = c("Total Effect (chr12_gain -> chronic bronchitis)",
             "Direct Effect (chr12_gain -> chronic bronchitis, adjusted for CLL)", 
             "Pathway A (chr12_gain -> CLL)",
             "Pathway B (CLL -> chronic bronchitis)"),
  Hazard_Ratio = c(round(total_hr, 3),
                   round(cb_direct_hr, 3),
                   round(cll_hr, 3),
                   round(cll_to_cb_hr, 3)),
  CI_Lower = c(round(total_ci[1], 3),
               round(cb_direct_ci[1], 3), 
               round(cll_ci[1], 3),
               round(cll_to_cb_ci[1], 3)),
  CI_Upper = c(round(total_ci[2], 3),
               round(cb_direct_ci[2], 3),
               round(cll_ci[2], 3), 
               round(cll_to_cb_ci[2], 3)),
  P_Value = c(ifelse(total_p < 0.001, "< 0.001", round(total_p, 3)),
              ifelse(cb_direct_p < 0.001, "< 0.001", round(cb_direct_p, 3)),
              ifelse(cll_p < 0.001, "< 0.001", round(cll_p, 3)),
              ifelse(cll_to_cb_p < 0.001, "< 0.001", round(cll_to_cb_p, 3))),
  Significance = c(ifelse(total_p < 0.05, "Significant", "Not Significant"),
                   ifelse(cb_direct_p < 0.05, "Significant", "Not Significant"),
                   ifelse(cll_p < 0.05, "Significant", "Not Significant"), 
                   ifelse(cll_to_cb_p < 0.05, "Significant", "Not Significant"))
)

print(survival_results, row.names = FALSE)

# Save results
write.csv(survival_results, "~/Desktop/UK BB/chr12_cll_cb_survival_mediation_results_matched_filters.csv", row.names = FALSE)

# ==========================================
# BINARY MEDIATION ANALYSIS (for comparison)
# ==========================================

cat("\n=== BINARY MEDIATION ANALYSIS (for comparison) ===\n")

# Step 1: Mediator model (chr12_gain -> CLL)
cat("Step 1: chr12_gain -> CLL (logistic model, adjusted for age, sex, baseline disease burden)\n")
mediator_model <- glm(cll_status ~ chr12_gain + age + sex + n_baseline_diseases, 
                      data = analysis_data, 
                      family = binomial())
print(summary(mediator_model))

# Step 2: Outcome model (chr12_gain + CLL -> chronic bronchitis)
cat("\nStep 2: chr12_gain + CLL -> chronic bronchitis (logistic model, adjusted for age, sex, baseline disease burden)\n")
outcome_model <- glm(cb_status ~ chr12_gain + cll_status + age + sex + n_baseline_diseases, 
                     data = analysis_data, 
                     family = binomial())
print(summary(outcome_model))

# Step 3: Mediation analysis
cat("\nStep 3: Binary mediation analysis\n")
set.seed(123)
mediation_results <- mediate(
  model.m = mediator_model,      # mediator model
  model.y = outcome_model,       # outcome model
  treat = "chr12_gain",          # treatment variable
  mediator = "cll_status",       # mediator variable
  robustSE = TRUE,               # robust standard errors
  sims = 1000                    # bootstrap simulations
)

# View results
summary(mediation_results)

# Plot results
plot(mediation_results)

# Save plot as PDF (square format with title and axis label)
pdf("~/Desktop/UK BB/chr12_cll_cb_mediation_plot_matched_filters.pdf", width = 8, height = 8)
plot(mediation_results, 
     main = "Mediation of chr12 gain Association with Chronic Bronchitis through CLL Development",
     xlab = "Absolute Risk Increase (percentage points)")
dev.off()

cat("Plot saved as: ~/Desktop/UK BB/chr12_cll_cb_mediation_plot_matched_filters.pdf\n")

# ==========================================
# CALCULATE BASELINE RISKS AND CREATE RISK TABLE
# ==========================================

# Calculate baseline chronic bronchitis risk in No_mCA group
baseline_cb_risk <- analysis_data %>%
  filter(No_mCA == 1) %>%
  summarise(baseline_risk = mean(cb_status) * 100) %>%
  pull(baseline_risk)

# Calculate actual chr12_gain chronic bronchitis risk
chr12_cb_risk <- analysis_data %>%
  filter(chr12_gain == 1) %>%
  summarise(chr12_risk = mean(cb_status) * 100) %>%
  pull(chr12_risk)

# Extract mediation effects (convert to percentage points)
total_effect_pct <- mediation_results$tau.coef * 100
acme_effect_pct <- mediation_results$d0.coef * 100  # Average ACME
ade_effect_pct <- mediation_results$z0.coef * 100   # Average ADE

# Create risk table
risk_table <- data.frame(
  Group = c("No_mCA (baseline)", 
            "chr12_gain (observed)", 
            "chr12_gain (through CLL pathway)",
            "chr12_gain (direct pathway)"),
  CB_Risk_Percent = c(
    round(baseline_cb_risk, 2),
    round(chr12_cb_risk, 2),
    round(baseline_cb_risk + acme_effect_pct, 2),
    round(baseline_cb_risk + ade_effect_pct, 2)
  ),
  Risk_Difference = c(
    "Reference",
    paste0("+", round(chr12_cb_risk - baseline_cb_risk, 2), "%"),
    paste0("+", round(acme_effect_pct, 2), "% (through CLL)"),
    paste0("+", round(ade_effect_pct, 2), "% (direct)")
  ),
  P_Value = c(
    "-",
    "-",
    ifelse(mediation_results$d0.p < 0.001, "< 0.001", round(mediation_results$d0.p, 3)),
    ifelse(mediation_results$z0.p < 0.001, "< 0.001", round(mediation_results$z0.p, 3))
  )
)

cat("\n=== CHRONIC BRONCHITIS RISK TABLE (chr12 gain → CLL → CB) ===\n")
print(risk_table, row.names = FALSE)

cat("\n=== SUMMARY ===\n")
cat("Baseline chronic bronchitis risk (No_mCA):", round(baseline_cb_risk, 2), "%\n")
cat("Observed chr12_gain chronic bronchitis risk:", round(chr12_cb_risk, 2), "%\n")
cat("Total risk difference:", round(chr12_cb_risk - baseline_cb_risk, 2), "percentage points\n")
cat("Mediated through CLL:", round(acme_effect_pct, 2), "percentage points\n")
cat("Direct effect:", round(ade_effect_pct, 2), "percentage points\n")
if (!is.null(mediation_results$n0.coef)) {
  cat("Proportion mediated:", round(mediation_results$n0.coef * 100, 1), "% of total effect\n")
}

# Save risk table
write.csv(risk_table, "~/Desktop/UK BB/chr12_cll_cb_risk_table_matched_filters.csv", row.names = FALSE)
cat("\nRisk table saved as: ~/Desktop/UK BB/chr12_cll_cb_risk_table_matched_filters.csv\n")

# ==========================================
# MEDIATION INTERPRETATION
# ==========================================

cat("\n=== MEDIATION INTERPRETATION ===\n")

if(cll_p < 0.05 && cll_to_cb_p < 0.05) {
  if(cb_direct_p >= 0.05) {
    cat("RESULT: Full mediation - chr12_gain -> chronic bronchitis association is fully mediated through CLL\n")
  } else {
    cat("RESULT: Partial mediation - chr12_gain -> chronic bronchitis has both direct and CLL-mediated pathways\n")
  }
} else if(cll_p < 0.05 && cll_to_cb_p >= 0.05) {
  cat("RESULT: No mediation - CLL does not significantly affect chronic bronchitis risk\n")
} else if(cll_p >= 0.05) {
  cat("RESULT: No mediation - chr12_gain does not significantly affect CLL risk\n") 
} else {
  cat("RESULT: Complex pattern - see individual pathway significance\n")
}

cat("chr12_gain -> CLL pathway:", ifelse(cll_p < 0.05, "SIGNIFICANT", "NOT SIGNIFICANT"), "\n")
cat("CLL -> chronic bronchitis pathway:", ifelse(cll_to_cb_p < 0.05, "SIGNIFICANT", "NOT SIGNIFICANT"), "\n")
cat("chr12_gain -> chronic bronchitis total effect:", ifelse(total_p < 0.05, "SIGNIFICANT", "NOT SIGNIFICANT"), "\n")
cat("chr12_gain -> chronic bronchitis direct effect:", ifelse(cb_direct_p < 0.05, "SIGNIFICANT", "NOT SIGNIFICANT"), "\n")

cat("\n=== FILES SAVED ===\n")
cat("1. chr12_cll_cb_mediation_plot_matched_filters.pdf - Mediation plot\n")
cat("2. chr12_cll_cb_risk_table_matched_filters.csv - Risk comparison table\n") 
cat("3. chr12_cll_cb_survival_mediation_results_matched_filters.csv - Cox model results\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("\n=== KEY DIFFERENCES FROM ORIGINAL CODE ===\n")
cat("1. Filtered to analyzed samples only (ukb4777.samples_analyzed.xlsx)\n")
cat("2. Stricter date filtering: excluded dates <= 1901-01-01\n")
cat("3. Excluded patients with only 'unknown' CNV status\n")
cat("4. INCIDENT disease definition: excluded pre-existing diagnosis CODES\n")
cat("   (not just patients with any prior diagnosis)\n")
cat("5. Adjusted for baseline disease burden (in addition to age and sex)\n")
cat("6. This matches the competing risk analysis methodology exactly\n")
