# ==========================================
# COMPLETE ANALYSIS: CHR9 NEUTRAL → MPN → CARDIOVASCULAR DISEASE
# WITH BUILT-IN MODEL DIAGNOSTICS AND FIGURE GENERATION
# ==========================================

library(dplyr)
library(mediation)
library(survival)
library(survminer)
library(ggplot2)
library(readxl)
library(data.table)
library(lubridate)
library(tidyr)
library(stringr)

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate

cat("Loading UKBB data for chr9 neutral → MPN → CVD analysis...\n")

censor_date <- as.Date("2023-07-14")

# ===========================
# Load analyzed samples
# ===========================
analyzed <- read_excel("~/Desktop/UK BB/ukb4777.samples_analyzed.xlsx") %>%
  mutate(ID = as.character(ID))

cat("Total analyzed samples:", nrow(analyzed), "\n")

# ===========================
# Load mCA/cohort data
# ===========================
mca_info_raw <- read_excel("~/Desktop/UK BB/Age sex dob collection date.xlsx") %>%
  mutate(
    ID = as.character(ID),
    mCA_date = mdy(`Date of mCA assessment`)
  ) %>%
  filter(ID %in% analyzed$ID)

cat("\n=== Checking mCA assessment dates ===\n")
cat("Total rows before filtering:", nrow(mca_info_raw), "\n")
cat("Missing mCA dates:", sum(is.na(mca_info_raw$mCA_date)), "\n")
cat("Dates = 1900-01-01:", sum(mca_info_raw$mCA_date == as.Date("1900-01-01"), na.rm = TRUE), "\n")
cat("Dates <= 1901-01-01:", sum(mca_info_raw$mCA_date <= as.Date("1901-01-01"), na.rm = TRUE), "\n")

mca_info <- mca_info_raw %>%
  filter(!is.na(mCA_date),
         mCA_date != as.Date("1900-01-01"),
         mCA_date > as.Date("1901-01-01")) %>%
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
  distinct(ID, .keep_all = TRUE)

cat("\nAfter filtering:\n")
cat("Valid participants:", nrow(mca_info), "\n")
cat("Date range:", as.character(min(mca_info$mCA_date)), "to", as.character(max(mca_info$mCA_date)), "\n")

# ===========================
# Load CNV data
# ===========================
cnv_all <- read_excel("~/Desktop/UK BB/CNV_data_19808.xlsx") %>%
  mutate(ID = as.character(ID)) %>%
  filter(ID %in% analyzed$ID)

cnv_summary <- cnv_all %>%
  group_by(ID) %>%
  summarise(
    has_known_cnv = any(COPY_CHANGE %in% c("gain", "loss", "neutral")),
    has_only_unknown = all(COPY_CHANGE == "unknown"),
    .groups = "drop"
  )

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
  filter(mCA_status %in% c("mCA", "no_mCA")) %>%
  select(ID, mCA_status)

cat("\n=== mCA Status Distribution ===\n")
print(table(cnv_status$mCA_status))

mca_info <- mca_info %>%
  inner_join(cnv_status, by = "ID")

cat("\nFinal cohort after excluding 'unknown' CNV patients:", nrow(mca_info), "\n")

# ===========================
# Identify chr9 neutral carriers
# ===========================
chr9_neutral_ids <- cnv_all %>%
  filter(CHR == 9, COPY_CHANGE == "neutral") %>%
  distinct(ID) %>%
  pull(ID)

mca_info <- mca_info %>%
  mutate(
    chr9_neutral = ifelse(ID %in% chr9_neutral_ids, 1, 0),
    No_mCA = ifelse(mCA_status == "no_mCA", 1, 0)
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
# DIAGNOSTIC: Check MPN codes in database
# ===========================
cat("\n=== DIAGNOSTIC: Checking MPN codes in database ===\n")
mpn_codes <- diag_data %>%
  filter(grepl("^D45|^D473|^D474", diagnosis)) %>%
  group_by(diagnosis) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n))

cat("MPN codes found in database:\n")
print(mpn_codes)
cat("\n")

# ===========================
# Calculate baseline disease burden
# ===========================
cat("\n=== Calculating baseline disease burden ===\n")

setDT(diag_data)
mca_info_dt <- as.data.table(mca_info)
setkey(diag_data, patient_id)
setkey(mca_info_dt, ID)

diag_data_filtered <- diag_data[patient_id %in% mca_info_dt$ID]
diag_with_mca <- diag_data_filtered[mca_info_dt[, .(ID, mCA_date)], 
                                    on = .(patient_id = ID), 
                                    nomatch = 0]

baseline_disease_burden <- diag_with_mca[admission_date <= mCA_date, 
                                         .(n_baseline_diseases = uniqueN(diagnosis),
                                           first_disease_date = min(admission_date)),
                                         by = patient_id
]

baseline_disease_burden <- as.data.frame(baseline_disease_burden)

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
cat("Mean baseline diseases - chr9_neutral:", 
    round(mean(mca_info$n_baseline_diseases[mca_info$chr9_neutral == 1], na.rm = TRUE), 2), "\n")
cat("Mean baseline diseases - No_mCA:", 
    round(mean(mca_info$n_baseline_diseases[mca_info$No_mCA == 1], na.rm = TRUE), 2), "\n\n")

# Compare baseline disease burden between groups
baseline_comparison <- mca_info %>%
  filter(chr9_neutral == 1 | No_mCA == 1) %>%
  mutate(group = ifelse(chr9_neutral == 1, "chr9_neutral", "No_mCA")) %>%
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
# MPN - INCIDENT CASES
# ===========================
cat("\n=== Identifying INCIDENT MPN cases ===\n")

# MPN ICD-10 codes
# D45 - Polycythemia vera
# D473 - Essential thrombocythemia
# D474 - Primary myelofibrosis
mpn_pattern <- "^D45$|^D473$|^D474$"

cat("\n=== MPN ICD-10 Codes ===\n")
cat("D45 - Polycythemia vera\n")
cat("D473 - Essential thrombocythemia\n")
cat("D474 - Primary myelofibrosis\n\n")

# Show what MPN codes will be matched
mpn_codes_matched <- diag_data %>%
  filter(grepl(mpn_pattern, diagnosis)) %>%
  group_by(diagnosis) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n))

cat("MPN codes that will be captured by pattern:\n")
print(mpn_codes_matched)
cat("\n")

pre_existing_mpn <- diag_data %>%
  filter(grepl(mpn_pattern, diagnosis)) %>%
  inner_join(mca_info %>% select(ID, mCA_date), 
             by = c("patient_id" = "ID"),
             relationship = "many-to-one") %>%
  filter(admission_date <= mCA_date) %>%
  distinct(patient_id, diagnosis)

cat("Patients with pre-existing MPN codes:", n_distinct(pre_existing_mpn$patient_id), "\n")

post_mca_mpn <- diag_data %>%
  filter(grepl(mpn_pattern, diagnosis)) %>%
  inner_join(mca_info %>% select(ID, mCA_date), 
             by = c("patient_id" = "ID"),
             relationship = "many-to-one") %>%
  filter(admission_date > mCA_date)

cat("Post-mCA MPN diagnoses:", nrow(post_mca_mpn), "\n")

incident_mpn <- post_mca_mpn %>%
  anti_join(pre_existing_mpn, by = c("patient_id", "diagnosis")) %>%
  group_by(patient_id) %>%
  summarise(first_mpn_date = min(admission_date), .groups = "drop")

cat("INCIDENT MPN cases (new codes only):", nrow(incident_mpn), "patients\n")

# ===========================
# CARDIOVASCULAR DISEASE - INCIDENT CASES
# ===========================
cat("\n=== Identifying INCIDENT CVD cases ===\n")

# CVD ICD-10 codes
# I50 - Heart failure
# I63 - Cerebral infarction (stroke)
# I70 - Atherosclerosis
# I702-I707 - Atherosclerosis of specific arteries
# I10 - Essential hypertension
cvd_pattern <- "^I50$|^I63$|^I70$|^I10$|^I702$|^I703$|^I704$|^I705$|^I706$|^I707$"

cat("\n=== CVD ICD-10 Codes ===\n")
cat("I50 - Heart failure\n")
cat("I63 - Cerebral infarction (stroke)\n")
cat("I70 - Atherosclerosis\n")
cat("I702-I707 - Atherosclerosis of specific arteries\n")
cat("I10 - Essential hypertension\n\n")

pre_existing_cvd <- diag_data %>%
  filter(grepl(cvd_pattern, diagnosis)) %>%
  inner_join(mca_info %>% select(ID, mCA_date), 
             by = c("patient_id" = "ID"),
             relationship = "many-to-one") %>%
  filter(admission_date <= mCA_date) %>%
  distinct(patient_id, diagnosis)

cat("Patients with pre-existing CVD codes:", n_distinct(pre_existing_cvd$patient_id), "\n")

post_mca_cvd <- diag_data %>%
  filter(grepl(cvd_pattern, diagnosis)) %>%
  inner_join(mca_info %>% select(ID, mCA_date), 
             by = c("patient_id" = "ID"),
             relationship = "many-to-one") %>%
  filter(admission_date > mCA_date)

cat("Post-mCA CVD diagnoses:", nrow(post_mca_cvd), "\n")

incident_cvd <- post_mca_cvd %>%
  anti_join(pre_existing_cvd, by = c("patient_id", "diagnosis")) %>%
  group_by(patient_id) %>%
  summarise(first_cvd_date = min(admission_date), .groups = "drop")

cat("INCIDENT CVD cases (new codes only):", nrow(incident_cvd), "patients\n")

# ===========================
# Combine all data
# ===========================
cohort <- mca_info %>%
  left_join(incident_mpn, by = c("ID" = "patient_id")) %>%
  left_join(incident_cvd, by = c("ID" = "patient_id")) %>%
  mutate(
    mpn_status = ifelse(!is.na(first_mpn_date), 1, 0),
    cvd_status = ifelse(!is.na(first_cvd_date), 1, 0)
  )

cat("\n=== Cohort Summary ===\n")
cat("Total cohort:", nrow(cohort), "participants\n")
cat("chr9_neutral cases:", sum(cohort$chr9_neutral), "\n")
cat("No_mCA cases:", sum(cohort$No_mCA), "\n")
cat("INCIDENT MPN cases:", sum(cohort$mpn_status), "\n")
cat("INCIDENT CVD cases:", sum(cohort$cvd_status), "\n\n")

# ===========================
# Check if we have enough MPN cases to proceed
# ===========================
if (sum(cohort$mpn_status) == 0) {
  cat("\n")
  cat(rep("!", 80), "\n", sep = "")
  cat("CRITICAL ERROR: No MPN cases found!\n")
  cat("Cannot proceed with mediation analysis.\n")
  cat("Please check ICD-10 coding in your database.\n")
  cat(rep("!", 80), "\n", sep = "")
  stop("No MPN cases found - analysis cannot proceed")
}

# ===========================
# Create analysis dataset
# ===========================
analysis_data <- cohort %>%
  filter(chr9_neutral == 1 | No_mCA == 1) %>%
  filter(complete.cases(chr9_neutral, mpn_status, cvd_status, age, sex, n_baseline_diseases)) %>%
  mutate(
    time_to_mpn = ifelse(mpn_status == 1, 
                         as.numeric(first_mpn_date - mCA_date)/365.25,
                         as.numeric(censor_date - mCA_date)/365.25),
    mpn_event = mpn_status,
    
    time_to_cvd = ifelse(cvd_status == 1,
                         as.numeric(first_cvd_date - mCA_date)/365.25,
                         as.numeric(censor_date - mCA_date)/365.25),
    cvd_event = cvd_status,
    
    time_to_mpn = pmax(time_to_mpn, 0.001),
    time_to_cvd = pmax(time_to_cvd, 0.001),
    
    # Create factor for plotting
    chr9_neutral_factor = factor(chr9_neutral, levels = c(0, 1), labels = c("No mCA", "9=")),
    mpn_factor = factor(mpn_status, levels = c(0, 1), labels = c("No MPN", "MPN"))
  )

cat("=== Analysis Dataset ===\n")
cat("Total participants:", nrow(analysis_data), "\n")
cat("chr9_neutral carriers:", sum(analysis_data$chr9_neutral), "\n")
cat("No_mCA controls:", sum(analysis_data$No_mCA), "\n")
cat("INCIDENT MPN events:", sum(analysis_data$mpn_event), "\n")
cat("INCIDENT CVD events:", sum(analysis_data$cvd_event), "\n")
cat("Median follow-up time:", round(median(analysis_data$time_to_cvd), 2), "years\n\n")

# ===========================
# POWER WARNING CHECK
# ===========================
cat("\n=== POWER CHECK ===\n")
if (sum(analysis_data$chr9_neutral) < 100) {
  cat("⚠️  WARNING: Small exposure group (n=", sum(analysis_data$chr9_neutral), ")\n")
  cat("   Results may be underpowered and unstable.\n")
  cat("   Consider combining with other chromosomal abnormalities or\n")
  cat("   interpreting results as exploratory.\n\n")
}

mpn_in_exposed <- sum(analysis_data$mpn_event[analysis_data$chr9_neutral == 1])
cat("MPN cases in chr9_neutral group:", mpn_in_exposed, "\n")
if (mpn_in_exposed < 5) {
  cat("⚠️  WARNING: Very few MPN cases in exposed group.\n")
  cat("   Pathway A estimates may be unstable.\n\n")
}

# ==========================================
# BUILT-IN DIAGNOSTICS BEFORE MEDIATION
# ==========================================

cat("\n")
cat(rep("=", 80), "\n", sep = "")
cat("MODEL DIAGNOSTICS - CHR9 NEUTRAL → MPN → CVD\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("=== 1. CONTINGENCY TABLES (Checking for Separation) ===\n\n")

cat("A. chr9_neutral vs MPN events:\n")
tbl1 <- table(analysis_data$chr9_neutral, analysis_data$mpn_event,
              dnn = c("chr9_neutral", "MPN_event"))
print(tbl1)

if (any(tbl1 == 0)) {
  cat("❌ WARNING: ZERO cells detected! Complete separation!\n")
} else if (any(tbl1 > 0 & tbl1 < 5)) {
  cat("⚠️  WARNING: Small cells (<5) detected. Quasi-separation possible.\n")
} else {
  cat("✓ OK: No separation issues\n")
}

cat("\nB. chr9_neutral vs CVD events:\n")
tbl2 <- table(analysis_data$chr9_neutral, analysis_data$cvd_event,
              dnn = c("chr9_neutral", "CVD_event"))
print(tbl2)

if (any(tbl2 == 0)) {
  cat("❌ WARNING: ZERO cells detected! Complete separation!\n")
} else if (any(tbl2 > 0 & tbl2 < 5)) {
  cat("⚠️  WARNING: Small cells (<5) detected. Quasi-separation possible.\n")
} else {
  cat("✓ OK: No separation issues\n")
}

cat("\nC. MPN vs CVD (among chr9_neutral=1 only) - PATH B CHECK:\n")
exposed <- analysis_data[analysis_data$chr9_neutral == 1, ]
tbl3 <- table(exposed$mpn_event, exposed$cvd_event,
              dnn = c("MPN_event", "CVD_event"))
print(tbl3)

if (any(tbl3 == 0)) {
  cat("❌ WARNING: ZERO cells among exposed! Path B is unstable!\n")
} else if (any(tbl3 > 0 & tbl3 < 5)) {
  cat("⚠️  WARNING: Small cells (<5) among exposed. Path B may be unstable.\n")
} else {
  cat("✓ OK: No separation issues for Path B\n")
}

cat("\n")

# ==========================================
# COX MODELS (METHOD 1: SEQUENTIAL COX MODELS)
# ==========================================

cat("=== 2. COX MODEL RESULTS (METHOD 1: SEQUENTIAL COX MODELS) ===\n\n")

# Pathway A: chr9_neutral -> MPN
cat("Step 1: chr9_neutral → MPN (Cox model, adjusted for age, sex, baseline disease burden)\n")
mpn_cox <- coxph(Surv(time_to_mpn, mpn_event) ~ chr9_neutral + age + sex + n_baseline_diseases, 
                 data = analysis_data)
print(summary(mpn_cox))

# Total Effect
cat("\nTotal Effect: chr9_neutral → CVD (adjusted for age, sex, baseline disease burden, without MPN adjustment)\n")
cvd_total_cox <- coxph(Surv(time_to_cvd, cvd_event) ~ chr9_neutral + age + sex + n_baseline_diseases, 
                       data = analysis_data)
print(summary(cvd_total_cox))

# Direct Effect
cat("\nStep 2: chr9_neutral + MPN → CVD (Cox model, adjusted for age, sex, baseline disease burden)\n")
cvd_cox <- coxph(Surv(time_to_cvd, cvd_event) ~ chr9_neutral + mpn_status + age + sex + n_baseline_diseases, 
                 data = analysis_data)
print(summary(cvd_cox))

# Extract HRs - with error handling for NA coefficients
mpn_hr <- tryCatch(exp(coef(mpn_cox)["chr9_neutral"]), error = function(e) NA)
mpn_ci <- tryCatch(exp(confint(mpn_cox)["chr9_neutral", ]), error = function(e) c(NA, NA))
mpn_p <- tryCatch(summary(mpn_cox)$coefficients["chr9_neutral", "Pr(>|z|)"], error = function(e) NA)

cvd_direct_hr <- tryCatch(exp(coef(cvd_cox)["chr9_neutral"]), error = function(e) NA)
cvd_direct_ci <- tryCatch(exp(confint(cvd_cox)["chr9_neutral", ]), error = function(e) c(NA, NA))
cvd_direct_p <- tryCatch(summary(cvd_cox)$coefficients["chr9_neutral", "Pr(>|z|)"], error = function(e) NA)

mpn_to_cvd_hr <- tryCatch(exp(coef(cvd_cox)["mpn_status"]), error = function(e) NA)
mpn_to_cvd_ci <- tryCatch(exp(confint(cvd_cox)["mpn_status", ]), error = function(e) c(NA, NA))
mpn_to_cvd_p <- tryCatch(summary(cvd_cox)$coefficients["mpn_status", "Pr(>|z|)"], error = function(e) NA)

total_hr <- tryCatch(exp(coef(cvd_total_cox)["chr9_neutral"]), error = function(e) NA)
total_ci <- tryCatch(exp(confint(cvd_total_cox)["chr9_neutral", ]), error = function(e) c(NA, NA))
total_p <- tryCatch(summary(cvd_total_cox)$coefficients["chr9_neutral", "Pr(>|z|)"], error = function(e) NA)

# ==========================================
# LOGISTIC REGRESSION DIAGNOSTICS
# ==========================================

cat("\n=== 3. LOGISTIC REGRESSION DIAGNOSTICS ===\n\n")

cat("Mediator Model (chr9_neutral → MPN):\n")
mediator_model <- glm(mpn_status ~ chr9_neutral + age + sex + n_baseline_diseases, 
                      data = analysis_data, 
                      family = binomial())
print(summary(mediator_model))

cat("\nIterations:", mediator_model$iter, "\n")
if (mediator_model$iter > 10) {
  cat("❌ WARNING: >10 iterations. Possible convergence issues!\n")
} else {
  cat("✓ OK: Model converged normally\n")
}

se_med <- summary(mediator_model)$coefficients[, "Std. Error"]
if (any(se_med > 2)) {
  cat("❌ WARNING: Large SE detected (max =", round(max(se_med), 2), "). Quasi-separation!\n")
} else {
  cat("✓ OK: Standard errors reasonable\n")
}

cat("\nOutcome Model (chr9_neutral + MPN → CVD):\n")
outcome_model <- glm(cvd_status ~ chr9_neutral + mpn_status + age + sex + n_baseline_diseases, 
                     data = analysis_data, 
                     family = binomial())
print(summary(outcome_model))

cat("\nIterations:", outcome_model$iter, "\n")
if (outcome_model$iter > 10) {
  cat("❌ WARNING: >10 iterations. Possible convergence issues!\n")
} else {
  cat("✓ OK: Model converged normally\n")
}

se_out <- summary(outcome_model)$coefficients[, "Std. Error"]
if (any(se_out > 2)) {
  cat("❌ WARNING: Large SE detected (max =", round(max(se_out), 2), "). Quasi-separation!\n")
} else {
  cat("✓ OK: Standard errors reasonable\n")
}

# ==========================================
# MEDIATION ANALYSIS (BINARY)
# ==========================================

cat("\n=== 4. BINARY MEDIATION ANALYSIS ===\n\n")

# Check if mediation is feasible
mediation_feasible <- TRUE

# Check for NA coefficients in mediator model
if (any(is.na(coef(mediator_model)))) {
  cat("❌ WARNING: NA coefficients in mediator model. Mediation analysis may fail.\n")
  mediation_feasible <- FALSE
}

# Check for NA coefficients in outcome model
if (any(is.na(coef(outcome_model)))) {
  cat("❌ WARNING: NA coefficients in outcome model. Mediation analysis may fail.\n")
  mediation_feasible <- FALSE
}

if (mediation_feasible) {
  set.seed(123)
  mediation_results <- tryCatch({
    mediate(
      model.m = mediator_model,
      model.y = outcome_model,
      treat = "chr9_neutral",
      mediator = "mpn_status",
      robustSE = TRUE,
      sims = 1000
    )
  }, error = function(e) {
    cat("❌ Mediation analysis failed with error:", conditionMessage(e), "\n")
    NULL
  })
  
  if (!is.null(mediation_results)) {
    summary(mediation_results)
  }
} else {
  mediation_results <- NULL
  cat("Skipping mediation analysis due to model issues.\n")
}

# ==========================================
# SUMMARY RESULTS TABLE
# ==========================================

cat("\n")
cat(rep("=", 80), "\n", sep = "")
cat("TIME-TO-EVENT MEDIATION RESULTS SUMMARY\n")
cat(rep("=", 80), "\n\n", sep = "")

survival_results <- data.frame(
  Effect = c("Total Effect (9= → CVD)",
             "Direct Effect (9= → CVD, adjusted for MPN)", 
             "Pathway A (9= → MPN)",
             "Pathway B (MPN → CVD)"),
  Hazard_Ratio = c(round(total_hr, 3),
                   round(cvd_direct_hr, 3),
                   round(mpn_hr, 3),
                   round(mpn_to_cvd_hr, 3)),
  CI_Lower = c(round(total_ci[1], 3),
               round(cvd_direct_ci[1], 3), 
               round(mpn_ci[1], 3),
               round(mpn_to_cvd_ci[1], 3)),
  CI_Upper = c(round(total_ci[2], 3),
               round(cvd_direct_ci[2], 3),
               round(mpn_ci[2], 3), 
               round(mpn_to_cvd_ci[2], 3)),
  P_Value = c(ifelse(is.na(total_p), "NA", ifelse(total_p < 0.001, "< 0.001", round(total_p, 3))),
              ifelse(is.na(cvd_direct_p), "NA", ifelse(cvd_direct_p < 0.001, "< 0.001", round(cvd_direct_p, 3))),
              ifelse(is.na(mpn_p), "NA", ifelse(mpn_p < 0.001, "< 0.001", round(mpn_p, 3))),
              ifelse(is.na(mpn_to_cvd_p), "NA", ifelse(mpn_to_cvd_p < 0.001, "< 0.001", round(mpn_to_cvd_p, 3)))),
  Significance = c(ifelse(is.na(total_p), "NA", ifelse(total_p < 0.05, "Significant", "Not Significant")),
                   ifelse(is.na(cvd_direct_p), "NA", ifelse(cvd_direct_p < 0.05, "Significant", "Not Significant")),
                   ifelse(is.na(mpn_p), "NA", ifelse(mpn_p < 0.05, "Significant", "Not Significant")), 
                   ifelse(is.na(mpn_to_cvd_p), "NA", ifelse(mpn_to_cvd_p < 0.05, "Significant", "Not Significant")))
)

print(survival_results, row.names = FALSE)

# Save Cox model results
write.csv(survival_results, "~/Desktop/UK BB/chr9_neutral_mpn_cvd_survival_mediation_results.csv", row.names = FALSE)

# ==========================================
# CALCULATE BASELINE RISKS AND CREATE RISK TABLE
# ==========================================

# Calculate baseline CVD risk in No_mCA group
baseline_cvd_risk <- analysis_data %>%
  filter(No_mCA == 1) %>%
  summarise(baseline_risk = mean(cvd_status) * 100) %>%
  pull(baseline_risk)

# Calculate actual chr9_neutral CVD risk
chr9_neutral_cvd_risk <- analysis_data %>%
  filter(chr9_neutral == 1) %>%
  summarise(chr9_neutral_risk = mean(cvd_status) * 100) %>%
  pull(chr9_neutral_risk)

# Extract mediation effects if available
if (!is.null(mediation_results)) {
  total_effect_pct <- mediation_results$tau.coef * 100
  acme_effect_pct <- mediation_results$d0 * 100
  ade_effect_pct <- mediation_results$z0 * 100
  
  risk_table <- data.frame(
    Group = c("No_mCA (baseline)", 
              "chr9_neutral (observed)", 
              "chr9_neutral (through MPN pathway)",
              "chr9_neutral (direct pathway)"),
    CVD_Risk_Percent = c(
      round(baseline_cvd_risk, 2),
      round(chr9_neutral_cvd_risk, 2),
      round(baseline_cvd_risk + acme_effect_pct, 2),
      round(baseline_cvd_risk + ade_effect_pct, 2)
    ),
    Risk_Difference = c(
      "Reference",
      paste0(round(chr9_neutral_cvd_risk - baseline_cvd_risk, 2), "%"),
      paste0(round(acme_effect_pct, 2), "% (through MPN)"),
      paste0(round(ade_effect_pct, 2), "% (direct)")
    ),
    P_Value = c(
      "-",
      "-",
      ifelse(mediation_results$d0.p < 0.001, "< 0.001", round(mediation_results$d0.p, 3)),
      ifelse(mediation_results$z0.p < 0.001, "< 0.001", round(mediation_results$z0.p, 3))
    )
  )
} else {
  risk_table <- data.frame(
    Group = c("No_mCA (baseline)", "chr9_neutral (observed)"),
    CVD_Risk_Percent = c(round(baseline_cvd_risk, 2), round(chr9_neutral_cvd_risk, 2)),
    Risk_Difference = c("Reference", paste0(round(chr9_neutral_cvd_risk - baseline_cvd_risk, 2), "%")),
    P_Value = c("-", "-")
  )
  acme_effect_pct <- NA
  ade_effect_pct <- NA
}

cat("\n=== CVD RISK TABLE (chr9 neutral → MPN → CVD) ===\n")
print(risk_table, row.names = FALSE)

cat("\n=== SUMMARY ===\n")
cat("Baseline CVD risk (No_mCA):", round(baseline_cvd_risk, 2), "%\n")
cat("Observed chr9_neutral CVD risk:", round(chr9_neutral_cvd_risk, 2), "%\n")
cat("Total risk difference:", round(chr9_neutral_cvd_risk - baseline_cvd_risk, 2), "percentage points\n")
if (!is.na(acme_effect_pct)) {
  cat("Mediated through MPN:", round(acme_effect_pct, 2), "percentage points\n")
  cat("Direct effect:", round(ade_effect_pct, 2), "percentage points\n")
  if (!is.null(mediation_results) && !is.null(mediation_results$n0)) {
    cat("Proportion mediated:", round(mediation_results$n0 * 100, 1), "% of total effect\n")
  }
}

# Save risk table
write.csv(risk_table, "~/Desktop/UK BB/chr9_neutral_mpn_cvd_risk_table.csv", row.names = FALSE)

# ==========================================
# FIGURE 1: MEDIATION PLOT (Binary mediation results)
# ==========================================

cat("\n=== GENERATING FIGURES ===\n\n")

if (!is.null(mediation_results)) {
  pdf("~/Desktop/UK BB/chr9_neutral_mpn_cvd_mediation_plot.pdf", width = 8, height = 8)
  plot(mediation_results, 
       main = "Mediation of chr9 Neutral Association with CVD\nthrough MPN Development",
       xlab = "Absolute Risk Increase (percentage points)")
  dev.off()
  cat("Figure 1 saved: chr9_neutral_mpn_cvd_mediation_plot.pdf\n")
} else {
  cat("Figure 1 skipped: Mediation analysis did not complete\n")
}

# ==========================================
# FIGURE 2: KAPLAN-MEIER CURVE - chr9_neutral → MPN (Pathway A)
# ==========================================

surv_mpn <- survfit(Surv(time_to_mpn, mpn_event) ~ chr9_neutral_factor, data = analysis_data)

km_mpn <- ggsurvplot(
  surv_mpn,
  data = analysis_data,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  palette = c("#1f78b4", "#e31a1c"),
  title = "Pathway A: chr9 Neutral to MPN Development",
  xlab = "Years from mCA assessment",
  ylab = "MPN-Free Survival Probability",
  legend.title = "Group",
  legend.labs = c("No mCA", "9="),
  break.time.by = 2,
  xlim = c(0, 16),
  ggtheme = theme_minimal(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.3)
    ),
  risk.table.height = 0.25
)

pdf("~/Desktop/UK BB/chr9_neutral_mpn_cvd_KM_pathway_A.pdf", width = 7, height = 9)
print(km_mpn)
dev.off()
cat("Figure 2 saved: chr9_neutral_mpn_cvd_KM_pathway_A.pdf\n")

# ==========================================
# FIGURE 3: KAPLAN-MEIER CURVE - chr9_neutral → CVD (Total Effect)
# ==========================================

surv_cvd_total <- survfit(Surv(time_to_cvd, cvd_event) ~ chr9_neutral_factor, data = analysis_data)

km_cvd_total <- ggsurvplot(
  surv_cvd_total,
  data = analysis_data,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  palette = c("#1f78b4", "#e31a1c"),
  title = "Total Effect: chr9 Neutral to Cardiovascular Disease",
  xlab = "Years from mCA assessment",
  ylab = "CVD-Free Survival Probability",
  legend.title = "Group",
  legend.labs = c("No mCA", "9="),
  break.time.by = 2,
  xlim = c(0, 16),
  ggtheme = theme_minimal(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.3)
    ),
  risk.table.height = 0.25
)

pdf("~/Desktop/UK BB/chr9_neutral_mpn_cvd_KM_total_effect.pdf", width = 7, height = 9)
print(km_cvd_total)
dev.off()
cat("Figure 3 saved: chr9_neutral_mpn_cvd_KM_total_effect.pdf\n")

# ==========================================
# FIGURE 4: KAPLAN-MEIER CURVE - MPN → CVD (Pathway B)
# ==========================================

surv_cvd_mpn <- survfit(Surv(time_to_cvd, cvd_event) ~ mpn_factor, data = analysis_data)

km_cvd_mpn <- ggsurvplot(
  surv_cvd_mpn,
  data = analysis_data,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  palette = c("#1f78b4", "#e31a1c"),
  title = "Pathway B: MPN to Cardiovascular Disease",
  xlab = "Years from mCA assessment",
  ylab = "CVD-Free Survival Probability",
  legend.title = "MPN Status",
  legend.labs = c("No MPN", "MPN"),
  break.time.by = 2,
  xlim = c(0, 16),
  ggtheme = theme_minimal(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.3)
    ),
  risk.table.height = 0.25
)

pdf("~/Desktop/UK BB/chr9_neutral_mpn_cvd_KM_pathway_B.pdf", width = 7, height = 9)
print(km_cvd_mpn)
dev.off()
cat("Figure 4 saved: chr9_neutral_mpn_cvd_KM_pathway_B.pdf\n")

# ==========================================
# FIGURE 5: FOREST PLOT - All Pathway HRs
# ==========================================

forest_data <- data.frame(
  Pathway = c("Total Effect\n(9= -> CVD)", 
              "Direct Effect\n(9= -> CVD | MPN)", 
              "Pathway A\n(9= -> MPN)", 
              "Pathway B\n(MPN -> CVD)"),
  HR = c(total_hr, cvd_direct_hr, mpn_hr, mpn_to_cvd_hr),
  CI_Lower = c(total_ci[1], cvd_direct_ci[1], mpn_ci[1], mpn_to_cvd_ci[1]),
  CI_Upper = c(total_ci[2], cvd_direct_ci[2], mpn_ci[2], mpn_to_cvd_ci[2]),
  P_Value = c(total_p, cvd_direct_p, mpn_p, mpn_to_cvd_p)
)

# Remove NA rows for plotting
forest_data_plot <- forest_data %>% filter(!is.na(HR))

if (nrow(forest_data_plot) > 0) {
  forest_data_plot$Significant <- ifelse(forest_data_plot$P_Value < 0.05, "Yes", "No")
  
  forest_data_plot$Pathway <- factor(forest_data_plot$Pathway, 
                                     levels = rev(forest_data_plot$Pathway))
  
  forest_plot <- ggplot(forest_data_plot, aes(x = HR, y = Pathway)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper, color = Significant), 
                   height = 0.2, linewidth = 1) +
    geom_point(aes(color = Significant), size = 4) +
    scale_color_manual(values = c("Yes" = "#E41A1C", "No" = "#377EB8"),
                       name = "P < 0.05") +
    scale_x_log10(breaks = c(0.5, 1, 2, 5, 10, 20, 50, 100),
                  labels = c("0.5", "1", "2", "5", "10", "20", "50", "100")) +
    labs(
      title = "Mediation Analysis: chr9 Neutral to MPN to CVD",
      subtitle = "Hazard Ratios with 95% Confidence Intervals (log scale)",
      x = "Hazard Ratio (log scale)",
      y = ""
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      axis.text.x = element_text(size = 10),
      legend.position = "bottom"
    ) +
    geom_text(aes(label = sprintf("HR=%.2f", HR)), 
              hjust = -0.3, vjust = -0.5, size = 3.5)
  
  ggsave("~/Desktop/UK BB/chr9_neutral_mpn_cvd_forest_plot.pdf", 
         forest_plot, width = 7, height = 7)
  cat("Figure 5 saved: chr9_neutral_mpn_cvd_forest_plot.pdf\n")
} else {
  cat("Figure 5 skipped: No valid HR estimates\n")
}

# ==========================================
# FIGURE 6: PATHWAY DIAGRAM (Simplified visual)
# ==========================================

pathway_diagram <- ggplot() +
  annotate("rect", xmin = 0, xmax = 2, ymin = 0.8, ymax = 1.2, 
           fill = "#FFF2CC", color = "#D6B656", linewidth = 1) +
  annotate("rect", xmin = 3, xmax = 5, ymin = 0.8, ymax = 1.2, 
           fill = "#DAE8FC", color = "#6C8EBF", linewidth = 1) +
  annotate("rect", xmin = 6, xmax = 8, ymin = 0.8, ymax = 1.2, 
           fill = "#D5E8D4", color = "#82B366", linewidth = 1) +
  
  annotate("text", x = 1, y = 1, label = "chr9 Neutral\n(9=)", size = 4.5, fontface = "bold") +
  annotate("text", x = 4, y = 1, label = "MPN", size = 4.5, fontface = "bold") +
  annotate("text", x = 7, y = 1, label = "Cardiovascular\nDisease", size = 4, fontface = "bold") +
  
  annotate("segment", x = 2, xend = 2.9, y = 1, yend = 1, 
           arrow = arrow(length = unit(0.3, "cm")), linewidth = 1.2, color = "#E41A1C") +
  annotate("text", x = 2.45, y = 1.15, 
           label = ifelse(is.na(mpn_hr), "Pathway A\nHR=NA", 
                          sprintf("Pathway A\nHR=%.2f%s", mpn_hr, ifelse(is.na(mpn_p) | mpn_p >= 0.05, " (NS)", "***"))), 
           size = 3.5, color = "#E41A1C", fontface = "bold") +
  
  annotate("segment", x = 5, xend = 5.9, y = 1, yend = 1, 
           arrow = arrow(length = unit(0.3, "cm")), linewidth = 1.2, color = "#E41A1C") +
  annotate("text", x = 5.45, y = 1.15, 
           label = ifelse(is.na(mpn_to_cvd_hr), "Pathway B\nHR=NA",
                          sprintf("Pathway B\nHR=%.2f%s", mpn_to_cvd_hr, ifelse(is.na(mpn_to_cvd_p) | mpn_to_cvd_p >= 0.05, " (NS)", "***"))), 
           size = 3.5, color = "#E41A1C", fontface = "bold") +
  
  annotate("curve", x = 2, xend = 5.9, y = 0.75, yend = 0.75,
           curvature = -0.3, arrow = arrow(length = unit(0.3, "cm")), 
           linewidth = 1, color = "#377EB8", linetype = "dashed") +
  annotate("text", x = 4, y = 0.45, 
           label = ifelse(is.na(cvd_direct_hr), "Direct Effect\nHR=NA",
                          sprintf("Direct Effect\nHR=%.2f%s", cvd_direct_hr, 
                                  ifelse(is.na(cvd_direct_p) | cvd_direct_p >= 0.05, " (NS)", "***"))), 
           size = 3.5, color = "#377EB8") +
  
  annotate("text", x = 4, y = 0.15, 
           label = ifelse(is.na(total_hr), "Total Effect: HR=NA",
                          sprintf("Total Effect: HR=%.2f%s\n*** p<0.05, NS = Not Significant", 
                                  total_hr, ifelse(is.na(total_p) | total_p >= 0.05, " (NS)", "***"))), 
           size = 3.5, color = "gray30") +
  
  labs(title = "Mediation Pathway: chr9 Neutral to MPN to CVD",
       subtitle = ifelse(is.na(total_p) | total_p >= 0.05, "Null Total Effect", "Significant Mediation")) +
  
  theme_void() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  coord_cartesian(xlim = c(-0.5, 8.5), ylim = c(0, 1.5))

ggsave("~/Desktop/UK BB/chr9_neutral_mpn_cvd_pathway_diagram.pdf", 
       pathway_diagram, width = 10, height = 5)
cat("Figure 6 saved: chr9_neutral_mpn_cvd_pathway_diagram.pdf\n")

# ==========================================
# MEDIATION INTERPRETATION
# ==========================================

cat("\n")
cat(rep("=", 80), "\n", sep = "")
cat("MEDIATION INTERPRETATION\n")
cat(rep("=", 80), "\n\n", sep = "")

if(is.na(mpn_p) || is.na(mpn_to_cvd_p)) {
  cat("RESULT: Unable to fully assess mediation due to model estimation issues\n")
  cat("  -> This is likely due to very small sample sizes or rare events\n")
} else if(mpn_p < 0.05 && mpn_to_cvd_p < 0.05) {
  if(is.na(total_p) || total_p >= 0.05) {
    cat("RESULT: NULL TOTAL EFFECT despite significant individual pathways\n")
    cat("  -> chr9_neutral predisposes to MPN (Pathway A significant)\n")
    cat("  -> MPN significantly increases CVD risk (Pathway B significant)\n")
    cat("  -> BUT: chr9_neutral does NOT significantly increase overall CVD risk\n")
    cat("  -> INTERPRETATION: Low MPN incidence among chr9_neutral carriers limits\n")
    cat("                    the population-level impact of the mediated pathway\n")
  } else if(is.na(cvd_direct_p) || cvd_direct_p >= 0.05) {
    cat("RESULT: Full mediation - chr9_neutral -> CVD association is fully mediated through MPN\n")
  } else {
    cat("RESULT: Partial mediation - chr9_neutral -> CVD has both direct and MPN-mediated pathways\n")
  }
} else if(!is.na(mpn_p) && mpn_p < 0.05 && (is.na(mpn_to_cvd_p) || mpn_to_cvd_p >= 0.05)) {
  cat("RESULT: No mediation - MPN does not significantly affect CVD risk\n")
} else if(is.na(mpn_p) || mpn_p >= 0.05) {
  cat("RESULT: No mediation - chr9_neutral does not significantly affect MPN risk\n") 
} else {
  cat("RESULT: Complex pattern - see individual pathway significance\n")
}

cat("\nchr9_neutral -> MPN pathway:", 
    ifelse(is.na(mpn_p), "NA", ifelse(mpn_p < 0.05, "SIGNIFICANT", "NOT SIGNIFICANT")), "\n")
cat("MPN -> CVD pathway:", 
    ifelse(is.na(mpn_to_cvd_p), "NA", ifelse(mpn_to_cvd_p < 0.05, "SIGNIFICANT", "NOT SIGNIFICANT")), "\n")
cat("chr9_neutral -> CVD total effect:", 
    ifelse(is.na(total_p), "NA", ifelse(total_p < 0.05, "SIGNIFICANT", "NOT SIGNIFICANT")), "\n")
cat("chr9_neutral -> CVD direct effect:", 
    ifelse(is.na(cvd_direct_p), "NA", ifelse(cvd_direct_p < 0.05, "SIGNIFICANT", "NOT SIGNIFICANT")), "\n")

# ==========================================
# FINAL ASSESSMENT
# ==========================================

cat("\n")
cat(rep("=", 80), "\n", sep = "")
cat("OVERALL ASSESSMENT\n")
cat(rep("=", 80), "\n\n", sep = "")

warnings_count <- 0

if (any(tbl1 == 0) || any(tbl2 == 0) || any(tbl3 == 0)) warnings_count <- warnings_count + 1
if (mediator_model$iter > 10 || outcome_model$iter > 10) warnings_count <- warnings_count + 1
if (any(se_med > 2) || any(se_out > 2)) warnings_count <- warnings_count + 1
if (!is.na(mpn_hr) && mpn_hr > 50) warnings_count <- warnings_count + 0.5
if (!is.na(total_hr) && total_hr > 50) warnings_count <- warnings_count + 0.5
if (sum(analysis_data$chr9_neutral) < 100) warnings_count <- warnings_count + 1

if (warnings_count == 0) {
  cat("✓ EXCELLENT: No issues detected. Models fit well.\n")
  cat("  -> Results can be reported with confidence\n")
} else if (warnings_count < 2) {
  cat("⚠️  CAUTION: Minor issues detected.\n")
  cat("  -> Results valid but note:\n")
  cat("     • Present absolute risks alongside HRs\n")
  cat("     • Acknowledge wide CIs if present\n")
} else {
  cat("❌ CONCERN: Multiple issues detected.\n")
  cat("  -> RECOMMENDATIONS:\n")
  cat("     1. chr9_neutral sample size (n=", sum(analysis_data$chr9_neutral), ")\n")
  cat("     2. Consider this analysis EXPLORATORY if sample size is small\n")
  cat("     3. Results may be underpowered\n")
  cat("     4. Consider combining with other chromosomal abnormalities\n")
  cat("     5. MPN events in exposed:", sum(analysis_data$mpn_event[analysis_data$chr9_neutral == 1]), "\n")
}

# ==========================================
# FILES SAVED SUMMARY
# ==========================================

cat("\n")
cat(rep("=", 80), "\n", sep = "")
cat("FILES SAVED\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("=== CSV FILES ===\n")
cat("1. chr9_neutral_mpn_cvd_survival_mediation_results.csv - Cox model results\n")
cat("2. chr9_neutral_mpn_cvd_risk_table.csv - Risk comparison table\n")

cat("\n=== PDF FIGURES ===\n")
cat("1. chr9_neutral_mpn_cvd_mediation_plot.pdf - Binary mediation analysis plot\n")
cat("2. chr9_neutral_mpn_cvd_KM_pathway_A.pdf - Kaplan-Meier: chr9 neutral -> MPN\n")
cat("3. chr9_neutral_mpn_cvd_KM_total_effect.pdf - Kaplan-Meier: chr9 neutral -> CVD\n")
cat("4. chr9_neutral_mpn_cvd_KM_pathway_B.pdf - Kaplan-Meier: MPN -> CVD\n")
cat("5. chr9_neutral_mpn_cvd_forest_plot.pdf - Forest plot of all pathway HRs\n")
cat("6. chr9_neutral_mpn_cvd_pathway_diagram.pdf - Visual pathway diagram\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
