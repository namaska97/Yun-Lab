# ==========================================
# COMPLETE ANALYSIS: CHR13 LOSS → CLL → ACTINIC KERATOSIS
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

cat("Loading UKBB data for chr13 loss → CLL → Actinic Keratosis analysis...\n")

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
# Identify chr13 loss carriers
# ===========================
chr13_loss_ids <- cnv_all %>%
  filter(CHR == 13, COPY_CHANGE == "loss") %>%
  distinct(ID) %>%
  pull(ID)

mca_info <- mca_info %>%
  mutate(
    chr13_loss = ifelse(ID %in% chr13_loss_ids, 1, 0),
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
# DIAGNOSTIC: Check CLL codes in database
# ===========================
cat("\n=== DIAGNOSTIC: Checking CLL codes in database ===\n")
cll_codes <- diag_data %>%
  filter(grepl("^C91", diagnosis)) %>%
  group_by(diagnosis) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n))

cat("C91 codes found in database:\n")
print(cll_codes)
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
cat("Mean baseline diseases - chr13_loss:", 
    round(mean(mca_info$n_baseline_diseases[mca_info$chr13_loss == 1], na.rm = TRUE), 2), "\n")
cat("Mean baseline diseases - No_mCA:", 
    round(mean(mca_info$n_baseline_diseases[mca_info$No_mCA == 1], na.rm = TRUE), 2), "\n\n")

# Compare baseline disease burden between groups
baseline_comparison <- mca_info %>%
  filter(chr13_loss == 1 | No_mCA == 1) %>%
  mutate(group = ifelse(chr13_loss == 1, "chr13_loss", "No_mCA")) %>%
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
# CLL - INCIDENT CASES
# ===========================
cat("\n=== Identifying INCIDENT CLL cases ===\n")

# CLL ICD-10 codes
# C911 - Chronic lymphocytic leukemia of B-cell type
cll_pattern <- "^C911$|^C9110$|^C9111$|^C9112$"

cat("\n=== CLL ICD-10 Codes ===\n")
cat("C911 - Chronic lymphocytic leukemia of B-cell type\n")
cat("C9110 - CLL not having achieved remission\n")
cat("C9111 - CLL in remission\n")
cat("C9112 - CLL in relapse\n\n")

# Show what CLL codes will be matched
cll_codes_matched <- diag_data %>%
  filter(grepl(cll_pattern, diagnosis)) %>%
  group_by(diagnosis) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n))

cat("CLL codes that will be captured by pattern:\n")
print(cll_codes_matched)
cat("\n")

pre_existing_cll <- diag_data %>%
  filter(grepl(cll_pattern, diagnosis)) %>%
  inner_join(mca_info %>% select(ID, mCA_date), 
             by = c("patient_id" = "ID"),
             relationship = "many-to-one") %>%
  filter(admission_date <= mCA_date) %>%
  distinct(patient_id, diagnosis)

cat("Patients with pre-existing CLL codes:", n_distinct(pre_existing_cll$patient_id), "\n")

post_mca_cll <- diag_data %>%
  filter(grepl(cll_pattern, diagnosis)) %>%
  inner_join(mca_info %>% select(ID, mCA_date), 
             by = c("patient_id" = "ID"),
             relationship = "many-to-one") %>%
  filter(admission_date > mCA_date)

cat("Post-mCA CLL diagnoses:", nrow(post_mca_cll), "\n")

incident_cll <- post_mca_cll %>%
  anti_join(pre_existing_cll, by = c("patient_id", "diagnosis")) %>%
  group_by(patient_id) %>%
  summarise(first_cll_date = min(admission_date), .groups = "drop")

cat("INCIDENT CLL cases (new codes only):", nrow(incident_cll), "patients\n")

# ===========================
# ACTINIC KERATOSIS - INCIDENT CASES
# ===========================
cat("\n=== Identifying INCIDENT Actinic Keratosis cases ===\n")

# Actinic Keratosis ICD-10 codes
# L57 - Skin changes due to chronic exposure to nonionizing radiation
# L570 - Actinic keratosis
ak_pattern <- "^L57$|^L570$"

cat("\n=== Actinic Keratosis ICD-10 Codes ===\n")
cat("L57 - Skin changes due to chronic exposure to nonionizing radiation\n")
cat("L570 - Actinic keratosis\n\n")

pre_existing_ak <- diag_data %>%
  filter(grepl(ak_pattern, diagnosis)) %>%
  inner_join(mca_info %>% select(ID, mCA_date), 
             by = c("patient_id" = "ID"),
             relationship = "many-to-one") %>%
  filter(admission_date <= mCA_date) %>%
  distinct(patient_id, diagnosis)

cat("Patients with pre-existing Actinic Keratosis codes:", n_distinct(pre_existing_ak$patient_id), "\n")

post_mca_ak <- diag_data %>%
  filter(grepl(ak_pattern, diagnosis)) %>%
  inner_join(mca_info %>% select(ID, mCA_date), 
             by = c("patient_id" = "ID"),
             relationship = "many-to-one") %>%
  filter(admission_date > mCA_date)

cat("Post-mCA Actinic Keratosis diagnoses:", nrow(post_mca_ak), "\n")

incident_ak <- post_mca_ak %>%
  anti_join(pre_existing_ak, by = c("patient_id", "diagnosis")) %>%
  group_by(patient_id) %>%
  summarise(first_ak_date = min(admission_date), .groups = "drop")

cat("INCIDENT Actinic Keratosis cases (new codes only):", nrow(incident_ak), "patients\n")

# ===========================
# Combine all data
# ===========================
cohort <- mca_info %>%
  left_join(incident_cll, by = c("ID" = "patient_id")) %>%
  left_join(incident_ak, by = c("ID" = "patient_id")) %>%
  mutate(
    cll_status = ifelse(!is.na(first_cll_date), 1, 0),
    ak_status = ifelse(!is.na(first_ak_date), 1, 0)
  )

cat("\n=== Cohort Summary ===\n")
cat("Total cohort:", nrow(cohort), "participants\n")
cat("chr13_loss cases:", sum(cohort$chr13_loss), "\n")
cat("No_mCA cases:", sum(cohort$No_mCA), "\n")
cat("INCIDENT CLL cases:", sum(cohort$cll_status), "\n")
cat("INCIDENT Actinic Keratosis cases:", sum(cohort$ak_status), "\n\n")

# ===========================
# Check if we have enough CLL cases to proceed
# ===========================
if (sum(cohort$cll_status) == 0) {
  cat("\n")
  cat(rep("!", 80), "\n", sep = "")
  cat("CRITICAL ERROR: No CLL cases found!\n")
  cat("Cannot proceed with mediation analysis.\n")
  cat("Please check ICD-10 coding in your database.\n")
  cat(rep("!", 80), "\n", sep = "")
  stop("No CLL cases found - analysis cannot proceed")
}

# ===========================
# Create analysis dataset
# ===========================
analysis_data <- cohort %>%
  filter(chr13_loss == 1 | No_mCA == 1) %>%
  filter(complete.cases(chr13_loss, cll_status, ak_status, age, sex, n_baseline_diseases)) %>%
  mutate(
    time_to_cll = ifelse(cll_status == 1, 
                         as.numeric(first_cll_date - mCA_date)/365.25,
                         as.numeric(censor_date - mCA_date)/365.25),
    cll_event = cll_status,
    
    time_to_ak = ifelse(ak_status == 1,
                        as.numeric(first_ak_date - mCA_date)/365.25,
                        as.numeric(censor_date - mCA_date)/365.25),
    ak_event = ak_status,
    
    time_to_cll = pmax(time_to_cll, 0.001),
    time_to_ak = pmax(time_to_ak, 0.001),
    
    # Create factor for plotting
    chr13_loss_factor = factor(chr13_loss, levels = c(0, 1), labels = c("No mCA", "13-")),
    cll_factor = factor(cll_status, levels = c(0, 1), labels = c("No CLL", "CLL"))
  )

cat("=== Analysis Dataset ===\n")
cat("Total participants:", nrow(analysis_data), "\n")
cat("chr13_loss carriers:", sum(analysis_data$chr13_loss), "\n")
cat("No_mCA controls:", sum(analysis_data$No_mCA), "\n")
cat("INCIDENT CLL events:", sum(analysis_data$cll_event), "\n")
cat("INCIDENT Actinic Keratosis events:", sum(analysis_data$ak_event), "\n")
cat("Median follow-up time:", round(median(analysis_data$time_to_ak), 2), "years\n\n")

# ===========================
# POWER WARNING CHECK
# ===========================
cat("\n=== POWER CHECK ===\n")
if (sum(analysis_data$chr13_loss) < 100) {
  cat("⚠️  WARNING: Small exposure group (n=", sum(analysis_data$chr13_loss), ")\n")
  cat("   Results may be underpowered and unstable.\n")
  cat("   Consider combining with other chromosomal abnormalities or\n")
  cat("   interpreting results as exploratory.\n\n")
}

cll_in_exposed <- sum(analysis_data$cll_event[analysis_data$chr13_loss == 1])
cat("CLL cases in chr13_loss group:", cll_in_exposed, "\n")
if (cll_in_exposed < 5) {
  cat("⚠️  WARNING: Very few CLL cases in exposed group.\n")
  cat("   Pathway A estimates may be unstable.\n\n")
}

# ==========================================
# BUILT-IN DIAGNOSTICS BEFORE MEDIATION
# ==========================================

cat("\n")
cat(rep("=", 80), "\n", sep = "")
cat("MODEL DIAGNOSTICS - CHR13 LOSS → CLL → ACTINIC KERATOSIS\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("=== 1. CONTINGENCY TABLES (Checking for Separation) ===\n\n")

cat("A. chr13_loss vs CLL events:\n")
tbl1 <- table(analysis_data$chr13_loss, analysis_data$cll_event,
              dnn = c("chr13_loss", "CLL_event"))
print(tbl1)

if (any(tbl1 == 0)) {
  cat("❌ WARNING: ZERO cells detected! Complete separation!\n")
} else if (any(tbl1 > 0 & tbl1 < 5)) {
  cat("⚠️  WARNING: Small cells (<5) detected. Quasi-separation possible.\n")
} else {
  cat("✓ OK: No separation issues\n")
}

cat("\nB. chr13_loss vs Actinic Keratosis events:\n")
tbl2 <- table(analysis_data$chr13_loss, analysis_data$ak_event,
              dnn = c("chr13_loss", "AK_event"))
print(tbl2)

if (any(tbl2 == 0)) {
  cat("❌ WARNING: ZERO cells detected! Complete separation!\n")
} else if (any(tbl2 > 0 & tbl2 < 5)) {
  cat("⚠️  WARNING: Small cells (<5) detected. Quasi-separation possible.\n")
} else {
  cat("✓ OK: No separation issues\n")
}

cat("\nC. CLL vs Actinic Keratosis (among chr13_loss=1 only) - PATH B CHECK:\n")
exposed <- analysis_data[analysis_data$chr13_loss == 1, ]
tbl3 <- table(exposed$cll_event, exposed$ak_event,
              dnn = c("CLL_event", "AK_event"))
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

# Pathway A: chr13_loss -> CLL
cat("Step 1: chr13_loss → CLL (Cox model, adjusted for age, sex, baseline disease burden)\n")
cll_cox <- coxph(Surv(time_to_cll, cll_event) ~ chr13_loss + age + sex + n_baseline_diseases, 
                 data = analysis_data)
print(summary(cll_cox))

# Total Effect
cat("\nTotal Effect: chr13_loss → Actinic Keratosis (adjusted for age, sex, baseline disease burden, without CLL adjustment)\n")
ak_total_cox <- coxph(Surv(time_to_ak, ak_event) ~ chr13_loss + age + sex + n_baseline_diseases, 
                      data = analysis_data)
print(summary(ak_total_cox))

# Direct Effect
cat("\nStep 2: chr13_loss + CLL → Actinic Keratosis (Cox model, adjusted for age, sex, baseline disease burden)\n")
ak_cox <- coxph(Surv(time_to_ak, ak_event) ~ chr13_loss + cll_status + age + sex + n_baseline_diseases, 
                data = analysis_data)
print(summary(ak_cox))

# Extract HRs - with error handling for NA coefficients
cll_hr <- tryCatch(exp(coef(cll_cox)["chr13_loss"]), error = function(e) NA)
cll_ci <- tryCatch(exp(confint(cll_cox)["chr13_loss", ]), error = function(e) c(NA, NA))
cll_p <- tryCatch(summary(cll_cox)$coefficients["chr13_loss", "Pr(>|z|)"], error = function(e) NA)

ak_direct_hr <- tryCatch(exp(coef(ak_cox)["chr13_loss"]), error = function(e) NA)
ak_direct_ci <- tryCatch(exp(confint(ak_cox)["chr13_loss", ]), error = function(e) c(NA, NA))
ak_direct_p <- tryCatch(summary(ak_cox)$coefficients["chr13_loss", "Pr(>|z|)"], error = function(e) NA)

cll_to_ak_hr <- tryCatch(exp(coef(ak_cox)["cll_status"]), error = function(e) NA)
cll_to_ak_ci <- tryCatch(exp(confint(ak_cox)["cll_status", ]), error = function(e) c(NA, NA))
cll_to_ak_p <- tryCatch(summary(ak_cox)$coefficients["cll_status", "Pr(>|z|)"], error = function(e) NA)

total_hr <- tryCatch(exp(coef(ak_total_cox)["chr13_loss"]), error = function(e) NA)
total_ci <- tryCatch(exp(confint(ak_total_cox)["chr13_loss", ]), error = function(e) c(NA, NA))
total_p <- tryCatch(summary(ak_total_cox)$coefficients["chr13_loss", "Pr(>|z|)"], error = function(e) NA)

# ==========================================
# LOGISTIC REGRESSION DIAGNOSTICS
# ==========================================

cat("\n=== 3. LOGISTIC REGRESSION DIAGNOSTICS ===\n\n")

cat("Mediator Model (chr13_loss → CLL):\n")
mediator_model <- glm(cll_status ~ chr13_loss + age + sex + n_baseline_diseases, 
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

cat("\nOutcome Model (chr13_loss + CLL → Actinic Keratosis):\n")
outcome_model <- glm(ak_status ~ chr13_loss + cll_status + age + sex + n_baseline_diseases, 
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
      treat = "chr13_loss",
      mediator = "cll_status",
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
  Effect = c("Total Effect (13- → AK)",
             "Direct Effect (13- → AK, adjusted for CLL)", 
             "Pathway A (13- → CLL)",
             "Pathway B (CLL → AK)"),
  Hazard_Ratio = c(round(total_hr, 3),
                   round(ak_direct_hr, 3),
                   round(cll_hr, 3),
                   round(cll_to_ak_hr, 3)),
  CI_Lower = c(round(total_ci[1], 3),
               round(ak_direct_ci[1], 3), 
               round(cll_ci[1], 3),
               round(cll_to_ak_ci[1], 3)),
  CI_Upper = c(round(total_ci[2], 3),
               round(ak_direct_ci[2], 3),
               round(cll_ci[2], 3), 
               round(cll_to_ak_ci[2], 3)),
  P_Value = c(ifelse(is.na(total_p), "NA", ifelse(total_p < 0.001, "< 0.001", round(total_p, 3))),
              ifelse(is.na(ak_direct_p), "NA", ifelse(ak_direct_p < 0.001, "< 0.001", round(ak_direct_p, 3))),
              ifelse(is.na(cll_p), "NA", ifelse(cll_p < 0.001, "< 0.001", round(cll_p, 3))),
              ifelse(is.na(cll_to_ak_p), "NA", ifelse(cll_to_ak_p < 0.001, "< 0.001", round(cll_to_ak_p, 3)))),
  Significance = c(ifelse(is.na(total_p), "NA", ifelse(total_p < 0.05, "Significant", "Not Significant")),
                   ifelse(is.na(ak_direct_p), "NA", ifelse(ak_direct_p < 0.05, "Significant", "Not Significant")),
                   ifelse(is.na(cll_p), "NA", ifelse(cll_p < 0.05, "Significant", "Not Significant")), 
                   ifelse(is.na(cll_to_ak_p), "NA", ifelse(cll_to_ak_p < 0.05, "Significant", "Not Significant")))
)

print(survival_results, row.names = FALSE)

# Save Cox model results
write.csv(survival_results, "~/Desktop/UK BB/chr13_loss_cll_ak_survival_mediation_results.csv", row.names = FALSE)

# ==========================================
# CALCULATE BASELINE RISKS AND CREATE RISK TABLE
# ==========================================

# Calculate baseline AK risk in No_mCA group
baseline_ak_risk <- analysis_data %>%
  filter(No_mCA == 1) %>%
  summarise(baseline_risk = mean(ak_status) * 100) %>%
  pull(baseline_risk)

# Calculate actual chr13_loss AK risk
chr13_loss_ak_risk <- analysis_data %>%
  filter(chr13_loss == 1) %>%
  summarise(chr13_loss_risk = mean(ak_status) * 100) %>%
  pull(chr13_loss_risk)

# Extract mediation effects if available
if (!is.null(mediation_results)) {
  total_effect_pct <- mediation_results$tau.coef * 100
  acme_effect_pct <- mediation_results$d0 * 100
  ade_effect_pct <- mediation_results$z0 * 100
  
  risk_table <- data.frame(
    Group = c("No_mCA (baseline)", 
              "chr13_loss (observed)", 
              "chr13_loss (through CLL pathway)",
              "chr13_loss (direct pathway)"),
    AK_Risk_Percent = c(
      round(baseline_ak_risk, 2),
      round(chr13_loss_ak_risk, 2),
      round(baseline_ak_risk + acme_effect_pct, 2),
      round(baseline_ak_risk + ade_effect_pct, 2)
    ),
    Risk_Difference = c(
      "Reference",
      paste0(round(chr13_loss_ak_risk - baseline_ak_risk, 2), "%"),
      paste0(round(acme_effect_pct, 2), "% (through CLL)"),
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
    Group = c("No_mCA (baseline)", "chr13_loss (observed)"),
    AK_Risk_Percent = c(round(baseline_ak_risk, 2), round(chr13_loss_ak_risk, 2)),
    Risk_Difference = c("Reference", paste0(round(chr13_loss_ak_risk - baseline_ak_risk, 2), "%")),
    P_Value = c("-", "-")
  )
  acme_effect_pct <- NA
  ade_effect_pct <- NA
}

cat("\n=== ACTINIC KERATOSIS RISK TABLE (chr13 loss → CLL → AK) ===\n")
print(risk_table, row.names = FALSE)

cat("\n=== SUMMARY ===\n")
cat("Baseline Actinic Keratosis risk (No_mCA):", round(baseline_ak_risk, 2), "%\n")
cat("Observed chr13_loss Actinic Keratosis risk:", round(chr13_loss_ak_risk, 2), "%\n")
cat("Total risk difference:", round(chr13_loss_ak_risk - baseline_ak_risk, 2), "percentage points\n")
if (!is.na(acme_effect_pct)) {
  cat("Mediated through CLL:", round(acme_effect_pct, 2), "percentage points\n")
  cat("Direct effect:", round(ade_effect_pct, 2), "percentage points\n")
  if (!is.null(mediation_results) && !is.null(mediation_results$n0)) {
    cat("Proportion mediated:", round(mediation_results$n0 * 100, 1), "% of total effect\n")
  }
}

# Save risk table
write.csv(risk_table, "~/Desktop/UK BB/chr13_loss_cll_ak_risk_table.csv", row.names = FALSE)

# ==========================================
# FIGURE 1: MEDIATION PLOT (Binary mediation results)
# ==========================================

cat("\n=== GENERATING FIGURES ===\n\n")

if (!is.null(mediation_results)) {
  pdf("~/Desktop/UK BB/chr13_loss_cll_ak_mediation_plot.pdf", width = 8, height = 8)
  plot(mediation_results, 
       main = "Mediation of chr13 Loss Association with Actinic Keratosis\nthrough CLL Development",
       xlab = "Absolute Risk Increase (percentage points)")
  dev.off()
  cat("Figure 1 saved: chr13_loss_cll_ak_mediation_plot.pdf\n")
} else {
  cat("Figure 1 skipped: Mediation analysis did not complete\n")
}

# ==========================================
# FIGURE 2: KAPLAN-MEIER CURVE - chr13_loss → CLL (Pathway A)
# ==========================================

surv_cll <- survfit(Surv(time_to_cll, cll_event) ~ chr13_loss_factor, data = analysis_data)

km_cll <- ggsurvplot(
  surv_cll,
  data = analysis_data,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  palette = c("#1f78b4", "#e31a1c"),
  title = "Pathway A: chr13 Loss to CLL Development",
  xlab = "Years from mCA assessment",
  ylab = "CLL-Free Survival Probability",
  legend.title = "Group",
  legend.labs = c("No mCA", "13-"),
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

pdf("~/Desktop/UK BB/chr13_loss_cll_ak_KM_pathway_A.pdf", width = 7, height = 9)
print(km_cll)
dev.off()
cat("Figure 2 saved: chr13_loss_cll_ak_KM_pathway_A.pdf\n")

# ==========================================
# FIGURE 3: KAPLAN-MEIER CURVE - chr13_loss → AK (Total Effect)
# ==========================================

surv_ak_total <- survfit(Surv(time_to_ak, ak_event) ~ chr13_loss_factor, data = analysis_data)

km_ak_total <- ggsurvplot(
  surv_ak_total,
  data = analysis_data,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  palette = c("#1f78b4", "#e31a1c"),
  title = "Total Effect: chr13 Loss to Actinic Keratosis",
  xlab = "Years from mCA assessment",
  ylab = "Actinic Keratosis-Free Survival Probability",
  legend.title = "Group",
  legend.labs = c("No mCA", "13-"),
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

pdf("~/Desktop/UK BB/chr13_loss_cll_ak_KM_total_effect.pdf", width = 7, height = 9)
print(km_ak_total)
dev.off()
cat("Figure 3 saved: chr13_loss_cll_ak_KM_total_effect.pdf\n")

# ==========================================
# FIGURE 4: KAPLAN-MEIER CURVE - CLL → AK (Pathway B)
# ==========================================

surv_ak_cll <- survfit(Surv(time_to_ak, ak_event) ~ cll_factor, data = analysis_data)

km_ak_cll <- ggsurvplot(
  surv_ak_cll,
  data = analysis_data,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.col = "strata",
  palette = c("#1f78b4", "#e31a1c"),
  title = "Pathway B: CLL to Actinic Keratosis",
  xlab = "Years from mCA assessment",
  ylab = "Actinic Keratosis-Free Survival Probability",
  legend.title = "CLL Status",
  legend.labs = c("No CLL", "CLL"),
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

pdf("~/Desktop/UK BB/chr13_loss_cll_ak_KM_pathway_B.pdf", width = 7, height = 9)
print(km_ak_cll)
dev.off()
cat("Figure 4 saved: chr13_loss_cll_ak_KM_pathway_B.pdf\n")

# ==========================================
# FIGURE 5: FOREST PLOT - All Pathway HRs
# ==========================================

forest_data <- data.frame(
  Pathway = c("Total Effect\n(13- -> AK)", 
              "Direct Effect\n(13- -> AK | CLL)", 
              "Pathway A\n(13- -> CLL)", 
              "Pathway B\n(CLL -> AK)"),
  HR = c(total_hr, ak_direct_hr, cll_hr, cll_to_ak_hr),
  CI_Lower = c(total_ci[1], ak_direct_ci[1], cll_ci[1], cll_to_ak_ci[1]),
  CI_Upper = c(total_ci[2], ak_direct_ci[2], cll_ci[2], cll_to_ak_ci[2]),
  P_Value = c(total_p, ak_direct_p, cll_p, cll_to_ak_p)
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
      title = "Mediation Analysis: chr13 Loss to CLL to Actinic Keratosis",
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
  
  ggsave("~/Desktop/UK BB/chr13_loss_cll_ak_forest_plot.pdf", 
         forest_plot, width = 7, height = 7)
  cat("Figure 5 saved: chr13_loss_cll_ak_forest_plot.pdf\n")
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
  
  annotate("text", x = 1, y = 1, label = "chr13 Loss\n(13-)", size = 4.5, fontface = "bold") +
  annotate("text", x = 4, y = 1, label = "CLL", size = 4.5, fontface = "bold") +
  annotate("text", x = 7, y = 1, label = "Actinic\nKeratosis", size = 4, fontface = "bold") +
  
  annotate("segment", x = 2, xend = 2.9, y = 1, yend = 1, 
           arrow = arrow(length = unit(0.3, "cm")), linewidth = 1.2, color = "#E41A1C") +
  annotate("text", x = 2.45, y = 1.15, 
           label = ifelse(is.na(cll_hr), "Pathway A\nHR=NA", 
                          sprintf("Pathway A\nHR=%.2f%s", cll_hr, ifelse(is.na(cll_p) | cll_p >= 0.05, " (NS)", "***"))), 
           size = 3.5, color = "#E41A1C", fontface = "bold") +
  
  annotate("segment", x = 5, xend = 5.9, y = 1, yend = 1, 
           arrow = arrow(length = unit(0.3, "cm")), linewidth = 1.2, color = "#E41A1C") +
  annotate("text", x = 5.45, y = 1.15, 
           label = ifelse(is.na(cll_to_ak_hr), "Pathway B\nHR=NA",
                          sprintf("Pathway B\nHR=%.2f%s", cll_to_ak_hr, ifelse(is.na(cll_to_ak_p) | cll_to_ak_p >= 0.05, " (NS)", "***"))), 
           size = 3.5, color = "#E41A1C", fontface = "bold") +
  
  annotate("curve", x = 2, xend = 5.9, y = 0.75, yend = 0.75,
           curvature = -0.3, arrow = arrow(length = unit(0.3, "cm")), 
           linewidth = 1, color = "#377EB8", linetype = "dashed") +
  annotate("text", x = 4, y = 0.45, 
           label = ifelse(is.na(ak_direct_hr), "Direct Effect\nHR=NA",
                          sprintf("Direct Effect\nHR=%.2f%s", ak_direct_hr, 
                                  ifelse(is.na(ak_direct_p) | ak_direct_p >= 0.05, " (NS)", "***"))), 
           size = 3.5, color = "#377EB8") +
  
  annotate("text", x = 4, y = 0.15, 
           label = ifelse(is.na(total_hr), "Total Effect: HR=NA",
                          sprintf("Total Effect: HR=%.2f%s\n*** p<0.05, NS = Not Significant", 
                                  total_hr, ifelse(is.na(total_p) | total_p >= 0.05, " (NS)", "***"))), 
           size = 3.5, color = "gray30") +
  
  labs(title = "Mediation Pathway: chr13 Loss to CLL to Actinic Keratosis",
       subtitle = ifelse(is.na(total_p) | total_p >= 0.05, "Null Total Effect", "Significant Mediation")) +
  
  theme_void() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  coord_cartesian(xlim = c(-0.5, 8.5), ylim = c(0, 1.5))

ggsave("~/Desktop/UK BB/chr13_loss_cll_ak_pathway_diagram.pdf", 
       pathway_diagram, width = 10, height = 5)
cat("Figure 6 saved: chr13_loss_cll_ak_pathway_diagram.pdf\n")

# ==========================================
# MEDIATION INTERPRETATION
# ==========================================

cat("\n")
cat(rep("=", 80), "\n", sep = "")
cat("MEDIATION INTERPRETATION\n")
cat(rep("=", 80), "\n\n", sep = "")

if(is.na(cll_p) || is.na(cll_to_ak_p)) {
  cat("RESULT: Unable to fully assess mediation due to model estimation issues\n")
  cat("  -> This is likely due to very small sample sizes or rare events\n")
} else if(cll_p < 0.05 && cll_to_ak_p < 0.05) {
  if(is.na(total_p) || total_p >= 0.05) {
    cat("RESULT: NULL TOTAL EFFECT despite significant individual pathways\n")
    cat("  -> chr13_loss predisposes to CLL (Pathway A significant)\n")
    cat("  -> CLL significantly increases Actinic Keratosis risk (Pathway B significant)\n")
    cat("  -> BUT: chr13_loss does NOT significantly increase overall AK risk\n")
    cat("  -> INTERPRETATION: Low CLL incidence among chr13_loss carriers limits\n")
    cat("                    the population-level impact of the mediated pathway\n")
  } else if(is.na(ak_direct_p) || ak_direct_p >= 0.05) {
    cat("RESULT: Full mediation - chr13_loss -> AK association is fully mediated through CLL\n")
  } else {
    cat("RESULT: Partial mediation - chr13_loss -> AK has both direct and CLL-mediated pathways\n")
  }
} else if(!is.na(cll_p) && cll_p < 0.05 && (is.na(cll_to_ak_p) || cll_to_ak_p >= 0.05)) {
  cat("RESULT: No mediation - CLL does not significantly affect Actinic Keratosis risk\n")
} else if(is.na(cll_p) || cll_p >= 0.05) {
  cat("RESULT: No mediation - chr13_loss does not significantly affect CLL risk\n") 
} else {
  cat("RESULT: Complex pattern - see individual pathway significance\n")
}

cat("\nchr13_loss -> CLL pathway:", 
    ifelse(is.na(cll_p), "NA", ifelse(cll_p < 0.05, "SIGNIFICANT", "NOT SIGNIFICANT")), "\n")
cat("CLL -> Actinic Keratosis pathway:", 
    ifelse(is.na(cll_to_ak_p), "NA", ifelse(cll_to_ak_p < 0.05, "SIGNIFICANT", "NOT SIGNIFICANT")), "\n")
cat("chr13_loss -> Actinic Keratosis total effect:", 
    ifelse(is.na(total_p), "NA", ifelse(total_p < 0.05, "SIGNIFICANT", "NOT SIGNIFICANT")), "\n")
cat("chr13_loss -> Actinic Keratosis direct effect:", 
    ifelse(is.na(ak_direct_p), "NA", ifelse(ak_direct_p < 0.05, "SIGNIFICANT", "NOT SIGNIFICANT")), "\n")

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
if (!is.na(cll_hr) && cll_hr > 50) warnings_count <- warnings_count + 0.5
if (!is.na(total_hr) && total_hr > 50) warnings_count <- warnings_count + 0.5
if (sum(analysis_data$chr13_loss) < 100) warnings_count <- warnings_count + 1

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
  cat("     1. chr13_loss sample size (n=", sum(analysis_data$chr13_loss), ")\n")
  cat("     2. Consider this analysis EXPLORATORY if sample size is small\n")
  cat("     3. Results may be underpowered\n")
  cat("     4. Consider combining with other chromosomal abnormalities\n")
  cat("     5. CLL events in exposed:", sum(analysis_data$cll_event[analysis_data$chr13_loss == 1]), "\n")
}

# ==========================================
# FILES SAVED SUMMARY
# ==========================================

cat("\n")
cat(rep("=", 80), "\n", sep = "")
cat("FILES SAVED\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("=== CSV FILES ===\n")
cat("1. chr13_loss_cll_ak_survival_mediation_results.csv - Cox model results\n")
cat("2. chr13_loss_cll_ak_risk_table.csv - Risk comparison table\n")

cat("\n=== PDF FIGURES ===\n")
cat("1. chr13_loss_cll_ak_mediation_plot.pdf - Binary mediation analysis plot\n")
cat("2. chr13_loss_cll_ak_KM_pathway_A.pdf - Kaplan-Meier: chr13 loss -> CLL\n")
cat("3. chr13_loss_cll_ak_KM_total_effect.pdf - Kaplan-Meier: chr13 loss -> Actinic Keratosis\n")
cat("4. chr13_loss_cll_ak_KM_pathway_B.pdf - Kaplan-Meier: CLL -> Actinic Keratosis\n")
cat("5. chr13_loss_cll_ak_forest_plot.pdf - Forest plot of all pathway HRs\n")
cat("6. chr13_loss_cll_ak_pathway_diagram.pdf - Visual pathway diagram\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
