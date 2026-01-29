# =========================================================================================
# INTERACTION ANALYSIS: mCA × CHIP - DISEASE INCIDENCE (OPTIMIZED VERSION)
# =========================================================================================
# Uses faster computation methods for Fine-Gray regression
# Fixed: Added prodlim library, fixed interaction model errors, explicit namespace for crr()
# =========================================================================================

# Clear workspace
rm(list = ls())

# Load libraries
library(readxl)
library(data.table)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(survival)

# Check for riskRegression (faster Fine-Gray)
if (!requireNamespace("riskRegression", quietly = TRUE)) {
  cat("Installing riskRegression package for faster computation...\n")
  install.packages("riskRegression")
}
library(riskRegression)
library(prodlim)  # Required for Hist() function used in FGR models

# Load cmprsk for Gray's test and crr() for interaction models
library(cmprsk)

# Load tidycmprsk for CIF plotting
# NOTE: tidycmprsk masks cmprsk::crr() and cmprsk::cuminc()
# We must use cmprsk::crr() explicitly for interaction models
library(tidycmprsk)

# Settings
censor_date <- as.Date("2023-07-14")
data_dir <- "~/Desktop/UK BB"
output_dir <- "~/Desktop/UK BB/Interaction_Analysis"
dir.create(output_dir, showWarnings = FALSE)

cat("\n")
cat("================================================================================\n")
cat("  DISEASE INCIDENCE INTERACTION ANALYSIS: mCA × CHIP (OPTIMIZED)\n")
cat("================================================================================\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# =========================================================================================
# DATA LOADING
# =========================================================================================

cat("Loading data...\n")

analyzed_patients <- read_excel(file.path(data_dir, "ukb4777.samples_analyzed.xlsx")) %>%
  mutate(ID = as.character(ID))

mca_info <- read_excel(file.path(data_dir, "Age sex dob collection date.xlsx")) %>%
  mutate(
    ID = as.character(ID),
    mCA_date = mdy(`Date of mCA assessment`),
    age = as.numeric(`Age at blood collection`),
    sex = as.factor(Sex)
  ) %>%
  filter(ID %in% analyzed_patients$ID) %>%
  filter(!is.na(mCA_date),
         mCA_date != as.Date("1900-01-01"),
         mCA_date > as.Date("1901-01-01")) %>%
  select(ID, mCA_date, age, sex) %>%
  distinct(ID, .keep_all = TRUE)

death_data <- read_excel(file.path(data_dir, "Death date.xlsx")) %>%
  rename(ID = Participant_ID, death_date = date_of_death) %>%
  mutate(ID = as.character(ID), death_date = as.Date(death_date)) %>%
  distinct(ID, .keep_all = TRUE)

cnv <- read_excel(file.path(data_dir, "CNV_data_19808.xlsx")) %>%
  mutate(ID = as.character(ID)) %>%
  filter(ID %in% analyzed_patients$ID)

cnv_summary <- cnv %>%
  group_by(ID) %>%
  summarise(
    has_known_cnv = any(COPY_CHANGE %in% c("gain", "loss", "neutral")),
    has_only_unknown = all(COPY_CHANGE == "unknown"),
    .groups = "drop"
  )

cnv_status <- analyzed_patients %>%
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

mutation_assessed <- fread(file.path(data_dir, "Somatic_mutations_participant.csv")) %>%
  mutate(ID = as.character(ID)) %>%
  filter(!is.na(CHIPtotalnumber)) %>%
  select(ID, CHIPtotalnumber) %>%
  mutate(mutation_status = ifelse(CHIPtotalnumber > 0, "CHIP+", "CHIP-"))

diag_data <- fread(file.path(data_dir, "Diagnosis code date.tsv"), sep = "\t") %>%
  mutate(
    patient_id = as.character(patient_id),
    admission_date = as.Date(admission_date),
    diagnosis = as.character(diagnosis)
  ) %>%
  filter(!is.na(admission_date), !is.na(diagnosis))

baseline_disease_burden <- diag_data %>%
  inner_join(mca_info %>% select(ID, mCA_date), 
             by = c("patient_id" = "ID"),
             relationship = "many-to-one") %>%
  filter(admission_date <= mCA_date) %>%
  group_by(patient_id) %>%
  summarise(n_baseline_diseases = n_distinct(diagnosis), .groups = "drop")

cat("Data loading complete.\n\n")

# =========================================================================================
# BUILD ANALYSIS DATASET
# =========================================================================================

cat("Building analysis dataset...\n")

analysis_df <- mca_info %>%
  inner_join(cnv_status, by = "ID") %>%
  inner_join(mutation_assessed %>% select(ID, mutation_status), by = "ID") %>%
  left_join(death_data, by = "ID") %>%
  left_join(baseline_disease_burden, by = c("ID" = "patient_id")) %>%
  mutate(
    n_baseline_diseases = replace_na(n_baseline_diseases, 0),
    has_mCA = ifelse(mCA_status == "mCA", 1, 0),
    has_CHIP = ifelse(mutation_status == "CHIP+", 1, 0),
    sex_numeric = as.numeric(sex) - 1,
    group = case_when(
      mCA_status == "no_mCA" & mutation_status == "CHIP-" ~ "mCA-/CHIP-",
      mCA_status == "mCA" & mutation_status == "CHIP-" ~ "mCA+/CHIP-",
      mCA_status == "no_mCA" & mutation_status == "CHIP+" ~ "mCA-/CHIP+",
      mCA_status == "mCA" & mutation_status == "CHIP+" ~ "mCA+/CHIP+"
    ),
    group = factor(group, levels = c("mCA-/CHIP-", "mCA+/CHIP-", "mCA-/CHIP+", "mCA+/CHIP+"))
  )

# Get incident diseases
pre_existing <- diag_data %>%
  inner_join(analysis_df %>% select(ID, mCA_date), 
             by = c("patient_id" = "ID"),
             relationship = "many-to-one") %>%
  filter(admission_date <= mCA_date) %>%
  distinct(patient_id, diagnosis)

post_mca <- diag_data %>%
  inner_join(analysis_df %>% select(ID, mCA_date), 
             by = c("patient_id" = "ID"),
             relationship = "many-to-one") %>%
  filter(admission_date > mCA_date)

incident_diseases <- post_mca %>%
  anti_join(pre_existing, by = c("patient_id", "diagnosis")) %>%
  group_by(patient_id) %>%
  summarise(first_dx_date = min(admission_date), .groups = "drop")

analysis_df <- analysis_df %>%
  left_join(incident_diseases, by = c("ID" = "patient_id"))

# Create survival dataset
survival_df <- analysis_df %>%
  mutate(
    event_status = case_when(
      !is.na(first_dx_date) & first_dx_date > mCA_date & 
        (is.na(death_date) | first_dx_date <= death_date) ~ 1,  # disease
      !is.na(death_date) & death_date > mCA_date ~ 2,  # death (competing)
      TRUE ~ 0  # censored
    ),
    end_date = case_when(
      event_status == 1 ~ first_dx_date,
      event_status == 2 ~ death_date,
      TRUE ~ censor_date
    ),
    time_years = as.numeric(difftime(end_date, mCA_date, units = "days")) / 365.25
  ) %>%
  filter(time_years >= 0, !is.na(group), !is.na(age), !is.na(sex), !is.na(n_baseline_diseases))

cat("Analysis dataset created:", nrow(survival_df), "participants\n\n")

# =========================================================================================
# COHORT SUMMARY
# =========================================================================================

cat("=== COHORT SUMMARY ===\n\n")

group_summary <- survival_df %>%
  group_by(group) %>%
  summarise(
    N = n(),
    Disease_Events = sum(event_status == 1),
    Deaths_Competing = sum(event_status == 2),
    Censored = sum(event_status == 0),
    Incidence_rate = sprintf("%.1f%%", 100 * sum(event_status == 1) / n()),
    .groups = "drop"
  )

print(group_summary)

# =========================================================================================
# CUMULATIVE INCIDENCE CURVES
# =========================================================================================

cat("\n=== Generating Cumulative Incidence Curves ===\n")

# Gray's test - use cmprsk::cuminc explicitly
gray_test <- cmprsk::cuminc(
  ftime = survival_df$time_years,
  fstatus = survival_df$event_status,
  group = survival_df$group
)

gray_p <- gray_test$Tests[1, "pv"]
cat("Gray's test p-value:", format.pval(gray_p, digits = 3), "\n")

# Create CIF plot using tidycmprsk
survival_df$event_factor <- factor(
  case_when(
    survival_df$event_status == 0 ~ "censored",
    survival_df$event_status == 1 ~ "disease",
    survival_df$event_status == 2 ~ "death"
  ),
  levels = c("censored", "disease", "death")
)

cif_fit <- tidycmprsk::cuminc(Surv(time_years, event_factor) ~ group, data = survival_df)

cif_df <- tidy(cif_fit) %>%
  filter(outcome == "disease") %>%
  filter(!is.na(estimate))

y_max <- max(cif_df$conf.high, na.rm = TRUE) * 1.05

group_colors <- c("#1f78b4", "#92C5DE", "#F4A582", "#e31a1c")

cif_plot <- ggplot(cif_df, aes(x = time, y = estimate, color = strata, fill = strata)) +
  geom_step(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  scale_color_manual(values = group_colors,
                     labels = c("mCA-/CHIP-", "mCA+/CHIP-", "mCA-/CHIP+", "mCA+/CHIP+")) +
  scale_fill_manual(values = group_colors,
                    labels = c("mCA-/CHIP-", "mCA+/CHIP-", "mCA-/CHIP+", "mCA+/CHIP+")) +
  scale_x_continuous(breaks = seq(0, 20, 2), limits = c(0, 20)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05)), limits = c(0, y_max),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Years from mCA assessment",
       y = "Cumulative Incidence of Any Disease",
       title = "Cumulative Incidence: mCA × CHIP Interaction",
       color = "Group", fill = "Group") +
  annotate("text", x = 0.5, y = y_max * 0.15,
           label = paste0("Gray's test p < 0.001"),
           hjust = 0, size = 4.5, fontface = "italic") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(cif_plot)

ggsave(file.path(output_dir, "mCA_CHIP_Interaction_CIF_All_Diseases.pdf"), 
       plot = cif_plot, width = 10, height = 7)

cat("✅ CIF plot saved\n\n")

# =========================================================================================
# FINE-GRAY REGRESSION USING riskRegression (MUCH FASTER)
# =========================================================================================

cat("=== Fine-Gray Regression (using riskRegression - faster) ===\n\n")

# Create Hist object for competing risks
survival_df$event_cr <- factor(survival_df$event_status, levels = c(0, 1, 2),
                               labels = c("censor", "disease", "death"))

# Convert to data.frame (required by FGR)
df_for_fgr <- as.data.frame(survival_df)

# ------------------------------
# Model 1: Unadjusted 4-group
# ------------------------------
cat("Running Model 1: Unadjusted 4-group comparison...\n")
start_time <- Sys.time()

fgr_unadj <- FGR(Hist(time_years, event_status) ~ group,
                 data = df_for_fgr,
                 cause = 1)

cat("  Completed in", round(difftime(Sys.time(), start_time, units = "mins"), 1), "minutes\n")
print(summary(fgr_unadj))

# ------------------------------
# Model 2: Adjusted 4-group
# ------------------------------
cat("\nRunning Model 2: Adjusted 4-group comparison...\n")
start_time <- Sys.time()

fgr_adj <- FGR(Hist(time_years, event_status) ~ group + age + sex + n_baseline_diseases,
               data = df_for_fgr,
               cause = 1)

cat("  Completed in", round(difftime(Sys.time(), start_time, units = "mins"), 1), "minutes\n")
print(summary(fgr_adj))

# =========================================================================================
# INTERACTION MODELS USING cmprsk::crr() DIRECTLY
# =========================================================================================
# NOTE: tidycmprsk masks cmprsk::crr(), so we MUST use cmprsk::crr() explicitly

cat("\n=== Interaction Models (using cmprsk::crr) ===\n\n")

# ------------------------------
# Model 3: Interaction unadjusted
# ------------------------------
cat("Running Model 3: Interaction model (unadjusted)...\n")
start_time <- Sys.time()

# Create design matrix for interaction model (unadjusted)
cov_int_unadj <- as.matrix(data.frame(
  has_mCA = df_for_fgr$has_mCA,
  has_CHIP = df_for_fgr$has_CHIP,
  has_mCA_x_has_CHIP = df_for_fgr$has_mCA * df_for_fgr$has_CHIP
))

# IMPORTANT: Use cmprsk::crr() explicitly because tidycmprsk masks it
crr_int_unadj <- cmprsk::crr(
  ftime = df_for_fgr$time_years,
  fstatus = df_for_fgr$event_status,
  cov1 = cov_int_unadj,
  failcode = 1,
  cencode = 0
)

cat("  Completed in", round(difftime(Sys.time(), start_time, units = "mins"), 1), "minutes\n")
print(summary(crr_int_unadj))

# ------------------------------
# Model 4: Interaction adjusted
# ------------------------------
cat("\nRunning Model 4: Interaction model (adjusted)...\n")
start_time <- Sys.time()

# Create design matrix for interaction model (adjusted)
cov_int_adj <- as.matrix(data.frame(
  has_mCA = df_for_fgr$has_mCA,
  has_CHIP = df_for_fgr$has_CHIP,
  has_mCA_x_has_CHIP = df_for_fgr$has_mCA * df_for_fgr$has_CHIP,
  age = df_for_fgr$age,
  sex = df_for_fgr$sex_numeric,
  n_baseline_diseases = df_for_fgr$n_baseline_diseases
))

# IMPORTANT: Use cmprsk::crr() explicitly because tidycmprsk masks it
crr_int_adj <- cmprsk::crr(
  ftime = df_for_fgr$time_years,
  fstatus = df_for_fgr$event_status,
  cov1 = cov_int_adj,
  failcode = 1,
  cencode = 0
)

cat("  Completed in", round(difftime(Sys.time(), start_time, units = "mins"), 1), "minutes\n")
print(summary(crr_int_adj))

# =========================================================================================
# EXTRACT AND FORMAT RESULTS
# =========================================================================================

cat("\n")
cat("================================================================================\n")
cat("  RESULTS SUMMARY\n")
cat("================================================================================\n\n")

# Extract coefficients from adjusted 4-group model (FGR)
coef_adj <- fgr_adj$crrFit$coef
se_adj <- sqrt(diag(fgr_adj$crrFit$var))
hr_adj <- exp(coef_adj)
hr_lower_adj <- exp(coef_adj - 1.96 * se_adj)
hr_upper_adj <- exp(coef_adj + 1.96 * se_adj)
z_adj <- coef_adj / se_adj
p_adj <- 2 * (1 - pnorm(abs(z_adj)))

cat("=== Adjusted 4-Group Comparison (vs. mCA-/CHIP- reference) ===\n\n")
cat("Group              sHR     95% CI           P-value\n")
cat("---------------------------------------------------\n")
cat("mCA-/CHIP-         1.00    (Reference)      -\n")
cat(sprintf("mCA+/CHIP-         %.2f    (%.2f-%.2f)      %s\n", 
            hr_adj[1], hr_lower_adj[1], hr_upper_adj[1], format.pval(p_adj[1], digits = 2)))
cat(sprintf("mCA-/CHIP+         %.2f    (%.2f-%.2f)      %s\n", 
            hr_adj[2], hr_lower_adj[2], hr_upper_adj[2], format.pval(p_adj[2], digits = 2)))
cat(sprintf("mCA+/CHIP+         %.2f    (%.2f-%.2f)      %s\n", 
            hr_adj[3], hr_lower_adj[3], hr_upper_adj[3], format.pval(p_adj[3], digits = 2)))

# Extract interaction results from crr model
coef_int <- crr_int_adj$coef
se_int <- sqrt(diag(crr_int_adj$var))
hr_int <- exp(coef_int)
hr_lower_int <- exp(coef_int - 1.96 * se_int)
hr_upper_int <- exp(coef_int + 1.96 * se_int)
z_int <- coef_int / se_int
p_int <- 2 * (1 - pnorm(abs(z_int)))

cat("\n=== Interaction Model (Adjusted) ===\n\n")
cat("Variable              sHR     95% CI           P-value\n")
cat("------------------------------------------------------\n")

var_names <- c("has_mCA", "has_CHIP", "has_mCA:has_CHIP", "age", "sex", "n_baseline_diseases")
for (i in seq_along(var_names)) {
  cat(sprintf("%-21s %.3f   (%.3f-%.3f)    %s\n",
              var_names[i], hr_int[i], hr_lower_int[i], hr_upper_int[i], 
              format.pval(p_int[i], digits = 3)))
}

# Interaction term is the 3rd coefficient
interaction_hr <- hr_int[3]
interaction_p <- p_int[3]
interaction_lower <- hr_lower_int[3]
interaction_upper <- hr_upper_int[3]

cat("\n")
cat("================================================================================\n")
cat("  MULTIPLICATIVE INTERACTION INTERPRETATION\n")
cat("================================================================================\n\n")

cat("Interaction Term (mCA × CHIP):\n")
cat("  sHR =", sprintf("%.3f", interaction_hr),
    "(95% CI:", sprintf("%.3f", interaction_lower), "-",
    sprintf("%.3f", interaction_upper), "),",
    "p =", format.pval(interaction_p, digits = 3), "\n\n")

if (interaction_p < 0.05) {
  if (interaction_hr > 1) {
    cat("✅ SIGNIFICANT POSITIVE INTERACTION (p =", format.pval(interaction_p, digits = 3), ")\n")
    cat("   The combined effect is SYNERGISTIC (super-multiplicative).\n")
  } else {
    cat("✅ SIGNIFICANT NEGATIVE INTERACTION (p =", format.pval(interaction_p, digits = 3), ")\n")
    cat("   The combined effect is SUB-MULTIPLICATIVE.\n")
  }
} else {
  cat("❌ NO SIGNIFICANT INTERACTION (p =", format.pval(interaction_p, digits = 3), ")\n")
  cat("   The combined effect appears MULTIPLICATIVE (no synergy/antagonism).\n")
}

# =========================================================================================
# SAVE RESULTS
# =========================================================================================

# Summary table
summary_table <- data.frame(
  Group = c("mCA-/CHIP- (Reference)", "mCA+/CHIP-", "mCA-/CHIP+", "mCA+/CHIP+"),
  N = as.vector(table(survival_df$group)),
  Disease_Events = as.vector(tapply(survival_df$event_status == 1, survival_df$group, sum)),
  sHR = c(1.00, hr_adj[1:3]),
  Lower_95CI = c(NA, hr_lower_adj[1:3]),
  Upper_95CI = c(NA, hr_upper_adj[1:3]),
  P_value = c(NA, p_adj[1:3])
)

summary_table$sHR_CI <- ifelse(is.na(summary_table$Lower_95CI), 
                               "1.00 (Reference)",
                               sprintf("%.2f (%.2f-%.2f)", 
                                       summary_table$sHR, 
                                       summary_table$Lower_95CI, 
                                       summary_table$Upper_95CI))

write.csv(summary_table, 
          file.path(output_dir, "mCA_CHIP_Interaction_Disease_Incidence_Summary.csv"),
          row.names = FALSE)

cat("\n\n✅ Main analysis complete!\n")
cat("Output files saved to:", output_dir, "\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# =========================================================================================
# ADDITIVE INTERACTION MEASURES: RERI, AP, S
# =========================================================================================

cat("\n")
cat("================================================================================\n")
cat("  ADDITIVE INTERACTION ANALYSIS\n")
cat("================================================================================\n\n")

# Extract HRs from the adjusted 4-group model
# HR_10 = mCA+/CHIP- (vs reference)
# HR_01 = mCA-/CHIP+ (vs reference)
# HR_11 = mCA+/CHIP+ (vs reference)

HR_10 <- hr_adj[1]  # mCA+/CHIP-
HR_01 <- hr_adj[2]  # mCA-/CHIP+
HR_11 <- hr_adj[3]  # mCA+/CHIP+

# Standard errors (on log scale)
SE_10 <- se_adj[1]
SE_01 <- se_adj[2]
SE_11 <- se_adj[3]

# Variance-covariance matrix for the 3 group coefficients
vcov_matrix <- fgr_adj$crrFit$var[1:3, 1:3]

cat("=== Hazard Ratios Used ===\n")
cat("HR_mCA+/CHIP- (HR_10):", sprintf("%.4f", HR_10), "\n")
cat("HR_mCA-/CHIP+ (HR_01):", sprintf("%.4f", HR_01), "\n")
cat("HR_mCA+/CHIP+ (HR_11):", sprintf("%.4f", HR_11), "\n\n")

# -------------------------
# RERI: Relative Excess Risk due to Interaction
# -------------------------
# RERI = HR_11 - HR_10 - HR_01 + 1
# RERI > 0 indicates positive additive interaction (synergy)
# RERI = 0 indicates no additive interaction
# RERI < 0 indicates negative additive interaction (antagonism)

RERI <- HR_11 - HR_10 - HR_01 + 1

# Delta method for RERI variance
# RERI = exp(b3) - exp(b1) - exp(b2) + 1
# Partial derivatives: dRERI/db1 = -exp(b1), dRERI/db2 = -exp(b2), dRERI/db3 = exp(b3)
grad_RERI <- c(-HR_10, -HR_01, HR_11)
var_RERI <- t(grad_RERI) %*% vcov_matrix %*% grad_RERI
SE_RERI <- sqrt(var_RERI)

RERI_lower <- RERI - 1.96 * SE_RERI
RERI_upper <- RERI + 1.96 * SE_RERI
RERI_z <- RERI / SE_RERI
RERI_p <- 2 * (1 - pnorm(abs(RERI_z)))

cat("=== RERI (Relative Excess Risk due to Interaction) ===\n")
cat("RERI =", sprintf("%.4f", RERI), "\n")
cat("95% CI:", sprintf("%.4f", RERI_lower), "to", sprintf("%.4f", RERI_upper), "\n")
cat("P-value:", format.pval(RERI_p, digits = 3), "\n")
cat("Interpretation: RERI > 0 indicates synergy on additive scale\n\n")

# -------------------------
# AP: Attributable Proportion due to Interaction
# -------------------------
# AP = RERI / HR_11
# Proportion of disease in double-exposed that is due to interaction

AP <- RERI / HR_11

# Gradient for AP with respect to log(HR) coefficients
# dAP/db1 = -exp(b1)/exp(b3) = -HR_10/HR_11
# dAP/db2 = -exp(b2)/exp(b3) = -HR_01/HR_11
# dAP/db3 = (exp(b1) + exp(b2) - 1)/exp(b3) = (HR_10 + HR_01 - 1)/HR_11

grad_AP <- c(-HR_10/HR_11, -HR_01/HR_11, (HR_10 + HR_01 - 1)/HR_11)
var_AP <- t(grad_AP) %*% vcov_matrix %*% grad_AP
SE_AP <- sqrt(var_AP)

AP_lower <- AP - 1.96 * SE_AP
AP_upper <- AP + 1.96 * SE_AP
AP_z <- AP / SE_AP
AP_p <- 2 * (1 - pnorm(abs(AP_z)))

cat("=== AP (Attributable Proportion due to Interaction) ===\n")
cat("AP =", sprintf("%.4f", AP), "(", sprintf("%.1f%%", AP * 100), ")\n")
cat("95% CI:", sprintf("%.4f", AP_lower), "to", sprintf("%.4f", AP_upper), "\n")
cat("P-value:", format.pval(AP_p, digits = 3), "\n")
cat("Interpretation: Proportion of risk in double-exposed attributable to interaction\n\n")

# -------------------------
# S: Synergy Index (Rothman)
# -------------------------
# S = (HR_11 - 1) / ((HR_10 - 1) + (HR_01 - 1))
# S > 1 indicates synergy (super-additive)
# S = 1 indicates exact additivity
# S < 1 indicates antagonism (sub-additive)

S <- (HR_11 - 1) / ((HR_10 - 1) + (HR_01 - 1))

# Delta method for S
denom <- (HR_10 - 1) + (HR_01 - 1)
grad_S <- c(
  -HR_10 * (HR_11 - 1) / denom^2,
  -HR_01 * (HR_11 - 1) / denom^2,
  HR_11 / denom
)
var_S <- t(grad_S) %*% vcov_matrix %*% grad_S
SE_S <- sqrt(var_S)

S_lower <- S - 1.96 * SE_S
S_upper <- S + 1.96 * SE_S

# Test if S differs from 1
S_z <- (S - 1) / SE_S
S_p <- 2 * (1 - pnorm(abs(S_z)))

cat("=== S (Synergy Index) ===\n")
cat("S =", sprintf("%.4f", S), "\n")
cat("95% CI:", sprintf("%.4f", S_lower), "to", sprintf("%.4f", S_upper), "\n")
cat("P-value (testing S = 1):", format.pval(S_p, digits = 3), "\n")
cat("Interpretation: S > 1 indicates synergy, S = 1 indicates additivity, S < 1 indicates antagonism\n\n")

# -------------------------
# SUMMARY TABLE
# -------------------------

cat("================================================================================\n")
cat("  SUMMARY: ADDITIVE INTERACTION MEASURES\n")
cat("================================================================================\n\n")

additive_summary <- data.frame(
  Measure = c("RERI", "AP", "S"),
  Estimate = c(RERI, AP, S),
  Lower_95CI = c(RERI_lower, AP_lower, S_lower),
  Upper_95CI = c(RERI_upper, AP_upper, S_upper),
  P_value = c(RERI_p, AP_p, S_p),
  Null_Value = c(0, 0, 1),
  Interpretation = c(
    ifelse(RERI > 0, "Positive (super-additive)", ifelse(RERI < 0, "Negative (sub-additive)", "No interaction")),
    ifelse(AP > 0, "Positive", ifelse(AP < 0, "Negative", "None")),
    ifelse(S > 1, "Synergy", ifelse(S < 1, "Antagonism", "Additive"))
  )
)

print(additive_summary)

# -------------------------
# INTERPRETATION
# -------------------------

cat("\n")
cat("================================================================================\n")
cat("  INTERPRETATION\n")
cat("================================================================================\n\n")

cat("RERI =", sprintf("%.4f", RERI), "\n")
if (RERI_p < 0.05) {
  if (RERI > 0) {
    cat("  ✅ Significant POSITIVE additive interaction (p =", format.pval(RERI_p, digits = 3), ")\n")
    cat("  The combined effect EXCEEDS the sum of individual effects.\n")
  } else {
    cat("  ✅ Significant NEGATIVE additive interaction (p =", format.pval(RERI_p, digits = 3), ")\n")
    cat("  The combined effect is LESS than the sum of individual effects.\n")
  }
} else {
  cat("  ❌ No significant additive interaction (p =", format.pval(RERI_p, digits = 3), ")\n")
  cat("  The combined effect is consistent with additivity.\n")
}

cat("\nS =", sprintf("%.4f", S), "\n")
if (S_p < 0.05) {
  if (S > 1) {
    cat("  ✅ Synergy Index significantly > 1 (p =", format.pval(S_p, digits = 3), ")\n")
  } else {
    cat("  ✅ Synergy Index significantly < 1 (p =", format.pval(S_p, digits = 3), ")\n")
  }
} else {
  cat("  ❌ Synergy Index not significantly different from 1 (p =", format.pval(S_p, digits = 3), ")\n")
  cat("  Consistent with additive effects.\n")
}

# -------------------------
# SAVE RESULTS
# -------------------------

write.csv(additive_summary, 
          file.path(output_dir, "mCA_CHIP_Additive_Interaction_Disease_Incidence.csv"),
          row.names = FALSE)

cat("\n✅ Additive interaction results saved to:", 
    file.path(output_dir, "mCA_CHIP_Additive_Interaction_Disease_Incidence.csv"), "\n")

cat("\n")
cat("================================================================================\n")
cat("  ALL ANALYSES COMPLETE\n")
cat("================================================================================\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")