# =========================================================================================
# MORTALITY ANALYSIS: mCA SUBTYPES vs NO mCA (OVERALL MORTALITY)
# =========================================================================================
# Compares each mCA subtype against no mCA controls for all-cause mortality
# Uses exact same filtering as the original competing risk analysis
# =========================================================================================

library(readxl)
library(data.table)
library(dplyr)
library(tidyr)
library(lubridate)
library(cmprsk)
library(tidycmprsk)

# =========================================================================================
# DATA FILTERING (matching methods section and incidence analysis code)
# =========================================================================================
cat("\n")
cat("================================================================================\n")
cat("  DATA FILTERING\n")
cat("================================================================================\n")

# ===========================
# Load analyzed samples
# ===========================
cat("Loading analyzed samples...\n")
analyzed <- read_excel("~/Desktop/UK BB/ukb4777.samples_analyzed.xlsx") %>%
  mutate(ID = as.character(ID))

cat("Total analyzed samples:", nrow(analyzed), "\n\n")

# ===========================
# Load sample info with mCA dates, age, and sex
# ===========================
cat("Loading participant data...\n")
mca_info_raw <- read_excel("~/Desktop/UK BB/Age sex dob collection date.xlsx") %>%
  mutate(ID = as.character(ID),
         mCA_date = mdy(`Date of mCA assessment`),
         age = as.numeric(`Age at blood collection`),
         sex = as.factor(Sex)) %>%
  filter(ID %in% analyzed$ID)

cat("Participant records (analyzed samples only):", nrow(mca_info_raw), "\n")

# Filter out invalid dates
cat("\nStep 1: Filtering invalid mCA assessment dates...\n")
cat("  Missing mCA dates:", sum(is.na(mca_info_raw$mCA_date)), "\n")
cat("  Dates = 1900-01-01:", sum(mca_info_raw$mCA_date == as.Date("1900-01-01"), na.rm = TRUE), "\n")
cat("  Dates <= 1901-01-01:", sum(mca_info_raw$mCA_date <= as.Date("1901-01-01"), na.rm = TRUE), "\n")

mca_info <- mca_info_raw %>%
  filter(!is.na(mCA_date),
         mCA_date != as.Date("1900-01-01"),
         mCA_date > as.Date("1901-01-01")) %>%
  select(ID, mCA_date, age, sex) %>%
  distinct(ID, .keep_all = TRUE)

cat("  Valid participants after date filtering:", nrow(mca_info), "\n")
cat("  Date range:", as.character(min(mca_info$mCA_date)), "to", as.character(max(mca_info$mCA_date)), "\n\n")

# ===========================
# Load CNV data (filter AFTER date filtering to match previous code)
# ===========================
cat("Step 2: Loading CNV data...\n")
cnv_all <- read_excel("~/Desktop/UK BB/CNV_data_19808.xlsx") %>%
  mutate(ID = as.character(ID)) %>%
  filter(ID %in% analyzed$ID)

cat("  Total CNV records loaded:", nrow(cnv_all), "\n\n")

# ===========================
# Determine mCA status (using only date-filtered patients)
# ===========================
cat("Step 3: Determining mCA status from date-filtered cohort...\n")

# Filter CNV to only patients who passed date filtering
cnv_filtered_for_status <- cnv_all %>%
  filter(ID %in% mca_info$ID)

# Summarize CNV status for each patient
cnv_summary_status <- cnv_filtered_for_status %>%
  group_by(ID) %>%
  summarise(
    has_known_cnv = any(COPY_CHANGE %in% c("gain", "loss", "neutral")),
    has_only_unknown = all(COPY_CHANGE == "unknown"),
    .groups = "drop"
  )

# Get all date-filtered IDs
date_filtered_ids <- mca_info$ID

# Determine mCA status
cnv_status <- data.frame(ID = date_filtered_ids, stringsAsFactors = FALSE) %>%
  left_join(cnv_summary_status, by = "ID") %>%
  mutate(
    has_known_cnv = ifelse(is.na(has_known_cnv), FALSE, has_known_cnv),
    has_only_unknown = ifelse(is.na(has_only_unknown), FALSE, has_only_unknown),
    mCA_status = case_when(
      has_known_cnv ~ "mCA",
      !has_known_cnv & !has_only_unknown ~ "no_mCA",
      has_only_unknown ~ "unknown_only"
    )
  ) %>%
  filter(mCA_status %in% c("mCA", "no_mCA")) %>%  # Exclude patients with only "unknown"
  select(ID, mCA_status)

cat("  mCA positive:", sum(cnv_status$mCA_status == "mCA"), "\n")
cat("  mCA negative:", sum(cnv_status$mCA_status == "no_mCA"), "\n")
cat("  Total with definitive mCA results:", nrow(cnv_status), "\n\n")

# Merge mCA status with mca_info
mca_info <- mca_info %>%
  inner_join(cnv_status, by = "ID")

cat("  Final cohort after excluding 'unknown' CNV patients:", nrow(mca_info), "\n")
cat("  Expected: 452,594\n")
if (nrow(mca_info) != 452594) {
  warning("⚠️  Participant count differs from expected 452,594! Got: ", nrow(mca_info))
}
cat("\n")

# ===========================
# Define mCA subtypes (from filtered cohort)
# ===========================
cat("Step 4: Defining mCA subtypes...\n")

# Filter CNV data to only filtered cohort
cnv_filtered <- cnv_all %>%
  filter(ID %in% mca_info$ID)

# Patients with only 'unknown' CNVs (should be none after filtering)
unknown_only_ids <- cnv_filtered %>%
  group_by(ID) %>%
  summarise(only_unknown = all(COPY_CHANGE == "unknown")) %>%
  filter(only_unknown) %>%
  pull(ID)

# complex_mCA: >= 3 unique chromosomes
complex_mca <- cnv_filtered %>%
  filter(COPY_CHANGE %in% c("gain", "loss", "neutral"), !ID %in% unknown_only_ids) %>%
  distinct(ID, CHR) %>%
  group_by(ID) %>%
  summarise(n_chr = n_distinct(CHR)) %>%
  filter(n_chr >= 3) %>%
  pull(ID)

cat("  complex_mCA:", length(complex_mca), "\n")

# low_vaf: all CELL_FRAC < 0.1
low_vaf <- cnv_filtered %>%
  filter(COPY_CHANGE %in% c("gain", "loss", "neutral")) %>%
  group_by(ID) %>%
  filter(all(CELL_FRAC < 0.1)) %>%
  distinct(ID) %>%
  pull(ID)

cat("  low_vaf:", length(low_vaf), "\n")

# high_vaf: any CELL_FRAC >= 0.1
high_vaf <- cnv_filtered %>%
  filter(COPY_CHANGE %in% c("gain", "loss", "neutral")) %>%
  group_by(ID) %>%
  filter(any(CELL_FRAC >= 0.1)) %>%
  distinct(ID) %>%
  pull(ID)

cat("  high_vaf:", length(high_vaf), "\n")

# gain_only: has gain, and no loss or neutral
gain_only <- cnv_filtered %>%
  filter(COPY_CHANGE %in% c("gain", "loss", "neutral")) %>%
  group_by(ID) %>%
  summarise(type_set = list(unique(COPY_CHANGE)), .groups = "drop") %>%
  rowwise() %>%
  filter("gain" %in% type_set && !any(type_set %in% c("loss", "neutral"))) %>%
  ungroup() %>%
  pull(ID)

cat("  gain_only:", length(gain_only), "\n")

# loss_only: has loss, and no gain or neutral
loss_only <- cnv_filtered %>%
  filter(COPY_CHANGE %in% c("gain", "loss", "neutral")) %>%
  group_by(ID) %>%
  summarise(type_set = list(unique(COPY_CHANGE)), .groups = "drop") %>%
  rowwise() %>%
  filter("loss" %in% type_set && !any(type_set %in% c("gain", "neutral"))) %>%
  ungroup() %>%
  pull(ID)

cat("  loss_only:", length(loss_only), "\n")

# neutral_only: has neutral, and no gain or loss
neutral_only <- cnv_filtered %>%
  filter(COPY_CHANGE %in% c("gain", "loss", "neutral")) %>%
  group_by(ID) %>%
  summarise(type_set = list(unique(COPY_CHANGE)), .groups = "drop") %>%
  rowwise() %>%
  filter("neutral" %in% type_set && !any(type_set %in% c("gain", "loss"))) %>%
  ungroup() %>%
  pull(ID)

cat("  neutral_only:", length(neutral_only), "\n\n")

# ===========================
# Analyze mCA composition
# ===========================
cat("Step 4b: Analyzing mCA composition (total mCAs and types)...\n")

# Count total mCAs and their types from filtered cohort
mca_composition <- cnv_filtered %>%
  filter(COPY_CHANGE %in% c("gain", "loss", "neutral")) %>%
  group_by(COPY_CHANGE) %>%
  summarise(n_mCAs = n(), .groups = "drop") %>%
  mutate(
    percentage = round(100 * n_mCAs / sum(n_mCAs), 1),
    type_label = case_when(
      COPY_CHANGE == "gain" ~ "Gain ('+')",
      COPY_CHANGE == "loss" ~ "Loss ('-')",
      COPY_CHANGE == "neutral" ~ "CNLOH ('=')"
    )
  )

total_mCAs <- sum(mca_composition$n_mCAs)

cat("  Total mCAs identified:", total_mCAs, "\n")
cat("  Composition:\n")
for (i in 1:nrow(mca_composition)) {
  cat("    ", mca_composition$type_label[i], ": ", 
      mca_composition$n_mCAs[i], " (", 
      mca_composition$percentage[i], "%)\n", sep = "")
}
cat("\n")

# ===========================
# Load death data
# ===========================
cat("Step 5: Loading death data...\n")
death_info <- read_excel("~/Desktop/UK BB/Death date.xlsx") %>%
  rename(ID = Participant_ID, death_date = date_of_death) %>%
  mutate(ID = as.character(ID),
         death_date = as.Date(death_date)) %>%
  distinct(ID, .keep_all = TRUE)

mca_info <- mca_info %>%
  left_join(death_info, by = "ID")

cat("  Deaths recorded:", sum(!is.na(mca_info$death_date)), "\n\n")

# ===========================
# Load diagnosis data for baseline disease burden
# ===========================
cat("Step 6: Loading diagnosis data and calculating baseline disease burden...\n")
diag_data <- fread("~/Desktop/UK BB/Diagnosis code date.tsv", sep = "\t") %>%
  mutate(patient_id = as.character(patient_id),
         admission_date = as.Date(admission_date),
         diagnosis = as.character(diagnosis)) %>%
  filter(!is.na(admission_date), !is.na(diagnosis))

# Calculate baseline disease burden
baseline_disease_burden <- diag_data %>%
  inner_join(mca_info %>% select(ID, mCA_date), 
             by = c("patient_id" = "ID"),
             relationship = "many-to-one") %>%
  filter(admission_date <= mCA_date) %>%
  group_by(patient_id) %>%
  summarise(
    n_baseline_diseases = n_distinct(diagnosis),
    .groups = "drop"
  )

mca_info <- mca_info %>%
  left_join(baseline_disease_burden, by = c("ID" = "patient_id")) %>%
  mutate(n_baseline_diseases = replace_na(n_baseline_diseases, 0))

cat("  Baseline disease burden calculated\n\n")

# ===========================
# Create survival dataset
# ===========================
cat("Step 7: Creating survival dataset...\n")
censor_date <- as.Date("2023-07-14")

survival_df <- mca_info %>%
  mutate(
    end_date = if_else(!is.na(death_date), death_date, censor_date),
    time_years = as.numeric(difftime(end_date, mCA_date, units = "days")) / 365.25
  ) %>%
  filter(time_years >= 0) %>%  # Exclude deaths before mCA assessment
  select(ID, mCA_status, age, sex, n_baseline_diseases, mCA_date, death_date, time_years)

cat("  Final analytical cohort:", nrow(survival_df), "participants\n")
cat("  Expected: 452,594\n")
if (abs(nrow(survival_df) - 452594) > 100) {
  warning("⚠️  Participant count differs from expected 452,594!")
}
cat("  mCA participants:", sum(survival_df$mCA_status == "mCA"), "\n")
cat("  No mCA participants:", sum(survival_df$mCA_status == "no_mCA"), "\n")
cat("================================================================================\n\n")

# =========================================================================================
# QUALITY CHECKS BEFORE ANALYSIS
# =========================================================================================
cat("================================================================================\n")
cat("  QUALITY CHECKS - VERIFYING COHORT COMPOSITION\n")
cat("================================================================================\n\n")

cat("📊 Overall Cohort Summary:\n")
cat("  Total participants:", nrow(survival_df), "\n")
cat("  Expected:", 452594, "\n")
cat("  Difference:", nrow(survival_df) - 452594, "\n\n")

cat("📊 mCA Status Distribution:\n")
mca_status_table <- table(survival_df$mCA_status)
print(mca_status_table)
cat("  Percent mCA:", round(100 * mca_status_table["mCA"] / nrow(survival_df), 2), "%\n\n")

cat("📊 Mortality Summary:\n")
cat("  Total deaths:", sum(!is.na(survival_df$death_date)), "\n")
cat("  Deaths in mCA:", sum(!is.na(survival_df$death_date) & survival_df$mCA_status == "mCA"), "\n")
cat("  Deaths in no_mCA:", sum(!is.na(survival_df$death_date) & survival_df$mCA_status == "no_mCA"), "\n")
cat("  Overall mortality rate:", round(100 * sum(!is.na(survival_df$death_date)) / nrow(survival_df), 2), "%\n\n")

cat("📊 Follow-up Time Summary:\n")
cat("  Mean follow-up (years):", round(mean(survival_df$time_years), 2), "\n")
cat("  Median follow-up (years):", round(median(survival_df$time_years), 2), "\n")
cat("  Min follow-up (years):", round(min(survival_df$time_years), 2), "\n")
cat("  Max follow-up (years):", round(max(survival_df$time_years), 2), "\n\n")

cat("📊 Age Summary:\n")
cat("  Mean age:", round(mean(survival_df$age, na.rm = TRUE), 2), "\n")
cat("  Median age:", round(median(survival_df$age, na.rm = TRUE), 2), "\n")
age_q1 <- round(quantile(survival_df$age, 0.25, na.rm = TRUE), 0)
age_q3 <- round(quantile(survival_df$age, 0.75, na.rm = TRUE), 0)
cat("  IQR:", age_q1, "-", age_q3, "\n")
cat("  Age range:", round(min(survival_df$age, na.rm = TRUE), 2), "-", 
    round(max(survival_df$age, na.rm = TRUE), 2), "\n\n")

cat("📊 Sex Distribution:\n")
print(table(survival_df$sex))
cat("\n")

cat("📊 Baseline Disease Burden Summary:\n")
cat("  Mean:", round(mean(survival_df$n_baseline_diseases, na.rm = TRUE), 2), "\n")
cat("  Median:", round(median(survival_df$n_baseline_diseases, na.rm = TRUE), 2), "\n")
cat("  Range:", min(survival_df$n_baseline_diseases, na.rm = TRUE), "-", 
    max(survival_df$n_baseline_diseases, na.rm = TRUE), "\n\n")

cat("📊 mCA Composition Summary:\n")
cat("  Total mCAs identified:", total_mCAs, "\n")
cat("  Composition breakdown:\n")
for (i in 1:nrow(mca_composition)) {
  cat("    ", mca_composition$type_label[i], ": ", 
      mca_composition$n_mCAs[i], " (", 
      mca_composition$percentage[i], "%)\n", sep = "")
}
cat("\n")

cat("================================================================================\n\n")

# =========================================================================================
# MORTALITY ANALYSIS: mCA SUBTYPES vs NO mCA (OVERALL MORTALITY)
# =========================================================================================
output_dir <- "~/Desktop/UK BB/Mortality_Subtype_Analysis"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("================================================================================\n")
cat("  MORTALITY ANALYSIS: mCA SUBTYPES vs NO mCA (OVERALL MORTALITY)\n")
cat("================================================================================\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Output directory:", output_dir, "\n")
cat("Participants:", nrow(survival_df), "\n")
cat("Deaths:", sum(!is.na(survival_df$death_date)), "\n")
cat("Censored:", sum(is.na(survival_df$death_date)), "\n")
cat("================================================================================\n\n")

# =========================
# Define mCA subtypes for analysis
# =========================
cat("⚙️  Setting up mCA subtypes for analysis...\n")

# Get control IDs (no_mCA)
control_ids <- survival_df %>%
  filter(mCA_status == "no_mCA") %>%
  pull(ID)

# Get all mCA IDs
all_mca_ids <- survival_df %>%
  filter(mCA_status == "mCA") %>%
  pull(ID)

# Define subtypes with names
mca_subtypes <- list(
  list(name = "mCA", label = "All mCA", ids = all_mca_ids),
  list(name = "complex_mCA", label = "Complex mCA (≥3 chr)", ids = complex_mca),
  list(name = "high_vaf", label = "High VAF (≥10%)", ids = high_vaf),
  list(name = "low_vaf", label = "Low VAF (<10%)", ids = low_vaf),
  list(name = "gain_only", label = "Gain only", ids = gain_only),
  list(name = "loss_only", label = "Loss only", ids = loss_only),
  list(name = "neutral_only", label = "Neutral only", ids = neutral_only)
)

cat("✅ mCA subtypes defined:\n")
for (subtype in mca_subtypes) {
  cat("   ", subtype$label, ":", length(subtype$ids), "patients\n")
}
cat("   Controls (no_mCA):", length(control_ids), "patients\n\n")

# =========================
# QUALITY CHECK: Verify subtype overlaps and coverage
# =========================
cat("================================================================================\n")
cat("  QUALITY CHECK: mCA SUBTYPE VERIFICATION\n")
cat("================================================================================\n\n")

cat("📊 Subtype Overlaps Check:\n")
cat("  NOTE: Subtypes can overlap (e.g., a patient can be both complex and high VAF)\n\n")

# Check that all subtypes are subsets of all mCA
for (i in 2:length(mca_subtypes)) {
  subtype <- mca_subtypes[[i]]
  not_in_mca <- sum(!subtype$ids %in% all_mca_ids)
  if (not_in_mca > 0) {
    warning("⚠️  ", subtype$name, " has ", not_in_mca, " patients NOT in all mCA!")
  } else {
    cat("  ✅ ", subtype$name, ": all patients are in mCA cohort\n")
  }
}

cat("\n📊 Subtype Coverage:\n")
cat("  Total mCA patients:", length(all_mca_ids), "\n")
cat("  In complex_mCA:", length(complex_mca), "(", 
    round(100 * length(complex_mca) / length(all_mca_ids), 1), "%)\n")
cat("  In high_vaf:", length(high_vaf), "(", 
    round(100 * length(high_vaf) / length(all_mca_ids), 1), "%)\n")
cat("  In low_vaf:", length(low_vaf), "(", 
    round(100 * length(low_vaf) / length(all_mca_ids), 1), "%)\n")
cat("  In gain_only:", length(gain_only), "(", 
    round(100 * length(gain_only) / length(all_mca_ids), 1), "%)\n")
cat("  In loss_only:", length(loss_only), "(", 
    round(100 * length(loss_only) / length(all_mca_ids), 1), "%)\n")
cat("  In neutral_only:", length(neutral_only), "(", 
    round(100 * length(neutral_only) / length(all_mca_ids), 1), "%)\n")

# Check VAF coverage
vaf_overlap <- length(intersect(high_vaf, low_vaf))
vaf_coverage <- length(union(high_vaf, low_vaf))
cat("\n📊 VAF Category Check:\n")
cat("  High VAF + Low VAF overlap:", vaf_overlap, 
    ifelse(vaf_overlap > 0, " ⚠️  SHOULD BE 0!", " ✅"), "\n")
cat("  High VAF + Low VAF coverage:", vaf_coverage, "of", length(all_mca_ids), 
    "(", round(100 * vaf_coverage / length(all_mca_ids), 1), "%)\n")

# Check copy change exclusivity
copy_change_ids <- unique(c(gain_only, loss_only, neutral_only))
cat("\n📊 Copy Change Category Check:\n")
cat("  Patients in gain/loss/neutral_only:", length(copy_change_ids), "of", length(all_mca_ids),
    "(", round(100 * length(copy_change_ids) / length(all_mca_ids), 1), "%)\n")
cat("  Gain ∩ Loss overlap:", length(intersect(gain_only, loss_only)), 
    ifelse(length(intersect(gain_only, loss_only)) > 0, " ⚠️  SHOULD BE 0!", " ✅"), "\n")
cat("  Gain ∩ Neutral overlap:", length(intersect(gain_only, neutral_only)), 
    ifelse(length(intersect(gain_only, neutral_only)) > 0, " ⚠️  SHOULD BE 0!", " ✅"), "\n")
cat("  Loss ∩ Neutral overlap:", length(intersect(loss_only, neutral_only)), 
    ifelse(length(intersect(loss_only, neutral_only)) > 0, " ⚠️  SHOULD BE 0!", " ✅"), "\n")

cat("\n📊 Mortality Rates by Subtype (Unadjusted):\n")
for (subtype in mca_subtypes) {
  subtype_deaths <- sum(survival_df$ID %in% subtype$ids & !is.na(survival_df$death_date))
  mortality_rate <- round(100 * subtype_deaths / length(subtype$ids), 2)
  cat("  ", subtype$label, ": ", subtype_deaths, "/", length(subtype$ids), 
      " (", mortality_rate, "%)\n", sep = "")
}
control_deaths <- sum(survival_df$ID %in% control_ids & !is.na(survival_df$death_date))
control_mortality_rate <- round(100 * control_deaths / length(control_ids), 2)
cat("  no_mCA: ", control_deaths, "/", length(control_ids), 
    " (", control_mortality_rate, "%)\n", sep = "")

cat("\n================================================================================\n\n")

# =========================
# Initialize results
# =========================
all_results <- list()
results_file <- file.path(output_dir, "mortality_subtype_analysis_INCREMENTAL.csv")

# =========================
# Sequential analysis loop
# =========================
cat("================================================================================\n")
cat("  STARTING SEQUENTIAL ANALYSIS\n")
cat("================================================================================\n\n")

start_time <- Sys.time()

for (i in seq_along(mca_subtypes)) {
  subtype <- mca_subtypes[[i]]
  
  cat("--------------------------------------------------------------------------------\n")
  cat("[", i, "/", length(mca_subtypes), "] Processing:", subtype$name, "-", subtype$label, "\n")
  cat("--------------------------------------------------------------------------------\n")
  
  # Step 1: Create analysis dataset (subtype vs no_mCA)
  cat("  Step 1/5: Creating analysis dataset (", subtype$name, " vs no_mCA)...\n")
  
  # Combine subtype cases with controls
  analysis_ids <- c(subtype$ids, control_ids)
  
  sdf <- survival_df %>%
    filter(ID %in% analysis_ids) %>%
    mutate(
      group = ifelse(ID %in% subtype$ids, subtype$name, "no_mCA"),
      event = ifelse(!is.na(death_date), 1, 0),
      status = ifelse(event == 1, 1, 0)  # For overall mortality: 1 = death, 0 = censored
    )
  
  n_subtype <- sum(sdf$group == subtype$name)
  n_control <- sum(sdf$group == "no_mCA")
  n_events_subtype <- sum(sdf$event == 1 & sdf$group == subtype$name)
  n_events_control <- sum(sdf$event == 1 & sdf$group == "no_mCA")
  n_total_events <- sum(sdf$event == 1)
  
  cat("            ", subtype$name, ":", n_subtype, "(deaths:", n_events_subtype, ")\n")
  cat("            no_mCA:", n_control, "(deaths:", n_events_control, ")\n")
  cat("            Total events:", n_total_events, "\n")
  
  # Quality check: verify numbers
  cat("            QUALITY CHECK:\n")
  cat("              Expected subtype N:", length(subtype$ids), 
      ifelse(n_subtype == length(subtype$ids), " ✅", " ⚠️ MISMATCH!"), "\n")
  cat("              Expected control N:", length(control_ids), 
      ifelse(n_control == length(control_ids), " ✅", " ⚠️ MISMATCH!"), "\n")
  cat("              Mortality rate (subtype):", round(100 * n_events_subtype / n_subtype, 2), "%\n")
  cat("              Mortality rate (control):", round(100 * n_events_control / n_control, 2), "%\n")
  cat("              Mean follow-up (subtype):", round(mean(sdf$time_years[sdf$group == subtype$name]), 2), "years\n")
  cat("              Mean follow-up (control):", round(mean(sdf$time_years[sdf$group == "no_mCA"]), 2), "years\n")
  
  # Initialize result
  result <- data.frame(
    Subtype = subtype$name,
    Label = subtype$label,
    N_Subtype = n_subtype,
    N_Control = n_control,
    Events_Subtype = n_events_subtype,
    Events_Control = n_events_control,
    N_Events = n_total_events,
    HR_Unadjusted = NA,
    HR_Unadj_Lower = NA,
    HR_Unadj_Upper = NA,
    HR_Unadj_P = NA,
    HR_Adjusted = NA,
    HR_Adj_Lower = NA,
    HR_Adj_Upper = NA,
    HR_Adj_P = NA,
    Status = "insufficient_events",
    stringsAsFactors = FALSE
  )
  
  if (n_total_events >= 10 && n_events_subtype > 0 && n_events_control > 0) {
    
    # Step 2: Unadjusted Cox model
    cat("  Step 2/5: Computing unadjusted HR (Cox)...\n")
    
    tryCatch({
      library(survival)
      sdf$group_numeric <- as.numeric(sdf$group == subtype$name)
      cox_unadj <- coxph(Surv(time_years, event) ~ group_numeric, data = sdf)
      
      coef_summary <- summary(cox_unadj)$coefficients
      conf_int <- summary(cox_unadj)$conf.int
      
      result$HR_Unadjusted <- exp(coef_summary[1, "coef"])
      result$HR_Unadj_Lower <- conf_int[1, "lower .95"]
      result$HR_Unadj_Upper <- conf_int[1, "upper .95"]
      result$HR_Unadj_P <- coef_summary[1, "Pr(>|z|)"]
      
      cat("            Unadjusted HR:", round(result$HR_Unadjusted, 3), 
          "95% CI: [", round(result$HR_Unadj_Lower, 3), "-", round(result$HR_Unadj_Upper, 3), "]",
          "p =", format.pval(result$HR_Unadj_P, digits = 3), "\n")
    }, error = function(e) {
      cat("            ⚠️  Unadjusted HR calculation failed:", e$message, "\n")
    })
    
    # Step 3: Adjusted Cox model
    cat("  Step 3/5: Computing adjusted HR (age, sex, baseline disease burden)...\n")
    
    tryCatch({
      sdf$sex_numeric <- as.numeric(as.factor(sdf$sex)) - 1
      cox_adj <- coxph(Surv(time_years, event) ~ group_numeric + age + sex_numeric + n_baseline_diseases, 
                       data = sdf)
      
      coef_summary <- summary(cox_adj)$coefficients
      conf_int <- summary(cox_adj)$conf.int
      
      result$HR_Adjusted <- exp(coef_summary[1, "coef"])
      result$HR_Adj_Lower <- conf_int[1, "lower .95"]
      result$HR_Adj_Upper <- conf_int[1, "upper .95"]
      result$HR_Adj_P <- coef_summary[1, "Pr(>|z|)"]
      
      cat("            Adjusted HR:  ", round(result$HR_Adjusted, 3), 
          "95% CI: [", round(result$HR_Adj_Lower, 3), "-", round(result$HR_Adj_Upper, 3), "]",
          "p =", format.pval(result$HR_Adj_P, digits = 3), "\n")
    }, error = function(e) {
      cat("            ⚠️  Adjusted HR calculation failed:", e$message, "\n")
    })
    
    result$Status <- "success"
    cat("  ✅ Analysis complete\n")
  } else {
    cat("  ⚠️  Skipped (insufficient events)\n")
  }
  
  # Save incrementally
  all_results[[i]] <- result
  results_temp <- do.call(rbind, all_results)
  write.csv(results_temp, results_file, row.names = FALSE)
  
  cat("  💾 Saved to:", results_file, "\n")
  
  # Time estimate
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  avg_per_subtype <- elapsed / i
  remaining <- (length(mca_subtypes) - i) * avg_per_subtype
  cat("  ⏱️  Progress:", round(i/length(mca_subtypes)*100, 1), "% |",
      "Elapsed:", round(elapsed, 1), "min |",
      "Est. remaining:", round(remaining, 1), "min\n")
  cat("\n")
}

# =========================
# Final processing
# =========================
cat("================================================================================\n")
cat("  FINALIZING RESULTS\n")
cat("================================================================================\n\n")

cat("📊 Applying FDR correction...\n")

results_df <- do.call(rbind, all_results)

valid_unadj <- !is.na(results_df$HR_Unadj_P)
valid_adj <- !is.na(results_df$HR_Adj_P)

results_df$HR_Unadj_Q <- NA
results_df$HR_Adj_Q <- NA

if (sum(valid_unadj) > 0) {
  results_df$HR_Unadj_Q[valid_unadj] <- p.adjust(results_df$HR_Unadj_P[valid_unadj], method = "BH")
  cat("   Unadjusted FDR: applied to", sum(valid_unadj), "tests\n")
}

if (sum(valid_adj) > 0) {
  results_df$HR_Adj_Q[valid_adj] <- p.adjust(results_df$HR_Adj_P[valid_adj], method = "BH")
  cat("   Adjusted FDR: applied to", sum(valid_adj), "tests\n")
}

cat("\n📋 Formatting final table...\n")

results_final <- results_df %>%
  mutate(
    HR_Unadj_CI = ifelse(!is.na(HR_Unadjusted),
                         paste0(sprintf("%.3f", HR_Unadjusted), " (",
                                sprintf("%.3f", HR_Unadj_Lower), "-",
                                sprintf("%.3f", HR_Unadj_Upper), ")"), NA),
    HR_Adj_CI = ifelse(!is.na(HR_Adjusted),
                       paste0(sprintf("%.3f", HR_Adjusted), " (",
                              sprintf("%.3f", HR_Adj_Lower), "-",
                              sprintf("%.3f", HR_Adj_Upper), ")"), NA)
  ) %>%
  select(Subtype, Label, N_Subtype, N_Control, Events_Subtype, Events_Control, N_Events,
         HR_Unadjusted, HR_Unadj_Lower, HR_Unadj_Upper, HR_Unadj_P, HR_Unadj_Q, HR_Unadj_CI,
         HR_Adjusted, HR_Adj_Lower, HR_Adj_Upper, HR_Adj_P, HR_Adj_Q, HR_Adj_CI,
         Status)

final_file <- file.path(output_dir, "mortality_subtype_analysis_FINAL.csv")
write.csv(results_final, final_file, row.names = FALSE)

cat("\n")
cat("================================================================================\n")
cat("  ANALYSIS COMPLETE!\n")
cat("================================================================================\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Total time:", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1), "minutes\n")
cat("Final results saved to:", final_file, "\n")
cat("================================================================================\n\n")

# =========================
# Final Quality Check Summary
# =========================
cat("================================================================================\n")
cat("  FINAL QUALITY CHECK SUMMARY\n")
cat("================================================================================\n\n")

cat("📊 Analysis Completion Status:\n")
success_count <- sum(results_final$Status == "success")
insufficient_count <- sum(results_final$Status == "insufficient_events")
cat("  Successful analyses:", success_count, "/", nrow(results_final), "\n")
cat("  Insufficient events:", insufficient_count, "/", nrow(results_final), "\n\n")

cat("📊 Statistical Significance Summary (p < 0.05):\n")
cat("  Unadjusted HR significant:", sum(results_final$HR_Unadj_P < 0.05, na.rm = TRUE), "/", 
    sum(!is.na(results_final$HR_Unadj_P)), "\n")
cat("  Adjusted HR significant:", sum(results_final$HR_Adj_P < 0.05, na.rm = TRUE), "/", 
    sum(!is.na(results_final$HR_Adj_P)), "\n\n")

cat("📊 FDR Significance Summary (q < 0.05):\n")
cat("  Unadjusted HR FDR significant:", sum(results_final$HR_Unadj_Q < 0.05, na.rm = TRUE), "/", 
    sum(!is.na(results_final$HR_Unadj_Q)), "\n")
cat("  Adjusted HR FDR significant:", sum(results_final$HR_Adj_Q < 0.05, na.rm = TRUE), "/", 
    sum(!is.na(results_final$HR_Adj_Q)), "\n\n")

cat("📊 Hazard Ratio Range:\n")
cat("  Unadjusted HR range:", round(min(results_final$HR_Unadjusted, na.rm = TRUE), 3), "-",
    round(max(results_final$HR_Unadjusted, na.rm = TRUE), 3), "\n")
cat("  Adjusted HR range:", round(min(results_final$HR_Adjusted, na.rm = TRUE), 3), "-",
    round(max(results_final$HR_Adjusted, na.rm = TRUE), 3), "\n\n")

cat("📊 Top Results by Adjusted HR (highest to lowest):\n")
results_ordered <- results_final %>%
  filter(!is.na(HR_Adjusted)) %>%
  arrange(desc(HR_Adjusted)) %>%
  select(Subtype, Label, HR_Adjusted, HR_Adj_P, HR_Adj_Q)
print(results_ordered, row.names = FALSE, digits = 3)

cat("\n================================================================================\n\n")

cat("📊 FINAL RESULTS TABLE:\n\n")
print(results_final, width = 200)

cat("\n✅ All done! Check the output directory for your results.\n")
