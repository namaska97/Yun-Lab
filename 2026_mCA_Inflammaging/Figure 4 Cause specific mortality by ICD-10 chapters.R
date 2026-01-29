# =========================================================================================
# MORTALITY ANALYSIS WITH PROPER DATA FILTERING
# =========================================================================================
# This version includes the data filtering steps from the methods section
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
# Load CNV data and determine mCA status
# ===========================
cat("Step 2: Loading CNV data and determining mCA status...\n")
cnv <- read_excel("~/Desktop/UK BB/CNV_data_19808.xlsx") %>%
  mutate(ID = as.character(ID)) %>%
  filter(ID %in% analyzed$ID)

# Summarize CNV status for each patient
cnv_summary <- cnv %>%
  group_by(ID) %>%
  summarise(
    has_known_cnv = any(COPY_CHANGE %in% c("gain", "loss", "neutral")),
    has_only_unknown = all(COPY_CHANGE == "unknown"),
    .groups = "drop"
  )

# Determine mCA status
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
  filter(mCA_status %in% c("mCA", "no_mCA")) %>%  # Exclude patients with only "unknown"
  select(ID, mCA_status)

cat("  mCA positive:", sum(cnv_status$mCA_status == "mCA"), "\n")
cat("  mCA negative:", sum(cnv_status$mCA_status == "no_mCA"), "\n")
cat("  Total with definitive mCA results:", nrow(cnv_status), "\n\n")

# Merge mCA status with mca_info
mca_info <- mca_info %>%
  inner_join(cnv_status, by = "ID") %>%
  rename(group = mCA_status) %>%
  mutate(group = ifelse(group == "mCA", "mCA", "No mCA"))

cat("  Final cohort after excluding 'unknown' CNV patients:", nrow(mca_info), "\n\n")

# ===========================
# Load death data
# ===========================
cat("Step 3: Loading death data...\n")
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
cat("Step 4: Loading diagnosis data and calculating baseline disease burden...\n")
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
cat("Step 5: Creating survival dataset...\n")
censor_date <- as.Date("2023-07-14")

survival_df <- mca_info %>%
  mutate(
    end_date = if_else(!is.na(death_date), death_date, censor_date),
    time_years = as.numeric(difftime(end_date, mCA_date, units = "days")) / 365.25
  ) %>%
  filter(time_years >= 0) %>%  # Exclude deaths before mCA assessment
  select(ID, group, age, sex, n_baseline_diseases, mCA_date, death_date, time_years)

cat("  Final analytical cohort:", nrow(survival_df), "participants\n")
cat("  Expected: 452,594\n")
if (abs(nrow(survival_df) - 452594) > 100) {
  warning("⚠️  Participant count differs from expected 452,594!")
}
cat("  mCA participants:", sum(survival_df$group == "mCA"), "\n")
cat("  No mCA participants:", sum(survival_df$group == "No mCA"), "\n")
cat("================================================================================\n\n")

# ===========================
# Load death cause data (ICD-10 codes)
# ===========================
cat("Step 6: Loading death cause data...\n")
death_cause <- read_excel("~/Desktop/UK BB/Death cause.xlsx") %>%
  rename(ID = Participant_ID, icd10_code = `ICD-10 code`) %>%
  mutate(ID = as.character(ID),
         icd10_code = as.character(icd10_code)) %>%
  filter(!is.na(icd10_code))

cat("  Death cause records loaded:", nrow(death_cause), "\n")
cat("  Unique participants with death causes:", n_distinct(death_cause$ID), "\n\n")

# =========================================================================================
# MORTALITY ANALYSIS
# =========================================================================================
output_dir <- "~/Desktop/UK BB/Mortality_CIF_Plots_with_HR"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("================================================================================\n")
cat("  MORTALITY ANALYSIS - 13 CATEGORIES (SEQUENTIAL WITH PROGRESS TRACKING)\n")
cat("================================================================================\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Output directory:", output_dir, "\n")
cat("Participants:", nrow(survival_df), "\n")
cat("Deaths:", sum(!is.na(survival_df$death_date)), "\n")
cat("Censored:", sum(is.na(survival_df$death_date)), "\n")
cat("================================================================================\n\n")

# =========================
# Initialize
# =========================
all_results <- list()
results_file <- file.path(output_dir, "mortality_analysis_with_HR_INCREMENTAL.csv")

# =========================
# Define 13 categories
# =========================
cat("⚙️  Setting up analysis categories...\n")

all_death_chapters <- death_cause %>%
  mutate(chapter = substr(icd10_code, 1, 1)) %>%
  distinct(chapter) %>%
  arrange(chapter) %>%
  pull(chapter)

# Exclude A, B, C (in combined categories), D (split), and others
excluded_chapters <- c("", "A", "B", "C", "D", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z")
standard_chapters <- setdiff(all_death_chapters, excluded_chapters)

analysis_categories <- list(
  list(name = "AB", label = "A+B", type = "multi_chapter", chapters = c("A", "B")),
  list(name = "C_D00-D48", label = "C+D00-D48 (Neoplasms)", type = "neoplasm", chapters = NULL),
  list(name = "D50-D89", label = "D50-D89 (Blood)", type = "blood", chapters = NULL)
)

for (ch in standard_chapters) {
  analysis_categories[[length(analysis_categories) + 1]] <- list(
    name = ch, label = ch, type = "single_chapter", chapters = ch
  )
}

cat("✅ Categories defined:", length(analysis_categories), "total\n")
cat("   Combined: AB, C+D00-D48, D50-D89\n")
cat("   Individual:", paste(standard_chapters, collapse = ", "), "\n\n")

# =========================
# Sequential analysis loop
# =========================
cat("================================================================================\n")
cat("  STARTING SEQUENTIAL ANALYSIS\n")
cat("================================================================================\n\n")

start_time <- Sys.time()

for (i in seq_along(analysis_categories)) {
  category <- analysis_categories[[i]]
  
  cat("--------------------------------------------------------------------------------\n")
  cat("[", i, "/", length(analysis_categories), "] Processing:", category$name, "-", category$label, "\n")
  cat("--------------------------------------------------------------------------------\n")
  
  cat_name <- category$name
  cat_type <- category$type
  cat_chapters <- category$chapters
  
  # Step 1: Identify deaths
  cat("  Step 1/5: Identifying deaths in this category...\n")
  
  if (cat_type == "single_chapter") {
    deaths_in_category <- death_cause %>%
      filter(substr(icd10_code, 1, 1) == cat_chapters) %>%
      distinct(ID)
  } else if (cat_type == "multi_chapter") {
    deaths_in_category <- death_cause %>%
      filter(substr(icd10_code, 1, 1) %in% cat_chapters) %>%
      distinct(ID)
  } else if (cat_type == "neoplasm") {
    deaths_in_category <- death_cause %>%
      filter(substr(icd10_code, 1, 1) == "C" | 
               (substr(icd10_code, 1, 1) == "D" & 
                  substr(icd10_code, 2, 2) %in% c("0", "1", "2", "3", "4") &
                  as.numeric(substr(icd10_code, 2, 3)) <= 48)) %>%
      distinct(ID)
  } else if (cat_type == "blood") {
    deaths_in_category <- death_cause %>%
      filter(substr(icd10_code, 1, 1) == "D" & 
               substr(icd10_code, 2, 2) %in% c("5", "6", "7", "8") &
               as.numeric(substr(icd10_code, 2, 3)) >= 50 &
               as.numeric(substr(icd10_code, 2, 3)) <= 89) %>%
      distinct(ID)
  }
  
  # Step 2: Create dataset
  cat("  Step 2/5: Creating analysis dataset...\n")
  
  sdf <- survival_df %>%
    mutate(event = ifelse(ID %in% deaths_in_category$ID, 1, 0),
           status = case_when(!is.na(death_date) & event == 1 ~ 1,
                              !is.na(death_date) & event == 0 ~ 2,
                              TRUE ~ 0))
  
  n_events <- sum(sdf$status == 1)
  events_mca <- sum(sdf$status == 1 & sdf$group == "mCA")
  events_no_mca <- sum(sdf$status == 1 & sdf$group == "No mCA")
  n_mca <- sum(sdf$group == "mCA")
  n_no_mca <- sum(sdf$group == "No mCA")
  
  cat("            Total events:", n_events, "(mCA:", events_mca, ", No mCA:", events_no_mca, ")\n")
  
  # Initialize result
  result <- data.frame(Category = cat_name, Status = "insufficient_events",
                       N_mCA = n_mca, N_No_mCA = n_no_mca,
                       Events_mCA = events_mca, Events_No_mCA = events_no_mca,
                       N_Events = n_events, Gray_P = NA, HR_Unadjusted = NA,
                       HR_Unadj_Lower = NA, HR_Unadj_Upper = NA,
                       HR_Adjusted = NA, HR_Adj_Lower = NA, HR_Adj_Upper = NA,
                       CRR_Adj_P = NA, stringsAsFactors = FALSE)
  
  if (n_events >= 10 && n_distinct(sdf$group[sdf$status == 1]) == 2) {
    
    # Step 3: Gray's test
    cat("  Step 3/5: Running Gray's test...\n")
    
    gray_test <- tryCatch({
      cmprsk::cuminc(ftime = sdf$time_years, fstatus = sdf$status, group = sdf$group)
    }, error = function(e) NULL)
    
    if (!is.null(gray_test) && !is.null(gray_test$Tests) && nrow(gray_test$Tests) >= 1) {
      result$Gray_P <- gray_test$Tests[1, "pv"]
      cat("            Gray's p-value:", format.pval(result$Gray_P, digits = 3), "\n")
    }
    
    # Step 4: Unadjusted HR
    cat("  Step 4/5: Computing unadjusted HR...\n")
    
    tryCatch({
      sdf$group_numeric <- as.numeric(sdf$group == "mCA")
      crr_unadj <- cmprsk::crr(ftime = sdf$time_years, fstatus = sdf$status,
                               cov1 = matrix(sdf$group_numeric, ncol = 1),
                               failcode = 1, cencode = 0)
      if (!is.null(crr_unadj$coef) && !is.na(crr_unadj$coef[1])) {
        result$HR_Unadjusted <- exp(crr_unadj$coef[1])
        se <- sqrt(crr_unadj$var[1,1])
        result$HR_Unadj_Lower <- exp(crr_unadj$coef[1] - 1.96 * se)
        result$HR_Unadj_Upper <- exp(crr_unadj$coef[1] + 1.96 * se)
        cat("            Unadjusted HR:", round(result$HR_Unadjusted, 3), 
            "95% CI: [", round(result$HR_Unadj_Lower, 3), "-", round(result$HR_Unadj_Upper, 3), "]\n")
      }
    }, error = function(e) {
      cat("            ⚠️  Unadjusted HR calculation failed\n")
    })
    
    # Step 5: Adjusted HR
    cat("  Step 5/5: Computing adjusted HR (age, sex, baseline disease burden)...\n")
    
    tryCatch({
      sdf$sex_numeric <- as.numeric(as.factor(sdf$sex)) - 1
      cov_mat <- cbind(group = sdf$group_numeric, age = sdf$age,
                       sex = sdf$sex_numeric, n_baseline_diseases = sdf$n_baseline_diseases)
      crr_adj <- cmprsk::crr(ftime = sdf$time_years, fstatus = sdf$status,
                             cov1 = cov_mat, failcode = 1, cencode = 0)
      if (!is.null(crr_adj$coef) && !is.na(crr_adj$coef[1])) {
        result$HR_Adjusted <- exp(crr_adj$coef[1])
        z_score <- crr_adj$coef[1] / sqrt(crr_adj$var[1,1])
        result$CRR_Adj_P <- 2 * (1 - pnorm(abs(z_score)))
        se <- sqrt(crr_adj$var[1,1])
        result$HR_Adj_Lower <- exp(crr_adj$coef[1] - 1.96 * se)
        result$HR_Adj_Upper <- exp(crr_adj$coef[1] + 1.96 * se)
        cat("            Adjusted HR:  ", round(result$HR_Adjusted, 3), 
            "95% CI: [", round(result$HR_Adj_Lower, 3), "-", round(result$HR_Adj_Upper, 3), "]",
            "p =", format.pval(result$CRR_Adj_P, digits = 3), "\n")
      }
    }, error = function(e) {
      cat("            ⚠️  Adjusted HR calculation failed\n")
    })
    
    result$Status <- "success"
    cat("  ✅ Analysis complete\n")
  } else {
    cat("  ⚠️  Skipped (insufficient events or groups)\n")
  }
  
  # Save incrementally
  all_results[[i]] <- result
  results_temp <- do.call(rbind, all_results)
  write.csv(results_temp, results_file, row.names = FALSE)
  
  cat("  💾 Saved to:", results_file, "\n")
  
  # Time estimate
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  avg_per_category <- elapsed / i
  remaining <- (length(analysis_categories) - i) * avg_per_category
  cat("  ⏱️  Progress:", round(i/length(analysis_categories)*100, 1), "% |",
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

valid_gray <- !is.na(results_df$Gray_P)
valid_crr <- !is.na(results_df$CRR_Adj_P)

results_df$Gray_Q <- NA
results_df$CRR_Adj_Q <- NA

if (sum(valid_gray) > 0) {
  results_df$Gray_Q[valid_gray] <- p.adjust(results_df$Gray_P[valid_gray], method = "BH")
  cat("   Gray's test FDR: applied to", sum(valid_gray), "tests\n")
}

if (sum(valid_crr) > 0) {
  results_df$CRR_Adj_Q[valid_crr] <- p.adjust(results_df$CRR_Adj_P[valid_crr], method = "BH")
  cat("   CRR adjusted FDR: applied to", sum(valid_crr), "tests\n")
}

cat("\n📋 Formatting final table...\n")

results_final <- results_df %>%
  mutate(HR_Unadj_CI = ifelse(!is.na(HR_Unadjusted),
                              paste0(sprintf("%.3f", HR_Unadjusted), " (",
                                     sprintf("%.3f", HR_Unadj_Lower), "-",
                                     sprintf("%.3f", HR_Unadj_Upper), ")"), NA),
         HR_Adj_CI = ifelse(!is.na(HR_Adjusted),
                            paste0(sprintf("%.3f", HR_Adjusted), " (",
                                   sprintf("%.3f", HR_Adj_Lower), "-",
                                   sprintf("%.3f", HR_Adj_Upper), ")"), NA))

final_file <- file.path(output_dir, "mortality_analysis_with_HR_FINAL.csv")
write.csv(results_final, final_file, row.names = FALSE)

cat("\n")
cat("================================================================================\n")
cat("  ANALYSIS COMPLETE!\n")
cat("================================================================================\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Total time:", round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1), "minutes\n")
cat("Final results saved to:", final_file, "\n")
cat("================================================================================\n\n")

cat("📊 FINAL RESULTS TABLE:\n\n")
print(results_final)

cat("\n✅ All done! Check the output directory for your results.\n")
