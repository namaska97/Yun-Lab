# ===========================
# Clear workspace
# ===========================
rm(list = ls())

# ===========================
# Libraries
# ===========================
library(readxl)
library(data.table)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(cmprsk)
library(tidycmprsk)
library(survival)

# Try to load Excel writing package
excel_package <- NULL
if (requireNamespace("openxlsx", quietly = TRUE)) {
  library(openxlsx)
  excel_package <- "openxlsx"
} else if (requireNamespace("writexl", quietly = TRUE)) {
  library(writexl)
  excel_package <- "writexl"
} else {
  warning("Neither 'openxlsx' nor 'writexl' packages are available. Results will be saved as CSV files instead.")
  excel_package <- "csv"
}

# ===========================
# Config
# ===========================
censor_date <- as.Date("2023-07-14")
output_dir <- "~/Desktop/UK BB/KM_Curves_ICD10"
dir.create(output_dir, showWarnings = FALSE)

cat("\n=== STARTING ANALYSIS ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ===========================
# Load analyzed samples
# ===========================
cat("Loading analyzed samples...\n")
analyzed <- read_excel("~/Desktop/UK BB/ukb4777.samples_analyzed.xlsx") %>%
  mutate(ID = as.character(ID))

cat("Total analyzed samples:", nrow(analyzed), "\n\n")

# ===========================
# Load sample info with age and sex
# ===========================
cat("Loading participant data...\n")
mca_info_raw <- read_excel("~/Desktop/UK BB/Age sex dob collection date.xlsx") %>%
  mutate(ID = as.character(ID),
         mCA_date = mdy(`Date of mCA assessment`),
         age = as.numeric(`Age at blood collection`),
         sex = as.factor(Sex)) %>%
  filter(ID %in% analyzed$ID)  # Filter to analyzed samples only

cat("Participant records (analyzed samples only):", nrow(mca_info_raw), "\n")

# Filter out invalid dates
mca_info <- mca_info_raw %>%
  filter(!is.na(mCA_date),
         mCA_date != as.Date("1900-01-01"),
         mCA_date > as.Date("1901-01-01")) %>%
  select(ID, mCA_date, age, sex) %>%
  distinct(ID, .keep_all = TRUE)

cat("Valid participants after date filtering:", nrow(mca_info), "\n\n")

# ===========================
# Load CNV data and determine mCA status (SAME AS SUMMARY TABLE CODE)
# ===========================
cat("Loading CNV data and determining mCA status...\n")
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

cat("\n=== mCA Status Distribution ===\n")
print(table(cnv_status$mCA_status))
cat("Total patients with definitive mCA results:", nrow(cnv_status), "\n\n")

# Merge mCA status with mca_info (this filters to only patients with definitive mCA status)
mca_info <- mca_info %>%
  inner_join(cnv_status, by = "ID") %>%
  rename(mca_group = mCA_status) %>%
  mutate(mca_group = ifelse(mca_group == "mCA", "mCA", "No mCA"))

cat("Final cohort after excluding 'unknown' CNV patients:", nrow(mca_info), "\n")
cat("mCA participants:", sum(mca_info$mca_group == "mCA"), "\n")
cat("No mCA participants:", sum(mca_info$mca_group == "No mCA"), "\n\n")

# ===========================
# Load death data
# ===========================
cat("Loading death data...\n")
death_info <- read_excel("~/Desktop/UK BB/Death date.xlsx") %>%
  rename(ID = Participant_ID, death_date = date_of_death) %>%
  mutate(ID = as.character(ID),
         death_date = as.Date(death_date)) %>%
  distinct(ID, .keep_all = TRUE)

mca_info <- mca_info %>%
  left_join(death_info, by = "ID")

cat("Deaths recorded:", sum(!is.na(mca_info$death_date)), "\n\n")

# ===========================
# Load diagnosis data
# ===========================
cat("Loading diagnosis data...\n")
diag_data <- fread("~/Desktop/UK BB/Diagnosis code date.tsv", sep = "\t") %>%
  mutate(patient_id = as.character(patient_id),
         admission_date = as.Date(admission_date),
         diagnosis = as.character(diagnosis)) %>%
  filter(!is.na(admission_date), !is.na(diagnosis))

cat("Total diagnosis records:", nrow(diag_data), "\n")
cat("Unique patients in diagnosis data:", n_distinct(diag_data$patient_id), "\n\n")

# ===========================
# Calculate baseline disease burden (OPTIMIZED)
# ===========================
cat("=== Calculating baseline disease burden ===\n")

# Convert to data.table for speed
setDT(diag_data)
setDT(mca_info)

# Set keys for fast joining
setkey(diag_data, patient_id)
setkey(mca_info, ID)

# Filter diagnosis data to only include patients in mca_info FIRST
cat("Filtering diagnosis data to study participants...\n")
diag_data_filtered <- diag_data[patient_id %in% mca_info$ID]
cat("Filtered diagnosis records:", nrow(diag_data_filtered), "\n")

# Join with mCA dates
cat("Joining with mCA dates...\n")
diag_with_mca <- diag_data_filtered[mca_info[, .(ID, mCA_date)], 
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
mca_info <- as.data.frame(mca_info)

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
cat("Mean baseline diseases - mCA:", round(mean(mca_info$n_baseline_diseases[mca_info$mca_group == "mCA"]), 2), "\n")
cat("Mean baseline diseases - No mCA:", round(mean(mca_info$n_baseline_diseases[mca_info$mca_group == "No mCA"]), 2), "\n\n")

# ===========================
# Function to get incident diseases (OPTIMIZED)
# ===========================
get_incident_diseases <- function(diag_data, mca_info, filter_condition = NULL) {
  cat("  Getting incident diseases...")
  
  # Convert to data.table
  setDT(diag_data)
  dt_diag <- copy(diag_data)
  
  # Filter if condition provided
  if (!is.null(filter_condition)) {
    dt_diag <- dt_diag[eval(parse(text = filter_condition))]
    cat(" filtered to", format(nrow(dt_diag), big.mark = ","), "records...")
  }
  
  # Filter to study participants FIRST
  dt_diag <- dt_diag[patient_id %in% mca_info$ID]
  
  # Join with mCA dates
  mca_lookup <- data.table(ID = mca_info$ID, mCA_date = mca_info$mCA_date)
  setkey(mca_lookup, ID)
  setkey(dt_diag, patient_id)
  
  dt_diag_mca <- dt_diag[mca_lookup, on = .(patient_id = ID), nomatch = 0]
  
  # Pre-existing diseases
  pre_existing <- dt_diag_mca[admission_date <= mCA_date, 
                              .(patient_id, diagnosis)]
  pre_existing <- unique(pre_existing)
  
  # Post-mCA diseases
  post_mca <- dt_diag_mca[admission_date > mCA_date]
  
  # Exclude pre-existing
  setkey(post_mca, patient_id, diagnosis)
  setkey(pre_existing, patient_id, diagnosis)
  incident <- post_mca[!pre_existing]
  
  # Get first diagnosis date
  result <- incident[, .(first_dx_date = min(admission_date)), by = patient_id]
  
  cat(" found", format(nrow(result), big.mark = ","), "incident cases\n")
  return(as.data.frame(result))
}

# ===========================
# Competing Risk Analysis Function
# ===========================
perform_competing_risk_analysis <- function(df, disease_name) {
  df <- df %>%
    mutate(
      event_status = case_when(
        !is.na(first_dx_date) & first_dx_date > mCA_date & 
          (is.na(death_date) | first_dx_date <= death_date) ~ "disease",
        !is.na(death_date) & death_date > mCA_date ~ "death",
        TRUE ~ "censored"
      ),
      event_code = case_when(
        event_status == "disease" ~ 1,
        event_status == "death" ~ 2,
        TRUE ~ 0
      ),
      end_date = case_when(
        event_status == "disease" ~ first_dx_date,
        event_status == "death" ~ death_date,
        TRUE ~ censor_date
      ),
      time_years = as.numeric(difftime(end_date, mCA_date, units = "days")) / 365.25,
      group = factor(mca_group, levels = c("No mCA", "mCA"))
    ) %>%
    filter(time_years >= 0, !is.na(group), !is.na(age), !is.na(sex), !is.na(n_baseline_diseases))
  
  # Summary statistics
  n_mca <- sum(df$group == "mCA")
  n_no_mca <- sum(df$group == "No mCA")
  events_mca <- sum(df$group == "mCA" & df$event_code == 1)
  events_no_mca <- sum(df$group == "No mCA" & df$event_code == 1)
  
  # Initialize results
  hr_unadj <- hr_unadj_lower <- hr_unadj_upper <- gray_p <- NA
  hr_adj <- hr_adj_lower <- hr_adj_upper <- crr_p <- NA
  
  total_events <- sum(df$event_code == 1)
  groups_with_events <- n_distinct(df$group[df$event_code == 1])
  
  cat("  ", disease_name, "- Events:", total_events, 
      "(mCA:", events_mca, "No mCA:", events_no_mca, ")\n")
  
  if (nrow(df) > 0 && total_events > 10 && groups_with_events == 2) {
    
    # Gray's test (unadjusted)
    tryCatch({
      cfit_test <- cmprsk::cuminc(
        ftime = df$time_years,
        fstatus = df$event_code,
        group = df$group
      )
      
      if (!is.null(cfit_test$Tests) && nrow(cfit_test$Tests) > 0) {
        gray_p <- cfit_test$Tests[1, "pv"]
      }
    }, error = function(e) {
      cat("    Gray's test failed:", e$message, "\n")
    })
    
    # Unadjusted Fine-Gray regression
    tryCatch({
      df$group_numeric <- as.numeric(df$group == "mCA")
      
      crr_fit_unadj <- cmprsk::crr(
        ftime = df$time_years,
        fstatus = df$event_code,
        cov1 = matrix(df$group_numeric, ncol = 1),
        failcode = 1,
        cencode = 0
      )
      
      if (!is.null(crr_fit_unadj$coef) && !is.na(crr_fit_unadj$coef[1])) {
        hr_unadj <- exp(crr_fit_unadj$coef[1])
        se <- sqrt(crr_fit_unadj$var[1,1])
        hr_unadj_lower <- exp(crr_fit_unadj$coef[1] - 1.96 * se)
        hr_unadj_upper <- exp(crr_fit_unadj$coef[1] + 1.96 * se)
      }
      
    }, error = function(e) {
      cat("    Unadjusted CRR failed:", e$message, "\n")
    })
    
    # Adjusted Fine-Gray regression (age, sex, baseline disease burden)
    tryCatch({
      df$sex_numeric <- as.numeric(as.factor(df$sex)) - 1
      
      # Covariates: [group, age, sex, n_baseline_diseases]
      cov_matrix <- cbind(
        group = df$group_numeric,
        age = df$age,
        sex = df$sex_numeric,
        n_baseline_diseases = df$n_baseline_diseases
      )
      
      crr_fit_adj <- cmprsk::crr(
        ftime = df$time_years,
        fstatus = df$event_code,
        cov1 = cov_matrix,
        failcode = 1,
        cencode = 0
      )
      
      if (!is.null(crr_fit_adj$coef) && !is.na(crr_fit_adj$coef[1])) {
        hr_adj <- exp(crr_fit_adj$coef[1])
        
        # Get p-value
        z_score <- crr_fit_adj$coef[1] / sqrt(crr_fit_adj$var[1,1])
        crr_p <- 2 * (1 - pnorm(abs(z_score)))
        
        # Get CI
        se <- sqrt(crr_fit_adj$var[1,1])
        hr_adj_lower <- exp(crr_fit_adj$coef[1] - 1.96 * se)
        hr_adj_upper <- exp(crr_fit_adj$coef[1] + 1.96 * se)
        
        cat("    Adjusted HR:", round(hr_adj, 3), 
            "95% CI: [", round(hr_adj_lower, 3), "-", round(hr_adj_upper, 3), "]",
            "p =", format.pval(crr_p, digits = 3), "\n")
      }
      
    }, error = function(e) {
      cat("    Adjusted CRR failed:", e$message, "\n")
    })
  } else {
    cat("    Skipping analysis (insufficient events or groups)\n")
  }
  
  return(data.frame(
    Disease_Category = disease_name,
    N_mCA = n_mca,
    N_No_mCA = n_no_mca,
    Events_mCA = events_mca,
    Events_No_mCA = events_no_mca,
    HR_Unadjusted = hr_unadj,
    HR_Unadj_95CI_Lower = hr_unadj_lower,
    HR_Unadj_95CI_Upper = hr_unadj_upper,
    HR_Adjusted = hr_adj,
    HR_Adj_95CI_Lower = hr_adj_lower,
    HR_Adj_95CI_Upper = hr_adj_upper,
    Grays_Test_P = gray_p,
    CRR_Adjusted_P = crr_p,
    stringsAsFactors = FALSE
  ))
}

# ===========================
# Analyze disease categories (SKIP ALL ICD-10)
# ===========================
cat("=== STARTING COMPETING RISK ANALYSES ===\n")
cat("Note: Skipping ALL ICD-10 analysis for faster results\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

results_df <- data.frame()

# A+B Combined (Infectious diseases)
cat("[1/13] Analyzing ICD-10 A+B (Infectious diseases)\n")
incident_ab <- get_incident_diseases(diag_data, mca_info, 
                                     filter_condition = "substr(diagnosis, 1, 1) %in% c('A', 'B')")
df_ab <- mca_info %>% left_join(incident_ab, by = c("ID" = "patient_id"))
results_df <- rbind(results_df, perform_competing_risk_analysis(df_ab, "ICD10_A_B_Infectious"))

# C+D00-D48 (All Neoplasms)
cat("\n[2/13] Analyzing ICD-10 C+D00-D48 (Neoplasms)\n")
incident_neoplasm <- get_incident_diseases(diag_data, mca_info,
                                           filter_condition = "substr(diagnosis, 1, 1) == 'C' | (substr(diagnosis, 1, 1) == 'D' & substr(diagnosis, 2, 2) %in% c('0', '1', '2', '3', '4') & as.numeric(substr(diagnosis, 2, 3)) <= 48)")
df_neoplasm <- mca_info %>% left_join(incident_neoplasm, by = c("ID" = "patient_id"))
results_df <- rbind(results_df, perform_competing_risk_analysis(df_neoplasm, "ICD10_C_D00_D48_Neoplasms"))

# D50-D89 (Blood Diseases)
cat("\n[3/13] Analyzing ICD-10 D50-D89 (Blood disorders)\n")
incident_blood <- get_incident_diseases(diag_data, mca_info,
                                        filter_condition = "substr(diagnosis, 1, 1) == 'D' & substr(diagnosis, 2, 2) %in% c('5', '6', '7', '8') & as.numeric(substr(diagnosis, 2, 3)) >= 50 & as.numeric(substr(diagnosis, 2, 3)) <= 89")
df_blood <- mca_info %>% left_join(incident_blood, by = c("ID" = "patient_id"))
results_df <- rbind(results_df, perform_competing_risk_analysis(df_blood, "ICD10_D50_D89_Blood"))

# Individual chapters E, F, G, H, I, J, K, L, M, N (SEQUENTIAL)
individual_chapters <- c("E", "F", "G", "H", "I", "J", "K", "L", "M", "N")
chapter_names <- c(
  "E" = "Endocrine/Metabolic",
  "F" = "Mental/Behavioral",
  "G" = "Nervous System",
  "H" = "Eye/Ear",
  "I" = "Circulatory",
  "J" = "Respiratory",
  "K" = "Digestive",
  "L" = "Skin",
  "M" = "Musculoskeletal",
  "N" = "Genitourinary"
)

cat("\n[4-13] Analyzing Individual ICD-10 Chapters (Sequential)\n\n")

for (i in seq_along(individual_chapters)) {
  ch <- individual_chapters[i]
  ch_name <- chapter_names[ch]
  
  cat(paste0("[", i+3, "/13] Chapter ", ch, " (", ch_name, ")\n"))
  
  # Get incident diseases for this chapter
  incident_ch <- get_incident_diseases(diag_data, mca_info,
                                       filter_condition = paste0("substr(diagnosis, 1, 1) == '", ch, "'"))
  
  df_ch <- mca_info %>% left_join(incident_ch, by = c("ID" = "patient_id"))
  
  result <- perform_competing_risk_analysis(df_ch, paste0("ICD10_Chapter_", ch))
  results_df <- rbind(results_df, result)
  cat("\n")
}

cat("=== All analyses complete ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ===========================
# FDR correction
# ===========================
cat("=== Applying FDR Correction ===\n")

results_df$Grays_Test_Q <- NA
results_df$CRR_Adjusted_Q <- NA

# FDR for Gray's test
valid_gray_p <- !is.na(results_df$Grays_Test_P)
if (sum(valid_gray_p) > 0) {
  results_df$Grays_Test_Q[valid_gray_p] <- p.adjust(results_df$Grays_Test_P[valid_gray_p], method = "BH")
  cat("Gray's test FDR correction applied to", sum(valid_gray_p), "tests\n")
}

# FDR for adjusted CRR
valid_crr_p <- !is.na(results_df$CRR_Adjusted_P)
if (sum(valid_crr_p) > 0) {
  results_df$CRR_Adjusted_Q[valid_crr_p] <- p.adjust(results_df$CRR_Adjusted_P[valid_crr_p], method = "BH")
  cat("CRR adjusted FDR correction applied to", sum(valid_crr_p), "tests\n")
}

# ===========================
# Format final results
# ===========================
cat("\n=== Formatting results ===\n")

results_final <- results_df %>%
  mutate(
    HR_Unadj_95CI = ifelse(
      !is.na(HR_Unadjusted) & !is.na(HR_Unadj_95CI_Lower) & !is.na(HR_Unadj_95CI_Upper),
      paste0(sprintf("%.3f", HR_Unadjusted), " (",
             sprintf("%.3f", HR_Unadj_95CI_Lower), "-",
             sprintf("%.3f", HR_Unadj_95CI_Upper), ")"),
      NA
    ),
    HR_Adj_95CI = ifelse(
      !is.na(HR_Adjusted) & !is.na(HR_Adj_95CI_Lower) & !is.na(HR_Adj_95CI_Upper),
      paste0(sprintf("%.3f", HR_Adjusted), " (",
             sprintf("%.3f", HR_Adj_95CI_Lower), "-",
             sprintf("%.3f", HR_Adj_95CI_Upper), ")"),
      NA
    ),
    Grays_Test_P = ifelse(!is.na(Grays_Test_P), sprintf("%.6f", Grays_Test_P), NA),
    Grays_Test_Q = ifelse(!is.na(Grays_Test_Q), sprintf("%.6f", Grays_Test_Q), NA),
    CRR_Adjusted_P = ifelse(!is.na(CRR_Adjusted_P), sprintf("%.6f", CRR_Adjusted_P), NA),
    CRR_Adjusted_Q = ifelse(!is.na(CRR_Adjusted_Q), sprintf("%.6f", CRR_Adjusted_Q), NA)
  ) %>%
  select(
    Disease_Category, N_mCA, N_No_mCA, Events_mCA, Events_No_mCA,
    HR_Unadjusted, HR_Unadj_95CI_Lower, HR_Unadj_95CI_Upper, HR_Unadj_95CI,
    HR_Adjusted, HR_Adj_95CI_Lower, HR_Adj_95CI_Upper, HR_Adj_95CI,
    Grays_Test_P, Grays_Test_Q, CRR_Adjusted_P, CRR_Adjusted_Q
  ) %>%
  arrange(Disease_Category)

# Print results
cat("\n=== FINAL RESULTS ===\n\n")
print(results_final, row.names = FALSE)

# ===========================
# Save results
# ===========================
cat("\n=== Saving results ===\n")

excel_path <- file.path(output_dir, "Competing_Risk_Analysis_ICD10_Chapters_Only.xlsx")

if (excel_package == "openxlsx") {
  wb <- createWorkbook()
  addWorksheet(wb, "Summary_Results")
  writeData(wb, "Summary_Results", results_final)
  addWorksheet(wb, "Detailed_Results")
  writeData(wb, "Detailed_Results", results_df)
  saveWorkbook(wb, excel_path, overwrite = TRUE)
  cat("✅ Results saved to Excel:", excel_path, "\n")
} else if (excel_package == "writexl") {
  write_xlsx(list(
    Summary_Results = results_final,
    Detailed_Results = results_df
  ), excel_path)
  cat("✅ Results saved to Excel:", excel_path, "\n")
} else {
  csv_path_summary <- file.path(output_dir, "Competing_Risk_Analysis_ICD10_Chapters_Summary.csv")
  csv_path_detailed <- file.path(output_dir, "Competing_Risk_Analysis_ICD10_Chapters_Detailed.csv")
  write.csv(results_final, csv_path_summary, row.names = FALSE)
  write.csv(results_df, csv_path_detailed, row.names = FALSE)
  cat("✅ Results saved to CSV files\n")
}

cat("\n🎉 ANALYSIS COMPLETE! 🎉\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
