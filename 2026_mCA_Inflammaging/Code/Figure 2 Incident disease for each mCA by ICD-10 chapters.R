# ===============================
# FULL ICD-10 CHAPTERS COMPETING RISK ANALYSIS - UPDATED
# All 13 ICD chapters, all mCA subtypes (gain/loss/neutral)
# Fine-Gray model with death as competing risk
# Adjusted for age, sex, and baseline disease burden
# With robust interim saving and R file output
# ===============================

# FIRST: Stop any running parallel clusters
try(stopCluster(cl), silent = TRUE)
rm(list = ls())

# -------------------------------
# Libraries
# -------------------------------
library(readxl)
library(data.table)
library(dplyr)
library(lubridate)
library(survival)
library(tidyr)
library(writexl)
library(stringr)
library(tibble)
library(purrr)
library(doParallel)
library(foreach)
library(cmprsk)

# -------------------------------
# Config
# -------------------------------
censor_date <- as.Date("2023-07-14")
output_dir <- "~/Desktop/UK BB/HR_mCA_ICD10_chapters_competing_risk_FULL_UPDATED"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Thresholds for unadjusted model (1 covariate)
min_events_unadj <- 10
min_events_per_arm_unadj <- 5

# Thresholds for adjusted model (4 covariates: exposed, age, sex, n_baseline_diseases)
min_events_adj <- 30
min_events_per_arm_adj <- 5

# Progress tracking
analysis_start_time <- Sys.time()

# -------------------------------
# Load data with UPDATED FILTERING (matching second code)
# -------------------------------
cat("=== LOADING DATA WITH UPDATED FILTERING ===\n")
cat("Start time:", format(analysis_start_time), "\n\n")

if (!exists("cohort")) {
  
  # Step 1: Load analyzed samples (NEW - from second code)
  cat("Step 1: Loading analyzed samples...\n")
  analyzed <- read_excel("~/Desktop/UK BB/ukb4777.samples_analyzed.xlsx") %>%
    mutate(ID = as.character(ID))
  cat("✅ Total analyzed samples:", nrow(analyzed), "\n\n")
  
  # Step 2: Load participant info with age and sex
  cat("Step 2: Loading participant data...\n")
  mca_info_raw <- read_excel("~/Desktop/UK BB/Age sex dob collection date.xlsx") %>%
    mutate(
      ID = as.character(ID),
      mCA_date = mdy(`Date of mCA assessment`),
      age = as.numeric(`Age at blood collection`)
    ) %>%
    filter(ID %in% analyzed$ID)  # Filter to analyzed samples only
  
  cat("Participant records (analyzed samples only):", nrow(mca_info_raw), "\n")
  
  # Filter out invalid dates
  mca_info <- mca_info_raw %>%
    filter(!is.na(mCA_date),
           mCA_date != as.Date("1900-01-01"),
           mCA_date > as.Date("1901-01-01"),
           mCA_date <= censor_date) %>%
    filter(!is.na(age)) %>%
    mutate(
      sex = case_when(
        Sex == 0 ~ "Female",
        Sex == 1 ~ "Male",
        TRUE ~ NA_character_
      ),
      sex = factor(sex, levels = c("Female", "Male"))
    ) %>%
    filter(!is.na(sex)) %>%
    select(ID, mCA_date, age, sex) %>%
    distinct(ID, .keep_all = TRUE)
  
  cat("Valid participants after date filtering:", nrow(mca_info), "\n\n")
  
  # Step 3: Load CNV data and determine mCA status (UPDATED - matching second code)
  cat("Step 3: Loading CNV data and determining mCA status...\n")
  cnv_all <- read_excel("~/Desktop/UK BB/CNV_data_19808.xlsx") %>%
    mutate(ID = as.character(ID)) %>%
    filter(ID %in% analyzed$ID)
  
  # Summarize CNV status for each patient (NEW - from second code)
  cnv_summary <- cnv_all %>%
    group_by(ID) %>%
    summarise(
      has_known_cnv = any(COPY_CHANGE %in% c("gain", "loss", "neutral")),
      has_only_unknown = all(COPY_CHANGE == "unknown"),
      .groups = "drop"
    )
  
  # Determine mCA status (NEW - exclude "unknown only" patients)
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
  
  # Step 4: Create mCA subtype indicators
  cat("Step 4: Creating mCA subtype indicators...\n")
  cnv_data <- cnv_all %>%
    filter(COPY_CHANGE %in% c("gain", "loss", "neutral")) %>%
    mutate(mca_subtype = paste0("chr", CHR, "_", COPY_CHANGE)) %>%
    distinct(ID, mca_subtype)
  
  all_subtypes <- sort(unique(cnv_data$mca_subtype))
  cat("mCA subtypes found:", length(all_subtypes), "\n")
  
  mca_wide <- cnv_data %>%
    mutate(value = 1L) %>%
    tidyr::pivot_wider(id_cols = ID,
                       names_from = mca_subtype,
                       values_from = value,
                       values_fill = 0L)
  
  # Step 5: Build cohort with mCA status and No_mCA indicator
  cat("Step 5: Building final cohort...\n")
  
  # Filter mca_info to only include patients with definitive mCA status
  mca_info <- mca_info %>%
    filter(ID %in% cnv_status$ID)
  
  # Create No_mCA indicator based on cnv_status
  no_mca_ids <- cnv_status %>%
    filter(mCA_status == "no_mCA") %>%
    pull(ID)
  
  cohort <- mca_info %>%
    left_join(mca_wide, by = "ID") %>%
    mutate(across(any_of(all_subtypes), ~replace_na(., 0L))) %>%
    mutate(No_mCA = ifelse(ID %in% no_mca_ids, 1L, 0L))
  
  # Ensure all subtype columns exist
  for (st in all_subtypes) {
    if (!st %in% names(cohort)) {
      cohort[[st]] <- 0L
    }
  }
  
  cat("✅ Final cohort size:", nrow(cohort), "participants\n")
  cat("   - mCA participants:", sum(cohort$No_mCA == 0L), "\n")
  cat("   - No mCA participants:", sum(cohort$No_mCA == 1L), "\n\n")
}

# Step 6: Load diagnosis data
if (!exists("diag_data")) {
  cat("Step 6: Loading diagnosis data...\n")
  diag_data <- fread("~/Desktop/UK BB/Diagnosis code date.tsv", sep = "\t") %>%
    mutate(
      patient_id     = as.character(patient_id),
      admission_date = as.Date(admission_date),
      diagnosis      = str_replace_all(toupper(str_trim(as.character(diagnosis))), "\\.", "")
    ) %>%
    filter(!is.na(admission_date), !is.na(diagnosis), !is.na(patient_id)) %>%
    select(patient_id, admission_date, diagnosis)
  cat("✅ Diagnosis data loaded:", format(nrow(diag_data), big.mark = ","), "records\n\n")
}

# Step 7: Load death data
if (!exists("death_data")) {
  cat("Step 7: Loading death data...\n")
  death_data <- read_excel("~/Desktop/UK BB/Death date.xlsx") %>%
    mutate(
      patient_id = as.character(Participant_ID),
      death_date = as.Date(date_of_death)
    ) %>%
    filter(!is.na(death_date), !is.na(patient_id)) %>%
    select(patient_id, death_date) %>%
    distinct(patient_id, .keep_all = TRUE)
  cat("✅ Death data loaded:", nrow(death_data), "records\n\n")
}

# -------------------------------
# Calculate baseline disease burden (NEW - from second code)
# -------------------------------
cat("=== Calculating baseline disease burden ===\n")

# Convert to data.table for speed
dt_diag <- as.data.table(diag_data)
dt_cohort <- as.data.table(cohort)

# Set keys for fast joining
setkey(dt_diag, patient_id)
setkey(dt_cohort, ID)

# Filter diagnosis data to only include patients in cohort FIRST
cat("Filtering diagnosis data to study participants...\n")
diag_data_filtered <- dt_diag[patient_id %in% dt_cohort$ID]
cat("Filtered diagnosis records:", format(nrow(diag_data_filtered), big.mark = ","), "\n")

# Join with mCA dates
cat("Joining with mCA dates...\n")
diag_with_mca <- diag_data_filtered[dt_cohort[, .(ID, mCA_date)], 
                                    on = .(patient_id = ID), 
                                    nomatch = 0]

# Calculate baseline burden
cat("Calculating baseline disease counts...\n")
baseline_disease_burden <- diag_with_mca[admission_date <= mCA_date, 
                                         .(n_baseline_diseases = uniqueN(diagnosis)),
                                         by = patient_id]

cat("Baseline calculations complete for", nrow(baseline_disease_burden), "patients\n")

# Add baseline disease burden to cohort
cohort <- cohort %>%
  left_join(as.data.frame(baseline_disease_burden), by = c("ID" = "patient_id")) %>%
  mutate(n_baseline_diseases = replace_na(n_baseline_diseases, 0L))

cat("Mean baseline diseases - mCA:", round(mean(cohort$n_baseline_diseases[cohort$No_mCA == 0L]), 2), "\n")
cat("Mean baseline diseases - No mCA:", round(mean(cohort$n_baseline_diseases[cohort$No_mCA == 1L]), 2), "\n\n")

# -------------------------------
# ICD-10 CHAPTERS: 13 disease categories
# Note: A+B combined, C+D00-D49 combined, D50-D89 combined
#       E through N as individual chapters
#       No separate A, B, C, D analysis; No "all combined" category
# -------------------------------
icd_table <- tribble(
  ~Disease, ~Category, ~ICD10_pattern,
  "A+B: Infectious and parasitic diseases (A00-B99)", "ICD10_Chapters", "^[AB]",
  "C+D0-D4: Neoplasms (C00-D49)", "ICD10_Chapters", "^C|^D[0-4]",
  "D5-D8: Blood disorders (D50-D89)", "ICD10_Chapters", "^D[5-8]",
  "E: Endocrine and metabolic diseases", "ICD10_Chapters", "^E",
  "F: Mental and behavioral disorders", "ICD10_Chapters", "^F",
  "G: Nervous system diseases", "ICD10_Chapters", "^G",
  "H: Eye and ear diseases", "ICD10_Chapters", "^H",
  "I: Circulatory system diseases", "ICD10_Chapters", "^I",
  "J: Respiratory system diseases", "ICD10_Chapters", "^J",
  "K: Digestive system diseases", "ICD10_Chapters", "^K",
  "L: Skin and subcutaneous diseases", "ICD10_Chapters", "^L",
  "M: Musculoskeletal diseases", "ICD10_Chapters", "^M",
  "N: Genitourinary system diseases", "ICD10_Chapters", "^N"
) %>%
  mutate(idx = row_number())

cat("\n=== FULL ANALYSIS SETUP ===\n")
cat("Cohort size:", nrow(cohort), "\n")
cat("Strict No_mCA count:", sum(cohort$No_mCA == 1L, na.rm = TRUE), "\n")
cat("mCA subtypes (gain/loss/neutral):", length(all_subtypes), "\n")
cat("ICD chapters (all 13):", nrow(icd_table), "\n")
cat("Total models to run:", nrow(icd_table) * length(all_subtypes), "\n")

# Estimate runtime
models_per_hour <- 96 / 49.1
total_models <- nrow(icd_table) * length(all_subtypes)
estimated_hours <- total_models / models_per_hour
cat("Estimated runtime:", round(estimated_hours, 1), "hours (", round(estimated_hours/24, 1), "days)\n\n")

# -------------------------------
# Helper functions
# -------------------------------
# UPDATED: Get first INCIDENT diagnosis after mCA
# Only excludes specific diagnosis CODES that existed before mCA
# Allows new codes in the same chapter
first_dx_after_mCA_pattern <- function(icd_pattern, diag_tbl, cohort_tbl) {
  
  # Convert to data.table for speed
  dt_diag <- as.data.table(diag_tbl)
  dt_cohort <- as.data.table(cohort_tbl)
  
  # Filter to cohort and matching ICD pattern
  dt_diag_filtered <- dt_diag[patient_id %in% dt_cohort$ID]
  dt_diag_filtered <- dt_diag_filtered[grepl(icd_pattern, diagnosis, perl = TRUE)]
  
  # Join with mCA dates
  dt_diag_filtered <- dt_diag_filtered[dt_cohort[, .(ID, mCA_date)], 
                                       on = .(patient_id = ID), 
                                       nomatch = 0]
  
  # Get pre-existing (patient_id, diagnosis) pairs - SPECIFIC CODES before mCA
  pre_existing <- dt_diag_filtered[admission_date <= mCA_date, .(patient_id, diagnosis)]
  pre_existing <- unique(pre_existing)
  
  # Get post-mCA diagnoses
  post_mca <- dt_diag_filtered[admission_date > mCA_date]
  
  # Anti-join: remove only those specific diagnosis CODES that existed before
  setkey(post_mca, patient_id, diagnosis)
  setkey(pre_existing, patient_id, diagnosis)
  incident <- post_mca[!pre_existing]
  
  # Get first date from remaining truly NEW codes
  if (nrow(incident) == 0) {
    return(data.frame(patient_id = character(), first_dx_date = as.Date(character())))
  }
  
  result <- incident[, .(first_dx_date = min(admission_date)), by = patient_id]
  
  return(as.data.frame(result))
}

first_death_after_mCA <- function(death_tbl, cohort_ids) {
  death_tbl %>%
    filter(patient_id %in% cohort_ids) %>%
    group_by(patient_id) %>%
    summarise(death_date = suppressWarnings(min(death_date, na.rm = TRUE)), .groups = "drop") %>%
    filter(!is.na(death_date))
}

# UPDATED: Fine-Gray model with baseline disease burden adjustment
safe_finegray <- function(dat, formula_str) {
  tryCatch({
    cov_matrix <- model.matrix(as.formula(paste0("~", gsub(".*~", "", formula_str))), data = dat)[,-1, drop=FALSE]
    
    if (nrow(dat) < 10 || ncol(cov_matrix) == 0) return(NULL)
    
    fit <- crr(ftime = dat$time_years, 
               fstatus = dat$event_type, 
               cov1 = cov_matrix,
               failcode = 1)
    
    hr <- exp(fit$coef)
    se <- sqrt(diag(fit$var))
    ci_lower <- exp(fit$coef - 1.96 * se)
    ci_upper <- exp(fit$coef + 1.96 * se)
    pval <- 2 * (1 - pnorm(abs(fit$coef / se)))
    
    list(
      hr = hr[1], 
      ci = c(ci_lower[1], ci_upper[1]), 
      p = pval[1]
    )
  }, error = function(e) {
    return(NULL)
  })
}

create_na_row <- function(label, category, icd_pattern, subtype, reason) {
  tibble(
    Disease = label, Category = category, ICD10_pattern = icd_pattern,
    mCA_subtype = subtype, n_total = NA_integer_, events_total = NA_integer_,
    deaths_total = NA_integer_, events_exposed = NA_integer_, events_unexposed = NA_integer_,
    deaths_exposed = NA_integer_, deaths_unexposed = NA_integer_,
    model_type = "Competing_Risk_vs_No_mCA", HR_unadj = NA_real_, CI_lower_unadj = NA_real_,
    CI_upper_unadj = NA_real_, p_unadj = NA_real_, HR_adj = NA_real_,
    CI_lower_adj = NA_real_, CI_upper_adj = NA_real_, p_adj = NA_real_,
    exclusion_reason = reason
  )
}

# -------------------------------
# Get death events for cohort
# -------------------------------
death_events <- first_death_after_mCA(death_data, cohort_ids = cohort$ID)
cat("Deaths after mCA:", nrow(death_events), "\n\n")

# -------------------------------
# ROBUST INTERIM SAVING FUNCTIONS
# -------------------------------
save_interim_results <- function(results_so_far, chapter_num, output_dir) {
  if (length(results_so_far) > 0) {
    interim_df <- bind_rows(results_so_far[!sapply(results_so_far, is.null)])
    if (nrow(interim_df) > 0) {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      
      # Save as Excel
      interim_excel <- file.path(output_dir, paste0("INTERIM_Chapter_", chapter_num, "_", timestamp, ".xlsx"))
      write_xlsx(interim_df, path = interim_excel)
      
      # Save as RDS (R binary format - faster to load)
      interim_rds <- file.path(output_dir, paste0("INTERIM_Chapter_", chapter_num, "_", timestamp, ".rds"))
      saveRDS(interim_df, file = interim_rds)
      
      cat("💾 Interim results saved:\n")
      cat("   Excel:", interim_excel, "\n")
      cat("   RDS:  ", interim_rds, "\n")
    }
  }
}

# Function to save checkpoint (allows resuming)
save_checkpoint <- function(results_list, completed_chapters, output_dir) {
  checkpoint_file <- file.path(output_dir, "CHECKPOINT_progress.rds")
  checkpoint_data <- list(
    results_list = results_list,
    completed_chapters = completed_chapters,
    timestamp = Sys.time()
  )
  saveRDS(checkpoint_data, file = checkpoint_file)
  cat("🔖 Checkpoint saved:", checkpoint_file, "\n")
}

# Function to load checkpoint (if exists)
load_checkpoint <- function(output_dir) {
  checkpoint_file <- file.path(output_dir, "CHECKPOINT_progress.rds")
  if (file.exists(checkpoint_file)) {
    checkpoint_data <- readRDS(checkpoint_file)
    cat("📂 Found checkpoint from:", format(checkpoint_data$timestamp), "\n")
    cat("   Completed chapters:", checkpoint_data$completed_chapters, "\n")
    return(checkpoint_data)
  }
  return(NULL)
}

# -------------------------------
# OPTIMIZED PARALLEL SETUP
# -------------------------------
n_cores <- min(8, parallel::detectCores() - 1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)

cat("🚀 Starting FULL analysis with", n_cores, "cores...\n")
cat("This is a large-scale analysis - consider running overnight or over a weekend\n\n")

# Export objects to cluster (UPDATED - includes baseline disease burden)
clusterExport(cl, c("cohort", "diag_data", "death_events", "icd_table", "all_subtypes"), envir = environment())
export_syms <- c("first_dx_after_mCA_pattern", "safe_finegray", "create_na_row",
                 "censor_date", "min_events_unadj", "min_events_per_arm_unadj",
                 "min_events_adj", "min_events_per_arm_adj")

start_time <- Sys.time()

# Check for existing checkpoint
checkpoint <- load_checkpoint(output_dir)
if (!is.null(checkpoint)) {
  results_list <- checkpoint$results_list
  start_chapter <- checkpoint$completed_chapters + 1
  cat("📌 Resuming from chapter", start_chapter, "\n\n")
} else {
  results_list <- vector("list", nrow(icd_table))
  start_chapter <- 1
}

completed_models <- (start_chapter - 1) * length(all_subtypes)
total_models <- nrow(icd_table) * length(all_subtypes)

# MAIN ANALYSIS LOOP
for (i in start_chapter:nrow(icd_table)) {
  cat("\n", rep("=", 60), "\n")
  cat("Processing ICD chapter", i, "of", nrow(icd_table), ":", icd_table$Disease[i], "\n")
  cat("Progress: ", round(100 * (i-1) / nrow(icd_table), 1), "% of chapters completed\n")
  chapter_start <- Sys.time()
  
  # Process this chapter in parallel across subtypes
  chapter_results <- foreach(
    j = seq_along(all_subtypes),
    .packages = c("dplyr", "tidyr", "stringr", "tibble", "lubridate", "purrr", "data.table", "cmprsk"),
    .export   = export_syms,
    .errorhandling = "pass"
  ) %dopar% {
    
    label       <- icd_table$Disease[i]
    icd_pattern <- icd_table$ICD10_pattern[i]
    category    <- icd_table$Category[i]
    subtype     <- all_subtypes[j]
    
    # Get disease events using pattern matching (excludes only specific pre-existing codes)
    ch_diag <- first_dx_after_mCA_pattern(icd_pattern, diag_tbl = diag_data, cohort_tbl = cohort)
    
    # Build disease dataset (post mCA only) with competing risk setup
    df_disease <- cohort %>%
      left_join(ch_diag, by = c("ID" = "patient_id")) %>%
      left_join(death_events, by = c("ID" = "patient_id")) %>%
      filter(is.na(first_dx_date) | first_dx_date > mCA_date) %>%
      filter(is.na(death_date) | death_date > mCA_date) %>%
      mutate(
        disease_after_mca = !is.na(first_dx_date) & first_dx_date > mCA_date,
        death_after_mca = !is.na(death_date) & death_date > mCA_date,
        
        first_disease_time = ifelse(disease_after_mca, as.numeric(difftime(first_dx_date, mCA_date, units = "days"))/365.25, Inf),
        first_death_time = ifelse(death_after_mca, as.numeric(difftime(death_date, mCA_date, units = "days"))/365.25, Inf),
        
        event_type = case_when(
          first_disease_time < first_death_time ~ 1L,
          first_death_time < first_disease_time ~ 2L,
          TRUE ~ 0L
        ),
        
        time_years = case_when(
          event_type == 1L ~ first_disease_time,
          event_type == 2L ~ first_death_time,
          TRUE ~ as.numeric(difftime(censor_date, mCA_date, units = "days"))/365.25
        )
      ) %>%
      filter(is.finite(time_years), time_years > 0, !is.na(age), !is.na(sex), !is.na(n_baseline_diseases))
    
    # Subset to exposed (specific mCA subtype) vs unexposed (No_mCA)
    df_sub <- df_disease %>%
      mutate(exposed = as.integer(.data[[subtype]] == 1L)) %>%
      filter(exposed == 1L | No_mCA == 1L)
    
    # Check if we have both groups
    if (nrow(df_sub) == 0 || !any(df_sub$exposed == 1L) || !any(df_sub$exposed == 0L)) {
      return(create_na_row(label, category, icd_pattern, subtype, "Insufficient_data_for_subtype"))
    }
    
    # Calculate event counts
    events_total <- sum(df_sub$event_type == 1L)
    e1 <- sum(df_sub$event_type[df_sub$exposed == 1L] == 1L)  # events in exposed
    e0 <- sum(df_sub$event_type[df_sub$exposed == 0L] == 1L)  # events in unexposed
    d1 <- sum(df_sub$event_type[df_sub$exposed == 1L] == 2L)  # deaths in exposed
    d0 <- sum(df_sub$event_type[df_sub$exposed == 0L] == 2L)  # deaths in unexposed
    
    # Initialize results
    fit_unadj <- NULL
    fit_adj <- NULL
    exclusion_reason <- character(0)
    
    # Check thresholds for UNADJUSTED model (>=10 total events, >=5 per arm)
    can_run_unadj <- (events_total >= min_events_unadj) && 
      (e1 >= min_events_per_arm_unadj) && 
      (e0 >= min_events_per_arm_unadj)
    
    if (can_run_unadj) {
      fit_unadj <- safe_finegray(df_sub, "~ exposed")
      if (is.null(fit_unadj)) {
        exclusion_reason <- c(exclusion_reason, "Unadj_model_failed")
      }
    } else {
      exclusion_reason <- c(exclusion_reason, 
                            paste0("Unadj_threshold_fail_total=", events_total, "_e1=", e1, "_e0=", e0))
    }
    
    # Check thresholds for ADJUSTED model (>=30 total events, >=5 per arm)
    can_run_adj <- (events_total >= min_events_adj) && 
      (e1 >= min_events_per_arm_adj) && 
      (e0 >= min_events_per_arm_adj)
    
    if (can_run_adj) {
      fit_adj <- safe_finegray(df_sub, "~ exposed + age + sex + n_baseline_diseases")
      if (is.null(fit_adj)) {
        exclusion_reason <- c(exclusion_reason, "Adj_model_failed")
      }
    } else {
      exclusion_reason <- c(exclusion_reason, 
                            paste0("Adj_threshold_fail_total=", events_total, "_e1=", e1, "_e0=", e0))
    }
    
    # If both models failed, return NA row
    if (is.null(fit_unadj) && is.null(fit_adj)) {
      return(create_na_row(label, category, icd_pattern, subtype, 
                           paste(exclusion_reason, collapse = "; ")))
    }
    
    # Build result tibble (allow partial results)
    tibble(
      Disease = label, Category = category, ICD10_pattern = icd_pattern,
      mCA_subtype = subtype, n_total = nrow(df_sub), events_total = events_total,
      deaths_total = sum(df_sub$event_type == 2L), events_exposed = e1, events_unexposed = e0,
      deaths_exposed = d1, deaths_unexposed = d0, model_type = "Competing_Risk_vs_No_mCA",
      HR_unadj = if (!is.null(fit_unadj)) unname(fit_unadj$hr) else NA_real_,
      CI_lower_unadj = if (!is.null(fit_unadj)) unname(fit_unadj$ci[1]) else NA_real_,
      CI_upper_unadj = if (!is.null(fit_unadj)) unname(fit_unadj$ci[2]) else NA_real_,
      p_unadj = if (!is.null(fit_unadj)) unname(fit_unadj$p) else NA_real_,
      HR_adj = if (!is.null(fit_adj)) unname(fit_adj$hr) else NA_real_,
      CI_lower_adj = if (!is.null(fit_adj)) unname(fit_adj$ci[1]) else NA_real_,
      CI_upper_adj = if (!is.null(fit_adj)) unname(fit_adj$ci[2]) else NA_real_,
      p_adj = if (!is.null(fit_adj)) unname(fit_adj$p) else NA_real_,
      exclusion_reason = if (length(exclusion_reason) > 0) paste(exclusion_reason, collapse = "; ") else "Analysis_successful"
    )
  }
  
  results_list[[i]] <- bind_rows(chapter_results)
  completed_models <- completed_models + length(all_subtypes)
  
  chapter_time <- round(difftime(Sys.time(), chapter_start, units = "mins"), 2)
  elapsed_total <- round(difftime(Sys.time(), start_time, units = "hours"), 2)
  
  cat("✅ Chapter", i, "completed in", chapter_time, "minutes\n")
  cat("📊 Models completed:", completed_models, "of", total_models, "(", round(100*completed_models/total_models, 1), "%)\n")
  cat("⏱️  Total elapsed time:", elapsed_total, "hours\n")
  
  # ROBUST INTERIM SAVING: Save after EVERY chapter
  save_interim_results(results_list[1:i], i, output_dir)
  save_checkpoint(results_list, i, output_dir)
  
  # Estimate remaining time
  if (i > start_chapter) {
    chapters_done <- i - start_chapter + 1
    avg_time_per_chapter <- as.numeric(difftime(Sys.time(), start_time, units = "hours")) / chapters_done
    remaining_chapters <- nrow(icd_table) - i
    estimated_remaining <- remaining_chapters * avg_time_per_chapter
    cat("🕐 Estimated time remaining:", round(estimated_remaining, 1), "hours\n")
  }
}

stopCluster(cl)
end_time <- Sys.time()
total_runtime <- round(difftime(end_time, start_time, units = "hours"), 2)

cat("\n🎉 FULL ANALYSIS COMPLETED!\n")
cat("Total runtime:", total_runtime, "hours (", round(total_runtime/24, 1), "days)\n\n")

# -------------------------------
# Combine and save final results
# -------------------------------
if (length(results_list) > 0) {
  results_df <- bind_rows(results_list) %>%
    arrange(Category, Disease, mCA_subtype)
  
  # Apply FDR correction
  results_df$q_unadj <- NA_real_
  results_df$q_adj   <- NA_real_
  
  idx_u <- which(!is.na(results_df$p_unadj))
  idx_a <- which(!is.na(results_df$p_adj))
  
  if (length(idx_u) > 0) results_df$q_unadj[idx_u] <- p.adjust(results_df$p_unadj[idx_u], method = "fdr")
  if (length(idx_a) > 0) results_df$q_adj[idx_a]   <- p.adjust(results_df$p_adj[idx_a],   method = "fdr")
  
  # Save final results - multiple formats
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
  
  # 1. Excel file
  output_excel <- file.path(output_dir, paste0("FINAL_HR_Results_ICD10_All13Chapters_FullmCA_", timestamp, ".xlsx"))
  write_xlsx(results_df, path = output_excel)
  cat("✅ Final results saved to Excel:", output_excel, "\n")
  
  # 2. RDS file (R binary - preserves all data types)
  output_rds <- file.path(output_dir, paste0("FINAL_HR_Results_ICD10_All13Chapters_FullmCA_", timestamp, ".rds"))
  saveRDS(results_df, file = output_rds)
  cat("✅ Final results saved to RDS:", output_rds, "\n")
  
  # 3. RData file (can save multiple objects)
  output_rdata <- file.path(output_dir, paste0("FINAL_HR_Results_ICD10_All13Chapters_FullmCA_", timestamp, ".RData"))
  analysis_metadata <- list(
    analysis_date = Sys.time(),
    censor_date = censor_date,
    cohort_size = nrow(cohort),
    n_mca_subtypes = length(all_subtypes),
    n_icd_chapters = nrow(icd_table),
    total_runtime_hours = total_runtime,
    thresholds = list(
      unadjusted = list(
        min_events_total = min_events_unadj,
        min_events_per_arm = min_events_per_arm_unadj
      ),
      adjusted = list(
        min_events_total = min_events_adj,
        min_events_per_arm = min_events_per_arm_adj
      )
    ),
    adjustment_variables = c("age", "sex", "n_baseline_diseases")
  )
  save(results_df, analysis_metadata, all_subtypes, icd_table, 
       file = output_rdata)
  cat("✅ Final results saved to RData:", output_rdata, "\n\n")
  
  # Clean up checkpoint file after successful completion
  checkpoint_file <- file.path(output_dir, "CHECKPOINT_progress.rds")
  if (file.exists(checkpoint_file)) {
    file.remove(checkpoint_file)
    cat("🧹 Checkpoint file removed (analysis completed successfully)\n\n")
  }
  
  # Summary statistics
  total_tests <- nrow(results_df)
  successful_unadj <- sum(!is.na(results_df$p_unadj))
  successful_adj <- sum(!is.na(results_df$p_adj))
  significant_unadj_05 <- sum(results_df$p_unadj < 0.05, na.rm = TRUE)
  significant_adj_05 <- sum(results_df$p_adj < 0.05, na.rm = TRUE)
  significant_adj_01 <- sum(results_df$p_adj < 0.01, na.rm = TRUE)
  significant_fdr_unadj <- sum(results_df$q_unadj < 0.05, na.rm = TRUE)
  significant_fdr_adj <- sum(results_df$q_adj < 0.05, na.rm = TRUE)
  
  cat("=== FULL ANALYSIS SUMMARY ===\n")
  cat("Total tests performed:", total_tests, "\n\n")
  cat("UNADJUSTED MODEL (threshold: >=", min_events_unadj, " total events, >=", min_events_per_arm_unadj, " per arm)\n", sep = "")
  cat("  Successful analyses:", successful_unadj, "(", round(successful_unadj/total_tests*100, 1), "%)\n")
  cat("  Significant at p<0.05:", significant_unadj_05, "\n")
  cat("  Significant after FDR:", significant_fdr_unadj, "\n\n")
  cat("ADJUSTED MODEL (threshold: >=", min_events_adj, " total events, >=", min_events_per_arm_adj, " per arm)\n", sep = "")
  cat("  Successful analyses:", successful_adj, "(", round(successful_adj/total_tests*100, 1), "%)\n")
  cat("  Significant at p<0.05:", significant_adj_05, "\n")
  cat("  Significant at p<0.01:", significant_adj_01, "\n")
  cat("  Significant after FDR:", significant_fdr_adj, "\n\n")
  
  if (successful_adj > 0) {
    cat("=== TOP 10 MOST SIGNIFICANT RESULTS ===\n")
    top_results <- results_df %>%
      filter(!is.na(p_adj)) %>%
      arrange(p_adj) %>%
      select(Disease, mCA_subtype, HR_adj, CI_lower_adj, CI_upper_adj, p_adj, q_adj) %>%
      head(10)
    print(top_results)
    
    cat("\n=== RESULTS BY DISEASE CATEGORY ===\n")
    category_summary <- results_df %>%
      filter(!is.na(p_adj)) %>%
      group_by(Disease) %>%
      summarise(
        total_tests = n(),
        successful = sum(!is.na(p_adj)),
        sig_p05 = sum(p_adj < 0.05, na.rm = TRUE),
        sig_fdr = sum(q_adj < 0.05, na.rm = TRUE),
        median_hr = median(HR_adj, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(sig_fdr))
    print(category_summary)
    
    cat("\n=== TOTAL EVENT SUMMARY ===\n")
    event_summary <- results_df %>%
      filter(!is.na(p_adj)) %>%
      summarise(
        total_disease_events = sum(events_total, na.rm = TRUE),
        total_deaths = sum(deaths_total, na.rm = TRUE),
        .groups = "drop"
      )
    print(event_summary)
  }
  
  cat("\n🎯 Analysis complete! Check the output directory for all results.\n")
  cat("📁 Output directory:", output_dir, "\n")
} else {
  cat("⚠️ No results generated\n")
}

analysis_end_time <- Sys.time()
cat("\n=== FULL ICD-10 ANALYSIS COMPLETE ===\n")
cat("Analysis started:", format(analysis_start_time), "\n")
cat("Analysis ended:  ", format(analysis_end_time), "\n")
cat("Total duration: ", round(difftime(analysis_end_time, analysis_start_time, units = "hours"), 2), "hours\n")
cat("\n=== OUTPUT FILES ===\n")
cat("1. Excel: ", output_excel, "\n")
cat("2. RDS:   ", output_rds, "\n")
cat("3. RData: ", output_rdata, "\n")
