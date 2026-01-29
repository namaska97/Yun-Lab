# ==============================
# Cox Regression with LASSO Feature Selection
# Adjusted for Age, Sex, and Baseline Disease Burden
# ==============================

# ===========================
# Clear workspace
# ===========================
rm(list = ls())

# ==============================
# Libraries
# ==============================
library(readxl)
library(dplyr)
library(ggplot2)
library(cowplot)
library(GenomicRanges)
library(grid)
library(openxlsx)
library(data.table)
library(tidyr)
library(survival)
library(lubridate)
library(tibble)
library(glmnet)

# ==============================
# Config
# ==============================
ukbb_dir <- "~/Desktop/UK BB/"
censor_date <- as.Date("2023-07-14")
cyto_path <- file.path(ukbb_dir, "cytoBand.txt")
cnv_path <- file.path(ukbb_dir, "CNV_data_19808.xlsx")
diag_path <- file.path(ukbb_dir, "Diagnosis code date.tsv")
age_path <- file.path(ukbb_dir, "Age sex dob collection date.xlsx")
analyzed_path <- file.path(ukbb_dir, "ukb4777.samples_analyzed.xlsx")
death_path <- file.path(ukbb_dir, "Death date.xlsx")

cat("\n=== STARTING LASSO-COX ANALYSIS ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# ==============================
# Load cytoband data
# ==============================
cat("Loading cytoband data...\n")
cyto_df <- read.table(cyto_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(cyto_df) <- c("chrom", "start", "end", "name", "gieStain")
cyto_df$chrom <- gsub("^chr", "", cyto_df$chrom)

hg19IdeogramCyto <- GRanges(
  seqnames = paste0("chr", cyto_df$chrom),
  ranges = IRanges(start = cyto_df$start, end = cyto_df$end),
  name = cyto_df$name,
  gieStain = cyto_df$gieStain
)

# ==============================
# Load analyzed samples (SAME AS COMPETING RISK CODE)
# ==============================
cat("Loading analyzed samples...\n")
analyzed <- read_excel(analyzed_path) %>%
  mutate(ID = as.character(ID))

cat("Total analyzed samples:", nrow(analyzed), "\n\n")

# ==============================
# Load sample info with age and sex (SAME FILTERING)
# ==============================
cat("Loading participant data...\n")
mca_info_raw <- read_excel(age_path) %>%
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

# ==============================
# Load CNV data and determine mCA status (SAME AS COMPETING RISK CODE)
# ==============================
cat("Loading CNV data and determining mCA status...\n")
cnv_all <- read_excel(cnv_path) %>%
  mutate(ID = as.character(ID)) %>%
  filter(ID %in% analyzed$ID)

# ==============================
# Exclude homogeneous mCAs from LASSO-Cox analysis
# ==============================
# Based on homogeneity analysis, the following mCAs were classified as homogeneous:
# chr12 Gain (12+), chr14 Gain (14+), chr15 Gain (15+), chr18 Gain (18+), chrX Loss (X-)
# These are excluded because cytoband-level analysis is uninformative for 
# whole-chromosome or whole-arm events with uniform genomic boundaries.

homogeneous_mCAs <- list(
  list(chr = "12", type = "gain"),
  list(chr = "14", type = "gain"),
  list(chr = "15", type = "gain"),
  list(chr = "18", type = "gain"),
  list(chr = "X", type = "loss")
)

cat("\n=== Excluding Homogeneous mCAs ===\n")
cat("Homogeneous mCAs to exclude: 12+, 14+, 15+, 18+, X-\n")

# Count before exclusion
n_before <- nrow(cnv_all)

# Remove homogeneous mCAs
for (hm in homogeneous_mCAs) {
  cnv_all <- cnv_all %>%
    filter(!(CHR == hm$chr & COPY_CHANGE == hm$type))
}

n_after <- nrow(cnv_all)
cat("CNV records before exclusion:", n_before, "\n")
cat("CNV records after exclusion:", n_after, "\n")
cat("Records removed:", n_before - n_after, "\n\n")

# Summarize CNV status for each patient
cnv_summary <- cnv_all %>%
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

# Merge mCA status with mca_info
mca_info <- mca_info %>%
  inner_join(cnv_status, by = "ID")

cat("Final cohort after excluding 'unknown' CNV patients:", nrow(mca_info), "\n")
cat("mCA participants:", sum(mca_info$mCA_status == "mCA"), "\n")
cat("No mCA participants:", sum(mca_info$mCA_status == "no_mCA"), "\n\n")

# Get IDs for each group
mca_ids <- mca_info %>% filter(mCA_status == "mCA") %>% pull(ID)
no_mca_ids <- mca_info %>% filter(mCA_status == "no_mCA") %>% pull(ID)

# ==============================
# Load death data
# ==============================
cat("Loading death data...\n")
death_info <- read_excel(death_path) %>%
  dplyr::rename(ID = Participant_ID, death_date = date_of_death) %>%
  mutate(ID = as.character(ID),
         death_date = as.Date(death_date)) %>%
  distinct(ID, .keep_all = TRUE)

mca_info <- mca_info %>%
  left_join(death_info, by = "ID")

cat("Deaths recorded:", sum(!is.na(mca_info$death_date)), "\n\n")

# ==============================
# Load diagnosis data
# ==============================
cat("Loading diagnosis data...\n")
diag_data <- fread(diag_path, sep = "\t") %>%
  mutate(ID = as.character(patient_id),
         diagnosis = as.character(diagnosis),
         admission_date = as.Date(admission_date)) %>%
  filter(!is.na(admission_date), !is.na(diagnosis))

cat("Total diagnosis records:", nrow(diag_data), "\n")
cat("Unique patients in diagnosis data:", n_distinct(diag_data$ID), "\n\n")

# ==============================
# Calculate baseline disease burden (SAME AS COMPETING RISK CODE)
# ==============================
cat("=== Calculating baseline disease burden ===\n")

# Convert to data.table for speed
setDT(diag_data)
dt_mca_info <- as.data.table(mca_info)

# Set keys for fast joining
setkey(diag_data, ID)
setkey(dt_mca_info, ID)

# Filter diagnosis data to only include patients in mca_info FIRST
cat("Filtering diagnosis data to study participants...\n")
diag_data_filtered <- diag_data[ID %in% mca_info$ID]
cat("Filtered diagnosis records:", nrow(diag_data_filtered), "\n")

# Join with mCA dates
cat("Joining with mCA dates...\n")
diag_with_mca <- diag_data_filtered[dt_mca_info[, .(ID, mCA_date)], 
                                    on = .(ID = ID), 
                                    nomatch = 0]

# Calculate baseline burden
cat("Calculating baseline disease counts...\n")
baseline_disease_burden <- diag_with_mca[admission_date <= mCA_date, 
                                         .(n_baseline_diseases = uniqueN(diagnosis),
                                           first_disease_date = min(admission_date)),
                                         by = ID
]

cat("Baseline calculations complete for", nrow(baseline_disease_burden), "patients\n")

# Convert back to data.frame for compatibility
baseline_disease_burden <- as.data.frame(baseline_disease_burden)
mca_info <- as.data.frame(mca_info)

# Add to mca_info
mca_info <- mca_info %>%
  left_join(baseline_disease_burden, by = "ID") %>%
  mutate(
    n_baseline_diseases = replace_na(n_baseline_diseases, 0),
    time_from_first_disease = if_else(
      !is.na(first_disease_date),
      as.numeric(difftime(mCA_date, first_disease_date, units = "days")) / 365.25,
      NA_real_
    )
  )

cat("Baseline disease burden calculated\n")
cat("Mean baseline diseases - mCA:", round(mean(mca_info$n_baseline_diseases[mca_info$mCA_status == "mCA"]), 2), "\n")
cat("Mean baseline diseases - No mCA:", round(mean(mca_info$n_baseline_diseases[mca_info$mCA_status == "no_mCA"]), 2), "\n\n")

# ==============================
# Define disease categories (SAME AS COMPETING RISK CODE)
# ==============================
disease_categories <- list(
  "ICD10_A_B_Infectious" = list(
    filter_expr = "substr(diagnosis, 1, 1) %in% c('A', 'B')",
    description = "Infectious diseases"
  ),
  "ICD10_C_D00_D48_Neoplasms" = list(
    filter_expr = "substr(diagnosis, 1, 1) == 'C' | (substr(diagnosis, 1, 1) == 'D' & as.numeric(substr(diagnosis, 2, 3)) >= 0 & as.numeric(substr(diagnosis, 2, 3)) <= 48)",
    description = "Neoplasms"
  ),
  "ICD10_D50_D89_Blood" = list(
    filter_expr = "substr(diagnosis, 1, 1) == 'D' & as.numeric(substr(diagnosis, 2, 3)) >= 50 & as.numeric(substr(diagnosis, 2, 3)) <= 89",
    description = "Blood disorders"
  ),
  "ICD10_Chapter_E" = list(
    filter_expr = "substr(diagnosis, 1, 1) == 'E'",
    description = "Endocrine/Metabolic"
  ),
  "ICD10_Chapter_F" = list(
    filter_expr = "substr(diagnosis, 1, 1) == 'F'",
    description = "Mental/Behavioral"
  ),
  "ICD10_Chapter_G" = list(
    filter_expr = "substr(diagnosis, 1, 1) == 'G'",
    description = "Nervous System"
  ),
  "ICD10_Chapter_H" = list(
    filter_expr = "substr(diagnosis, 1, 1) == 'H'",
    description = "Eye/Ear"
  ),
  "ICD10_Chapter_I" = list(
    filter_expr = "substr(diagnosis, 1, 1) == 'I'",
    description = "Circulatory"
  ),
  "ICD10_Chapter_J" = list(
    filter_expr = "substr(diagnosis, 1, 1) == 'J'",
    description = "Respiratory"
  ),
  "ICD10_Chapter_K" = list(
    filter_expr = "substr(diagnosis, 1, 1) == 'K'",
    description = "Digestive"
  ),
  "ICD10_Chapter_L" = list(
    filter_expr = "substr(diagnosis, 1, 1) == 'L'",
    description = "Skin"
  ),
  "ICD10_Chapter_M" = list(
    filter_expr = "substr(diagnosis, 1, 1) == 'M'",
    description = "Musculoskeletal"
  ),
  "ICD10_Chapter_N" = list(
    filter_expr = "substr(diagnosis, 1, 1) == 'N'",
    description = "Genitourinary"
  )
)

# ==============================
# Function to get incident diseases
# ==============================
get_incident_diseases <- function(diag_dt, mca_info_df, filter_expr) {
  dt_diag <- copy(diag_dt)
  
  # Apply filter
  dt_diag <- dt_diag[eval(parse(text = filter_expr))]
  
  # Filter to study participants
  dt_diag <- dt_diag[ID %in% mca_info_df$ID]
  
  # Join with mCA dates
  mca_lookup <- data.table(ID = mca_info_df$ID, mCA_date = mca_info_df$mCA_date)
  setkey(mca_lookup, ID)
  setkey(dt_diag, ID)
  
  dt_diag_mca <- dt_diag[mca_lookup, on = .(ID = ID), nomatch = 0]
  
  # Pre-existing diseases (at or before mCA date)
  pre_existing <- dt_diag_mca[admission_date <= mCA_date, .(ID, diagnosis)]
  pre_existing <- unique(pre_existing)
  
  # Post-mCA diseases
  post_mca <- dt_diag_mca[admission_date > mCA_date]
  
  # Exclude pre-existing (incident only)
  setkey(post_mca, ID, diagnosis)
  setkey(pre_existing, ID, diagnosis)
  incident <- post_mca[!pre_existing]
  
  # Get first diagnosis date per patient
  result <- incident[, .(first_dx_date = min(admission_date)), by = ID]
  
  return(as.data.frame(result))
}

# ==============================
# Chromosomes (excluding homogeneous mCA chromosomes for specific CNV types)
# ==============================
# Note: Homogeneous mCAs (12+, 14+, 15+, 18+, X-) are already excluded from cnv_all
# The chromosome loop will naturally skip these since no CNV data exists for them
chromosomes <- c(as.character(1:22), "X", "Y")

cat("=== Chromosomes to analyze ===\n")
cat("All chromosomes:", paste(chromosomes, collapse = ", "), "\n")
cat("Note: Homogeneous mCAs (12+, 14+, 15+, 18+, X-) were excluded from CNV data\n\n")

# ==============================
# Main loop over disease categories and chromosomes
# ==============================
all_results <- list()

cat("=== STARTING LASSO-COX ANALYSES ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

for (cat_name in names(disease_categories)) {
  cat_info <- disease_categories[[cat_name]]
  cat("\n========================================\n")
  cat("Processing:", cat_name, "(", cat_info$description, ")\n")
  cat("========================================\n")
  
  # Get incident diseases for this category
  incident_dx <- get_incident_diseases(diag_data, mca_info, cat_info$filter_expr)
  cat("Incident cases found:", nrow(incident_dx), "\n")
  
  # Filter to incident cohort (exclude prevalent cases)
  incident_ids <- mca_info %>%
    left_join(incident_dx, by = "ID") %>%
    filter(is.na(first_dx_date) | first_dx_date > mCA_date) %>%
    pull(ID)
  
  cat("Incident cohort size:", length(incident_ids), "\n")
  
  for (chr in chromosomes) {
    cat("  Processing chromosome:", chr, "... ")
    
    # Get cytoband info for this chromosome
    ideo_chr <- subset(hg19IdeogramCyto, seqnames == paste0("chr", chr))
    
    # Get CNV data for this chromosome (only from mCA patients)
    cnv_chr <- cnv_all %>% 
      filter(CHR == chr, ID %in% mca_ids, COPY_CHANGE %in% c("gain", "loss", "neutral"))
    
    if (nrow(cnv_chr) == 0) {
      cat("no CNVs\n")
      next
    }
    
    # Create GenomicRanges for CNVs
    gr_cnv <- GRanges(
      seqnames = paste0("chr", chr),
      ranges = IRanges(start = cnv_chr$START_MB * 1e6,
                       end = cnv_chr$END_MB * 1e6),
      ID = cnv_chr$ID,
      cnv_type = cnv_chr$COPY_CHANGE
    )
    
    # Find overlaps with cytobands
    hits <- findOverlaps(gr_cnv, ideo_chr)
    band_map <- data.frame(
      ID = gr_cnv$ID[queryHits(hits)],
      cytoband = ideo_chr$name[subjectHits(hits)],
      cnv_type = gr_cnv$cnv_type[queryHits(hits)]
    ) %>% distinct()
    
    if (nrow(band_map) == 0) {
      cat("no band mappings\n")
      next
    }
    
    # Create binary feature matrix for cytobands
    band_mat <- band_map %>%
      mutate(val = 1) %>%
      unite(col = "band_type", cytoband, cnv_type, sep = "_") %>%
      pivot_wider(names_from = band_type, values_from = val, values_fill = 0)
    
    # Build analysis cohort: mCA patients with CNVs on this chr + no_mCA patients
    # All must be in incident cohort
    cnv_chr_ids <- unique(cnv_chr$ID)
    
    group_df <- tibble(ID = c(cnv_chr_ids, no_mca_ids)) %>%
      filter(ID %in% incident_ids) %>%
      mutate(group = ifelse(ID %in% cnv_chr_ids, paste0("has_chr", chr, "_CNV"), "no_mCA"))
    
    # Build survival data with covariates
    surv_df <- mca_info %>%
      filter(ID %in% group_df$ID) %>%
      left_join(group_df, by = "ID") %>%
      left_join(incident_dx, by = "ID") %>%
      mutate(
        # Event = disease diagnosis after mCA date (and before death if applicable)
        event = ifelse(!is.na(first_dx_date) & first_dx_date > mCA_date & 
                         (is.na(death_date) | first_dx_date <= death_date), 1, 0),
        # End date: disease date, death date, or censor date (whichever comes first)
        end_date = case_when(
          event == 1 ~ first_dx_date,
          !is.na(death_date) & death_date > mCA_date & death_date < censor_date ~ death_date,
          TRUE ~ censor_date
        ),
        time_years = as.numeric(difftime(end_date, mCA_date, units = "days")) / 365.25,
        sex_numeric = as.numeric(sex) - 1  # Convert sex to numeric for model
      ) %>%
      filter(time_years >= 0, !is.na(age), !is.na(sex), !is.na(n_baseline_diseases))
    
    # Filter band matrix to patients in survival data
    band_mat_filtered <- band_mat %>% filter(ID %in% surv_df$ID)
    
    # Check event thresholds (SAME AS COMPETING RISK CODE)
    total_events <- sum(surv_df$event)
    events_mca <- sum(surv_df$event == 1 & surv_df$group == paste0("has_chr", chr, "_CNV"))
    events_no_mca <- sum(surv_df$event == 1 & surv_df$group == "no_mCA")
    groups_with_events <- (events_mca > 0) + (events_no_mca > 0)
    
    if (!(nrow(surv_df) > 0 && total_events > 10 && groups_with_events == 2)) {
      cat("insufficient events (total:", total_events, ", mCA:", events_mca, ", no_mCA:", events_no_mca, ")\n")
      next
    }
    
    # Prepare matrices for LASSO
    X_bands <- as.matrix(column_to_rownames(band_mat_filtered, "ID"))
    
    # Add no_mCA patients with all zeros for band features
    no_mca_in_surv <- surv_df %>% 
      filter(group == "no_mCA", !ID %in% rownames(X_bands)) %>% 
      pull(ID)
    
    if (length(no_mca_in_surv) > 0) {
      zero_mat <- matrix(0, nrow = length(no_mca_in_surv), ncol = ncol(X_bands))
      rownames(zero_mat) <- no_mca_in_surv
      colnames(zero_mat) <- colnames(X_bands)
      X_bands <- rbind(X_bands, zero_mat)
    }
    
    # Align survival data with X matrix
    Y <- surv_df %>% 
      filter(ID %in% rownames(X_bands)) %>% 
      arrange(match(ID, rownames(X_bands)))
    
    X_bands <- X_bands[Y$ID, , drop = FALSE]
    
    # Create survival object
    surv_obj <- Surv(Y$time_years, Y$event)
    
    # Run LASSO with cross-validation
    set.seed(42)
    
    tryCatch({
      lasso_fit <- cv.glmnet(X_bands, surv_obj, family = "cox", alpha = 1)
      coefs <- coef(lasso_fit, s = "lambda.min")
      selected_bands <- rownames(coefs)[as.numeric(coefs) != 0]
      
      if (length(selected_bands) == 0) {
        cat("no bands selected by LASSO\n")
        next
      }
      
      cat(length(selected_bands), "bands selected... ")
      
      # Run adjusted Cox regression for each selected band
      cox_results <- list()
      
      for (band in selected_bands) {
        Y_band <- Y %>% mutate(marker = X_bands[, band])
        
        # Check thresholds (consistent with competing risk code)
        n_with_marker <- sum(Y_band$marker == 1)
        n_without_marker <- sum(Y_band$marker == 0)
        events_with_marker <- sum(Y_band$marker == 1 & Y_band$event == 1)
        events_without_marker <- sum(Y_band$marker == 0 & Y_band$event == 1)
        
        # Need variation in marker and events in both groups
        if (length(unique(Y_band$marker)) < 2) next
        if (events_with_marker == 0 || events_without_marker == 0) next
        if ((events_with_marker + events_without_marker) <= 10) next
        
        # Unadjusted Cox model
        fit_unadj <- coxph(Surv(time_years, event) ~ marker, data = Y_band)
        sumfit_unadj <- summary(fit_unadj)
        
        # Adjusted Cox model (age + sex + baseline disease burden)
        fit_adj <- coxph(Surv(time_years, event) ~ marker + age + sex_numeric + n_baseline_diseases, 
                         data = Y_band)
        sumfit_adj <- summary(fit_adj)
        
        cox_results[[band]] <- tibble(
          cytoband = band,
          # Unadjusted results
          HR_unadj = sumfit_unadj$conf.int["marker", "exp(coef)"],
          HR_unadj_lower = sumfit_unadj$conf.int["marker", "lower .95"],
          HR_unadj_upper = sumfit_unadj$conf.int["marker", "upper .95"],
          p_unadj = sumfit_unadj$coefficients["marker", "Pr(>|z|)"],
          # Adjusted results
          HR_adj = sumfit_adj$conf.int["marker", "exp(coef)"],
          HR_adj_lower = sumfit_adj$conf.int["marker", "lower .95"],
          HR_adj_upper = sumfit_adj$conf.int["marker", "upper .95"],
          p_adj = sumfit_adj$coefficients["marker", "Pr(>|z|)"],
          # Log2 HR
          log2HR_unadj = log2(sumfit_unadj$conf.int["marker", "exp(coef)"][[1]]),
          log2HR_adj = log2(sumfit_adj$conf.int["marker", "exp(coef)"][[1]]),
          # Sample info
          n_with_marker = n_with_marker,
          n_without_marker = n_without_marker,
          n_events_with_marker = events_with_marker,
          n_events_without_marker = events_without_marker,
          n_total = nrow(Y_band),
          n_events_total = sum(Y_band$event),
          # Metadata
          icd_category = cat_name,
          chromosome = chr
        )
      }
      
      if (length(cox_results) > 0) {
        cox_df <- bind_rows(cox_results) %>%
          mutate(
            q_unadj = p.adjust(p_unadj, method = "BH"),
            q_adj = p.adjust(p_adj, method = "BH")
          )
        all_results[[paste0(cat_name, "_chr", chr)]] <- cox_df
        cat(nrow(cox_df), "results\n")
      } else {
        cat("no valid results\n")
      }
      
    }, error = function(e) {
      cat("ERROR:", e$message, "\n")
    })
  }
}

# ==============================
# Combine and save results
# ==============================
cat("\n=== FINALIZING RESULTS ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

if (length(all_results) > 0) {
  final_df <- bind_rows(all_results)
  
  # Apply global FDR correction
  final_df <- final_df %>%
    mutate(
      q_unadj_global = p.adjust(p_unadj, method = "BH"),
      q_adj_global = p.adjust(p_adj, method = "BH")
    )
  
  # Format for output
  final_df_formatted <- final_df %>%
    mutate(
      HR_unadj_95CI = paste0(sprintf("%.3f", HR_unadj), " (",
                             sprintf("%.3f", HR_unadj_lower), "-",
                             sprintf("%.3f", HR_unadj_upper), ")"),
      HR_adj_95CI = paste0(sprintf("%.3f", HR_adj), " (",
                           sprintf("%.3f", HR_adj_lower), "-",
                           sprintf("%.3f", HR_adj_upper), ")")
    ) %>%
    arrange(icd_category, chromosome, p_adj)
  
  # Summary statistics
  cat("Total significant associations (FDR < 0.05, adjusted):", 
      sum(final_df$q_adj_global < 0.05, na.rm = TRUE), "\n")
  cat("Total significant associations (FDR < 0.10, adjusted):", 
      sum(final_df$q_adj_global < 0.10, na.rm = TRUE), "\n")
  cat("Total results:", nrow(final_df), "\n\n")
  
  # Save to Excel
  output_path <- file.path(ukbb_dir, "Cox_LASSO_Adjusted_AllChromosomes_ICD10_Categories.xlsx")
  
  wb <- createWorkbook()
  
  # Summary sheet
  addWorksheet(wb, "Summary")
  writeData(wb, "Summary", final_df_formatted %>% 
              select(icd_category, chromosome, cytoband, 
                     HR_unadj_95CI, p_unadj, q_unadj_global,
                     HR_adj_95CI, p_adj, q_adj_global,
                     n_with_marker, n_without_marker,
                     n_events_with_marker, n_events_without_marker,
                     n_total, n_events_total))
  
  # Detailed sheet
  addWorksheet(wb, "Detailed")
  writeData(wb, "Detailed", final_df)
  
  # Results by category
  for (cat_name in unique(final_df$icd_category)) {
    sheet_name <- substr(cat_name, 1, 31)  # Excel sheet name limit
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, final_df %>% filter(icd_category == cat_name))
  }
  
  saveWorkbook(wb, output_path, overwrite = TRUE)
  cat("Results saved to:", output_path, "\n")
  
} else {
  cat("No results to save.\n")
}

cat("\n🎉 ANALYSIS COMPLETE! 🎉\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
