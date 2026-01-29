# ============================
# Clear workspace
# ============================
rm(list = ls())

# ============================
# Load required libraries
# ============================
library(readxl)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(circlize)
library(lubridate)

# ============================
# Config
# ============================
ukbb_dir <- "~/Desktop/UK BB/"

# ============================
# STEP 1: Apply EXACT filtering from competing risk analysis
# ============================

# Load analyzed samples
analyzed <- read_excel(file.path(ukbb_dir, "ukb4777.samples_analyzed.xlsx")) %>%
  mutate(ID = as.character(ID))

cat("\n=== Step 1: Load analyzed samples ===\n")
cat("Analyzed samples:", nrow(analyzed), "\n")

# Load mCA info with date filtering
mca_info_raw <- read_excel(file.path(ukbb_dir, "Age sex dob collection date.xlsx")) %>%
  mutate(ID = as.character(ID),
         mCA_date = mdy(`Date of mCA assessment`)) %>%
  filter(ID %in% analyzed$ID)

cat("\n=== Step 2: Date filtering ===\n")
cat("Total rows before filtering:", nrow(mca_info_raw), "\n")
cat("Missing mCA dates:", sum(is.na(mca_info_raw$mCA_date)), "\n")
cat("Dates = 1900-01-01:", sum(mca_info_raw$mCA_date == as.Date("1900-01-01"), na.rm = TRUE), "\n")

# Filter out invalid dates
mca_info <- mca_info_raw %>%
  filter(!is.na(mCA_date),
         mCA_date != as.Date("1900-01-01"),
         mCA_date > as.Date("1901-01-01")) %>%
  select(ID, mCA_date) %>%
  distinct(ID, .keep_all = TRUE)

cat("After date filtering:", nrow(mca_info), "\n")

# Load CNV data
cna_df <- read_excel(file.path(ukbb_dir, "CNV_data_19808.xlsx")) %>%
  mutate(ID = as.character(ID)) %>%
  filter(ID %in% analyzed$ID)

cat("\n=== Step 3: CNV data ===\n")
cat("Total CNV records:", nrow(cna_df), "\n")

# Clean column names
colnames(cna_df) <- trimws(colnames(cna_df))

# Summarize CNV status
cnv_summary_status <- cna_df %>%
  group_by(ID) %>%
  summarise(
    has_known_cnv = any(COPY_CHANGE %in% c("gain", "loss", "neutral")),
    has_only_unknown = all(COPY_CHANGE == "unknown"),
    .groups = "drop"
  )

# Determine mCA status
cnv_status <- analyzed %>%
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
  filter(mCA_status %in% c("mCA", "no_mCA")) %>%
  select(ID, mCA_status)

cat("\n=== Step 4: mCA status ===\n")
print(table(cnv_status$mCA_status))

# Create final filtered cohort
mca_info_filtered <- mca_info %>%
  inner_join(cnv_status, by = "ID")

cat("\n=== Step 5: Final filtered cohort ===\n")
cat("Total patients:", nrow(mca_info_filtered), "\n")
cat("Distribution:\n")
print(table(mca_info_filtered$mCA_status))

# Verify we got 452,594 patients
if (nrow(mca_info_filtered) == 452594) {
  cat("\n✅ SUCCESS: Got exactly 452,594 patients!\n")
} else {
  cat("\n⚠️  WARNING: Expected 452,594 but got", nrow(mca_info_filtered), "\n")
}

# ============================
# STEP 2: Filter CNV data to mCA patients only
# ============================

# Get IDs of mCA patients in filtered cohort
mca_patient_ids <- mca_info_filtered %>%
  filter(mCA_status == "mCA") %>%
  pull(ID)

cat("\n=== Step 6: Filtering CNV data for heatmap ===\n")
cat("mCA patients in filtered cohort:", length(mca_patient_ids), "\n")

# Filter CNA data to only filtered mCA patients
cna_df_filtered <- cna_df %>%
  filter(ID %in% mca_patient_ids)

cat("Total CNV records for heatmap:", nrow(cna_df_filtered), "\n")

# ============================
# STEP 3: Prepare data for heatmap
# ============================

# Remove rows with "unknown" CNA status and prepare for heatmap
cna_df_clean <- cna_df_filtered %>%
  filter(tolower(COPY_CHANGE) != "unknown") %>%
  mutate(
    COPY_CHANGE = tolower(COPY_CHANGE),
    CHR = ifelse(CHR %in% c("23", "24"), c("X", "Y")[as.integer(CHR) - 22], as.character(CHR)),
    CHR = paste0("chr", CHR)
  )

cat("CNV records after removing 'unknown':", nrow(cna_df_clean), "\n")
cat("Unique patients in heatmap:", n_distinct(cna_df_clean$ID), "\n")

# Define CNA priority mapping
priority_map <- c("gain" = 3, "loss" = 2, "neutral" = 1)
cna_df_clean <- cna_df_clean %>%
  mutate(priority = priority_map[COPY_CHANGE])

# Select highest-priority CNA per chromosome–patient pair
cna_summary <- cna_df_clean %>%
  group_by(CHR, ID) %>%
  slice_max(order_by = priority, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(CHR, ID, CNA = COPY_CHANGE)

# Ensure chromosome order
chr_levels <- paste0("chr", c(1:22, "X", "Y"))
cna_summary$CHR <- factor(cna_summary$CHR, levels = chr_levels)

# Generate all combinations (to fill missing with "none")
all_combinations <- expand.grid(
  CHR = chr_levels,
  ID = unique(cna_summary$ID),
  stringsAsFactors = FALSE
)

# Merge and fill missing with "none"
cna_summary_full <- left_join(all_combinations, cna_summary, by = c("CHR", "ID")) %>%
  mutate(CNA = ifelse(is.na(CNA), "none", CNA))

# Pivot to wide matrix
cna_matrix_wide <- cna_summary_full %>%
  pivot_wider(names_from = ID, values_from = CNA)

# Convert to matrix
chr_labels <- cna_matrix_wide$CHR
cna_matrix <- as.matrix(cna_matrix_wide[, -1])
rownames(cna_matrix) <- chr_labels

# Remove patients with only "none"
cna_matrix <- cna_matrix[, colSums(cna_matrix != "none") > 0]

cat("\n=== Heatmap matrix dimensions ===\n")
cat("Chromosomes (rows):", nrow(cna_matrix), "\n")
cat("Patients (columns):", ncol(cna_matrix), "\n")

# ============================
# STEP 4: Custom sorting
# ============================

# 1. Order chromosomes (rows) by total CNAs (non-"none") in descending order
chromosome_order <- order(rowSums(cna_matrix != "none"), decreasing = TRUE)
cna_matrix <- cna_matrix[chromosome_order, ]

# 2. Order patients (columns) by aggregate CNA priority (row-wise sort)
# Map CNA values to numeric priorities
cna_priority_matrix <- matrix(priority_map[as.vector(cna_matrix)], 
                              nrow = nrow(cna_matrix), 
                              dimnames = dimnames(cna_matrix))

# Replace "none" (NA) with 0 for proper sorting
cna_priority_matrix[is.na(cna_priority_matrix)] <- 0

# Custom column sort: patients with gain > loss > neutral (by row)
patient_order <- do.call(order, as.data.frame(-t(cna_priority_matrix)))
cna_matrix <- cna_matrix[, patient_order]

# ============================
# STEP 5: Prepare MCF data for annotation
# ============================

# Subset MCF data for annotation
mcf_data <- cna_df_filtered %>%
  filter(!is.na(CELL_FRAC)) %>%
  mutate(
    CHR = ifelse(CHR %in% c("23", "24"), c("X", "Y")[as.integer(CHR) - 22], as.character(CHR)),
    CHR = paste0("chr", CHR)
  ) %>%
  filter(CHR %in% rownames(cna_matrix)) %>%
  filter(ID %in% colnames(cna_matrix)) %>%
  mutate(CHR = factor(CHR, levels = rownames(cna_matrix)))

cat("\n=== MCF annotation data ===\n")
cat("Records with cell fraction data:", nrow(mcf_data), "\n")

# Prepare numeric list of MCF values by chromosome
violin_values <- lapply(levels(mcf_data$CHR), function(chr) {
  vals <- mcf_data$CELL_FRAC[mcf_data$CHR == chr]
  if(length(vals) == 0 || all(is.na(vals))) vals <- 0
  vals
})
names(violin_values) <- levels(mcf_data$CHR)

# ============================
# STEP 6: Plot setup
# ============================

# Define CNA colors
cna_colors <- c(
  "gain" = "firebrick",
  "loss" = "steelblue",
  "neutral" = "orange",
  "none" = "gray90"
)

# Top annotation: CNA count per patient
top_anno <- HeatmapAnnotation(
  CNA_Count = anno_barplot(
    colSums(cna_matrix != "none"),
    border = FALSE,
    gp = gpar(fill = "skyblue", col = NA),
    height = unit(3, "cm")
  ),
  annotation_name_side = "left"
)

# Right annotation: Combine patient count barplot and MCF boxplot per chromosome
right_anno <- rowAnnotation(
  Patient_Count = anno_barplot(
    rowSums(cna_matrix != "none"),
    border = FALSE,
    gp = gpar(fill = "forestgreen", col = NA),
    width = unit(3, "cm"),
    axis_param = list(side = "top"),
    bar_width = 0.75
  ),
  MCF_Distribution = anno_boxplot(
    violin_values,
    axis_param = list(at = seq(0, 1, 0.2), labels = seq(0, 1, 0.2), side = "bottom"),
    border = TRUE,
    gp = gpar(fill = "lightblue", col = "black"),
    width = unit(6, "cm"),
    ylim = c(0, 1)
  ),
  annotation_name_side = "top",
  gap = unit(5, "mm")
)

# ============================
# STEP 7: Save to PDF
# ============================

output_file <- file.path(ukbb_dir, "cna_oncoplot_grouped_sorted_19808_with_MCF_filtered.pdf")

pdf(output_file, width = 25, height = 12)

Heatmap(
  cna_matrix,
  name = NULL,
  col = cna_colors,
  na_col = "gray90",
  top_annotation = top_anno,
  right_annotation = right_anno,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 45,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_heatmap_legend = FALSE
)

dev.off()

cat("\n✅ Complex heatmap saved to:", output_file, "\n")
cat("   Based on", length(mca_patient_ids), "mCA patients from filtered cohort of 452,594\n")
cat("   Matrix dimensions:", nrow(cna_matrix), "chromosomes ×", ncol(cna_matrix), "patients\n")
