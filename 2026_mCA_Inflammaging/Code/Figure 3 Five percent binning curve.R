# ===============================
# 5% BINNING DOSE-RESPONSE CURVES
# FOR SELECTED mCA-DISEASE COMBINATIONS
# RED COLOR SCHEME
# FIXED X-AXIS MAX AT 60%
# ===============================

message("=== GENERATING 5% BINNING CURVES FOR SELECTED COMBINATIONS ===")
message("Start time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

# Load the saved workspace
workspace_file <- "~/Desktop/UK BB/dose_response_results/spline_analysis_workspace_competing_risk_20251126.RData"
message("Loading workspace: ", workspace_file)
load(workspace_file)
message("✅ Workspace loaded successfully")

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(survival)
  library(grid)
  library(gridExtra)
})

# ===============================
# DEFINE SELECTED COMBINATIONS
# Same as spline curves
# ===============================

selected_combinations <- tribble(
  ~subtype, ~disease,
  # (1) chrY_loss combinations
  "chrY_loss", "Blood disorders",
  "chrY_loss", "Cancer",
  "chrY_loss", "Respiratory system diseases",
  "chrY_loss", "Mental and behavioral disorders",
  
  # (2) chr9_neutral combinations
  "chr9_neutral", "Blood disorders",
  "chr9_neutral", "Cancer",
  "chr9_neutral", "Circulatory system diseases",
  "chr9_neutral", "Respiratory system diseases",
  "chr9_neutral", "Infectious and parasitic diseases",
  
  # (3) chr12_gain combinations
  "chr12_gain", "Blood disorders",
  "chr12_gain", "Infectious and parasitic diseases",
  "chr12_gain", "Cancer",
  "chr12_gain", "Respiratory system diseases",
  
  # (4) chr13_loss combinations
  "chr13_loss", "Blood disorders",
  "chr13_loss", "Cancer",
  "chr13_loss", "Infectious and parasitic diseases",
  "chr13_loss", "Respiratory system diseases",
  "chr13_loss", "Skin and subcutaneous diseases"
)

message("\nSelected combinations to analyze: ", nrow(selected_combinations))
print(selected_combinations)

# ===============================
# LOOK UP FDR-ADJUSTED Q-VALUES FROM SAVED RESULTS
# ===============================

message("\n=== LOOKING UP FDR-ADJUSTED Q-VALUES ===")

# Check if successful_results exists in workspace
if (!exists("successful_results")) {
  stop("❌ successful_results not found in workspace!")
}

message("Found ", nrow(successful_results), " successful results with FDR-adjusted q-values")

# Create lookup function for q-values
get_fdr_qvalue <- function(disease_name, subtype_name, results_df) {
  match_row <- results_df %>%
    dplyr::filter(Disease == disease_name, mCA_subtype == subtype_name)
  
  if (nrow(match_row) == 1) {
    return(list(
      wald_p = match_row$spline_wald_p,
      wald_q = match_row$spline_wald_q
    ))
  } else {
    return(list(wald_p = NA, wald_q = NA))
  }
}

# ===============================
# 5% BINNING CURVE CREATION FUNCTION
# ===============================

create_5percent_binning_curve <- function(disease_name, subtype_name, disease_data) {
  
  message("\n--- Creating 5% binning curve: ", subtype_name, " → ", disease_name, " ---")
  
  tryCatch({
    df_disease <- disease_data$data
    
    # Check if subtype exists
    if (!subtype_name %in% names(df_disease)) {
      message("❌ Subtype column not found: ", subtype_name)
      return(NULL)
    }
    
    # Create subtype dataset
    df_sub <- df_disease %>%
      mutate(exposed = as.integer(.data[[subtype_name]] == 1L)) %>%
      dplyr::filter(exposed == 1L | No_mCA == 1L)
    
    # Get fraction values
    frac_col_name <- paste0(subtype_name, "_frac")
    if (frac_col_name %in% names(df_sub)) {
      df_sub <- df_sub %>%
        mutate(frac = ifelse(exposed == 1L, pmin(pmax(.data[[frac_col_name]], 0), 1), 0))
    } else {
      message("❌ Fraction column not found: ", frac_col_name)
      return(NULL)
    }
    
    # Filter to exposed individuals
    exposed_data <- df_sub %>% dplyr::filter(exposed == 1L)
    unexposed_data <- df_sub %>% dplyr::filter(exposed == 0L)
    
    message("  Exposed individuals: ", nrow(exposed_data))
    message("  Unexposed individuals: ", nrow(unexposed_data))
    message("  Events in exposed: ", sum(exposed_data$event))
    message("  Fraction range: ", round(min(exposed_data$frac, na.rm = TRUE), 4), 
            " to ", round(max(exposed_data$frac, na.rm = TRUE), 4))
    
    # Check minimum requirements
    if (nrow(exposed_data) < 50) {
      message("❌ Insufficient exposed individuals")
      return(NULL)
    }
    
    # ===============================
    # CREATE 5% BINS
    # ===============================
    
    # Create 5% bins (0-5%, 5-10%, 10-15%, ..., 95-100%)
    exposed_data <- exposed_data %>%
      mutate(
        frac_percent = frac * 100,
        bin_5pct = floor(frac_percent / 5)  # 0, 1, 2, ... 19 (each represents 5% range)
      )
    
    # Cap at 19 for 100% fractions (95-100% bin)
    exposed_data$bin_5pct[exposed_data$bin_5pct >= 20] <- 19
    
    # Summarize by 5% bins
    bin_summary <- exposed_data %>%
      group_by(bin_5pct) %>%
      summarise(
        n = n(),
        events = sum(event),
        frac_mean = mean(frac, na.rm = TRUE),
        frac_min = min(frac, na.rm = TRUE),
        frac_max = max(frac, na.rm = TRUE),
        .groups = "drop"
      )
    
    message("  Total 5% bins with data: ", nrow(bin_summary))
    message("  Bins with ≥10 individuals: ", sum(bin_summary$n >= 10))
    message("  Bins with ≥3 events: ", sum(bin_summary$events >= 3))
    
    # ===============================
    # CALCULATE HR FOR EACH BIN VS UNEXPOSED
    # ===============================
    
    # Filter to bins with sufficient data
    valid_bins <- bin_summary %>%
      dplyr::filter(n >= 10, events >= 3)
    
    if (nrow(valid_bins) < 2) {
      message("❌ Insufficient valid bins for plotting")
      return(NULL)
    }
    
    message("  Valid bins for HR calculation: ", nrow(valid_bins))
    
    bin_hrs <- list()
    
    for (i in 1:nrow(valid_bins)) {
      bin_num <- valid_bins$bin_5pct[i]
      bin_data <- exposed_data %>% dplyr::filter(bin_5pct == bin_num)
      
      # Combine with unexposed for comparison
      combined_data <- bind_rows(
        unexposed_data %>% mutate(is_bin = 0L),
        bin_data %>% mutate(is_bin = 1L)
      ) %>%
        mutate(is_bin = factor(is_bin, levels = c(0, 1)))
      
      # Fit Cox model adjusted for age, sex, baseline disease burden
      fit <- tryCatch({
        coxph(Surv(time_years, event) ~ is_bin + age + sex + n_baseline_diseases, 
              data = combined_data)
      }, error = function(e) NULL)
      
      if (!is.null(fit) && "is_bin1" %in% names(coef(fit))) {
        coef_val <- coef(fit)["is_bin1"]
        se_val <- summary(fit)$coefficients["is_bin1", "se(coef)"]
        p_val <- summary(fit)$coefficients["is_bin1", "Pr(>|z|)"]
        
        bin_hrs[[i]] <- data.frame(
          bin_5pct = bin_num,
          bin_start = bin_num * 5,        # e.g., 0, 5, 10, ...
          bin_end = (bin_num + 1) * 5,    # e.g., 5, 10, 15, ...
          bin_center = bin_num * 5 + 2.5, # Center of bin (e.g., 2.5% for 0-5% bin)
          n = valid_bins$n[i],
          events = valid_bins$events[i],
          hr = exp(coef_val),
          hr_lower = exp(coef_val - 1.96 * se_val),
          hr_upper = exp(coef_val + 1.96 * se_val),
          log_hr = coef_val,
          se = se_val,
          p_value = p_val
        )
      }
    }
    
    if (length(bin_hrs) == 0) {
      message("❌ No valid HR estimates calculated")
      return(NULL)
    }
    
    # Combine all bin results
    bin_results <- do.call(rbind, bin_hrs)
    
    message("  ✅ HR calculated for ", nrow(bin_results), " bins")
    message("  HR range: ", round(min(bin_results$hr), 2), " to ", round(max(bin_results$hr), 2))
    message("  MCF coverage: ", min(bin_results$bin_start), "% to ", max(bin_results$bin_end), "%")
    
    # Get FDR-adjusted q-value
    qvals <- get_fdr_qvalue(disease_name, subtype_name, successful_results)
    
    message("  Wald test p-value: ", format(qvals$wald_p, scientific = TRUE, digits = 3))
    message("  FDR-adjusted q-value: ", format(qvals$wald_q, scientific = TRUE, digits = 3))
    
    return(list(
      disease = disease_name,
      subtype = subtype_name,
      bin_results = bin_results,
      total_exposed = nrow(exposed_data),
      total_events = sum(exposed_data$event),
      n_bins = nrow(bin_results),
      max_mcf_covered = max(bin_results$bin_end),
      wald_p = qvals$wald_p,
      wald_q = qvals$wald_q
    ))
    
  }, error = function(e) {
    message("❌ Error creating 5% binning curve: ", e$message)
    return(NULL)
  })
}

# ===============================
# 5% BINNING CURVE PLOTTING FUNCTION
# RED COLOR SCHEME
# FIXED X-AXIS MAX AT 60%
# STRAIGHT LINES CONNECTING POINTS + SHADED CI RIBBON
# ===============================

plot_5percent_binning_curve <- function(binning_result, show_ci = TRUE) {
  
  if (is.null(binning_result)) return(NULL)
  
  bin_results <- binning_result$bin_results
  custom_red <- "#e31a1c"  # Red color
  x_max <- 60  # Fixed x-axis max at 60%
  
  # Sort by bin_center to ensure proper line connection
  bin_results <- bin_results %>% arrange(bin_center)
  
  # Format FDR-adjusted q-value for display
  wald_q <- binning_result$wald_q
  if (!is.na(wald_q)) {
    if (wald_q < 0.001) {
      q_label <- paste0("Q = ", format(wald_q, scientific = TRUE, digits = 2))
    } else {
      q_label <- paste0("Q = ", format(round(wald_q, 3), nsmall = 3))
    }
  } else {
    q_label <- "Q = NA"
  }
  
  # Create base plot
  p <- ggplot(bin_results, aes(x = bin_center, y = hr)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    
    # Add shaded ribbon for 95% CI (instead of error bars)
    {if (show_ci) geom_ribbon(aes(ymin = hr_lower, ymax = hr_upper), 
                              alpha = 0.2, fill = custom_red)} +
    
    # Add straight lines connecting HR points (instead of loess smooth)
    geom_line(color = custom_red, linewidth = 1.2) +
    
    # Add points sized by sample size
    geom_point(aes(size = n), color = custom_red, alpha = 0.8) +
    
    # Add FDR-adjusted q-value annotation in top-left corner
    annotate("text", 
             x = 1,  # Fixed position near left edge
             y = max(bin_results$hr_upper, na.rm = TRUE) * 0.95,
             label = q_label,
             hjust = 0, vjust = 1,
             size = 4, fontface = "bold", color = "black") +
    
    scale_size_continuous(name = "Sample Size", range = c(2, 6), guide = "legend") +
    
    # Fixed x-axis max at 60%
    scale_x_continuous(limits = c(0, x_max), breaks = seq(0, x_max, by = 10)) +
    
    labs(
      title = paste0("Dose-Response: ", binning_result$subtype, " → ", binning_result$disease),
      subtitle = paste0("5% binning, n=", binning_result$total_exposed, 
                        ", events=", binning_result$total_events,
                        ", bins=", binning_result$n_bins),
      x = "Mosaic Cell Fraction (%)",
      y = "Hazard Ratio",
      caption = "Adjusted for age, sex, and baseline disease burden | 95% CI shaded | Q = FDR-adjusted"
    ) +
    
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10),
      legend.position = "right",
      
      # CLEAN AXES
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_line(color = "black", linewidth = 0.5),
      axis.ticks.x = element_line(color = "black", linewidth = 0.3),
      axis.ticks.y = element_line(color = "black", linewidth = 0.3),
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.title.x = element_text(color = "black", face = "bold"),
      axis.title.y = element_text(color = "black", face = "bold"),
      
      # NO GRID LINES
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # WHITE BACKGROUND
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      
      # SQUARE ASPECT RATIO
      aspect.ratio = 1
    )
  
  return(p)
}

# Function to save individual PDF files
save_pdf <- function(plot_obj, filename, output_dir, width = 8, height = 8) {
  pdf_file <- file.path(output_dir, paste0(filename, ".pdf"))
  pdf(pdf_file, width = width, height = height)
  print(plot_obj)
  dev.off()
  return(pdf_file)
}

# ===============================
# GENERATE 5% BINNING CURVES FOR SELECTED COMBINATIONS
# ===============================

message("\n=== GENERATING 5% BINNING CURVES ===")

# Create output directory
binning_output_dir <- file.path(output_dir, "binning_curves_5pct_selected")
dir.create(binning_output_dir, showWarnings = FALSE, recursive = TRUE)
message("Output directory: ", binning_output_dir)

# Store results
all_binning_results <- list()
all_binning_plots <- list()
all_pdf_files <- list()

for (i in 1:nrow(selected_combinations)) {
  
  subtype_name <- selected_combinations$subtype[i]
  disease_name <- selected_combinations$disease[i]
  
  message("\n", rep("=", 60))
  message("Processing ", i, "/", nrow(selected_combinations), ": ", subtype_name, " → ", disease_name)
  message(rep("=", 60))
  
  # Check if disease dataset exists
  if (!disease_name %in% names(disease_datasets)) {
    message("❌ Disease dataset not found: ", disease_name)
    next
  }
  
  disease_data <- disease_datasets[[disease_name]]
  
  # Create 5% binning curve
  binning_result <- create_5percent_binning_curve(
    disease_name = disease_name,
    subtype_name = subtype_name,
    disease_data = disease_data
  )
  
  if (!is.null(binning_result)) {
    all_binning_results[[i]] <- binning_result
    
    # Create plot
    message("  🎨 Creating plot...")
    plot_obj <- plot_5percent_binning_curve(binning_result, show_ci = TRUE)
    
    if (!is.null(plot_obj)) {
      all_binning_plots[[i]] <- plot_obj
      
      # Save individual PDF
      filename <- paste0(sprintf("%02d", i), "_", 
                         gsub(" ", "_", disease_name), "_", 
                         subtype_name, "_5pct_binning")
      pdf_file <- save_pdf(plot_obj, filename, binning_output_dir, 8, 8)
      all_pdf_files[[i]] <- pdf_file
      
      message("  ✅ Plot created and saved")
      message("  📄 ", basename(pdf_file))
    }
  } else {
    message("  ❌ Failed to create 5% binning curve")
  }
}

# ===============================
# CREATE COMBINED PDFs BY mCA SUBTYPE
# ===============================

message("\n=== CREATING COMBINED PDFs BY mCA SUBTYPE ===")

# Group by subtype
subtypes_to_combine <- unique(selected_combinations$subtype)

for (subtype in subtypes_to_combine) {
  
  message("\nCreating combined PDF for: ", subtype)
  
  # Get indices for this subtype
  subtype_indices <- which(selected_combinations$subtype == subtype)
  subtype_plots <- all_binning_plots[subtype_indices]
  subtype_plots <- subtype_plots[!sapply(subtype_plots, is.null)]
  
  if (length(subtype_plots) > 0) {
    
    combined_pdf <- file.path(binning_output_dir, paste0(subtype, "_all_diseases_5pct_binning_combined.pdf"))
    
    pdf(combined_pdf, width = 12, height = 10, onefile = TRUE)
    
    # Title page
    grid.newpage()
    grid.text(paste0("5% Binning Dose-Response Curves: ", subtype), 
              x = 0.5, y = 0.6, gp = gpar(fontsize = 24, fontface = "bold", col = "black"))
    grid.text("Model-Free Approach with 95% Confidence Intervals", 
              x = 0.5, y = 0.5, gp = gpar(fontsize = 16, col = "#e31a1c"))
    grid.text("Adjusted for age, sex, and baseline disease burden", 
              x = 0.5, y = 0.45, gp = gpar(fontsize = 14, col = "gray50"))
    grid.text("Q-values are FDR-adjusted (Benjamini-Hochberg)",
              x = 0.5, y = 0.40, gp = gpar(fontsize = 12, col = "gray50"))
    grid.text(paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), 
              x = 0.5, y = 0.30, gp = gpar(fontsize = 12, col = "gray50"))
    
    # Individual plots
    for (plot_obj in subtype_plots) {
      print(plot_obj)
    }
    
    dev.off()
    
    message("  ✅ Combined PDF saved: ", basename(combined_pdf))
  }
}

# ===============================
# CREATE MASTER COMBINED PDF
# ===============================

message("\n=== CREATING MASTER COMBINED PDF ===")

# Get all successful plots
successful_plots <- all_binning_plots[!sapply(all_binning_plots, is.null)]

if (length(successful_plots) > 0) {
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  master_pdf <- file.path(binning_output_dir, paste0("ALL_5pct_binning_curves_selected_", timestamp, ".pdf"))
  
  pdf(master_pdf, width = 12, height = 10, onefile = TRUE)
  
  # Title page
  grid.newpage()
  grid.text("5% Binning Dose-Response Analysis", 
            x = 0.5, y = 0.70, gp = gpar(fontsize = 22, fontface = "bold", col = "black"))
  grid.text("Selected mCA-Disease Combinations", 
            x = 0.5, y = 0.63, gp = gpar(fontsize = 18, col = "black"))
  grid.text("Model-Free Approach with 95% Confidence Intervals", 
            x = 0.5, y = 0.55, gp = gpar(fontsize = 16, col = "#e31a1c"))
  grid.text("Adjusted for age, sex, and baseline disease burden", 
            x = 0.5, y = 0.50, gp = gpar(fontsize = 14, col = "gray50"))
  grid.text("Q-values are FDR-adjusted (Benjamini-Hochberg)",
            x = 0.5, y = 0.45, gp = gpar(fontsize = 12, col = "gray50"))
  
  # List combinations
  grid.text("Combinations analyzed:", 
            x = 0.5, y = 0.35, gp = gpar(fontsize = 14, fontface = "bold"))
  
  y_pos <- 0.30
  for (subtype in subtypes_to_combine) {
    grid.text(paste0("• ", subtype), 
              x = 0.5, y = y_pos, gp = gpar(fontsize = 12, col = "#e31a1c"))
    y_pos <- y_pos - 0.04
  }
  
  grid.text(paste0("Total curves: ", length(successful_plots)), 
            x = 0.5, y = y_pos - 0.02, gp = gpar(fontsize = 12, col = "gray50"))
  grid.text(paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), 
            x = 0.5, y = 0.08, gp = gpar(fontsize = 10, col = "gray50"))
  
  # Individual plots
  for (plot_obj in successful_plots) {
    print(plot_obj)
  }
  
  dev.off()
  
  message("✅ Master combined PDF saved: ", basename(master_pdf))
}

# ===============================
# FINAL SUMMARY
# ===============================

message("\n", rep("=", 60))
message("🎯 5% BINNING CURVE GENERATION COMPLETE")
message(rep("=", 60))

successful_count <- sum(!sapply(all_binning_results, is.null))
message("\nSuccessful curves: ", successful_count, "/", nrow(selected_combinations))

message("\n📊 Summary with FDR-adjusted q-values and MCF coverage:")
for (i in 1:length(all_binning_results)) {
  if (!is.null(all_binning_results[[i]])) {
    result <- all_binning_results[[i]]
    q_formatted <- if (!is.na(result$wald_q)) {
      if (result$wald_q < 0.001) format(result$wald_q, scientific = TRUE, digits = 2)
      else format(round(result$wald_q, 4), nsmall = 4)
    } else "NA"
    
    sig_marker <- if (!is.na(result$wald_q) && result$wald_q < 0.05) " *" else ""
    message("  ", result$subtype, " → ", result$disease, 
            ": Q = ", q_formatted, sig_marker,
            " (", result$n_bins, " bins, covers 0-", result$max_mcf_covered, "%)")
  }
}

message("\n  (* = significant at FDR < 0.05)")

message("\n📊 Summary by mCA subtype:")
for (subtype in subtypes_to_combine) {
  subtype_indices <- which(selected_combinations$subtype == subtype)
  subtype_results <- all_binning_results[subtype_indices]
  subtype_results_valid <- subtype_results[!sapply(subtype_results, is.null)]
  success_count <- length(subtype_results_valid)
  total_count <- length(subtype_indices)
  
  if (success_count > 0) {
    max_coverage <- max(sapply(subtype_results_valid, function(x) x$max_mcf_covered))
    avg_bins <- round(mean(sapply(subtype_results_valid, function(x) x$n_bins)), 1)
    message("  ", subtype, ": ", success_count, "/", total_count, " curves",
            " (avg ", avg_bins, " bins, max coverage ", max_coverage, "%)")
  } else {
    message("  ", subtype, ": ", success_count, "/", total_count, " curves")
  }
}

message("\n📁 Output directory: ", binning_output_dir)
message("\n📄 Files created:")

# List individual PDFs
individual_pdfs <- list.files(binning_output_dir, pattern = "^\\d+_.*\\.pdf$", full.names = FALSE)
if (length(individual_pdfs) > 0) {
  message("\n  Individual curve PDFs:")
  for (pdf in individual_pdfs) {
    message("    • ", pdf)
  }
}

# List combined PDFs
combined_pdfs <- list.files(binning_output_dir, pattern = "(combined|ALL_)", full.names = FALSE)
if (length(combined_pdfs) > 0) {
  message("\n  Combined PDFs:")
  for (pdf in combined_pdfs) {
    message("    • ", pdf)
  }
}

message("\n✅ Analysis complete!")
message("End time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
