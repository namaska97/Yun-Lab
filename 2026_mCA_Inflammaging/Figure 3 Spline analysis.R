# ===============================
# SPLINE CURVE GENERATION FOR SELECTED mCA-DISEASE COMBINATIONS
# WITH FDR-ADJUSTED Q-VALUES DISPLAYED
# Uses saved workspace from successful spline analysis
# ===============================

message("=== GENERATING SPLINE CURVES FOR SELECTED COMBINATIONS ===")
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
  library(splines)
  library(grid)
  library(gridExtra)
})

# ===============================
# DEFINE SELECTED COMBINATIONS
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

# Preview q-values for selected combinations
message("\nQ-values for selected combinations:")
for (i in 1:nrow(selected_combinations)) {
  qvals <- get_fdr_qvalue(selected_combinations$disease[i], 
                          selected_combinations$subtype[i], 
                          successful_results)
  if (!is.na(qvals$wald_q)) {
    message("  ", selected_combinations$subtype[i], " → ", selected_combinations$disease[i], 
            ": Q = ", format(qvals$wald_q, scientific = TRUE, digits = 2))
  }
}

# ===============================
# SPLINE CURVE CREATION FUNCTION
# ===============================

create_spline_curve <- function(disease_name, subtype_name, disease_data, n_points = 100) {
  
  message("\n--- Creating spline curve: ", subtype_name, " → ", disease_name, " ---")
  
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
    
    # Check fraction variation
    frac_range <- diff(range(exposed_data$frac, na.rm = TRUE))
    if (frac_range < 0.01) {
      message("❌ Insufficient fraction variation")
      return(NULL)
    }
    
    # Create spline basis functions
    spline_df <- 3
    
    if (length(unique(exposed_data$frac)) < spline_df + 1) {
      message("❌ Insufficient unique fraction values for spline")
      return(NULL)
    }
    
    # Create spline basis matrix for exposed individuals
    spline_basis <- ns(exposed_data$frac, df = spline_df)
    
    # Create full dataset with spline terms
    spline_data <- df_sub
    
    # Initialize spline columns
    for (i in 1:spline_df) {
      spline_data[[paste0("spline_", i)]] <- 0
    }
    
    # Fill in spline values for exposed individuals
    exposed_indices <- which(spline_data$exposed == 1L)
    for (i in 1:spline_df) {
      spline_data[[paste0("spline_", i)]][exposed_indices] <- spline_basis[, i]
    }
    
    # Fit the spline model (adjusted for age, sex, and baseline disease burden)
    spline_formula <- as.formula(paste0(
      "Surv(time_years, event) ~ exposed + ",
      paste0("spline_", 1:spline_df, collapse = " + "),
      " + age + sex + n_baseline_diseases"
    ))
    
    spline_fit <- coxph(spline_formula, data = spline_data)
    
    message("  ✅ Spline model fitted successfully")
    
    # ===============================
    # GET FDR-ADJUSTED Q-VALUE FROM SAVED RESULTS
    # ===============================
    qvals <- get_fdr_qvalue(disease_name, subtype_name, successful_results)
    wald_p <- qvals$wald_p
    wald_q <- qvals$wald_q
    
    message("  Wald test p-value: ", format(wald_p, scientific = TRUE, digits = 3))
    message("  FDR-adjusted q-value: ", format(wald_q, scientific = TRUE, digits = 3))
    
    # Create prediction data points across the fraction range
    frac_min <- min(exposed_data$frac, na.rm = TRUE)
    frac_max <- max(exposed_data$frac, na.rm = TRUE)
    pred_fractions <- seq(frac_min, frac_max, length.out = n_points)
    
    # Create spline basis for prediction points
    pred_spline_basis <- predict(ns(exposed_data$frac, df = spline_df), pred_fractions)
    
    # Use representative values for covariates
    mean_age <- mean(spline_data$age, na.rm = TRUE)
    ref_sex <- names(sort(table(spline_data$sex), decreasing = TRUE))[1]
    mean_baseline <- mean(spline_data$n_baseline_diseases, na.rm = TRUE)
    
    # Create prediction dataset
    pred_data <- data.frame(
      exposed = rep(1L, n_points),
      age = rep(mean_age, n_points),
      sex = factor(rep(ref_sex, n_points), levels = levels(spline_data$sex)),
      n_baseline_diseases = rep(mean_baseline, n_points)
    )
    
    # Add spline terms
    for (i in 1:spline_df) {
      pred_data[[paste0("spline_", i)]] <- pred_spline_basis[, i]
    }
    
    # Predict log hazard ratios
    pred_log_hr <- predict(spline_fit, newdata = pred_data, type = "lp", se.fit = TRUE)
    
    # Calculate baseline (unexposed) prediction
    baseline_data <- data.frame(
      exposed = 0L,
      age = mean_age,
      sex = factor(ref_sex, levels = levels(spline_data$sex)),
      n_baseline_diseases = mean_baseline
    )
    
    # Add zero spline terms for baseline
    for (i in 1:spline_df) {
      baseline_data[[paste0("spline_", i)]] <- 0
    }
    
    baseline_log_hr <- predict(spline_fit, newdata = baseline_data, type = "lp")
    
    # Calculate relative log HR (exposed vs unexposed)
    relative_log_hr <- pred_log_hr$fit - as.numeric(baseline_log_hr)
    relative_se <- pred_log_hr$se.fit
    
    # Convert to HR scale
    pred_hr <- exp(relative_log_hr)
    pred_hr_lower <- exp(relative_log_hr - 1.96 * relative_se)
    pred_hr_upper <- exp(relative_log_hr + 1.96 * relative_se)
    
    # Create results dataframe
    curve_data <- data.frame(
      fraction = pred_fractions,
      fraction_percent = pred_fractions * 100,
      hr = pred_hr,
      hr_lower = pred_hr_lower,
      hr_upper = pred_hr_upper
    )
    
    message("  ✅ Spline curve calculated")
    message("  HR range: ", round(min(pred_hr), 2), " to ", round(max(pred_hr), 2))
    
    # Create observed data points for overlay
    n_groups <- 8
    exposed_data$frac_group <- cut(exposed_data$frac, 
                                   breaks = n_groups, 
                                   include.lowest = TRUE, 
                                   labels = FALSE)
    
    # Calculate observed stats for each group
    observed_points <- exposed_data %>%
      dplyr::filter(!is.na(frac_group)) %>%
      group_by(frac_group) %>%
      summarise(
        n = n(),
        events = sum(event),
        frac_mean = mean(frac, na.rm = TRUE),
        frac_median = median(frac, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::filter(n >= 15, events >= 3)
    
    message("  Groups for observed points: ", nrow(observed_points))
    
    # Calculate HR for each group vs unexposed
    if (nrow(observed_points) > 0) {
      group_hrs <- list()
      
      for (i in 1:nrow(observed_points)) {
        group_num <- observed_points$frac_group[i]
        group_data <- exposed_data %>% dplyr::filter(frac_group == group_num)
        
        combined_data <- bind_rows(
          unexposed_data %>% mutate(is_group = 0L),
          group_data %>% mutate(is_group = 1L)
        ) %>%
          mutate(is_group = factor(is_group, levels = c(0, 1)))
        
        fit <- tryCatch({
          coxph(Surv(time_years, event) ~ is_group + age + sex + n_baseline_diseases, 
                data = combined_data)
        }, error = function(e) NULL)
        
        if (!is.null(fit) && "is_group1" %in% names(coef(fit))) {
          coef_val <- coef(fit)["is_group1"]
          se_val <- summary(fit)$coefficients["is_group1", "se(coef)"]
          
          group_hrs[[i]] <- data.frame(
            frac_group = group_num,
            hr_obs = exp(coef_val),
            hr_obs_lower = exp(coef_val - 1.96 * se_val),
            hr_obs_upper = exp(coef_val + 1.96 * se_val)
          )
        }
      }
      
      if (length(group_hrs) > 0) {
        group_hr_df <- do.call(rbind, group_hrs)
        observed_points <- observed_points %>%
          left_join(group_hr_df, by = "frac_group") %>%
          dplyr::filter(!is.na(hr_obs))
        
        message("  Observed points with valid HRs: ", nrow(observed_points))
      } else {
        observed_points <- NULL
      }
    } else {
      observed_points <- NULL
    }
    
    # Return results including FDR-adjusted q-value
    return(list(
      disease = disease_name,
      subtype = subtype_name,
      curve_data = curve_data,
      observed_points = observed_points,
      total_exposed = nrow(exposed_data),
      total_events = sum(exposed_data$event),
      model_fit = spline_fit,
      wald_p = wald_p,
      wald_q = wald_q
    ))
    
  }, error = function(e) {
    message("❌ Error creating spline curve: ", e$message)
    return(NULL)
  })
}

# ===============================
# SPLINE CURVE PLOTTING FUNCTION
# Now displays FDR-adjusted q-value
# ===============================

plot_spline_curve <- function(spline_result, show_ci = TRUE, show_points = TRUE) {
  
  if (is.null(spline_result)) return(NULL)
  
  curve_data <- spline_result$curve_data
  observed_points <- spline_result$observed_points
  custom_blue <- "#1f78b4"
  
  # Format FDR-adjusted q-value for display
  wald_q <- spline_result$wald_q
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
  p <- ggplot(curve_data, aes(x = fraction_percent, y = hr)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    
    # Add confidence band
    {if (show_ci) geom_ribbon(aes(ymin = hr_lower, ymax = hr_upper), 
                              alpha = 0.2, fill = custom_blue)} +
    
    # Add main spline curve
    geom_line(color = custom_blue, linewidth = 1.2) +
    
    # Add observed data points if available
    {if (show_points && !is.null(observed_points) && nrow(observed_points) > 0) 
      geom_point(data = observed_points, 
                 aes(x = frac_median * 100, y = hr_obs, size = n, alpha = pmin(events, 15)/15),
                 color = custom_blue)} +
    
    # Add FDR-adjusted q-value annotation in top-left corner (fixed position for 0-60% x-axis)
    annotate("text", 
             x = 1,  # Fixed position near left edge
             y = max(curve_data$hr_upper, na.rm = TRUE) * 0.95,
             label = q_label,
             hjust = 0, vjust = 1,
             size = 4, fontface = "bold", color = "black") +
    
    scale_size_continuous(name = "Sample Size", range = c(1, 4), guide = "legend") +
    scale_alpha_continuous(name = "Events", range = c(0.4, 0.9), guide = "none") +
    
    # Fixed x-axis max at 60%
    scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, by = 10)) +
    
    labs(
      title = paste0("Dose-Response: ", spline_result$subtype, " → ", spline_result$disease),
      subtitle = paste0("Natural splines (df=3), n=", spline_result$total_exposed, 
                        ", events=", spline_result$total_events),
      x = "Mosaic Cell Fraction (%)",
      y = "Hazard Ratio",
      caption = "Adjusted for age, sex, and baseline disease burden | 95% CI shown | Q = FDR-adjusted"
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
# GENERATE CURVES FOR SELECTED COMBINATIONS
# ===============================

message("\n=== GENERATING SPLINE CURVES ===")

# Create output directory
spline_output_dir <- file.path(output_dir, "spline_curves_selected")
dir.create(spline_output_dir, showWarnings = FALSE, recursive = TRUE)
message("Output directory: ", spline_output_dir)

# Store results organized by subtype
all_spline_results <- list()
all_spline_plots <- list()
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
  
  # Create spline curve
  spline_result <- create_spline_curve(
    disease_name = disease_name,
    subtype_name = subtype_name,
    disease_data = disease_data,
    n_points = 100
  )
  
  if (!is.null(spline_result)) {
    all_spline_results[[i]] <- spline_result
    
    # Create plot
    message("  🎨 Creating plot...")
    plot_obj <- plot_spline_curve(spline_result, show_ci = TRUE, show_points = TRUE)
    
    if (!is.null(plot_obj)) {
      all_spline_plots[[i]] <- plot_obj
      
      # Save individual PDF
      filename <- paste0(sprintf("%02d", i), "_", 
                         gsub(" ", "_", disease_name), "_", 
                         subtype_name, "_spline")
      pdf_file <- save_pdf(plot_obj, filename, spline_output_dir, 8, 8)
      all_pdf_files[[i]] <- pdf_file
      
      message("  ✅ Plot created and saved")
      message("  📄 ", basename(pdf_file))
    }
  } else {
    message("  ❌ Failed to create spline curve")
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
  subtype_plots <- all_spline_plots[subtype_indices]
  subtype_plots <- subtype_plots[!sapply(subtype_plots, is.null)]
  
  if (length(subtype_plots) > 0) {
    
    combined_pdf <- file.path(spline_output_dir, paste0(subtype, "_all_diseases_combined.pdf"))
    
    pdf(combined_pdf, width = 12, height = 10, onefile = TRUE)
    
    # Title page
    grid.newpage()
    grid.text(paste0("Dose-Response Curves: ", subtype), 
              x = 0.5, y = 0.6, gp = gpar(fontsize = 24, fontface = "bold", col = "black"))
    grid.text("Natural Splines (df=3) with 95% Confidence Bands", 
              x = 0.5, y = 0.5, gp = gpar(fontsize = 16, col = "#1f78b4"))
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
successful_plots <- all_spline_plots[!sapply(all_spline_plots, is.null)]

if (length(successful_plots) > 0) {
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  master_pdf <- file.path(spline_output_dir, paste0("ALL_spline_curves_selected_", timestamp, ".pdf"))
  
  pdf(master_pdf, width = 12, height = 10, onefile = TRUE)
  
  # Title page
  grid.newpage()
  grid.text("Dose-Response Analysis: Selected mCA-Disease Combinations", 
            x = 0.5, y = 0.7, gp = gpar(fontsize = 22, fontface = "bold", col = "black"))
  grid.text("Natural Splines (df=3) with 95% Confidence Bands", 
            x = 0.5, y = 0.6, gp = gpar(fontsize = 16, col = "#1f78b4"))
  grid.text("Adjusted for age, sex, and baseline disease burden", 
            x = 0.5, y = 0.55, gp = gpar(fontsize = 14, col = "gray50"))
  grid.text("Q-values are FDR-adjusted (Benjamini-Hochberg)",
            x = 0.5, y = 0.50, gp = gpar(fontsize = 12, col = "gray50"))
  
  # List combinations
  grid.text("Combinations analyzed:", 
            x = 0.5, y = 0.40, gp = gpar(fontsize = 14, fontface = "bold"))
  
  y_pos <- 0.35
  for (subtype in subtypes_to_combine) {
    grid.text(paste0("• ", subtype), 
              x = 0.5, y = y_pos, gp = gpar(fontsize = 12, col = "#1f78b4"))
    y_pos <- y_pos - 0.04
  }
  
  grid.text(paste0("Total curves: ", length(successful_plots)), 
            x = 0.5, y = y_pos - 0.02, gp = gpar(fontsize = 12, col = "gray50"))
  grid.text(paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), 
            x = 0.5, y = 0.1, gp = gpar(fontsize = 10, col = "gray50"))
  
  # Individual plots
  for (plot_obj in successful_plots) {
    print(plot_obj)
  }
  
  dev.off()
  
  message("✅ Master combined PDF saved: ", basename(master_pdf))
}

# ===============================
# FINAL SUMMARY WITH FDR Q-VALUES
# ===============================

message("\n", rep("=", 60))
message("🎯 SPLINE CURVE GENERATION COMPLETE")
message(rep("=", 60))

successful_count <- sum(!sapply(all_spline_results, is.null))
message("\nSuccessful curves: ", successful_count, "/", nrow(selected_combinations))

message("\n📊 Summary with FDR-adjusted q-values:")
for (i in 1:length(all_spline_results)) {
  if (!is.null(all_spline_results[[i]])) {
    result <- all_spline_results[[i]]
    q_formatted <- if (!is.na(result$wald_q)) {
      if (result$wald_q < 0.001) format(result$wald_q, scientific = TRUE, digits = 2)
      else format(round(result$wald_q, 4), nsmall = 4)
    } else "NA"
    
    sig_marker <- if (!is.na(result$wald_q) && result$wald_q < 0.05) " *" else ""
    message("  ", result$subtype, " → ", result$disease, ": Q = ", q_formatted, sig_marker)
  }
}

message("\n  (* = significant at FDR < 0.05)")

message("\n📊 Summary by mCA subtype:")
for (subtype in subtypes_to_combine) {
  subtype_indices <- which(selected_combinations$subtype == subtype)
  subtype_results <- all_spline_results[subtype_indices]
  success_count <- sum(!sapply(subtype_results, is.null))
  total_count <- length(subtype_indices)
  message("  ", subtype, ": ", success_count, "/", total_count, " curves created")
}

message("\n📁 Output directory: ", spline_output_dir)
message("\n📄 Files created:")

# List individual PDFs
individual_pdfs <- list.files(spline_output_dir, pattern = "^\\d+_.*\\.pdf$", full.names = FALSE)
if (length(individual_pdfs) > 0) {
  message("\n  Individual curve PDFs:")
  for (pdf in individual_pdfs) {
    message("    • ", pdf)
  }
}

# List combined PDFs
combined_pdfs <- list.files(spline_output_dir, pattern = "(combined|ALL_)", full.names = FALSE)
if (length(combined_pdfs) > 0) {
  message("\n  Combined PDFs:")
  for (pdf in combined_pdfs) {
    message("    • ", pdf)
  }
}

message("\n✅ Analysis complete!")
message("End time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
