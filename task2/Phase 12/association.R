# ==============================================================================
# Phase 12: Metabolite-Diabetes Phenotype Association
# ==============================================================================

# Install packages if missing: install.packages(c("dplyr", "ggplot2", "tidyr"))
library(dplyr)
library(ggplot2)
library(tidyr)

# ------------------------------------------------------------------------------
# 0. Load and Merge Data
# ------------------------------------------------------------------------------
cat("Loading datasets...\n")
metabo_data <- read.csv("../Phase 11/outputs/Qatari_metabolomics_Cleaned_Scaled.csv", check.names = FALSE)

# Attempt to load covariates from Task 1. 
# (Assuming your covariates file has columns like IID, Sex, PC1, PC2, etc.)

covariates <- read.table("../../task1/Phase 4/outputs/covariates.txt", header = TRUE, stringsAsFactors = FALSE)
# Merge on ID (main_id in metabo matches IID in covariates)
merged_data <- merge(metabo_data, covariates, by.x = "main_id", by.y = "IID", all.x = TRUE)
cat("Covariates successfully merged.\n")

# Ensure Diabetes is a factor for plotting
merged_data$Diabetes_Label <- factor(merged_data$Diabetes, levels = c(0, 1), labels = c("Non-Diabetic", "Diabetic"))

# Get the list of the 136 metabolites (Everything after Diabetes, before covariates)
# We know the clean data had main_id, Diabetes, then 136 features.
metabolite_cols <- colnames(metabo_data)[3:ncol(metabo_data)]

# ------------------------------------------------------------------------------
# 1. Exploratory Summary (Model-Free First Look)
# ------------------------------------------------------------------------------
cat("Calculating exploratory statistics...\n")
# Calculate Mean and SD for every metabolite grouped by Diabetes status
summary_stats <- merged_data %>%
  select(Diabetes_Label, all_of(metabolite_cols)) %>%
  pivot_longer(cols = -Diabetes_Label, names_to = "Metabolite", values_to = "Level") %>%
  group_by(Metabolite, Diabetes_Label) %>%
  summarise(Mean = mean(Level, na.rm = TRUE),
            SD = sd(Level, na.rm = TRUE),
            .groups = 'drop')

write.csv(summary_stats, "outputs/Deliverable_Metabolite_Means_by_Diabetes.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 2 & 3. Formal Association Test & Bonferroni Correction
# ------------------------------------------------------------------------------
cat("Running Multivariable Regressions for 136 metabolites...\n")

results_list <- list()

for (met in metabolite_cols) {
  # Build the regression formula dynamically based on available covariates
  # Wrapping metabolite in backticks handles any weird characters in chemical names
  if("PC1" %in% colnames(merged_data)) {
    form_str <- paste0("`", met, "` ~ Diabetes + Sex + PC1 + PC2 + PC3 + PC4 + PC5")
  } else {
    form_str <- paste0("`", met, "` ~ Diabetes")
  }
  
  fit <- lm(as.formula(form_str), data = merged_data)
  coefs <- summary(fit)$coefficients
  
  # Extract the row corresponding to the 'Diabetes' effect
  if("Diabetes" %in% rownames(coefs)) {
    results_list[[met]] <- data.frame(
      Metabolite = met,
      Beta = coefs["Diabetes", "Estimate"],
      SE = coefs["Diabetes", "Std. Error"],
      P_value = coefs["Diabetes", "Pr(>|t|)"]
    )
  }
}

# Combine all results into one master table
assoc_results <- bind_rows(results_list)

# Apply Bonferroni Correction
num_tests <- nrow(assoc_results)
bonferroni_thresh <- 0.05 / num_tests

assoc_results <- assoc_results %>%
  mutate(
    Bonferroni_Significant = ifelse(P_value < bonferroni_thresh, "Y", "N"),
    Neg_Log10_P = -log10(P_value)
  )

# Sort by most significant
assoc_results <- assoc_results[order(assoc_results$P_value), ]

# Export Full Association Table
write.csv(assoc_results, "outputs/Deliverable_Full_Association_Table.csv", row.names = FALSE)

# Export ONLY the significant ones (This goes to Phase 13)
significant_metabolites <- assoc_results %>% filter(Bonferroni_Significant == "Y") %>% pull(Metabolite)
write.table(significant_metabolites, "outputs/Deliverable_Significant_Metabolites_Phase13.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# ------------------------------------------------------------------------------
# 4. Data Visualization: Volcano Plot & Top Bar Plot
# ------------------------------------------------------------------------------
cat("Generating plots...\n")

# A. Volcano Plot
volcano_plot <- ggplot(assoc_results, aes(x = Beta, y = Neg_Log10_P, color = Bonferroni_Significant)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_hline(yintercept = -log10(bonferroni_thresh), linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", alpha = 0.5) +
  scale_color_manual(values = c("N" = "darkgray", "Y" = "dodgerblue")) +
  theme_minimal() +
  labs(title = "Volcano Plot: Metabolite vs Diabetes Status",
       subtitle = paste("Red line = Bonferroni threshold ( p <", signif(bonferroni_thresh, 3), ")"),
       x = "Effect Size (Beta)",
       y = "-log10(p-value)") +
  theme(legend.position = "bottom")

ggsave("outputs/Deliverable_Volcano_Plot.png", plot = volcano_plot, width = 8, height = 6, dpi = 300)

# B. Exploratory Bar Plot (Plotting ONLY the Top 10 to keep it readable)
# (Plotting 136 bars on a single chart is an anti-pattern in data science)
top_10_metabolites <- head(assoc_results$Metabolite, 10)
top_10_data <- summary_stats %>% filter(Metabolite %in% top_10_metabolites)

bar_plot <- ggplot(top_10_data, aes(x = Metabolite, y = Mean, fill = Diabetes_Label)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), 
                position = position_dodge(0.8), width = 0.25, alpha = 0.6) +
  coord_flip() + # Flip coordinates for readability of chemical names
  theme_minimal() +
  labs(title = "Top 10 Most Associated Metabolites",
       subtitle = "Mean ± Standard Deviation (Z-score scaled)",
       x = "", y = "Normalized Abundance") +
  theme(legend.title = element_blank(), legend.position = "top")

ggsave("Deliverable_Exploratory_BarPlot.png", plot = bar_plot, width = 8, height = 6, dpi = 300)

# ------------------------------------------------------------------------------
# 5. Summary Printout
# ------------------------------------------------------------------------------
cat("\n=== PHASE 12 ASSOCIATION SUMMARY ===\n")
cat("Total Metabolites Tested:    ", num_tests, "\n")
cat("Bonferroni Threshold:        ", signif(bonferroni_thresh, 3), "\n")
cat("Significant Hits Found:      ", length(significant_metabolites), "\n")
if(length(significant_metabolites) > 0) {
  cat("Top Hit:                     ", assoc_results$Metabolite[1], "(p =", signif(assoc_results$P_value[1], 3), ")\n")
}
cat("====================================\n")