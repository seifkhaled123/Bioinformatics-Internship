# ==============================================================================
# Metabolite PCA and Visualization
# ==============================================================================

# Install ggplot2 if you don't have it: install.packages("ggplot2")
library(ggplot2)

cat("Loading cleaned metabolite data...\n")
clean_data <- read.csv("outputs/Qatari_metabolomics_Cleaned_Scaled.csv", check.names = FALSE)

# 1. Isolate the Metadata and the Math Matrix
# From Phase 11, we know column 1 is main_id and column 2 is Diabetes
metadata <- clean_data[, c("main_id", "Diabetes")]

# Make Diabetes a readable label for our plot legend
metadata$Diabetes_Label <- factor(metadata$Diabetes, levels = c(0, 1), labels = c("Non-Diabetic", "Diabetic"))

# Extract just the 136 metabolite columns for the math
metabolite_matrix <- clean_data[, 3:ncol(clean_data)]

# 2. Run the PCA Algorithm
cat("Calculating Principal Components...\n")
# Note: scale. = FALSE because we ALREADY Z-score scaled the data in Phase 11!
pca_result <- prcomp(metabolite_matrix, center = FALSE, scale. = FALSE)

# Extract the percentage of variance explained by PC1 and PC2
pca_summary <- summary(pca_result)
var_explained_pc1 <- pca_summary$importance[2, 1] * 100
var_explained_pc2 <- pca_summary$importance[2, 2] * 100

# 3. Create a Plotting Dataframe
# Map the X and Y coordinates (PC1 and PC2) back to the patient IDs
pca_plot_data <- data.frame(
  main_id = metadata$main_id,
  Diabetes_Label = metadata$Diabetes_Label,
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2]
)

# 4. Generate the PC1 vs PC2 Scatter Plot
cat("Generating Plot...\n")
pca_plot <- ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = Diabetes_Label)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "PCA of Patient Metabolic Profiles",
    subtitle = "Are diabetic patients clustering apart from non-diabetic patients?",
    x = paste0("PC1 (", round(var_explained_pc1, 1), "% Variance)"),
    y = paste0("PC2 (", round(var_explained_pc2, 1), "% Variance)")
  ) +
  scale_color_manual(values = c("Non-Diabetic" = "darkgray", "Diabetic" = "dodgerblue")) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold")
  )

# 5. Save the Plot
ggsave("Metabolite_PCA_Plot.png", plot = pca_plot, width = 8, height = 6, dpi = 300)

cat("Success! PCA Plot saved as 'Metabolite_PCA_Plot.png'.\n")