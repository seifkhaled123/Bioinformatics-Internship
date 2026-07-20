# ==============================================================================
# Phase 11: Metabolite Data Quality Control Pipeline (Tibble-Proof Version)
# ==============================================================================

library(readxl)

# ------------------------------------------------------------------------------
# 0. Load the Data & Map IDs (Strictly as base Data Frames)
# ------------------------------------------------------------------------------
# Read the excel files and IMMEDIATELY force them out of 'tibble' format
raw_metabolites <- as.data.frame(read_excel("../data/Qatari_metabolomics.xlsx"))
mapping_dict <- as.data.frame(read_excel("../data/mapping.xlsx"))

# Map the IDs using base R merge to avoid any tibble conversion
mapped_data <- merge(raw_metabolites, mapping_dict, by = "mapped_id", all.x = TRUE)

# Isolate the Metadata
metadata <- data.frame(main_id = mapped_data$main_id, Diabetes = mapped_data$Diabetes)
# Now setting row names will work perfectly
rownames(metadata) <- metadata$main_id 

# Isolate the Pure Metabolite Matrix by dropping non-numeric metadata columns
cols_to_drop <- c("mapped_id", "main_id", "Diabetes")
metabolite_matrix <- mapped_data[, !(names(mapped_data) %in% cols_to_drop)]

# Force format and set row names so the math functions track the patients properly
metabolite_matrix <- as.data.frame(metabolite_matrix)
rownames(metabolite_matrix) <- metadata$main_id

# Record initial dimensions for your summary deliverable
num_samples_raw <- nrow(metabolite_matrix)
num_metabolites_raw <- ncol(metabolite_matrix)

# ------------------------------------------------------------------------------
# 1. Assess Metabolite-level Missingness (Clean the Columns)
# ------------------------------------------------------------------------------
# colMeans calculates the exact percentage of missing values per column
metabolite_missing_prop <- colMeans(is.na(metabolite_matrix))

# Apply the 20% (0.20) mathematical gate
metabolites_to_keep <- names(metabolite_missing_prop[metabolite_missing_prop <= 0.20])
metabolites_removed_missing <- names(metabolite_missing_prop[metabolite_missing_prop > 0.20])

# Slice the matrix
qc1_matrix <- metabolite_matrix[, metabolites_to_keep]

# ------------------------------------------------------------------------------
# 2. Assess Sample-level Missingness (Clean the Rows)
# ------------------------------------------------------------------------------
# rowMeans calculates the exact percentage of missing values per row
sample_missing_prop <- rowMeans(is.na(qc1_matrix))

# Apply the 20% (0.20) mathematical gate
samples_to_keep <- rownames(qc1_matrix)[sample_missing_prop <= 0.20]
samples_removed_missing <- rownames(qc1_matrix)[sample_missing_prop > 0.20]

# Slice the matrix
qc2_matrix <- qc1_matrix[samples_to_keep, ]

# ------------------------------------------------------------------------------
# 3. Remove Uninformative Metabolites (Drop Zero-Variance Features)
# ------------------------------------------------------------------------------
# Calculate the variance for every column
metabolite_variances <- apply(qc2_matrix, 2, var, na.rm = TRUE)

# Threshold set to 1e-6 (near-zero variance)
variance_threshold <- 1e-6
metabolites_passed_var <- names(metabolite_variances[metabolite_variances > variance_threshold])
metabolites_removed_var <- names(metabolite_variances[metabolite_variances <= variance_threshold])

# Slice the matrix one last time
qc3_matrix <- qc2_matrix[, metabolites_passed_var]

# ------------------------------------------------------------------------------
# 4. Normalize / Standardize (Z-Score Scaling)
# ------------------------------------------------------------------------------
# Standardize the pure matrix
qc_final_matrix <- as.data.frame(scale(qc3_matrix))

# Record final dimensions
num_samples_qc <- nrow(qc_final_matrix)
num_metabolites_qc <- ncol(qc_final_matrix)

# ------------------------------------------------------------------------------
# 5. Export Deliverables
# ------------------------------------------------------------------------------

# Deliverables A & B: Missingness Tables
write.csv(data.frame(Metabolite = names(metabolite_missing_prop), Missing_Proportion = metabolite_missing_prop), 
          "outputs/Deliverable_Metabolite_Missingness.csv", row.names = FALSE)
write.csv(data.frame(Sample = names(sample_missing_prop), Missing_Proportion = sample_missing_prop), 
          "outputs/Deliverable_Sample_Missingness.csv", row.names = FALSE)

# Deliverables C & D: Lists of Metabolites
all_removed <- c(metabolites_removed_missing, metabolites_removed_var)
if(length(all_removed) == 0) all_removed <- "None"
write.table(all_removed, "outputs/Deliverable_Metabolites_Removed.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(metabolites_passed_var, "outputs/Deliverable_Metabolites_Passed.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Re-attach Metadata to the Clean Data!
final_metadata <- metadata[rownames(qc_final_matrix), ]
final_output <- cbind(final_metadata, qc_final_matrix)

# Save the master dataset
write.csv(final_output, "outputs/Qatari_metabolomics_Cleaned_Scaled.csv", row.names = FALSE)

# Deliverable E: Summary Note (Printed to your console)
cat("\n=== PHASE 11 QC SUMMARY ===\n")
cat("Starting Dimensions:     ", num_samples_raw, "samples |", num_metabolites_raw, "metabolites\n")
cat("Parameters Used:          >20% missingness dropped, Variance <= 1e-6 dropped\n")
cat("--------------------------------------------------\n")
cat("Samples Dropped:         ", length(samples_removed_missing), "\n")
cat("Metabs Dropped (Missing):", length(metabolites_removed_missing), "\n")
cat("Metabs Dropped (Low Var):", length(metabolites_removed_var), "\n")
cat("--------------------------------------------------\n")
cat("Final Clean Dimensions:  ", num_samples_qc, "samples |", num_metabolites_qc, "metabolites\n")
cat("==================================================\n")