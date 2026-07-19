source("../functions.R")
library(corrplot)

run_plink(
    "--bfile ../data/dataQC --make-rel square --out outputs/sample_kinship"
)

# 1. Load the Data
kin_data <- read.table("outputs/sample_kinship.rel", header = FALSE)
kin_matrix <- as.matrix(kin_data)

# 2. Deliverable A: Plot the Kinship Heatmap
# Most of it will be white (unrelated), but related pairs will show up as colored dots off the diagonal.
corrplot(kin_matrix, 
         method = "color", 
         tl.pos = "n", 
         is.corr = FALSE,   # <--- Change 'iscale' to 'is.corr' here!
         col = colorRampPalette(c("white", "orange", "red"))(200),
         title = "Genomic Kinship Matrix (Sample Relatedness)",
         mar = c(0,0,2,0))

# 3. Deliverable B: Count the Relatives
# We only check the "upper triangle" of the matrix so we don't count the same pair of people twice, 
# and we ignore the diagonal (because a person is 100% related to themselves).
kin_matrix[lower.tri(kin_matrix, diag = TRUE)] <- NA
vals <- kin_matrix[!is.na(kin_matrix)]

# Filter for pairs with relatedness > 0.05 (Task requirement)
related_vals <- vals[vals > 0.05]

# Standard Biological Thresholds for Kinship
mz_twins   <- sum(related_vals > 0.354)
first_deg  <- sum(related_vals > 0.177 & related_vals <= 0.354)
second_deg <- sum(related_vals > 0.088 & related_vals <= 0.177)
third_deg  <- sum(related_vals > 0.05  & related_vals <= 0.088)

# Print out the Table for your Deliverable
cat("\n--- RELATEDNESS REPORT ---\n")
cat("Identical Twins / Duplicates :", mz_twins, "\n")
cat("1st Degree (Siblings/Parent):", first_deg, "\n")
cat("2nd Degree (Half-Sib/Uncle) :", second_deg, "\n")
cat("3rd Degree (Cousins)        :", third_deg, "\n")
cat("----------------------------\n")
cat("Total non-unrelated pairs identified in the cohort:", length(related_vals), "\n")