# Load the output files from Phase 4 (Linear) and Phase 5 (Logistic)
# PLINK output files are space/tab separated and have a header row
pc1_results <- read.table("../Phase\ 4/outputs/GWAS_PC1.assoc.linear", header = TRUE, stringsAsFactors = FALSE)
pc2_results <- read.table("../Phase\ 4/outputs/GWAS_PC2.assoc.linear", header = TRUE, stringsAsFactors = FALSE)
sex_results <- read.table("outputs/GWAS_Sex.assoc.logistic", header = TRUE, stringsAsFactors = FALSE)

# Clean out any NA values just in case
pc1_results <- na.omit(pc1_results)
pc2_results <- na.omit(pc2_results)
sex_results <- na.omit(sex_results)

# Sort by P-value (lowest to highest) and grab the top 5
top_pc1 <- head(pc1_results[order(pc1_results$P), ], 5)
top_pc2 <- head(pc2_results[order(pc2_results$P), ], 5)
top_sex <- head(sex_results[order(sex_results$P), ], 5)

cat("\n==================================================\n")
cat("PHASE 4: TOP 5 SNPs DRIVING ANCESTRY (PC1)\n")
cat("==================================================\n")
print(top_pc1[, c("CHR", "SNP", "BP", "BETA", "P")])

cat("\n==================================================\n")
cat("PHASE 5: TOP 5 SNPs DRIVING BIOLOGICAL SEX\n")
cat("==================================================\n")
print(top_sex[, c("CHR", "SNP", "BP", "OR", "P")])