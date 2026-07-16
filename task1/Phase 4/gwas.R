source("../functions.R")

# ====================================================================
# TASK 1: Prepare the Covariate File
# ====================================================================
# 1. Load the fam file to extract Biological Sex (Column 5)
fam <- read.table("../data/dataQC.fam", header = FALSE)
colnames(fam) <- c("FID", "IID", "PID", "MID", "Sex", "Pheno")

# 2. Load the PCs from Phase 2
pcs <- read_plink_out("../Phase\ 2/outputs/PCA.eigenvec", header = FALSE)
colnames(pcs) <- c("FID", "IID", paste0("PC", 1:10))

# 3. Merge them together by Sample ID (FID and IID)
covariates <- merge(fam[, c("FID", "IID", "Sex")], pcs, by = c("FID", "IID"))

# 4. Save the master covariate file (space-delimited, no quotes)
write.table(covariates, "outputs/covariates.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
cat("Task 1 Complete: covariates.txt generated successfully.\n")

# ====================================================================
# TASK 2: Run Association Testing on PC1 and PC2
# ====================================================================
# We tell PLINK to use the eigenvec file as the phenotype file. 
# --mpheno 1 tells PLINK to use the 1st column after FID/IID (which is PC1).
# --mpheno 2 tells PLINK to use the 2nd column after FID/IID (which is PC2).

run_plink("--bfile ../data/dataQC --pheno ../'Phase\ 2'/outputs/PCA.eigenvec --mpheno 1 --linear --allow-no-sex --out outputs/GWAS_PC1", step_name = "GWAS on PC1")

run_plink("--bfile ../data/dataQC --pheno ../'Phase\ 2'/outputs/PCA.eigenvec --mpheno 2 --linear --allow-no-sex --out outputs/GWAS_PC2", step_name = "GWAS on PC2")

cat("Task 2 Complete: Linear association models finished for PC1 and PC2.\n")