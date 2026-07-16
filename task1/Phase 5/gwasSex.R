source("../functions.R")

# ====================================================================
# TASK 1: Confirm Sex is Coded Correctly
# ====================================================================
cat("Checking Sex coding in covariates file...\n")
covariates <- read.table("../Phase\ 4/outputs/covariates.txt", header = TRUE)

# Print a summary table to the console. 
# You should see 1 (Male) and 2 (Female). If you see 0, those are missing.
print(table(covariates$Sex, useNA = "ifany"))

# ====================================================================
# TASK 2: Run Logistic Regression
# ====================================================================
# We use the exact same PLINK architecture, but swap --linear for --logistic.
# --pheno-name Sex tells PLINK to look at the 'Sex' column in our text file.

run_plink("--bfile ../data/dataQC --pheno ../'Phase\ 4'/outputs/covariates.txt --pheno-name Sex --logistic --allow-no-sex --out outputs/GWAS_Sex", step_name = "Logistic GWAS for Sex")

cat("Phase 5 Complete: Logistic association model finished. Check outputs/GWAS_Sex.assoc.logistic\n")