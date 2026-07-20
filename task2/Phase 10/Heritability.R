cat("Loading Data...\n")

# 1. Load the Kinship Matrix (GRM) from Phase 9
grm_data <- read.table("../Phase 9/outputs/sample_kinship.rel", header = FALSE)
grm_matrix <- as.matrix(grm_data)

# 2. Load the sample IDs to attach to the matrix
# PLINK generates a .rel.id file alongside the .rel matrix. 
# We need to give our matrix row and column names so the model knows who is who.
grm_ids <- read.table("../Phase 9/outputs/sample_kinship.rel.id", header = FALSE, comment.char="")
# Usually PLINK2 .id files have two columns (FID, IID). We use IID (column 2) or column 1 if there's only one.
sample_names <- as.character(grm_ids$V1) 
if(ncol(grm_ids) > 1) sample_names <- as.character(grm_ids$V2)

rownames(grm_matrix) <- sample_names
colnames(grm_matrix) <- sample_names

# 3. Load a Phenotype to test (We will use your PC1 covariate as a baseline test)
pheno <- read.table("../../task1/Phase 4/outputs/covariates.txt", header = TRUE)

# sommer requires an explicit ID column to match against the GRM row names
pheno$ID <- as.character(pheno$IID) 

# Ensure the phenotype data only includes people actually in our GRM
pheno <- pheno[pheno$ID %in% sample_names, ]

# This breaks the perfect mathematical collinearity caused by the duplicate samples!
diag(grm_matrix) <- diag(grm_matrix) + 0.001    

cat("Fitting Mixed-Effects Model (This involves complex matrix algebra, give it a second)...\n")

# 4. Fit the Variance Component Model
# Formula: PC1 ~ 1 (Intercept only)
# Random Effect: vsr(ID, Gu=grm_matrix) tells the model to use our Kinship matrix to explain variance!
fit <- mmer(fixed = PC1 ~ 1, 
            random = ~vsr(ID, Gu = grm_matrix), 
            rcov = ~units, 
            data = pheno, 
            tolParInv = 1e-3, 
            verbose = FALSE)

# 5. Extract the Heritability (h^2) Manually (Bypassing vpredict)
# Extract Genetic Variance (Vg) - The first component in the sigma list
Vg <- as.numeric(fit$sigma[[1]])

# Extract Environmental/Residual Variance (Ve) - The last component in the sigma list
Ve <- as.numeric(fit$sigma[[length(fit$sigma)]])

# Calculate Heritability manually
h2 <- Vg / (Vg + Ve)

cat("\n--- HERITABILITY REPORT ---\n")
cat("Sample Size (N) :", nrow(pheno), "\n")
cat("Target Trait    : PC1 (Genetic Ancestry Baseline)\n")
cat("Genetic Variance (Vg) :", round(Vg, 4), "\n")
cat("Residual Variance (Ve):", round(Ve, 4), "\n")
cat("Heritability (h2)     :", round(h2, 4), "\n")
cat("---------------------------\n")