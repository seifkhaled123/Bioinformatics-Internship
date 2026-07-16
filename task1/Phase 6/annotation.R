library(qqman)
library(biomaRt)

# 1. Load your Phase 4 and Phase 5 Results
pc1_results <- read.table("../Phase\ 4/outputs/GWAS_PC1.assoc.linear", header = TRUE, stringsAsFactors = FALSE)
sex_results <- read.table("../Phase\ 5/outputs/GWAS_Sex.assoc.logistic", header = TRUE, stringsAsFactors = FALSE)

# Clean out NAs
pc1_results <- na.omit(pc1_results)
sex_results <- na.omit(sex_results)

# ====================================================================
# TASK 4: Produce Manhattan Plots
# ====================================================================
cat("Drawing Manhattan Plots...\n")

# Phase 4 Plot (PC1 / Ancestry)
jpeg("outputs/Manhattan_PC1.jpg", width = 1000, height = 500)
manhattan(pc1_results, chr="CHR", bp="BP", snp="SNP", p="P", 
          main="GWAS Phase 4: Ancestry Markers (PC1)", 
          suggestiveline = -log10(1e-5), genomewideline = -log10(7.3e-7),
          annotatePval = 1e-23) # This tells qqman to label the top skyscrapers!
dev.off()

# Phase 5 Plot (Biological Sex - The "Noise" check)
jpeg("outputs/Manhattan_Sex.jpg", width = 1000, height = 500)
manhattan(sex_results, chr="CHR", bp="BP", snp="SNP", p="P", 
          main="GWAS Phase 5: Biological Sex (Autosomal Noise)", 
          suggestiveline = -log10(1e-5), genomewideline = -log10(7.3e-7))
dev.off()

# ====================================================================
# TASKS 1, 2 & 3: biomaRt Annotation of Top SNPs
# ====================================================================
cat("Connecting to Ensembl Database to fetch Gene Names...\n")

# Grab the top 5 SNPs from your PC1 run
top_snps <- head(pc1_results[order(pc1_results$P), "SNP"], 5)

# Connect to the human genome database
ensembl <- useEnsembl(biomart = "snps", dataset = "hsapiens_snp")

# Query the database using your rsIDs
annotation_data <- getBM(attributes = c('refsnp_id', 'chr_name', 'chrom_start', 'associated_gene', 'consequence_type_tv'), 
                         filters = 'snp_filter', 
                         values = top_snps, 
                         mart = ensembl)

cat("\n==================================================\n")
cat("PHASE 6: ANNOTATED GENE TABLE\n")
cat("==================================================\n")
print(annotation_data)