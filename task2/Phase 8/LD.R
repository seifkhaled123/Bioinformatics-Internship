source("../functions.R")
library(glue)
library(corrplot)

# 1. Load the data
path <- "../data/GWAS_PC1.assoc.linear"
pc1_values <- read.table(path, header = TRUE, stringsAsFactors = FALSE)

# 2. Sort by P-value (decreasing = FALSE to get the smallest decimals first)
top_pc1 <- head(pc1_values[order(pc1_values$P, decreasing = FALSE), ], 1)

# 3. Extract the exact Coordinates (Chromosome and Base Pair)
mx_chr <- top_pc1$CHR
mx_bp <- top_pc1$BP

run_plink(
    glue("--bfile ../data/dataQC --chr {mx_chr} --from-bp {mx_bp - 500000} --to-bp {mx_bp + 500000} --r2 square --out outputs/LD_chr{mx_chr}_bp{mx_bp}")
)

ld_data <- read.table("outputs/LD_chr11_bp124159136.ld", header = FALSE)

ld_matrix <- as.matrix(ld_data)

# We only plot the 'upper' triangle since the matrix is symmetrical
# tl.pos="n" removes text labels (since plotting hundreds of SNP names overlaps into a messy blur)
corrplot(ld_matrix, 
         method = "color", 
         type = "upper", 
         diag = FALSE,
         tl.pos = "n", 
         cl.lim = c(0, 1), 
         col = colorRampPalette(c("white", "yellow", "red"))(200),
         title = "Linkage Disequilibrium (LD) Heatmap",
         mar = c(0,0,2,0))
