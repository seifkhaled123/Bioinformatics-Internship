library(tidyverse)
library(qqman)
library(mclust)
library(plotly)
library(htmlwidgets)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggrepel)

plink   <- "plink1.9"
input   <- "data/Qatari156_filtered_pruned"
results <- "Results"

dir.create(results, showWarnings = FALSE)

run_plink <- function(cmd, step_name = NULL) {
  full_cmd <- paste(plink, cmd)
  status <- system(full_cmd)
  if (status != 0) {
    stop(sprintf(
      "PLINK step failed%s (exit status %d).\nCommand: %s",
      if (!is.null(step_name)) paste0(" [", step_name, "]") else "",
      status, full_cmd
    ))
  }
  invisible(status)
}

read_plink_out <- function(path, ...) {
  if (!file.exists(path)) {
    stop(sprintf("Expected PLINK output file not found: %s", path))
  }
  read.table(path, ...)
}

save_hist <- function(x, file, title, xlab) {
  jpeg(file, 900, 650)
  hist(x, col = "steelblue", breaks = 40, main = title, xlab = xlab)
  dev.off()
}

cat("========== QC ==========\n")

run_plink(
  paste("--bfile", input, "--freq", "--missing", "--out Results/qc"),
  step_name = "freq/missing"
)

frq   <- read_plink_out("Results/qc.frq", header = TRUE)
lmiss <- read_plink_out("Results/qc.lmiss", header = TRUE)
imiss <- read_plink_out("Results/qc.imiss", header = TRUE)

cat("Minimum MAF :", min(frq$MAF), "\n")
cat("Maximum MAF :", max(frq$MAF), "\n")

save_hist(frq$MAF,       "Results/MAF_Histogram.jpg",      "Minor Allele Frequency", "MAF")
save_hist(lmiss$F_MISS,  "Results/SNP_Missingness.jpg",    "SNP Missingness",        "Missing Rate")
save_hist(imiss$F_MISS,  "Results/Sample_Missingness.jpg", "Sample Missingness",     "Missing Rate")

qc <- data.frame()
maf_list  <- c(0.01, 0.03, 0.05)
geno_list <- c(0.10, 0.05, 0.02)
hwe_list  <- c(1e-4, 1e-5, 1e-6)

for (i in maf_list) {
  out <- paste0("Results/maf_", i)
  run_plink(paste("--bfile", input, "--maf", i, "--make-bed", "--out", out),
            step_name = paste("MAF filter", i))
  bim <- read_plink_out(paste0(out, ".bim"))
  qc <- rbind(qc, data.frame(Filter = "MAF", Threshold = i,
                             Remaining = nrow(bim), Reason = "Remove rare variants"))
}

for (i in geno_list) {
  out <- paste0("Results/geno_", i)
  run_plink(paste("--bfile", input, "--geno", i, "--make-bed", "--out", out),
            step_name = paste("GENO filter", i))
  bim <- read_plink_out(paste0(out, ".bim"))
  qc <- rbind(qc, data.frame(Filter = "GENO", Threshold = i,
                             Remaining = nrow(bim), Reason = "Remove SNPs with missing calls"))
}

for (i in hwe_list) {
  out <- paste0("Results/hwe_", i)
  run_plink(paste("--bfile", input, "--hwe", i, "--make-bed", "--out", out),
            step_name = paste("HWE filter", i))
  bim <- read_plink_out(paste0(out, ".bim"))
  qc <- rbind(qc, data.frame(Filter = "HWE", Threshold = i,
                             Remaining = nrow(bim), Reason = "Remove SNPs violating HWE"))
}

run_plink(
  paste("--bfile", input, "--maf 0.05", "--geno 0.05", "--mind 0.05",
        "--hwe 1e-6", "--make-bed", "--out Results/Cleaned"),
  step_name = "final QC"
)

write.csv(qc, "Results/QC_Report.csv", row.names = FALSE)
cat("QC Finished\n")

fam_before <- read_plink_out(paste0(input, ".fam"))
bim_before <- read_plink_out(paste0(input, ".bim"))
fam_after  <- read_plink_out("Results/Cleaned.fam")
bim_after  <- read_plink_out("Results/Cleaned.bim")

summary_qc <- data.frame(
  Metric = c("Samples", "SNPs"),
  Before = c(nrow(fam_before), nrow(bim_before)),
  After  = c(nrow(fam_after), nrow(bim_after))
)
write.csv(summary_qc, "Results/QC_Summary.csv", row.names = FALSE)

cat("\n========== PCA ==========\n")
run_plink("--bfile Results/Cleaned --pca 10 --out Results/PCA", step_name = "PCA")

pcs <- read_plink_out("Results/PCA.eigenvec", header = FALSE)
eig <- read_plink_out("Results/PCA.eigenval", header = FALSE)
colnames(pcs) <- c("FID", "IID", paste0("PC", 1:10))

var_exp <- eig$V1 / sum(eig$V1) * 100
jpeg("Results/ScreePlot.jpg", 900, 650)
plot(var_exp, type = "b", pch = 19,
     xlab = "Principal Component", ylab = "Variance Explained (%)",
     main = "Scree Plot")
dev.off()

p1 <- ggplot(pcs, aes(PC1, PC2)) +
  geom_point(size = 2, color = "steelblue") +
  theme_minimal() +
  labs(title = "Population Structure")
ggsave("Results/PCA_PC1_PC2.jpg", p1, width = 7, height = 5)

cat("\n========== Clustering ==========\n")

wss <- sapply(2:8, function(k) {
  kmeans(pcs[, c("PC1", "PC2")], centers = k, nstart = 25)$tot.withinss
})

jpeg("Results/ElbowPlot.jpg", 900, 650)
plot(2:8, wss, type = "b", pch = 19, xlab = "Clusters", ylab = "Within SS",
     main = "Elbow Method")
dev.off()

mc2 <- Mclust(pcs[, c("PC1", "PC2")])
pcs$Cluster2 <- factor(mc2$classification)

p2 <- ggplot(pcs, aes(PC1, PC2, color = Cluster2)) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "Clusters using PC1 & PC2")
ggsave("Results/Clusters_2PC.jpg", p2, width = 7, height = 5)

mc3 <- Mclust(pcs[, c("PC1", "PC2", "PC3")])
pcs$Cluster3 <- factor(mc3$classification)

fig <- plot_ly(pcs, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Cluster3,
               type = "scatter3d", mode = "markers")
saveWidget(fig, "Results/Clusters_3PC.html", selfcontained = TRUE)

comparison <- data.frame(
  PCs      = c("PC1+PC2", "PC1+PC2+PC3"),
  Clusters = c(mc2$G, mc3$G),
  Max_BIC  = c(max(mc2$BIC), max(mc3$BIC))
)
write.csv(comparison, "Results/Cluster_Comparison.csv", row.names = FALSE)
write.csv(pcs, "Results/PCA_with_Clusters.csv", row.names = FALSE)

cat("PCA & Clustering Finished\n")

cat("\n========== Linear GWAS ==========\n")

write.table(pcs[, c("FID", "IID", "PC1")], "Results/PC1.pheno",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(pcs[, c("FID", "IID", "PC2")], "Results/PC2.pheno",
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

write.table(pcs[, c("FID", "IID", "PC2")], "Results/covar_for_PC1.txt",
            quote = FALSE, row.names = FALSE, sep = "\t")

write.table(pcs[, c("FID", "IID", "PC1")], "Results/covar_for_PC2.txt",
            quote = FALSE, row.names = FALSE, sep = "\t")

run_plink(
  "--bfile Results/Cleaned --pheno Results/PC1.pheno --linear --covar Results/covar_for_PC1.txt --allow-no-sex --out Results/GWAS_PC1",
  step_name = "linear GWAS on PC1"
)
run_plink(
  "--bfile Results/Cleaned --pheno Results/PC2.pheno --linear --covar Results/covar_for_PC2.txt --allow-no-sex --out Results/GWAS_PC2",
  step_name = "linear GWAS on PC2"
)

cat("\n========== Logistic GWAS ==========\n")

run_plink("--bfile Results/Cleaned --check-sex --out Results/SexCheck",
          step_name = "check-sex")
sexcheck <- read_plink_out("Results/SexCheck.sexcheck", header = TRUE)
n_problem <- sum(sexcheck$STATUS == "PROBLEM")
cat("Sex-check flagged", n_problem, "sample(s) as PROBLEM (see Results/SexCheck.sexcheck)\n")

fam <- read_plink_out("Results/Cleaned.fam", header = FALSE)
sex <- fam[, c(1, 2, 5)]
write.table(sex, "Results/Sex.pheno", quote = FALSE, row.names = FALSE,
            col.names = FALSE, sep = "\t")

run_plink(
  "--bfile Results/Cleaned --pheno Results/Sex.pheno --logistic --allow-no-sex --out Results/GWAS_Sex",
  step_name = "logistic GWAS on Sex"
)

plot_gwas <- function(file, prefix) {
  gwas <- read_plink_out(file, header = TRUE, stringsAsFactors = FALSE)
  
  gwas$P   <- as.numeric(as.character(gwas$P))
  gwas$CHR <- as.numeric(as.character(gwas$CHR))
  gwas$BP  <- as.numeric(as.character(gwas$BP))
  
  if ("TEST" %in% colnames(gwas)) {
    gwas <- subset(gwas, TEST == "ADD")
  }
  
  gwas <- subset(gwas, !is.na(P) & !is.na(CHR) & !is.na(BP))
  
  if (nrow(gwas) == 0) {
    warning(paste0(prefix, ": no valid rows to plot after cleaning; skipping Manhattan/QQ."))
    return(NA_real_)
  }
  
  jpeg(paste0("Results/", prefix, "_Manhattan.jpg"), 1000, 700)
  manhattan(gwas, chr = "CHR", bp = "BP", p = "P", snp = "SNP", main = prefix)
  dev.off()
  
  jpeg(paste0("Results/", prefix, "_QQ.jpg"), 800, 600)
  qq(gwas$P, main = paste(prefix, "QQ Plot"))
  dev.off()
  
  chisq  <- qchisq(1 - gwas$P, 1)
  lambda <- median(chisq, na.rm = TRUE) / 0.456
  cat(prefix, " Lambda GC =", round(lambda, 3), "\n")
  
  return(lambda)
}

lambda_pc1 <- plot_gwas("Results/GWAS_PC1.assoc.linear", "PC1")
lambda_pc2 <- plot_gwas("Results/GWAS_PC2.assoc.linear", "PC2")
lambda_sex <- plot_gwas("Results/GWAS_Sex.assoc.logistic", "Sex")

lambda_report <- data.frame(
  Analysis = c("PC1", "PC2", "Sex"),
  Lambda   = c(lambda_pc1, lambda_pc2, lambda_sex)
)
write.csv(lambda_report, "Results/Lambda_Report.csv", row.names = FALSE)

pc1 <- read_plink_out("Results/GWAS_PC1.assoc.linear", header = TRUE, stringsAsFactors = FALSE)
pc2 <- read_plink_out("Results/GWAS_PC2.assoc.linear", header = TRUE, stringsAsFactors = FALSE)
sex <- read_plink_out("Results/GWAS_Sex.assoc.logistic", header = TRUE, stringsAsFactors = FALSE)

if ("TEST" %in% colnames(pc1)) pc1 <- subset(pc1, TEST == "ADD")
if ("TEST" %in% colnames(pc2)) pc2 <- subset(pc2, TEST == "ADD")
if ("TEST" %in% colnames(sex)) sex <- subset(sex, TEST == "ADD")

top_pc1 <- pc1 %>% filter(!is.na(P)) %>% arrange(as.numeric(P)) %>% head(20)
top_pc2 <- pc2 %>% filter(!is.na(P)) %>% arrange(as.numeric(P)) %>% head(20)
top_sex <- sex %>% filter(!is.na(P)) %>% arrange(as.numeric(P)) %>% head(20)

write.csv(top_pc1, "Results/Top_PC1.csv", row.names = FALSE)
write.csv(top_pc2, "Results/Top_PC2.csv", row.names = FALSE)
write.csv(top_sex, "Results/Top_Sex.csv", row.names = FALSE)

cat("GWAS Finished\n")

cat("\n========== SNP Annotation ==========\n")

annotate_gwas <- function(gwas_file, prefix) {
  gwas <- read_plink_out(gwas_file, header = TRUE, stringsAsFactors = FALSE)
  
  gwas$P   <- as.numeric(as.character(gwas$P))
  gwas$CHR <- as.numeric(as.character(gwas$CHR))
  gwas$BP  <- as.numeric(as.character(gwas$BP))
  
  if ("TEST" %in% colnames(gwas)) {
    gwas <- subset(gwas, TEST == "ADD")
  }
  
  gwas <- subset(gwas, !is.na(P) & !is.na(CHR) & !is.na(BP))
  sig  <- subset(gwas, P < 1e-5)
  
  if (nrow(sig) == 0) {
    cat(prefix, " : No significant SNPs (P < 1e-5)\n")
    return(NULL)
  }
  
  mart <- tryCatch(
    useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp"),
    error = function(e) {
      warning(paste0(prefix, ": could not connect to Ensembl SNP mart (",
                     conditionMessage(e), "). Annotation skipped."))
      NULL
    }
  )
  if (is.null(mart)) return(NULL)
  
  ann <- tryCatch(
    getBM(
      attributes = c("refsnp_id", "chr_name", "chrom_start",
                     "ensembl_gene_stable_id", "distance_to_transcript",
                     "consequence_type_tv"),
      filters = "snp_filter",
      values  = sig$SNP,
      mart    = mart
    ),
    error = function(e) {
      warning(paste0(prefix, ": biomaRt query failed (", conditionMessage(e), ")"))
      NULL
    }
  )
  
  if (is.null(ann) || nrow(ann) == 0) {
    cat(prefix, " : SNPs found but no annotations retrieved from biomaRt",
        "(check that SNP IDs are dbSNP rsIDs).\n")
    return(NULL)
  }
  
  gene_mart <- tryCatch(
    useMart("ensembl", dataset = "hsapiens_gene_ensembl"),
    error = function(e) NULL
  )
  if (!is.null(gene_mart) && any(!is.na(ann$ensembl_gene_stable_id) & ann$ensembl_gene_stable_id != "")) {
    gene_map <- tryCatch(
      getBM(
        attributes = c("ensembl_gene_id", "external_gene_name"),
        filters    = "ensembl_gene_id",
        values     = unique(ann$ensembl_gene_stable_id),
        mart       = gene_mart
      ),
      error = function(e) NULL
    )
    if (!is.null(gene_map)) {
      ann <- merge(ann, gene_map, by.x = "ensembl_gene_stable_id",
                   by.y = "ensembl_gene_id", all.x = TRUE)
    }
  }
  if (!"external_gene_name" %in% colnames(ann)) ann$external_gene_name <- NA_character_
  
  write.csv(ann, paste0("Results/", prefix, "_Annotation.csv"), row.names = FALSE)
  
  gwas <- merge(gwas, ann, by.x = "SNP", by.y = "refsnp_id", all.x = TRUE)
  
  p <- ggplot(gwas, aes(BP, -log10(P))) +
    geom_point(aes(color = factor(CHR)), size = 1) +
    theme_minimal() +
    theme(legend.position = "none") +
    geom_text_repel(
      data = subset(gwas, P < 1e-5 & !is.na(external_gene_name) & external_gene_name != ""),
      aes(label = external_gene_name), size = 3, max.overlaps = 20
    ) +
    labs(title = paste(prefix, "GWAS Annotated Manhattan"))
  
  ggsave(paste0("Results/", prefix, "_Gene_Manhattan.jpg"), p, width = 10, height = 5)
  
  return(unique(ann$external_gene_name))
}

genes1 <- annotate_gwas("Results/GWAS_PC1.assoc.linear", "PC1")
genes2 <- annotate_gwas("Results/GWAS_PC2.assoc.linear", "PC2")
genes3 <- annotate_gwas("Results/GWAS_Sex.assoc.logistic", "Sex")

cat("\n========== Pathway Analysis ==========\n")

genes <- unique(c(genes1, genes2, genes3))
genes <- genes[!is.na(genes) & genes != ""]

if (length(genes) > 0) {
  ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  if (nrow(ids) > 0) {
    ego <- enrichGO(gene = ids$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                    ont = "BP", pAdjustMethod = "BH")
    write.csv(as.data.frame(ego), "Results/GO_Enrichment.csv", row.names = FALSE)
    
    jpeg("Results/GO_Dotplot.jpg", 900, 700)
    print(dotplot(ego, showCategory = 10))
    dev.off()
    
    ekegg <- enrichKEGG(gene = ids$ENTREZID, organism = "hsa")
    write.csv(as.data.frame(ekegg), "Results/KEGG_Enrichment.csv", row.names = FALSE)
    
    jpeg("Results/KEGG_Barplot.jpg", 900, 700)
    print(barplot(ekegg, showCategory = 10))
    dev.off()
    
    cat("\nPipeline Completed Successfully!\n")
  } else {
    cat("\nCould not map input Gene Symbols to ENTREZ IDs. Pathway Analysis skipped.\n")
  }
} else {
  cat("\nNo significant genes found across any GWAS runs. Pathway Analysis skipped.\n")
}