# Reproducible Task 1 pipeline. Run from task1: Rscript professional_pipeline.R 1
require_pkgs <- function(x) { miss <- x[!vapply(x, requireNamespace, logical(1), quietly = TRUE)]; if (length(miss)) stop("Install: ", paste(miss, collapse = ", ")) }
require_pkgs(c("ggplot2", "dplyr", "mclust", "cluster", "qqman"))
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

task1_dir <- normalizePath(".")
out_dir <- function(phase) { x <- file.path(task1_dir, paste("Phase", phase), "outputs"); dir.create(x, recursive = TRUE, showWarnings = FALSE); x }
plink <- Sys.getenv("PLINK_BIN", "plink1.9")
run_plink <- function(args, label) { status <- system2(plink, shQuote(args)); if (status != 0) stop("PLINK failed at ", label) }
theme_report <- function() theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold", colour = "#17324D"), plot.subtitle = element_text(colour = "#4E6476"), panel.grid.minor = element_blank(), legend.title = element_blank())
save_fig <- function(p, filename, w = 8, h = 5) ggsave(filename, p, width = w, height = h, dpi = 320, bg = "white")
read_assoc <- function(path) { d <- read.table(path, header = TRUE, check.names = FALSE); if ("TEST" %in% names(d)) d <- d[d$TEST == "ADD", ]; d[is.finite(d$P) & d$P > 0, ] }
lambda_gc <- function(p) median(qchisq(1 - p, 1)) / qchisq(.5, 1)

phase1 <- function() {
  input <- file.path(task1_dir, "data", "Qatari156_filtered_pruned")
  out <- out_dir(1); prefix <- file.path(out, "qc")
  run_plink(c("--bfile", input, "--freq", "--missing", "--out", prefix), "QC metrics")
  frq <- read.table(paste0(prefix, ".frq"), header = TRUE); lmiss <- read.table(paste0(prefix, ".lmiss"), header = TRUE); imiss <- read.table(paste0(prefix, ".imiss"), header = TRUE)
  base_n <- nrow(read.table(paste0(input, ".bim"))); base_samples <- nrow(read.table(paste0(input, ".fam")))
  strategies <- list(Standard = c("--maf", ".05", "--geno", ".05", "--mind", ".05", "--hwe", "1e-6"), Strict_missingness = c("--maf", ".05", "--geno", ".001", "--mind", ".05", "--hwe", "1e-6"), Strict_MAF = c("--maf", ".20", "--geno", ".05", "--mind", ".05", "--hwe", "1e-6"))
  counts <- lapply(names(strategies), function(name) { p <- file.path(out, paste0("qc_", name)); run_plink(c("--bfile", input, strategies[[name]], "--make-bed", "--out", p), name); data.frame(Strategy = name, Samples = nrow(read.table(paste0(p, ".fam"))), SNPs = nrow(read.table(paste0(p, ".bim")))) }) |> bind_rows()
  write.csv(bind_rows(data.frame(Strategy = "Before_QC", Samples = base_samples, SNPs = base_n), counts), file.path(out, "qc_filter_waterfall.csv"), row.names = FALSE)
  plots <- list(MAF = data.frame(value = frq$MAF, type = "Minor allele frequency"), SNP_missingness = data.frame(value = lmiss$F_MISS, type = "SNP missingness"), Sample_missingness = data.frame(value = imiss$F_MISS, type = "Sample missingness"))
  for (name in names(plots)) { p <- ggplot(plots[[name]], aes(value)) + geom_histogram(bins = 40, fill = "#0072B2", colour = "white") + labs(title = gsub("_", " ", name), x = "Proportion", y = "Count") + theme_report(); save_fig(p, file.path(out, paste0(tolower(name), ".png"))) }
  write.csv(data.frame(Metric = c("MAF_min", "MAF_max", "SNP_missingness_max", "sample_missingness_max"), Value = c(min(frq$MAF), max(frq$MAF), max(lmiss$F_MISS), max(imiss$F_MISS))), file.path(out, "qc_summary.csv"), row.names = FALSE)
}

phase2 <- function() {
  input <- file.path(out_dir(1), "qc_Standard"); out <- out_dir(2); prune <- file.path(out, "pca_pruning")
  run_plink(c("--bfile", input, "--indep-pairwise", "200", "50", ".2", "--out", prune), "LD pruning for PCA")
  run_plink(c("--bfile", input, "--extract", paste0(prune, ".prune.in"), "--pca", "10", "--out", file.path(out, "PCA")), "PCA")
  pcs <- read.table(file.path(out, "PCA.eigenvec"), header = FALSE); names(pcs) <- c("FID", "IID", paste0("PC", 1:10)); eig <- scan(file.path(out, "PCA.eigenval")); variance <- 100 * eig / sum(eig)
  write.csv(data.frame(PC = seq_along(variance), Eigenvalue = eig, Variance_explained_percent = variance), file.path(out, "pca_variance_explained.csv"), row.names = FALSE); write.csv(pcs, file.path(out, "pca_scores.csv"), row.names = FALSE)
  p1 <- ggplot(pcs, aes(PC1, PC2)) + geom_point(size = 2.4, alpha = .8, colour = "#0072B2") + labs(title = "Genetic population structure", x = sprintf("PC1 (%.1f%%)", variance[1]), y = sprintf("PC2 (%.1f%%)", variance[2])) + theme_report(); save_fig(p1, file.path(out, "pca_pc1_pc2.png"))
  p2 <- ggplot(data.frame(PC = 1:10, Variance = variance[1:10]), aes(PC, Variance)) + geom_line(colour = "#0072B2") + geom_point(size = 2.5, colour = "#0072B2") + scale_x_continuous(breaks = 1:10) + labs(title = "PCA scree plot", x = "Principal component", y = "Variance explained (%)") + theme_report(); save_fig(p2, file.path(out, "pca_scree_plot.png"))
}

phase3 <- function() {
  out <- out_dir(3); pcs <- read.csv(file.path(out_dir(2), "pca_scores.csv")); x <- scale(pcs[, paste0("PC", 1:3)])
  candidates <- 2:8; diagnostics <- lapply(candidates, function(k) { fit <- kmeans(x, centers = k, nstart = 100); sil <- cluster::silhouette(fit$cluster, dist(x)); data.frame(k = k, WSS = fit$tot.withinss, Mean_silhouette = mean(sil[, "sil_width"])) }) |> bind_rows(); best_k <- diagnostics$k[which.max(diagnostics$Mean_silhouette)]
  fit2 <- kmeans(x[, 1:2], centers = best_k, nstart = 100); fit3 <- kmeans(x, centers = best_k, nstart = 100); pcs$Cluster_2PC <- factor(fit2$cluster); pcs$Cluster_3PC <- factor(fit3$cluster)
  write.csv(diagnostics, file.path(out, "cluster_diagnostics.csv"), row.names = FALSE); write.csv(pcs, file.path(out, "cluster_assignments.csv"), row.names = FALSE); write.csv(as.data.frame(table(pcs$Cluster_3PC)), file.path(out, "cluster_sizes.csv"), row.names = FALSE)
  diag_plot <- ggplot(diagnostics, aes(k, Mean_silhouette)) + geom_line(colour = "#0072B2") + geom_point(size = 2.5, colour = "#0072B2") + labs(title = "Cluster-number selection", subtitle = paste("Selected k =", best_k, "by mean silhouette width"), x = "k", y = "Mean silhouette width") + theme_report(); save_fig(diag_plot, file.path(out, "cluster_silhouette_selection.png"))
  cluster_plot <- ggplot(pcs, aes(PC1, PC2, colour = Cluster_3PC)) + geom_point(size = 2.5, alpha = .85) + labs(title = "PCA clusters", subtitle = "Unsupervised clusters; ancestry labels require external reference validation", x = "PC1", y = "PC2") + theme_report(); save_fig(cluster_plot, file.path(out, "pca_clusters_2d.png"))
}

write_gwas_figures <- function(d, out, label) {
  threshold <- .05 / nrow(d); jpeg(file.path(out, paste0(label, "_manhattan.jpg")), 1800, 900, res = 180); qqman::manhattan(d, main = paste(label, "(descriptive PC association)"), suggestiveline = -log10(1e-5), genomewideline = -log10(threshold)); dev.off(); jpeg(file.path(out, paste0(label, "_qq.jpg")), 1200, 900, res = 180); qqman::qq(d$P, main = paste(label, "QQ plot")); dev.off(); data.frame(Analysis = label, Tests = nrow(d), Bonferroni_threshold = threshold, Lambda_GC = lambda_gc(d$P), Significant_SNPs = sum(d$P < threshold))
}

phase4 <- function() {
  out <- out_dir(4); input <- file.path(out_dir(1), "qc_Standard"); pcs <- read.csv(file.path(out_dir(2), "pca_scores.csv")); fam_master <- read.table(paste0(input, ".fam")); names(fam_master) <- c("FID", "IID", "PID", "MID", "Sex", "Phenotype"); covariates_master <- merge(fam_master[, c("FID", "IID", "Sex")], pcs, by = c("FID", "IID")); write.table(covariates_master, file.path(out, "covariates.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t");
  for (trait in c("PC1", "PC2")) { other <- setdiff(paste0("PC", 1:5), trait); pheno <- pcs[, c("FID", "IID", trait)]; cov <- pcs[, c("FID", "IID", other)]; write.table(pheno, file.path(out, paste0(trait, ".pheno")), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t"); write.table(cov, file.path(out, paste0(trait, ".covar")), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"); run_plink(c("--bfile", input, "--pheno", file.path(out, paste0(trait, ".pheno")), "--covar", file.path(out, paste0(trait, ".covar")), "--linear", "hide-covar", "--out", file.path(out, paste0("GWAS_", trait))), paste("GWAS", trait)) }
  results <- lapply(c("PC1", "PC2"), function(pc) { d <- read_assoc(file.path(out, paste0("GWAS_", pc, ".assoc.linear"))); write.csv(head(d[order(d$P), ], 50), file.path(out, paste0(pc, "_top_hits.csv")), row.names = FALSE); write_gwas_figures(d, out, pc) }) |> bind_rows(); write.csv(results, file.path(out, "gwas_summary.csv"), row.names = FALSE)
}

phase5 <- function() {
  out <- out_dir(5); input <- file.path(out_dir(1), "qc_Standard"); fam <- read.table(paste0(input, ".fam")); pcs <- read.csv(file.path(out_dir(2), "pca_scores.csv")); fam_master <- read.table(paste0(input, ".fam")); names(fam_master) <- c("FID", "IID", "PID", "MID", "Sex", "Phenotype"); covariates_master <- merge(fam_master[, c("FID", "IID", "Sex")], pcs, by = c("FID", "IID")); write.table(covariates_master, file.path(out, "covariates.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"); names(fam) <- c("FID", "IID", "PID", "MID", "Sex", "Phenotype"); pheno <- fam[, c("FID", "IID", "Sex")]; cov <- pcs[, c("FID", "IID", paste0("PC", 1:5))]; write.table(pheno, file.path(out, "sex.pheno"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t"); write.table(cov, file.path(out, "sex.covar"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t"); run_plink(c("--bfile", input, "--pheno", file.path(out, "sex.pheno"), "--covar", file.path(out, "sex.covar"), "--logistic", "hide-covar", "--out", file.path(out, "GWAS_Sex")), "sex GWAS")
  d <- read_assoc(file.path(out, "GWAS_Sex.assoc.logistic")); write.csv(head(d[order(d$P), ], 50), file.path(out, "sex_top_hits.csv"), row.names = FALSE); write.csv(write_gwas_figures(d, out, "Sex"), file.path(out, "sex_gwas_summary.csv"), row.names = FALSE); write.csv(as.data.frame(table(Sex_code = fam$Sex)), file.path(out, "sex_coding_check.csv"), row.names = FALSE)
}

phase6 <- function() {
  out <- out_dir(6); gwas <- bind_rows(lapply(c("PC1", "PC2"), function(pc) mutate(read_assoc(file.path(out_dir(4), paste0("GWAS_", pc, ".assoc.linear"))), Analysis = pc))); threshold <- .05 / nrow(filter(gwas, Analysis == "PC1")); leads <- gwas |> filter(P < threshold) |> group_by(Analysis, CHR) |> slice_min(P, n = 1, with_ties = FALSE) |> ungroup() |> select(Analysis, SNP, CHR, BP, P, BETA); write.csv(leads, file.path(out, "annotation_input_lead_loci.csv"), row.names = FALSE); writeLines(c("Genome build must be confirmed before annotation (GRCh37 vs GRCh38).", "This file deliberately contains no hard-coded annotations.", "Use the optional annotation script only after confirming the build."), file.path(out, "annotation_status.txt"))
}

phase7 <- function() {
  out <- out_dir(7); input <- file.path(out_dir(6), "annotated_lead_genes.csv"); if (!file.exists(input)) { writeLines("No enrichment run: first create annotated_lead_genes.csv from confirmed-build Phase 6 annotations. No hard-coded gene set is used.", file.path(out, "enrichment_status.txt")); return(invisible()) }; require_pkgs(c("clusterProfiler", "org.Hs.eg.db")); genes <- read.csv(input)$Gene; ids <- clusterProfiler::bitr(unique(genes), "SYMBOL", "ENTREZID", org.Hs.eg.db); ego <- clusterProfiler::enrichGO(ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", readable = TRUE); write.csv(as.data.frame(ego), file.path(out, "go_enrichment.csv"), row.names = FALSE); save_fig(clusterProfiler::dotplot(ego, showCategory = 15) + ggtitle("GO enrichment of observed GWAS lead-locus genes"), file.path(out, "go_enrichment_dotplot.png"), 9, 7)
}

phase <- Sys.getenv("BIOINF_PHASE", unset = commandArgs(trailingOnly = TRUE)[1])
switch(phase, `1` = phase1(), `2` = phase2(), `3` = phase3(), `4` = phase4(), `5` = phase5(), `6` = phase6(), `7` = phase7(), stop("Usage: Rscript professional_pipeline.R {1|2|3|4|5|6|7}"))
