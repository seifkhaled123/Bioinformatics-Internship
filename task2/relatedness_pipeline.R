# Task 2 Phases 8--10: LD, KING/GRM, and baseline metabolite heritability.
# Run from task2: Rscript relatedness_pipeline.R {8|9|10}
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
task2_dir <- normalizePath("."); root <- normalizePath("..")
out_dir <- function(phase) { x <- file.path(task2_dir, paste("Phase", phase), "outputs"); dir.create(x, recursive = TRUE, showWarnings = FALSE); x }
plink1 <- Sys.getenv("PLINK_BIN", "plink1.9")
run1 <- function(args) { if (system2(plink1, shQuote(args)) != 0) stop("PLINK 1.9 command failed.") }

phase8 <- function() {
  out <- out_dir(8); association <- read.table(file.path(root, "task1", "Phase 4", "outputs", "GWAS_PC1.assoc.linear"), header = TRUE)
  association <- association[association$TEST == "ADD" & is.finite(association$P), ]; lead <- association[which.min(association$P), ]
  start <- max(1, lead$BP - 500000); end <- lead$BP + 500000; prefix <- file.path(out, sprintf("LD_chr%s_bp%s", lead$CHR, lead$BP))
  run1(c("--bfile", file.path(task2_dir, "data", "dataQC"), "--chr", lead$CHR, "--from-bp", start, "--to-bp", end, "--r2", "square", "--out", prefix))
  ld <- as.matrix(read.table(paste0(prefix, ".ld"), header = FALSE)); write.csv(data.frame(Lead_SNP = lead$SNP, Chromosome = lead$CHR, Position = lead$BP, Window_start = start, Window_end = end, SNPs_in_window = nrow(ld)), file.path(out, "ld_region_summary.csv"), row.names = FALSE); write.csv(ld, file.path(out, "ld_r2_matrix.csv"), row.names = FALSE)
  png(file.path(out, "ld_heatmap.png"), width = 2200, height = 1800, res = 250); image(seq_len(nrow(ld)), seq_len(ncol(ld)), ld, col = colorRampPalette(c("white", "#FEE08B", "#D73027"))(200), axes = FALSE, xlab = "SNP order in 1 Mb window", ylab = "SNP order in 1 Mb window", main = sprintf("LD (r²) around %s, chr%s:%s", lead$SNP, lead$CHR, lead$BP)); box(); dev.off()
}

phase9 <- function() {
  out <- out_dir(9); input <- file.path(task2_dir, "data", "dataQC")
  # GRM is generated with installed PLINK 1.9. KING robust output deliberately requires PLINK 2.
  run1(c("--bfile", input, "--make-rel", "square", "--out", file.path(out, "grm")))
  plink2 <- Sys.getenv("PLINK2_BIN", "plink2")
  if (!nzchar(Sys.which(plink2))) { writeLines(c("KING-robust kinship has not been substituted with a GRM.", "Install PLINK 2 and rerun with PLINK2_BIN=/path/to/plink2.", "Required command: --make-king-table."), file.path(out, "KING_REQUIRED.txt")); stop("PLINK 2 is required for the KING-robust deliverable; GRM was written, but no relatedness categories were reported.") }
  status <- system2(plink2, shQuote(c("--bfile", input, "--make-king-table", "--out", file.path(out, "king"))))
  if (status != 0) stop("PLINK 2 KING calculation failed.")
  king_path <- file.path(out, "king.kin0"); king <- read.table(king_path, header = TRUE, check.names = FALSE, comment.char = "")
  kinship_col <- grep("KINSHIP", names(king), value = TRUE)[1]; ids <- grep("^IID[12]$", names(king), value = TRUE)
  king$Relatedness <- cut(king[[kinship_col]], c(-Inf, .0442, .0884, .177, .354, Inf), right = FALSE, labels = c("Unrelated", "3rd-degree", "2nd-degree", "1st-degree", "Duplicate/MZ"))
  summary <- as.data.frame(table(king$Relatedness)); names(summary) <- c("Relatedness", "Pairs"); write.csv(summary, file.path(out, "king_relatedness_summary.csv"), row.names = FALSE); related <- king[king[[kinship_col]] >= .0442, ]; write.csv(related, file.path(out, "king_related_pairs.csv"), row.names = FALSE); write.csv(data.frame(Individuals_in_nonunrelated_pairs = unique(unlist(related[, ids]))), file.path(out, "individuals_in_related_pairs.csv"), row.names = FALSE)
}

phase10 <- function() {
  if (!requireNamespace("sommer", quietly = TRUE)) stop("Install sommer before Phase 10.")
  out <- out_dir(10); grm <- as.matrix(read.table(file.path(out_dir(9), "grm.rel"))); ids <- read.table(file.path(out_dir(9), "grm.rel.id"), stringsAsFactors = FALSE)[, 2]; rownames(grm) <- colnames(grm) <- ids
  metabolomics <- read.csv(file.path(out_dir(11), "Qatari_metabolomics_Cleaned_Imputed.csv"), check.names = FALSE); trait <- Sys.getenv("HERITABILITY_TRAIT", "Glucose"); if (!trait %in% names(metabolomics)) stop("Trait not found: ", trait)
  pheno <- data.frame(ID = metabolomics$main_id, trait = metabolomics[[trait]]); pheno <- pheno[pheno$ID %in% ids, ]; grm <- grm[pheno$ID, pheno$ID]
  fit <- sommer::mmer(fixed = trait ~ 1, random = ~sommer::vsr(ID, Gu = grm), rcov = ~units, data = pheno, verbose = FALSE); vc <- unlist(fit$sigma); h2 <- vc[1] / sum(vc); h2_delta <- tryCatch(sommer::vpredict(fit, h2 ~ V1/(V1 + V2)), error = function(e) c(Estimate = NA_real_, SE = NA_real_)); h2_delta <- as.numeric(h2_delta); h2_se <- if (length(h2_delta) >= 2) h2_delta[2] else NA_real_; boundary <- isTRUE(vc[1] <= 1e-10); note <- if (boundary) "Genetic variance reached the non-negative boundary; h2 is estimated as 0 and a delta-method SE is not available." else "Variance components converged; h2 SE is from sommer::vpredict."; write.csv(data.frame(Trait = trait, N = nrow(pheno), Genetic_variance = vc[1], Residual_variance = vc[length(vc)], Heritability = h2, Heritability_SE_delta_method = h2_se, Boundary_estimate = boundary, Model_converged = fit$convergence, Note = note), file.path(out, "heritability_result.csv"), row.names = FALSE); write.csv(data.frame(Model_converged = fit$convergence, AIC = fit$AIC, Genetic_variance_boundary = boundary, Participants = nrow(pheno)), file.path(out, "heritability_diagnostics.csv"), row.names = FALSE); png(file.path(out, "variance_components.png"), width = 1600, height = 1100, res = 180); barplot(c(Genetic = vc[1], Residual = vc[length(vc)]), col = c("#2A9D8F", "#457B9D"), ylab = "Estimated variance", main = sprintf("Variance components for %s", trait)); dev.off()
}

phase <- Sys.getenv("BIOINF_PHASE", unset = commandArgs(trailingOnly = TRUE)[1]); switch(phase, `8` = phase8(), `9` = phase9(), `10` = phase10(), stop("Usage: Rscript relatedness_pipeline.R {8|9|10}"))
