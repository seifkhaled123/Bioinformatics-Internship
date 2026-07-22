# GENESIS GRM-adjusted mQTL scan for Task 2 Phase 15.
mqtl_lambda <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  if (!length(p)) return(NA_real_)
  median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
}

mqtl_qq <- function(p, trait, model, max_points = 2500L) {
  p <- sort(p[is.finite(p) & p > 0 & p <= 1])
  if (!length(p)) return(data.frame())
  take <- unique(round(seq(1, length(p), length.out = min(max_points, length(p)))))
  data.frame(Expected = -log10(take / (length(p) + 1)), Observed = -log10(p[take]), Trait = trait, Model = model)
}

run_mqtl_phase15 <- function() {
  required <- c("GENESIS", "SeqArray", "SeqVarTools", "Biobase", "BiocParallel", "ggplot2")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required packages before Phase 15: ", paste(missing, collapse = ", "))
  out <- out_dir(15); full_dir <- file.path(out, "full_results"); dir.create(full_dir, recursive = TRUE, showWarnings = FALSE)
  metabolites <- read.csv(file.path(task2_dir, "Phase 11", "outputs", "Qatari_metabolomics_Cleaned_Imputed.csv"), check.names = FALSE)
  covariates <- read.table(file.path(task2_dir, "..", "task1", "Phase 4", "outputs", "covariates.txt"), header = TRUE, check.names = FALSE)
  selected_path <- file.path(task2_dir, "Phase 12", "outputs", "bonferroni_significant_metabolites.txt")
  traits <- if (file.exists(selected_path)) scan(selected_path, what = character(), quiet = TRUE) else character()
  traits <- intersect(traits, names(metabolites)); if (!length(traits)) stop("No selected metabolites were found for Phase 15.")
  gds_path <- file.path(out, "dataQC.gds")
  if (!file.exists(gds_path)) SeqArray::seqBED2GDS(file.path(task2_dir, "data", "dataQC.bed"), file.path(task2_dir, "data", "dataQC.fam"), file.path(task2_dir, "data", "dataQC.bim"), gds_path, verbose = FALSE)
  gds <- SeqArray::seqOpen(gds_path); on.exit(SeqArray::seqClose(gds), add = TRUE)
  ids <- SeqArray::seqGetData(gds, "sample.id"); variant_id <- SeqArray::seqGetData(gds, "variant.id"); chromosome <- SeqArray::seqGetData(gds, "chromosome")
  autosomal <- variant_id[as.character(chromosome) %in% as.character(1:22)]
  SeqArray::seqSetFilter(gds, variant.id = autosomal, verbose = FALSE)
  annotation <- covariates[match(ids, covariates$IID), c("IID", "Sex", paste0("PC", 1:5))]
  annotation$sample.id <- ids; annotation$Sex <- factor(annotation$Sex)
  if (anyNA(annotation$IID) || !identical(annotation$sample.id, ids)) stop("Genotype, metabolite, and covariate IDs do not align.")
  grm <- as.matrix(read.table(file.path(out_dir(9), "grm.rel"))); grm_ids <- read.table(file.path(out_dir(9), "grm.rel.id"), stringsAsFactors = FALSE)[, 2]
  rownames(grm) <- colnames(grm) <- grm_ids; grm <- grm[ids, ids]; grm <- (grm + t(grm)) / 2; grm_ridge <- max(0, 1e-6 - min(eigen(grm, symmetric = TRUE, only.values = TRUE)$values)); if (grm_ridge > 0) diag(grm) <- diag(grm) + grm_ridge; write.csv(data.frame(GRM_diagonal_ridge = grm_ridge), file.path(out, "grm_stabilization.csv"), row.names = FALSE)
  seq_data <- SeqVarTools::SeqVarData(gds, sampleData = Biobase::AnnotatedDataFrame(annotation)); covar_names <- c("Sex", paste0("PC", 1:5))
  write.csv(data.frame(Metabolite = traits, Phenotype_column = traits, Covariates = paste(covar_names, collapse = ","), Mixed_model = "GENESIS score test with Phase 9 GRM", Status = "completed"), file.path(out, "mqtl_analysis_manifest.csv"), row.names = FALSE)
  write.csv(annotation, file.path(out, "mqtl_covariates.csv"), row.names = FALSE)
  summaries <- list(); top_hits <- list(); qq <- list()
  for (i in seq_along(traits)) {
    SeqArray::seqResetFilter(gds, verbose = FALSE)
    SeqArray::seqSetFilter(gds, variant.id = autosomal, verbose = FALSE)
    trait <- traits[i]; annot <- annotation; annot$outcome <- metabolites[[trait]][match(ids, metabolites$main_id)]
    if (anyNA(annot$outcome)) stop("Missing outcome values remain for ", trait)
    seq_data <- SeqVarTools::SeqVarData(gds, sampleData = Biobase::AnnotatedDataFrame(annot))
    naive_null <- GENESIS::fitNullModel(seq_data, outcome = "outcome", family = "gaussian", verbose = FALSE)
    adjusted_null <- GENESIS::fitNullModel(seq_data, outcome = "outcome", covars = covar_names, cov.mat = grm, family = "gaussian", verbose = FALSE)
    naive <- GENESIS::assocTestSingle(SeqVarTools::SeqVarBlockIterator(seq_data, variantBlock = 10000, verbose = FALSE), naive_null, verbose = FALSE, BPPARAM = BiocParallel::SerialParam())
    SeqArray::seqResetFilter(gds, verbose = FALSE)
    SeqArray::seqSetFilter(gds, variant.id = autosomal, verbose = FALSE)
    adjusted <- GENESIS::assocTestSingle(SeqVarTools::SeqVarBlockIterator(seq_data, variantBlock = 10000, verbose = FALSE), adjusted_null, verbose = FALSE, BPPARAM = BiocParallel::SerialParam())
    naive$SNP <- variant_id[naive$variant.id]; adjusted$SNP <- variant_id[adjusted$variant.id]
    naive$Trait <- trait; naive$Model <- "Naive"; adjusted$Trait <- trait; adjusted$Model <- "GRM + sex + PCs"
    naive$FDR_within_trait <- p.adjust(naive$Score.pval, method = "BH"); adjusted$FDR_within_trait <- p.adjust(adjusted$Score.pval, method = "BH")
    write.csv(naive, gzfile(file.path(full_dir, paste0(make.names(trait), "_naive.csv.gz"))), row.names = FALSE)
    write.csv(adjusted, gzfile(file.path(full_dir, paste0(make.names(trait), "_grm_adjusted.csv.gz"))), row.names = FALSE)
    top_hits[[i]] <- rbind(utils::head(naive[order(naive$Score.pval), ], 10), utils::head(adjusted[order(adjusted$Score.pval), ], 10))
    qq[[length(qq) + 1]] <- mqtl_qq(naive$Score.pval, trait, "Naive"); qq[[length(qq) + 1]] <- mqtl_qq(adjusted$Score.pval, trait, "GRM + sex + PCs")
    summaries[[i]] <- data.frame(Trait = trait, Variants_tested = nrow(adjusted), Naive_lambda = mqtl_lambda(naive$Score.pval), Adjusted_lambda = mqtl_lambda(adjusted$Score.pval), Naive_min_P = min(naive$Score.pval, na.rm = TRUE), Adjusted_min_P = min(adjusted$Score.pval, na.rm = TRUE), Naive_within_trait_FDR_hits = sum(naive$FDR_within_trait < 0.05, na.rm = TRUE), Adjusted_within_trait_FDR_hits = sum(adjusted$FDR_within_trait < 0.05, na.rm = TRUE), Null_model_genetic_variance = unname(adjusted_null$varComp[1]))
  }
  summary_table <- dplyr::bind_rows(summaries); top_table <- dplyr::bind_rows(top_hits); qq_table <- dplyr::bind_rows(qq)
  write.csv(summary_table, file.path(out, "mqtl_model_comparison.csv"), row.names = FALSE); write.csv(top_table, file.path(out, "mqtl_top_hits.csv"), row.names = FALSE); write.csv(qq_table, file.path(out, "mqtl_qq_points.csv"), row.names = FALSE)
  p <- ggplot2::ggplot(qq_table, ggplot2::aes(Expected, Observed, colour = Model)) + ggplot2::geom_point(alpha = 0.16, size = 0.35) + ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) + ggplot2::facet_wrap(~Trait, scales = "free") + ggplot2::scale_colour_manual(values = c("Naive" = "#9E9E9E", "GRM + sex + PCs" = "#0072B2")) + ggplot2::labs(title = "mQTL calibration: naive vs GRM-adjusted models", x = "Expected -log10(P)", y = "Observed -log10(P)") + ggplot2::theme_bw(base_size = 10) + ggplot2::theme(legend.position = "bottom")
  ggplot2::ggsave(file.path(out, "mqtl_qq_comparison.png"), p, width = 14, height = 11, dpi = 260)
  best <- summary_table[order(summary_table$Adjusted_min_P), ]
  p2 <- ggplot2::ggplot(best, ggplot2::aes(x = stats::reorder(Trait, Adjusted_min_P), y = -log10(Adjusted_min_P))) + ggplot2::geom_col(fill = "#0072B2") + ggplot2::coord_flip() + ggplot2::labs(title = "Strongest GRM-adjusted mQTL signal by metabolite", x = NULL, y = "-log10(minimum P)") + ggplot2::theme_bw(base_size = 11)
  ggplot2::ggsave(file.path(out, "mqtl_strongest_signal_by_trait.png"), p2, width = 9, height = 10, dpi = 260)
  writeLines(c("Each compressed result file contains all autosomal variants for one metabolite and model.", "Naive models include no covariates or GRM; adjusted models use sex, PC1-PC5, and the Phase 9 genomic relationship matrix.", "FDR is calculated within each metabolite scan. Interpret hits only after confirming genome build for annotation."), file.path(out, "README_RESULTS.txt"))
  write_session_info(file.path(out, "sessionInfo.txt"))
}
