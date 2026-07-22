# Phase 15--16 preparation. These phases intentionally refuse to fabricate
# mixed-model or database annotations when their required inputs are absent.
suppressPackageStartupMessages({ library(dplyr); library(readxl) })
task2_dir <- normalizePath("."); out_dir <- function(x) { d <- file.path(task2_dir, paste("Phase", x), "outputs"); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }

phase15 <- function() {
  source("professional_helpers.R")
  source("mqtl_pipeline.R")
  run_mqtl_phase15()
}

phase16 <- function() {
  out <- out_dir(16)
  metabolites <- read.csv(file.path(task2_dir, "Phase 11", "outputs", "Qatari_metabolomics_Cleaned_Imputed.csv"), check.names = FALSE)
  selected_path <- file.path(task2_dir, "Phase 12", "outputs", "bonferroni_significant_metabolites.txt")
  selected <- if (file.exists(selected_path)) scan(selected_path, what = character(), quiet = TRUE) else character()
  template <- data.frame(Metabolite = names(metabolites)[-(1:2)], Carried_forward = names(metabolites)[-(1:2)] %in% selected, HMDB_ID = NA_character_, KEGG_Compound_ID = NA_character_, PubChem_CID = NA_character_, ChEBI_ID = NA_character_, Chemical_class = NA_character_, Chemical_superclass = NA_character_, Pathway = NA_character_, Evidence_source = NA_character_)
  write.csv(template, file.path(out, "metabolite_annotation_template.csv"), row.names = FALSE)
  writeLines(c("The supplied mapping.xlsx maps participant IDs only; it is not a metabolite identifier cross-reference.", "Fill this template from a versioned HMDB/KEGG/ChEBI source before enrichment.", "Use the QC-passed 136 metabolites as the enrichment background and the carried-forward set as foreground; retain database version and query date."), file.path(out, "ANNOTATION_AND_ENRICHMENT_REQUIREMENTS.txt"))
}

phase <- Sys.getenv("BIOINF_PHASE", unset = commandArgs(trailingOnly = TRUE)[1]); switch(phase, `15` = phase15(), `16` = phase16(), stop("Usage: Rscript future_phases.R {15|16}"))
