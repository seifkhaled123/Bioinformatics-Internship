# Confirmed-build annotation and observed-gene enrichment for Phases 6 and 7.

read_genome_build <- function() {
  path <- file.path(task1_dir, "data", "genome_build.tsv")
  if (!file.exists(path)) stop("Missing ", path, ". Confirm and record the source assembly before annotation.")
  build <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("UCSC_name", "NCBI_assembly")
  if (nrow(build) != 1 || !all(required %in% names(build))) {
    stop("genome_build.tsv must contain exactly one row and columns: ", paste(required, collapse = ", "))
  }
  if (!build$UCSC_name %in% c("hg18", "hg19", "hg38")) stop("Unsupported source assembly: ", build$UCSC_name)
  build
}

first_or_na <- function(x, name) {
  value <- x[[name]]
  if (is.null(value) || !length(value)) NA_character_ else as.character(value[[1]])
}

fetch_vep_annotations <- function(ids) {
  require_pkgs(c("httr2", "jsonlite"))
  endpoint <- Sys.getenv("ENSEMBL_REST_URL", "https://rest.ensembl.org")
  response <- httr2::request(paste0(endpoint, "/vep/human/id")) |>
    httr2::req_url_query(symbol = 1, canonical = 1, pick = 1) |>
    httr2::req_headers(Accept = "application/json") |>
    httr2::req_body_json(list(ids = unname(ids))) |>
    httr2::req_retry(max_tries = 6, retry_on_failure = TRUE) |>
    httr2::req_perform()
  records <- httr2::resp_body_json(response, simplifyVector = FALSE)
  if (!length(records)) stop("Ensembl VEP returned no annotations.")
  bind_rows(lapply(records, function(x) {
    tc <- x$transcript_consequences
    chosen <- if (is.null(tc) || !length(tc)) list() else tc[[1]]
    data.frame(
      SNP = first_or_na(x, "id"),
      Annotation_assembly = first_or_na(x, "assembly_name"),
      Annotation_CHR = first_or_na(x, "seq_region_name"),
      Annotation_BP = suppressWarnings(as.numeric(first_or_na(x, "start"))),
      Gene = first_or_na(chosen, "gene_symbol"),
      Ensembl_gene_id = first_or_na(chosen, "gene_id"),
      Transcript = first_or_na(chosen, "transcript_id"),
      Consequence = first_or_na(x, "most_severe_consequence"),
      Impact = first_or_na(chosen, "impact"),
      Gene_biotype = first_or_na(chosen, "biotype"),
      stringsAsFactors = FALSE
    )
  }))
}
fetch_nearest_gene <- function(chromosome, position, max_distance = 1e6) {
  endpoint <- Sys.getenv("ENSEMBL_REST_URL", "https://rest.ensembl.org")
  region_start <- max(1, position - max_distance)
  region_end <- position + max_distance
  url <- paste0(endpoint, "/overlap/region/human/", chromosome, ":", region_start, "-", region_end)
  response <- httr2::request(url) |>
    httr2::req_url_query(feature = "gene", biotype = "protein_coding") |>
    httr2::req_headers(Accept = "application/json") |>
    httr2::req_retry(max_tries = 6, retry_on_failure = TRUE) |>
    httr2::req_perform()
  records <- httr2::resp_body_json(response, simplifyVector = FALSE)
  if (!length(records)) {
    return(data.frame(Gene = NA_character_, Ensembl_gene_id = NA_character_, Gene_biotype = NA_character_, Distance_to_gene_bp = NA_real_))
  }
  genes <- bind_rows(lapply(records, function(x) {
    gene_start <- as.numeric(first_or_na(x, "start"))
    gene_end <- as.numeric(first_or_na(x, "end"))
    distance <- if (position < gene_start) gene_start - position else if (position > gene_end) position - gene_end else 0
    data.frame(
      Gene = first_or_na(x, "external_name"),
      Ensembl_gene_id = first_or_na(x, "id"),
      Gene_biotype = first_or_na(x, "biotype"),
      Distance_to_gene_bp = distance,
      stringsAsFactors = FALSE
    )
  })) |>
    filter(!is.na(Gene), nzchar(Gene)) |>
    arrange(Distance_to_gene_bp, Gene)
  if (!nrow(genes)) {
    return(data.frame(Gene = NA_character_, Ensembl_gene_id = NA_character_, Gene_biotype = NA_character_, Distance_to_gene_bp = NA_real_))
  }
  genes[1, , drop = FALSE]
}


ensure_plot_associations <- function() {
  linear_files <- file.path(out_dir(4), paste0("GWAS_", c("PC1", "PC2"), ".assoc.linear"))
  logistic_file <- file.path(out_dir(5), "GWAS_Sex.assoc.logistic")
  if (all(file.exists(c(linear_files, logistic_file)))) return(invisible())

  prefixes <- c(file.path(out_dir(1), "qc_Standard"), file.path(task1_dir, "data", "Qatari156_filtered_pruned"), file.path(task1_dir, "data", "dataQC"))
  available <- prefixes[file.exists(paste0(prefixes, ".bed")) & file.exists(paste0(prefixes, ".bim")) & file.exists(paste0(prefixes, ".fam"))]
  if (!length(available)) stop("No PLINK genotype prefix is available to recreate Manhattan plot inputs.")
  input <- available[[1]]
  pcs <- read.csv(file.path(out_dir(2), "pca_scores.csv"))

  if (!all(file.exists(linear_files))) {
    for (trait in c("PC1", "PC2")) {
      other <- setdiff(paste0("PC", 1:5), trait)
      pheno_path <- file.path(out_dir(4), paste0(trait, ".pheno"))
      covar_path <- file.path(out_dir(4), paste0(trait, ".covar"))
      write.table(pcs[, c("FID", "IID", trait)], pheno_path, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
      write.table(pcs[, c("FID", "IID", other)], covar_path, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
      run_plink(c("--bfile", input, "--pheno", pheno_path, "--covar", covar_path, "--linear", "hide-covar", "--out", file.path(out_dir(4), paste0("GWAS_", trait))), paste("GWAS", trait, "for labeled Manhattan plot"))
    }
  }

  if (!file.exists(logistic_file)) {
    fam <- read.table(paste0(input, ".fam"))
    names(fam) <- c("FID", "IID", "PID", "MID", "Sex", "Phenotype")
    pheno_path <- file.path(out_dir(5), "sex.pheno")
    covar_path <- file.path(out_dir(5), "sex.covar")
    write.table(fam[, c("FID", "IID", "Sex")], pheno_path, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    write.table(pcs[, c("FID", "IID", paste0("PC", 1:5))], covar_path, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    run_plink(c("--bfile", input, "--pheno", pheno_path, "--covar", covar_path, "--logistic", "hide-covar", "--out", file.path(out_dir(5), "GWAS_Sex")), "sex GWAS for labeled Manhattan plot")
  }
  invisible()
}

plot_gene_labeled_manhattan <- function(association, labels, analysis, model, filename) {
  require_pkgs("ggrepel")
  d <- association |>
    filter(is.finite(P), P > 0, CHR >= 1, CHR <= 23) |>
    arrange(CHR, BP)
  threshold <- .05 / nrow(d)
  chromosome_layout <- d |>
    group_by(CHR) |>
    summarise(Chromosome_length = as.numeric(max(BP)), .groups = "drop") |>
    arrange(CHR) |>
    mutate(Offset = lag(cumsum(Chromosome_length), default = 0), Center = Offset + Chromosome_length / 2)
  d <- d |>
    left_join(chromosome_layout[, c("CHR", "Offset")], by = "CHR") |>
    mutate(Cumulative_BP = as.numeric(BP) + Offset, Minus_log10_P = -log10(P), Chromosome_group = factor(CHR %% 2))
  label_data <- d |>
    inner_join(labels[, c("SNP", "Gene"), drop = FALSE], by = "SNP") |>
    filter(!is.na(Gene), nzchar(Gene))
  has_labels <- nrow(label_data) > 0
  subtitle <- if (has_labels) {
    sprintf("%d Bonferroni-significant lead loci labeled; dashed line is p = %.3g", nrow(label_data), threshold)
  } else {
    sprintf("No Bonferroni-significant loci; no gene labels are valid (threshold p = %.3g)", threshold)
  }
  p <- ggplot(d, aes(Cumulative_BP, Minus_log10_P, colour = Chromosome_group)) +
    geom_point(size = .75, alpha = .72) +
    geom_hline(yintercept = -log10(threshold), linetype = "dashed", colour = "#C0392B", linewidth = .55) +
    scale_colour_manual(values = c("#2962A3", "#55A6A6"), guide = "none") +
    scale_x_continuous(breaks = chromosome_layout$Center, labels = ifelse(chromosome_layout$CHR == 23, "X", chromosome_layout$CHR), expand = expansion(mult = c(.01, .02))) +
    labs(title = paste0(analysis, " ", model, " GWAS — gene-labeled Manhattan plot"), subtitle = subtitle, x = "Chromosome", y = expression(-log[10](italic(P)))) +
    theme_report() +
    theme(panel.grid.major.x = element_blank())
  if (has_labels) {
    p <- p + ggrepel::geom_text_repel(
      data = label_data,
      aes(Cumulative_BP, Minus_log10_P, label = Gene),
      inherit.aes = FALSE, seed = 42, size = 3, box.padding = .35,
      point.padding = .2, min.segment.length = 0, max.overlaps = Inf,
      segment.colour = "#6B7280", colour = "#172B4D"
    )
  }
  save_fig(p, filename, 12, 7)
  data.frame(
    Analysis = analysis, Model = model, Tests = nrow(d),
    Bonferroni_threshold = threshold,
    Significant_variants = sum(d$P < threshold),
    Gene_labels = nrow(label_data), File = basename(filename)
  )
}

phase6 <- function() {
  out <- out_dir(6)
  build <- read_genome_build()
  ensure_plot_associations()
  association_files <- file.path(out_dir(4), paste0("GWAS_", c("PC1", "PC2"), ".assoc.linear"))
  lead_path <- file.path(out, "annotation_input_lead_loci.csv")
  if (all(file.exists(association_files))) {
    gwas <- bind_rows(lapply(c("PC1", "PC2"), function(pc) {
      mutate(read_assoc(file.path(out_dir(4), paste0("GWAS_", pc, ".assoc.linear"))), Analysis = pc)
    }))
    threshold <- .05 / nrow(filter(gwas, Analysis == "PC1"))
    leads <- gwas |>
      filter(P < threshold) |>
      group_by(Analysis, CHR) |>
      slice_min(P, n = 1, with_ties = FALSE) |>
      ungroup() |>
      select(Analysis, SNP, CHR, BP, P, BETA)
    write.csv(leads, lead_path, row.names = FALSE)
  } else if (file.exists(lead_path)) {
    leads <- read.csv(lead_path, stringsAsFactors = FALSE)
    required <- c("Analysis", "SNP", "CHR", "BP", "P", "BETA")
    if (!all(required %in% names(leads))) stop("Existing lead-locus file is missing required columns.")
    leads <- leads[, required]
  } else {
    stop("Phase 4 association files and the derived Phase 6 lead-locus input are both missing.")
  }

  annotations <- fetch_vep_annotations(leads$SNP)
  annotated <- leads |>
    rename(Source_CHR = CHR, Source_BP = BP) |>
    mutate(Source_assembly = build$UCSC_name) |>
    left_join(annotations, by = "SNP")
  annotated$Mapping_type <- ifelse(is.na(annotated$Gene) | !nzchar(annotated$Gene), "unmapped", "VEP_consequence_gene")
  annotated$Distance_to_gene_bp <- ifelse(annotated$Mapping_type == "VEP_consequence_gene", 0, NA_real_)
  missing_gene <- which(is.na(annotated$Gene) | !nzchar(annotated$Gene))
  for (i in missing_gene) {
    nearest <- fetch_nearest_gene(annotated$Annotation_CHR[i], annotated$Annotation_BP[i])
    if (!is.na(nearest$Gene[1])) {
      annotated$Gene[i] <- nearest$Gene[1]
      annotated$Ensembl_gene_id[i] <- nearest$Ensembl_gene_id[1]
      annotated$Gene_biotype[i] <- nearest$Gene_biotype[1]
      annotated$Mapping_type[i] <- "nearest_protein_coding_gene"
      annotated$Distance_to_gene_bp[i] <- nearest$Distance_to_gene_bp[1]
    }
  }
  if (anyNA(annotated$Annotation_BP)) {
    stop("VEP did not return a current placement for: ", paste(annotated$SNP[is.na(annotated$Annotation_BP)], collapse = ", "))
  }
  write.csv(annotated, file.path(out, "annotated_lead_genes.csv"), row.names = FALSE, na = "")

  plot_summary <- bind_rows(lapply(c("PC1", "PC2"), function(pc) {
    association <- read_assoc(file.path(out_dir(4), paste0("GWAS_", pc, ".assoc.linear")))
    labels <- annotated |> filter(Analysis == pc)
    plot_gene_labeled_manhattan(
      association, labels, pc, "linear",
      file.path(out, paste0(pc, "_gene_labeled_manhattan.png"))
    )
  }))
  sex_association <- read_assoc(file.path(out_dir(5), "GWAS_Sex.assoc.logistic"))
  plot_summary <- bind_rows(
    plot_summary,
    plot_gene_labeled_manhattan(
      sex_association, data.frame(SNP = character(), Gene = character()),
      "Sex", "logistic", file.path(out, "Sex_gene_labeled_manhattan.png")
    )
  )
  write.csv(plot_summary, file.path(out, "manhattan_plot_summary.csv"), row.names = FALSE)

  summary <- annotated |>
    group_by(Analysis) |>
    summarise(
      Lead_loci = n(),
      Loci_with_gene_symbol = sum(!is.na(Gene) & nzchar(Gene)),
      Unique_gene_symbols = n_distinct(Gene[!is.na(Gene) & nzchar(Gene)]),
      Intergenic_loci = sum(Consequence == "intergenic_variant", na.rm = TRUE),
      .groups = "drop"
    )
  write.csv(summary, file.path(out, "annotation_summary.csv"), row.names = FALSE)
  writeLines(c(
    sprintf("Completed: %d lead loci annotated through Ensembl VEP.", nrow(annotated)),
    sprintf("Source coordinates: %s (%s), documented in data/genome_build.tsv.", build$UCSC_name, build$NCBI_assembly),
    sprintf("Annotation coordinates/consequences: %s.", paste(unique(na.omit(annotated$Annotation_assembly)), collapse = ", ")),
    "The rsID is used to remap old source coordinates to the current assembly; source and annotation coordinates are both retained.",
    "Final PC1/PC2 linear and sex-logistic Manhattan plots are written here; the null sex scan has no valid gene labels.",
    "PC1/PC2 associations are descriptive because the PCs were calculated from the same genotypes."
  ), file.path(out, "annotation_status.txt"))
}

phase7 <- function() {
  out <- out_dir(7)
  input <- file.path(out_dir(6), "annotated_lead_genes.csv")
  if (!file.exists(input)) stop("Run Phase 6 first: annotated_lead_genes.csv is missing.")
  require_pkgs(c("clusterProfiler", "org.Hs.eg.db"))
  annotations <- read.csv(input, stringsAsFactors = FALSE)
  if (!all(c("Analysis", "Gene") %in% names(annotations))) stop("Phase 6 annotation output lacks Analysis/Gene columns.")
  annotations <- annotations |> filter(!is.na(Gene), nzchar(Gene))
  if (!nrow(annotations)) stop("Phase 6 produced no mappable gene symbols.")

  mapping <- suppressMessages(clusterProfiler::bitr(
    unique(annotations$Gene), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db
  ))
  gene_map <- annotations |>
    distinct(Analysis, SNP, Gene) |>
    left_join(mapping, by = c("Gene" = "SYMBOL"))
  write.csv(gene_map, file.path(out, "enrichment_gene_mapping.csv"), row.names = FALSE, na = "")

  result_list <- list()
  summary_list <- list()
  for (analysis_name in unique(gene_map$Analysis)) {
    entrez <- gene_map |>
      filter(Analysis == analysis_name, !is.na(ENTREZID)) |>
      pull(ENTREZID) |>
      unique()
    if (length(entrez) < 3) {
      summary_list[[analysis_name]] <- data.frame(
        Analysis = analysis_name,
        Annotated_symbols = sum(gene_map$Analysis == analysis_name),
        Mapped_Entrez_genes = length(entrez),
        Tested_GO_terms = 0,
        BH_significant_terms = 0
      )
      next
    }
    ego <- clusterProfiler::enrichGO(
      gene = entrez, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
      pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
      minGSSize = 3, readable = TRUE
    )
    table <- as.data.frame(ego)
    if (nrow(table)) {
      table$Analysis <- analysis_name
      result_list[[analysis_name]] <- table
    }
    summary_list[[analysis_name]] <- data.frame(
      Analysis = analysis_name,
      Annotated_symbols = n_distinct(gene_map$Gene[gene_map$Analysis == analysis_name]),
      Mapped_Entrez_genes = length(entrez),
      Tested_GO_terms = nrow(table),
      BH_significant_terms = if (nrow(table)) sum(table$p.adjust < .05, na.rm = TRUE) else 0
    )
  }
  results <- bind_rows(result_list)
  summary <- bind_rows(summary_list)
  significant <- if (nrow(results)) filter(results, p.adjust < .05) else results
  write.csv(results, file.path(out, "go_enrichment_all.csv"), row.names = FALSE)
  write.csv(significant, file.path(out, "go_enrichment_significant.csv"), row.names = FALSE)
  write.csv(summary, file.path(out, "enrichment_summary.csv"), row.names = FALSE)
  pathway_table <- significant |>
    transmute(
      Analysis,
      Pathway_ID = ID,
      Pathway_Name = Description,
      Gene_Count = Count,
      Adjusted_P_Value = p.adjust
    )
  write.csv(pathway_table, file.path(out, "enriched_pathways.csv"), row.names = FALSE)
  interpretation <- if (nrow(pathway_table)) {
    sprintf(
      "%d GO Biological Process pathways passed BH-adjusted p < 0.05. Because the input loci are descriptive associations with PCs derived from the same genotypes, biological plausibility should be treated as exploratory rather than as independent phenotype validation.",
      nrow(pathway_table)
    )
  } else {
    "No GO Biological Process pathway passed BH-adjusted p < 0.05 for either PC1 or PC2. Therefore, there is no statistically supported enriched pathway to interpret as biologically plausible for population structure. The leading terms shown in the dot plot are exploratory and non-significant; moreover, the PC associations are descriptive because the PCs were calculated from the same genotype matrix."
  }
  writeLines(interpretation, file.path(out, "pathway_interpretation.txt"))

  if (nrow(results)) {
    plot_data <- results |>
      group_by(Analysis) |>
      slice_min(p.adjust, n = 10, with_ties = FALSE) |>
      ungroup() |>
      mutate(Description = reorder(Description, -log10(pmax(p.adjust, .Machine$double.xmin))))
    p <- ggplot(plot_data, aes(-log10(pmax(p.adjust, .Machine$double.xmin)), Description, size = Count, colour = p.adjust)) +
      geom_point(alpha = .85) +
      facet_wrap(~Analysis, scales = "free_y") +
      scale_colour_viridis_c(direction = -1) +
      labs(
        title = "GO biological-process enrichment",
        subtitle = "Top terms per descriptive PC axis; colour is BH-adjusted p-value",
        x = expression(-log[10]("BH-adjusted p-value")), y = NULL
      ) +
      theme_report()
    save_fig(p, file.path(out, "go_enrichment_dotplot.png"), 10, 8)
  }
  total_significant <- sum(summary$BH_significant_terms)
  writeLines(c(
    sprintf("Completed GO Biological Process over-representation analysis for %d PC axis/axes.", nrow(summary)),
    sprintf("%d terms pass BH-adjusted p < 0.05.", total_significant),
    "Gene lists come only from Phase 6 Ensembl annotations; no hard-coded pathway genes are used.",
    "Interpretation is exploratory: PC association is circular/descriptive and per-axis gene lists are small."
  ), file.path(out, "enrichment_status.txt"))
}
