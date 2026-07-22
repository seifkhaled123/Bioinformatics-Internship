# Reproducible analysis runbook

The legacy phase scripts are retained as entry points, but now dispatch to the professional pipelines. Run any phase script directly, or use the commands below from the corresponding task directory.

## Task 1

```bash
cd task1
Rscript professional_pipeline.R 1  # QC and filter waterfall
Rscript professional_pipeline.R 2  # LD-pruned PCA
Rscript professional_pipeline.R 3  # silhouette-selected clustering
Rscript professional_pipeline.R 4  # covariate-adjusted PC association
Rscript professional_pipeline.R 5  # PC-adjusted sex association
Rscript professional_pipeline.R 6  # lead-locus annotation input
Rscript professional_pipeline.R 7  # enrichment only from observed annotations
```

PC association is explicitly descriptive: PCs are constructed from the same genotypes. Do not present its lead SNPs as independent trait discoveries. Phase 6 deliberately requires a confirmed genome build before annotation; Phase 7 never uses a hard-coded gene list.

## Task 2

```bash
cd task2
Rscript relatedness_pipeline.R 8
Rscript relatedness_pipeline.R 9
HERITABILITY_TRAIT=Glucose Rscript relatedness_pipeline.R 10
Rscript professional_pipeline.R 11
Rscript professional_pipeline.R 12
Rscript professional_pipeline.R 13
Rscript professional_pipeline.R 14
```

Phase 9 requires PLINK 2 for the KING-robust table. Set `PLINK2_BIN` to its executable path. The GRM is generated separately and is never interpreted using KING thresholds. Phase 10 should run only after duplicate handling is decided from the KING output; report heritability with a model-derived standard error or confidence interval.

Phase 13 uses a fixed stratified holdout, training-only scaling, elastic-net logistic regression, and random forest. Its 30-person test set is small: report it as internal validation, not as a final clinical-performance estimate.

Phase 14 is an exploratory ridge-regularised partial-correlation network. The exported Cytoscape edge/node tables include weights, signs, degree, and component membership. Interpret hubs only after stability analysis or external replication.

## Output convention

Each phase writes human-readable CSV/TXT tables, 320-dpi PNG figures, model/configuration files, and (for Task 2) `sessionInfo.txt` into that phase's `outputs/` directory. Existing historical outputs are preserved; files with the new descriptive names are the current deliverables.
