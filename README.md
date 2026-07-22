# Bioinformatics internship: integrated genomic and metabolomic analysis

## Executive summary

This repository contains two connected analyses of 156 participants:

1. **Task 1:** genotype quality control, ancestry structure, clustering, and descriptive genome-wide association analyses.
2. **Task 2:** metabolite quality control, diabetes associations, prediction, metabolite network analysis, genetic relatedness/heritability, and mQTL analysis.

The original phase scripts are retained as entry points, but dispatch to reproducible professional pipelines. Results are written inside each phase's `outputs/` directory; the links and figures below point to the current deliverables.

### Headline findings

| Area | Result | Interpretation |
|---|---|---|
| Genotype QC | 156 participants; 67,735 starting variants; 55,226 after strict missingness and 31,694 after strict MAF | The raw genotype data have low missingness; strict filters are supplied as sensitivity-ready subsets. |
| Population structure | PC1 explains 27.35% and PC2 13.68% of LD-pruned variation | Genetic structure is substantial and must be accounted for in association models. |
| Relatedness | 12,084 unrelated pairs, 6 third-degree pairs, no first-/second-degree or duplicate pairs | No close-relative removal was required. |
| Diabetes–metabolite associations | 60 FDR and 33 Bonferroni-significant metabolites | Strongest association: lower Anhydroglucitol_1_5 in diabetes (beta −1.38, P = 2.79e−20). |
| Diabetes prediction | Both elastic net and random forest: held-out balanced accuracy 0.955, AUC 1.00 | Promising internal result, but the 30-person holdout is too small for clinical claims. |
| Metabolite network | 328 ridge-regularised partial-correlation edges among 136 metabolites | Exploratory network; hub/edge claims need stability analysis or replication. |
| Glucose heritability | h² = 0 with a converged boundary estimate | This sample/GRM model finds no detectable additive genetic component for Glucose. |
| mQTL scan | 33 traits × 66,078 autosomal variants in naïve and GRM-adjusted models | One within-trait FDR-significant adjusted result (Glutamate); treat as discovery-stage evidence. |

## Data, software, and reproducibility

The study includes PLINK genotype data, participant mapping information, metabolomics data, and diabetes/sex/PC covariates. The analyses use R, PLINK 1.9, and PLINK 2 available on the system PATH (or through `PLINK2_BIN`). The mQTL analysis uses GENESIS, SeqArray, SeqVarTools, and the genomic relationship matrix (GRM).

Run commands and software assumptions are in [RUNBOOK.md](RUNBOOK.md). Each Task 2 output folder also contains `sessionInfo.txt` where relevant.

### Important interpretation principles

- Association does not establish causation.
- Population PCs were built from the same genotypes. PC-GWAS is therefore **descriptive of genetic structure**, not an independent disease-trait discovery study.
- Genome build information was not confirmed in the supplied material. Variant-to-gene and metabolite-database annotations are deliberately not fabricated.
- Multiple-testing results are explicitly labelled as Bonferroni, BH/FDR, or within-trait FDR; these are not interchangeable.
- The small sample size (n = 156) limits power and external generalisability.

---

# Task 1 — genotype and population-structure analysis

## Phase 1 — genotype quality control

**Question:** Are the genotype data suitable for downstream analysis, and what is the effect of reasonable QC thresholds?

**Method and decisions**

- Calculated SNP and participant missingness plus allele frequencies.
- Preserved the baseline dataset, then generated strict sensitivity subsets rather than silently deleting data.
- Standard QC: 156 participants and 67,735 variants retained.
- Strict missingness: 55,226 variants retained.
- Strict missingness plus MAF filtering: 31,694 variants retained.

The observed maximum SNP missingness was 0.00641 and maximum sample missingness was 0.01088, so no participant failed the chosen thresholds.

**Outputs:** [filter waterfall](task1/Phase%201/outputs/qc_filter_waterfall.csv), [QC summary](task1/Phase%201/outputs/qc_summary.csv), [MAF distribution](task1/Phase%201/outputs/maf.png), [sample missingness](task1/Phase%201/outputs/sample_missingness.png), and [SNP missingness](task1/Phase%201/outputs/snp_missingness.png).

![MAF distribution](task1/Phase%201/outputs/maf.png)

## Phase 2 — LD-pruned principal-component analysis

**Question:** What genetic structure is present, and which covariates should control it?

**Method and decisions**

- Performed explicit LD pruning before PCA to avoid regions of correlated variants dominating the components.
- Retained PC1–PC10; PC1–PC5 are used as covariates in later adjusted analyses.
- PC1 explains 27.35%, PC2 13.68%, and PC3 9.51% of the LD-pruned genetic variation.

**Outputs:** [variance table](task1/Phase%202/outputs/pca_variance_explained.csv), [sample PC scores](task1/Phase%202/outputs/pca_scores.csv), [scree plot](task1/Phase%202/outputs/pca_scree_plot.png), and [PC1/PC2 plot](task1/Phase%202/outputs/pca_pc1_pc2.png).

![PCA](task1/Phase%202/outputs/pca_pc1_pc2.png)

## Phase 3 — unsupervised genetic clustering

**Question:** Does the PCA structure contain discrete clusters?

**Method and decisions**

- Evaluated k = 2–8 with within-cluster sum of squares and mean silhouette width.
- The maximum silhouette width was at k = 5 (0.582), narrowly above k = 4 (0.581); k = 5 was selected transparently rather than chosen by visual appearance alone.
- Cluster sizes were 76, 21, 5, 38, and 16. The five-person cluster warrants cautious interpretation.

**Outputs:** [diagnostics](task1/Phase%203/outputs/cluster_diagnostics.csv), [assignments](task1/Phase%203/outputs/cluster_assignments.csv), [silhouette selection](task1/Phase%203/outputs/cluster_silhouette_selection.png), and [cluster plot](task1/Phase%203/outputs/pca_clusters_2d.png).

![Genetic clusters](task1/Phase%203/outputs/pca_clusters_2d.png)

## Phase 4 — PC association analysis (descriptive)

**Question:** Which variants track the major axes of genetic variation?

**Method and decisions**

- Ran PLINK linear association for PC1 and PC2.
- Adjusted each model for the other leading PCs to reduce direct confounding by the same broad structure.
- This is not independent GWAS discovery: the phenotype (PC) is constructed from the genotype data.

**Result:** PC1 had lambda GC 4.61 and 2,485 Bonferroni-significant variants; PC2 had lambda GC 2.55 and 118. The high inflation is consistent with the intentionally structure-derived phenotype, so these results are used for structure/LD follow-up only. The leading PC1 variant was `rs7355960` at chr3:180,566,740 (P = 3.83e−28).

**Outputs:** [model summary](task1/Phase%204/outputs/gwas_summary.csv), [PC1 top hits](task1/Phase%204/outputs/PC1_top_hits.csv), [PC1 Manhattan](task1/Phase%204/outputs/PC1_manhattan.jpg), [PC1 QQ](task1/Phase%204/outputs/PC1_qq.jpg), and equivalent PC2 files.

![PC1 Manhattan plot](task1/Phase%204/outputs/PC1_manhattan.jpg)

## Phase 5 — sex association scan

**Question:** After controlling for structure, are there autosomal genotype associations with recorded sex?

**Method and decisions**

- Verified sex coding first.
- Used logistic association adjusted for PC1–PC5.
- Kept this analysis separate from the PC-derived analyses because its phenotype is externally recorded.

**Result:** lambda GC = 0.989 and no Bonferroni-significant variants among 67,718 tests. This well-calibrated null result is preferable to overinterpreting nominal hits.

**Outputs:** [summary](task1/Phase%205/outputs/sex_gwas_summary.csv), [coding check](task1/Phase%205/outputs/sex_coding_check.csv), [Manhattan plot](task1/Phase%205/outputs/Sex_manhattan.jpg), and [QQ plot](task1/Phase%205/outputs/Sex_qq.jpg).

## Phases 6–7 — annotation and enrichment safeguards

**Decision:** No gene or pathway enrichment result is reported until the genome build and actual variant-to-gene mapping are confirmed. Phase 6 writes the reproducible lead-locus input and Phase 7 refuses hard-coded gene lists.

**Outputs:** [annotation input](task1/Phase%206/outputs/annotation_input_lead_loci.csv), [annotation status](task1/Phase%206/outputs/annotation_status.txt), and [enrichment status](task1/Phase%207/outputs/enrichment_status.txt).

---

# Task 2 — integrated metabolomic and genetic analysis

## Phase 8 — LD around the lead locus

**Question:** What local LD structure surrounds the strongest descriptive PC1 hit?

**Method:** Calculated the full r² matrix in a ±500 kb window around `rs7355960` (chr3:180,566,740), containing 19 variants.

**Outputs:** [region summary](task2/Phase%208/outputs/ld_region_summary.csv), [r² matrix](task2/Phase%208/outputs/ld_r2_matrix.csv), and [LD heatmap](task2/Phase%208/outputs/ld_heatmap.png).

![LD heatmap](task2/Phase%208/outputs/ld_heatmap.png)

## Phase 9 — relatedness and GRM

**Question:** Are close relatives or duplicate samples present, and what covariance matrix is appropriate for mixed models?

**Method and decisions**

- Used PLINK 2 `--make-king-table` for KING-robust pairwise relatedness.
- Used a separate PLINK 1.9 GRM for mixed models. A GRM is not labelled or thresholded as KING.
- Applied degree cut-offs of 0.0442, 0.0884, 0.177, and 0.354.

**Result:** Of 12,090 pairs, 12,084 are unrelated and 6 are third-degree; there are no closer pairs or duplicates. Therefore no duplicate/removal decision was needed before downstream analyses.

**Outputs:** [KING category summary](task2/Phase%209/outputs/king_relatedness_summary.csv), [related pairs](task2/Phase%209/outputs/king_related_pairs.csv), [affected participant list](task2/Phase%209/outputs/individuals_in_related_pairs.csv).

## Phase 10 — metabolite heritability

**Question:** What fraction of Glucose variance is explained by the GRM in this sample?

**Method:** Fitted a `sommer` variance-component mixed model with the Phase 9 GRM. The trait can be changed with `HERITABILITY_TRAIT`.

**Result:** For Glucose (n = 156), genetic variance was estimated at the non-negative boundary (0), residual variance was 0.863, and h² = 0. The model converged. Because this is a boundary estimate, the delta-method SE is not available; this is a valid null/boundary result, not evidence that heritability is universally zero.

**Outputs:** [heritability result](task2/Phase%2010/outputs/heritability_result.csv), [diagnostics](task2/Phase%2010/outputs/heritability_diagnostics.csv), and [variance components](task2/Phase%2010/outputs/variance_components.png).

## Phase 11 — metabolite QC and preprocessing

**Question:** Are metabolite measurements complete and prepared consistently?

**Method and decisions**

- Checked sample and metabolite missingness before filtering.
- Used threshold-first QC followed by median imputation only if needed.
- Saved clean, imputed, and scaled analysis inputs with the QC metrics.

**Result:** all 156 participants and all 136 metabolites passed; no values required imputation.

**Outputs:** [QC summary](task2/Phase%2011/outputs/qc_summary.csv), [metabolite metrics](task2/Phase%2011/outputs/metabolite_qc_metrics.csv), [cleaned data](task2/Phase%2011/outputs/Qatari_metabolomics_Cleaned_Imputed.csv), and [missingness figure](task2/Phase%2011/outputs/missingness_distribution.png).

## Phase 12 — diabetes–metabolite associations

**Question:** Which metabolites differ with diabetes status after accounting for sex and population structure?

**Method:** One linear model per metabolite with diabetes, sex, and PC1–PC5 as predictors. Effects have confidence intervals and BH/FDR plus Bonferroni correction.

**Results:** 60 metabolites were FDR-significant and 33 Bonferroni-significant. The five strongest reported effects were:

| Metabolite | Direction in diabetes | Beta | P value |
|---|---:|---:|---:|
| Anhydroglucitol_1_5 | lower | −1.375 | 2.79e−20 |
| Mannose | higher | 1.230 | 1.35e−15 |
| Citrulline | lower | −1.082 | 1.14e−11 |
| Enyl_palmitoyl_GPC | lower | −0.998 | 5.70e−10 |
| Palmitoylcholine | lower | −0.954 | 2.77e−09 |

**Outputs:** [complete association table](task2/Phase%2012/outputs/metabolite_diabetes_associations.csv), [model specification](task2/Phase%2012/outputs/model_specification.csv), [forest plot](task2/Phase%2012/outputs/top_associations_forest.png), [volcano plot](task2/Phase%2012/outputs/volcano_plot.png), and [top-metabolite plot](task2/Phase%2012/outputs/top_metabolites_by_diabetes.png).

![Diabetes associations](task2/Phase%2012/outputs/top_associations_forest.png)

## Phase 13 — diabetes classification

**Question:** How well can the metabolite panel discriminate diabetes status internally?

**Method and decisions**

- Used a fixed stratified holdout split (training-only scaling) to avoid data leakage.
- Compared elastic-net logistic regression with random forest.
- Saved held-out predictions, ROC curve, and feature importance.

**Result:** both models had accuracy 0.967, sensitivity 1.00, specificity 0.909, balanced accuracy 0.955, and AUC 1.00 on the 30-person holdout. This is an internal result only; it may be optimistic and requires cross-validation/external validation.

**Outputs:** [metrics](task2/Phase%2013/outputs/classification_metrics.csv), [held-out predictions](task2/Phase%2013/outputs/test_set_predictions.csv), [ROC](task2/Phase%2013/outputs/test_set_roc.png), and [feature importance](task2/Phase%2013/outputs/metabolite_importance.png).

## Phase 14 — metabolite partial-correlation network

**Question:** Which conditional metabolite relationships remain after accounting for all other measured metabolites?
