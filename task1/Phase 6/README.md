# Phase 6: Gene and Functional Annotation

Phase 6 selects the strongest Bonferroni-significant locus on each chromosome for each
descriptive PC association scan, then annotates those observed rsIDs with Ensembl VEP.

The original PLINK coordinates are NCBI36/hg18. This was confirmed by eight sampled rsID
and position matches against the UCSC hg18 dbSNP130 track and is recorded in
`../data/genome_build.tsv`. Ensembl remaps each stable rsID to GRCh38 for current gene and
consequence annotation. Both source and current coordinates are retained in the output.
If VEP provides no named consequence gene, the script finds the nearest GRCh38
protein-coding gene within 1 Mb and records the mapping type and distance.

Run from `task1`:

```bash
Rscript professional_pipeline.R 6
```

Current result: all 45 lead loci were mapped to named genes. Of these, 27 use a direct VEP
The final Manhattan deliverables include gene-labeled PC1 and PC2 linear plots plus a
sex-logistic plot. The logistic scan has no Bonferroni-significant variants, so its figure
correctly has no gene labels and states that null result.
consequence gene and 18 use an explicitly labelled nearest protein-coding gene.

Deliverables:

- `outputs/annotation_input_lead_loci.csv`: reproducible lead-locus selection.
- `outputs/annotated_lead_genes.csv`: source/current coordinates, gene, mapping type,
  distance, transcript, consequence, impact, and biotype.
- `outputs/annotation_summary.csv`: annotation coverage by PC axis.
- `outputs/PC1_gene_labeled_manhattan.png`: labeled PC1 linear-GWAS Manhattan plot.
- `outputs/PC2_gene_labeled_manhattan.png`: labeled PC2 linear-GWAS Manhattan plot.
- `outputs/Sex_gene_labeled_manhattan.png`: final logistic-GWAS Manhattan plot; no labels
  because there are no significant loci.
- `outputs/manhattan_plot_summary.csv`: tests, thresholds, significant variants, and label counts.
- `outputs/annotation_status.txt`: assemblies, method, and interpretation safeguard.

These are associations with PCs calculated from the same genotype matrix. They describe
population structure and are not independent phenotype discoveries.
