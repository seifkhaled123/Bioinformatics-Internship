# Phase 7: Pathway Enrichment Analysis

Phase 7 performs GO Biological Process over-representation analysis separately for the
PC1 and PC2 gene lists produced by Phase 6. Gene symbols are mapped to Entrez IDs through
`org.Hs.eg.db`; no hard-coded gene list or pathway result is used.

Run from `task1` after Phase 6:

```bash
Rscript professional_pipeline.R 7
```

Current result: all 23 PC1 and 22 PC2 gene symbols mapped to Entrez IDs. The analysis
tested 363 GO terms for PC1 and 389 for PC2. No term passed BH-adjusted p < 0.05 (the
smallest adjusted p-value was approximately 0.128), so no enriched pathway is claimed.

Deliverables:

- `outputs/enriched_pathways.csv`: requested pathway name, gene count, and adjusted
  p-value table; currently header-only because no pathway is significant.
- `outputs/enrichment_gene_mapping.csv`: auditable symbol-to-Entrez mapping.
- `outputs/go_enrichment_all.csv`: all tested terms.
- `outputs/go_enrichment_significant.csv`: terms with BH-adjusted p-value below 0.05;
  currently header-only because none pass.
- `outputs/enrichment_summary.csv`: genes, tested terms, and significant-term counts.
- `outputs/go_enrichment_dotplot.png`: top terms per PC axis, clearly shown with adjusted
  p-values even though they are not significant.
- `outputs/enrichment_status.txt`: concise completion and interpretation status.
- `outputs/pathway_interpretation.txt`: standalone biological-plausibility paragraph.

The analysis is exploratory because each input gene list is small and comes from
descriptive associations with PCs derived from the same genotypes. A plausible-looking
pathway is not treated as a validated phenotype association.
