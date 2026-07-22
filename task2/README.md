# Task 2: relatedness and metabolomics

Use the scripts and interpretation rules in [`../RUNBOOK.md`](../RUNBOOK.md).

The current source-of-truth entry points are the existing phase `.R` files, which dispatch to:

- `professional_pipeline.R` for Phases 11--14;
- `relatedness_pipeline.R` for Phases 8--10; and
- `future_phases.R` for Phases 15--16.

Important safeguards: the GRM is not treated as a KING table; PC-GWAS results remain descriptive; the classifier is only internally validated; and metabolite database annotations are not inferred from free-text names.
