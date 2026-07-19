## Phase 9: Kinship Matrix and Relatedness Structure

**Objective:** Compute a Genomic Relationship Matrix (GRM) to identify cryptic relatedness, duplicate samples, and family structures within the cohort prior to running mixed-model association analyses.

### Methodology
Pairwise genetic relatedness was estimated using `PLINK2`'s `--make-rel square` function, which calculates the kinship coefficient based on shared allele frequencies across all QC-passing autosomes. The resulting $156 \times 156$ matrix was parsed in R to visualize the relatedness structure and quantify biological relationships based on standard kinship thresholds.

### Deliverable 1: Kinship Heatmap
*(Insert/Attach your Rplots.pdf image here)*

### Deliverable 2: Relatedness Report

| Relatedness Category | Mathematical Threshold (Kinship) | Number of Pairs |
| :--- | :--- | :--- |
| **Identical Twins / Duplicates** | $> 0.354$ | 2 |
| **1st Degree** (Siblings/Parent-Child) | $0.177 < r \le 0.354$ | 45 |
| **2nd Degree** (Half-Sib/Uncle-Aunt) | $0.088 < r \le 0.177$ | 87 |
| **3rd Degree** (First Cousins) | $0.050 < r \le 0.088$ | 143 |

**Total non-unrelated pairs identified in the cohort:** 277

**Note on Sample Exclusion:** We identified 2 pairs of identical twins or duplicate sample entries. To preserve the integrity of downstream predictive models and avoid artificial deflation of standard errors, one individual from each of these duplicate pairs **must be excluded** prior to Phase 11. The remaining related individuals (1st-3rd degree) do not need to be hard-filtered; instead, their covariance will be explicitly modeled by including this GRM as a random effect in our upcoming Mixed-Model GWAS.