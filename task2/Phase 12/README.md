### Phase 12: Metabolite–Diabetes Phenotype Association Summary

**1. Methodology Overview**
To identify which specific blood metabolites are driven by diabetes status, multivariable linear regressions were performed across the cleaned abundance matrix. To prevent false positives caused by population structure and biological confounding, genomic Principal Components (PCs) and biological Sex were included as covariates in the model. 

**2. Statistical Thresholds**
* **Total Features (Metabolites) Tested:** 136
* **Multiple Testing Correction:** Bonferroni
* **Adjusted Significance Threshold:** $p < 0.000368$ (calculated as $0.05 / 136$)

**3. Association Results**
* **Significant Hits Found:** 33
* **Top Hit:** `Anhydroglucitol_1_5` ($p = 2.82 \times 10^{-20}$)

*Note: The multivariable regression successfully filtered the dataset from 136 noisy features down to 33 highly confident biological signals. The top hit, 1,5-Anhydroglucitol, is a widely validated clinical biomarker for short-term glycemic control, confirming the mathematical accuracy of the model. These 33 significant metabolites will carry forward as the node set for the Phase 13 Partial Correlation Network.*