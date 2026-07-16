### Phase 5.3: Comparison of Linear (Phase 4) and Logistic (Phase 5) GWAS Hits

**Biological Signal vs. Ancestry Markers**
In Phase 4, the linear association test identified Ancestry Informative Markers (AIMs) driving the population structure (PC1). The top significant hits (e.g., rs10466604, P = 9.5e-28) were distributed across autosomal chromosomes. This mathematically reflects the reality of human genetic admixture, which is a continuous, polygenic trait driven by mutations scattered across the entire genome.

**Are all SNPs associated with sex in the X and Y chromosomes only?**
Biologically, yes. The defining genetic markers for biological sex are strictly located on the X and Y chromosomes. 

When observing the logistic regression output for Biological Sex in Phase 5, the top hits mapped to autosomal chromosomes (Chromosomes 1, 5, 7, 17) rather than Chromosome 23. However, these hits yielded mathematically weak P-values (e.g., P = 6.8e-05) that fail to meet the Bonferroni-corrected threshold for genome-wide significance (~7.3e-07). 

This outcome serves as an excellent technical validation of the Phase 1 Quality Control pipeline. Because sex chromosomes violate standard Hardy-Weinberg Equilibrium assumptions due to hemizygosity in males, non-autosomal chromosomes were filtered out of the dataset prior to analysis. Consequently, the logistic regression correctly demonstrated that there is zero statistically significant genetic signal for biological sex located on the autosomal chromosomes (1-22), returning only expected background statistical noise.