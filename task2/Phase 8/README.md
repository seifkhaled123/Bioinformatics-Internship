## Phase 8: Linkage Disequilibrium (LD) Structure Analysis

**Objective:** Characterize the Linkage Disequilibrium (LD) structure and feature correlation between SNPs in the highly significant genomic region identified during the GWAS, and visualize it as a heatmap.

### Methodology
To avoid the computational bottleneck of a genome-wide $N \times N$ correlation matrix, we isolated a manageable 1 Megabase (MB) sliding window centered around our most statistically significant signal (lowest P-value / tallest peak from the Phase 6 Manhattan Plot). Pairwise Pearson correlation squared ($r^2$) was computed for all SNPs within this bounded region using `PLINK`. The resulting mathematical matrix was visualized as a heatmap using R.

### Deliverable 1: LD Heatmap
Located at "outputs/Rplots.PDF"

### Deliverable 2: Observation Note on Haplotype Structure
**Analysis of LD Blocks and Regional Haplotypes:**

Based on the visual observation of the generated heatmap, the feature correlation ($r^2$) is not randomly distributed. Instead, we observe distinct, highly rigid "blocks" of intense correlation (deep red/warm color clusters approaching $r^2 \approx 1.0$) clinging closely to the diagonal axis.

* **LD Blocks:** These distinct geometric clusters represent strong Linkage Disequilibrium blocks. Mathematically, this indicates severe multicollinearity between these adjacent genetic features.
* **Haplotype Structure:** Biologically, these blocks correspond to highly conserved haplotype structures within the study population. Because these adjacent SNPs are inherited together as a single physical package, they form unbroken segments of DNA that have survived generations without being separated.
* **Recombination Hotspots:** The sharp transitions from high correlation (red) to near-zero correlation (white/cool colors) at the edges of these blocks indicate historical recombination hotspots. These are the exact physical coordinates where the DNA string is routinely "cut and shuffled" during meiosis, effectively breaking the feature correlation between adjacent genetic zones.

**Conclusion:** The presence of these clearly defined LD blocks confirms that the strong association signals observed in our GWAS are likely inherited together as part of a single, dominant haplotype in this region, rather than acting as independent mutations.