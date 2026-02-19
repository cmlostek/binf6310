# BINF 6310 Final Project Proposal

**Title:** Bayesian Estimation of Allele Frequencies Across Human Populations: Detecting Signatures of Natural Selection Using the Beta-Binomial Model

---

## Background & Significance

Natural selection leaves detectable signatures in the genome by shifting allele frequencies in ways that deviate from neutral evolutionary expectations. Identifying these loci is fundamental to understanding human adaptation, disease susceptibility, and population history. Classical approaches (e.g., FST-based outlier tests) treat allele frequencies as point estimates, discarding the uncertainty inherent in finite sample sizes. A Bayesian framework, using the Beta distribution as a conjugate prior over allele frequencies, offers a principled alternative: posterior distributions explicitly capture uncertainty and enable coherent comparison across populations with different sample sizes.

The 1000 Genomes Project Phase 3 provides whole-genome sequencing data for 2,504 individuals across 26 populations and 5 continental groups, making it an ideal resource for population-level allele frequency analysis.

---

## Biological Question

**Can a Bayesian Beta-Binomial model identify SNP loci whose posterior allele frequency distributions differ significantly across continental populations, recovering known signals of natural selection?**

---

## Specific Aims

### Aim 1 — Bayesian Allele Frequency Estimation
Fit a Beta-Binomial model to allele count data for a curated set of SNPs across five continental population groups (AFR, AMR, EAS, EUR, SAS). For each SNP and population, define a Beta(α, β) prior and update with observed reference/alternate allele counts to obtain a posterior distribution. Compare uninformative (Beta(1,1)) and weakly informative priors derived from genome-wide allele frequency distributions.

### Aim 2 — Population Differentiation via Posterior Comparison
Quantify allele frequency differentiation between populations by comparing posterior Beta distributions using the Bhattacharyya coefficient and Kullback-Leibler divergence. Rank SNPs by the degree of inter-population differentiation to identify candidate loci under selection.

### Aim 3 — Validation Against Known Selected Loci
Assess whether the model recovers well-documented selection signals at biologically validated loci, including:
- *LCT* (rs4988235) — lactase persistence, strong selection in Europeans
- *DARC* (rs2814778) — Duffy antigen, near fixation in Africans
- *EDAR* (rs3827760) — hair/sweat gland morphology, high frequency in East Asians
- *SLC24A5* (rs1426654) — skin pigmentation, strong sweep in Europeans

---

## Methods

### Data
- **Dataset:** 1000 Genomes Project Phase 3 (release 20130502)
- **FTP:** https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/
- **Format:** Per-chromosome VCF files (`.vcf.gz`) with tabix indices (`.tbi`). Population-level allele frequencies are encoded directly in the VCF INFO fields — no separate `.frq` files exist in this release.
- **Scope:** A targeted panel of ~500–1,000 SNPs including the validation loci above and a random background set for comparison

**Relevant VCF INFO fields:**

| Field | Description |
|---|---|
| `AF` | Overall allele frequency |
| `AFR_AF` | African allele frequency |
| `AMR_AF` | Admixed American allele frequency |
| `EAS_AF` | East Asian allele frequency |
| `EUR_AF` | European allele frequency |
| `SAS_AF` | South Asian allele frequency |

**Files needed for validation loci:**

| Gene | Variant | Chromosome | File |
|---|---|---|---|
| *LCT* | rs4988235 | chr2 | `ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz` |
| *DARC* | rs2814778 | chr1 | `ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz` |
| *EDAR* | rs3827760 | chr2 | `ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz` |
| *SLC24A5* | rs1426654 | chr15 | `ALL.chr15.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz` |

Since the VCFs are tabix-indexed, specific variants can be queried without downloading entire chromosome files. Allele frequencies for target SNPs can also be retrieved directly via the **Ensembl REST API** in R using `httr`, requiring no file downloads at all.

### Statistical Model

**Likelihood:**
$$x_i \sim \text{Binomial}(n_i, p_i)$$

Where $x_i$ = alternate allele count and $n_i$ = total allele count in population $i$.

**Prior:**
$$p_i \sim \text{Beta}(\alpha, \beta)$$

**Posterior** (analytically derived via conjugacy):
$$p_i \mid x_i \sim \text{Beta}(\alpha + x_i, \; \beta + n_i - x_i)$$

### Analysis Steps
1. Download and parse allele count data per population in R
2. Fit posterior Beta distributions for each SNP × population combination
3. Visualize posterior distributions and credible intervals across populations
4. Compute pairwise Bhattacharyya coefficients and rank SNPs by differentiation
5. Compare model rankings against known selection loci for validation

### Software & Packages
- `vcfR` — VCF parsing (if using raw genotype data)
- `ggplot2` / `patchwork` — visualization of posterior distributions
- `dplyr` / `tidyr` — data manipulation
- Base R — Beta distribution functions (`dbeta`, `pbeta`, `rbeta`)
- Optional: `rstan` for extended hierarchical models (HPC if needed)

---

## Expected Outcomes

- Posterior allele frequency estimates with credible intervals for each population
- A ranked list of SNPs showing the greatest inter-population differentiation
- Successful recovery of known selection signals, validating the Bayesian approach
- Comparison of uninformative vs. informative prior sensitivity

---

## Data Access

**1000 Genomes Project Phase 3 (release 20130502)**
- FTP: https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/
- Population-level allele frequencies (AFR_AF, AMR_AF, EAS_AF, EUR_AF, SAS_AF) are stored in the INFO fields of per-chromosome VCF files — there are no separate pre-computed `.frq` files in this release
- Population code reference: AFR (African), AMR (Admixed American), EAS (East Asian), EUR (European), SAS (South Asian)
- Sample panel file: `integrated_call_samples_v3.20130502.ALL.panel` (maps individual IDs to populations)

Since all VCFs are tabix-indexed, targeted queries for specific variants or genomic regions are possible without downloading full chromosome files (~200 MB–1.2 GB each). For the validation SNPs only, the Ensembl REST API can return population allele frequencies directly in R with no downloads required. The HPC would only be needed for genome-wide analysis.

---

## Timeline

| Week | Milestone |
|------|-----------|
| 1 | Data download, parsing, and QC in R |
| 2 | Implement Beta-Binomial model; Aim 1 complete |
| 3 | Population differentiation metrics; Aim 2 complete |
| 4 | Validation against known loci; Aim 3 complete |
| 5 | Figures, write-up, final submission |
