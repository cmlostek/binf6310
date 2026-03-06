---
output:
  pdf_document: default
  html_document: default
---
# BINF 6310 Final Project Proposal

**Title:** Bayesian Estimation of Allele Frequencies Across Human Populations: Detecting Signatures of Natural Selection Using the Beta-Binomial Model

---

## Background & Significance

Natural selection leaves detectable signatures in the genome by shifting allele frequencies in ways that deviate from neutral evolutionary expectations. Over the course of human prehistory, populations adapting to distinct environments accumulated differences in the frequency of functionally important variants — changes that are now reflected in patterns of disease susceptibility, metabolic traits, and morphological variation. Identifying the genomic loci underlying these adaptations is a central problem in human evolutionary genetics and has direct implications for understanding population-specific disease risk.

Classical approaches to detecting selection, such as FST-based outlier tests, treat allele frequencies as fixed point estimates. This ignores the uncertainty inherent in finite sample sizes: a variant observed in 3 out of 10 sampled chromosomes and one observed in 300 out of 1,000 carry very different amounts of information, yet both yield a point estimate of 0.30. This uncertainty is particularly important when comparing populations with unequal or modest sample sizes, where naive frequency comparisons can produce misleading rankings of candidate loci.

A Bayesian framework offers a principled solution. By modeling allele frequency as a random variable with a prior distribution updated by observed data, we obtain a full posterior distribution over the frequency rather than a single point estimate. This posterior explicitly encodes uncertainty, naturally shrinks extreme estimates from small samples toward the prior, and enables statistically coherent comparisons across populations regardless of sample size differences. For this project, we adopt the Beta-Binomial model, in which the Beta distribution serves as the conjugate prior for the binomial likelihood of observing allele counts. This conjugacy yields analytically tractable posteriors and makes the approach computationally accessible in R without requiring Markov chain Monte Carlo sampling.

---

## Biological Question

Can a Bayesian Beta-Binomial model identify SNP loci whose posterior allele frequency distributions differ significantly across continental populations, recovering known signals of natural selection?

---

## Specific Aims

### Aim 1 — Bayesian Allele Frequency Estimation

For a curated set of SNPs spanning five continental population groups (AFR, AMR, EAS, EUR, SAS), we will fit a Beta-Binomial model to allele count data extracted from the 1000 Genomes Project. Each SNP-population combination will be assigned a Beta(a, B) prior, which will then be updated with the observed reference and alternate allele counts to produce an analytic posterior distribution. We will evaluate both an uninformative uniform prior (Beta(1,1)) and a weakly informative prior derived from genome-wide allele frequency distributions, allowing us to assess prior sensitivity as a secondary outcome.

### Aim 2 — Population Differentiation via Posterior Comparison

Rather than comparing point estimates of allele frequency, we will quantify population differentiation by directly comparing the posterior Beta distributions across continental groups. We will compute pairwise Bhattacharyya coefficients and Kullback-Leibler divergences between distributions, which measure the overlap and information distance between populations, respectively. SNPs will be ranked by the degree of inter-population differentiation identified through these metrics, producing a list of candidate loci showing unusually divergent frequency distributions consistent with population-specific selection.

### Aim 3 — Validation Against Known Selected Loci

As a validation step, we will test whether our model successfully recovers well-documented selection signals at four biologically characterized loci: *LCT* (rs4988235), which underlies lactase persistence and shows strong positive selection in European populations; *DARC* (rs2814778), the Duffy antigen variant that has reached near-fixation in African populations due to selection pressure from *Plasmodium vivax* malaria; *EDAR* (rs3827760), associated with hair and sweat gland morphology and at high frequency in East Asian populations; and *SLC24A5* (rs1426654), a major determinant of skin pigmentation under a strong selective sweep in Europeans. If the Bayesian differentiation metrics rank these loci highly relative to a neutral background set, it provides direct evidence that the model is capturing genuine selection signals rather than sampling noise.

---

## Statistical Framework

The model rests on a straightforward conjugate update. The observed alternate allele count $x_i$ in population $i$, drawn from $n_i$ total alleles sampled, follows a Binomial likelihood:

$$x_i \sim \text{Binomial}(n_i, p_i)$$

The true allele frequency $p_i$ is treated as unknown and given a Beta prior:

$$p_i \sim \text{Beta}(\alpha, \beta)$$

Because the Beta is the conjugate prior for the Binomial likelihood, the posterior distribution is also Beta-distributed and can be derived analytically:

$$p_i \mid x_i \sim \text{Beta}(\alpha + x_i, \; \beta + n_i - x_i)$$

The practical significance of this conjugacy is substantial: for each SNP-population pair we obtain a complete posterior distribution — including the posterior mean, variance, and highest density credible intervals — in closed form, with no numerical approximation required. The posterior mean $(\alpha + x_i) / (\alpha + \beta + n_i)$ illustrates the shrinkage property directly: estimates from small samples are pulled toward the prior, while estimates from large samples dominate the prior and converge toward the maximum likelihood estimate. This shrinkage is precisely what makes the Bayesian approach more reliable than point estimates when comparing populations of varying sample size.

Differentiation between populations is then quantified at the level of entire posterior distributions rather than point estimates. The Bhattacharyya coefficient measures the overlap between two Beta distributions, with values near zero indicating nearly non-overlapping distributions and values near one indicating near-identical distributions. Complementarily, the Kullback-Leibler divergence measures the information lost when one distribution is used to approximate another. Together, these metrics provide a richer characterization of population differentiation than a single FST value and are directly interpretable in probabilistic terms.

---

## Data

The 1000 Genomes Project Phase 3 (release 20130502) provides whole-genome sequencing data for 2,504 individuals across 26 populations grouped into five continental superpopulations: AFR (African), AMR (Admixed American), EAS (East Asian), EUR (European), and SAS (South Asian). This dataset is ideally suited to the present analysis for several reasons. First, the scale of sequencing — covering the full genome at sufficient depth to call SNPs reliably — provides a rich catalog of common variation with population-level allele frequency information. Second, the breadth of population coverage, spanning five continental groups with meaningful within-continent diversity, enables the kind of multi-population comparison required to identify population-specific selection signals. Third, the data are publicly accessible and well-annotated, with population-level allele frequencies pre-computed and stored directly in the INFO fields of per-chromosome VCF files (as `AFR_AF`, `AMR_AF`, `EAS_AF`, `EUR_AF`, and `SAS_AF`), enabling extraction of allele counts without individual-level genotype processing.

For this project, we will work with a targeted panel of approximately 500–1,000 SNPs: the four validation loci described above, supplemented by a randomly drawn background set to provide a neutral comparison. Because all VCF files are tabix-indexed, specific variants can be queried directly without downloading full per-chromosome files, which range from roughly 200 MB to 1.2 GB. For the four validation loci specifically, allele frequencies can be retrieved programmatically from the Ensembl REST API in R using the `httr` package, requiring no file downloads at all. This combination of local VCF queries and REST API access makes the data pipeline lightweight and reproducible on a standard workstation, with HPC resources reserved only if we extend to genome-wide analysis.

The analysis pipeline will be implemented entirely in R. VCF parsing, where needed, will use the `vcfR` package; allele count data will be manipulated with `dplyr` and `tidyr`; posterior Beta distributions will be computed using base R functions (`dbeta`, `pbeta`, `rbeta`); and posterior distributions along with credible intervals will be visualized using `ggplot2` and `patchwork`. An optional extension using `rstan` would allow us to fit a hierarchical model that partially pools information across populations, but this is not required for the core aims.

---

## Expected Outcomes

The primary deliverable is a ranked list of SNPs ordered by the degree of inter-population differentiation in their posterior allele frequency distributions, alongside visualizations of the posterior Beta distributions for each SNP across the five continental groups. We expect the four validation loci to rank near the top of this list, confirming that the model recovers genuine selection signals. A secondary analysis will compare results under the uninformative and weakly informative priors to characterize how sensitive the rankings are to prior specification, which has practical relevance for applying the framework to future studies with different background frequency assumptions.

---