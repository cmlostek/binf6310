# =============================================================================
# BINF 6310 — Final Project
# Bayesian Estimation of Allele Frequencies Across Human Populations:
# Detecting Signatures of Natural Selection Using the Beta-Binomial Model
#
# Author : Cole M.
# Date   : May 2026
# Data   : 1000 Genomes Project Phase 3, via Ensembl REST API
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Packages
# -----------------------------------------------------------------------------
required_pkgs <- c("httr", "jsonlite", "dplyr", "tidyr",
                   "ggplot2", "patchwork", "knitr", "scales", "kableExtra")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[,"Package"])]
if (length(new_pkgs)) install.packages(new_pkgs, repos = "https://cloud.r-project.org")

library(httr); library(jsonlite); library(dplyr); library(tidyr)
library(ggplot2); library(patchwork); library(knitr); library(scales)
library(kableExtra)

# -----------------------------------------------------------------------------
# 1. Constants
# -----------------------------------------------------------------------------

# 1000 Genomes Phase 3 haplotype counts (2 × sampled individuals)
N_HAP <- c(AFR = 1322, AMR = 694, EAS = 1008, EUR = 1006, SAS = 978)

# Superpopulation color palette
POP_COLORS <- c(AFR = "#E41A1C", AMR = "#FF7F00",
                EAS = "#4DAF4A", EUR = "#377EB8", SAS = "#984EA3")

# One representative subpopulation per superpopulation for LD queries
LD_POP_CODES <- c(
  AFR = "1000GENOMES:phase_3:YRI",
  AMR = "1000GENOMES:phase_3:MXL",
  EAS = "1000GENOMES:phase_3:CHB",
  EUR = "1000GENOMES:phase_3:CEU",
  SAS = "1000GENOMES:phase_3:GIH"
)

# Validation SNP panel — 12 loci with documented selection histories
VALIDATION_SNPS <- data.frame(
  rsid         = c("rs4988235","rs2814778","rs3827760","rs1426654",
                   "rs334",    "rs17822931","rs16891982","rs12913832",
                   "rs1805007","rs1229984", "rs1800414", "rs12203592"),
  gene         = c("LCT",  "DARC",  "EDAR",   "SLC24A5",
                   "HBB",  "ABCC11","SLC45A2","HERC2",
                   "MC1R", "ADH1B", "OCA2",   "IRF4"),
  focal_allele = c("A","T","A","A",
                   "T","A","C","A",
                   "T","G","A","C"),
  sel_type     = c("positive","positive","positive","positive",
                   "balancing","positive","positive","positive",
                   "positive","positive","positive","positive"),
  selection_pop = c("EUR","AFR","EAS","EUR/SAS",
                    "AFR","EAS","EUR","EUR",
                    "EUR","EAS","EAS","EUR"),
  ld_query_pop  = c("EUR","AFR","EAS","EUR",
                    "AFR","EAS","EUR","EUR",
                    "EUR","EAS","EAS","EUR"),
  stringsAsFactors = FALSE
)

# Background panel — 10 SNPs without documented strong directional selection
BACKGROUND_RSIDS <- c(
  "rs1801133", "rs4680",   "rs6265",   "rs1800497",
  "rs1801282", "rs7412",   "rs429358", "rs1042522",
  "rs9939609", "rs1695"
)

# -----------------------------------------------------------------------------
# 2. API Helper Functions
# -----------------------------------------------------------------------------

# Fetch allele + phenotype data for a single variant from Ensembl REST API
fetch_variant <- function(rsid, pause_sec = 0.35) {
  url <- paste0("https://rest.ensembl.org/variation/human/", rsid,
                "?pops=1;phenotypes=1")
  Sys.sleep(pause_sec)
  resp <- tryCatch(
    GET(url, add_headers(`Content-Type` = "application/json"), timeout(30)),
    error = function(e) NULL
  )
  if (is.null(resp) || status_code(resp) != 200) return(NULL)
  content(resp, as = "parsed", type = "application/json")
}

# Extract per-superpopulation allele counts from the API response
parse_variant <- function(api_data, rsid, focal_allele = NULL) {
  if (is.null(api_data) || is.null(api_data$populations)) return(NULL)
  target <- if (!is.null(focal_allele)) focal_allele else api_data$minor_allele
  if (is.null(target)) return(NULL)

  pops_df   <- bind_rows(lapply(api_data$populations, as.data.frame))
  superpops <- c("AFR","AMR","EAS","EUR","SAS")

  bind_rows(lapply(superpops, function(sp) {
    entries   <- pops_df[pops_df$population == paste0("1000GENOMES:phase_3:", sp), ]
    if (nrow(entries) == 0) return(NULL)
    n         <- sum(entries$allele_count, na.rm = TRUE)
    focal_row <- entries[entries$allele == target, ]
    alt_count <- if (nrow(focal_row) > 0) focal_row$allele_count[1] else 0L
    data.frame(rsid = rsid, population = sp,
               alt_count = as.integer(alt_count), n = as.integer(n),
               freq_obs  = if (n > 0) alt_count / n else NA_real_,
               stringsAsFactors = FALSE)
  }))
}

# Extract phenotype annotations from the API response
parse_phenotypes <- function(api_data, rsid, gene_label) {
  if (is.null(api_data) || is.null(api_data$phenotypes)) return(NULL)
  phenos <- bind_rows(lapply(api_data$phenotypes, function(p) {
    data.frame(
      rsid        = rsid,
      gene        = gene_label,
      trait       = if (!is.null(p$trait))       p$trait       else NA_character_,
      source      = if (!is.null(p$source))      p$source      else NA_character_,
      risk_allele = if (!is.null(p$risk_allele)) p$risk_allele else NA_character_,
      stringsAsFactors = FALSE
    )
  }))
  phenos[!is.na(phenos$trait) & phenos$trait != "", ]
}

# -----------------------------------------------------------------------------
# 3. Data Retrieval
# -----------------------------------------------------------------------------

val_data_list <- list()
pheno_list    <- list()

cat("Fetching validation SNPs...\n")
for (i in seq_len(nrow(VALIDATION_SNPS))) {
  rs  <- VALIDATION_SNPS$rsid[i]
  fa  <- VALIDATION_SNPS$focal_allele[i]
  cat("  ", rs, "(", VALIDATION_SNPS$gene[i], ")\n")
  api <- fetch_variant(rs)

  parsed <- parse_variant(api, rs, fa)
  if (!is.null(parsed)) {
    parsed$gene  <- VALIDATION_SNPS$gene[i]
    parsed$panel <- "validation"
    val_data_list[[i]] <- parsed
  }
  pheno <- parse_phenotypes(api, rs, VALIDATION_SNPS$gene[i])
  if (!is.null(pheno)) pheno_list[[i]] <- pheno
}
val_data   <- bind_rows(val_data_list)
pheno_data <- bind_rows(pheno_list)

bg_data_list <- list()
cat("Fetching background SNPs...\n")
for (i in seq_along(BACKGROUND_RSIDS)) {
  rs <- BACKGROUND_RSIDS[i]
  cat("  ", rs, "\n")
  api    <- fetch_variant(rs)
  parsed <- parse_variant(api, rs, focal_allele = NULL)
  if (!is.null(parsed)) {
    parsed$gene  <- rs
    parsed$panel <- "background"
    bg_data_list[[i]] <- parsed
  }
}
bg_data  <- bind_rows(bg_data_list)
all_data <- bind_rows(val_data, bg_data)

cat("Retrieved — validation SNPs:", n_distinct(val_data$rsid),
    "| background SNPs:", n_distinct(bg_data$rsid),
    "| phenotype records:", nrow(pheno_data), "\n")

# -----------------------------------------------------------------------------
# 4. Bayesian Beta-Binomial Posterior Computation
# -----------------------------------------------------------------------------
# Model: x_i ~ Binomial(n_i, p_i), p_i ~ Beta(alpha, beta)
# Posterior: p_i | x_i ~ Beta(alpha + x_i, beta + n_i - x_i)

compute_posterior <- function(data, alpha_prior = 1, beta_prior = 1) {
  data %>%
    mutate(
      alpha_post = alpha_prior + alt_count,
      beta_post  = beta_prior  + (n - alt_count),
      post_mean  = alpha_post / (alpha_post + beta_post),
      post_sd    = sqrt((alpha_post * beta_post) /
                        ((alpha_post + beta_post)^2 * (alpha_post + beta_post + 1))),
      ci_lo      = qbeta(0.025, alpha_post, beta_post),
      ci_hi      = qbeta(0.975, alpha_post, beta_post)
    )
}

posterior_data <- compute_posterior(all_data)
val_post       <- posterior_data %>% filter(panel == "validation")

# Wide-format observed frequency table (for inspection)
freq_wide <- val_data %>%
  select(rsid, gene, population, freq_obs) %>%
  pivot_wider(names_from = population, values_from = freq_obs) %>%
  left_join(VALIDATION_SNPS %>% select(rsid, sel_type, selection_pop), by = "rsid") %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
print(freq_wide)

# -----------------------------------------------------------------------------
# 5. Prior Sensitivity Analysis
# -----------------------------------------------------------------------------

compute_posterior_multi <- function(data, priors_list) {
  bind_rows(lapply(seq_along(priors_list), function(i) {
    ap <- priors_list[[i]][1]
    bp <- priors_list[[i]][2]
    compute_posterior(data, ap, bp) %>%
      mutate(prior_label = names(priors_list)[i])
  }))
}

priors_list <- list(
  "Beta(1,1) uniform"      = c(1,   1  ),
  "Beta(2,2) weak"         = c(2,   2  ),
  "Beta(0.5,0.5) Jeffreys" = c(0.5, 0.5)
)

focal_snps_sens <- c("rs4988235", "rs2814778", "rs3827760", "rs334")
focal_names     <- c("LCT",       "DARC",      "EDAR",      "HBB")

sens_data <- all_data %>%
  filter(rsid %in% focal_snps_sens) %>%
  compute_posterior_multi(priors_list)

sens_plot_data <- sens_data %>%
  filter(population %in% c("AFR","EUR","EAS")) %>%
  left_join(data.frame(rsid = focal_snps_sens, gene2 = focal_names), by = "rsid")

fig_sensitivity <- ggplot(
  sens_plot_data,
  aes(x = post_mean, xmin = ci_lo, xmax = ci_hi,
      y = interaction(population, prior_label),
      color = prior_label, shape = prior_label)
) +
  geom_point(size = 2.5) +
  geom_errorbarh(height = 0.25, linewidth = 0.6) +
  facet_wrap(~ gene2, scales = "free_x", ncol = 2) +
  scale_color_manual(values = c("#1f77b4","#ff7f0e","#2ca02c"), name = "Prior") +
  scale_shape_manual(values = c(16, 17, 15), name = "Prior") +
  labs(title = "Prior sensitivity: posterior mean and 95% credible interval",
       x = "Posterior mean allele frequency", y = NULL) +
  theme_classic(base_size = 9) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))

print(fig_sensitivity)

# -----------------------------------------------------------------------------
# 6. Posterior Distribution Visualization
# -----------------------------------------------------------------------------

x_seq <- seq(0.001, 0.999, length.out = 400)

plot_posterior_snp <- function(rsid_sel, gene_name) {
  snp_data <- val_post %>% filter(rsid == rsid_sel)
  dens_df  <- snp_data %>%
    rowwise() %>%
    mutate(dens = list(data.frame(p = x_seq,
                                  density = dbeta(x_seq, alpha_post, beta_post)))) %>%
    unnest(dens)

  ggplot(dens_df, aes(x = p, y = density, color = population, fill = population)) +
    geom_line(linewidth = 0.9) +
    geom_area(alpha = 0.12, position = "identity") +
    scale_color_manual(values = POP_COLORS) +
    scale_fill_manual(values  = POP_COLORS) +
    labs(title = paste0(gene_name, "  (", rsid_sel, ")"),
         x = "Allele frequency p", y = "Posterior density",
         color = "Population", fill = "Population") +
    theme_classic(base_size = 9) +
    theme(legend.position = "right", plot.title = element_text(size = 9, face = "bold"))
}

plots <- mapply(plot_posterior_snp,
                VALIDATION_SNPS$rsid, VALIDATION_SNPS$gene,
                SIMPLIFY = FALSE)

fig_posteriors <- wrap_plots(plots, ncol = 3) +
  plot_annotation(
    title    = "Posterior Beta distributions — all 12 validation SNPs",
    subtitle = "Prior: Beta(1,1). Each curve = one 1000 Genomes superpopulation.",
    theme    = theme(plot.title = element_text(face = "bold", size = 12))
  )

print(fig_posteriors)

# -----------------------------------------------------------------------------
# 7. Population Differentiation Metrics
# -----------------------------------------------------------------------------

# Bhattacharyya coefficient (numerical integration)
bhattacharyya_coef <- function(a1, b1, a2, b2) {
  if (any(is.na(c(a1, b1, a2, b2)))) return(NA_real_)
  tryCatch(
    integrate(function(x) sqrt(dbeta(x, a1, b1) * dbeta(x, a2, b2)),
              0, 1, subdivisions = 200L, rel.tol = 1e-6)$value,
    error = function(e) NA_real_
  )
}

# Closed-form KL divergence between two Beta distributions
kl_div_beta <- function(a1, b1, a2, b2) {
  if (any(is.na(c(a1, b1, a2, b2)))) return(NA_real_)
  tryCatch(
    lbeta(a2, b2) - lbeta(a1, b1) +
      (a1 - a2) * digamma(a1) +
      (b1 - b2) * digamma(b1) +
      (a2 - a1 + b2 - b1) * digamma(a1 + b1),
    error = function(e) NA_real_
  )
}

# Symmetrized KL divergence
kl_sym <- function(a1, b1, a2, b2) {
  (kl_div_beta(a1, b1, a2, b2) + kl_div_beta(a2, b2, a1, b1)) / 2
}

# Compute all pairwise metrics for one SNP's posterior data
pairwise_metrics <- function(snp_post_data) {
  pops  <- snp_post_data$population
  pairs <- combn(length(pops), 2, simplify = FALSE)
  bind_rows(lapply(pairs, function(idx) {
    i  <- idx[1]; j <- idx[2]
    bc <- bhattacharyya_coef(
      snp_post_data$alpha_post[i], snp_post_data$beta_post[i],
      snp_post_data$alpha_post[j], snp_post_data$beta_post[j]
    )
    kl <- kl_sym(
      snp_post_data$alpha_post[i], snp_post_data$beta_post[i],
      snp_post_data$alpha_post[j], snp_post_data$beta_post[j]
    )
    data.frame(
      pop1 = pops[i], pop2 = pops[j],
      BD   = if (!is.na(bc) && bc > 0) -log(bc) else NA_real_,
      KL   = kl
    )
  }))
}

cat("Computing pairwise differentiation metrics...\n")
snp_scores <- posterior_data %>%
  group_by(rsid, gene, panel) %>%
  group_split() %>%
  lapply(function(df) {
    pw       <- pairwise_metrics(df)
    pw$rsid  <- df$rsid[1]
    pw$gene  <- df$gene[1]
    pw$panel <- df$panel[1]
    pw
  }) %>%
  bind_rows() %>%
  group_by(rsid, gene, panel) %>%
  summarise(
    max_BD  = max(BD,  na.rm = TRUE),
    mean_BD = mean(BD, na.rm = TRUE),
    max_KL  = max(KL,  na.rm = TRUE),
    mean_KL = mean(KL, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(max_BD))

print(snp_scores)

# -----------------------------------------------------------------------------
# 8. Lollipop Ranking Plot
# -----------------------------------------------------------------------------

rank_df <- snp_scores %>%
  mutate(label = ifelse(panel == "validation", gene, rsid)) %>%
  arrange(max_BD) %>%
  mutate(label = factor(label, levels = label))

fig_rank <- ggplot(rank_df, aes(x = max_BD, y = label, color = panel, shape = panel)) +
  geom_segment(aes(x = 0, xend = max_BD, y = label, yend = label),
               color = "grey75", linewidth = 0.5) +
  geom_point(size = 3.5) +
  scale_color_manual(values = c(validation = "#D7191C", background = "#2C7BB6"),
                     labels = c("Validation (selection)", "Background (neutral)"),
                     name = NULL) +
  scale_shape_manual(values = c(validation = 17, background = 16), name = NULL) +
  labs(title = "Population differentiation: max pairwise Bhattacharyya distance",
       subtitle = "Validation SNPs (red) vs. background (blue)",
       x = "Max pairwise Bhattacharyya distance (BD)", y = NULL) +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_text(size = 9), legend.position = "bottom")

print(fig_rank)

# -----------------------------------------------------------------------------
# 9. Statistical Testing: Wilcoxon Rank-Sum Test
# -----------------------------------------------------------------------------

val_scores <- snp_scores %>% filter(panel == "validation") %>% pull(max_BD)
bg_scores  <- snp_scores %>% filter(panel == "background")  %>% pull(max_BD)

# Background distribution parameters for Z-score computation
bg_mean <- mean(bg_scores, na.rm = TRUE)
bg_sd   <- sd(bg_scores,   na.rm = TRUE)

# Z-score standardization relative to background
snp_scores_z <- snp_scores %>%
  mutate(
    z_score    = (max_BD - bg_mean) / bg_sd,
    is_outlier = z_score > 2,
    label      = ifelse(panel == "validation", gene, rsid)
  )

# One-sided Wilcoxon rank-sum test
wt <- wilcox.test(val_scores, bg_scores, alternative = "greater")
cat(sprintf("\nWilcoxon rank-sum (validation > background):\n  W = %.0f\n  p = %.6f\n",
            wt$statistic, wt$p.value))

# Boxplot
fig_boxplot <- ggplot(snp_scores, aes(x = panel, y = max_BD, color = panel)) +
  geom_boxplot(outlier.shape = NA, width = 0.35, linewidth = 0.8, fill = NA) +
  geom_jitter(aes(shape = panel), width = 0.12, size = 3.5, alpha = 0.85) +
  geom_text(
    data    = snp_scores %>% filter(panel == "validation"),
    aes(label = gene), hjust = -0.25, size = 3, color = "black"
  ) +
  scale_color_manual(values = c(validation = "#D7191C", background = "#2C7BB6")) +
  scale_shape_manual(values = c(validation = 17, background = 16)) +
  scale_x_discrete(labels = c(
    validation = "Selection loci\n(n = 12)",
    background = "Background\n(n = 10)"
  )) +
  labs(
    title    = "Max pairwise BD: selection loci vs. neutral background",
    subtitle = sprintf("One-sided Wilcoxon: W = %.0f, p = %.4f", wt$statistic, wt$p.value),
    x = NULL, y = "Max pairwise Bhattacharyya distance"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")

print(fig_boxplot)

# Z-score table
z_table <- snp_scores_z %>%
  filter(panel == "validation") %>%
  select(label, max_BD, mean_BD, z_score, is_outlier) %>%
  arrange(desc(z_score)) %>%
  mutate(across(c(max_BD, mean_BD, z_score), ~ round(.x, 3)))

cat("\nZ-score table (validation SNPs vs. background):\n")
print(z_table)

# KL vs BD comparison plot
combined_metrics <- snp_scores_z %>%
  filter(panel == "validation") %>%
  select(label, max_BD, max_KL, z_score) %>%
  arrange(desc(max_BD)) %>%
  mutate(label = factor(label, levels = rev(label)))

p_bd <- ggplot(combined_metrics, aes(x = max_BD, y = label)) +
  geom_col(fill = "#D7191C", alpha = 0.8) +
  labs(x = "Max BD", y = NULL, title = "Bhattacharyya distance") +
  theme_classic(base_size = 9)

p_kl <- ggplot(combined_metrics, aes(x = max_KL, y = label)) +
  geom_col(fill = "#2C7BB6", alpha = 0.8) +
  labs(x = "Max KL divergence", y = NULL, title = "KL divergence") +
  theme_classic(base_size = 9) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

fig_metrics <- (p_bd | p_kl) +
  plot_annotation(
    title = "Differentiation metrics: BD vs. KL divergence (validation SNPs)",
    theme = theme(plot.title = element_text(face = "bold", size = 11))
  )

print(fig_metrics)

# -----------------------------------------------------------------------------
# 10. Haplotype-Based LD Analysis
# -----------------------------------------------------------------------------

fetch_ld <- function(rsid, population, window_kb = 300, r2_min = 0.1) {
  pop_full <- LD_POP_CODES[population]
  url <- paste0("https://rest.ensembl.org/ld/human/", rsid, "/",
                URLencode(pop_full, reserved = TRUE),
                "?window_size=", window_kb, "&r2=", r2_min)
  Sys.sleep(0.5)
  resp <- tryCatch(
    GET(url, add_headers(`Content-Type` = "application/json"), timeout(60)),
    error = function(e) NULL
  )
  if (is.null(resp) || status_code(resp) != 200) return(NULL)
  result <- content(resp, as = "parsed", type = "application/json")
  if (length(result) == 0) return(NULL)
  bind_rows(lapply(result, function(x) {
    data.frame(
      focal_snp   = rsid,
      partner_snp = if (!is.null(x$variation2)) x$variation2 else NA_character_,
      r2          = as.numeric(x$r2),
      d_prime     = as.numeric(x$d_prime),
      population  = population,
      stringsAsFactors = FALSE
    )
  }))
}

cat("Fetching LD data...\n")
ld_results <- list()

for (i in seq_len(nrow(VALIDATION_SNPS))) {
  rs  <- VALIDATION_SNPS$rsid[i]
  pop <- VALIDATION_SNPS$ld_query_pop[i]
  cat("  LD:", rs, "in", pop, "\n")
  ld <- fetch_ld(rs, pop, window_kb = 300, r2_min = 0.1)
  if (!is.null(ld)) {
    ld$gene     <- VALIDATION_SNPS$gene[i]
    ld$sel_type <- VALIDATION_SNPS$sel_type[i]
    ld_results[[length(ld_results) + 1]] <- ld
  }
}

for (rs in BACKGROUND_RSIDS) {
  cat("  LD (bg):", rs, "\n")
  ld <- fetch_ld(rs, "EUR", window_kb = 300, r2_min = 0.1)
  if (!is.null(ld)) {
    ld$gene     <- rs
    ld$sel_type <- "background"
    ld_results[[length(ld_results) + 1]] <- ld
  }
}

ld_all <- bind_rows(ld_results)
cat("Total LD variant pairs retrieved:", nrow(ld_all), "\n")

# Haplotype extent scores
hap_scores <- ld_all %>%
  group_by(gene, focal_snp, sel_type) %>%
  summarise(
    n_partners  = n(),
    frac_r2_hi  = mean(r2 >= 0.8, na.rm = TRUE),
    median_r2   = median(r2, na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  arrange(desc(frac_r2_hi))

cat("\nHaplotype extent scores:\n")
print(hap_scores)

# LD density plot
if (nrow(ld_all) > 0) {
  fig_ld_density <- ld_all %>%
    ggplot(aes(x = r2, fill = sel_type, color = sel_type)) +
    geom_density(alpha = 0.3, adjust = 1.3) +
    facet_wrap(~ gene, ncol = 4) +
    scale_fill_manual(
      values = c(background = "steelblue", positive = "#D7191C", balancing = "#FF7F00"),
      name   = "Selection type"
    ) +
    scale_color_manual(
      values = c(background = "steelblue", positive = "#D7191C", balancing = "#FF7F00"),
      name   = "Selection type"
    ) +
    labs(
      title    = "r² distribution within 300 kb of each focal SNP",
      subtitle = "Right-skewed r² = extended haplotype = selective sweep",
      x = "r² with focal SNP", y = "Density"
    ) +
    theme_classic(base_size = 8) +
    theme(strip.text = element_text(face = "bold", size = 7), legend.position = "bottom")

  print(fig_ld_density)
}

# Haplotype extent score barplot
fig_hes <- hap_scores %>%
  arrange(frac_r2_hi) %>%
  mutate(gene = factor(gene, levels = gene)) %>%
  ggplot(aes(x = frac_r2_hi, y = gene, fill = sel_type)) +
  geom_col(alpha = 0.85) +
  scale_fill_manual(
    values = c(background = "steelblue", positive = "#D7191C", balancing = "#FF7F00"),
    name   = "Selection type"
  ) +
  labs(
    title = "Haplotype extent scores — fraction of 300 kb partners with r² ≥ 0.8",
    x = "Fraction r² ≥ 0.8", y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(legend.position = "right")

print(fig_hes)

# Concordance: BD vs. haplotype extent
concordance_df <- snp_scores %>%
  filter(panel == "validation") %>%
  select(gene, max_BD) %>%
  left_join(hap_scores %>% select(gene, frac_r2_hi, sel_type), by = "gene")

fig_concordance <- ggplot(
  concordance_df,
  aes(x = max_BD, y = frac_r2_hi, color = sel_type, label = gene)
) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text(hjust = -0.15, size = 3, color = "black") +
  scale_color_manual(
    values = c(positive = "#D7191C", balancing = "#FF7F00"),
    name   = "Selection type"
  ) +
  labs(
    title    = "BD vs. haplotype extent: concordance of two independent metrics",
    subtitle = "Higher BD + higher HES = stronger selection signal",
    x = "Max pairwise Bhattacharyya distance", y = "Haplotype extent score (frac r² ≥ 0.8)"
  ) +
  theme_classic(base_size = 11)

print(fig_concordance)

# -----------------------------------------------------------------------------
# 11. Session Info
# -----------------------------------------------------------------------------
sessionInfo()
