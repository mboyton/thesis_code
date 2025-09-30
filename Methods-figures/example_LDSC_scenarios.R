# Load libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)  # For combining plots

set.seed(123)

# Simulate LD scores for 500 SNPs
n_snps <- 500
ld_scores <- runif(n_snps, 0, 300)

# -------------------------------
# Univariate LDSC Simulations
# -------------------------------

# Scenario 1: Polygenicity - test statistics increase with LD
chi2_polygenic <- 1 + 0.01 * ld_scores + rnorm(n_snps, mean = 0, sd = 0.5)

# Scenario 2: Confounding (e.g. drift) - uniform inflation
chi2_drift <- 1.3 + rnorm(n_snps, mean = 0, sd = 0.5)

# Bin LD scores
ld_bins <- cut(ld_scores, breaks = seq(0, 300, by = 30), include.lowest = TRUE)

df_uni <- data.frame(ld_score = ld_scores,
                     ld_bin = ld_bins,
                     chi2_polygenic = chi2_polygenic,
                     chi2_drift = chi2_drift)

# Compute mean chi2 per LD bin
mean_chi2 <- df_uni %>%
 group_by(ld_bin) %>%
 summarise(
  mean_ld_score = mean(ld_score),
  mean_chi2_polygenic = mean(chi2_polygenic),
  mean_chi2_drift = mean(chi2_drift)
 )

# Reshape for plotting
mean_chi2_long <- mean_chi2 %>%
 pivot_longer(cols = c(mean_chi2_polygenic, mean_chi2_drift),
              names_to = "scenario", values_to = "mean_chi2") %>%
 mutate(scenario = recode(scenario,
                          mean_chi2_polygenic = "Polygenic",
                          mean_chi2_drift = "Drift/Confounding"))

# -------------------------------
# Cross-trait LDSC Simulations
# -------------------------------

# Scenario 1: Genetic correlation
z1_corr <- rnorm(n_snps, mean = 0.01 * ld_scores, sd = 1)
z2_corr <- z1_corr + rnorm(n_snps, mean = 0, sd = 1)
zprod_corr <- z1_corr * z2_corr

# Scenario 2: Shared confounding only
z1_conf <- rnorm(n_snps, mean = 0.3, sd = 1)
z2_conf <- rnorm(n_snps, mean = 0.3, sd = 1)
zprod_conf <- z1_conf * z2_conf

df_biv <- data.frame(ld_score = ld_scores,
                     ld_bin = ld_bins,
                     zprod_corr = zprod_corr,
                     zprod_conf = zprod_conf)

# Compute mean zprod per LD bin
mean_zprod <- df_biv %>%
 group_by(ld_bin) %>%
 summarise(
  mean_ld_score = mean(ld_score),
  mean_zprod_corr = mean(zprod_corr),
  mean_zprod_conf = mean(zprod_conf)
 )

# Reshape for plotting
mean_zprod_long <- mean_zprod %>%
 pivot_longer(cols = c(mean_zprod_corr, mean_zprod_conf),
              names_to = "scenario", values_to = "mean_zprod") %>%
 mutate(scenario = recode(scenario,
                          mean_zprod_corr = "Genetic Correlation",
                          mean_zprod_conf = "Confounding Only"))

# -------------------------------
# Plot A: Univariate LDSC
# -------------------------------

p1b <- ggplot(mean_chi2_long, aes(x = mean_ld_score, y = mean_chi2, color = scenario)) +
 geom_point(size = 3) +
 geom_smooth(method = "lm", se = FALSE) +
 scale_color_manual(values = c("Polygenic" = "#5DADE2",
                               "Drift/Confounding" = "#EC7063")) +
 labs(
  x = "Mean LD Score (per bin)",
  y = expression("Mean " ~ chi^2),
  color = "Scenario",
  title = expression("Figure 1B. Univariate LDSC: LD Score vs " ~ chi^2),
  subtitle = "Blue = Polygenic (slope > 0), Red = Confounding (flat)"
 ) +
 theme_minimal() +
 theme(legend.title = element_text(hjust = 0.5))

# -------------------------------
# Plot B : Cross-trait LDSC
# -------------------------------

p1c <- ggplot(mean_zprod_long, aes(x = mean_ld_score, y = mean_zprod, color = scenario)) +
 geom_point(size = 3) +
 geom_smooth(method = "lm", se = FALSE) +
 scale_color_manual(values = c("Genetic Correlation" = "#5DADE2",
                               "Confounding Only" = "#EC7063")) +
 labs(
  x = "Mean LD Score (per bin)",
  y = bquote("Mean " ~ Z[1] %.% Z[2]),
  color = "Scenario",
  title = bquote("Figure 1C. Cross-trait LDSC: LD Score vs " ~ Z[1] %.% Z[2]),
  subtitle = "Blue = Genetic Correlation (slope > 0), Red = Confounding (flat)"
 ) +
 theme_minimal() +
 theme(legend.title = element_text(hjust = 0.5))

# -------------------------------
# Display side-by-side
# -------------------------------
p1b + p1c
