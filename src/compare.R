# =============================================================================
# STAT 405 — Week 4: Model Comparison, Trend Posterior, Final Figures
# =============================================================================

library(cmdstanr)
library(dplyr)
library(ggplot2)
library(bayesplot)
library(loo)
library(tidyr)
library(sf)

dir.create("figures", showWarnings = FALSE)
dir.create("output",  showWarnings = FALSE)

# =============================================================================
# 1. LOAD DATA (same prep as always — must match what models were fit on)
# =============================================================================

nfdb <- st_read("data/NFDB_point/NFDB_point_20250519.shp")

df <- nfdb %>%
  st_drop_geometry() %>%
  filter(!is.na(SIZE_HA), SIZE_HA > 0,
         !is.na(MONTH), !is.na(YEAR),
         YEAR >= 1959, YEAR <= 2024,
         !is.na(SRC_AGENCY)) %>%
  mutate(
    log_size     = log(SIZE_HA),
    month_c      = as.numeric(MONTH) - mean(as.numeric(MONTH), na.rm = TRUE),
    year_c       = YEAR - mean(YEAR, na.rm = TRUE),
    province_fac = factor(SRC_AGENCY)
  ) %>%
  filter(is.finite(log_size))

set.seed(405)
df <- df %>%
  slice_sample(n = 10000) %>%
  mutate(province_fac = droplevels(province_fac))

province_levels <- levels(df$province_fac)
y_obs           <- df$log_size
year_mean       <- mean(nfdb %>% st_drop_geometry() %>%
                          filter(YEAR >= 1959, YEAR <= 2024) %>%
                          pull(YEAR), na.rm = TRUE)

cat("N =", nrow(df), " P =", length(province_levels), "\n")
cat("year_mean used for centering:", round(year_mean, 1), "\n")

# =============================================================================
# 2. LOAD FITTED MODELS AND LOO OBJECTS
# =============================================================================

fit1      <- readRDS("output/fit1.rds")
fit2_nuts <- readRDS("output/fit2_nuts.rds")
loo1      <- readRDS("output/loo1.rds")
loo2      <- readRDS("output/loo2.rds")

# =============================================================================
# 3. LOO-CV MODEL COMPARISON
# =============================================================================

cat("\n===== LOO-CV: Model 1 vs. Model 2 =====\n")
loo_comp <- loo_compare(loo1, loo2)
print(loo_comp)

# loo_compare ranks models best-first.
# elpd_diff > 0 for Model 2 means it has better predictive accuracy.
# se_diff: if |elpd_diff| > 2 * se_diff, the difference is substantial.

cat("\nInterpretation:\n")
best <- rownames(loo_comp)[1]
cat(" Best model by LOO-CV:", best, "\n")
cat(" elpd_diff:", round(loo_comp[2, "elpd_diff"], 1),
    " (SE:", round(loo_comp[2, "se_diff"], 1), ")\n")
if (abs(loo_comp[2, "elpd_diff"]) > 2 * loo_comp[2, "se_diff"]) {
  cat(" → Difference is statistically meaningful (|diff| > 2 SE)\n")
} else {
  cat(" → Difference is within 2 SE — models are roughly equivalent predictively\n")
}

# Visualise LOO comparison as a simple bar chart
loo_df <- data.frame(
  model  = c("Model 1\n(fixed province effects)",
             "Model 2\n(hierarchical)"),
  elpd   = c(loo1$estimates["elpd_loo", "Estimate"],
             loo2$estimates["elpd_loo", "Estimate"]),
  se     = c(loo1$estimates["elpd_loo", "SE"],
             loo2$estimates["elpd_loo", "SE"])
)

p_loo <- ggplot(loo_df, aes(x = model, y = elpd, ymin = elpd - se, ymax = elpd + se)) +
  geom_col(fill = "steelblue", alpha = 0.8, width = 0.5) +
  geom_errorbar(width = 0.15, colour = "grey30") +
  labs(
    title    = "LOO-CV Comparison: Model 1 vs. Model 2",
    subtitle = "Higher ELPD = better out-of-sample predictive accuracy (±1 SE shown)",
    x        = NULL,
    y        = "Expected Log Predictive Density (ELPD)"
  ) +
  theme_minimal(base_size = 13)

ggsave("figures/loo_comparison.png", p_loo, width = 7, height = 5, dpi = 150)
print(p_loo)

# =============================================================================
# 4. POSTERIOR PREDICTIVE CHECK COMPARISON (Model 1 vs Model 2 side by side)
# =============================================================================

y_rep1 <- fit1$draws("y_rep",      format = "matrix")[1:200, ]
y_rep2 <- fit2_nuts$draws("y_rep", format = "matrix")[1:200, ]

p_ppc1 <- ppc_dens_overlay(y = y_obs, yrep = y_rep1) +
  labs(title = "Model 1 PPC", x = "log(SIZE_HA)") +
  theme_minimal(base_size = 12) +
  xlim(-15, 15)

p_ppc2 <- ppc_dens_overlay(y = y_obs, yrep = y_rep2) +
  labs(title = "Model 2 PPC", x = "log(SIZE_HA)") +
  theme_minimal(base_size = 12) +
  xlim(-15, 15)

# Save side by side using patchwork
library(patchwork)
p_ppc_both <- p_ppc1 + p_ppc2 +
  plot_annotation(
    title    = "Posterior Predictive Checks: Model 1 vs. Model 2",
    subtitle = "Both models miss the bimodal peaks — log-normal likelihood limitation"
  )
ggsave("figures/ppc_comparison.png", p_ppc_both, width = 13, height = 5, dpi = 150)
print(p_ppc_both)

# =============================================================================
# 5. TREND POSTERIOR PLOT
# "What does the model say about expected fire size across time?"
# This is the key substantive figure for the write-up.
# =============================================================================

# Use Model 2 (the preferred model) draws
draws_df2 <- fit2_nuts$draws(format = "df")

# For trend plot, we evaluate mu at: average month, average province (mu_alpha)
# mu = beta_month * 0 + beta_year * year_c + mu_alpha
# => mu = mu_alpha + beta_year * year_c

years_plot <- seq(1960, 2024, by = 1)
year_c_plot <- years_plot - year_mean

# For each posterior draw, compute the expected log-size at each year
# then exponentiate to get median fire size in hectares
set.seed(42)
n_draws_trend <- 500
draw_idx <- sample(nrow(draws_df2), n_draws_trend)

trend_list <- lapply(draw_idx, function(i) {
  mu_alpha   <- draws_df2$mu_alpha[i]
  beta_year  <- draws_df2$beta_year[i]
  mu_t       <- mu_alpha + beta_year * year_c_plot
  data.frame(
    year        = years_plot,
    log_size_mu = mu_t,
    ha_median   = exp(mu_t)   # median of lognormal = exp(mu)
  )
})

trend_df <- bind_rows(trend_list, .id = "draw")

# Summarise across draws
trend_summary <- trend_df %>%
  group_by(year) %>%
  summarise(
    median_ha = median(ha_median),
    lo90      = quantile(ha_median, 0.05),
    hi90      = quantile(ha_median, 0.95),
    lo50      = quantile(ha_median, 0.25),
    hi50      = quantile(ha_median, 0.75),
    .groups   = "drop"
  )

p_trend <- ggplot(trend_summary, aes(x = year)) +
  geom_ribbon(aes(ymin = lo90, ymax = hi90), fill = "steelblue", alpha = 0.2) +
  geom_ribbon(aes(ymin = lo50, ymax = hi50), fill = "steelblue", alpha = 0.35) +
  geom_line(aes(y = median_ha), colour = "steelblue", linewidth = 1) +
  scale_y_log10(labels = scales::comma) +
  labs(
    title    = "Posterior Trend: Expected Median Fire Size Over Time",
    subtitle = "Model 2 (hierarchical), averaged over provinces and mid-season\n50% and 90% credible bands shown",
    x        = "Year",
    y        = "Expected Median Fire Size (ha, log scale)"
  ) +
  theme_minimal(base_size = 13)

ggsave("figures/trend_posterior.png", p_trend, width = 9, height = 5, dpi = 150)
print(p_trend)

# Print key numbers for write-up
pred_1960 <- trend_summary %>% filter(year == 1960)
pred_2024 <- trend_summary %>% filter(year == 2024)

cat("\n--- Trend posterior summary ---\n")
cat("Expected median fire size in 1960:",
    round(pred_1960$median_ha, 1), "ha",
    "[", round(pred_1960$lo90, 1), ",", round(pred_1960$hi90, 1), "] 90% CI\n")
cat("Expected median fire size in 2024:",
    round(pred_2024$median_ha, 1), "ha",
    "[", round(pred_2024$lo90, 1), ",", round(pred_2024$hi90, 1), "] 90% CI\n")
cat("Implied multiplicative change 1960→2024:",
    round(pred_2024$median_ha / pred_1960$median_ha, 2), "x\n")

# =============================================================================
# 6. BETA_YEAR POSTERIOR — explicit display of the trend coefficient
# =============================================================================

beta_year_draws2 <- draws_df2$beta_year

cat("\n--- beta_year (Model 2) ---\n")
cat("Posterior mean:", round(mean(beta_year_draws2), 5), "\n")
cat("90% CI: [", round(quantile(beta_year_draws2, 0.05), 5), ",",
    round(quantile(beta_year_draws2, 0.95), 5), "]\n")
cat("P(beta_year > 0):", round(mean(beta_year_draws2 > 0), 3), "\n")
cat("Implied % change per decade:",
    round((mean(exp(10 * beta_year_draws2)) - 1) * 100, 1), "%\n")

p_beta_year <- ggplot(data.frame(beta_year = beta_year_draws2), aes(x = beta_year)) +
  geom_density(fill = "steelblue", alpha = 0.5, colour = "steelblue") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = mean(beta_year_draws2), colour = "steelblue", linewidth = 1) +
  labs(
    title    = "Posterior Distribution of beta_year (Model 2)",
    subtitle = "Dashed = 0 (no trend); solid = posterior mean",
    x        = "beta_year (log-scale change in fire size per year)",
    y        = "Density"
  ) +
  theme_minimal(base_size = 13)

ggsave("figures/beta_year_posterior.png", p_beta_year, width = 7, height = 4, dpi = 150)
print(p_beta_year)

# =============================================================================
# 7. PARETO K DIAGNOSTIC PLOT (influential observations)
# =============================================================================

pk_df <- data.frame(
  obs    = seq_along(loo2$diagnostics$pareto_k),
  k      = loo2$diagnostics$pareto_k,
  log_size = y_obs
)

p_pk <- ggplot(pk_df, aes(x = obs, y = k, colour = k > 0.7)) +
  geom_point(size = 0.8, alpha = 0.6) +
  geom_hline(yintercept = 0.7, linetype = "dashed", colour = "firebrick") +
  scale_colour_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick"),
                      labels = c("k ≤ 0.7 (ok)", "k > 0.7 (influential)")) +
  labs(
    title    = "LOO Pareto k Diagnostics (Model 2)",
    subtitle = "Points above dashed line are influential observations",
    x        = "Observation index",
    y        = "Pareto k",
    colour   = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

ggsave("figures/pareto_k.png", p_pk, width = 8, height = 4, dpi = 150)
print(p_pk)

cat("\n✓ Week 4 complete. All figures saved to figures/\n")
cat("  Key outputs for write-up:\n")
cat("  - figures/loo_comparison.png\n")
cat("  - figures/trend_posterior.png\n")
cat("  - figures/beta_year_posterior.png\n")
cat("  - figures/ppc_comparison.png\n")
cat("  - figures/pareto_k.png\n")