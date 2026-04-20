library(cmdstanr)   # install_cmdstan() once if not already done
library(dplyr)
library(ggplot2)
library(bayesplot)  # install.packages("bayesplot")
library(loo)        # install.packages("loo")
library(sf)

dir.create("figures", showWarnings = FALSE)
dir.create("output",  showWarnings = FALSE)   # for fitted model objects

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
  slice_sample(n = 10000) %>%
  mutate(province_fac = droplevels(province_fac))   # drop empty levels

province_levels <- levels(df$province_fac)
cat("Provinces in model:", paste(province_levels, collapse = ", "), "\n")
cat("N =", nrow(df), " P =", length(province_levels), "\n")

# ----- Stan data list -----
stan_data <- list(
  N        = nrow(df),
  P        = length(province_levels),
  log_size = df$log_size,
  month    = df$month_c,
  year_c   = df$year_c,
  province = as.integer(df$province_fac)
)

mod1 <- cmdstan_model("src/model1_lognormal.stan")

fit1 <- mod1$sample(
  data            = stan_data,
  seed            = 405,
  chains          = 4,
  parallel_chains = 2,         # use all 4 cores; reduce if needed
  iter_warmup     = 1000,
  iter_sampling   = 2000,
  refresh         = 10,
  adapt_delta     = 0.9        # slightly above default; helps with province effects
)

# Save fitted object so you don't have to re-run
fit1$save_object("output/fit1.rds")
fit1 <- readRDS("output/fit1.rds")

cat("\n===== MCMC Diagnostics =====\n")
fit1$diagnostic_summary()   # prints divergences, max treedepth, E-BFMI

# R-hat and ESS for scalar parameters
params_scalar <- c("beta_month", "beta_year", "sigma")
summary_scalar <- fit1$summary(params_scalar)
print(summary_scalar)

# Rule of thumb: Rhat < 1.01, ESS_bulk > 400, ESS_tail > 400
if (any(summary_scalar$rhat > 1.01, na.rm = TRUE)) {
  warning("Some Rhat > 1.01 — chains may not have converged. Consider longer warmup.")
} else {
  cat("✓ All Rhat < 1.01\n")
}

# 3a. Trace plots for scalar parameters
draws_arr <- fit1$draws(format = "array")

p_trace <- mcmc_trace(draws_arr, pars = params_scalar) +
  labs(title = "Model 1 — Trace Plots (scalar parameters)") +
  theme_minimal()
ggsave("figures/m1_trace.png", p_trace, width = 10, height = 7, dpi = 150)

# 3b. Pairs plot (check for funnel / correlation between key params)
p_pairs <- mcmc_pairs(draws_arr,
                      pars = c("beta_0", "beta_month", "beta_year", "sigma"),
                      off_diag_args = list(size = 0.3, alpha = 0.3))
ggsave("figures/m1_pairs.png", p_pairs, width = 8, height = 8, dpi = 150)

# =============================================================================
# 4. POSTERIOR SUMMARIES
# =============================================================================

# 4a. Scalar coefficients with 90% credible intervals
draws_df <- fit1$draws(format = "df")

coef_summary <- fit1$summary(
  variables = params_scalar,
  mean, sd,
  ~quantile(.x, probs = c(0.05, 0.95))
)
cat("\n===== Posterior Summary (scalar parameters) =====\n")
print(coef_summary)

# Interpret year trend on original scale
beta_year_draws <- draws_df$beta_year
cat("\nbeta_year posterior mean:", round(mean(beta_year_draws), 5))
cat("\n  → exp(beta_year) =", round(mean(exp(beta_year_draws)), 5),
    " (multiplicative change in fire size per year)\n")
cat("  → Implied % change per decade:",
    round((mean(exp(10 * beta_year_draws)) - 1) * 100, 1), "%\n")

# 4b. Province effects
alpha_summary <- fit1$summary(
  paste0("alpha[", seq_along(province_levels), "]"),
  mean, sd, ~quantile(.x, probs = c(0.05, 0.95))
) %>%
  mutate(province = province_levels)

p_alpha <- ggplot(alpha_summary,
                  aes(x = reorder(province, mean), y = mean,
                      ymin = `5%`, ymax = `95%`)) +
  geom_pointrange(colour = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  coord_flip() +
  labs(
    title    = "Model 1 — Province Fixed Effects (alpha)",
    subtitle = "Posterior mean ± 90% CI on log scale\n(Model 2 will shrink these toward a common mean)",
    x        = NULL,
    y        = "Effect on log(SIZE_HA)"
  ) +
  theme_minimal(base_size = 13)

ggsave("figures/m1_province_effects.png", p_alpha, width = 8, height = 6, dpi = 150)
print(p_alpha)

# =============================================================================
# 5. POSTERIOR PREDICTIVE CHECKS
# =============================================================================

y_rep_mat <- fit1$draws("y_rep", format = "matrix")[1:200, ]
y_obs <- df$log_size          # <-- single source of truth

set.seed(42)

p_ppc_dens <- ppc_dens_overlay(
  y    = y_obs,
  yrep = y_rep_mat[idx, ]
) +
  labs(
    title    = "Model 1 — Posterior Predictive Check (density overlay)",
    subtitle = "Dark line = observed; light lines = 200 simulated datasets",
    x        = "log(SIZE_HA)"
  ) +
  theme_minimal(base_size = 13)

ggsave("figures/m1_ppc_density.png", p_ppc_dens, width = 8, height = 5, dpi = 150)
print(p_ppc_dens)

# 5b. Stat checks: mean and SD of log-size
p_ppc_mean <- ppc_stat(df$log_size, y_rep_mat, stat = "mean") +
  labs(title = "PPC — Mean of log(SIZE_HA)") + theme_minimal()

p_ppc_sd   <- ppc_stat(df$log_size, y_rep_mat, stat = "sd") +
  labs(title = "PPC — SD of log(SIZE_HA)") + theme_minimal()

ggsave("figures/m1_ppc_mean.png", p_ppc_mean, width = 6, height = 4, dpi = 150)
ggsave("figures/m1_ppc_sd.png",   p_ppc_sd,   width = 6, height = 4, dpi = 150)

# =============================================================================
# 6. LOO-CV (save for Model 1 vs Model 2 comparison in Week 4)
# =============================================================================

# Replace the entire Section 6 with this:
loo1 <- fit1$loo(cores = 1)
saveRDS(loo1, "output/loo1.rds")
#log_lik_1 <- fit1$draws("log_lik", format = "matrix")
#loo1       <- loo(log_lik_1, cores = 1)
#saveRDS(loo1, "output/loo1.rds")

#cat("\n===== LOO-CV (Model 1) =====\n")
print(loo1)

# Check for highly influential observations (Pareto k > 0.7)
bad_k <- sum(loo1$diagnostics$pareto_k > 0.7)
cat("Observations with Pareto k > 0.7:", bad_k, "\n")
if (bad_k > 0) {
  cat("  → Consider moment matching or removing extreme outliers\n")
}

# =============================================================================
# 7. MARGINAL POSTERIOR PLOTS (for write-up)
# =============================================================================

p_post <- mcmc_areas(
  draws_arr,
  pars = params_scalar,
  prob = 0.9
) +
  labs(title = "Model 1 — Marginal Posteriors (90% CI shaded)") +
  theme_minimal(base_size = 13)

ggsave("figures/m1_marginal_posteriors.png", p_post, width = 8, height = 6, dpi = 150)
print(p_post)

cat("\n✓ Week 2 complete. Figures saved to figures/, model to output/fit1.rds\n")
cat("  Next: fit Model 2 (hierarchical) in fit_model2.R\n")

