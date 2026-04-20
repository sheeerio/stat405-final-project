// =============================================================================
// Model 2: Hierarchical log-normal regression
// log(SIZE_HA_i) ~ Normal(mu_i, sigma^2)
// mu_i = beta_0 + beta_month * MONTH_i + beta_year * YEAR_centered_i
//        + alpha[province_i]
//
// KEY DIFFERENCE FROM MODEL 1:
//   alpha[j] ~ Normal(mu_alpha, sigma_alpha)   <- hierarchical prior
//
// mu_alpha and sigma_alpha are estimated from data (hyperparameters).
// This partial-pools province effects: provinces with little data are
// shrunk toward the global mean; data-rich provinces stay near their MLE.
// =============================================================================

data {
  int<lower=1> N;                              // number of fires
  int<lower=1> P;                              // number of provinces
  vector[N] log_size;                          // log(SIZE_HA) — response
  vector[N] month;                             // MONTH centred
  vector[N] year_c;                            // YEAR centred
  array[N] int<lower=1, upper=P> province;// province index (1..P)
}

parameters {
  real beta_month;           // linear month effect
  real beta_year;            // linear year trend

  // Hyperparameters for province effects
  real mu_alpha;             // mean province effect
  real<lower=0> sigma_alpha; // SD of province effects — controls pooling degree

  // Province-level effects (non-centred parameterisation for sampler efficiency)
  vector[P] z_alpha;         // standard normal offsets

  real<lower=0> sigma;       // residual SD
}

transformed parameters {
  // Non-centred parameterisation: alpha = mu_alpha + sigma_alpha * z
  // This avoids the funnel geometry that causes divergences in centred form
  vector[P] alpha = mu_alpha + sigma_alpha * z_alpha;
}

model {
  // ----- hyperpriors -----
  mu_alpha    ~ normal(0, 2);        // prior on mean province effect
  sigma_alpha ~ exponential(1);      // prior on pooling SD; exp(1) puts
                                     // most mass below 2 log-units

  // ----- priors -----
  beta_month ~ normal(0, 1);
  beta_year  ~ normal(0, 0.1);
  z_alpha    ~ std_normal();         // implies alpha ~ Normal(mu_alpha, sigma_alpha)
  sigma      ~ exponential(1);

  // ----- likelihood (vectorised) -----
  vector[N] mu = beta_month * month + beta_year * year_c
                 + alpha[province];
  log_size ~ normal(mu, sigma);
}

generated quantities {
  vector[N] log_lik;   // for LOO-CV
  vector[N] y_rep;     // for PPCs

  {
    vector[N] mu = beta_0 + beta_month * month + beta_year * year_c
                   + alpha[province];
    for (i in 1:N) {
      log_lik[i] = normal_lpdf(log_size[i] | mu[i], sigma);
      y_rep[i]   = normal_rng(mu[i], sigma);
    }
  }
}
