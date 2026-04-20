// =============================================================================
// Model 1: Simple log-normal regression (baseline)
// log(SIZE_HA_i) ~ Normal(mu_i, sigma^2)
// mu_i = beta_0 + beta_month * MONTH_i + beta_year * YEAR_centered_i
//        + alpha[province_i]
//
// Province effects are FIXED (unpooled) here — contrast with Model 2
// which will give them a hierarchical prior.
// =============================================================================

data {
  int<lower=1> N;                        // number of fires
  int<lower=1> P;                        // number of provinces
  vector[N] log_size;                    // log(SIZE_HA) — response
  vector[N] month;                       // MONTH (1–12, centered below)
  vector[N] year_c;                      // YEAR centered on mean year
  array[N] int<lower=1, upper=P> province; // province index (1..P)
}

parameters {
  real beta_month;           // linear month effect
  real beta_year;            // linear year trend
  vector[P] alpha;           // province fixed effects (one per province)
  real<lower=0> sigma;       // residual SD on log scale
}

model {
  beta_month ~ normal(0, 1);
  beta_year  ~ normal(0, 0.1);
  alpha      ~ normal(0, 2);
  sigma      ~ exponential(1);

  // vectorized — much faster than a for loop
  vector[N] mu = beta_month * month + beta_year * year_c
                 + alpha[province];
  log_size ~ normal(mu, sigma);
}

generated quantities {
  vector[N] mu = beta_0 + beta_month * month + beta_year * year_c
                 + alpha[province];
  vector[N] log_lik;
  vector[N] y_rep;
  for (i in 1:N) {
    log_lik[i] = normal_lpdf(log_size[i] | mu[i], sigma);
    y_rep[i]   = normal_rng(mu[i], sigma);
  }
}
