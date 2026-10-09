// int_fit_nb_rate_hill_bidir_splitnH.stan
// Bidirectional ATN blocking-probability fit, SPLIT-nH variant: NB
// rate-thinning likelihood x Hill decay on mosquito-level oocyst counts,
// with separate nH and s_half on each side of infection.  Same shared
// peak (b_max = 1 - r_min).  Compare against int_fit_nb_rate_hill_bidir.stan
// (which shares nH between sides) using LOO.
//
// Decay shape (s = absolute time of exposure relative to infection):
//   side = 0 (pre):  kernel = s_half_pre^nH_pre   / (s_half_pre^nH_pre   + s^nH_pre)
//   side = 1 (post): kernel = s_half_post^nH_post / (s_half_post^nH_post + s^nH_post)
// Continuity at s = 0 holds independently of nH: both kernels equal 1 at
// s = 0, so b(0-) = b(0+) = 1 - r_min on either side.
//
// Likelihood: identical to int_fit_nb_rate_hill_bidir.stan -- NB
// rate-thinning with shared dispersion phi and a per-experiment log-rate
// OLRE.  log_lik aggregates over all mosquitoes within each row n, so it
// is directly comparable to the shared-nH bidir fit via loo_compare().

data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> M;
  array[N] int<lower=1, upper=K> exp_id;
  vector<lower=0>[N] s;
  array[N] int<lower=0, upper=1> side;     // 0 = pre, 1 = post
  array[M] int<lower=1, upper=N> row_id;
  array[M] int<lower=0, upper=1> arm;
  array[M] int<lower=0> y;

  // Prior hyperparameters -- canonical values live in priors_atn_main.R.
  // Drivers source that file and pass the values via prior_stan_data().
  // Do NOT hardcode prior constants in this file.
  real<lower=0> r_min_alpha;
  real<lower=0> r_min_beta;
  real log_s_half_mean;
  real<lower=0> log_s_half_sd;
  real log_nH_mean;
  real<lower=0> log_nH_sd;
  real mu_L_mean;
  real<lower=0> mu_L_sd;
  real<lower=0> sigma_exp_sd;
  real<lower=0> reciprocal_phi_sd;
}

parameters {
  real mu_L;
  real<lower=0, upper=1> r_min;
  real log_s_half_pre;
  real log_s_half_post;
  real log_nH_pre;
  real log_nH_post;
  real<lower=0> sigma_exp;
  vector[K] z_exp;
  real<lower=0> reciprocal_phi;
}

transformed parameters {
  real<lower=0> s_half_pre  = exp(log_s_half_pre);
  real<lower=0> s_half_post = exp(log_s_half_post);
  real<lower=0> nH_pre      = exp(log_nH_pre);
  real<lower=0> nH_post     = exp(log_nH_post);
  vector[K] beta_exp        = sigma_exp * z_exp;
  real<lower=0> phi         = inv(reciprocal_phi);
}

model {
  mu_L             ~ normal(mu_L_mean, mu_L_sd);
  r_min            ~ beta(r_min_alpha, r_min_beta);
  log_s_half_pre   ~ normal(log_s_half_mean, log_s_half_sd);
  log_s_half_post  ~ normal(log_s_half_mean, log_s_half_sd);
  log_nH_pre       ~ normal(log_nH_mean,     log_nH_sd);
  log_nH_post      ~ normal(log_nH_mean,     log_nH_sd);
  sigma_exp        ~ normal(0, sigma_exp_sd);
  z_exp            ~ std_normal();
  reciprocal_phi   ~ normal(0, reciprocal_phi_sd);

  real shn_pre  = pow(s_half_pre,  nH_pre);
  real shn_post = pow(s_half_post, nH_post);
  vector[N] r_s;
  for (n in 1:N) {
    real kern;
    if (side[n] == 0) {
      kern = shn_pre  / (shn_pre  + pow(s[n], nH_pre));
    } else {
      kern = shn_post / (shn_post + pow(s[n], nH_post));
    }
    r_s[n] = 1 - (1 - r_min) * kern;
  }

  for (m in 1:M) {
    int n     = row_id[m];
    real L_c  = exp(mu_L + beta_exp[exp_id[n]]);
    real mu_m = (arm[m] == 0) ? L_c : r_s[n] * L_c;
    target += neg_binomial_2_lpmf(y[m] | mu_m, phi);
  }
}

generated quantities {
  real b_max = 1 - r_min;

  vector[N] r_vals;
  vector[N] L_c_row;
  array[M] int y_rep;
  vector[N] log_lik = rep_vector(0, N);

  {
    real shn_pre  = pow(s_half_pre,  nH_pre);
    real shn_post = pow(s_half_post, nH_post);
    for (n in 1:N) {
      real kern;
      if (side[n] == 0) {
        kern = shn_pre  / (shn_pre  + pow(s[n], nH_pre));
      } else {
        kern = shn_post / (shn_post + pow(s[n], nH_post));
      }
      r_vals[n]  = 1 - (1 - r_min) * kern;
      L_c_row[n] = exp(mu_L + beta_exp[exp_id[n]]);
    }
  }

  for (m in 1:M) {
    int n     = row_id[m];
    real mu_m = (arm[m] == 0) ? L_c_row[n] : r_vals[n] * L_c_row[n];
    y_rep[m]    = neg_binomial_2_rng(mu_m, phi);
    log_lik[n] += neg_binomial_2_lpmf(y[m] | mu_m, phi);
  }
}
