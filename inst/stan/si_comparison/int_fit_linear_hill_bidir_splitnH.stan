// int_fit_linear_hill_bidir_splitnH.stan
// Bidirectional ATN blocking-probability fit, SPLIT-nH variant: linear
// (TBA) binomial likelihood on per-row prevalence x Hill decay, with
// separate nH and s_half on each side of infection.  Same shared peak
// (b_max = 1 - r_min).  Compare against int_fit_linear_hill_bidir.stan
// (which shares nH between sides) using LOO.
//
// Likelihood: identical to int_fit_linear_hill_bidir.stan (binomial on
// per-row prevalence with logit-scale OLRE).
//
// Identifiability caveat: the post-infection rows in the current dataset
// cover only a few unique s values, so (s_half_post, nH_post) may show a
// strong banana correlation in the joint posterior.  This is expected --
// the goal is to see whether freeing nH_post improves out-of-sample fit
// despite the wider posterior.

data {
  int<lower=1> N;
  int<lower=1> K;
  array[N] int<lower=1, upper=K> exp_id;
  vector<lower=0>[N] s;
  array[N] int<lower=0, upper=1> side;     // 0 = pre, 1 = post
  array[N] int<lower=0> n_ctl;
  array[N] int<lower=0> pos_ctl;
  array[N] int<lower=0> n_trt;
  array[N] int<lower=0> pos_trt;

  // Prior hyperparameters -- canonical values live in priors_atn_main.R.
  // Drivers source that file and pass the values via prior_stan_data().
  // Do NOT hardcode prior constants in this file.
  real<lower=0> r_min_alpha;
  real<lower=0> r_min_beta;
  real log_s_half_mean;
  real<lower=0> log_s_half_sd;
  real log_nH_mean;
  real<lower=0> log_nH_sd;
  real mu_p_mean;
  real<lower=0> mu_p_sd;
  real<lower=0> sigma_exp_sd;
}

parameters {
  real mu_p;
  real<lower=0, upper=1> r_min;
  real log_s_half_pre;
  real log_s_half_post;
  real log_nH_pre;
  real log_nH_post;
  real<lower=0> sigma_exp;
  vector[K] z_exp;
}

transformed parameters {
  real<lower=0> s_half_pre  = exp(log_s_half_pre);
  real<lower=0> s_half_post = exp(log_s_half_post);
  real<lower=0> nH_pre      = exp(log_nH_pre);
  real<lower=0> nH_post     = exp(log_nH_post);
  vector[K] beta_exp        = sigma_exp * z_exp;
}

model {
  mu_p             ~ normal(mu_p_mean, mu_p_sd);
  r_min            ~ beta(r_min_alpha, r_min_beta);
  log_s_half_pre   ~ normal(log_s_half_mean, log_s_half_sd);
  log_s_half_post  ~ normal(log_s_half_mean, log_s_half_sd);
  log_nH_pre       ~ normal(log_nH_mean,     log_nH_sd);
  log_nH_post      ~ normal(log_nH_mean,     log_nH_sd);
  sigma_exp        ~ normal(0, sigma_exp_sd);
  z_exp            ~ std_normal();

  real shn_pre  = pow(s_half_pre,  nH_pre);
  real shn_post = pow(s_half_post, nH_post);

  for (n in 1:N) {
    real p_ctl_n = inv_logit(mu_p + beta_exp[exp_id[n]]);
    real kern;
    if (side[n] == 0) {
      kern = shn_pre  / (shn_pre  + pow(s[n], nH_pre));
    } else {
      kern = shn_post / (shn_post + pow(s[n], nH_post));
    }
    real r_s     = 1 - (1 - r_min) * kern;
    real p_trt_n = r_s * p_ctl_n;
    target += binomial_lpmf(pos_ctl[n] | n_ctl[n], p_ctl_n);
    target += binomial_lpmf(pos_trt[n] | n_trt[n], p_trt_n);
  }
}

generated quantities {
  real b_max = 1 - r_min;

  array[N] int pos_ctl_rep;
  array[N] int pos_trt_rep;
  vector[N] p_ctl;
  vector[N] p_trt;
  vector[N] r_vals;
  vector[N] log_lik;

  {
    real shn_pre  = pow(s_half_pre,  nH_pre);
    real shn_post = pow(s_half_post, nH_post);
    for (n in 1:N) {
      real p_ctl_n = inv_logit(mu_p + beta_exp[exp_id[n]]);
      real kern;
      if (side[n] == 0) {
        kern = shn_pre  / (shn_pre  + pow(s[n], nH_pre));
      } else {
        kern = shn_post / (shn_post + pow(s[n], nH_post));
      }
      real r_s     = 1 - (1 - r_min) * kern;
      real p_trt_n = r_s * p_ctl_n;

      p_ctl[n]       = p_ctl_n;
      p_trt[n]       = p_trt_n;
      r_vals[n]      = r_s;
      pos_ctl_rep[n] = binomial_rng(n_ctl[n], p_ctl_n);
      pos_trt_rep[n] = binomial_rng(n_trt[n], p_trt_n);
      log_lik[n]     = binomial_lpmf(pos_ctl[n] | n_ctl[n], p_ctl_n)
                     + binomial_lpmf(pos_trt[n] | n_trt[n], p_trt_n);
    }
  }
}
