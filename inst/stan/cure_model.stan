data {
  int<lower=1> N;
  vector<lower=0>[N] tiempo;
  int<lower=0, upper=1> evento[N];
  vector[N] arm;
}

parameters {
  // Non-centered parameterization for stability
  real beta_cure_intercept;
  real beta_surv_intercept;
  real beta_cure_arm_raw;
  real beta_surv_arm_raw;
  real<lower=0> alpha;
}

transformed parameters {
  // Scaling the raw parameters to get interpretable effects
  real beta_cure_arm = 2.5 * beta_cure_arm_raw;
  real beta_surv_arm = 2.5 * beta_surv_arm_raw;
}

model {
  // Priors on base parameters
  beta_cure_intercept ~ student_t(4, 0, 2.5);
  beta_surv_intercept ~ student_t(4, 0, 2.5);
  alpha ~ gamma(1, 1);
  beta_cure_arm_raw ~ std_normal();
  beta_surv_arm_raw ~ std_normal();

  // Likelihood with corrected logic
  for (i in 1:N) {
    // cure_logit is the log-odds of being CURED
    real cure_logit = beta_cure_intercept + beta_cure_arm * arm[i];

    // surv_log_scale for the Weibull
    real surv_log_scale = exp(beta_surv_intercept + beta_surv_arm * arm[i]);

    if (evento[i] == 1) {
      // --- EVENT ---
      // log(Prob of NOT being cured) + log(Weibull PDF at time t)
      // bernoulli_logit_lpmf(0 | cure_logit) is log(1 - sigmoid(cure_logit))
      target += bernoulli_logit_lpmf(0 | cure_logit) +
                weibull_lpdf(tiempo[i] | alpha, surv_log_scale);
    } else {
      // --- CENSORED ---
      // log( Prob(Cured) + Prob(Not Cured) * S_weibull(t) )
      target += log_sum_exp(
        bernoulli_logit_lpmf(1 | cure_logit), // log(Prob of being CURED)
        bernoulli_logit_lpmf(0 | cure_logit) + weibull_lccdf(tiempo[i] | alpha, surv_log_scale) // log(Prob of NOT CURED) + log(Survival of non-cured)
      );
    }
  }
}


