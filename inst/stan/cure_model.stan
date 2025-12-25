// Stan Model: Weibull AFT Mixture Cure Model (Optimized Shape Version)

data {
  int<lower=1> N;
  vector<lower=0>[N] tiempo;
  int<lower=0, upper=1> evento[N];
  vector[N] arm;
  int<lower=1, upper=3> cure_prior_type;
  int<lower=0, upper=1> use_historical_prior;
  vector[2] historical_prior_params; 
  int<lower=0, upper=1> shared_shape; // 1 = shared, 0 = independent
}

parameters {
  real beta_cure_intercept;
  real beta_cure_arm_raw;
  real beta_surv_intercept;
  real beta_surv_arm_raw;
  real<lower=0> alpha_control;
  
  // MODIFICATION: Vector of size 0 if shared_shape is 1
  vector[1 - shared_shape] beta_alpha_arm_raw;
}

transformed parameters {
  real beta_cure_arm = 2.5 * beta_cure_arm_raw;
  real beta_surv_arm = 2.5 * beta_surv_arm_raw;
  
  // MODIFICATION: Define beta_alpha_arm as a constant 0 if shared
  real beta_alpha_arm;
  if (shared_shape == 1) {
    beta_alpha_arm = 0.0;
  } else {
    beta_alpha_arm = beta_alpha_arm_raw[1];
  }
}

model {
  beta_surv_intercept ~ student_t(4, 0, 2.5);
  beta_surv_arm_raw ~ std_normal();
  alpha_control ~ gamma(1, 1);

  // MODIFICATION: Only apply prior if the parameter exists
  if (shared_shape == 0) {
    beta_alpha_arm_raw ~ normal(0, 1);
  }

  beta_cure_intercept ~ student_t(4, 0, 2.5);
    
  if (use_historical_prior == 1) {
    beta_cure_arm_raw ~ normal(
        historical_prior_params[1] / 2.5, 
        historical_prior_params[2] / 2.5
    );
  } else {
    if (cure_prior_type == 1) {
      beta_cure_arm_raw ~ std_normal();
    } else if (cure_prior_type == 2) {
      beta_cure_arm_raw ~ double_exponential(0, 0.08); 
    } else {
      beta_cure_arm_raw ~ double_exponential(0, 0.04);
    }
  }

  for (i in 1:N) {
    real cure_logit = beta_cure_intercept + beta_cure_arm * arm[i];
    real surv_log_scale = exp(beta_surv_intercept + beta_surv_arm * arm[i]);

    // MODIFICATION: Logic is now clean because beta_alpha_arm is 0 if shared
    real alpha_i = alpha_control * exp(beta_alpha_arm * arm[i]);

    if (evento[i] == 1) {
      target += bernoulli_logit_lpmf(0 | cure_logit) +
                weibull_lpdf(tiempo[i] | alpha_i, surv_log_scale);
    } else {
      target += log_sum_exp(
        bernoulli_logit_lpmf(1 | cure_logit), 
        bernoulli_logit_lpmf(0 | cure_logit) + 
          weibull_lccdf(tiempo[i] | alpha_i, surv_log_scale)
      );
    }
  }
}
