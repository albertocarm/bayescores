// Stan Model: Weibull AFT Mixture Cure Model
// Modified for 'Tail Assumption' logic

data {
  int<lower=1> N;
  vector<lower=0>[N] tiempo;
  int<lower=0, upper=1> evento[N];
  vector[N] arm;
  // cure_prior_type mapping:
  // 1: neutral (Data-Driven / Agnostic)
  // 2: immature_skeptical (Safety Brake / Skeptical)
  // 3: biologically_null (Structural AFT / Null)
  // 4: supportive (Phase II Evidence / Mild)
  // 5: optimistic (Strong Evidence)
  int<lower=1, upper=5> cure_prior_type; 
  int<lower=0, upper=1> use_historical_prior;
  vector[2] historical_prior_params; 
  int<lower=0, upper=1> shared_shape; 
}

parameters {
  real beta_cure_intercept;
  real beta_cure_arm_raw; 
  real beta_surv_intercept;
  real beta_surv_arm_raw;
  real<lower=0> alpha_control;
  vector[1 - shared_shape] beta_alpha_arm_raw;
}

transformed parameters {
  // Scaling factors kept exactly as original
  real beta_cure_arm = 2.5 * beta_cure_arm_raw;
  real beta_surv_arm = 2.5 * beta_surv_arm_raw;
    
  real beta_alpha_arm;
  if (shared_shape == 1) {
    beta_alpha_arm = 0.0;
  } else {
    beta_alpha_arm = beta_alpha_arm_raw[1];
  }
}

model {
  // Priors for survival and shape parameters (Unchanged)
  beta_surv_intercept ~ student_t(4, 0, 2.5);
  beta_surv_arm_raw ~ std_normal();
  alpha_control ~ gamma(1, 1);

  if (shared_shape == 0) {
    beta_alpha_arm_raw ~ normal(0, 1);
  }

  beta_cure_intercept ~ student_t(4, 0, 2.5);
     
  // Logic for Cure Arm Prior (Modified labels and order)
  if (use_historical_prior == 1) {
    beta_cure_arm_raw ~ normal(
        historical_prior_params[1] / 2.5, 
        historical_prior_params[2] / 2.5
    );
  } else {
    if (cure_prior_type == 1) {
      // 1. NEUTRAL (Agostic / Data-Driven)
      // Corresponds to old 'unknown'
      beta_cure_arm_raw ~ std_normal();
      
    } else if (cure_prior_type == 2) {
      // 2. IMMATURE_SKEPTICAL (Safety Brake)
      // Corresponds to old 'unlikely'. 
      // Using Laplace with scale 0.08 (strong shrinkage for immature data)
      beta_cure_arm_raw ~ double_exponential(0, 0.08); 
      
    } else if (cure_prior_type == 3) {
      // 3. BIOLOGICALLY_NULL (Structural AFT)
      // Corresponds to old 'very_unlikely'.
      // Using Laplace with scale 0.04 (extreme shrinkage/collapse)
      beta_cure_arm_raw ~ double_exponential(0, 0.04);
      
    } else if (cure_prior_type == 4) {
      // 4. SUPPORTIVE (Mild Evidence)
      // Corresponds to old 'mild_optimistic' (previously type 5)
      // Logic swapped to match the new vector order.
      if (beta_cure_arm_raw < 0) {
        target += negative_infinity();
      } else {
        // Weibull parameters for mild/supportive prior
        beta_cure_arm_raw ~ weibull(2.0, 0.2828);
      }
      
    } else if (cure_prior_type == 5) {
      // 5. OPTIMISTIC (Strong Evidence)
      // Corresponds to old 'optimistic' (previously type 4)
      if (beta_cure_arm_raw < 0) {
        target += negative_infinity();
      } else {
        // Weibull parameters for optimistic/strong prior
        beta_cure_arm_raw ~ weibull(1.5, 0.4);
      }
    }
  }

  // Likelihood (Unchanged)
  for (i in 1:N) {
    real cure_logit = beta_cure_intercept + beta_cure_arm * arm[i];
    real surv_log_scale = exp(beta_surv_intercept + beta_surv_arm * arm[i]);
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


