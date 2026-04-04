// Stan Model: Weibull AFT Mixture Cure Model
// Modificado: Dinámica condicional + Parametrización No Centrada (NCP) para evitar divergencias.

data {
  int<lower=1> N;
  vector<lower=0>[N] tiempo;
  int<lower=0, upper=1> evento[N];
  vector[N] arm;
  
  int<lower=1, upper=5> cure_prior_type;
  int<lower=0, upper=1> use_historical_prior;
  vector[2] historical_prior_params;
  int<lower=0, upper=1> shared_shape;
}

parameters {
  real beta_cure_intercept;
  real beta_cure_arm;       
  real beta_surv_intercept;
  
  // Muestreo en el espacio base limpio para evitar divergencias (NCP)
  real beta_surv_arm_raw;       
  
  real<lower=0> alpha_control;
  vector[1 - shared_shape] beta_alpha_arm_raw;
}

transformed parameters {
  real beta_alpha_arm;
  real beta_surv_arm;
  
  if (shared_shape == 1) {
    beta_alpha_arm = 0.0;
  } else {
    beta_alpha_arm = beta_alpha_arm_raw[1];
  }
  
  // DINÁMICA CONDICIONAL CON NCP:
  if (cure_prior_type == 2 || cure_prior_type == 3) {
    // Escenarios Escépticos: OR está hundida.
    // Multiplicador de 2.5 expande el parámetro raw para permitir la hipertrofia del TR.
    beta_surv_arm = 2.5 * beta_surv_arm_raw; 
  } else {
    // Escenarios Neutros o Direccionales: OR libre.
    // Multiplicador de 1.0 restringe la expansión para evitar inestabilidad.
    beta_surv_arm = 1.0 * beta_surv_arm_raw; 
  }
}

model {
  // 1. Priors for Survival (AFT) and Shape Parameters
  beta_surv_intercept ~ student_t(4, 0, 2.5);
  
  // El MCMC siempre ve una distribución dócil y fácil de explorar sin divergencias
  beta_surv_arm_raw ~ std_normal(); 
  
  alpha_control ~ gamma(1, 1);
  if (shared_shape == 0) {
    beta_alpha_arm_raw ~ normal(0, 1);
  }
  
  // 2. Priors for Cure Parameters
  beta_cure_intercept ~ student_t(4, 0, 2.5);

  if (use_historical_prior == 1) {
    beta_cure_arm ~ normal(historical_prior_params[1], historical_prior_params[2]);
  } else {
    if (cure_prior_type == 1) {
      beta_cure_arm ~ std_normal();
    } else if (cure_prior_type == 2) {
      beta_cure_arm ~ double_exponential(0, 0.2); 
    } else if (cure_prior_type == 3) {
      beta_cure_arm ~ double_exponential(0, 0.1);
    } else if (cure_prior_type == 4) {
      if (beta_cure_arm < 0) {
        target += negative_infinity();
      } else {
        beta_cure_arm ~ weibull(2.0, 0.707);
      }
    } else if (cure_prior_type == 5) {
      if (beta_cure_arm < 0) {
        target += negative_infinity();
      } else {
        beta_cure_arm ~ weibull(1.5, 1.0);
      }
    }
  }

  // 3. Likelihood
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


