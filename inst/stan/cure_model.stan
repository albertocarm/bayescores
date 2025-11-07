// Stan Model: Weibull AFT Mixture Cure Model
//
// Model structure:
// 1. A logistic regression component models the probability of being "cured".
// 2. A Weibull AFT survival component models the time-to-event for the
//    "susceptible" (non-cured) population.
//
// This implementation includes a conditional skeptical prior on the
// adjuvant treatment effect (the log-OR for cure).

data {
  int<lower=1> N;                   // Number of subjects
  vector<lower=0>[N] tiempo;        // Observed time (event or censoring)
  int<lower=0, upper=1> evento[N];  // Event status (1=event, 0=censored)
  vector[N] arm;                    // Treatment arm (0=control, 1=treatment)

  // Flag for conditional priors on the adjuvant cure effect (log-OR)
  // 1 = (TRUE)  Standard, weakly-informative prior (default)
  // 0 = (FALSE) Skeptical shrinkage prior (Laplace) on the arm effect
  int<lower=0, upper=1> suspect_cure;
}

parameters {
  // --- Cure Component Parameters ---
  // Log-odds of cure for control arm (e.g., surgery alone)
  real beta_cure_intercept;
  
  // Non-centered effect of 'arm' on cure log-odds
  real beta_cure_arm_raw;

  // --- Survival Component Parameters (for non-cured population) ---
  // Log-scale of survival for control arm
  real beta_surv_intercept;
  
  // Non-centered effect of 'arm' on survival log-scale
  real beta_surv_arm_raw;

  // --- Shared Weibull Shape Parameter ---
  real<lower=0> alpha;
}

transformed parameters {
  // Re-scale non-centered parameters (NCP) for interpretability
  // and improved sampler stability.

  // beta_cure_arm: log-OR of cure for treatment vs. control
  // Effective prior (if suspect_cure=1): normal(0, 2.5)
  real beta_cure_arm = 2.5 * beta_cure_arm_raw;

  // beta_surv_arm: log-Time Ratio for treatment vs. control
  // Effective prior: normal(0, 2.5)
  real beta_surv_arm = 2.5 * beta_surv_arm_raw;
}

model {
  // --- Priors ---

  // Priors for survival component (non-cured population)
  // These are unconditional and identical to the original model.
  beta_surv_intercept ~ student_t(4, 0, 2.5);
  beta_surv_arm_raw ~ std_normal();
  alpha ~ gamma(1, 1);

  // Prior for cure intercept (control arm)
  // This is unconditional, reflecting the baseline (e.g., surgical) cure rate.
  // It is identical to the original model.
  beta_cure_intercept ~ student_t(4, 0, 2.5);
  
  // Conditional prior for the adjuvant cure *effect*
  if (suspect_cure == 1) {
    // Standard, weakly-informative prior (identical to original)
    beta_cure_arm_raw ~ std_normal();
  
  } else {
    // Skeptical (L1 shrinkage) prior on the arm effect.
    // This is equivalent to beta_cure_arm ~ double_exponential(0, 1)
    // on the transformed scale (since 0.4 = 1 / 2.5).
    beta_cure_arm_raw ~ double_exponential(0, 0.4);
  }

  // --- Likelihood ---
  // This section is identical to the original model.
  for (i in 1:N) {
    // Log-odds of being cured
    real cure_logit = beta_cure_intercept + beta_cure_arm * arm[i];

    // Weibull scale parameter (on log scale)
    // Note: Stan's weibull_lpdf uses the *scale*, not log-scale
    real surv_log_scale = exp(beta_surv_intercept + beta_surv_arm * arm[i]);

    if (evento[i] == 1) {
      // --- EVENT ---
      // Subject must be "susceptible" (not cured).
      // log(Prob(Not Cured)) + log(Weibull PDF at time t)
      target += bernoulli_logit_lpmf(0 | cure_logit) +
                weibull_lpdf(tiempo[i] | alpha, surv_log_scale);
    } else {
      // --- CENSORED ---
      // Subject is *either* "cured" *or* "susceptible but alive".
      // log( Prob(Cured) + Prob(Not Cured) * S_weibull(t) )
      target += log_sum_exp(
        bernoulli_logit_lpmf(1 | cure_logit), // log(Prob(Cured))
        bernoulli_logit_lpmf(0 | cure_logit) + // log(Prob(Not Cured))
          weibull_lccdf(tiempo[i] | alpha, surv_log_scale) // log(S(t))
      );
    }
  }
}

