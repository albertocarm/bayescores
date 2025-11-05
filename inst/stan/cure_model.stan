// Stan Model: cure_model.stan
// Implements a Weibull AFT (Accelerated Failure Time) mixture cure model.
//
// This model assumes the population is a mix of:
// 1. "Cured" patients (probability pi) who will never experience the event.
// 2. "Susceptible" patients (probability 1-pi) whose event time follows
//    a Weibull survival distribution.
//
// The model includes conditional priors to handle non-cure scenarios
// robustly by collapsing to a standard AFT model.

data {
  int<lower=1> N;                   // Number of subjects
  vector<lower=0>[N] tiempo;        // Observed time (event or censoring)
  int<lower=0, upper=1> evento[N];  // Event status (1=event, 0=censored)
  vector[N] arm;                    // Treatment arm (0=control, 1=treatment)

  // --- Conditional Prior Flag ---
  // Passed from R (as integer 1 or 0)
  // 1 = (TRUE)  Suspect cure: Use standard, weakly informative priors.
  // 0 = (FALSE) No cure:     Use skeptical priors to improve identifiability.
  int<lower=0, upper=1> suspect_cure;
}

parameters {
  // --- Cure Component Parameters ---
  real beta_cure_intercept;     // Log-odds of cure for arm=0
  real beta_cure_arm_raw;       // Non-centered effect of 'arm' on cure log-odds

  // --- Survival Component Parameters ---
  real beta_surv_intercept;     // Log-scale of survival for arm=0 (non-cured group)
  real beta_surv_arm_raw;       // Non-centered effect of 'arm' on survival log-scale

  // --- Weibull Shape Parameter ---
  real<lower=0> alpha;          // Weibull shape (common to all subjects)
}

transformed parameters {
  // Re-scale the non-centered parameters to get interpretable effects.
  // This non-centered parameterization (NCP) helps Stan's sampler
  // avoid divergent transitions, especially when the scale (2.5) is large.
  
  // Effective prior: beta_cure_arm ~ normal(0, 2.5)
  real beta_cure_arm = 2.5 * beta_cure_arm_raw;
  
  // Effective prior: beta_surv_arm ~ normal(0, 2.5)
  real beta_surv_arm = 2.5 * beta_surv_arm_raw;
}

model {
  // --- Priors for Survival Component (Always applied) ---
  beta_surv_intercept ~ student_t(4, 0, 2.5); // Robust prior for intercept
  beta_surv_arm_raw ~ std_normal();            // Prior on non-centered param
  alpha ~ gamma(1, 1);                         // Standard prior for shape > 0

  // --- Priors for Cure Component (Conditional) ---
  if (suspect_cure == 1) {
    // --- CASE 1: Cure is Suspected (Standard Model) ---
    // Use weakly informative, robust priors centered at 0 (log-odds).
    beta_cure_intercept ~ student_t(4, 0, 2.5);
    beta_cure_arm_raw ~ std_normal(); // (Effective prior is normal(0, 2.5))

  } else {
    // --- CASE 2: No Cure is Suspected (Skeptical Model) ---
    // Use skeptical priors to force the cure probability (pi) to ~0.
    // This resolves identifiability issues and collapses the mixture
    // model into a standard AFT model.

    // 1. Force intercept log-odds to be very negative (e.g., -10).
    //    inv_logit(-10) is ~0.000045, effectively forcing pi -> 0.
    beta_cure_intercept ~ normal(-10, 1);

    // 2. Force arm effect on cure to be zero.
    //    We apply a Laplace (L1) shrinkage prior.
    //    We want: beta_cure_arm ~ double_exponential(0, 1)
    //    On the raw param: 2.5 * beta_cure_arm_raw ~ double_exponential(0, 1)
    //    This is equivalent to: beta_cure_arm_raw ~ double_exponential(0, 1 / 2.5)
    beta_cure_arm_raw ~ double_exponential(0, 0.4);
  }

  // --- Likelihood ---
  // This logic is identical for both models.
  for (i in 1:N) {
    // 1. Calculate the log-odds of being CURED for subject i
    real cure_logit = beta_cure_intercept + beta_cure_arm * arm[i];

    // 2. Calculate the Weibull log-scale for subject i
    //    (Note: Stan's weibull uses scale, not log(scale))
    real surv_log_scale = exp(beta_surv_intercept + beta_surv_arm * arm[i]);

    if (evento[i] == 1) {
      // --- EVENT OBSERVED ---
      // This person *cannot* be in the "cured" group.
      // log-likelihood = log(Prob of NOT Cured) + log(Weibull PDF at time t)
      
      target += bernoulli_logit_lpmf(0 | cure_logit) + // log(1 - pi_i)
                weibull_lpdf(tiempo[i] | alpha, surv_log_scale);

    } else {
      // --- CENSORED ---
      // This person is *either* (A) in the "cured" group
      // *or* (B) in the "susceptible" group but hasn't had the event yet.
      // log-likelihood = log( Prob(A) + Prob(B) )
      // log-likelihood = log( pi_i + (1 - pi_i) * S_weibull(t) )
      
      target += log_sum_exp(
        bernoulli_logit_lpmf(1 | cure_logit), // log(pi_i)
        bernoulli_logit_lpmf(0 | cure_logit) + // log(1 - pi_i)
          weibull_lccdf(tiempo[i] | alpha, surv_log_scale) // log(S_weibull(t))
      );
    }
  }
}
