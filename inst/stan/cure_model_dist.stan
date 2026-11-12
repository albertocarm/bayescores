// Stan Model: AFT Mixture Cure Model with Alternative Latency Distributions
// dist: 1 exponential, 2 log-normal, 3 log-logistic, 4 gamma, 5 generalized gamma

functions {
  // Regularized incomplete gamma at x = a * exp(lx), Temme uniform expansion
  real cure_temme(real a, real lx, int is_upper) {
    real mu = expm1(lx);
    real s2 = 1 / a;
    real r;
    real eta;
    real z;
    real c0;
    real c1;
    real c2;
    real c3;
    real corr;
    if (lx > -0.01 && lx < 0.01) {
      r = 0.5 + mu * (-0.3333333333333333 + mu * (0.25 + mu * (-0.2
          + mu * (0.16666666666666666 + mu * (-0.14285714285714285 + mu * 0.125)))));
    } else {
      r = (1 - lx / mu) / mu;
    }
    eta = mu * sqrt(2 * r);
    z = eta * sqrt(a);
    if (z > 38) {
      if (is_upper == 1) {
        return 0.0;
      }
      return 1.0;
    }
    if (z < -38) {
      if (is_upper == 1) {
        return 1.0;
      }
      return 0.0;
    }
    if (eta > -0.02 && eta < 0.02) {
      c0 = -0.3333333333333333 + eta * (0.08333333333333333 + eta * (-0.014814814814814815
           + eta * (0.0011574074074074073 + eta * (0.0003527336860670194
           + eta * (-0.0001787551440329218 + eta * (0.00003919263178522438
           - eta * 0.000002185448510679992))))));
    } else {
      c0 = 1 / mu - 1 / eta;
    }
    if (eta > -0.1 && eta < 0.1) {
      c1 = -0.001851851851851852 + eta * (-0.003472222222222222 + eta * (0.0026455026455026454
           + eta * (-0.0009902263374485596 + eta * (0.00020576131687242798
           + eta * (-0.0000004018775720164609 + eta * (-0.000018098550334489977
           + eta * 0.00000764916091608111))))));
    } else {
      c1 = 1 / eta ^ 3 - 1 / mu ^ 3 - 1 / mu ^ 2 - 1 / (12 * mu);
    }
    if (eta > -0.25 && eta < 0.25) {
      c2 = 0.004133597883597883 + eta * (-0.0026813271604938273 + eta * (0.0007716049382716049
           + eta * (0.0000020093878600823047 + eta * (-0.0001073665322636516
           + eta * (0.000052923448829120125 + eta * (-0.000012760635188618728
           + eta * 0.00000003423578734096138))))));
    } else {
      c2 = -3 / eta ^ 5 + 3 / mu ^ 5 + 5 / mu ^ 4 + 25 / (12 * mu ^ 3) + 1 / (12 * mu ^ 2)
           + 1 / (288 * mu);
    }
    if (eta > -0.4 && eta < 0.4) {
      c3 = 0.0006494341563786008 + eta * (0.00022947209362139917 + eta * (-0.0004691894943952557
           + eta * (0.00026772063206283885 + eta * (-0.00007561801671883977
           + eta * (-0.00000023965051138672968 + eta * (0.000011082654115347302
           - eta * 0.0000056749528269915965))))));
    } else {
      c3 = 15 / eta ^ 7 - 15 / mu ^ 7 - 35 / mu ^ 6 - 105 / (4 * mu ^ 5) - 77 / (12 * mu ^ 4)
           - 49 / (288 * mu ^ 3) - 1 / (288 * mu ^ 2) + 139 / (51840 * mu);
    }
    corr = exp(-0.5 * z * z) * 0.3989422804014327 * sqrt(s2)
           * (c0 + s2 * (c1 + s2 * (c2 + s2 * c3)));
    if (is_upper == 1) {
      return 0.5 * erfc(z * 0.7071067811865476) + corr;
    }
    return 0.5 * erfc(-z * 0.7071067811865476) - corr;
  }

  // Regularized incomplete gamma at x = a * exp(lx), any shape
  real cure_inc_gamma(real a, real lx, int is_upper) {
    real aa = a;
    int m = 0;
    real logx;
    real x;
    real term;
    real acc = 0;
    real out;
    while (aa < 20) {
      aa += 1;
      m += 1;
    }
    if (m == 0) {
      return cure_temme(a, lx, is_upper);
    }
    logx = log(a) + lx;
    out = cure_temme(aa, logx - log(aa), is_upper);
    if (is_upper == 1 && out <= 0) {
      return 0.0;
    }
    if (is_upper == 0 && out >= 1) {
      return 1.0;
    }
    x = exp(logx);
    term = exp(a * logx - x - lgamma(a + 1));
    for (j in 1:m) {
      acc += term;
      term *= x / (a + j);
    }
    if (is_upper == 1) {
      return out - acc;
    }
    return out + acc;
  }

  // Log density of the standardized generalized gamma
  real cure_gengamma_logdens(real w, real q) {
    real q2 = q * q;
    real u = q * w;
    real k;
    real c;
    if (u > 700) {
      u = 700;
    }
    if (u > -0.01 && u < 0.01) {
      k = 0.5 + u * (0.16666666666666666 + u * (0.041666666666666664
          + u * (0.008333333333333333 + u * (0.001388888888888889 + u * 0.0001984126984126984))));
    } else {
      k = (expm1(u) - u) / (u * u);
    }
    if (q2 < 0.05) {
      c = -0.9189385332046727 - q2 / 12 + q2 ^ 3 / 360 - q2 ^ 5 / 1260;
    } else {
      c = 0.5 * log(q2) - (log(q2) + 1) / q2 - lgamma(1 / q2);
    }
    return c - w * w * k;
  }

  // Survival function of the standardized generalized gamma
  real cure_gengamma_surv(real w, real q) {
    real q2 = q * q;
    real u = q * w;
    if (q2 < 1e-200) {
      return 0.5 * erfc(w * 0.7071067811865476);
    }
    if (u > 700) {
      u = 700;
    }
    if (q > 0) {
      return cure_inc_gamma(1 / q2, u, 1);
    }
    return cure_inc_gamma(1 / q2, u, 0);
  }
}

data {
  int<lower=1> N;
  vector<lower=0>[N] tiempo;
  vector<lower=0, upper=1>[N] evento;
  vector[N] arm;

  int<lower=1, upper=3> cure_prior_type;
  int<lower=0, upper=1> use_historical_prior;
  vector[2] historical_prior_params;
  int<lower=0, upper=1> shared_shape;
  int<lower=1, upper=5> dist;
}

transformed data {
  vector[N] log_t = log(tiempo);
  int has_shape = (dist == 1) ? 0 : 1;
  int has_q = (dist == 5) ? 1 : 0;
}

parameters {
  real beta_cure_intercept;
  real beta_cure_arm;
  real beta_surv_intercept;

  real beta_surv_arm_raw;

  vector<lower=0>[has_shape] alpha_control_raw;
  vector[has_shape * (1 - shared_shape)] beta_alpha_arm_raw;
  vector[has_q] q_shape;
}

transformed parameters {
  real alpha_control;
  real beta_alpha_arm;
  real beta_surv_arm;

  if (has_shape == 1) {
    alpha_control = alpha_control_raw[1];
  } else {
    alpha_control = 1.0;
  }

  if (has_shape == 1 && shared_shape == 0) {
    beta_alpha_arm = beta_alpha_arm_raw[1];
  } else {
    beta_alpha_arm = 0.0;
  }

  beta_surv_arm = 2.5 * beta_surv_arm_raw;
}

model {
  real q = 0;
  if (has_q == 1) {
    q = q_shape[1];
  }

  // 1. Priors for Survival (AFT) and Shape Parameters
  beta_surv_intercept ~ student_t(4, 0, 2.5);

  beta_surv_arm_raw ~ std_normal();

  alpha_control_raw ~ gamma(1, 1);
  beta_alpha_arm_raw ~ normal(0, 1);
  q_shape ~ std_normal();

  // 2. Priors for Cure Parameters
  beta_cure_intercept ~ student_t(4, 0, 2.5);

  if (use_historical_prior == 1) {
    beta_cure_arm ~ normal(historical_prior_params[1], historical_prior_params[2]);
  } else {
    if (cure_prior_type == 1) {
      beta_cure_arm ~ normal(0, 2.5);
    } else if (cure_prior_type == 2) {
      beta_cure_arm ~ double_exponential(0, 0.2);
    } else if (cure_prior_type == 3) {
      beta_cure_arm ~ double_exponential(0, 0.1);
    }
  }

  // 3. Likelihood
  for (i in 1:N) {
    real cure_logit = beta_cure_intercept + beta_cure_arm * arm[i];
    real mu_i = beta_surv_intercept + beta_surv_arm * arm[i];
    real alpha_i = alpha_control * exp(beta_alpha_arm * arm[i]);
    real w = (log_t[i] - mu_i) * alpha_i;

    if (evento[i] == 1) {
      real lf;
      if (dist == 1) {
        lf = -mu_i - exp(w);
      } else if (dist == 2) {
        lf = log(alpha_i) - log_t[i] - 0.9189385332046727 - 0.5 * w * w;
      } else if (dist == 3) {
        lf = log(alpha_i) - log_t[i] + w - 2 * log1p_exp(w);
      } else if (dist == 4) {
        real logx = log(alpha_i) + log_t[i] - mu_i;
        lf = alpha_i * logx - exp(logx) - lgamma(alpha_i) - log_t[i];
      } else {
        lf = log(alpha_i) - log_t[i] + cure_gengamma_logdens(w, q);
      }
      target += bernoulli_logit_lpmf(0 | cure_logit) + lf;
    } else {
      real s;
      if (dist == 1) {
        s = exp(-exp(w));
      } else if (dist == 2) {
        s = 0.5 * erfc(w * 0.7071067811865476);
      } else if (dist == 3) {
        s = inv_logit(-w);
      } else if (dist == 4) {
        s = cure_inc_gamma(alpha_i, fmin(log_t[i] - mu_i, 700), 1);
      } else {
        s = cure_gengamma_surv(w, q);
      }
      target += log(inv_logit(cure_logit) + inv_logit(-cure_logit) * fmax(s, 0));
    }
  }
}
