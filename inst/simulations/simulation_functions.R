# Additional simulation scenarios of the BayeScores manuscript (Scenarios 6
# to 12): scenario definitions and data generation.

# Scenario definitions ---------------------------------------------------
# True parameter values as in Scenarios 1 (benefit) and 2 (null) of Table 1.
base <- list(n = 200, median_ctrl = 12, max_fu = 120, weibull_shape = 1.5,
             latency = "weibull", delay = 0, dropout_annual = 0, accrual = 0)
benefit <- list(cure_ctrl = 0.30, cure_exp = 0.75, tr = 2)
null    <- list(cure_ctrl = 0.40, cure_exp = 0.40, tr = 1)
sc <- function(...) utils::modifyList(base, c(...))

scenarios <- list(
  # Non-Weibull survival times of the uncured patients.
  sc06_lognormal_benefit   = sc(benefit, list(latency = "lognormal")),
  sc07_lognormal_null      = sc(null,    list(latency = "lognormal")),
  sc08_loglogistic_benefit = sc(benefit, list(latency = "loglogistic")),
  sc09_loglogistic_null    = sc(null,    list(latency = "loglogistic")),
  # Delayed treatment effect: no effect on event times during the first 3 months.
  sc10_delayed_effect      = sc(benefit, list(delay = 3)),
  # Alternative censoring mechanisms.
  sc11_random_dropout      = sc(benefit, list(dropout_annual = 0.05)),
  sc12_staggered_accrual   = sc(benefit, list(accrual = 24))
)

# Latency distributions share the control median and the standard deviation
# of the log time of the Weibull distribution with shape 1.5.
sd_log_time <- function(shape) pi / (sqrt(6) * shape)

latency_quantile <- function(u, s) {
  sdl <- sd_log_time(s$weibull_shape)
  switch(s$latency,
    weibull     = stats::qweibull(u, shape = s$weibull_shape,
                                  scale = s$median_ctrl / log(2)^(1 / s$weibull_shape)),
    lognormal   = stats::qlnorm(u, meanlog = log(s$median_ctrl), sdlog = sdl),
    loglogistic = s$median_ctrl * (u / (1 - u))^(sdl * sqrt(3) / pi),
    stop("unknown latency"))
}

# One trial. The same seed gives the same underlying uniforms in every scenario.
generate_trial <- function(s, seed) {
  set.seed(seed)
  n <- s$n
  u <- matrix(stats::runif(5 * n), ncol = 5)
  exper <- u[, 1] < 0.5
  cured <- u[, 2] < ifelse(exper, s$cure_exp, s$cure_ctrl)
  t_ctrl <- latency_quantile(u[, 3], s)
  t_exp <- if (s$delay > 0) {
    ifelse(t_ctrl <= s$delay, t_ctrl, s$delay + s$tr * (t_ctrl - s$delay))
  } else {
    s$tr * t_ctrl
  }
  latent <- ifelse(exper, t_exp, t_ctrl)
  latent[cured] <- Inf
  cens <- rep(s$max_fu, n)
  if (s$accrual > 0) cens <- s$max_fu - s$accrual * u[, 4]
  if (s$dropout_annual > 0) {
    cens <- pmin(cens, stats::qexp(u[, 5], rate = -log(1 - s$dropout_annual) / 12))
  }
  data.frame(
    time  = pmin(latent, cens),
    event = as.numeric(latent <= cens),
    arm   = factor(ifelse(exper, "Experimental", "Control"),
                   levels = c("Control", "Experimental"))
  )
}

# True values (NA when the fitted Weibull parameter has no counterpart).
true_values <- function(s) {
  weib <- s$latency == "weibull"
  c(alpha_control       = if (weib && s$delay == 0) s$weibull_shape else NA,
    beta_cure_arm       = stats::qlogis(s$cure_exp) - stats::qlogis(s$cure_ctrl),
    beta_cure_intercept = stats::qlogis(s$cure_ctrl),
    beta_surv_arm       = if (s$delay == 0) log(s$tr) else NA,
    beta_surv_intercept = if (weib) log(s$median_ctrl / log(2)^(1 / s$weibull_shape)) else NA)
}

pars <- c("alpha_control", "beta_cure_arm", "beta_cure_intercept",
          "beta_surv_arm", "beta_surv_intercept")
