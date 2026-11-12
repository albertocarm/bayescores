# Fit the Bayesian Cure Model (Engine)

Prepares data, calls the external Stan cure model, and stores MCMC
draws.

## Usage

``` r
fit_bayesian_cure_model(
  data,
  time_col = "time",
  event_col = "event",
  arm_col = "arm",
  tail_assumption = "neutral",
  shared_shape = FALSE,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  seed = 555,
  adapt_delta = 0.99,
  use_historical_prior = FALSE,
  historical_prior_params = c(0, 1),
  dist = "weibull",
  ...
)
```

## Arguments

- data:

  A data frame with time, event, and arm columns.

- time_col, event_col, arm_col:

  Character strings for column names.

- tail_assumption:

  Character string. Sets the prior belief strategy for the adjuvant cure
  effect. One of "neutral", "immature_skeptical", or
  "biologically_null".

- shared_shape:

  Logical. If TRUE, both treatment arms share the same Weibull shape
  parameter (proportional hazards). If FALSE (default), allows
  independent shapes for each arm.

- chains, iter, warmup, seed:

  Numeric arguments passed to \`rstan::stan\`.

- adapt_delta:

  Target acceptance rate for Stan's NUTS algorithm.

- use_historical_prior:

  Logical. If TRUE, overrides \`tail_assumption\` and uses an
  informative historical prior defined by \`historical_prior_params\`.
  Default is FALSE.

- historical_prior_params:

  Numeric vector of length 2 (Mean, SD). Used only if
  \`use_historical_prior = TRUE\`. Defines the Normal prior for the
  treatment effect log-OR. Default is c(0, 1).

- dist:

  Character string. Parametric distribution assumed for the event times
  of the non-cured patients. One of "weibull" (default), "exponential",
  "lognormal", "loglogistic", "gamma", or "gengamma" (generalized
  gamma). All of them are accelerated failure time models, so the time
  ratio keeps its interpretation. The Weibull distribution is the
  pre-specified default; the alternatives are intended for sensitivity
  analyses and can optionally be compared with \`compare_models\`.

- ...:

  Additional arguments passed to \`rstan::stan\`.

## Value

A custom S3 object of class \`bcm_fit\`, a list with elements: -
\`stan_fit\`: the original \`stanfit\` object - \`original_data\`: the
input data frame - \`column_map\`: list mapping \`time_col\`,
\`event_col\`, \`arm_col\` - \`posterior_draws\`: list of posterior
samples for each parameter - \`n_draws\`: the total number of
post-warmup posterior draws - \`dist\`: the distribution used for the
non-cured patients

## Details

For every distribution the location of the non-cured patients on the
log-time scale is \`beta_surv_intercept + beta_surv_arm \* arm\`, and
\`exp(beta_surv_arm)\` is the time ratio. The shape parameter
\`alpha_control\` (multiplied by \`exp(beta_alpha_arm)\` in the
experimental arm when \`shared_shape = FALSE\`) has the following
meaning:

- "weibull": Weibull shape; \`exp(location)\` is the Weibull scale.

- "exponential": no shape parameter (fixed at 1); \`exp(location)\` is
  the mean.

- "lognormal": inverse of the standard deviation of the log times;
  \`exp(location)\` is the median.

- "loglogistic": log-logistic shape; \`exp(location)\` is the median.

- "gamma": gamma shape; \`exp(location)\` is the mean.

- "gengamma": inverse of the scale parameter \`sigma\` of the
  generalized gamma distribution in the parameterization of Prentice
  (1974), with \`location\` as \`mu\` and the additional shape parameter
  \`q_shape\`, common to both arms, as \`Q\` (\`Q = 1\` is the Weibull
  and \`Q = 0\` the log-normal distribution).

The priors of the parameters shared with the Weibull model are the same
for every distribution; \`q_shape\` has a standard normal prior. With
"exponential" the argument \`shared_shape\` has no effect.

Priors of the treatment effects. With \`tail_assumption = "neutral"\`
(the default) the log odds ratio of cure (\`beta_cure_arm\`) and the log
time ratio (\`beta_surv_arm\`) have independent Normal(0, 2.5) priors. A
standard deviation of 2.5 is the scale of the weakly informative default
of Gelman et al. (2008): a shift of 5 on the logit scale moves a
probability from 0.01 to 0.5, so larger effects are implausible but not
excluded. The prior places 95% of its mass on odds ratios and time
ratios between 0.007 and 134 and is nearly flat over the range of
effects seen in clinical trials, so it stabilizes the computation
without pulling the estimates towards the null. With
"immature_skeptical" and "biologically_null" the log odds ratio has a
Laplace(0, 0.2) and a Laplace(0, 0.1) prior, respectively, and the log
time ratio keeps the Normal(0, 2.5) prior. The intercepts have
Student-t(4, 0, 2.5) priors and the shape parameter a Gamma(1, 1) prior.

## References

Gelman A, Jakulin A, Pittau MG, Su YS (2008). A weakly informative
default prior distribution for logistic and other regression models. The
Annals of Applied Statistics, 2, 1360-1383.

Prentice RL (1974). A log gamma model and its maximum likelihood
estimation. Biometrika, 61, 539-544.
