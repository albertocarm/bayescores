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
  cure_belief = "unknown",
  shared_shape = FALSE,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  seed = 555,
  adapt_delta = 0.99,
  use_historical_prior = FALSE,
  historical_prior_params = c(0, 1),
  ...
)
```

## Arguments

- data:

  A data frame with time, event, and arm columns.

- time_col, event_col, arm_col:

  Character strings for column names.

- cure_belief:

  Character string. Sets the prior belief for the adjuvant cure effect.
  One of "unknown", "unlikely", "very_unlikely", "optimistic" (Heavy
  Radial), or "mild_optimistic" (Standard Radial).

- shared_shape:

  Logical. If TRUE, both treatment arms share the same Weibull shape
  parameter (proportional hazards). If FALSE (default), allows
  independent shapes for each arm.

- chains, iter, warmup, seed:

  Numeric arguments passed to \`rstan::stan\`.

- adapt_delta:

  Target acceptance rate for Stan's NUTS algorithm.

- use_historical_prior:

  Logical. If TRUE, overrides \`cure_belief\` and uses an informative
  historical prior defined by \`historical_prior_params\`. Default is
  FALSE.

- historical_prior_params:

  Numeric vector of length 2 (Mean, SD). Used only if
  \`use_historical_prior = TRUE\`. Defines the Normal prior for the
  treatment effect log-OR. Default is c(0, 1).

- ...:

  Additional arguments passed to \`rstan::stan\`.

## Value

A custom S3 object of class \`bcm_fit\`, a list with elements: -
\`stan_fit\`: the original \`stanfit\` object - \`original_data\`: the
input data frame - \`column_map\`: list mapping \`time_col\`,
\`event_col\`, \`arm_col\` - \`posterior_draws\`: list of posterior
samples for each parameter - \`n_draws\`: the total number of
post-warmup posterior draws
