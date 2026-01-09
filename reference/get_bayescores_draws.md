# Get Original or Shrunk MCMC Draws for Subsequent Analyses (Updated)

Extracts or generates posterior draws for efficacy inputs. This version
incorporates advanced shrinkage methods, including a skeptical prior via
whitening, and strategically applies them based on the primary utility
endpoint and parameter correlation.

## Usage

``` r
get_bayescores_draws(
  fit,
  shrinkage_method = "none",
  shrinkage_target = "primary",
  calibration_args = list()
)
```

## Arguments

- fit:

  The fitted model object from \`fit_bayesian_cure_model\`.

- shrinkage_method:

  A character string specifying the shrinkage method. Options are
  "zwet", "sherry", "skeptical", or "none" (default).

- shrinkage_target:

  A character string controlling shrinkage application. Options are
  "primary" or "both_if_uncorrelated".

- calibration_args:

  An optional list to override default utility calibration, used to
  determine the primary endpoint.

## Value

A list containing two named vectors of posterior samples:

- `tr_posterior_samples`: Draws for the Time Ratio (already
  exponentiated).

- `cure_posterior_samples`: Draws for the absolute difference in cure
  rates.
