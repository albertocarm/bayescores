# Generate a Summary Table of Model Outcomes (Final Version 2.1)

Extracts data directly from the fit object to calculate and display key
model outcomes, now including the Shape Ratio (SR) to account for
flexible hazard dynamics.

## Usage

``` r
outcomes(
  fit,
  calibration_args = list(),
  digits = 2,
  shrinkage_method = "none",
  shrinkage_target = "primary",
  correlation_method = c("pearson", "spearman", "kendall")
)
```

## Arguments

- fit:

  The fitted model object from \`fit_bayesian_cure_model\`.

- calibration_args:

  An optional list to override default utility calibration.

- digits:

  The number of decimal places for rounding.

- shrinkage_method:

  A character string specifying the shrinkage method. Options are "none"
  (default), "zwet", "sherry", or "skeptical".

- shrinkage_target:

  A character string controlling the application of shrinkage. Options
  are "primary" or "both_if_uncorrelated".

- correlation_method:

  Character. The method to use: 'pearson' (default), 'spearman', or
  'kendall'. Used for the identifiability check.

## Value

A tibble with a summary of key model outcomes.
