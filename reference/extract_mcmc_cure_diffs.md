# Extract MCMC Samples for the Difference in Cure Rates

Extracts the full vector of posterior MCMC samples for the absolute
difference in cure rates between the treatment and control arms.

## Usage

``` r
extract_mcmc_cure_diffs(model)
```

## Arguments

- model:

  An object of class 'bcm_fit' returned by \`fit_bayesian_cure_model\`.

## Value

A numeric vector containing the posterior samples of the difference in
cure rates (treatment - control).
