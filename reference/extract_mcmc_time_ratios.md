# Extract MCMC Samples for the Time Ratio

Extracts the full vector of posterior MCMC samples for the Time Ratio
(TR), allowing for custom plotting and analysis.

## Usage

``` r
extract_mcmc_time_ratios(model)
```

## Arguments

- model:

  An object of class 'bcm_fit' returned by \`fit_bayesian_cure_model\`.

## Value

A numeric vector containing the posterior samples of the Time Ratio.
