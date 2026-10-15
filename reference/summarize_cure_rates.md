# Summarize Cure Rate Results

Calculates the posterior estimates of the cure rates for each treatment
arm and the absolute difference between them.

## Usage

``` r
summarize_cure_rates(model)
```

## Arguments

- model:

  An object of class 'bcm_fit' returned by \`fit_bayesian_cure_model\`.

## Value

A list containing the estimated cure rates for each arm (named
dynamically using the factor levels), their absolute difference, and an
\`interpretation\` string.
