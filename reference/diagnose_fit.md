# Diagnose Bayesian MCMC fit

Extracts effective sample size (n_eff) and Rhat diagnostics from a
fitted Bayesian model (e.g., Stan fit) and provides categorical
interpretations along with practical recommendations for improvement.

## Usage

``` r
diagnose_fit(
  fit,
  n_eff_threshold_high = 1000,
  n_eff_threshold_ok = 500,
  rhat_excellent = 1.01,
  rhat_ok = 1.05,
  include_lp = TRUE
)
```

## Arguments

- fit:

  A fitted Bayesian model object (e.g., stanfit).

- n_eff_threshold_high:

  Numeric. Threshold above which n_eff is considered "high" (default =
  1000).

- n_eff_threshold_ok:

  Numeric. Threshold above which n_eff is considered "moderate" (default
  = 500).

- rhat_excellent:

  Numeric. Threshold below which Rhat is considered "excellent" (default
  = 1.01).

- rhat_ok:

  Numeric. Threshold below which Rhat is considered "acceptable"
  (default = 1.05).

- include_lp:

  Logical. If FALSE, excludes lp from the output table (default = TRUE).

## Value

A data frame with columns

- n_eff effective sample size

- Rhat Gelman-Rubin convergence diagnostic

- Rhat_interpretation categorical interpretation of Rhat

- n_eff_interpretation categorical interpretation of n_eff

- recommendation practical suggestions if diagnostics are poor

## Examples

``` r
if (FALSE) { # \dontrun{
diagnose_fit(bayesian_fit_kn199$stan_fit)
diagnose_fit(bayesian_fit_kn199$stan_fit, include_lp = FALSE)
} # }
```
