# Simulate Survival Data from a Cure Model

Generates a dataset for survival analysis based on a two-arm (Control
vs. Experimental) clinical trial scenario. It uses a mixture cure model,
where a fraction of subjects are considered "cured". For the
"susceptible" subjects, their time-to-event follows a Weibull
distribution.

## Usage

``` r
simulate_weibull_cure_data(
  n_patients,
  cure_fraction_ctrl,
  cure_fraction_exp,
  max_follow_up,
  weibull_shape,
  median_survival_ctrl,
  time_ratio_exp,
  seed = 42
)
```

## Arguments

- n_patients:

  The total number of patients to simulate.

- cure_fraction_ctrl:

  The proportion (0-1) of "cured" patients in the Control arm.

- cure_fraction_exp:

  The proportion (0-1) of "cured" patients in the Experimental arm.

- max_follow_up:

  The maximum patient follow-up time (administrative censoring).

- weibull_shape:

  The shape parameter of the Weibull distribution.

- median_survival_ctrl:

  The median survival time for susceptible patients in the Control arm.

- time_ratio_exp:

  The ratio of the experimental arm's median survival time to the
  control arm's.

- seed:

  An optional integer to set the random seed for reproducibility.

## Value

A data.frame with the class \`survival_sim_data\`.

## Examples

``` r
if (FALSE) { # \dontrun{
  sim_data <- simulate_weibull_cure_data(n_patients = 100,
                                         cure_fraction_ctrl = 0.3,
                                         cure_fraction_exp = 0.6,
                                         max_follow_up = 60,
                                         weibull_shape = 1.5,
                                         median_survival_ctrl = 15,
                                         time_ratio_exp = 2.0)
  head(sim_data)
} # }
```
