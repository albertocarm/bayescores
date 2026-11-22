# Additional simulation scenarios

Scripts for the additional simulation scenarios of the BayeScores manuscript
(Scenarios 6 to 12). All datasets are analyzed with the default Weibull
mixture cure model of bayescores 1.2 (neutral Normal(0, 2.5) priors on the log
odds ratio of cure and on the log time ratio, common shape).

## Files

- `simulation_functions.R`: scenario definitions and data generation.
- `model_fit.R`: model fit with the Stan program of the package.
- `summary_functions.R`: performance metrics, as in Table 1 of the manuscript.
- `run_additional_scenarios.R`: runs the seven scenarios and writes the results.

## Run

Copy this folder to a working directory and run:

```
Rscript run_additional_scenarios.R 500 results
```

The first argument is the number of replicates per scenario and the second
the output folder.

## Design

200 patients randomized 1:1 and 500 replicates per scenario (seed = replicate
number). With benefit, the cure fractions are 0.30 (control) and 0.75
(experimental) and the time ratio is 2. Under the null, the cure fraction is
0.40 in both arms and the time ratio is 1. The median survival time of uncured
control patients is 12 months, with administrative censoring at 120 months.

| Scenario | Name in the code | Departure from the Weibull model |
|---|---|---|
| 6 | `sc06_lognormal_benefit` | Log-normal survival times, benefit |
| 7 | `sc07_lognormal_null` | Log-normal survival times, null effects |
| 8 | `sc08_loglogistic_benefit` | Log-logistic survival times, benefit |
| 9 | `sc09_loglogistic_null` | Log-logistic survival times, null effects |
| 10 | `sc10_delayed_effect` | No effect on event times during the first 3 months |
| 11 | `sc11_random_dropout` | Random loss to follow-up, 5% per year |
| 12 | `sc12_staggered_accrual` | Uniform accrual over 24 months, data cut-off at 120 months |

Log-normal and log-logistic survival times share the median and the standard
deviation of the log time of the Weibull distribution with shape 1.5.

MCMC: 2 chains of 2,000 iterations (1,000 warm-up), `adapt_delta` 0.99.

The results in the manuscript were obtained with R 4.3.3, rstan 2.32.5 and
bayescores 1.2.
