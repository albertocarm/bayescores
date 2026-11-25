# Simulation scenarios

Scripts for the simulation scenarios of the BayeScores manuscript (Scenarios 1
to 12). All datasets are analyzed with the default Weibull mixture cure model
of bayescores (neutral Normal(0, 2.5) priors on the log odds ratio of cure and
on the log time ratio, common shape).

## Files

- `simulation_functions.R`: scenario definitions and data generation.
- `model_fit.R`: model fit with the Stan program of the package.
- `summary_functions.R`: performance metrics, as in Table 1 of the manuscript.
- `run_scenarios.R`: runs the scenarios and writes the results.

## Run

Copy this folder to a working directory and run:

```
Rscript run_scenarios.R 500 results
```

The first argument is the number of replicates per scenario and the second
the output folder. A third argument selects scenarios (names separated by
commas) and a fourth splits the replicates between parallel processes:

```
Rscript run_scenarios.R 500 results sc05_high_censoring 1/2
Rscript run_scenarios.R 500 results sc05_high_censoring 2/2
```

One file is saved per replicate, so an interrupted run can be resumed.
Scenario 5 takes about 4 minutes per replicate and the other scenarios about
5 seconds.

## Design

500 replicates per scenario (seed = replicate number, so the same random
numbers are used in every scenario). Patients are randomized 1:1 and the
Weibull shape is 1.5.

Scenarios 1 to 4 and 6 to 12 have 200 patients, a median survival time of
uncured control patients of 12 months and administrative censoring at 120
months. Scenario 5 has 5,000 patients, a median survival time of uncured
control patients of 30.1 months (Weibull scale exp(3.65)) and administrative
censoring at 84 months, which leaves about 58% of the patients censored.

With benefit, the cure fractions are 0.30 (control) and 0.75 (experimental)
and the time ratio is 2. Under the null, the cure fraction is 0.40 in both
arms and the time ratio is 1.

| Scenario | Name in the code | Design |
|---|---|---|
| 1 | `sc01_benefit` | Benefit for the cure fraction and the survival time |
| 2 | `sc02_null` | Null effects |
| 3 | `sc03_cure_up_time_down` | Cure fractions 0.30 and 0.75, time ratio 0.5 |
| 4 | `sc04_detrimental` | Cure fractions 0.60 and 0.30, time ratio 0.5 |
| 5 | `sc05_high_censoring` | Benefit, with high censoring |
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

The results in the manuscript were obtained with R 4.3.3 and rstan 2.32.5.
The Stan program is the same in bayescores 1.2 and 1.2.1.
