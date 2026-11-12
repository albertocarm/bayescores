# Compare Bayesian Cure Models and Recommend a Distribution

Ranks two or more models fitted with \`fit_bayesian_cure_model\` to the
same data by their expected out-of-sample predictive accuracy and
reports which distribution for the non-cured patients is preferred.
Models whose predictive accuracy is practically the same as that of the
best one are flagged as practical ties.

## Usage

``` r
compare_models(..., criterion = c("loo", "waic"), tie_elpd = 4, tie_se = 2)

# S3 method for class 'bcm_comparison'
print(x, digits = 1, ...)
```

## Arguments

- ...:

  Two or more objects of class \`bcm_fit\`, for example
  \`compare_models(fit1, fit2, fit3)\`, or a single list of them. Models
  are labelled with the argument or list names and, failing that, with
  their distribution.

- criterion:

  Character string. "loo" (default) uses approximate leave-one-out
  cross-validation with Pareto smoothed importance sampling; "waic" uses
  the widely applicable information criterion.

- tie_elpd:

  Numeric. A model whose ELPD is lower than that of the best model by
  less than this amount (default 4) is a practical tie.

- tie_se:

  Numeric. A model that is not a practical tie is "clearly worse" when
  the ELPD difference is at least this many standard errors (default 2)
  and "probably worse" otherwise.

- x:

  An object of class \`bcm_comparison\`.

- digits:

  Number of decimal places shown in the report.

## Value

An object of class \`bcm_comparison\`, whose print method shows the
report. It is a list with elements:

- `comparison`: a data frame with one row per model, ordered from best
  to worst, with the difference in ELPD with respect to the best model
  and its standard error, the ELPD, the effective number of parameters,
  the information criterion, the stacking weight, for "loo" the largest
  Pareto k value, and the verdict.

- `best_model`: the name of the model with the highest ELPD.

- `practical_ties`: the names of the models that are practical ties with
  the best one.

- `recommendation`: the text of the recommendation.

- `criterion`, `tie_elpd`, `tie_se`: the settings used.

- `n_obs`: the number of observations.

- `estimates`: the list of objects returned by the \`loo\` package for
  each model.

## Details

Fitting is the slow step, so the comparison works on models that have
already been fitted, and possibly saved with \`saveRDS\`: nothing is
refitted.

The pointwise log-likelihood of every model is computed from its
posterior draws: an observed event contributes the density of the
non-cured patients weighted by the probability of not being cured, and a
censored observation contributes the mixture of the cure probability and
the survival function of the non-cured patients. The expected log
pointwise predictive density (ELPD) is then estimated with the \`loo\`
package (Vehtari, Gelman and Gabry, 2017). Stacking weights (Yao et al.,
2018) give the weight of each model in the combination with the best
predictive performance. With \`criterion = "loo"\`, Pareto k values
above 0.7 flag models for which the approximation is unreliable.

Each model receives a verdict:

- "best": the highest ELPD.

- "practical tie": the ELPD difference with the best model is below
  \`tie_elpd\`. Differences of this size are small: the models make very
  similar predictions and choosing one or the other has little
  consequence (Sivula et al., 2020).

- "probably worse": the difference is at least \`tie_elpd\` but smaller
  than \`tie_se\` standard errors, so it is not conclusive.

- "clearly worse": the difference is at least \`tie_elpd\` and at least
  \`tie_se\` standard errors.

Any threshold is a convention. Among practical ties the choice should
rest on other grounds, such as the pre-specified model or parsimony, and
the conclusions should be checked to agree across the tied models.

## References

Vehtari A, Gelman A, Gabry J (2017). Practical Bayesian model evaluation
using leave-one-out cross-validation and WAIC. Statistics and Computing,
27, 1413-1432.

Yao Y, Vehtari A, Simpson D, Gelman A (2018). Using stacking to average
Bayesian predictive distributions. Bayesian Analysis, 13, 917-1007.

Sivula T, Magnusson M, Matamoros AA, Vehtari A (2020). Uncertainty in
Bayesian leave-one-out cross-validation based model comparison.
arXiv:2008.10296.

## Examples

``` r
if (FALSE) { # \dontrun{
sim_data <- simulate_weibull_cure_data(
  n_patients = 300, cure_fraction_ctrl = 0.15, cure_fraction_exp = 0.40,
  max_follow_up = 36, weibull_shape = 1.2, median_survival_ctrl = 12,
  time_ratio_exp = 1.25
)
fit1 <- fit_bayesian_cure_model(sim_data, dist = "weibull")
fit2 <- fit_bayesian_cure_model(sim_data, dist = "lognormal")
fit3 <- fit_bayesian_cure_model(sim_data, dist = "loglogistic")
comparison <- compare_models(fit1, fit2, fit3)
comparison
comparison$comparison
} # }
```
