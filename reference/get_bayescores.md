# Calculate BayeScores Final Utility Score (Version 4, Corrected Toxicity Model)

This version corrects the toxicity calculation model as per user
feedback. Toxicity's impact is now proportional to the base efficacy
score.

## Usage

``` r
get_bayescores(
  efficacy_inputs,
  qol_scores,
  toxicity_scores,
  calibration_args = list(),
  correlation_method = "pearson"
)
```

## Arguments

- efficacy_inputs:

  A list with \`cure_posterior_samples\`, \`tr_posterior_samples\`.

- qol_scores:

  A numeric vector of posterior samples for QoL effect.

- toxicity_scores:

  A numeric vector of posterior samples for Toxicity effect.

- calibration_args:

  Optional list to override default calibration.

- correlation_method:

  Character. The method to use: 'pearson' (default), 'spearman', or
  'kendall'.

## Value

A list with the final utility vector and a summary data frame.
