# Summarizes the final utility score with advanced dynamic weighting.

This function calculates the final weighted utility score. It
incorporates a sophisticated weighting policy that adapts based on the
detected efficacy profile ('Cure' vs. 'Survival').

In the 'Survival' profile, the weight of the 'cure' component is scaled
by its observed magnitude, preventing a weak signal from diluting the
score.

In all cases, the weight of the 'QoL' component is also scaled by its
magnitude, preventing a neutral QoL from diluting the score.

Any "unused" weight from these gradual adjustments is redistributed
proportionally among the other active components.

## Usage

``` r
summarize_final_utility(
  efficacy_scores,
  toxicity_scores,
  qol_scores,
  cure_benefit_threshold = 0.02,
  prob_certainty_threshold = 0.95,
  tr_min_if_cure = 1,
  tr_min_if_no_cure = 1.15,
  tr_hopeful = 1.25,
  cure_min_relevant = 0.04,
  cure_hopeful = 0.12,
  tr_calibration_k = 0.35,
  weights_if_cure = c(tr = 0.15, cure = 0.6, tox = 0.15, qol = 0.1),
  weights_if_no_cure = c(tr = 0.7, cure = 0.1, tox = 0.1, qol = 0.1)
)
```

## Arguments

- efficacy_scores:

  A list with 'tr_posterior_samples' and 'cure_posterior_samples'.

- toxicity_scores:

  A numeric vector of toxicity adjustment scores (-1 to +1).

- qol_scores:

  A numeric vector of sampled QoL scores (-2 to +2).

- cure_benefit_threshold:

  Numeric. The threshold above which the cure fraction is considered
  beneficial. Default is 0.02.

- prob_certainty_threshold:

  Numeric. The probability threshold (e.g., 0.95) to determine if a
  benefit is "certain".

- tr_min_if_cure:

  Numeric. The minimum Time Ratio (TR) considered beneficial if a cure
  is present.

- tr_min_if_no_cure:

  Numeric. The minimum Time Ratio (TR) considered beneficial if no cure
  is present.

- tr_hopeful:

  Numeric. The Time Ratio threshold for a "hopeful" or highly desirable
  outcome.

- cure_min_relevant:

  Numeric. The minimum cure fraction considered clinically relevant.

- cure_hopeful:

  Numeric. The cure fraction threshold for a "hopeful" or highly
  desirable outcome.

- tr_calibration_k:

  Numeric. The calibration exponent (k) applied to the Time Ratio (TR)
  utility scores. Default is 0.35.

- weights_if_cure:

  A named numeric vector of weights for \`c(tr, cure, tox, qol)\` if a
  "Cure" profile is detected.

- weights_if_no_cure:

  A named numeric vector of weights for \`c(tr, cure, tox, qol)\` if a
  "Survival" (no cure) profile is detected.

## Value

A list containing the final utility vector (with floor at 0) and a
component summary.
