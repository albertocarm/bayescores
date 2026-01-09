# Calculate a complete toxicity analysis

Calculates weighted toxicity scores (WTS).

## Usage

``` r
calculate_toxicity_analysis(
  trial_data,
  n_simulations,
  unacceptable_rel_increase = 0.5,
  k_uncertainty = 5,
  soc_weights = c(`Gastrointestinal disorders` = 1.2,
    `Blood and lymphatic system disorders` = 1.6,
    `General, metabolic, and other disorders` = 1.3, `Dermatologic disorders` = 1.1,
    `Infections and infestations` = 1.6,
    `Respiratory, thoracic and mediastinal disorders` = 2.5)
)
```

## Arguments

- trial_data:

  A list containing the toxicity data.

- n_simulations:

  An integer for the number of simulations.

- unacceptable_rel_increase:

  A number for the unacceptable increase limit.

- k_uncertainty:

  A number for the uncertainty constant.

- soc_weights:

  The default weights for toxicity categories

## Value

A list containing the results.
