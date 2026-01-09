# Calculate the ASCO Net Health Benefit (NHB) Score

This function calculates the Net Health Benefit (NHB) based on the
American Society of Clinical Oncology (ASCO) Value Framework. It
differentiates between two clinical settings: "advanced" (metastatic)
and "adjuvant" (curative intent), as their scoring criteria differ
significantly.

## Usage

``` r
calculate_asco_nhb(
  setting,
  os_hr = NA,
  os_median_control = NA,
  os_median_treatment = NA,
  pfs_hr = NA,
  pfs_median_control = NA,
  pfs_median_treatment = NA,
  orr_treatment = NA,
  dfs_hr = NA,
  total_tox_points_control = 0,
  total_tox_points_treatment = 0,
  bonus_tail_of_curve = FALSE,
  bonus_qol = FALSE,
  bonus_palliation = FALSE,
  bonus_treatment_free_interval_points = 0
)
```

## Arguments

- setting:

  A string specifying the clinical context. Must be either "advanced" or
  "adjuvant".

- os_hr:

  Numeric. The Hazard Ratio for Overall Survival.

- os_median_control:

  Numeric. Median Overall Survival for the control arm (in months).

- os_median_treatment:

  Numeric. Median Overall Survival for the treatment arm (in months).

- pfs_hr:

  Numeric. The Hazard Ratio for Progression-Free Survival.

- pfs_median_control:

  Numeric. Median Progression-Free Survival for the control arm (in
  months).

- pfs_median_treatment:

  Numeric. Median Progression-Free Survival for the treatment arm (in
  months).

- orr_treatment:

  Numeric. The Objective Response Rate for the treatment arm (as a
  decimal, e.g., 0.40 for 40%).

- dfs_hr:

  Numeric. The Hazard Ratio for Disease-Free Survival (used in
  'adjuvant' setting).

- total_tox_points_control:

  A numeric value for the sum of toxicity points for the
  control/standard-of-care arm.

- total_tox_points_treatment:

  A numeric value for the sum of toxicity points for the
  new/experimental arm.

- bonus_tail_of_curve:

  Logical (TRUE/FALSE). Was a significant long-term survival benefit
  ("tail of the curve") demonstrated?

- bonus_qol:

  Logical (TRUE/FALSE). Was a statistically significant improvement in
  Quality of Life reported?

- bonus_palliation:

  Logical (TRUE/FALSE). Was significant symptom palliation shown? (for
  'advanced' setting).

- bonus_treatment_free_interval_points:

  Numeric. Points for treatment-free interval (for 'advanced' setting).

## Value

A named list with four elements: `clinical_benefit_score`,
`toxicity_adjustment`, `bonus_score`, and `net_health_benefit`.

## Examples

``` r
# Example 1: Advanced Disease (Metastatic)
# Hypothetical trial where the primary endpoint is OS improvement shown by HR.
calculate_asco_nhb(
  setting = "advanced",
  os_hr = 0.75,
  total_tox_points_control = 50,
  total_tox_points_treatment = 60,
  bonus_qol = TRUE,
  bonus_palliation = FALSE
)
#> $clinical_benefit_score
#> [1] 25
#> 
#> $toxicity_adjustment
#> [1] -4
#> 
#> $bonus_score
#> [1] 10
#> 
#> $net_health_benefit
#> [1] 31
#> 

# Example 2: Adjuvant Setting (Curative Intent)
# Based on the SPARTAN trial for apalutamide.
calculate_asco_nhb(
  setting = "adjuvant",
  os_hr = 0.78,
  total_tox_points_control = 37.0,
  total_tox_points_treatment = 43.7,
  bonus_tail_of_curve = TRUE,
  bonus_qol = FALSE
)
#> $clinical_benefit_score
#> [1] 22
#> 
#> $toxicity_adjustment
#> [1] -3.62
#> 
#> $bonus_score
#> [1] 20
#> 
#> $net_health_benefit
#> [1] 38.38
#> 
```
