# Interactively or non-interactively simulate a realistic clinical trial dataset

This function generates a flattened list object containing simulated
data for toxicity, patient numbers, and Quality of Life (QoL).

## Usage

``` r
simulate_trial_data(
  n_control = NULL,
  ratio_str = NULL,
  control_g1_4_pct = NULL,
  control_g3_4_pct = NULL,
  tox_ratio = NULL,
  qol_scenario = NULL,
  qol_strength = NULL
)
```

## Arguments

- n_control:

  Number of patients in the control arm. If NULL, prompts interactively.

- ratio_str:

  Character string indicating the randomization ratio (e.g., "2:1"). If
  NULL, prompts interactively.

- control_g1_4_pct:

  Overall percentage of any-grade toxicity (grades 1–4) in the control
  arm. If NULL, prompts interactively.

- control_g3_4_pct:

  Overall percentage of severe (grade 3–4) toxicity in the control arm.
  If NULL, prompts interactively.

- tox_ratio:

  Relative toxicity ratio (experimental vs control arm). If NULL,
  prompts interactively.

- qol_scenario:

  Integer from 1 to 5 indicating the predefined QoL scenario. If NULL,
  prompts interactively.

- qol_strength:

  Integer from 1 to 5 indicating the strength of evidence for the QoL
  effect. If NULL, prompts interactively.

## Value

A list with three main components:

- toxicity:

  A data frame with simulated adverse event data, including grades and
  patient IDs.

- N_patients:

  A named numeric vector with patient counts in the control and
  experimental arms.

- qol:

  A numeric vector of length 5 representing the final probability
  distribution for QoL outcomes.
