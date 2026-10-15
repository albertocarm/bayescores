# Toxicity and quality of life data from the SPOTLIGHT trial

This dataset contains simulated toxicity and quality of life (QoL) data
from the SPOTLIGHT randomized clinical trial (zolbetuximab in CLDN18.2+
gastric/GEJ cancer), including both any‑grade (G1–4) and severe (G3–4)
adverse events. Patient‑reported QoL scores are stored in the `qol`
column. The dataset provides the number of patients per trial arm.
Please refer to the accompanying documentation or vignette for details
on how QoL scores are defined and interpreted.

## Usage

``` r
spotlight_toxicity
```

## Format

A data frame with X rows and Y variables:

- EventName:

  Name of the adverse event.

- SystemOrganClass:

  System organ class for the adverse event.

- Incidence_G1_4_Control:

  Incidence (%) of any‑grade (G1–4) toxicity in the control arm.

- Incidence_G1_4_Experimental:

  Incidence (%) of any‑grade (G1–4) toxicity in the experimental arm.

- Incidence_G3_4_Control:

  Incidence (%) of severe (G3–4) toxicity in the control arm.

- Incidence_G3_4_Experimental:

  Incidence (%) of severe (G3–4) toxicity in the experimental arm.

- qol:

  Numeric QoL score associated with each adverse event.

- N_patients:

  Named numeric vector with the number of patients per trial arm.

## Source

<https://clinicaltrials.gov/ct2/show/NCT03504397>
