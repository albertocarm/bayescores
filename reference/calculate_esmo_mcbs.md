# Calculate the ESMO-MCBS v2.0 Score (Guided by Objective)

Calculates the ESMO-MCBS v2.0 score based on a primary study objective.
The function first requires the user to specify the study's main goal,
and then uses the corresponding arguments to apply the correct ESMO form
logic.

## Usage

``` r
calculate_esmo_mcbs(objective, ...)
```

## Arguments

- objective:

  A character string specifying the main goal of the clinical trial.
  This is a \*\*mandatory\*\* argument.

- ...:

  Additional arguments required for the chosen objective. See the
  documentation and examples for which arguments are needed for each
  objective.

## Value

A character string with the final ESMO-MCBS score.
