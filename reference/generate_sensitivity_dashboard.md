# Generate Comprehensive Sensitivity Analysis Dashboard

Orchestrates all simulations and generates the final 8-panel
publication-ready plot.

## Usage

``` r
generate_sensitivity_dashboard(calibration_args)
```

## Arguments

- calibration_args:

  A list containing the final calibration parameters for the utility
  model (e.g., efficacy targets).

## Value

A composite ggplot object created with the patchwork package.
