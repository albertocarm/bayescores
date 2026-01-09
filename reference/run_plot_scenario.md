# Run a Scenario for a Sensitivity Plot

Generates the data needed for a single sensitivity plot line chart.

## Usage

``` r
run_plot_scenario(
  x_var_name,
  x_range,
  line_var_name,
  line_levels,
  fixed_params,
  cal_args
)
```

## Arguments

- x_var_name:

  The parameter to vary on the x-axis.

- x_range:

  The range for the x-axis variable.

- line_var_name:

  The parameter to vary as different lines.

- line_levels:

  The levels for the line variable.

- fixed_params:

  A list of fixed parameters.

- cal_args:

  A list of calibration arguments.

## Value

A data frame with simulation results.
