# Plot Utility Score Gauge Visualization

Generates a gauge visualization from the output of
\`summarize_final_utility\`. This version correctly parses component
names, dynamically includes a penalty component if present, and displays
a clear summary.

## Usage

``` r
plot_utility_donut(final_utility_results, trial_name = NULL)
```

## Arguments

- final_utility_results:

  The list object returned by \`summarize_final_utility\`.

- trial_name:

  Optional. A character string for the trial name. If NULL, the function
  will prompt the user to enter a name.

## Value

A ggplot2 object representing the gauge plot.
