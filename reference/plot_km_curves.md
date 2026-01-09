# Plot Kaplan-Meier Survival Curves

Creates a Kaplan-Meier survival plot from a data frame, ensuring the
treatment arm is handled correctly as a factor. This function is robust
enough to handle data from simulations or external sources.

## Usage

``` r
plot_km_curves(data, time_col = "time", event_col = "event", arm_col = "arm")
```

## Arguments

- data:

  A data frame containing survival data.

- time_col:

  The name of the column with time-to-event data.

- event_col:

  The name of the column with the event indicator (1=event, 0=censored).

- arm_col:

  The name of the column with the treatment arm.

## Value

A ggplot object representing the Kaplan-Meier curves.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Example with simulated data
  sim_data <- simulate_weibull_cure_data(n_patients = 100)
  plot_km_curves(sim_data)

  # Example with a package dataset (e.g., MONALEESA2_1)
  # data(MONALEESA2_1)
  # plot_km_curves(MONALEESA2_1)
} # }
```
