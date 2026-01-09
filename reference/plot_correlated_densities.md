# Plot Correlated Posterior Densities

Creates a 2D density plot for two correlated parameters from a Stan
model fit. Axis limits are computed independently for X and Y using
quantiles for a tight data-driven window, with optional padding.
Ellipses mark specified probability levels under a bivariate normal
approximation.

## Usage

``` r
plot_correlated_densities(
  x,
  n_grid = 100,
  level_ellipses = c(0.5, 0.8, 0.95),
  quantile_range = c(0.001, 0.999),
  padding = 0.05,
  correlation_method = c("pearson", "spearman", "kendall")
)
```

## Arguments

- x:

  A list-like object that contains a `stanfit` at `x$stan_fit`. The fit
  must have parameters `beta_cure_arm` and `beta_surv_arm`.

- n_grid:

  Integer. Number of grid points per axis for the 2D kernel density
  estimation passed to
  [`MASS::kde2d`](https://rdrr.io/pkg/MASS/man/kde2d.html). Default is
  `100`.

- level_ellipses:

  Numeric vector of probability levels for confidence ellipses
  (bivariate normal contours). Default `c(0.5, 0.8, 0.95)`.

- quantile_range:

  Numeric vector of length 2 with lower and upper tail probabilities
  used to set axis limits independently. Defaults to the central 99.8%:
  `c(0.001, 0.999)`.

- padding:

  Numeric scalar (\>= 0). Fractional padding to expand each axis range
  beyond the selected quantiles. Default `0.05` (5%).

- correlation_method:

  Character. The method to use: 'pearson' (default), 'spearman', or
  'kendall'. Passed to [`stats::cor`](https://rdrr.io/r/stats/cor.html).

## Value

A `ggplot` object representing the joint posterior density.

## Details

This helper expects a list-like object `x` that contains a fitted Stan
object under `x$stan_fit`. The Stan fit must expose posterior draws for
two parameters named `beta_cure_arm` and `beta_surv_arm`. These are
interpreted as:

- `beta_cure_arm`: log(OR) for the cure component.

- `beta_surv_arm`: log(TR) for the survival component among the uncured.

The function:

1.  extracts posterior samples with
    [`rstan::extract()`](https://mc-stan.org/rstan/reference/stanfit-method-extract.html),

2.  computes Pearson, Spearman, or Kendall correlation,

3.  sets independent axis limits via
    [`stats::quantile()`](https://rdrr.io/r/stats/quantile.html),

4.  estimates a 2D KDE via
    [`MASS::kde2d()`](https://rdrr.io/pkg/MASS/man/kde2d.html),

5.  builds a `ggplot2` heatmap with contours and confidence ellipses.

## Examples

``` r
if (FALSE) { # \dontrun{
# Suppose 'fit' is a stanfit with parameters beta_cure_arm and beta_surv_arm:
obj <- list(stan_fit = fit)
p <- plot_correlated_densities(obj,
                               n_grid = 150,
                               correlation_method = "spearman",
                               quantile_range = c(0.005, 0.995))
print(p)
} # }
```
