# Reconstruct IPD for a Single Curve Using Guyot's Algorithm

This function takes digitized Kaplan-Meier (KM) curve data and a table
of the number of patients at risk at specific time points to reconstruct
individual patient data (IPD) for a single treatment arm or group. It
implements the algorithm described by Guyot et al. (2012).

## Usage

``` r
recon_one_curve_guyot(
  d_curve,
  nr_curve,
  totev = NULL,
  eps = 1e-06,
  ensure_monotone = TRUE,
  censor_placement = c("even", "uniform"),
  jitter_admin = 1e-06,
  check = TRUE
)
```

## Arguments

- d_curve:

  (\`data.frame\`)  
  The digitized KM curve points. Must contain 'time' and 'St' (survival
  probability) columns.

- nr_curve:

  (\`data.frame\`)  
  The number of patients at risk. Must contain 'time_tick' and 'nrisk'
  columns.

- totev:

  (\`integer\` or \`NULL\`)  
  Optional total number of events. If provided, the algorithm will
  adjust the event count in the last interval to match this total.

- eps:

  (\`numeric\`)  
  A small epsilon value for numerical comparisons.

- ensure_monotone:

  (\`logical\`)  
  If \`TRUE\`, forces the survival curve to be monotonically decreasing.

- censor_placement:

  (\`character\`)  
  Method for placing censored observations, "even" or "uniform".

- jitter_admin:

  (\`numeric\`)  
  A small time value to add to the last time point for administrative
  censoring.

- check:

  (\`logical\`)  
  If \`TRUE\`, performs a check to verify if the reconstructed number at
  risk matches the published numbers.

## Value

A list containing: \`ipd\` (the reconstructed IPD \`data.frame\`),
\`drops\` (details of survival probability drops), \`check\` (the
verification \`data.frame\`), and \`compat\` (a logical indicating if
the check was successful).

## References

Guyot, P., Ades, A. E., Ouwens, M. J., & Welton, N. J. (2012). Enhanced
secondary analysis of survival data: reconstructing the data from
published Kaplan-Meier survival curves. \*BMC Medical Research
Methodology\*, 12(1), 1-13.
