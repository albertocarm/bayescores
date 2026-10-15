# Reconstruct Individual Patient Data from Multiple Kaplan-Meier Curves

A wrapper function that applies the Guyot reconstruction algorithm to
multiple curves (e.g., different treatment arms) provided in a single
data frame. It splits the data by curve, applies
\`recon_one_curve_guyot\` to each, and then combines the results.

## Usage

``` r
reconstruct_ipd(
  plot_km,
  nrisk_tbl,
  totev = NULL,
  ensure_monotone = TRUE,
  censor_placement = c("even", "uniform"),
  jitter_admin = 1e-06,
  verbose = TRUE
)
```

## Arguments

- plot_km:

  (\`data.frame\`)  
  The digitized KM curve points for one or more curves. Must contain
  'time', 'St', and an optional 'curve' identifier column. If 'curve' is
  missing, all data is assumed to belong to a single curve.

- nrisk_tbl:

  (\`data.frame\`)  
  The number of patients at risk for one or more curves. Must contain
  'time_tick', 'nrisk', and an optional 'curve' identifier column.

- totev:

  (\`integer\`, named \`list\`, or \`NULL\`)  
  Optional total number of events for each curve. If a single integer,
  it's applied to all curves. If a named list, names should match the
  'curve' identifiers.

- ensure_monotone:

  (\`logical\`)  
  If \`TRUE\`, forces the survival curves to be monotonically
  decreasing.

- censor_placement:

  (\`character\`)  
  Method for placing censored observations, "even" or "uniform".

- jitter_admin:

  (\`numeric\`)  
  A small time value to add to the last time point for administrative
  censoring.

- verbose:

  (\`logical\`)  
  If \`TRUE\`, prints a summary of the verification check for each
  curve.

## Value

A list containing the combined results: \`ipd\` (\`data.frame\`),
\`drops\` (\`data.frame\`), \`check\` (\`data.frame\`), and \`compat\`
(\`logical\`).
