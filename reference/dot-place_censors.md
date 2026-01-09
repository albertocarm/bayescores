# Helper to Place Censored Observations

Distributes a given number of censored observations within a time
interval \`\[a, b\]\`.

## Usage

``` r
.place_censors(a, b, ci, eps = 1e-06, mode = c("even", "uniform"))
```

## Arguments

- a:

  (\`numeric\`)  
  Start of the interval.

- b:

  (\`numeric\`)  
  End of the interval.

- ci:

  (\`integer\`)  
  Number of censored observations to place.

- eps:

  (\`numeric\`)  
  A small epsilon value for numerical comparisons.

- mode:

  (\`character\`)  
  Method for placing censored observations, either "even" (evenly
  spaced) or "uniform" (from a uniform distribution).

## Value

A numeric vector of censored times.
