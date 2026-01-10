# Check for False Plateau or Structural Instability

Analyzes the estimated long-term survival fraction (plateau) to
distinguish between:

- Absence of plateau (below threshold).

- Robust plateau (high estimate, 0

- False/Unstable plateau (visual plateau but high probability of being
  an artifact).

## Usage

``` r
check_false_plateau(fit, threshold = 0.05)
```

## Arguments

- fit:

  An object containing the fitted Bayesian model (expected to have a
  `stan_fit` slot).

- threshold:

  A numeric value defining the boundary for clinical relevance. Defaults
  to 0.05 (5%).

## Value

Prints a context-aware diagnostic summary and invisibly returns a list
containing the calculated probabilities and mean rates for both arms.
