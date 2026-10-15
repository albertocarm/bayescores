# Create an AMIT plot (Adverse Event Dot Plot) from a clinical trial object.

This version is flexible and allows specifying the names of the elements
within the trial object.

## Usage

``` r
create_amit_plot(
  trial_object,
  grade_type = "any_grade",
  main_title = "",
  data_element = "datos",
  n_element = "N_pacientes",
  ...
)
```

## Arguments

- trial_object:

  The list object containing the trial data.

- grade_type:

  A string: "any_grade" (G1–4) or "severe_grade" (G3–4).

- main_title:

  The main title of the plot.

- data_element:

  The name (as a character string) of the element containing the
  toxicity data frame.

- n_element:

  The name (as a character string) of the element containing the vector
  with patient counts.

- ...:

  Additional arguments passed to the HH::AEdotplot function.
