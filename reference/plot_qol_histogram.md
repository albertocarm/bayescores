# Plot an Elegant Histogram of QoL Outcomes

This function takes a numeric vector of sampled Quality of Life (QoL)
scores and creates an elegant bar plot that functions as a histogram. It
visualizes the distribution of outcomes, showing the probability of each
category on the y-axis and using descriptive labels on the x-axis.

## Usage

``` r
plot_qol_histogram(qol_scores_vector)
```

## Arguments

- qol_scores_vector:

  A numeric vector of sampled QoL scores, where each element is one of
  \`\[-2, -1, 0, 1, 2\]\`. This is typically the output from
  \`sample_qol_scores()\`.

## Value

A \`ggplot\` object representing the customized histogram. This object
can be printed to display the plot or be further modified.

## Details

The function first calculates the proportion (probability) of each
unique score \`\[-2, -1, 0, 1, 2\]\`. It then maps these numeric scores
to meaningful labels (e.g., "Significant Deterioration") and uses a
diverging color palette to intuitively represent negative, neutral, and
positive outcomes. The bars are directly labeled with their
corresponding probabilities.

## Examples

``` r
# --- 1. First, generate a vector of sampled scores ---
# prob_vector <- c(0.05, 0.20, 0.50, 0.20, 0.05)
# scores <- sample_qol_scores(prob_vector, n_samples = 2000)

# --- 2. Now, plot the results with the new function ---
# final_plot <- plot_qol_histogram(scores)

# --- 3. Display the plot ---
# print(final_plot)
```
