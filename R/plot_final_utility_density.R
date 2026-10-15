
#' Print method for bayes_summary objects
#'
#' Displays a formatted summary of the disaggregated and final scores.
#'
#' @param x An object of class 'bayes_summary'.
#' @param ... further arguments passed to or from other methods.
#' @export
print.bayes_summary <- function(x, ...) {

  # Calculate cure weight for printing
  cure_weight <- 1 - x$tr_weight

  cat("--- Bayesian Utility Score Summary ---\n\n")

  # --- Component Scores Section ---
  cat("## Component Scores\n")
  cat(sprintf(
    "%-12s (Weight: %2.0f%%): Median = %5.2f, 95%% CrI = [%5.2f, %5.2f]\n",
    "Time Ratio",
    x$tr_weight * 100,
    x$tr_summary$median,
    x$tr_summary$ci_95[1],
    x$tr_summary$ci_95[2]
  ))
  cat(sprintf(
    "%-12s (Weight: %2.0f%%): Median = %5.2f, 95%% CrI = [%5.2f, %5.2f]\n\n",
    "Cure Rate",
    cure_weight * 100,
    x$cure_summary$median,
    x$cure_summary$ci_95[1],
    x$cure_summary$ci_95[2]
  ))

  # --- Final Score Section ---
  cat("## Final Composite Score\n")
  cat(sprintf("Median Score : %5.2f\n", x$final_summary$median))
  cat(sprintf("95%% CrI      : [%5.2f, %5.2f]\n", x$final_summary$ci_95[1], x$final_summary$ci_95[2]))
  cat(sprintf("5-Level Score: %5.2f\n", x$final_summary$rescaled_score))

  cat("\n------------------------------------\n")
  invisible(x)
}



#' Plot the Density of Final Utility Scores
#'
#' @description
#' Creates a gradient-filled density plot to visualize the posterior
#' distribution of final utility scores. This function uses a two-step process:
#' first, it calculates the density from the data, and then it uses ggplot2
#' to create the visualization.
#'
#' @param final_utility_results A list containing the model's output. This list
#'   must contain an element named \code{final_utility_vector}, which is a
#'   numeric vector of the posterior utility scores.
#'
#' @return A \code{ggplot} object representing the density plot.
#'
#' @importFrom stats density
#' @importFrom ggplot2 ggplot aes geom_segment geom_line scale_color_viridis_c
#' @importFrom ggplot2 labs scale_y_continuous coord_cartesian scale_x_continuous
#' @importFrom ggplot2 theme_classic theme element_blank element_text
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # 1. Create a sample data structure
#' # (simulating 1000 scores from a Beta distribution scaled to 0-100)
#' sample_scores <- rbeta(1000, shape1 = 30, shape2 = 10) * 100
#' sample_results <- list(final_utility_vector = sample_scores)
#'
#' # 2. Generate the plot
#' plot_final_utility_density(sample_results)
#' }
plot_final_utility_density <- function(final_utility_results) {

  # 1. DATA EXTRACTION AND VALIDATION
  final_scores <- as.numeric(final_utility_results$final_utility_vector)
  if (length(final_scores) == 0 || all(!is.finite(final_scores))) {
    stop("The 'final_utility_vector' is empty or does not contain finite values.")
  }

  # 2. DENSITY ESTIMATION (MANUAL)
  # Use the density() function from R base to calculate coordinates.
  # Specify 'from' and 'to' to cover the full 0-100 range.
  density_estimation <- stats::density(
    final_scores,
    from = 0,
    to = 100,
    n = 1000 # Controls the smoothness of the curve
  )

  # 3. DATA PREPARATION FOR PLOTTING
  # Create a data frame with the results, ready for ggplot2.
  df_den <- data.frame(
    x = density_estimation$x,
    y = density_estimation$y
  )
  # Artifact correction: Set any small negative density values to zero.
  df_den$y[df_den$y < 0] <- 0


  # 4. PLOT CREATION WITH GGPLOT2
  density_plot <- ggplot2::ggplot(df_den, ggplot2::aes(x = x, y = y)) +

    # Layer 1: Vertical segments to create the gradient fill.
    ggplot2::geom_segment(
      ggplot2::aes(xend = x, yend = 0, color = x)
    ) +

    # Layer 2: A solid line on top to outline the shape.
    ggplot2::geom_line(color = "black", linewidth = 1) +

    # Color Scale
    ggplot2::scale_color_viridis_c(option = "plasma") +

    # Labels and Titles
    ggplot2::labs(
      x     = "Final Utility Score",
      y     = NULL,
      title = "Posterior Distribution of Final Utility"
    ) +

    # Axes and Theming
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::coord_cartesian(xlim = c(0, 100)) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, 100, by = 10),
      expand = c(0, 0)
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.line.y  = ggplot2::element_blank(),
      plot.title   = ggplot2::element_text(hjust = 0.5, size = 14),
      legend.position = "none"
    )

  return(density_plot)
}
