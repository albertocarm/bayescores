
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




#' Plot Utility Score Density Distribution
#'
#' @description
#' Creates a smooth, filled density plot visualizing the posterior distribution of final utility scores.
#'
#' This function uses the `bde` package for bounded density estimation, which is appropriate for scores
#' that are constrained to a specific range (e.g., 0 to 100). The resulting plot is styled with a
#' vertical color gradient fill. It also includes a crucial fix to prevent the density line
#' from dropping below y=0, a common artifact in kernel density estimation near boundaries.
#'
#' @param final_utility_results A list object, typically the output of another function
#'   like `summarize_final_utility()`. This list **must** contain a numeric vector named
#'   `final_utility_vector` which holds the posterior samples of the utility scores.
#'
#' @return A `ggplot` object representing the density plot. This object can be further
#'   customized using standard `ggplot2` functions (e.g., adding more layers or modifying the theme).
#'
#' @export
#' @importFrom bde bde
#' @importFrom ggplot2 ggplot aes geom_segment geom_line scale_color_viridis_c labs scale_x_continuous scale_y_continuous coord_cartesian theme_classic theme element_blank element_text
#'
#' @examples
#' # \dontrun{
#' # --- Create a sample data structure ---
#' # This simulates the expected input for the function.
#' # In a real case, this would come from another part of your analysis.
#' set.seed(42) # for reproducibility
#' sample_scores <- rbeta(2000, shape1 = 30, shape2 = 15) * 100
#' results <- list(final_utility_vector = sample_scores)
#'
#' # --- Generate the plot ---
#' density_plot <- plot_final_utility_density(results)
#'
#' # --- Display the plot ---
#' print(density_plot)
#' # }
plot_final_utility_density <- function(final_utility_results) {

  # -------------------------------------------------------------------------- #
  # 1. DATA EXTRACTION AND VALIDATION
  # -------------------------------------------------------------------------- #

  # Extract the numeric vector of scores from the provided list.
  # Using as.numeric() ensures we are working with a plain numeric vector.
  final_scores <- as.numeric(final_utility_results$final_utility_vector)

  # Input validation: It's crucial to check that the input data is valid before proceeding.
  # This prevents cryptic errors later in the function.
  if (length(final_scores) == 0 || all(!is.finite(final_scores))) {
    stop("The 'final_utility_vector' is empty or does not contain finite values.")
  }

  # -------------------------------------------------------------------------- #
  # 2. DENSITY ESTIMATION
  # -------------------------------------------------------------------------- #

  # Estimate the density using the bde::bde function, which is designed for
  # data that has known upper and lower bounds.
  density_estimation <- bde::bde(
    data  = final_scores,
    estimator = "vitale",   # The 'vitale' estimator is well-suited for bounded variables.
    lower.limit = 0,        # Define the theoretical lower bound of the utility score.
    upper.limit = 100,      # Define the theoretical upper bound of the utility score.
    dataPointsCache = seq(0, 100, by = 0.1) # Use a fine grid for a smooth density curve.
  )

  # -------------------------------------------------------------------------- #
  # 3. DATA PREPARATION FOR PLOTTING
  # -------------------------------------------------------------------------- #

  # Create a data frame from the density estimation results.
  # ggplot2 works best with data frames.
  df_den <- data.frame(
    x = density_estimation@dataPointsCache * 100, # The x-values (scores)
    y = density_estimation@densityCache         # The y-values (estimated density)
  )

  # --- ARTIFACT CORRECTION ---
  # Kernel density estimators can sometimes compute small negative density values
  # near the boundaries of the data range. This is a mathematical artifact.
  # The following line corrects this by setting any negative density values to zero,
  # preventing the plot from incorrectly dipping below the horizontal axis.
  df_den$y[df_den$y < 0] <- 0
  # --- END OF CORRECTION ---

  # -------------------------------------------------------------------------- #
  # 4. PLOT CREATION WITH GGPLOT2
  # -------------------------------------------------------------------------- #

  # Build the plot layer by layer for full control over its appearance.
  density_plot <- ggplot2::ggplot(df_den, ggplot2::aes(x = x, y = y)) +

    # Layer 1: Vertical segments to create a filled area effect with a gradient.
    # Each segment runs from the y=0 axis up to the density curve (y).
    # The 'color = x' aesthetic maps the score value to a color, creating the gradient.
    ggplot2::geom_segment(
      ggplot2::aes(xend = x, yend = 0, color = x),
      show.legend = FALSE # The color legend is not needed for this visualization.
    ) +

    # Layer 2: A solid line on top to clearly outline the density shape.
    ggplot2::geom_line(color = "black", linewidth = 1) +

    # Color Scale: Apply the 'plasma' color palette from viridis.
    # It's perceptually uniform and good for accessibility.
    ggplot2::scale_color_viridis_c(option = "plasma", direction = 1) +

    # Labels and Titles
    ggplot2::labs(
      x     = "Final Utility Score",
      y     = NULL, # Y-axis label is removed for a cleaner look.
      title = "Posterior Distribution of Final Utility"
    ) +

    # X-Axis Customization: Set breaks every 10 units and remove padding.
    ggplot2::scale_x_continuous(
      breaks = seq(0, 100, by = 10),
      expand = c(0, 0)
    ) +

    # Y-Axis Customization: Remove padding so the plot touches the axis.
    ggplot2::scale_y_continuous(expand = c(0, 0)) +

    # Coordinate System: Enforce strict limits on the x-axis.
    ggplot2::coord_cartesian(xlim = c(0, 100)) +

    # Theming: Start with a classic theme and then customize it.
    ggplot2::theme_classic() +
    ggplot2::theme(
      # Remove all y-axis elements for a minimalist 'ridgeline plot' style.
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.line.y  = ggplot2::element_blank(),
      # Center the plot title and adjust its size.
      plot.title   = ggplot2::element_text(hjust = 0.5, size = 14)
    )

  # Return the completed ggplot object. It can now be printed or modified further.
  return(density_plot)
}
