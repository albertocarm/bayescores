#' @title Plot Utility Score Gauge Visualization
#' @description
#' Generates a gauge visualization from the output of `summarize_final_utility`.
#' This version correctly parses component names, dynamically includes a penalty
#' component if present, and displays a clear summary.
#'
#' @param final_utility_results The list object returned by `summarize_final_utility`.
#' @param trial_name Optional. A character string for the trial name. If NULL,
#'  the function will prompt the user to enter a name.
#'
#' @return A ggplot2 object representing the gauge plot.
#' @export
#'
plot_utility_donut <- function(final_utility_results, trial_name = NULL) {

  # --- 1. Load packages ---
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2")
  if (!requireNamespace("viridis", quietly = TRUE)) stop("Please install viridis")

  # --- 2. Prompt for trial name ---
  if (is.null(trial_name)) {
    trial_name <- readline(prompt = "Enter the trial name for the plot title: ")
  }

  # --- 3. Extract medians and look for penalty ---
  summary_data <- final_utility_results$component_summary

  # Extract the main components
  median_eff_score <- summary_data$Median[summary_data$Component == "Efficacy Score (Combined)"]
  median_tox_adj <- summary_data$Median[summary_data$Component == "Toxicity Contribution (points)"]
  median_qol_adj <- summary_data$Median[summary_data$Component == "QoL Contribution (points)"]
  median_final <- summary_data$Median[summary_data$Component == "FINAL UTILITY SCORE"]

  # Dynamic search for the penalty component
  penalty_row <- summary_data[grep("Penalty", summary_data$Component), ]
  penalty_text <- "" # Initialize as an empty string

  if (nrow(penalty_row) > 0) {
    # If one or more penalty rows are found, take the first one
    penalty_name <- penalty_row$Component[1]
    median_penalty <- penalty_row$Median[1]

    # Shorten the specific penalty name for a better fit in the plot
    penalty_name <- gsub(" \\(for TR<1\\)", "<1", penalty_name)

    # Format the text to be inserted into the plot
    penalty_text <- paste0("\n", penalty_name, ": ", round(median_penalty, 1))
  }

  # Safety checks in case a component is missing
  if (length(median_eff_score) == 0) median_eff_score <- 0
  if (length(median_tox_adj) == 0) median_tox_adj <- 0
  if (length(median_qol_adj) == 0) median_qol_adj <- 0
  if (length(median_final) == 0) median_final <- 0

  # --- 4. Prepare data for the plot ---
  plot_max_value <- 100
  median_final_clamped <- max(0, min(plot_max_value, median_final))
  final_score_color <- viridis::viridis(101)[round(median_final_clamped) + 1]

  # --- 5. Generate the plot with dynamic info text ---
  donut_plot <- ggplot2::ggplot() +
    ggplot2::geom_rect(ggplot2::aes(ymax = plot_max_value, ymin = 0, xmax = 4, xmin = 3), fill = "#F0F0F0") +
    ggplot2::geom_rect(ggplot2::aes(ymax = median_final_clamped, ymin = 0, xmax = 4, xmin = 3), fill = final_score_color) +
    # The central text now includes the 'penalty_text' variable
    ggplot2::annotate("text", x = 0, y = 0,
                      label = paste0(
                        "Efficacy Score: ", round(median_eff_score, 1),
                        penalty_text, # This variable will contain the penalty text or be empty
                        "\nToxicity Adj: ", round(median_tox_adj, 1),
                        "\nQoL Adj: ", round(median_qol_adj, 1),
                        "\n-----------------",
                        "\nFinal Score: ", round(median_final, 1)
                      ),
                      size = 5.0, lineheight = 1.1, hjust = 0.5, vjust = 0.5) +
    ggplot2::coord_polar(theta = "y", start = 0) +
    ggplot2::xlim(c(0, 5)) +
    ggplot2::ylim(c(0, plot_max_value)) +
    ggplot2::labs(title = paste("Clinical Utility Score:", trial_name)) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16, margin = ggplot2::margin(b = 5))
    )

  return(donut_plot)
}
