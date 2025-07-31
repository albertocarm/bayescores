#' @title Plot Utility Score Gauge Visualization
#' @description
#' Generates a gauge visualization from the output of `summarize_final_utility`.
#' This version correctly parses component names and displays a clear summary.
#'
#' @param final_utility_results The list object returned by `summarize_final_utility`.
#' @param trial_name Optional. A character string for the trial name. If NULL,
#'  the function will prompt the user to enter a name.
#'
#' @return A ggplot2 object representing the gauge plot.
#'
plot_utility_donut <- function(final_utility_results, trial_name = NULL) {

  # --- 1. Cargar paquetes ---
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2")
  if (!requireNamespace("viridis", quietly = TRUE)) stop("Please install viridis")

  # --- 2. Pedir nombre del ensayo ---
  if (is.null(trial_name)) {
    trial_name <- readline(prompt = "Enter the trial name for the plot title: ")
  }

  # --- 3. Extraer medianas usando los nombres EXACTOS del input ---
  summary_data <- final_utility_results$component_summary

  # Se extraen los 3 componentes principales y el score final
  median_eff_score <- summary_data$Median[summary_data$Component == "Efficacy Score (Combined)"]
  median_tox_adj <- summary_data$Median[summary_data$Component == "Toxicity Contribution (points)"]
  median_qol_adj <- summary_data$Median[summary_data$Component == "QoL Contribution (points)"]
  median_final <- summary_data$Median[summary_data$Component == "FINAL UTILITY SCORE"]

  # Comprobaciones de seguridad por si algún componente falta
  if (length(median_eff_score) == 0) median_eff_score <- 0
  if (length(median_tox_adj) == 0) median_tox_adj <- 0
  if (length(median_qol_adj) == 0) median_qol_adj <- 0
  if (length(median_final) == 0) median_final <- 0

  # --- 4. Preparar datos para el gráfico ---
  plot_max_value <- 100
  median_final_clamped <- max(0, min(plot_max_value, median_final))
  final_score_color <- viridis::viridis(101)[round(median_final_clamped) + 1]

  # --- 5. Generar el gráfico con el texto informativo actualizado ---
  donut_plot <- ggplot2::ggplot() +
    ggplot2::geom_rect(ggplot2::aes(ymax = plot_max_value, ymin = 0, xmax = 4, xmin = 3), fill = "#F0F0F0") +
    ggplot2::geom_rect(ggplot2::aes(ymax = median_final_clamped, ymin = 0, xmax = 4, xmin = 3), fill = final_score_color) +
    # Texto central actualizado para ser más informativo
    ggplot2::annotate("text", x = 0, y = 0,
                      label = paste0(
                        "Efficacy Score: ", round(median_eff_score, 1),
                        "\nToxicity Adj: ", round(median_tox_adj, 1),
                        "\nQoL Adj: ", round(median_qol_adj, 1),
                        "\n-----------------",
                        "\nFinal Score: ", round(median_final, 1)
                      ),
                      size = 5.5, lineheight = 1.1, hjust = 0.5, vjust = 0.5) +
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
