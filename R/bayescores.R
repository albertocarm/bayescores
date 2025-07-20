#' Calculate Bayesian scores for clinical benefit (BSCB)
#'
#' Transforms posterior MCMC samples of the Time Ratio (TR) and the cure rate
#' difference into utility score distributions (0–100) and a weighted composite
#' score.
#'
#' Scaling uses the Normal CDF:
#' \deqn{\text{TR\_score} = 100\,\Phi\!\left(\frac{t - \mu_{TR}}{\sigma_{TR}}\right), \qquad
#' \sigma_{TR} = \frac{\text{tr\_max}-\text{tr\_mean}}{z_{0.99}},\ z_{0.99}=\Phi^{-1}(0.99)}
#' and analogously for the cure rate difference. Thus `tr_mean` (or `cure_mean`)
#' maps to ≈50 and `tr_max` (or `cure_max`) to ≈99.
#'
#' @param tr_posterior_samples Numeric vector of MCMC samples for the Time Ratio.
#' @param cure_posterior_samples Numeric vector of MCMC samples for the cure rate difference.
#' @param tr_mean TR value mapping to utility 50. Default 1.25.
#' @param tr_max TR value mapping to utility ≈99. Default 1.50.
#' @param cure_mean Cure rate difference mapping to utility 50. Default 0.065.
#' @param cure_max Cure rate difference mapping to utility ≈99. Default 0.13.
#' @param tr_weight Weight (in [0,1]) of the TR score in the composite. Default 0.8.
#'
#' @return A list with components:
#' \describe{
#'   \item{tr_scores}{Numeric vector of TR utility scores (0–100).}
#'   \item{cure_scores}{Numeric vector of cure rate difference utility scores (0–100).}
#'   \item{final_scores}{Numeric vector of weighted composite scores (0–100).}
#' }
#'
#' @examples
#' \dontrun{
#' # Assume `bayesian_fit` is an object returned by your model fitting function,
#' # and the extractor helpers return numeric vectors of posterior samples:
#'
#' example_tr_mcmc   <- extract_mcmc_time_ratios(bayesian_fit)
#' example_cure_mcmc <- extract_mcmc_cure_diffs(bayesian_fit)
#'
#' scores <- calculate_bayescores(example_tr_mcmc, example_cure_mcmc)
#' summary(scores$final_scores)
#' }
#'
#' # (Optional) Self‑contained illustration (remove if you prefer only the real-workflow example):
#' # set.seed(123)
#' # tr_sim   <- rnorm(2000, mean = 1.30, sd = 0.08)
#' # cure_sim <- rnorm(2000, mean = 0.07,  sd = 0.02)
#' # simulate_scores <- calculate_bayescores(tr_sim, cure_sim)
#' # summary(simulate_scores$final_scores)
#'
#' @importFrom stats qnorm pnorm
#' @export
calculate_bayescores <- function(
    tr_posterior_samples,
    cure_posterior_samples,
    tr_mean = 1.25,
    tr_max = 1.50,
    cure_mean = 0.065,
    cure_max = 0.13,
    tr_weight = 0.8
) {
  # 1. Utility function parameters
  z_99   <- stats::qnorm(0.99)
  tr_sd   <- (tr_max   - tr_mean)   / z_99
  cure_sd <- (cure_max - cure_mean) / z_99

  # 2. Sub‑scores (0–100)
  tr_scores   <- 100 * stats::pnorm(tr_posterior_samples,   mean = tr_mean,   sd = tr_sd)
  cure_scores <- 100 * stats::pnorm(cure_posterior_samples, mean = cure_mean, sd = cure_sd)

  # 3. Weighted composite
  final_scores <- tr_weight * tr_scores + (1 - tr_weight) * cure_scores

  # 4. Output
  list(
    tr_scores    = tr_scores,
    cure_scores  = cure_scores,
    final_scores = final_scores
  )
}


#' Summarize BayeScores Results
#'
#' Calculates summary statistics for final BayeScores. This function returns an
#' object of class 'bayes_summary'. The summary is printed to the console
#' when the object is called.
#'
#' @param scores_list A list object returned by 'calculate_bayescores'.
#'
#' @return An object of class 'bayes_summary' containing the median, 95% CrI,
#'   and the rescaled score.
#' @export
#' @importFrom stats median quantile
summarize_bayescores <- function(scores_list) {
  # Extraer scores finales
  final_scores <- scores_list$final_scores

  # Calcular estadísticas
  summary_data <- list(
    median = stats::median(final_scores),
    ci_95 = stats::quantile(final_scores, probs = c(0.025, 0.975)),
    rescaled_score = (stats::median(final_scores) / 25) + 1
  )

  # Asignar la clase especial para usar el método print personalizado
  class(summary_data) <- c("bayes_summary", "list")

  return(summary_data)
}

#' Print method for bayes_summary objects
#' @param x An object of class 'bayes_summary'.
#' @param ... further arguments passed to or from other methods.
#' @export
#' @importFrom stats median quantile
print.bayes_summary <- function(x, ...) {
  cat("--- BayeScore Final Summary ---\n")
  cat(sprintf("Median Score : %.2f\n", x$median))
  cat(sprintf("95%% CrI      : [%.2f, %.2f]\n", x$ci_95[1], x$ci_95[2]))
  cat(sprintf("5-Level Score: %.2f\n", x$rescaled_score))
  cat("-----------------------------\n")
  invisible(x)
}


#' Plot BayeScores Density Distribution
#'
#' Creates a smooth, filled density plot of the posterior distribution of the final scores.
#'
#' @param scores_list A list object returned by 'calculate_bayescores'.
#'
#' @return A ggplot object representing the density plot.
#' @export
#' @importFrom bde bde
#' @importFrom ggplot2 ggplot aes geom_col geom_line scale_fill_viridis_c labs scale_x_continuous scale_y_continuous coord_cartesian theme_classic theme element_blank element_text
plot_bayes_dist <- function(scores_list) {

  final_scores <- as.numeric(scores_list$final_scores)

  density_estimation <-
    bde::bde(final_scores, estimator="vitale", lower.limit = 0, upper.limit = 100, dataPointsCache = seq(0, 100, 0.1/2))

  df_den <- data.frame(
    x = density_estimation@dataPointsCache*100,
    y = density_estimation@densityCache
  )

  dx <- diff(df_den$x)[1]

  density_plot <- ggplot2::ggplot(df_den, ggplot2::aes(x = x, y = y, fill = x)) +
    ggplot2::geom_col(
      width        = dx,
      alpha        = 0.9,
      show.legend = FALSE
    ) +
    ggplot2::geom_line(
      color     = "black",
      linewidth = 1
    ) +
    ggplot2::scale_fill_viridis_c(option = "plasma", direction = 1, guide = "none") +
    ggplot2::labs(
      x     = "BayeScores",
      y     = NULL,
      title = "Posterior distribution of BayeScores"
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, 100, by = 10),
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 100)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.y    = ggplot2::element_blank(),
      axis.ticks.y   = ggplot2::element_blank(),
      axis.title.y   = ggplot2::element_blank(),
      axis.line.y    = ggplot2::element_blank(),
      plot.title     = ggplot2::element_text(hjust = 0.5, size = 14)
    )

  return(density_plot)
}

#' Plot BayeScore Donut Visualization
#'
#' Generates a donut-shaped visualization of the BayeScore.
#'
#' @param scores_list A list object returned by `calculate_bayescores`.
#' @param tr_weight The weight assigned to the Time Ratio component in the final score. Default is 0.8.
#'
#' @return A ggplot2 object representing the donut plot.
#'
#' @importFrom dplyr %>% mutate lag
#' @importFrom stats setNames median
#' @importFrom viridis viridis
#' @importFrom ggplot2 ggplot annotate geom_rect aes scale_fill_manual coord_polar xlim ylim theme_void guides
#' @export
plot_bayes_donut <- function(scores_list, tr_weight = 0.8) {
  median_final_score <- stats::median(scores_list$final_scores)
  transformed_score <- (median_final_score / 25) + 1
  score_color <- viridis::viridis(101, option = "plasma")[round(median_final_score) + 1]

  total_area <- median_final_score / 100
  mean_weighted_tr <- mean(scores_list$tr_scores * tr_weight)
  mean_weighted_cure <- mean(scores_list$cure_scores * (1 - tr_weight))

  total_mean_contribution <- mean_weighted_tr + mean_weighted_cure

  if (total_mean_contribution > 0) {
    proportion_tr <- mean_weighted_tr / total_mean_contribution
    proportion_cure <- mean_weighted_cure / total_mean_contribution
  } else {
    proportion_tr <- 1
    proportion_cure <- 0
  }

  area_tr <- total_area * proportion_tr
  area_cure <- total_area * proportion_cure

  plot_data_gauge <- data.frame(component = c("Cure", "TR"), value = c(area_cure, area_tr)) %>%
    dplyr::mutate(
      ymax_pos = cumsum(value),
      ymin_pos = dplyr::lag(ymax_pos, default = 0),
      color_hex = c("#729fcf", score_color)
    )

  gauge_plot <- ggplot2::ggplot(plot_data_gauge) +
    ggplot2::annotate("rect", ymax = 1, ymin = 0, xmax = 4, xmin = 3, fill = "#D3D3D3", color = NA) +
    ggplot2::geom_rect(ggplot2::aes(ymax = ymax_pos, ymin = ymin_pos, xmax = 4, xmin = 3, fill = component)) +
    ggplot2::scale_fill_manual(values = stats::setNames(plot_data_gauge$color_hex, plot_data_gauge$component)) +
    ggplot2::annotate("text", x = 0.3, y = 0.3, label = paste0(round(median_final_score, 1), "%"), size = 7, fontface = "bold") +
    ggplot2::annotate("text", x = 0.3, y = 0.7, label = paste("5-level:", round(transformed_score, 1)), size = 7, fontface = "bold") +
    ggplot2::coord_polar(theta = "y", start = -pi / 2, clip = "off") +
    ggplot2::xlim(c(0, 4)) + ggplot2::ylim(c(0, 1)) +
    ggplot2::theme_void() +
    ggplot2::guides(fill = "none")

  return(gauge_plot)
}
