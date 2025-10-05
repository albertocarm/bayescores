


# --- S3 Methods for Created Objects ---
#' @title Plot Kaplan-Meier Survival Curves
#' @description Creates a Kaplan-Meier survival plot from a data frame, ensuring
#' the treatment arm is handled correctly as a factor. This function is robust
#' enough to handle data from simulations or external sources.
#'
#' @param data A data frame containing survival data.
#' @param time_col The name of the column with time-to-event data.
#' @param event_col The name of the column with the event indicator (1=event, 0=censored).
#' @param arm_col The name of the column with the treatment arm.
#'
#' @return A ggplot object representing the Kaplan-Meier curves.
#' @import ggplot2
#' @importFrom survival survfit Surv
#' @importFrom broom tidy
#' @export
#'
#' @examples
#' \dontrun{
#'   # Example with simulated data
#'   sim_data <- simulate_weibull_cure_data(n_patients = 100)
#'   plot_km_curves(sim_data)
#'
#'   # Example with a package dataset (e.g., MONALEESA2_1)
#'   # data(MONALEESA2_1)
#'   # plot_km_curves(MONALEESA2_1)
#' }
plot_km_curves <- function(data, time_col = "time", event_col = "event", arm_col = "arm") {

  # --- Robustness Check: Ensure 'arm' is a factor ---
  if (!is.factor(data[[arm_col]])) {
    message(paste0("Note: Converting column '", arm_col, "' to a factor for plotting."))
    data[[arm_col]] <- as.factor(data[[arm_col]])
  }

  # Create the survival formula dynamically
  formula_obj <- as.formula(paste("Surv(", time_col, ", ", event_col, ") ~ ", arm_col))

  # Fit the Kaplan-Meier model
  km_fit <- survival::survfit(formula_obj, data = data)

  # Tidy the results for ggplot
  km_df <- broom::tidy(km_fit)

  # Clean up strata names for the legend
  km_df$strata <- gsub(paste0(arm_col, "="), "", km_df$strata)

  # Create the plot
  final_plot <- ggplot(km_df, aes(x = time, y = estimate, color = strata)) +
    geom_step(linewidth = 1) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = strata), alpha = 0.2, linetype = 0) +
    scale_y_continuous(limits = c(0, 1), name = "Survival Probability") +
    scale_x_continuous(name = "Time") +
    labs(
      title = "Kaplan-Meier Survival Curves by Treatment Arm",
      color = "Arm",
      fill = "Arm"
    ) +
    theme_minimal(base_size = 14)

  return(final_plot)
}



#' @method plot bcm_fit
#' @importFrom rstan extract
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_step geom_line scale_y_continuous labs theme_minimal
#' @importFrom survival survfit Surv
#' @importFrom broom tidy
#' @importFrom stats as.formula
#' @export
plot.bcm_fit <- function(x, ...) {
  plot_data <- x$original_data
  cols <- x$column_map
  time_col <- cols$time_col
  event_col <- cols$event_col
  arm_col <- cols$arm_col
  posterior_samples <- rstan::extract(x$stan_fit)
  time_grid <- seq(0, max(plot_data[[time_col]]), length.out = 150)
  predicted_curves <- list()
  arm_levels <- levels(plot_data[[arm_col]])
  for (i in seq_along(arm_levels)) {
    b <- i - 1
    cure_logit <- posterior_samples$beta_cure_intercept + posterior_samples$beta_cure_arm * b
    cure_prob <- 1 / (1 + exp(-cure_logit))
    log_scale <- posterior_samples$beta_surv_intercept + posterior_samples$beta_surv_arm * b
    scale <- exp(log_scale)
    shape <- posterior_samples$alpha
    surv_matrix <- sapply(time_grid, function(t) {
      s_weibull <- exp(-(t / scale)^shape)
      s_total <- cure_prob + (1 - cure_prob) * s_weibull
      return(s_total)
    })
    predicted_curves[[arm_levels[i]]] <- data.frame(
      time = time_grid,
      survival = apply(surv_matrix, 2, mean),
      arm = arm_levels[i]
    )
  }
  predicted_df <- do.call(rbind, predicted_curves)
  formula_obj <- as.formula(paste("survival::Surv(", time_col, ", ", event_col, ") ~ ", arm_col))
  km_fit <- survival::survfit(formula_obj, data = plot_data)
  km_df <- broom::tidy(km_fit)
  km_df$strata <- gsub(paste0(arm_col, "="), "", km_df$strata)
  final_plot <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = km_df,
                         ggplot2::aes(x = time, ymin = conf.low, ymax = conf.high, group = strata),
                         alpha = 0.2, fill = "black") +
    ggplot2::geom_step(data = km_df,
                       ggplot2::aes(x = time, y = estimate, linetype = strata),
                       color = "black") +
    ggplot2::geom_line(data = predicted_df,
                       ggplot2::aes(x = time, y = survival, linetype = arm),
                       color = "red", linewidth = 1) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      title = "Kaplan-Meier vs. Bayesian Model Prediction",
      subtitle = "Black lines: Kaplan-Meier fit. Red lines: Bayesian model prediction.",
      x = "Time",
      y = "Survival Probability",
      linetype = "Arm"
    ) +
    ggplot2::theme_minimal(base_size = 14)
  print(final_plot)
}









#' Plot Posterior Densities
#'
#' A generic function to visualize posterior densities of key parameters
#' from a Bayesian model fit.
#'
#' @param x The fitted model object.
#' @param ... Additional arguments passed to specific methods.
#'
#' @return A ggplot object containing the density plots.
#' @export
plot_densities <- function(x, ...) {
  UseMethod("plot_densities")
}

#' @method plot_densities bcm_fit
#' @rdname plot_densities
#' @export
#' @importFrom rstan extract
#' @importFrom ggplot2 ggplot aes geom_segment geom_line geom_vline scale_color_gradient2 labs theme theme_minimal
#' @importFrom patchwork plot_layout wrap_plots
#' @importFrom stats density quantile
plot_densities.bcm_fit <- function(x, ...) {
  posterior_samples <- rstan::extract(x$stan_fit)

  # --- Plot 1: Time Ratio ---
  time_ratio_draws <- exp(posterior_samples$beta_surv_arm)
  density_surv <- density(time_ratio_draws, n = 4500)
  df_surv <- data.frame(x = density_surv$x, y = density_surv$y)

  plot_surv <- ggplot2::ggplot(df_surv, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_segment(ggplot2::aes(xend = x, yend = 0, colour = x)) +
    ggplot2::geom_line(color = "black") +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "darkgrey") +
    ggplot2::scale_color_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 1) +
    ggplot2::labs(
      title = "Post. Dens. of Time Ratio",
      subtitle = "Effect on Survival Time",
      x = "Time Ratio (TR)",
      y = "Density"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  # --- Plot 2: Cure Rate Difference ---
  prob_cure_ctrl <- 1 / (1 + exp(-posterior_samples$beta_cure_intercept))
  prob_cure_exp <- 1 / (1 + exp(-(posterior_samples$beta_cure_intercept + posterior_samples$beta_cure_arm)))
  cure_rate_diff <- prob_cure_exp - prob_cure_ctrl
  density_diff <- density(cure_rate_diff, n = 2048)
  df_diff <- data.frame(x = density_diff$x, y = density_diff$y)

  plot_diff <- ggplot2::ggplot(df_diff, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_segment(ggplot2::aes(xend = x, yend = 0, colour = x)) +
    ggplot2::geom_line(color = "black") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
    ggplot2::scale_color_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0) +
    ggplot2::labs(
      title = "Cure Rate Diff.",
      subtitle = "P(Cure|Exp) - P(Cure|Ctrl)",
      x = "Difference in Cure Probability",
      y = "Density"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  # --- Plot 3: Odds Ratio ---
  odds_ratio_draws <- exp(posterior_samples$beta_cure_arm)
  upper_limit <- quantile(odds_ratio_draws, probs = 0.95)

  # La corrección está aquí:
  density_cure <- density(odds_ratio_draws, n = 6000, from = 0, to = upper_limit)

  df_cure <- data.frame(x = density_cure$x, y = density_cure$y)

  plot_cure <- ggplot2::ggplot(df_cure, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_segment(ggplot2::aes(xend = x, yend = 0, colour = x)) +
    ggplot2::geom_line(color = "black") +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "darkgrey") +
    ggplot2::scale_color_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 1) +
    ggplot2::labs(
      title = "Post. Dens. of Cure (OR)",
      subtitle = "Effect on Odds of Being Cured",
      x = "Odds Ratio (OR)",
      y = "Density"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  # --- Combinar y mostrar Gráficos ---
  final_plot <- patchwork::wrap_plots(plot_surv, plot_cure, plot_diff)

  # Imprimir el objeto ggplot combinado
  print(final_plot)
}








#' @title Calculate the Time Ratio
#' @description Extracts the Time Ratio (or Acceleration Factor) for the survival
#' component of the cure model and its 95% credible interval.
#'
#' @param model An object of class 'bcm_fit' returned by
#'   `fit_bayesian_cure_model`.
#'
#' @return A named list containing the estimated `time_ratio` (posterior median),
#'   its 95% credible interval `ci_95`, and an `interpretation` string.
#' @importFrom rstan extract
#' @importFrom stats median quantile
#' @export
time_ratio <- function(model) {
  # Check if the input is the correct class
  if (!inherits(model, "bcm_fit")) {
    stop("Input must be an object of class 'bcm_fit'.")
  }

  # Extract posterior samples from the stanfit object inside the model list
  posterior <- rstan::extract(model$stan_fit)

  # Get the samples for the survival effect parameter
  beta_surv_samples <- posterior$beta_surv_arm

  # Calculate the Time Ratio by exponentiating the samples
  tr_samples <- exp(beta_surv_samples)

  # Calculate the posterior median as the point estimate
  estimate <- stats::median(tr_samples)

  # Calculate the 95% credible interval
  ci <- stats::quantile(tr_samples, probs = c(0.025, 0.975))

  # --- NEW: Add automatic interpretation ---
  if (ci[1] > 1) {
    interpretation <- "Significant evidence that the treatment extends the time to event for non-cured patients (Credible Interval is above 1)."
  } else if (ci[2] < 1) {
    interpretation <- "Significant evidence that the treatment shortens the time to event for non-cured patients (Credible Interval is below 1)."
  } else {
    interpretation <- "No significant evidence of an effect on the time to event for non-cured patients (Credible Interval includes 1)."
  }
  # ----------------------------------------

  # Return the results in a clean list
  result <- list(
    time_ratio = estimate,
    ci_95 = ci,
    interpretation = interpretation # Add interpretation to the output
  )

  return(result)
}


#' @title Plot Correlated Posterior Densities
#' @description
#' Creates a 2D density plot for two correlated parameters from a Stan model fit.
#' Axis limits are computed independently for X and Y using quantiles for a tight
#' data-driven window, with optional padding. Ellipses mark specified probability
#' levels under a bivariate normal approximation.
#'
#' @details
#' This helper expects a list-like object \code{x} that contains a fitted Stan
#' object under \code{x$stan_fit}. The Stan fit must expose posterior draws for
#' two parameters named \code{beta_cure_arm} and \code{beta_surv_arm}. These are
#' interpreted as:
#' \itemize{
#'   \item \code{beta_cure_arm}: log(OR) for the cure component.
#'   \item \code{beta_surv_arm}: log(TR) for the survival component among the uncured.
#' }
#' The function:
#' \enumerate{
#'   \item extracts posterior samples with \code{rstan::extract()},
#'   \item computes Pearson correlation,
#'   \item sets independent axis limits via \code{stats::quantile()},
#'   \item estimates a 2D KDE via \code{MASS::kde2d()},
#'   \item builds a \code{ggplot2} heatmap with contours and confidence ellipses.
#' }
#'
#' @param x A list-like object that contains a \code{stanfit} at \code{x$stan_fit}.
#'   The fit must have parameters \code{beta_cure_arm} and \code{beta_surv_arm}.
#' @param n_grid Integer. Number of grid points per axis for the 2D kernel density
#'   estimation passed to \code{MASS::kde2d}. Default is \code{100}.
#' @param level_ellipses Numeric vector of probability levels for confidence
#'   ellipses (bivariate normal contours). Default \code{c(0.5, 0.8, 0.95)}.
#' @param quantile_range Numeric vector of length 2 with lower and upper tail
#'   probabilities used to set axis limits independently. Defaults to the central
#'   99.8\%: \code{c(0.001, 0.999)}.
#' @param padding Numeric scalar (>= 0). Fractional padding to expand each axis
#'   range beyond the selected quantiles. Default \code{0.05} (5\%).
#'
#' @return A \code{ggplot} object representing the joint posterior density.
#'
#' @importFrom rstan extract
#' @importFrom MASS kde2d
#' @importFrom dplyr mutate
#' @importFrom tidyr expand_grid
#' @importFrom stats quantile cor
#' @importFrom rlang .data
#' @import ggplot2
#' @export
#'
#' @examples
#' \dontrun{
#' # Suppose 'fit' is a stanfit with parameters beta_cure_arm and beta_surv_arm:
#' obj <- list(stan_fit = fit)
#' p <- plot_correlated_densities(obj,
#'                                n_grid = 150,
#'                                level_ellipses = c(0.5, 0.9),
#'                                quantile_range = c(0.005, 0.995),
#'                                padding = 0.1)
#' print(p)
#' }
plot_correlated_densities <- function(x, n_grid = 100,
                                      level_ellipses = c(0.5, 0.8, 0.95),
                                      quantile_range = c(0.001, 0.999),
                                      padding = 0.05) {

  # --- 1) Extract posterior draws ---
  post <- rstan::extract(x$stan_fit)
  posterior_draws <- data.frame(
    log_or = as.numeric(post$beta_cure_arm),
    log_hr = as.numeric(post$beta_surv_arm)
  )

  # --- 2) Calculate correlation ---
  rho <- cor(posterior_draws$log_or, posterior_draws$log_hr)

  # --- 3) Define axis limits INDEPENDENTLY for X and Y ---
  xlims <- stats::quantile(posterior_draws$log_or, probs = quantile_range, na.rm = TRUE)
  ylims <- stats::quantile(posterior_draws$log_hr, probs = quantile_range, na.rm = TRUE)

  # Add padding to each axis
  xlims <- xlims + c(-1, 1) * diff(xlims) * padding
  ylims <- ylims + c(-1, 1) * diff(ylims) * padding

  # --- 4) 2D density estimation ---
  kde_lims <- c(xlims, ylims)
  kd <- MASS::kde2d(posterior_draws$log_or, posterior_draws$log_hr, n = n_grid, lims = kde_lims)
  density_df <- tidyr::expand_grid(x = kd$x, y = kd$y) |>
    dplyr::mutate(density = as.vector(t(kd$z)))

  # --- 5) Build the plot ---
  p <- ggplot(density_df, aes(x = x, y = y)) +
    geom_tile(aes(fill = density)) +
    geom_contour(aes(z = density), colour = "white", alpha = 0.6, bins = 10) +
    scale_fill_viridis_c(option = "plasma") +
    annotate("point", x = mean(posterior_draws$log_or), y = mean(posterior_draws$log_hr),
             colour = "red", size = 3) +
    coord_cartesian(xlim = xlims, ylim = ylims, expand = FALSE) +
    labs(
      title = "Joint Posterior Density",
      subtitle = paste0("Pearson correlation: ", round(rho, 3)),
      x = "log(OR) Cure", y = "log(TR) Survival", fill = "Density"
    ) +
    theme_minimal(base_size = 14)

  return(p)
}
