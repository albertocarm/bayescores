


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

  # --- Plot 3: Odds Ratio (CORREGIDO) ---
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






#' @method model_diagnostics bcm_fit
#' @importFrom rstan stan_trace
#' @export
model_diagnostics.bcm_fit <- function(x, ...) {
  stan_fit <- x$stan_fit
  params_to_check <- c("beta_cure_intercept", "beta_cure_arm",
                       "beta_surv_intercept", "beta_surv_arm", "alpha")
  trace_plot <- rstan::stan_trace(stan_fit, pars = params_to_check)
  print(trace_plot)
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


