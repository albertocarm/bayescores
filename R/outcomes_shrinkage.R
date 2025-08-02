# -------------------------------------------------------------------
# Internal Helper Functions (not exported to the user)
# These functions are used by both main functions to avoid code duplication.
# -------------------------------------------------------------------

#' @noRd
#' @title Conditional Distribution Calculator
#' @description Computes the parameters of the corrected posterior distribution for the SNR.
#' This code is an implementation of the method from Van Zwet et al. (2021).
#' @param z The observed z-score.
#' @param p A vector of mixture proportions for the prior.
#' @param tau A vector of standard deviations for the SNR prior components.
#' @return A data frame with the parameters (q, m, sigma) of the conditional mixture distribution.
conditional <- function(z, p, tau) {
  tau2 <- tau^2
  q_unscaled <- p * dnorm(z, 0, sqrt(tau2 + 1))
  if (all(q_unscaled == 0)) {
    q <- rep(1/length(p), length(p))
  } else {
    q <- q_unscaled / sum(q_unscaled)
  }
  m <- z * tau2 / (tau2 + 1)
  v <- tau2 / (tau2 + 1)
  sigma <- sqrt(v)
  data.frame(q, m, sigma)
}

#' @noRd
#' @title Shrinkage Sampler
#' @description Generates a new vector of "corrected" MCMC draws based on an empirical prior.
#' @param original_draws A numeric vector of the original MCMC draws (on the log scale).
#' @param method A character string, either "zwet" or "sherry", specifying the prior.
#' @return A numeric vector of the new, corrected MCMC draws.
sample_shrinkage_draws <- function(original_draws, method) {
  b <- mean(original_draws)
  s <- sd(original_draws)
  if (s == 0) return(original_draws) # No variation, no shrinkage
  z_obs <- b / s
  n_samples <- length(original_draws)

  if (method == "zwet") {
    p_prior <- c(0.32, 0.31, 0.30, 0.07)
    tau_prior <- c(0.61, 1.42, 2.16, 5.64)
  } else if (method == "sherry") {
    p_prior <- c(0.0137, 0.2015, 0.7848)
    tau_prior <- c(1.357, 1.368, 3.603)
  } else {
    stop("Internal error: Invalid shrinkage method passed to sampler.")
  }

  corrected_dist_params <- conditional(z = z_obs, p = p_prior, tau = tau_prior)
  q_corr <- corrected_dist_params$q
  m_corr <- corrected_dist_params$m
  sigma_corr <- corrected_dist_params$sigma
  n_components <- length(p_prior)

  corrected_snr_samples <- numeric(n_samples)
  for (i in 1:n_samples) {
    component <- sample(1:n_components, size = 1, prob = q_corr)
    corrected_snr_samples[i] <- rnorm(1, mean = m_corr[component], sd = sigma_corr[component])
  }
  s * corrected_snr_samples
}


# -------------------------------------------------------------------
# Main User-Facing Functions
# -------------------------------------------------------------------

#' Generate a Summary Table of Model Outcomes
#'
#' @description
#' Calculates key outcomes from a Bayesian model fit and can display
#' side-by-side comparisons with empirically shrunk estimates.
#'
#' @param fit The fitted model object from `fit_bayesian_cure_model`.
#' @param digits The number of decimal places for rounding.
#' @param shrinkage_method A character string specifying the shrinkage method for
#'   the corrected estimates. Options are "zwet", "sherry", or "none" (default).
#'
#' @return A tibble with a summary of key model outcomes.
#'
#' @importFrom rstan extract
#' @importFrom tibble tibble as_tibble
#' @importFrom stats quantile median sd rnorm dnorm
#' @importFrom tools toTitleCase
#' @export
outcomes <- function(fit,
                     digits = 2,
                     shrinkage_method = "none") {

  # --- Input Validation ---
  valid_methods <- c("zwet", "sherry", "none")
  if (!shrinkage_method %in% valid_methods) {
    stop(paste("Invalid shrinkage_method specified. Please use one of:",
               paste(valid_methods, collapse = ", ")))
  }

  posterior <- rstan::extract(fit$stan_fit)
  if (!"beta_cure_intercept" %in% names(posterior)) {
    stop("Required parameters not found in the model fit.")
  }

  # --- Prepare Draws for Summary ---
  original_log_tr <- posterior$beta_surv_arm
  original_log_or <- posterior$beta_cure_arm

  original_p_cure_ctrl <- 1 / (1 + exp(-posterior$beta_cure_intercept))
  original_p_cure_exp <- 1 / (1 + exp(-(posterior$beta_cure_intercept + original_log_or)))
  original_cure_diff <- original_p_cure_exp - original_p_cure_ctrl

  # --- Helper to format the output string ---
  summarize_draws <- function(draws, is_log_scale = FALSE, multiplier = 1) {
    summary_draws <- if (is_log_scale) exp(draws) else draws
    point_estimate <- median(summary_draws) * multiplier
    ci <- quantile(summary_draws, probs = c(0.025, 0.975)) * multiplier
    paste0(
      format(round(point_estimate, digits), nsmall = digits),
      " (", format(round(ci[1], digits), nsmall = digits), " - ", format(round(ci[2], digits), nsmall = digits), ")"
    )
  }

  # --- Build the Summary Table ---
  summary_df <- data.frame(
    Metric = c("Time Ratio (TR)", "Odds Ratio (OR) for Cure", "Long-Term Survival Rate (%) - Control",
               "Long-Term Survival Rate (%) - Experimental", "Absolute Difference in Survival Rate (%)"),
    `Result (95% CI)` = c(
      summarize_draws(original_log_tr, is_log_scale = TRUE),
      summarize_draws(original_log_or, is_log_scale = TRUE),
      summarize_draws(original_p_cure_ctrl, multiplier = 100),
      summarize_draws(original_p_cure_exp, multiplier = 100),
      summarize_draws(original_cure_diff, multiplier = 100)
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (shrinkage_method != "none") {
    # Generate corrected draws only for the summary table
    corrected_log_tr <- sample_shrinkage_draws(original_log_tr, method = shrinkage_method)
    corrected_log_or <- sample_shrinkage_draws(original_log_or, method = shrinkage_method)

    # Recalculate derived quantities using the corrected OR draws
    corrected_p_cure_exp <- 1 / (1 + exp(-(posterior$beta_cure_intercept + corrected_log_or)))
    corrected_cure_diff <- corrected_p_cure_exp - original_p_cure_ctrl

    # Update the rows for Experimental Survival and Difference with corrected values
    summary_df[4, "Result (95% CI)"] <- summarize_draws(corrected_p_cure_exp, multiplier = 100)
    summary_df[5, "Result (95% CI)"] <- summarize_draws(corrected_cure_diff, multiplier = 100)

    # Update the metric names to indicate they are corrected
    method_label <- tools::toTitleCase(shrinkage_method)
    summary_df$Metric[4] <- paste0(summary_df$Metric[4], " - Corrected (", method_label, ")")
    summary_df$Metric[5] <- paste0(summary_df$Metric[5], " - Corrected (", method_label, ")")

    # Create and insert the new rows for corrected TR and OR
    corrected_rows <- data.frame(
      Metric = c(paste0("Time Ratio (TR) - Corrected (", method_label, ")"),
                 paste0("Odds Ratio (OR) - Corrected (", method_label, ")")),
      `Result (95% CI)` = c(
        summarize_draws(corrected_log_tr, is_log_scale = TRUE),
        summarize_draws(corrected_log_or, is_log_scale = TRUE)
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    summary_df <- rbind(summary_df[1:2,], corrected_rows, summary_df[3:5,])
  }

  return(tibble::as_tibble(summary_df))
}


#' Get Original or Shrunk MCMC Draws for Subsequent Analyses
#'
#' @description
#' Extracts or generates posterior draws for efficacy inputs, which can be used
#' in other functions like `get_bayescores`. It can return the original MCMC draws
#' or generate new draws corrected with an empirical shrinkage method.
#'
#' @param fit The fitted model object from `fit_bayesian_cure_model`.
#' @param shrinkage_method A character string specifying the shrinkage method.
#'   Options are "zwet", "sherry", or "none" (default).
#'
#' @return A list containing two named vectors of posterior samples:
#'   \itemize{
#'     \item{\code{tr_posterior_samples}: Draws for the Time Ratio (already exponentiated).}
#'     \item{\code{cure_posterior_samples}: Draws for the absolute difference
#'       in cure rates.}
#'   }
#'
#' @importFrom rstan extract
#' @export
get_bayescores_draws <- function(fit, shrinkage_method = "none") {

  # --- Input Validation ---
  valid_methods <- c("zwet", "sherry", "none")
  if (!shrinkage_method %in% valid_methods) {
    stop(paste("Invalid shrinkage_method specified. Please use one of:",
               paste(valid_methods, collapse = ", ")))
  }

  posterior <- rstan::extract(fit$stan_fit)
  if (!"beta_cure_intercept" %in% names(posterior)) {
    stop("Required parameters not found in the model fit.")
  }

  # --- Determine final draws based on method ---
  if (shrinkage_method == "none") {
    # For "none", use the original draws directly from the fit object
    final_log_tr_draws <- posterior$beta_surv_arm
    final_log_or_draws <- posterior$beta_cure_arm
  } else {
    # For "zwet" or "sherry", generate new corrected draws
    final_log_tr_draws <- sample_shrinkage_draws(posterior$beta_surv_arm, method = shrinkage_method)
    final_log_or_draws <- sample_shrinkage_draws(posterior$beta_cure_arm, method = shrinkage_method)
  }

  # --- Calculate derived quantities from the final set of draws ---
  p_cure_ctrl_draws <- 1 / (1 + exp(-posterior$beta_cure_intercept))
  p_cure_exp_draws <- 1 / (1 + exp(-(posterior$beta_cure_intercept + final_log_or_draws)))
  cure_diff_draws <- p_cure_exp_draws - p_cure_ctrl_draws

  # --- Return the final list of draws ---
  return(
    list(
      tr_posterior_samples = exp(final_log_tr_draws),
      cure_posterior_samples = cure_diff_draws
    )
  )
}
