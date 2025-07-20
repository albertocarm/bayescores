

#' @title Summarize Cure Rate Results
#' @description Calculates the posterior estimates of the cure rates for each
#' treatment arm and the absolute difference between them.
#'
#' @param model An object of class 'bcm_fit' returned by
#'   `fit_bayesian_cure_model`.
#'
#' @return A list containing the estimated cure rates for each arm (named
#'   dynamically using the factor levels), their absolute difference, and an
#'   `interpretation` string.
#' @importFrom rstan extract
#' @importFrom stats median quantile
#' @export
summarize_cure_rates <- function(model) {
  # Check if the input is the correct class
  if (!inherits(model, "bcm_fit")) {
    stop("Input must be an object of class 'bcm_fit'.")
  }

  # --- NEW: Get factor levels from original data ---
  arm_col_name <- model$column_map$arm_col
  original_data <- model$original_data
  arm_levels <- levels(original_data[[arm_col_name]])
  if (length(arm_levels) != 2) {
    stop("This function currently supports only two treatment arms.")
  }
  level_arm1 <- arm_levels[1]
  level_arm2 <- arm_levels[2]
  # -----------------------------------------------

  # Helper function to convert log-odds to probability
  sigmoid <- function(x) { 1 / (1 + exp(-x)) }

  # Extract posterior samples
  posterior <- rstan::extract(model$stan_fit)

  # --- Calculate posterior distribution for each quantity ---
  p_cure_arm1_samples <- sigmoid(posterior$beta_cure_intercept)
  p_cure_arm2_samples <- sigmoid(posterior$beta_cure_intercept + posterior$beta_cure_arm)
  p_diff_samples <- p_cure_arm2_samples - p_cure_arm1_samples

  # --- Summarize each distribution ---
  res_arm1 <- list(estimate = stats::median(p_cure_arm1_samples), ci_95 = stats::quantile(p_cure_arm1_samples, probs = c(0.025, 0.975)))
  res_arm2 <- list(estimate = stats::median(p_cure_arm2_samples), ci_95 = stats::quantile(p_cure_arm2_samples, probs = c(0.025, 0.975)))
  res_diff <- list(estimate = stats::median(p_diff_samples), ci_95 = stats::quantile(p_diff_samples, probs = c(0.025, 0.975)))

  # --- NEW: Add automatic interpretation for the difference ---
  if (res_diff$ci_95[1] > 0) {
    interpretation <- paste("Credible evidence that the cure rate is higher in the '", level_arm2, "' arm compared to the '", level_arm1, "' arm.", sep="")
  } else if (res_diff$ci_95[2] < 0) {
    interpretation <- paste("Credible evidence that the cure rate is lower in the '", level_arm2, "' arm compared to the '", level_arm1, "' arm.", sep="")
  } else {
    interpretation <- "No credible evidence of a difference in cure rates between the arms (Credible Interval for the difference includes 0)."
  }
  # -------------------------------------------------------------

  # --- Combine into the final list using dynamic names ---
  final_result <- list()
  final_result[[paste0("cure_rate_arm_", level_arm1)]] <- res_arm1
  final_result[[paste0("cure_rate_arm_", level_arm2)]] <- res_arm2
  final_result$absolute_difference <- res_diff
  final_result$interpretation <- interpretation # Add interpretation to the output

  return(final_result)
}


#' @title Extract MCMC Samples for the Time Ratio
#' @description Extracts the full vector of posterior MCMC samples for the Time
#' Ratio (TR), allowing for custom plotting and analysis.
#'
#' @param model An object of class 'bcm_fit' returned by
#'   `fit_bayesian_cure_model`.
#'
#' @return A numeric vector containing the posterior samples of the Time Ratio.
#' @importFrom rstan extract
#' @export
extract_mcmc_time_ratios <- function(model) {
  # Input validation
  if (!inherits(model, "bcm_fit")) {
    stop("Input must be an object of class 'bcm_fit'.")
  }

  # Extract all posterior samples
  posterior <- rstan::extract(model$stan_fit)

  # Calculate the Time Ratio for each sample
  tr_samples <- exp(posterior$beta_surv_arm)

  return(tr_samples)
}



#' @title Extract MCMC Samples for the Difference in Cure Rates
#' @description Extracts the full vector of posterior MCMC samples for the
#' absolute difference in cure rates between the treatment and control arms.
#'
#' @param model An object of class 'bcm_fit' returned by
#'   `fit_bayesian_cure_model`.
#'
#' @return A numeric vector containing the posterior samples of the difference
#'   in cure rates (treatment - control).
#' @importFrom rstan extract
#' @export
extract_mcmc_cure_diffs <- function(model) {
  # Input validation
  if (!inherits(model, "bcm_fit")) {
    stop("Input must be an object of class 'bcm_fit'.")
  }

  # Helper function to convert log-odds to probability
  sigmoid <- function(x) { 1 / (1 + exp(-x)) }

  # Extract all posterior samples
  posterior <- rstan::extract(model$stan_fit)

  # Calculate cure rate posterior for the treatment arm (arm=1)
  p_cure_treatment_samples <- sigmoid(posterior$beta_cure_intercept + posterior$beta_cure_arm)

  # Calculate cure rate posterior for the control arm (arm=0)
  p_cure_control_samples <- sigmoid(posterior$beta_cure_intercept)

  # Calculate the posterior for the difference
  difference_samples <- p_cure_treatment_samples - p_cure_control_samples

  return(difference_samples)
}
