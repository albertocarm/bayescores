#' Generate a Summary Table with Optional Empirical Shrinkage
#'
#' @description
#' This function calculates key outcomes from a Bayesian model fit and provides
#' an option to include empirically shrunk estimates for the Time Ratio and
#' Odds Ratio based on the method by Van Zwet et al. (2021).
#'
#' @param fit The fitted model object from `fit_bayesian_cure_model`.
#' @param digits The number of decimal places for rounding.
#' @param empirical_shrinkage Logical. If TRUE, adds corrected estimates for
#'   TR and OR to the output table. Defaults to TRUE.
#'
#' @return A tibble with key model outcomes. If `empirical_shrinkage` is TRUE,
#'   includes additional rows with corrected estimates.
#'
#' @importFrom rstan extract
#' @importFrom tibble tibble
#' @importFrom stats quantile median sd rnorm
#' @export
outcomes_shrinkage <- function(fit, digits = 2, empirical_shrinkage = TRUE) {

  # ===================================================================
  # Helper Function 1: Van Zwet's `conditional` function
  # Computes the parameters of the corrected posterior distribution for the SNR.
  # This code is directly from the appendix of Van Zwet et al. (2021).
  # ===================================================================
  conditional <- function(z, p, tau) {
    tau2 <- tau^2
    q_unscaled <- p * dnorm(z, 0, sqrt(tau2 + 1))
    q <- q_unscaled / sum(q_unscaled) # conditional mixing probs
    m <- z * tau2 / (tau2 + 1)         # conditional means
    v <- tau2 / (tau2 + 1)             # conditional variances
    sigma <- sqrt(v)                   # conditional std devs
    data.frame(q, m, sigma)
  }

  # ===================================================================
  # Helper Function 2: The New Sampler
  # Takes original MCMC draws, applies the shrinkage, and returns a new
  # "corrected" vector of MCMC draws.
  # ===================================================================
  sample_shrinkage_draws <- function(original_draws) {
    # 1. Summarize the original posterior to get a single estimate (b) and error (s)
    b <- mean(original_draws)
    s <- sd(original_draws)
    z_obs <- b / s
    n_samples <- length(original_draws)

    # 2. Define the empirical prior parameters from Van Zwet et al. (Table 1)
    p_zwet <- c(0.32, 0.31, 0.30, 0.07)
    tau_zwet <- c(0.61, 1.42, 2.16, 5.64)

    # 3. Get the parameters of the corrected posterior distribution for the SNR
    corrected_dist_params <- conditional(z = z_obs, p = p_zwet, tau = tau_zwet)
    q_corr <- corrected_dist_params$q
    m_corr <- corrected_dist_params$m
    sigma_corr <- corrected_dist_params$sigma

    # 4. Generate new samples from this corrected mixture distribution
    corrected_snr_samples <- numeric(n_samples)
    for (i in 1:n_samples) {
      # 4a. Choose one of the 4 components based on the new probabilities (q_corr)
      component <- sample(1:4, size = 1, prob = q_corr)
      # 4b. Draw a sample from the normal distribution of that chosen component
      corrected_snr_samples[i] <- rnorm(1, mean = m_corr[component], sd = sigma_corr[component])
    }

    # 5. Scale the corrected SNR samples by 's' to get the final corrected draws for the parameter
    corrected_draws <- s * corrected_snr_samples
    return(corrected_draws)
  }

  # --- Main Function Logic ---

  # 1. Extract posterior samples from the fit object
  posterior <- rstan::extract(fit$stan_fit)
  if (!"beta_cure_intercept" %in% names(posterior)) {
    stop("Required parameters not found in the model fit.")
  }

  # 2. Calculate original (uncorrected) derived quantities
  # Note: these are on the log scale for TR and OR
  log_tr_draws <- posterior$beta_surv_arm
  log_or_draws <- posterior$beta_cure_arm

  # Other quantities that don't get corrected
  p_cure_ctrl_draws <- 1 / (1 + exp(-posterior$beta_cure_intercept))
  p_cure_exp_draws <- 1 / (1 + exp(-(posterior$beta_cure_intercept + posterior$beta_cure_arm)))
  cure_diff_draws <- p_cure_exp_draws - p_cure_ctrl_draws

  # 3. Generate corrected draws if requested
  if (empirical_shrinkage) {
    log_tr_draws_corrected <- sample_shrinkage_draws(log_tr_draws)
    log_or_draws_corrected <- sample_shrinkage_draws(log_or_draws)
  }

  # 4. Helper function to summarize any vector of draws into a formatted string
  summarize_draws <- function(draws, is_log_scale = FALSE, multiplier = 1) {
    # If draws are on log scale (like log-TR), exponentiate them first
    if (is_log_scale) {
      summary_draws <- exp(draws)
    } else {
      summary_draws <- draws
    }

    point_estimate <- median(summary_draws) * multiplier
    ci <- quantile(summary_draws, probs = c(0.025, 0.975)) * multiplier

    paste0(
      format(round(point_estimate, digits), nsmall = digits),
      " (", format(round(ci[1], digits), nsmall = digits), " - ", format(round(ci[2], digits), nsmall = digits), ")"
    )
  }

  # 5. Build the results table
  # Start with a base table of metrics
  summary_df <- data.frame(
    Metric = c("Time Ratio (TR)",
               "Odds Ratio (OR) for Cure",
               "Long-Term Survival Rate (%) - Control",
               "Long-Term Survival Rate (%) - Experimental",
               "Absolute Difference in Survival Rate (%)"),
    `Result (95% CI)` = c(
      summarize_draws(log_tr_draws, is_log_scale = TRUE),
      summarize_draws(log_or_draws, is_log_scale = TRUE),
      summarize_draws(p_cure_ctrl_draws, multiplier = 100),
      summarize_draws(p_cure_exp_draws, multiplier = 100),
      summarize_draws(cure_diff_draws, multiplier = 100)
    ),
    stringsAsFactors = FALSE
  )

  # If shrinkage was applied, add the corrected rows
  if (empirical_shrinkage) {
    corrected_rows <- data.frame(
      Metric = c("Time Ratio (TR) - Corrected",
                 "Odds Ratio (OR) - Corrected"),
      `Result (95% CI)` = c(
        summarize_draws(log_tr_draws_corrected, is_log_scale = TRUE),
        summarize_draws(log_or_draws_corrected, is_log_scale = TRUE)
      )
    )
    # Combine the original table with the new corrected rows
    summary_df <- rbind(summary_df[1:2,], corrected_rows, summary_df[3:5,])
  }

  # Convert to a tibble for nice printing
  summary_table <- tibble::as_tibble(summary_df)

  return(summary_table)
}





#' Generate a Summary Table or Return MCMC Draws with Optional Shrinkage
#'
#' @description
#' This function calculates key outcomes from a Bayesian model fit. It can either
#' return a summary table or a list of MCMC draws. It also provides an option
#' to include empirically shrunk estimates based on Van Zwet et al. (2021).
#'
#' @param fit The fitted model object from `fit_bayesian_cure_model`.
#' @param digits The number of decimal places for rounding in the summary table.
#' @param empirical_shrinkage Logical. If TRUE, corrected estimates/draws are generated.
#' @param return_draws Logical. If TRUE, the function returns a list of MCMC
#'   draws instead of a summary table. Defaults to FALSE.
#'
#' @return If `return_draws` is FALSE (default), a tibble with a summary of outcomes.
#'   If `return_draws` is TRUE, a list containing vectors of MCMC draws for
#'   various parameters, including corrected versions if requested.
#'
#' @importFrom rstan extract
#' @importFrom tibble tibble as_tibble
#' @importFrom stats quantile median sd rnorm
#' @importFrom stats dnorm
#' @export
outcomes_shrinkage <- function(fit, digits = 2, empirical_shrinkage = TRUE, return_draws = FALSE) {

  # ===================================================================
  # Helper Function 1: Van Zwet's `conditional` function
  # ===================================================================
  conditional <- function(z, p, tau) {
    tau2 <- tau^2
    q_unscaled <- p * dnorm(z, 0, sqrt(tau2 + 1))
    q <- q_unscaled / sum(q_unscaled)
    m <- z * tau2 / (tau2 + 1)
    v <- tau2 / (tau2 + 1)
    sigma <- sqrt(v)
    data.frame(q, m, sigma)
  }

  # ===================================================================
  # Helper Function 2: The New Sampler
  # ===================================================================
  sample_shrinkage_draws <- function(original_draws) {
    b <- mean(original_draws)
    s <- sd(original_draws)
    if (s == 0) return(original_draws) # No variation, no shrinkage
    z_obs <- b / s
    n_samples <- length(original_draws)
    p_zwet <- c(0.32, 0.31, 0.30, 0.07)
    tau_zwet <- c(0.61, 1.42, 2.16, 5.64)
    corrected_dist_params <- conditional(z = z_obs, p = p_zwet, tau = tau_zwet)
    q_corr <- corrected_dist_params$q
    m_corr <- corrected_dist_params$m
    sigma_corr <- corrected_dist_params$sigma
    corrected_snr_samples <- numeric(n_samples)
    for (i in 1:n_samples) {
      component <- sample(1:4, size = 1, prob = q_corr)
      corrected_snr_samples[i] <- rnorm(1, mean = m_corr[component], sd = sigma_corr[component])
    }
    corrected_draws <- s * corrected_snr_samples
    return(corrected_draws)
  }

  # --- Main Function Logic ---

  # 1. Extract posterior samples
  posterior <- rstan::extract(fit$stan_fit)
  if (!"beta_cure_intercept" %in% names(posterior)) {
    stop("Required parameters not found in the model fit.")
  }

  # 2. Calculate original (uncorrected) derived quantities
  log_tr_draws <- posterior$beta_surv_arm
  log_or_draws <- posterior$beta_cure_arm
  p_cure_ctrl_draws <- 1 / (1 + exp(-posterior$beta_cure_intercept))
  p_cure_exp_draws <- 1 / (1 + exp(-(posterior$beta_cure_intercept + log_or_draws)))
  cure_diff_draws <- p_cure_exp_draws - p_cure_ctrl_draws

  # Initialize list for returning draws
  draws_list <- list(
    log_tr_draws = log_tr_draws,
    log_or_draws = log_or_draws,
    cure_diff_draws = cure_diff_draws,
    p_cure_ctrl_draws = p_cure_ctrl_draws,
    p_cure_exp_draws = p_cure_exp_draws
  )

  # 3. Generate corrected draws if requested
  if (empirical_shrinkage) {
    log_tr_draws_corrected <- sample_shrinkage_draws(log_tr_draws)
    log_or_draws_corrected <- sample_shrinkage_draws(log_or_draws)

    # IMPORTANT: Recalculate derived quantities using the corrected draws
    p_cure_exp_draws_corrected <- 1 / (1 + exp(-(posterior$beta_cure_intercept + log_or_draws_corrected)))
    cure_diff_draws_corrected <- p_cure_exp_draws_corrected - p_cure_ctrl_draws

    # Add corrected draws to the list
    draws_list$log_tr_draws_corrected <- log_tr_draws_corrected
    draws_list$log_or_draws_corrected <- log_or_draws_corrected
    draws_list$cure_diff_draws_corrected <- cure_diff_draws_corrected
  }

  # --- Return either the draws or a summary table ---

  if (return_draws) {
    return(draws_list)
  } else {
    # Helper function to summarize draws into a formatted string
    summarize_draws <- function(draws, is_log_scale = FALSE, multiplier = 1) {
      summary_draws <- if (is_log_scale) exp(draws) else draws
      point_estimate <- median(summary_draws) * multiplier
      ci <- quantile(summary_draws, probs = c(0.025, 0.975)) * multiplier
      paste0(
        format(round(point_estimate, digits), nsmall = digits),
        " (", format(round(ci[1], digits), nsmall = digits), " - ", format(round(ci[2], digits), nsmall = digits), ")"
      )
    }

    # Build the results table
    summary_df <- data.frame(
      Metric = c("Time Ratio (TR)", "Odds Ratio (OR) for Cure", "Long-Term Survival Rate (%) - Control",
                 "Long-Term Survival Rate (%) - Experimental", "Absolute Difference in Survival Rate (%)"),
      `Result (95% CI)` = c(
        summarize_draws(draws_list$log_tr_draws, is_log_scale = TRUE),
        summarize_draws(draws_list$log_or_draws, is_log_scale = TRUE),
        summarize_draws(draws_list$p_cure_ctrl_draws, multiplier = 100),
        summarize_draws(draws_list$p_cure_exp_draws, multiplier = 100),
        summarize_draws(draws_list$cure_diff_draws, multiplier = 100)
      ),
      stringsAsFactors = FALSE
    )

    if (empirical_shrinkage) {
      corrected_rows <- data.frame(
        Metric = c("Time Ratio (TR) - Corrected", "Odds Ratio (OR) - Corrected"),
        `Result (95% CI)` = c(
          summarize_draws(draws_list$log_tr_draws_corrected, is_log_scale = TRUE),
          summarize_draws(draws_list$log_or_draws_corrected, is_log_scale = TRUE)
        )
      )
      summary_df <- rbind(summary_df[1:2,], corrected_rows, summary_df[3:5,])
    }

    return(tibble::as_tibble(summary_df))
  }
}




