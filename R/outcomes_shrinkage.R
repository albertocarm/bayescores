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
#' Generate a Summary Table of Model Outcomes (Final Version 2.1)
#'
#' @description
#' Extracts data directly from the fit object to calculate and display key model
#' outcomes, now including the Shape Ratio (SR) to account for flexible hazard dynamics.
#'
#' @param fit The fitted model object from `fit_bayesian_cure_model`.
#' @param calibration_args An optional list to override default utility calibration.
#' @param digits The number of decimal places for rounding.
#' @param shrinkage_method A character string specifying the shrinkage method.
#'    Options are "none" (default), "zwet", "sherry", or "skeptical".
#' @param shrinkage_target A character string controlling the application of
#'    shrinkage. Options are "primary" or "both_if_uncorrelated".
#' @param correlation_method Character. The method to use: 'pearson' (default),
#'    'spearman', or 'kendall'. Used for the identifiability check.
#'
#' @return A tibble with a summary of key model outcomes.
#'
#' @importFrom rstan extract
#' @importFrom tibble tibble as_tibble
#' @importFrom stats quantile median sd rnorm dnorm cov cor
#' @importFrom tools toTitleCase
#' @importFrom dplyr case_when
#' @export
#'
outcomes <- function(fit,
                     calibration_args = list(),
                     digits = 2,
                     shrinkage_method = "none",
                     shrinkage_target = "primary",
                     correlation_method = c("pearson", "spearman", "kendall")) {

  # --- 1) Input Validation and Data Extraction ---
  valid_methods <- c("zwet", "sherry", "skeptical", "none")
  if (!shrinkage_method %in% valid_methods) {
    stop(paste("Invalid shrinkage_method. Use one of:", paste(valid_methods, collapse = ", ")))
  }
  valid_targets <- c("primary", "both_if_uncorrelated")
  if (!shrinkage_target %in% valid_targets) {
    stop(paste("Invalid shrinkage_target. Use one of:", paste(valid_targets, collapse = ", ")))
  }

  correlation_method <- match.arg(correlation_method)

  posterior <- rstan::extract(fit$stan_fit)
  log_or <- as.numeric(posterior$beta_cure_arm)
  log_tr <- as.numeric(posterior$beta_surv_arm)
  log_sr <- as.numeric(posterior$beta_alpha_arm) # NUEVO: Shape Ratio

  if (is.null(log_tr) || is.null(log_or)) {
    stop("`fit` object must contain `beta_cure_arm` and `beta_surv_arm` draws.")
  }

  rho <- cor(log_or, log_tr, method = correlation_method)

  # --- 2) Helper functions ---
  summarize_draws <- function(draws, is_log_scale = FALSE, multiplier = 1) {
    vals <- if (is_log_scale) exp(draws) else draws
    med <- stats::median(vals)
    ci  <- stats::quantile(vals, probs = c(0.025, 0.975))
    paste0(sprintf("%.*f", digits, med * multiplier), " (",
           sprintf("%.*f", digits, ci[1] * multiplier), " - ",
           sprintf("%.*f", digits, ci[2] * multiplier), ")")
  }

  # Helper for skeptical shrinkage on whitened draws
  apply_skeptical_shrinkage <- function(log_or_draws, log_tr_draws) {
    theta <- cbind(as.numeric(log_or_draws), as.numeric(log_tr_draws))
    mu <- colMeans(theta)
    S  <- stats::cov(theta)
    eg <- eigen(S, symmetric = TRUE)

    lam <- base::pmax(eg$values, 1e-12)
    Q   <- eg$vectors

    W    <- Q %*% diag(1/sqrt(lam)) %*% t(Q)
    Winv <- Q %*% diag(sqrt(lam)) %*% t(Q)

    Z <- sweep(theta, 2, mu, "-") %*% W
    mz_sh <- as.numeric(W %*% mu) / 2
    Z_sh  <- Z / sqrt(2)

    theta_sh <- sweep(Z_sh %*% Winv, 2, as.numeric(Winv %*% mz_sh), "+")

    list(
      cor_log_or = as.numeric(theta_sh[, 1]),
      cor_log_tr = as.numeric(theta_sh[, 2])
    )
  }

  # --- 3) Base Table (Original Results) ---
  p_ctrl <- 1 / (1 + exp(-posterior$beta_cure_intercept))
  p_exp_orig <- 1 / (1 + exp(-(posterior$beta_cure_intercept + log_or)))
  diff_orig  <- p_exp_orig - p_ctrl

  summary_df <- data.frame(
    Metric = c("Time Ratio (TR)", "Odds Ratio (OR) for Cure", "Shape Ratio (SR)", # NUEVO: SR añadido
               "Long-Term Survival Rate (%) - Control",
               "Long-Term Survival Rate (%) - Experimental", "Absolute Difference in Survival Rate (%)"),
    `Result (95% CI)` = c(
      summarize_draws(log_tr, is_log_scale = TRUE),
      summarize_draws(log_or, is_log_scale = TRUE),
      summarize_draws(log_sr, is_log_scale = TRUE), # NUEVO: SR draws
      summarize_draws(p_ctrl, multiplier = 100),
      summarize_draws(p_exp_orig, multiplier = 100),
      summarize_draws(diff_orig, multiplier = 100)
    ),
    check.names = FALSE, stringsAsFactors = FALSE
  )

  # --- 4) Shrinkage Application (if requested) ---
  if (shrinkage_method != "none") {
    cor_log_or <- log_or
    cor_log_tr <- log_tr
    shrink_label <- ""

    defaults <- list(efficacy = list(cure_utility_target = list(effect_value = 0.15, utility_value = 40), tr_utility_target = list(effect_value = 1.25, utility_value = 35)))
    recursive_merge <- function(default, user) { for (name in names(user)) { if (is.list(user[[name]]) && is.list(default[[name]])) { default[[name]] <- recursive_merge(default[[name]], user[[name]]) } else { default[[name]] <- user[[name]] } }; return(default) }
    cal <- recursive_merge(defaults, calibration_args)
    .create_exp_utility_fn <- function(target, baseline = 0) {
      effect_gain <- target$effect_value - baseline
      if (effect_gain <= 0) return(function(x) rep(0, length(x)))
      lambda <- -log(1 - target$utility_value / 100) / effect_gain
      return(function(x) 100 * (1 - exp(-lambda * base::pmax(0, x))))
    }
    fn_u_cure <- .create_exp_utility_fn(cal$efficacy$cure_utility_target, 0)
    fn_u_tr <- .create_exp_utility_fn(cal$efficacy$tr_utility_target, 1.0)
    median_u_cure <- median(fn_u_cure(diff_orig))
    median_u_tr <- median(fn_u_tr(exp(log_tr) - 1.0))
    primary_endpoint <- if (median_u_cure >= median_u_tr) "OR" else "TR"

    method_name <- tools::toTitleCase(shrinkage_method)
    apply_to_or <- FALSE; apply_to_tr <- FALSE
    is_low_corr <- abs(rho) < 0.30

    if (shrinkage_target == "primary") {
      if (primary_endpoint == "OR") apply_to_or <- TRUE else apply_to_tr <- TRUE
      shrink_label <- paste0(method_name, ", Primary Only (", primary_endpoint, ")")
    } else if (shrinkage_target == "both_if_uncorrelated") {
      if (is_low_corr) {
        apply_to_or <- TRUE; apply_to_tr <- TRUE
        shrink_label <- paste0(method_name, ", Both (Low Corr.)")
      } else {
        if (primary_endpoint == "OR") apply_to_or <- TRUE else apply_to_tr <- TRUE
        shrink_label <- paste0(method_name, ", Primary Only (", primary_endpoint, ", High Corr.)")
      }
    }

    if (shrinkage_method == "skeptical") {
      if (apply_to_or && apply_to_tr) {
        shrunk_draws <- apply_skeptical_shrinkage(log_or, log_tr)
        cor_log_or <- shrunk_draws$cor_log_or
        cor_log_tr <- shrunk_draws$cor_log_tr
      } else if (apply_to_or) {
        temp_skeptical_sampler <- function(draws) {
          b <- mean(draws); s <- sd(draws); if(s==0) return(draws)
          shrunk_b <- b / 2
          return((draws - b) / sqrt(2) + shrunk_b)
        }
        cor_log_or <- temp_skeptical_sampler(log_or)
      } else {
        temp_skeptical_sampler <- function(draws) {
          b <- mean(draws); s <- sd(draws); if(s==0) return(draws)
          shrunk_b <- b / 2
          return((draws - b) / sqrt(2) + shrunk_b)
        }
        cor_log_tr <- temp_skeptical_sampler(log_tr)
      }
    } else {
      if (apply_to_or) cor_log_or <- sample_shrinkage_draws(log_or, method = shrinkage_method)
      if (apply_to_tr) cor_log_tr <- sample_shrinkage_draws(log_tr, method = shrinkage_method)
    }

    p_exp_cor <- 1 / (1 + exp(-(posterior$beta_cure_intercept + cor_log_or)))
    diff_cor  <- p_exp_cor - p_ctrl
    corrected_rows <- data.frame(
      Metric = c(
        paste0("Time Ratio (TR) - Corrected (", shrink_label, ")"),
        paste0("Odds Ratio (OR) - Corrected (", shrink_label, ")"),
        paste0("Long-Term Survival Rate (%) - Exp - Corrected (", shrink_label, ")"),
        paste0("Abs. Diff. in Survival (%) - Corrected (", shrink_label, ")")
      ),
      `Result (95% CI)` = c(
        summarize_draws(cor_log_tr, is_log_scale = TRUE),
        summarize_draws(cor_log_or, is_log_scale = TRUE),
        summarize_draws(p_exp_cor, multiplier = 100),
        summarize_draws(diff_cor, multiplier = 100)
      ),
      check.names = FALSE, stringsAsFactors = FALSE
    )

    # AJUSTE DE ÍNDICES: Ahora que SR es la fila 3, las tasas empiezan en la 4.
    summary_df <- rbind(
      summary_df[1,], if(!identical(cor_log_tr, log_tr)) corrected_rows[1,], # TR
      summary_df[2,], if(!identical(cor_log_or, log_or)) corrected_rows[2,], # OR
      summary_df[3,], # SR (Sin corrección aplicada)
      summary_df[4,], # Control Rate
      summary_df[5,], if(!identical(cor_log_or, log_or)) corrected_rows[3,], # Exp Rate
      summary_df[6,], if(!identical(cor_log_or, log_or)) corrected_rows[4,]  # Diff
    )
  }

  # --- 5) Print Identifiability Check and Return Tibble ---
  identifiability_level <- dplyr::case_when(
    rho > -0.15 ~ "Negligible (0 to -0.15): effects separable, estimation robust.",
    rho > -0.30 ~ "Low (-0.15 to -0.30): minor coupling, interpretation still reliable.",
    rho > -0.50 ~ "Moderate (-0.30 to -0.50): noticeable trade-off, joint summaries recommended.",
    rho > -0.70 ~ "Strong (-0.50 to -0.70): high ambiguity, marginal estimates fragile.",
    TRUE         ~ "Extreme (< -0.70): practical non-identifiability, requires longer follow-up or informative priors."
  )
  cat("\n---\n")
  cat("Posterior correlation as an identifiability check:\n")

  cor_symbol <- switch(correlation_method, "pearson" = "r", "spearman" = "rho", "kendall" = "tau")
  cat(paste0("    Correlation (", correlation_method, "): ", cor_symbol, " = ", round(rho, 3), "\n"))

  if (correlation_method %in% c("pearson", "spearman")) {
    rho2 <- round(rho^2 * 100, 1)
    cat(paste0("    Explained Variance (", cor_symbol, "^2): ", round(rho^2, 3),
               " (approx. ", rho2, "% information overlap)\n"))
    cat(paste0("    Note: Roughly ", rho2,
               "% of the cure variance is structurally indistinguishable from the time effect.\n"))
  }

  cat(paste0("    Interpretation: ", identifiability_level, "\n"))
  cat("---\n")
  tibble::as_tibble(summary_df)
}
