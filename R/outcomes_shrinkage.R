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
#'    Options are "none" (default) or "skeptical".
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
  valid_methods <- c("skeptical", "none")
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
  log_sr <- as.numeric(posterior$beta_alpha_arm)

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
    Metric = c("Time Ratio (TR)", "Odds Ratio (OR) for Cure", "Shape Ratio (SR)",
               "Long-Term Survival Rate (%) - Control",
               "Long-Term Survival Rate (%) - Experimental", "Absolute Difference in Survival Rate (%)"),
    `Result (95% CI)` = c(
      summarize_draws(log_tr, is_log_scale = TRUE),
      summarize_draws(log_or, is_log_scale = TRUE),
      summarize_draws(log_sr, is_log_scale = TRUE),
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

    summary_df <- rbind(
      summary_df[1,], if(!identical(cor_log_tr, log_tr)) corrected_rows[1,], # TR
      summary_df[2,], if(!identical(cor_log_or, log_or)) corrected_rows[2,], # OR
      summary_df[3,], # SR
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
