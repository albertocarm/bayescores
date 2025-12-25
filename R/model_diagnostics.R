#' Diagnose Bayesian MCMC fit
#'
#' Extracts effective sample size (n_eff) and Rhat diagnostics from a fitted
#' Bayesian model (e.g., Stan fit) and provides categorical interpretations
#' along with practical recommendations for improvement.
#'
#' @param fit A fitted Bayesian model object (e.g., stanfit).
#' @param n_eff_threshold_high Numeric. Threshold above which n_eff is considered "high" (default = 1000).
#' @param n_eff_threshold_ok Numeric. Threshold above which n_eff is considered "moderate" (default = 500).
#' @param rhat_excellent Numeric. Threshold below which Rhat is considered "excellent" (default = 1.01).
#' @param rhat_ok Numeric. Threshold below which Rhat is considered "acceptable" (default = 1.05).
#' @param include_lp Logical. If FALSE, excludes lp from the output table (default = TRUE).
#'
#' @return A data frame with columns
#' \itemize{
#'   \item n_eff  effective sample size
#'   \item Rhat  Gelman-Rubin convergence diagnostic
#'   \item Rhat_interpretation  categorical interpretation of Rhat
#'   \item n_eff_interpretation  categorical interpretation of n_eff
#'   \item recommendation  practical suggestions if diagnostics are poor
#' }
#'
#' @examples
#' \dontrun{
#' diagnose_fit(bayesian_fit_kn199$stan_fit)
#' diagnose_fit(bayesian_fit_kn199$stan_fit, include_lp = FALSE)
#' }
#'
#' @export
diagnose_fit <- function(fit,
                         n_eff_threshold_high = 1000,
                         n_eff_threshold_ok   = 500,
                         rhat_excellent       = 1.01,
                         rhat_ok              = 1.05,
                         include_lp           = TRUE) {
  # extract summary
  summ <- rstan::summary(fit)$summary[, c("n_eff", "Rhat")]
  df <- as.data.frame(summ)

  # Excluir lp__ si se pide
  if (!include_lp) {
    df <- df[rownames(df) != "lp__", , drop = FALSE]
  }

  #  Rhat interpretation
  Rhat_interpret <- function(x) {
    if (is.na(x)) return("Not available")
    if (x <= rhat_excellent) return("Excellent convergence")
    if (x <= rhat_ok)        return("Acceptable convergence")
    if (x <= 1.10)           return("Potential issues, check chains")
    return("Poor convergence, re-run / reparameterize")
  }

  #  n_eff interpretation
  n_eff_interpret <- function(x) {
    if (is.na(x)) return("Not available")
    if (x >= n_eff_threshold_high) return("High effective sample size")
    if (x >= n_eff_threshold_ok)   return("Moderate, acceptable")
    if (x >= 100)                  return("Low, estimates less precise")
    return("Very low, estimates unreliable")
  }

  # Recomendaciones
  recommend <- function(rhat, neff) {
    bad_rhat <- (!is.na(rhat)) && (rhat > rhat_ok)
    very_bad_rhat <- (!is.na(rhat)) && (rhat > 1.10)
    low_neff <- (!is.na(neff)) && (neff < n_eff_threshold_ok)
    very_low_neff <- (!is.na(neff)) && (neff < 100)

    if (very_bad_rhat || very_low_neff) {
      return(paste(
        "Strongly consider:",
        "increase iterations; set adapt_delta approx. 0.99;",
        "raise max_treedepth (15-20);",
        "non-centered reparameterization;",
        "rescale predictors; tighten priors;",
        "inspect divergences and BFMI."
      ))
    }
    if (bad_rhat || low_neff) {
      return(paste(
        "Consider:",
        "increase iterations; set adapt_delta >= 0.95;",
        "raise max_treedepth (12-15);",
        "check divergences;",
        "rescale parameters;",
        "try non-centered parameterization and informative priors."
      ))
    }
    return("No action needed.")
  }

  # apply
  df$Rhat_interpretation   <- vapply(df$Rhat,  Rhat_interpret,  character(1))
  df$n_eff_interpretation  <- vapply(df$n_eff, n_eff_interpret, character(1))
  df$recommendation        <- mapply(recommend, df$Rhat, df$n_eff, USE.NAMES = FALSE)

  # reorder
  df <- df[, c("n_eff", "Rhat", "Rhat_interpretation", "n_eff_interpretation", "recommendation")]
  return(df)
}


#' @title Model Diagnostics
#' @description Generic function to run diagnostics on fitted models.
#' @param x A fitted model object.
#' @param ... Additional arguments passed to methods.
#' @export
model_diagnostics <- function(x, ...) {
  UseMethod("model_diagnostics")
}




#' @method model_diagnostics bcm_fit
#' @importFrom rstan stan_trace
#' @export
model_diagnostics.bcm_fit <- function(x, ...) {
  stan_fit <- x$stan_fit

  # Define base parameters that are always present
  # Note: alpha has been renamed to alpha_control in the new model
  params_to_check <- c("beta_cure_intercept", "beta_cure_arm",
                       "beta_surv_intercept", "beta_surv_arm",
                       "alpha_control")

  # Only include the shape treatment effect if shared_shape is FALSE
  if (!isTRUE(x$shared_shape)) {
    params_to_check <- c(params_to_check, "beta_alpha_arm")
  }

  # Generate and print the trace plots for the selected parameters
  trace_plot <- rstan::stan_trace(stan_fit, pars = params_to_check)
  print(trace_plot)
}









