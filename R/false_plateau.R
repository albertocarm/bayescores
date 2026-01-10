#' Check for False Plateau or Structural Instability
#'
#' @description
#' Analyzes the estimated long-term survival fraction (plateau) to distinguish between:
#' \itemize{
#'   \item Absence of plateau (below threshold).
#'   \item Robust plateau (high estimate, 0% probability of being an artifact).
#'   \item False/Unstable plateau (visual plateau but high probability of being an artifact).
#' }
#'
#' @param fit An object containing the fitted Bayesian model (expected to have a \code{stan_fit} slot).
#' @param threshold A numeric value defining the boundary for clinical relevance.
#' Defaults to 0.05 (5\%).
#'
#' @return Prints a context-aware diagnostic summary and invisibly returns a list containing
#' the calculated probabilities and mean rates for both arms.
#'
#' @importFrom rstan extract
#' @importFrom stats plogis
#' @export
check_false_plateau <- function(fit, threshold = 0.05) {

  if (!"stan_fit" %in% names(fit)) stop("Object must contain a 'stan_fit' slot.")

  samples <- rstan::extract(fit$stan_fit)

  if (is.null(samples$beta_cure_intercept) || is.null(samples$beta_cure_arm)) {
    stop("Model missing required cure parameters.")
  }

  # Calculate Probabilities
  # We use 'Long-Term Fraction' instead of 'Cure' to be neutral
  rate_ctrl <- stats::plogis(samples$beta_cure_intercept)
  rate_exp  <- stats::plogis(samples$beta_cure_intercept + samples$beta_cure_arm)

  # Helper function to print dynamic status per arm
  print_arm_status <- function(name, rates, thresh) {
    mean_val <- mean(rates)
    prob_false <- mean(rates < thresh)

    cat(sprintf("%s:\n", name))
    cat(sprintf("  Est. Long-Term Fraction:  %.2f%%\n", mean_val * 100))

    # LOGIC:
    # 1. If the mean is below threshold -> No plateau exists to be judged.
    if (mean_val < thresh) {
      cat("  Status:                   NO PLATEAU DETECTED (Below threshold)\n")
    } else {
      # 2. If mean is above threshold -> Check if it's solid or false
      cat(sprintf("  Prob. of False Plateau:   %.1f%% (< %.1f%%)\n", prob_false * 100, thresh*100))

      if (prob_false < 0.01) {
        cat("  Status:                   ROBUST PLATEAU (Structurally Stable)\n")
      } else if (prob_false > 0.20) {
        cat("  Status:                   UNSTABLE / FALSE PLATEAU (Artifact likely)\n")
      } else {
        cat("  Status:                   BORDERLINE / AMBIGUOUS\n")
      }
    }
    cat("\n")
    return(list(mean = mean_val, prob_false = prob_false))
  }

  cat("--- PLATEAU STRUCTURAL CHECK ---\n")
  cat(sprintf("Clinical Threshold: %.1f%%\n\n", threshold * 100))

  res_ctrl <- print_arm_status("CONTROL ARM", rate_ctrl, threshold)
  res_exp  <- print_arm_status("EXPERIMENTAL ARM", rate_exp, threshold)

  invisible(list(control = res_ctrl, experimental = res_exp))
}
