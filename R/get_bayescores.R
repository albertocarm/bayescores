#' @title Calculate BayeScores Final Utility Score (Version 4, Corrected Toxicity Model)
#' @description
#' This version corrects the toxicity calculation model as per user feedback.
#' Toxicity's impact is now proportional to the base efficacy score.
#'
#' @param efficacy_inputs A list with `cure_posterior_samples`, `tr_posterior_samples`.
#' @param qol_scores A numeric vector of posterior samples for QoL effect.
#' @param toxicity_scores A numeric vector of posterior samples for Toxicity effect.
#' @param calibration_args Optional list to override default calibration.
#' @return A list with the final utility vector and a summary data frame.
#' @export
#' @importFrom stats median quantile
get_bayescores <- function(
    efficacy_inputs, qol_scores, toxicity_scores, calibration_args = list()
) {

  # --- Step 1 & 2: No changes ---
  defaults <- list(efficacy = list(cure_utility_target = list(effect_value = 0.15, utility_value = 40), tr_utility_target = list(effect_value = 1.25, utility_value = 35), tr_disutility_target = list(effect_value = 0.75, utility_value = 35), tr_baseline = 1.0))
  recursive_merge <- function(default, user) { for (name in names(user)) { if (is.list(user[[name]]) && is.list(default[[name]])) { default[[name]] <- recursive_merge(default[[name]], user[[name]]) } else { default[[name]] <- user[[name]] } }; return(default) }
  cal <- recursive_merge(defaults, calibration_args)
  cure_samples <- base::pmax(0, efficacy_inputs$cure_posterior_samples)
  tr_samples <- efficacy_inputs$tr_posterior_samples

  # --- Step 3 & 4: No changes ---
  .create_exp_utility_fn <- function(target, baseline = 0) {
    effect_gain <- target$effect_value - baseline
    if (effect_gain <= 0) stop(paste("Target effect_value", target$effect_value, "must be greater than baseline", baseline))
    lambda <- -log(1 - target$utility_value / 100) / effect_gain
    return(function(x) 100 * (1 - exp(-lambda * base::pmax(0, x))))
  }
  fn_u_cure <- .create_exp_utility_fn(cal$efficacy$cure_utility_target, 0)
  fn_u_tr <- .create_exp_utility_fn(cal$efficacy$tr_utility_target, cal$efficacy$tr_baseline)
  disutility_target_for_fn <- cal$efficacy$tr_disutility_target
  disutility_target_for_fn$effect_value <- cal$efficacy$tr_baseline - disutility_target_for_fn$effect_value
  fn_du_tr <- .create_exp_utility_fn(disutility_target_for_fn, 0)
  U_cure <- fn_u_cure(cure_samples)
  is_beneficial <- tr_samples >= cal$efficacy$tr_baseline
  U_tr <- fn_u_tr(tr_samples - cal$efficacy$tr_baseline) * is_beneficial
  DU_tr <- fn_du_tr(cal$efficacy$tr_baseline - tr_samples) * !is_beneficial
  bonus_term <- (pmin(U_tr, U_cure) / 100) * (100 - base::pmax(U_tr, U_cure))
  penalty_term <- (DU_tr / 100) * U_cure
  score_post_eficacia <- base::pmax(0, base::pmax(U_tr, U_cure) + bonus_term - penalty_term)

  # --- Step 5: Incorporate QoL and Toxicity ---
  qol_effect <- qol_scores / 2
  qol_power <- qol_effect * (score_post_eficacia / 100)
  score_post_qol <- score_post_eficacia + qol_power * (100 - score_post_eficacia)

  tox_effect <- toxicity_scores
  tox_power <- tox_effect * (score_post_eficacia / 100)
  final_utility_vector <- score_post_qol + tox_power * (100 - score_post_qol)

  # --- Step 6: Summarize Components ---
  .summarize <- function(vec) { c(Median = median(vec, na.rm=TRUE), Lower_95_CrI=quantile(vec, .025, na.rm=TRUE), Upper_95_CrI=quantile(vec, .975, na.rm=TRUE)) }
  summary_list <- list(
    "Utility Cure (0-100)" = .summarize(U_cure),
    "Utility TR (for TR>1)" = .summarize(U_tr),
    "Penalty TR (for TR<1)" = .summarize(-DU_tr),
    "Efficacy Score (Combined)" = .summarize(score_post_eficacia),
    "QoL Contribution (points)" = .summarize(qol_power * (100 - score_post_eficacia)),
    "Toxicity Contribution (points)" = .summarize(tox_power * (100 - score_post_qol)),
    "FINAL UTILITY SCORE" = .summarize(final_utility_vector)
  )
  component_summary <- as.data.frame(do.call(rbind, summary_list)); component_summary <- cbind(Component = rownames(component_summary), component_summary)
  rownames(component_summary) <- NULL

  # --- Step 7: Identifiability Check ---
  rho <- cor(
    efficacy_inputs$cure_posterior_samples,
    efficacy_inputs$tr_posterior_samples
  )

  identifiability_level <- dplyr::case_when(
    rho > -0.15 ~ "Negligible",
    rho > -0.30 ~ "Low",
    rho > -0.50 ~ "Moderate",
    rho > -0.70 ~ "Strong",
    TRUE        ~ "Extreme"
  )

  # --- Step 8: Return list with the new element ---
  return(list(
    final_utility_vector = final_utility_vector,
    component_summary = component_summary,
    identifiability_level = identifiability_level
  ))
}



#' Get Original or Shrunk MCMC Draws for Subsequent Analyses (Updated)
#'
#' @description
#' Extracts or generates posterior draws for efficacy inputs. This version
#' incorporates advanced shrinkage methods, including a skeptical prior via
#' whitening, and strategically applies them based on the primary utility
#' endpoint and parameter correlation.
#'
#' @param fit The fitted model object from `fit_bayesian_cure_model`.
#' @param shrinkage_method A character string specifying the shrinkage method.
#'   Options are "zwet", "sherry", "skeptical", or "none" (default).
#' @param shrinkage_target A character string controlling shrinkage application.
#'   Options are "primary" or "both_if_uncorrelated".
#' @param calibration_args An optional list to override default utility calibration,
#'   used to determine the primary endpoint.
#'
#' @return A list containing two named vectors of posterior samples:
#'   \itemize{
#'     \item{\code{tr_posterior_samples}: Draws for the Time Ratio (already exponentiated).}
#'     \item{\code{cure_posterior_samples}: Draws for the absolute difference
#'       in cure rates.}
#'   }
#'
#' @importFrom rstan extract
#' @importFrom stats median cov cor
#' @importFrom tools toTitleCase
#' @export
get_bayescores_draws <- function(fit,
                                 shrinkage_method = "none",
                                 shrinkage_target = "primary",
                                 calibration_args = list()) {

  # --- Input Validation ---
  valid_methods <- c("zwet", "sherry", "skeptical", "none")
  if (!shrinkage_method %in% valid_methods) {
    stop(paste("Invalid shrinkage_method specified. Use one of:", paste(valid_methods, collapse = ", ")))
  }
  valid_targets <- c("primary", "both_if_uncorrelated")
  if (!shrinkage_target %in% valid_targets) {
    stop(paste("Invalid shrinkage_target. Use one of:", paste(valid_targets, collapse = ", ")))
  }

  # --- Extract Draws ---
  posterior <- rstan::extract(fit$stan_fit)
  if (!"beta_cure_intercept" %in% names(posterior)) {
    stop("Required parameters not found in the model fit.")
  }
  log_or_draws <- as.numeric(posterior$beta_cure_arm)
  log_tr_draws <- as.numeric(posterior$beta_surv_arm)

  # --- Determine final draws based on shrinkage method ---
  final_log_or_draws <- log_or_draws
  final_log_tr_draws <- log_tr_draws

  if (shrinkage_method != "none") {
    rho <- cor(log_or_draws, log_tr_draws)

    # --- Determine Primary Endpoint via Utility ---
    p_ctrl <- 1 / (1 + exp(-posterior$beta_cure_intercept))
    p_exp_orig <- 1 / (1 + exp(-(posterior$beta_cure_intercept + log_or_draws)))
    diff_orig  <- p_exp_orig - p_ctrl

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
    median_u_tr <- median(fn_u_tr(exp(log_tr_draws) - 1.0))
    primary_endpoint <- if (median_u_cure >= median_u_tr) "OR" else "TR"

    # --- Determine which parameters to shrink ---
    apply_to_or <- FALSE; apply_to_tr <- FALSE
    is_low_corr <- abs(rho) < 0.30

    if (shrinkage_target == "primary") {
      if (primary_endpoint == "OR") apply_to_or <- TRUE else apply_to_tr <- TRUE
    } else if (shrinkage_target == "both_if_uncorrelated") {
      if (is_low_corr) {
        apply_to_or <- TRUE; apply_to_tr <- TRUE
      } else {
        if (primary_endpoint == "OR") apply_to_or <- TRUE else apply_to_tr <- TRUE
      }
    }

    # --- Apply Shrinkage ---
    if (shrinkage_method == "skeptical") {
      if (apply_to_or && apply_to_tr) {
        theta <- cbind(log_or_draws, log_tr_draws); mu <- colMeans(theta); S  <- stats::cov(theta)
        eg <- eigen(S, symmetric = TRUE); lam <- base::pmax(eg$values, 1e-12); Q <- eg$vectors
        Winv <- Q %*% diag(sqrt(lam)) %*% t(Q); W <- Q %*% diag(1/sqrt(lam)) %*% t(Q)
        mz <- as.numeric(W %*% mu); mz_sh <- mz / 2; mu_sh <- as.numeric(Winv %*% mz_sh)
        shift <- mu_sh - mu
        final_log_or_draws <- log_or_draws + shift[1]
        final_log_tr_draws <- log_tr_draws + shift[2]
      } else if (apply_to_or) {
        b <- mean(log_or_draws); shrunk_b <- b / 2
        final_log_or_draws <- log_or_draws - (b - shrunk_b)
      } else { # apply_to_tr
        b <- mean(log_tr_draws); shrunk_b <- b / 2
        final_log_tr_draws <- log_tr_draws - (b - shrunk_b)
      }
    } else { # Zwet or Sherry
      if (apply_to_or) final_log_or_draws <- sample_shrinkage_draws(log_or_draws, method = shrinkage_method)
      if (apply_to_tr) final_log_tr_draws <- sample_shrinkage_draws(log_tr_draws, method = shrinkage_method)
    }
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









