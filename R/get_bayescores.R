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

  # --- Step 1 & 2: No changes here ---
  defaults <- list(efficacy = list(cure_utility_target = list(effect_value = 0.15, utility_value = 40), tr_utility_target = list(effect_value = 1.25, utility_value = 35), tr_disutility_target = list(effect_value = 0.75, utility_value = 35), tr_baseline = 1.0))
  recursive_merge <- function(default, user) { for (name in names(user)) { if (is.list(user[[name]]) && is.list(default[[name]])) { default[[name]] <- recursive_merge(default[[name]], user[[name]]) } else { default[[name]] <- user[[name]] } }; return(default) }
  cal <- recursive_merge(defaults, calibration_args)
  cure_samples <- pmax(0, efficacy_inputs$cure_posterior_samples)
  tr_samples <- efficacy_inputs$tr_posterior_samples

  # --- Step 3 & 4: No changes here ---
  .create_exp_utility_fn <- function(target, baseline = 0) {
    effect_gain <- target$effect_value - baseline
    if (effect_gain <= 0) stop(paste("Target effect_value", target$effect_value, "must be greater than baseline", baseline))
    lambda <- -log(1 - target$utility_value / 100) / effect_gain
    return(function(x) 100 * (1 - exp(-lambda * pmax(0, x))))
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
  bonus_term <- (pmin(U_tr, U_cure) / 100) * (100 - pmax(U_tr, U_cure))
  penalty_term <- (DU_tr / 100) * U_cure
  score_post_eficacia <- pmax(0, pmax(U_tr, U_cure) + bonus_term - penalty_term)

  # --- Step 5: Incorporate QoL and Toxicity (CORRECTED LOGIC) ---
  # QoL adjustment remains the same
  qol_effect <- qol_scores / 2
  qol_power <- qol_effect * (score_post_eficacia / 100)
  score_post_qol <- score_post_eficacia + qol_power * (100 - score_post_eficacia)

  tox_effect <- toxicity_scores
  # -- CORRECTION: Toxicity power is now proportional to the base efficacy score --
  tox_power <- tox_effect * (score_post_eficacia / 100)

  # The final score is calculated sequentially, but tox_power's magnitude is based on efficacy
  final_utility_vector <- score_post_qol + tox_power * (100 - score_post_qol)

  # --- Step 6: Summarize Components (No changes here) ---
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
  return(list(final_utility_vector = final_utility_vector, component_summary = component_summary))
}
