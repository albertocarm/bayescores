
#' @title Calculate BayeScores Final Utility Score
#' @description
#' Combines the posterior distributions of efficacy (cure fraction and
#' tail-risk), quality of life (QoL), and toxicity into a single final
#' utility score (0-100) using exponential utility functions.
#'
#' @param efficacy_inputs A list containing the posterior samples for efficacy.
#'   Must include `cure_posterior_samples` and `tr_posterior_samples` as
#'   numeric vectors.
#' @param qol_scores A numeric vector of posterior samples for the Quality of
#'   Life effect.
#' @param toxicity_scores A numeric vector of posterior samples for the
#'   Toxicity effect.
#' @param calibration_args An optional list with calibration parameters to
#'   override the model's default values.
#' @return A list containing:
#'   \itemize{
#'     \item `final_utility_vector`: The posterior distribution vector of the final utility.
#'     \item `component_summary`: A data frame with the summary (median and 95% CrI)
#'       of each utility component.
#'   }
#' @export
#' @importFrom stats median quantile
#' @examples
#' \dontrun{
#'   efficacy_data <- list(
#'     cure_posterior_samples = rbeta(1000, 20, 80),
#'     tr_posterior_samples = rnorm(1000, 1.1, 0.2)
#'   )
#'   qol_data <- rnorm(1000, 5, 2)
#'   toxicity_data <- rnorm(1000, -3, 1.5)
#'
#'   results <- get_bayescores(
#'     efficacy_inputs = efficacy_data,
#'     qol_scores = qol_data,
#'     toxicity_scores = toxicity_data
#'   )
#'   print(results$component_summary)
#' }

get_bayescores <- function(
    efficacy_inputs, qol_scores, toxicity_scores, calibration_args = list()
) {
  defaults <- list(efficacy = list(cure_utility_target = list(effect_value = 0.15, utility_value = 40), tr_utility_target = list(effect_value = 1.25, utility_value = 35), tr_baseline = 1.0))
  recursive_merge <- function(default, user) { for (name in names(user)) { if (is.list(user[[name]]) && is.list(default[[name]])) { default[[name]] <- recursive_merge(default[[name]], user[[name]]) } else { default[[name]] <- user[[name]] } }; return(default) }
  cal <- recursive_merge(defaults, calibration_args)
  cure_samples <- pmax(0, efficacy_inputs$cure_posterior_samples)
  tr_samples <- efficacy_inputs$tr_posterior_samples
  .create_exp_utility_fn <- function(target, baseline = 0) {
    effect_gain <- target$effect_value - baseline
    if (effect_gain <= 0) stop(paste("Target effect_value", target$effect_value, "must be greater than baseline", baseline))
    lambda <- -log(1 - target$utility_value / 100) / effect_gain
    return(function(x) 100 * (1 - exp(-lambda * x)))
  }
  fn_u_cure <- .create_exp_utility_fn(cal$efficacy$cure_utility_target, 0)
  fn_u_tr <- .create_exp_utility_fn(cal$efficacy$tr_utility_target, cal$efficacy$tr_baseline)
  U_cure <- fn_u_cure(cure_samples); U_tr <- fn_u_tr(pmax(0, tr_samples - cal$efficacy$tr_baseline))
  score_post_eficacia <- pmax(U_tr, U_cure) + (pmin(U_tr, U_cure) / 100) * (100 - pmax(U_tr, U_cure))
  qol_effect <- qol_scores / 2
  qol_power <- qol_effect * (score_post_eficacia / 100)
  score_post_qol <- score_post_eficacia + qol_power * (100 - score_post_eficacia)
  tox_effect <- toxicity_scores
  tox_power <- tox_effect * (score_post_qol / 100)
  final_utility_vector <- score_post_qol + tox_power * (100 - score_post_qol)
  .summarize <- function(vec) { c(Median = median(vec), Lower_95_CrI=quantile(vec, .025), Upper_95_CrI=quantile(vec, .975)) }
  summary_list <- list("Utility TR (0-100)" = .summarize(U_tr), "Utility Cure (0-100)" = .summarize(U_cure), "Efficacy Score (Combined)" = .summarize(score_post_eficacia), "QoL Contribution (points)" = .summarize(qol_power * (100 - score_post_eficacia)), "Toxicity Contribution (points)" = .summarize(tox_power * (100 - score_post_qol)), "FINAL UTILITY SCORE" = .summarize(final_utility_vector))
  component_summary <- as.data.frame(do.call(rbind, summary_list)); component_summary <- cbind(Component = rownames(component_summary), component_summary)
  return(list(final_utility_vector = final_utility_vector, component_summary = component_summary))
}
