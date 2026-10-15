#' @title Calculate the ESMO-MCBS v2.0 Score (Guided by Objective)
#' @description Calculates the ESMO-MCBS v2.0 score based on a primary study
#' objective. The function first requires the user to specify the study's
#' main goal, and then uses the corresponding arguments to apply the
#' correct ESMO form logic.
#'
#' @param objective A character string specifying the main goal of the clinical
#'   trial. This is a **mandatory** argument.
#'
#' @param ... Additional arguments required for the chosen objective. See the
#'   documentation and examples for which arguments are needed for each objective.
#'
#' @return A character string with the final ESMO-MCBS score.

calculate_esmo_mcbs <- function(objective, ...) {

  args <- list(...)

  # Helper operator for setting default values if an argument is NULL
  `%||%` <- function(a, b) { if (is.null(a)) b else a }

  switch(objective,
         "non_curative_os" = {
           # Check for missing or NULL arguments
           required <- c("median_control_survival_months", "hr", "gain_months")
           missing_args <- required[sapply(required, function(p) !p %in% names(args) || is.null(args[[p]]))]
           if (length(missing_args) > 0) {
             stop(paste("For objective 'non_curative_os', the following required arguments are missing or NULL:",
                        paste(missing_args, collapse = ", ")))
           }

           qol <- args$qol_improved %||% FALSE
           tox <- args$fewer_toxicities_wellbeing %||% FALSE
           plat <- args$long_term_plateau_os %||% FALSE
           surv_inc <- args$survival_rate_increase

           if (args$median_control_survival_months < 12) {
             if ((args$hr <= 0.65 && args$gain_months >= 3) || (!is.null(surv_inc) && surv_inc >= 10)) { score <- 4 }
             else if (args$hr <= 0.65 && args$gain_months >= 2) { score <- 3 }
             else if ((args$hr <= 0.65 && args$gain_months >= 1.5) || (args$hr > 0.65 && args$hr <= 0.70 && args$gain_months >= 1.5)) { score <- 2 }
             else { score <- 1 }
           } else if (args$median_control_survival_months < 24) {
             if ((args$hr <= 0.70 && args$gain_months >= 5) || (!is.null(surv_inc) && surv_inc >= 10)) { score <- 4 }
             else if (args$hr <= 0.70 && args$gain_months >= 3) { score <- 3 }
             else if ((args$hr <= 0.70 && args$gain_months >= 1.5) || (args$hr > 0.70 && args$hr <= 0.75 && args$gain_months >= 1.5)) { score <- 2 }
             else { score <- 1 }
           } else if (args$median_control_survival_months < 36) {
             if ((args$hr <= 0.70 && args$gain_months >= 9) || (!is.null(surv_inc) && surv_inc >= 10)) { score <- 4 }
             else if (args$hr <= 0.70 && args$gain_months >= 6) { score <- 3 }
             else if ((args$hr <= 0.70 && args$gain_months >= 4) || (args$hr > 0.70 && args$hr <= 0.75 && args$gain_months >= 4)) { score <- 2 }
             else { score <- 1 }
           } else {
             if ((args$hr <= 0.70 && args$gain_months >= 12) || (!is.null(surv_inc) && surv_inc >= 10)) { score <- 4 }
             else if (args$hr <= 0.70 && args$gain_months >= 8) { score <- 3 }
             else if ((args$hr <= 0.70 && args$gain_months >= 6) || (args$hr > 0.70 && args$hr <= 0.75 && args$gain_months >= 6)) { score <- 2 }
             else { score <- 1 }
           }
           if (qol || tox) { score <- min(5, score + 1) }
           if (plat) { return(paste0("A/", score, " (Note: The plateau in the OS curve suggests also scoring with the objective 'adjuvant_rct'.)")) }
           return(as.character(score))
         },

         "non_curative_pfs" = {
           # Check for missing or NULL arguments
           required <- c("median_control_survival_months", "hr", "gain_months")
           missing_args <- required[sapply(required, function(p) !p %in% names(args) || is.null(args[[p]]))]
           if (length(missing_args) > 0) {
             stop(paste("For objective 'non_curative_pfs', the following required arguments are missing or NULL:",
                        paste(missing_args, collapse = ", ")))
           }

           qol <- args$qol_improved %||% FALSE
           tox <- args$fewer_toxicities_wellbeing %||% FALSE
           no_os <- args$mature_os_shows_no_advantage %||% FALSE
           inc_tox <- args$incremental_toxicity %||% FALSE
           cross <- args$early_crossover %||% FALSE

           if (args$median_control_survival_months < 6) { score <- ifelse(args$hr <= 0.65 && args$gain_months >= 1.5, 3, ifelse(args$hr <= 0.65, 2, 1)) }
           else if (args$median_control_survival_months < 12) { score <- ifelse(args$hr <= 0.65 && args$gain_months >= 3, 3, ifelse(args$hr <= 0.65, 2, 1)) }
           else { score <- ifelse(args$hr <= 0.65 && args$gain_months >= 5, 3, ifelse(args$hr <= 0.65, 2, 1)) }

           if (no_os && !qol) { score <- score - 1 }
           if (inc_tox) { score <- score - 1 }
           if (qol || tox || cross) { score <- score + 1 }

           return(as.character(max(1, min(4, score))))
         },

         "adjuvant_rct" = {
           if (is.null(args$survival_gain_percent) && (is.null(args$dfs_gain_percent) || is.null(args$hr))) {
             stop("For objective 'adjuvant_rct', provide 'survival_gain_percent' OR both 'dfs_gain_percent' and 'hr'.")
           }
           no_os <- args$mature_os_no_benefit %||% FALSE
           if (!is.null(args$survival_gain_percent)) { score <- ifelse(args$survival_gain_percent >= 5, "A", ifelse(args$survival_gain_percent > 3, "B", "C")) }
           else {
             if (args$hr <= 0.65 && args$dfs_gain_percent >= 3) score <- "A"
             else if ((args$hr <= 0.65 && args$dfs_gain_percent >= 1) || (args$hr > 0.65 && args$hr <= 0.75 && args$dfs_gain_percent >= 3)) score <- "B"
             else score <- "C"
           }
           if (no_os && score != "C") { score <- ifelse(score == "A", "B", "C") }
           return(score)
         },

         "single_arm_rare_disease" = {
           if (is.null(args$pfs_months) && is.null(args$orr_percent)) {
             stop("For objective 'single_arm_rare_disease', you must provide at least one of: 'pfs_months' or 'orr_percent'.")
           }
           is_grade_3 <- (!is.null(args$pfs_months) && args$pfs_months >= 6) || (!is.null(args$orr_percent) && args$orr_percent >= 60) || (!is.null(args$orr_percent) && args$orr_percent >= 20 && !is.null(args$dor_months) && args$dor_months >= 9)
           is_grade_2 <- (!is.null(args$pfs_months) && args$pfs_months >= 3 && args$pfs_months < 6) || (!is.null(args$orr_percent) && args$orr_percent >= 40 && args$orr_percent < 60) || (!is.null(args$orr_percent) && args$orr_percent >= 20 && args$orr_percent < 40 && !is.null(args$dor_months) && args$dor_months >= 6 && args$dor_months < 9)
           if (is_grade_3) { score <- 3 } else if (is_grade_2) { score <- 2 } else { score <- 1 }
           if (args$high_grade_toxicities_wellbeing %||% FALSE) { score <- score - 1 }
           if ((args$qol_improved %||% FALSE) || (args$is_confirmatory_phase4 %||% FALSE)) { score <- score + 1 }
           return(as.character(max(1, min(4, score))))
         },

         "non_inferiority_or_other" = {
           non_inf <- args$is_non_inferiority %||% FALSE
           qol <- args$qol_improved %||% FALSE
           fewer_tox <- args$fewer_toxicities_wellbeing %||% FALSE
           rr_inc <- args$rr_increase_percent
           inc_tox <- args$incremental_toxicity %||% FALSE

           if (non_inf && (qol || fewer_tox)) { score <- 4 }
           else if (qol) { score <- 3 }
           else if (!is.null(rr_inc) && rr_inc >= 20) { score <- 2 }
           else { score <- 1 }
           if (qol || fewer_tox) { score <- min(4, score + 1) }
           if (inc_tox) { score <- max(1, score - 1) }
           return(as.character(score))
         },

         "de_escalation" = {
           # Check for missing or NULL arguments
           if (is.null(args$observed_vs_target_diff_percent)) {
             stop("For objective 'de_escalation', the following required argument is missing or NULL: observed_vs_target_diff_percent")
           }
           multi <- args$is_multicentre %||% TRUE
           if (multi && args$observed_vs_target_diff_percent <= 2) { score <- "A" }
           else if (!multi && args$observed_vs_target_diff_percent <= 2) { score <- "B" }
           else { score <- "C" }
           return(score)
         },

         {
           valid_objectives <- c("non_curative_os", "non_curative_pfs", "adjuvant_rct",
                                 "single_arm_rare_disease", "non_inferiority_or_other",
                                 "de_escalation")
           stop(paste("Invalid 'objective' specified. Please use one of the following options:\n-",
                      paste(valid_objectives, collapse = "\n- ")))
         }
  )
}
