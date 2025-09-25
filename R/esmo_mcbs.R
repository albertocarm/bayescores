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
           # CORRECTED: Specific check for missing or NULL arguments
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
           # CORRECTED: Specific check for missing or NULL arguments
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
           # CORRECTED: Specific check for missing or NULL arguments
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




# ---------------------------------------------------------------------------
# MASTER TEMPLATE FOR calculate_esmo_mcbs FUNCTION
#
# HOW TO USE:
# 1. Choose your 'objective' from the list in the first argument.
# 2. Find the corresponding section below and fill in the required data from your paper.
# 3. Add any 'Optional Adjustments' that apply to your study.
# 4. Leave all other arguments as NULL or FALSE. The function will throw an
#    error if a required argument for your chosen objective is missing.
# ---------------------------------------------------------------------------

# calculate_esmo_mcbs(
#
#   # Step 1: Choose ONE of the following objectives (mandatory)
#   #--------------------------------------------------------------------------
#   objective = "non_curative_os", # "non_curative_os", "non_curative_pfs", "adjuvant_rct",
#   # "single_arm_rare_disease", "non_inferiority_or_other",
#   # "de_escalation"
#
#   # Step 2: Fill in the arguments corresponding to your chosen objective
#   #--------------------------------------------------------------------------
#
#   # --- A) For Non-Curative RCTs (objective: "non_curative_os" or "non_curative_pfs") ---
#   median_control_survival_months = NULL, # Median OS/PFS in the control arm (in months)
#   hr = NULL,                             # Hazard Ratio for the primary endpoint (OS/PFS)
#   gain_months = NULL,                    # Absolute gain in median OS/PFS (in months)
#
#   # --- B) For Adjuvant/Curative RCTs (objective: "adjuvant_rct") ---
#   survival_gain_percent = NULL,          # Use this for absolute N-year survival gain (e.g., 5 for 5%)
#   dfs_gain_percent = NULL,               # OR use this for N-year DFS gain...
#   # hr = NULL,                           # ...with the Hazard Ratio for DFS. (hr is shared with section A)
#
#   # --- C) For Single-Arm, Rare Disease Trials (objective: "single_arm_rare_disease") ---
#   pfs_months = NULL,                     # Median PFS observed in the study
#   orr_percent = NULL,                    # Overall Response Rate observed (e.g., 45 for 45%)
#   dor_months = NULL,                     # Duration of Response (in months)
#
#   # --- D) For De-escalation Trials (objective: "de_escalation") ---
#   observed_vs_target_diff_percent = NULL,# % difference between observed and pre-specified goal
#
#   # --- E) For Non-Inferiority or Other Endpoint Trials (objective: "non_inferiority_or_other") ---
#   rr_increase_percent = NULL,            # Absolute increase in Response Rate
#
#   # Step 3: Add any of the following optional adjustments that apply
#   #--------------------------------------------------------------------------
#
#   # General Adjustments
#   qol_improved = FALSE,                  # Set to TRUE if Quality of Life is improved
#   fewer_toxicities_wellbeing = FALSE,    # Set to TRUE if there are fewer grade 3-4 toxicities
#   incremental_toxicity = FALSE,          # Set to TRUE if the new treatment has incremental toxicity
#
#   # Adjustments for specific objectives
#   survival_rate_increase = NULL,         # For "non_curative_os": Absolute N-year survival increase (e.g., 10 for 10%)
#   long_term_plateau_os = FALSE,          # For "non_curative_os": Set to TRUE if there is a long-term plateau
#   mature_os_shows_no_advantage = FALSE,  # For "non_curative_pfs": Set to TRUE if mature OS shows no benefit
#   early_crossover = FALSE,               # For "non_curative_pfs": Set to TRUE if early crossover occurred
#   mature_os_no_benefit = FALSE,          # For "adjuvant_rct": Set to TRUE if mature OS shows no benefit
#   high_grade_toxicities_wellbeing = FALSE,# For "single_arm_rare_disease": Set to TRUE if G3-4 toxicity >= 30%
#   is_confirmatory_phase4 = FALSE,        # For "single_arm_rare_disease": Set to TRUE for confirmatory phase 4 data
#   is_non_inferiority = FALSE,            # For "non_inferiority_or_other": Set to TRUE if it's a non-inferiority study
#   is_multicentre = TRUE                  # For "de_escalation": Defaults to TRUE, set to FALSE if single-centre
# )
