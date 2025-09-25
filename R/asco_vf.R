#' @title Calculate the ASCO Net Health Benefit (NHB) Score
#' @description This function calculates the Net Health Benefit (NHB) based on the American
#' Society of Clinical Oncology (ASCO) Value Framework. It differentiates
#' between two clinical settings: "advanced" (metastatic) and "adjuvant"
#' (curative intent), as their scoring criteria differ significantly.
#'
#' @param setting A string specifying the clinical context. Must be either "advanced" or "adjuvant".
#' @param os_hr Numeric. The Hazard Ratio for Overall Survival.
#' @param os_median_control Numeric. Median Overall Survival for the control arm (in months).
#' @param os_median_treatment Numeric. Median Overall Survival for the treatment arm (in months).
#' @param pfs_hr Numeric. The Hazard Ratio for Progression-Free Survival.
#' @param pfs_median_control Numeric. Median Progression-Free Survival for the control arm (in months).
#' @param pfs_median_treatment Numeric. Median Progression-Free Survival for the treatment arm (in months).
#' @param orr_treatment Numeric. The Objective Response Rate for the treatment arm (as a decimal, e.g., 0.40 for 40\%).
#' @param dfs_hr Numeric. The Hazard Ratio for Disease-Free Survival (used in 'adjuvant' setting).
#' @param total_tox_points_control A numeric value for the sum of toxicity points for the control/standard-of-care arm.
#' @param total_tox_points_treatment A numeric value for the sum of toxicity points for the new/experimental arm.
#' @param bonus_tail_of_curve Logical (TRUE/FALSE). Was a significant long-term survival benefit ("tail of the curve") demonstrated?
#' @param bonus_qol Logical (TRUE/FALSE). Was a statistically significant improvement in Quality of Life reported?
#' @param bonus_palliation Logical (TRUE/FALSE). Was significant symptom palliation shown? (for 'advanced' setting).
#' @param bonus_treatment_free_interval_points Numeric. Points for treatment-free interval (for 'advanced' setting).
#'
#' @return A named list with four elements: \code{clinical_benefit_score}, \code{toxicity_adjustment},
#' \code{bonus_score}, and \code{net_health_benefit}.
#'
#' @export
#'
#' @examples
#' # Example 1: Advanced Disease (Metastatic)
#' # Hypothetical trial where the primary endpoint is OS improvement shown by HR.
#' calculate_asco_nhb(
#'   setting = "advanced",
#'   os_hr = 0.75,
#'   total_tox_points_control = 50,
#'   total_tox_points_treatment = 60,
#'   bonus_qol = TRUE,
#'   bonus_palliation = FALSE
#' )
#'
#' # Example 2: Adjuvant Setting (Curative Intent)
#' # Based on the SPARTAN trial for apalutamide.
#' calculate_asco_nhb(
#'   setting = "adjuvant",
#'   os_hr = 0.78,
#'   total_tox_points_control = 37.0,
#'   total_tox_points_treatment = 43.7,
#'   bonus_tail_of_curve = TRUE,
#'   bonus_qol = FALSE
#' )



calculate_asco_nhb <- function(
    setting,
    # Efficacy arguments
    os_hr = NA,
    os_median_control = NA,
    os_median_treatment = NA,
    pfs_hr = NA,
    pfs_median_control = NA,
    pfs_median_treatment = NA,
    orr_treatment = NA,
    dfs_hr = NA,
    # Toxicity arguments
    total_tox_points_control = 0,
    total_tox_points_treatment = 0,
    # Bonus points arguments
    bonus_tail_of_curve = FALSE,
    bonus_qol = FALSE,
    bonus_palliation = FALSE,
    bonus_treatment_free_interval_points = 0
) {

  # --- Input Validation ---
  if (missing(setting) || !setting %in% c("advanced", "adjuvant")) {
    stop("Invalid 'setting'. You must specify either 'advanced' or 'adjuvant'. This is crucial because the two settings use different clinical endpoints (e.g., OS/PFS vs. DFS), scoring hierarchies, and maximum possible scores according to the ASCO Framework.")
  }

  clinical_benefit_score <- 0

  # --- I. Calculate Clinical Benefit Score based on Setting ---
  if (setting == "advanced") {
    # Check for missing paired values
    if (!is.na(os_median_control) && is.na(os_median_treatment)) stop("Error: 'os_median_control' was provided, but 'os_median_treatment' is missing.")
    if (is.na(os_median_control) && !is.na(os_median_treatment)) stop("Error: 'os_median_treatment' was provided, but 'os_median_control' is missing.")
    if (!is.na(pfs_median_control) && is.na(pfs_median_treatment)) stop("Error: 'pfs_median_control' was provided, but 'pfs_median_treatment' is missing.")
    if (is.na(pfs_median_control) && !is.na(pfs_median_treatment)) stop("Error: 'pfs_median_treatment' was provided, but 'pfs_median_control' is missing.")

    # Hierarchy for Advanced Disease
    if (!is.na(os_hr)) {
      clinical_benefit_score <- (1 - os_hr) * 100
    } else if (!is.na(os_median_control) && !is.na(os_median_treatment)) {
      improvement_pct <- ((os_median_treatment - os_median_control) / os_median_control) * 100
      clinical_benefit_score <- improvement_pct
    } else if (!is.na(pfs_hr)) {
      clinical_benefit_score <- (1 - pfs_hr) * 100 * 0.8
    } else if (!is.na(pfs_median_control) && !is.na(pfs_median_treatment)) {
      improvement_pct <- ((pfs_median_treatment - pfs_median_control) / pfs_median_control) * 100
      clinical_benefit_score <- improvement_pct * 0.8
    } else if (!is.na(orr_treatment)) {
      clinical_benefit_score <- orr_treatment * 100 * 0.7
    } else {
      stop("Error: No valid efficacy data provided for 'advanced' setting. You must provide at least one of the following: os_hr, os_median values, pfs_hr, pfs_median values, or orr_treatment.")
    }

  } else if (setting == "adjuvant") {
    # Hierarchy for Adjuvant Disease
    if (!is.na(os_hr)) {
      clinical_benefit_score <- (1 - os_hr) * 100
    } else if (!is.na(dfs_hr)) {
      clinical_benefit_score <- (1 - dfs_hr) * 100
    } else {
      stop("Error: No valid efficacy data provided for 'adjuvant' setting. You must provide either 'os_hr' or 'dfs_hr'.")
    }
  }

  # --- II. Calculate Toxicity Adjustment ---
  toxicity_adjustment <- 0
  if (total_tox_points_control > 0) {
    # Formula: ((Control - Treatment) / Control) * 20
    # If treatment is more toxic, this results in a negative adjustment.
    toxicity_adjustment <- ((total_tox_points_control - total_tox_points_treatment) / total_tox_points_control) * 20
  }

  # --- III. Calculate Bonus Points ---
  bonus_score <- 0
  if (bonus_tail_of_curve) bonus_score <- bonus_score + 20
  if (bonus_qol) bonus_score <- bonus_score + 10

  # Palliation and TFI points are specific to advanced setting
  if (setting == "advanced") {
    if (bonus_palliation) bonus_score <- bonus_score + 10
    if (bonus_treatment_free_interval_points > 0) {
      bonus_score <- bonus_score + bonus_treatment_free_interval_points
    }
  }

  # --- IV. Calculate Final Net Health Benefit ---
  net_health_benefit <- clinical_benefit_score + toxicity_adjustment + bonus_score

  # --- V. Return Results ---
  return(
    list(
      clinical_benefit_score = round(clinical_benefit_score, 2),
      toxicity_adjustment = round(toxicity_adjustment, 2),
      bonus_score = round(bonus_score, 2),
      net_health_benefit = round(net_health_benefit, 2)
    )
  )
}
