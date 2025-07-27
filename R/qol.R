#' Interactively Generate a Quality of Life (QoL) Probability Vector
#'
#' @description
#' This is an interactive helper function that guides a user through a series of
#' questions to generate a probability vector for Quality of Life (QoL) outcomes.
#' The outcomes are discretized into five levels, from significant deterioration (-2)
#' to significant improvement (+2).
#'
#' The final vector is calculated by blending a predefined "scenario vector" with a
#' uniform "uncertainty vector". The weight of this blend is determined by the
#' user's stated confidence in the evidence, creating a final probability
#' distribution that reflects both the expected outcome and its certainty.
#'
#' @details
#' The function performs the following steps:
#' 1.  Prompts the user to select one of five predefined QoL scenarios (e.g., "Significant Improvement").
#' 2.  Prompts the user to rate the strength of evidence for that scenario on a 5-point scale.
#' 3.  Calculates a confidence factor `alpha` (from 0 to 1) based on the evidence strength.
#' 4.  Computes the final probability vector using the formula:
#'     `final_vector = alpha * scenario_vector + (1 - alpha) * uncertainty_vector`.
#'
#' @return A named numeric vector of length 5. The names (`P(-2)` to `P(+2)`)
#'   represent the QoL change level, and the values are the corresponding
#'   probabilities, which sum to 1.
#'
#' @export
#'
#' @examples
#' # The function is interactive, so the call should be wrapped in \dontrun{}
#' # for automated checking.
#' \dontrun{
#'
#' # Run the function
#' qol_vector <- generate_qol_vector()
#'
#' # The console will prompt you to enter information. For example:
#' #
#' # Select the main Quality of Life (QoL) scenario:
#' # 1: Significant Improvement
#' # 2: Stabilization / Probable Benefit
#' # 3: No Difference / Marginal Benefit
#' # 4: Deterioration
#' # 5: Insufficient Data / Unknown
#' # Enter a number (1-5): 2
#' #
#' # Select the strength/quality of the evidence for this scenario:
#' # 1: Very Low (almost speculative)
#' # 2: Low
#' # 3: Moderate (standard)
#' # 4: High
#' # 5: Very High (e.g., large RCT, low bias, p<0.001)
#' # Enter a number (1-5): 4
#'
#' # The resulting vector would be stored in `qol_vector`
#' # print(qol_vector)
#' # P(-2)  P(-1)   P(0)  P(+1)  P(+2)
#' # 0.0875 0.1250 0.2000 0.4250 0.1625
#'
#' }
generate_qol_vector <- function() {

  # -------------------------------------------------------------------------- #
  # 1. DEFINE BASE PROBABILITY VECTORS
  # -------------------------------------------------------------------------- #

  # These vectors are predefined probability distributions (archetypes) for
  # different clinical outcomes regarding Quality of Life.
  # The order of probabilities corresponds to outcomes: [-2, -1, 0, +1, +2].
  base_vectors <- list(
    "1" = c(0.05, 0.10, 0.15, 0.30, 0.40), # Scenario 1: Significant Improvement (mode is +1/+2)
    "2" = c(0.05, 0.10, 0.20, 0.50, 0.15), # Scenario 2: Stabilization (mode is +1)
    "3" = c(0.10, 0.15, 0.50, 0.15, 0.10), # Scenario 3: No Difference (mode is 0)
    "4" = c(0.40, 0.30, 0.15, 0.10, 0.05), # Scenario 4: Deterioration (mode is -2/-1)
    "5" = c(0.20, 0.20, 0.20, 0.20, 0.20)  # Scenario 5: Insufficient Data (uniform/max entropy)
  )

  # -------------------------------------------------------------------------- #
  # 2. PROMPT USER FOR THE MAIN QoL SCENARIO
  # -------------------------------------------------------------------------- #

  # Display the available options to the user.
  cat("Select the main Quality of Life (QoL) scenario:\n")
  cat("1: Significant Improvement\n")
  cat("2: Stabilization / Probable Benefit\n")
  cat("3: No Difference / Marginal Benefit\n")
  cat("4: Deterioration\n")
  cat("5: Insufficient Data / Unknown\n")

  # Capture user input from the console.
  scenario_num <- readline(prompt = "Enter a number (1-5): ")

  # Validate the input to ensure it's one of the defined scenarios.
  if (!(scenario_num %in% names(base_vectors))) {
    stop("Error: Invalid scenario number. Please enter a number from 1 to 5.")
  }

  # Retrieve the corresponding base vector based on user selection.
  selected_base_vector <- base_vectors[[scenario_num]]

  # -------------------------------------------------------------------------- #
  # 3. PROMPT USER FOR THE STRENGTH OF EVIDENCE
  # -------------------------------------------------------------------------- #

  # Display the evidence strength options. This will determine confidence.
  cat("\nSelect the strength/quality of the evidence for this scenario:\n")
  cat("1: Very Low (almost speculative)\n")
  cat("2: Low\n")
  cat("3: Moderate (standard)\n")
  cat("4: High\n")
  cat("5: Very High (e.g., large RCT, low bias, p<0.001)\n")

  # Capture and convert the user's input to an integer.
  strength_num <- as.integer(readline(prompt = "Enter a number (1-5): "))

  # Validate the input to ensure it's a number between 1 and 5.
  if (is.na(strength_num) || !(strength_num %in% 1:5)) {
    stop("Error: Invalid strength number. Please enter a number from 1 to 5.")
  }

  # -------------------------------------------------------------------------- #
  # 4. CALCULATE THE FINAL PROBABILITY VECTOR
  # -------------------------------------------------------------------------- #

  # Convert the 5-point strength scale into a "confidence factor" (alpha) from 0 to 1.
  # Strength 1 (Very Low)   -> alpha = (1-1)/4 = 0   (Total uncertainty)
  # Strength 3 (Moderate) -> alpha = (3-1)/4 = 0.5 (50% confidence)
  # Strength 5 (Very High)  -> alpha = (5-1)/4 = 1   (Full confidence in the scenario)
  alpha <- (strength_num - 1) / 4

  # Define the vector for total uncertainty (a uniform probability distribution).
  uniform_vector <- rep(0.2, 5)

  # The final vector is a weighted average (linear interpolation) between the
  # selected scenario vector and the uncertainty vector, weighted by alpha.
  final_vector <- alpha * selected_base_vector + (1 - alpha) * uniform_vector

  # Assign names to the elements of the vector for better readability and interpretation.
  names(final_vector) <- c("P(-2)", "P(-1)", "P(0)", "P(+1)", "P(+2)")

  # Return the final calculated vector.
  return(final_vector)
}



#' Sample Quality of Life (QoL) Scores from a Multinomial Distribution
#'
#' @description
#' This function generates a vector of random Quality of Life (QoL) scores by
#' sampling from a discrete probability distribution. It simulates outcomes
#' where each sample has a specific probability of falling into one of five
#' predefined QoL change categories.
#'
#' @details
#' This function is a robust wrapper around the base R `sample()` function. It's
#' specifically configured for sampling with replacement from the discrete set
#' of scores `[-2, -1, 0, 1, 2]`. It includes validation to ensure the input
#' probability vector is correctly formatted. This is useful for simulating
#' patient cohorts in health economic models or clinical trial simulations.
#'
#' @param prob_vector A numeric vector containing exactly 5 probabilities, which
#'   **must sum to 1**. Each element corresponds to the probability of the QoL
#'   scores `[-2, -1, 0, +1, +2]`, respectively.
#' @param n_samples An integer specifying the total number of samples to generate
#'   (i.e., the length of the output vector).
#'
#' @return A numeric vector of length `n_samples` containing the sampled QoL
#'   scores. Each element in the vector will be one of the values from
#'   `[-2, -1, 0, 1, 2]`.
#'
#' @export
#'
#' @examples
#' # --- 1. Define a probability vector ---
#' # This vector represents a "No Difference / Marginal Benefit" scenario.
#' # P(-2)=0.1, P(-1)=0.15, P(0)=0.5, P(+1)=0.15, P(+2)=0.1
#' my_probs <- c(0.1, 0.15, 0.5, 0.15, 0.1)
#'
#' # --- 2. Set the number of samples ---
#' # For example, simulate a cohort of 1000 patients.
#' num_patients <- 1000
#'
#' # --- 3. Generate the sampled scores ---
#' qol_draws <- sample_qol_scores(prob_vector = my_probs, n_samples = num_patients)
#'
#' # --- 4. View the results ---
#' # Look at the first few sampled scores
#' head(qol_draws)
#'
#' # Check the distribution of the generated samples.
#' # It should approximate the input probabilities.
#' table(qol_draws) / num_patients
#'
sample_qol_scores <- function(prob_vector, n_samples) {

  # -------------------------------------------------------------------------- #
  # 1. DEFINE MODEL PARAMETERS
  # -------------------------------------------------------------------------- #

  # Define the five possible discrete QoL outcome scores. These are fixed values
  # representing the different levels of change (e.g., -2 is significant deterioration).
  scores <- c(-2, -1, 0, 1, 2)

  # -------------------------------------------------------------------------- #
  # 2. INPUT VALIDATION
  # -------------------------------------------------------------------------- #

  # It's crucial to validate inputs to prevent errors and ensure correct behavior.
  # Check 1: Ensure the probability vector has exactly 5 elements, one for each score.
  if (length(prob_vector) != 5) {
    stop("Error: 'prob_vector' must contain exactly 5 probabilities.")
  }

  # Check 2: Ensure the probabilities sum to 1.
  # Use all.equal() for safe comparison of floating-point numbers.
  if (!isTRUE(all.equal(sum(prob_vector), 1))) {
    stop("Error: The probabilities in 'prob_vector' must sum to 1.")
  }

  # -------------------------------------------------------------------------- #
  # 3. PERFORM RANDOM SAMPLING
  # -------------------------------------------------------------------------- #

  # Use the base R sample() function to efficiently generate the random vector.
  # This function draws random samples from a given set of values based on
  # a vector of probabilities.
  sampled_vector <- sample(
    x       = scores,       # The set of values to sample from.
    size    = n_samples,    # The total number of items to choose.
    replace = TRUE,         # Set to TRUE to allow values to be chosen more than once.
    prob    = prob_vector   # A vector of weights/probabilities for each value in 'x'.
  )

  # -------------------------------------------------------------------------- #
  # 4. RETURN RESULT
  # -------------------------------------------------------------------------- #

  # Return the final vector containing the randomly sampled QoL scores.
  return(sampled_vector)
}




#' Summarizes the final utility score with advanced dynamic weighting.
#'
#' @description
#' This function calculates the final weighted utility score. It incorporates a
#' sophisticated weighting policy that adapts based on the detected efficacy
#' profile ('Cure' vs. 'Survival').
#'
#' In the 'Survival' profile, the weight of the 'cure' component is scaled by its
#' observed magnitude, preventing a weak signal from diluting the score.
#'
#' In all cases, the weight of the 'QoL' component is also scaled by its
#' magnitude, preventing a neutral QoL from diluting the score.
#'
#' Any "unused" weight from these gradual adjustments is redistributed
#' proportionally among the other active components.
#'
#' @param efficacy_scores A list with 'tr_posterior_samples' and 'cure_posterior_samples'.
#' @param toxicity_scores A numeric vector of toxicity adjustment scores (-1 to +1).
#' @param qol_scores A numeric vector of sampled QoL scores (-2 to +2).
#' @param cure_benefit_threshold Numeric. The threshold above which the cure fraction is considered beneficial. Default is 0.02.
#' @param prob_certainty_threshold Numeric. The probability threshold (e.g., 0.95) to determine if a benefit is "certain".
#' @param tr_min_if_cure Numeric. The minimum Time Ratio (TR) considered beneficial if a cure is present.
#' @param tr_min_if_no_cure Numeric. The minimum Time Ratio (TR) considered beneficial if no cure is present.
#' @param tr_hopeful Numeric. The Time Ratio threshold for a "hopeful" or highly desirable outcome.
#' @param cure_min_relevant Numeric. The minimum cure fraction considered clinically relevant.
#' @param cure_hopeful Numeric. The cure fraction threshold for a "hopeful" or highly desirable outcome.
#' @param tr_calibration_k Numeric. The calibration exponent (k) applied to the Time Ratio (TR) utility scores. Default is 0.35.
#' @param weights_if_cure A named numeric vector of weights for `c(tr, cure, tox, qol)` if a "Cure" profile is detected.
#' @param weights_if_no_cure A named numeric vector of weights for `c(tr, cure, tox, qol)` if a "Survival" (no cure) profile is detected.
#'
#' @return A list containing the final utility vector (with floor at 0) and a component summary.
#' @importFrom stats median quantile
#' @export
#'

summarize_final_utility <- function(
    efficacy_scores,
    toxicity_scores,
    qol_scores,
    cure_benefit_threshold = 0.02,
    prob_certainty_threshold = 0.95,
    tr_min_if_cure = 1.0,
    tr_min_if_no_cure = 1.15,
    tr_hopeful = 1.25,
    cure_min_relevant = 0.04,
    cure_hopeful = 0.12,
    tr_calibration_k = 0.35,
    weights_if_cure = c(tr = 0.15, cure = 0.60, tox = 0.15, qol = 0.10),
    weights_if_no_cure = c(tr = 0.7, cure = 0.10, tox = 0.10, qol = 0.10)
) {

  # ========================================================================== #
  # STEP A & B: EFFICACY PROFILE DETECTION
  # -------------------------------------------------------------------------- #
  # First, we determine if the evidence for a cure fraction is strong enough
  # to be considered a 'Cure' profile, or if it's more likely a 'Survival' profile.
  # This decision dictates the primary weighting policy.
  # ========================================================================== #
  if (!all(c("tr_posterior_samples", "cure_posterior_samples") %in% names(efficacy_scores))) {
    stop("The 'efficacy_scores' list must contain 'tr_posterior_samples' and 'cure_posterior_samples'.")
  }
  tr_posterior_samples <- efficacy_scores$tr_posterior_samples
  cure_posterior_samples <- efficacy_scores$cure_posterior_samples
  prob_benefit <- mean(cure_posterior_samples > cure_benefit_threshold)

  if (prob_benefit > prob_certainty_threshold) {
    profile <- "Cure"
    tr_min_relevant_to_use <- tr_min_if_cure
  } else {
    profile <- "Survival"
    tr_min_relevant_to_use <- tr_min_if_no_cure
  }
  cat(sprintf("Info: P(Cure Benefit > %.2f) = %.2f. Profile detected: '%s'.\n",
              cure_benefit_threshold, prob_benefit, profile))

  # ========================================================================== #
  # STEP C: CONVERT RAW EFFICACY TO 0-100 UTILITY SCORES
  # -------------------------------------------------------------------------- #
  # We use logistic functions to convert the raw efficacy metrics (Time Ratio
  # and Cure Fraction) into standardized 0-100 utility scores. This ensures
  # that a clinically meaningful difference maps to a relevant utility score.
  # ========================================================================== #
  solve_logistic_params <- function(x1, y1, x2, y2) {
    logit <- function(p) { log(p / (1 - p)) }; p1 <- y1 / 100; p2 <- y2 / 100
    logit_p1 <- logit(p1); logit_p2 <- logit(p2); s <- (x2 - x1) / (logit_p2 - logit_p1)
    x0 <- x1 - s * logit_p1; return(list(location = x0, scale = s))
  }
  utility_min <- 10;

  # Convert Time Ratio (TR) scores
  tr_params <- solve_logistic_params(tr_min_relevant_to_use, utility_min, tr_hopeful, 38)
  tr_scores_original <- 100 * stats::plogis(tr_posterior_samples, location = tr_params$location, scale = tr_params$scale)
  tr_scores <- 100 * (tr_scores_original / 100)^tr_calibration_k

  # Convert Cure Fraction scores
  cure_params <- solve_logistic_params(cure_min_relevant, 10, cure_hopeful, 38)
  cure_scores <- 100 * stats::plogis(cure_posterior_samples, location = cure_params$location, scale = cure_params$scale)

  # ========================================================================== #
  # STEP D: DYNAMIC WEIGHTING BASED ON PROFILE AND MAGNITUDE
  # -------------------------------------------------------------------------- #
  # This is the core logic. First, select the base weights based on the
  # detected profile. Then, dynamically adjust those weights based on the
  # observed magnitude of the Cure (if profile is 'Survival') and QoL effects.
  # ========================================================================== #

  # --- D.1: Select base weights and apply profile-specific logic ---
  weights_to_use <- if (profile == "Cure") weights_if_cure else weights_if_no_cure

  if (profile == "Survival") {
    # Gradual Weight Transfer for the CURE component
    # Instead of an arbitrary bonus 'k', we scale the cure component's weight
    # by its observed magnitude. If the cure signal is negligible, its weight
    # becomes negligible, and is redistributed to other components.

    # 1. Calculate the magnitude of the cure effect (0 to 1).
    # We scale it relative to the 'hopeful' cure threshold.
    median_cure_fraction <- stats::median(cure_posterior_samples)
    cure_magnitude_factor <- pmin(1, median_cure_fraction / cure_hopeful)

    # 2. Calculate the 'actual' cure weight and the 'unused' portion to redistribute.
    original_cure_weight <- weights_to_use["cure"]
    actual_cure_weight <- original_cure_weight * cure_magnitude_factor
    unused_cure_weight <- original_cure_weight - actual_cure_weight

    # 3. Redistribute the unused weight proportionally to TR, Tox, and QoL.
    sum_other_weights <- 1 - original_cure_weight
    if (sum_other_weights > 0) {
      new_weights <- weights_to_use
      new_weights["cure"] <- actual_cure_weight
      new_weights["tr"] <- weights_to_use["tr"] + (unused_cure_weight * (weights_to_use["tr"] / sum_other_weights))
      new_weights["tox"] <- weights_to_use["tox"] + (unused_cure_weight * (weights_to_use["tox"] / sum_other_weights))
      new_weights["qol"] <- weights_to_use["qol"] + (unused_cure_weight * (weights_to_use["qol"] / sum_other_weights))
      weights_to_use <- new_weights
    }

    cat(sprintf("Info: 'Survival' profile. Cure magnitude (vs hopeful) is %.2f. Applying gradual weighting.\n", cure_magnitude_factor))
    # In this new logic, final_cure_scores are simply the original scores.
    final_cure_scores <- cure_scores
    # END

  } else { # Profile is "Cure"
    # This logic remains: when a cure is detected, we adjust the TR/Cure balance
    # based on how dominant the cure signal is.
    median_cure_diff <- stats::median(cure_posterior_samples)
    cure_dominance_factor <- stats::plogis(median_cure_diff, location = 0.10, scale = 0.1)
    total_efficacy_weight <- weights_to_use["tr"] + weights_to_use["cure"]
    weights_to_use["cure"] <- total_efficacy_weight * cure_dominance_factor
    weights_to_use["tr"] <- total_efficacy_weight * (1 - cure_dominance_factor)
    final_cure_scores <- cure_scores
    cat(sprintf("Info: 'Cure' profile. Median cure diff=%.3f -> Dominance Factor=%.2f\n",
                median_cure_diff, cure_dominance_factor))
  }

  # --- D.2: Apply Gradual Weight Transfer for the QoL component ---
  # This logic is now applied sequentially AFTER the cure weight adjustments.
  median_qol_score <- stats::median(qol_scores)
  qol_magnitude_factor <- pmin(1, abs(median_qol_score) / 2)

  original_qol_weight <- weights_to_use["qol"]
  actual_qol_weight <- original_qol_weight * qol_magnitude_factor
  unused_qol_weight <- original_qol_weight - actual_qol_weight

  sum_other_weights_for_qol <- 1 - original_qol_weight
  if (sum_other_weights_for_qol > 0 && unused_qol_weight > 0) {
    final_weights <- weights_to_use
    final_weights["qol"] <- actual_qol_weight
    final_weights["tr"] <- weights_to_use["tr"] + (unused_qol_weight * (weights_to_use["tr"] / sum_other_weights_for_qol))
    final_weights["cure"] <- weights_to_use["cure"] + (unused_qol_weight * (weights_to_use["cure"] / sum_other_weights_for_qol))
    final_weights["tox"] <- weights_to_use["tox"] + (unused_qol_weight * (weights_to_use["tox"] / sum_other_weights_for_qol))
    weights_to_use <- final_weights
  }

  cat(sprintf("Info: QoL magnitude is %.2f. Applying gradual weighting.\n", qol_magnitude_factor))
  cat("Info: Final adjusted weights ->", paste(names(weights_to_use), round(weights_to_use, 3), collapse=", "), "\n")

  # ========================================================================== #
  # STEP E: CALCULATE FINAL UTILITY AND SUMMARIZE
  # -------------------------------------------------------------------------- #
  toxicity_contribution <- weights_to_use["tox"] * (toxicity_scores * 100)
  qol_contribution <- weights_to_use["qol"] * (qol_scores * 50)
  tr_contribution <- weights_to_use["tr"] * tr_scores
  cure_contribution <- weights_to_use["cure"] * final_cure_scores

  final_utility_vector <- tr_contribution + cure_contribution + toxicity_contribution + qol_contribution
  final_utility_vector <- pmax(0, final_utility_vector)

  # ... (The summary_df and return code does not change) ...
  summary_df <- data.frame(
    Component = c("1. -> Score TR (Calibrated)", "2. -> Score Cure (0-100)", "3. -> Score Toxicity (scaled)", "4. -> Score QoL (scaled)",
                  "Contribution TR (Weighted)", "Contribution Cure (Weighted)", "Contribution Toxicity (Weighted)", "Contribution QoL (Weighted)",
                  "FINAL UTILITY SCORE"),
    Median = c(stats::median(tr_scores), stats::median(cure_scores), stats::median(toxicity_scores * 100), stats::median(qol_scores * 50),
               stats::median(tr_contribution), stats::median(cure_contribution), stats::median(toxicity_contribution), stats::median(qol_contribution),
               stats::median(final_utility_vector)),
    Lower_95_CrI = c(stats::quantile(tr_scores, 0.025), stats::quantile(cure_scores, 0.025), stats::quantile(toxicity_scores * 100, 0.025), stats::quantile(qol_scores * 50, 0.025),
                     stats::quantile(tr_contribution, 0.025), stats::quantile(cure_contribution, 0.025), stats::quantile(toxicity_contribution, 0.025), stats::quantile(qol_contribution, 0.025),
                     stats::quantile(final_utility_vector, 0.025)),
    Upper_95_CrI = c(stats::quantile(tr_scores, 0.975), stats::quantile(cure_scores, 0.975), stats::quantile(toxicity_scores * 100, 0.975), stats::quantile(qol_scores * 50, 0.975),
                     stats::quantile(tr_contribution, 0.975), stats::quantile(cure_contribution, 0.975), stats::quantile(toxicity_contribution, 0.975), stats::quantile(qol_contribution, 0.975),
                     stats::quantile(final_utility_vector, 0.975))
  )
  return(list(
    final_utility_vector = final_utility_vector,
    component_summary = summary_df
  ))
}

#' Plot Utility Score Gauge Visualization
#'
#' Generates a gauge visualization from the output of `summarize_final_utility`.
#' It now interactively asks for a trial name to use as a plot title if one
#' is not provided as an argument.
#'
#' @param final_utility_results The list object returned by `summarize_final_utility`.
#' @param trial_name Optional. A character string for the trial name. If NULL,
#'  the function will prompt the user to enter a name.
#'
#' @return A ggplot2 object representing the gauge plot.
#' @importFrom graphics plot
#' @importFrom grDevices dev.off png
#' @importFrom stats median
#' @export
#'
plot_utility_donut <- function(final_utility_results, trial_name = NULL) {

  # --- 1. Cargar paquetes necesarios ---
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2")
  if (!requireNamespace("viridis", quietly = TRUE)) stop("Please install viridis")

  # --- 2. [NUEVO] Preguntar por el nombre del ensayo si no se proporciona ---
  if (is.null(trial_name)) {
    trial_name <- readline(prompt = "Enter the trial name for the plot title: ")
  }

  # --- 3. Extraer medianas del objeto de resultados ---
  summary_data <- final_utility_results$component_summary

  median_tr_contribution <- summary_data$Median[summary_data$Component == "Contribution TR (Weighted)"]
  median_cure_contribution <- summary_data$Median[summary_data$Component == "Contribution Cure (Weighted)"]
  median_tox_contribution <- summary_data$Median[summary_data$Component == "Contribution Toxicity (Weighted)"]
  median_qol_contribution <- summary_data$Median[summary_data$Component == "Contribution QoL (Weighted)"]
  median_final <- summary_data$Median[summary_data$Component == "FINAL UTILITY SCORE"]

  median_eff_total_contribution <- median_tr_contribution + median_cure_contribution

  if (length(median_eff_total_contribution) == 0) median_eff_total_contribution <- 0
  if (length(median_tox_contribution) == 0) median_tox_contribution <- 0
  if (length(median_qol_contribution) == 0) median_qol_contribution <- 0
  if (length(median_final) == 0) median_final <- 0

  # --- 4. Preparar datos para el gráfico ---
  plot_max_value <- 100
  median_final_clamped <- max(0, min(plot_max_value, median_final))
  final_score_color <- viridis::viridis(101)[round(median_final_clamped) + 1]

  # --- 5. Generar el gráfico ---
  donut_plot <- ggplot2::ggplot() +
    ggplot2::geom_rect(ggplot2::aes(ymax = plot_max_value, ymin = 0, xmax = 4, xmin = 3), fill = "#F0F0F0") +
    ggplot2::geom_rect(ggplot2::aes(ymax = median_final_clamped, ymin = 0, xmax = 4, xmin = 3), fill = final_score_color) +
    ggplot2::annotate("text", x = 0, y = 0,
                      label = paste0(
                        "Efficacy: ", ifelse(median_eff_total_contribution >= 0, "+", ""), round(median_eff_total_contribution, 1),
                        "\nToxicity: ", round(median_tox_contribution, 1),
                        "\nQoL: ", ifelse(median_qol_contribution >= 0, "+", ""), round(median_qol_contribution, 1),
                        "\n-----------------",
                        "\nFinal Score: ", round(median_final, 1)
                      ),
                      size = 5.5, lineheight = 1.1, hjust = 0.5, vjust = 0.5) +
    ggplot2::coord_polar(theta = "y", start = 0) +
    ggplot2::xlim(c(0, 5)) +
    ggplot2::ylim(c(0, plot_max_value)) +

    # --- 6. [NUEVO] Añadir título y ajustar su estilo ---
    ggplot2::labs(title = paste("Clinical Utility Score:", trial_name)) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16, margin = ggplot2::margin(b = 5))
    )

  return(donut_plot)
}



#' Plot an Elegant Histogram of QoL Outcomes
#'
#' @description
#' This function takes a numeric vector of sampled Quality of Life (QoL) scores
#' and creates an elegant bar plot that functions as a histogram. It visualizes
#' the distribution of outcomes, showing the probability of each category on the
#' y-axis and using descriptive labels on the x-axis.
#'
#' @details
#' The function first calculates the proportion (probability) of each unique score
#' `[-2, -1, 0, 1, 2]`. It then maps these numeric scores to meaningful labels
#' (e.g., "Significant Deterioration") and uses a diverging color palette to
#' intuitively represent negative, neutral, and positive outcomes. The bars are
#' directly labeled with their corresponding probabilities.
#'
#' @param qol_scores_vector A numeric vector of sampled QoL scores, where each
#'   element is one of `[-2, -1, 0, 1, 2]`. This is typically the output from
#'   `sample_qol_scores()`.
#'
#' @return A `ggplot` object representing the customized histogram. This object
#'   can be printed to display the plot or be further modified.
#'
#' @export
#' @importFrom dplyr count mutate
#' @importFrom ggplot2 ggplot aes geom_col geom_text labs scale_y_continuous scale_fill_manual theme_minimal theme element_text
#' @importFrom scales percent
#'
#' @examples
#' # --- 1. First, generate a vector of sampled scores ---
#' # prob_vector <- c(0.05, 0.20, 0.50, 0.20, 0.05)
#' # scores <- sample_qol_scores(prob_vector, n_samples = 2000)
#'
#' # --- 2. Now, plot the results with the new function ---
#' # final_plot <- plot_qol_histogram(scores)
#'
#' # --- 3. Display the plot ---
#' # print(final_plot)
#'
plot_qol_histogram <- function(qol_scores_vector) {

  # --- 1. Data Preparation ---
  # This step is key to transforming the raw score data into a format
  # suitable for a well-labeled plot.

  # Define the labels that correspond to the numeric scores.
  score_labels <- c(
    "-2" = "Sig. Deterioration",
    "-1" = "Slight Deterioration",
    "0"  = "No Change",
    "1"  = "Slight Improvement",
    "2"  = "Sig. Improvement"
  )

  # Create a data frame and calculate the probability of each score.
  plot_data <- data.frame(score = qol_scores_vector) |>
    dplyr::count(.data$score, name = "count") |>
    dplyr::mutate(
      probability = .data$count / sum(.data$count),
      # Convert the numeric score into an ordered factor with descriptive labels.
      # This ensures the x-axis is sorted correctly and displays the "true names".
      score_label = factor(
        .data$score,
        levels = names(score_labels),
        labels = score_labels
      )
    )

  # --- 2. Plot Creation ---
  # We build the plot layer by layer using ggplot2 for maximum control.

  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$score_label, y = .data$probability, fill = .data$score_label)) +
    # Use geom_col() because we have already calculated the y-values (probabilities).
    ggplot2::geom_col(show.legend = FALSE, width = 0.7) +

    # Add text labels on top of each bar to show the exact probability.
    ggplot2::geom_text(
      ggplot2::aes(label = scales::percent(.data$probability, accuracy = 0.1)),
      vjust = -0.4, # Vertically adjust to sit just above the bar.
      size = 3.5,
      color = "gray20"
    ) +

    # Define a custom, diverging color palette to be more intuitive.
    # Reds for negative, gray for neutral, blues for positive.
    ggplot2::scale_fill_manual(
      values = c(
        "Sig. Deterioration"   = "#d73027",
        "Slight Deterioration" = "#fc8d59",
        "No Change"            = "gray80",
        "Slight Improvement"   = "#91bfdb",
        "Sig. Improvement"     = "#4575b4"
      )
    ) +

    # Format the y-axis to display probabilities as percentages.
    ggplot2::scale_y_continuous(
      labels = scales::percent,
      limits = c(0, max(plot_data$probability) * 1.15), # Add space for labels
      expand = c(0, 0) # Make the bars start directly at the y=0 line
    ) +

    # Set all plot labels and title.
    ggplot2::labs(
      title = "Distribution of Quality of Life Outcomes",
      subtitle = "Based on the sampled scores from the trial",
      x = NULL, # Remove x-axis title as labels are self-explanatory
      y = "Probability"
    ) +

    # Apply a clean and minimal theme.
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "gray40", size = 11),
      panel.grid.major.x = ggplot2::element_blank(), # Remove vertical grid lines
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 0, vjust = 0.5) # Adjust if labels overlap
    )
}








#' Interactively or non-interactively simulate a realistic clinical trial dataset
#'
#' This function generates a flattened list object containing simulated data for
#' toxicity, patient numbers, and Quality of Life (QoL).
#'
#' @param n_control Number of patients in the control arm. If NULL, prompts interactively.
#' @param ratio_str Character string indicating the randomization ratio (e.g., "2:1"). If NULL, prompts interactively.
#' @param control_g1_4_pct Overall percentage of any-grade toxicity (grades 1–4) in the control arm. If NULL, prompts interactively.
#' @param control_g3_4_pct Overall percentage of severe (grade 3–4) toxicity in the control arm. If NULL, prompts interactively.
#' @param tox_ratio Relative toxicity ratio (experimental vs control arm). If NULL, prompts interactively.
#' @param qol_scenario Integer from 1 to 5 indicating the predefined QoL scenario. If NULL, prompts interactively.
#' @param qol_strength Integer from 1 to 5 indicating the strength of evidence for the QoL effect. If NULL, prompts interactively.
#'
#' @return A list with three main components:
#' \describe{
#'   \item{toxicity}{A data frame with simulated adverse event data, including grades and patient IDs.}
#'   \item{N_patients}{A named numeric vector with patient counts in the control and experimental arms.}
#'   \item{qol}{A numeric vector of length 5 representing the final probability distribution for QoL outcomes.}
#' }
#'
#' @importFrom dplyr sample_n %>% mutate lag
#' @importFrom stats rmultinom runif
#' @export


simulate_trial_data <- function(n_control = NULL, ratio_str = NULL,
                                control_g1_4_pct = NULL, control_g3_4_pct = NULL,
                                tox_ratio = NULL, qol_scenario = NULL, qol_strength = NULL) {

  # --- 0. Adverse Event Database ---
  ae_database <- data.frame(
    EventName = c("Fatigue", "Diarrhoea", "Nausea", "Anorexia", "Stomatitis", "Hand-foot skin reaction", "Hypertension", "Neutropenia", "Thrombocytopenia", "Anaemia", "Vomiting", "Alopecia", "Rash", "Alanine aminotransferase increased", "Aspartate aminotransferase increased", "Hyperbilirubinaemia", "Myalgia", "Arthralgia", "Headache", "Dyspnoea", "Cough", "Abdominal pain", "Constipation", "Insomnia", "Pyrexia", "Oedema peripheral", "Mucositis", "Proteinuria", "Hypokalaemia", "Cardiac arrest", "Heart failure"),
    SystemOrganClass = c("General, metabolic, and other disorders", "Gastrointestinal disorders", "Gastrointestinal disorders", "Gastrointestinal disorders", "Gastrointestinal disorders", "Dermatologic disorders", "General, metabolic, and other disorders", "Blood and lymphatic system disorders", "Blood and lymphatic system disorders", "Blood and lymphatic system disorders", "Gastrointestinal disorders", "Dermatologic disorders", "Dermatologic disorders", "General, metabolic, and other disorders", "General, metabolic, and other disorders", "General, metabolic, and other disorders", "General, metabolic, and other disorders", "General, metabolic, and other disorders", "General, metabolic, and other disorders", "Respiratory, thoracic and mediastinal disorders", "Respiratory, thoracic and mediastinal disorders", "Gastrointestinal disorders", "Gastrointestinal disorders", "General, metabolic, and other disorders", "General, metabolic, and other disorders", "General, metabolic, and other disorders", "Gastrointestinal disorders", "General, metabolic, and other disorders", "General, metabolic, and other disorders", "General, metabolic, and other disorders", "General, metabolic, and other disorders"),
    base_prevalence = c(10, 8, 8, 6, 5, 5, 5, 7, 6, 6, 5, 4, 4, 5, 5, 3, 3, 3, 3, 3, 3, 4, 3, 2, 2, 2, 3, 2, 2, 1, 1)
  )

  # --- 1. Check for Interactive Mode ---
  is_interactive <- is.null(n_control)

  if (is_interactive) {
    cat("--- Running in Interactive Mode ---\n\n")
    cat("--- Section 1: Patient Numbers ---\n")
    n_control <- as.integer(readline(prompt = "Enter number of patients in Control arm: "))
    ratio_str <- readline(prompt = "Enter randomization ratio (e.g., '2:1' for 2 Exp per 1 Ctrl): ")
    cat("\n--- Section 2: Toxicity Parameters ---\n")
    control_g1_4_pct <- as.numeric(readline(prompt = "Enter overall % of any-grade (G1-4) toxicity in Control arm (e.g., 85): "))
    control_g3_4_pct <- as.numeric(readline(prompt = "Enter overall % of severe (G3-4) toxicity in Control arm (e.g., 40): "))
    tox_ratio <- as.numeric(readline(prompt = "Enter toxicity ratio (e.g., 1.2 for 20% more AEs in Exp arm): "))

    # --- [FIX] Expanded and clearer interactive prompts ---
    cat("\n--- Section 3: Quality of Life (QoL) Assessment ---\n")
    cat("Select the main Quality of Life (QoL) scenario:\n")
    cat("1: Significant Improvement\n")
    cat("2: Stabilization / Probable Benefit\n")
    cat("3: No Difference / Marginal Benefit\n")
    cat("4: Deterioration\n")
    cat("5: Insufficient Data / Unknown\n")
    qol_scenario <- readline(prompt = "Enter a number (1-5): ")

    cat("\nSelect the strength/quality of the evidence for this scenario:\n")
    cat("1: Very Low\n")
    cat("2: Low\n")
    cat("3: Moderate\n")
    cat("4: High\n")
    cat("5: Very High\n")
    qol_strength <- as.integer(readline(prompt = "Enter a number (1-5): "))
  }

  # --- 2. Process Inputs and Generate Data ---
  parts <- strsplit(ratio_str, ":")[[1]]
  if (length(parts) != 2 || is.na(as.numeric(parts[1])) || is.na(as.numeric(parts[2]))) {
    stop("Invalid ratio format. Please use the format 'X:Y', e.g., '2:1' or '1:1'.")
  }
  ratio <- as.numeric(parts[1]) / as.numeric(parts[2])
  n_experimental <- round(n_control * ratio)
  N_patients <- c(n_experimental, n_control)
  names(N_patients) <- c("Experimental", "Control")

  base_vectors <- list("1" = c(0.05, 0.1, 0.15, 0.3, 0.4), "2" = c(0.05, 0.1, 0.2, 0.5, 0.15), "3" = c(0.1, 0.15, 0.5, 0.15, 0.1), "4" = c(0.4, 0.3, 0.15, 0.1, 0.05), "5" = c(0.2, 0.2, 0.2, 0.2, 0.2))
  alpha <- (qol_strength - 1) / 4
  qol_vector <- alpha * base_vectors[[as.character(qol_scenario)]] + (1 - alpha) * rep(0.2, 5)
  names(qol_vector) <- c("P(-2)", "P(-1)", "P(0)", "P(+1)", "P(+2)")

  n_events <- sample(25:nrow(ae_database), 1)
  trial_aes <- dplyr::sample_n(ae_database, n_events, weight = ae_database$base_prevalence, replace = FALSE)

  total_events_g1_4_ctrl <- round(n_control * (control_g1_4_pct / 100))
  event_counts_g1_4_ctrl <- stats::rmultinom(1, total_events_g1_4_ctrl, prob = trial_aes$base_prevalence)
  total_events_g3_4_ctrl <- round(n_control * (control_g3_4_pct / 100))
  event_counts_g3_4_ctrl <- stats::rmultinom(1, total_events_g3_4_ctrl, prob = trial_aes$base_prevalence)
  event_counts_g3_4_ctrl <- pmin(event_counts_g1_4_ctrl, event_counts_g3_4_ctrl)
  Incidence_G1_4_Control <- round(event_counts_g1_4_ctrl / n_control * 100)
  Incidence_G3_4_Control <- round(event_counts_g3_4_ctrl / n_control * 100)
  Incidence_G1_4_Experimental <- round(Incidence_G1_4_Control * tox_ratio * stats::runif(n_events, 0.85, 1.15))
  Incidence_G3_4_Experimental <- round(Incidence_G3_4_Control * tox_ratio * stats::runif(n_events, 0.90, 1.25))
  Incidence_G1_4_Experimental[Incidence_G1_4_Experimental > 100] <- sample(90:100, 1)
  Incidence_G3_4_Experimental <- pmin(Incidence_G1_4_Experimental, Incidence_G3_4_Experimental)

  toxicity <- data.frame(
    EventName = trial_aes$EventName,
    SystemOrganClass = trial_aes$SystemOrganClass,
    Incidence_G1_4_Experimental,
    Incidence_G3_4_Experimental,
    Incidence_G1_4_Control,
    Incidence_G3_4_Control
  )

  # --- 3. Assemble the Final Flattened List Object ---
  final_list <- list(
    toxicity = toxicity,
    N_patients = N_patients,
    qol = qol_vector
  )

  if(is_interactive) cat("\n--- Simulation Complete! ---\n")
  return(final_list)
}
