

################################################################################
#
#        SCRIPT FOR CLINICAL TOXICiTY ANALYSIS
#
################################################################################
## --------------------------------------------------------------------------##
## PART 1: THE ANALYSIS ENGINE (CORE FUNCTIONS) ----
## --------------------------------------------------------------------------##
# These functions perform the calculations and do not need to be edited.

################################################################################
#
#       DEFINITIVE SCRIPT WITH CLINICALLY-ADJUSTED WEIGHTS
#
# This script uses updated SOC weights to better reflect the severity of
# hematologic toxicities and provides a more clinically-aligned result.
#
################################################################################

## --------------------------------------------------------------------------##
## PART 1: THE ANALYSIS ENGINE (WITH ADJUSTED WEIGHTS) ----
## --------------------------------------------------------------------------##

#' Calculate Toxicity Scores and Adjustment Distribution
#'
#' Merges the calculation of summary toxicity scores and the generation of the
#' posterior adjustment distribution into a single, streamlined function.
#'
#' @param trial_data A list object from simulation, containing 'toxicity' and 'N_patients'.
#' @param n_simulations The number of simulations for the posterior distribution vector.
#'
#' @return A list containing two elements:
#' \item{summary_scores}{A data frame with the total weighted toxicity scores.}
#' \item{adjustment_vector}{A numeric vector of the posterior distribution (-1 to +1).}
#' @importFrom stats rbeta
#' @export
#'
calculate_toxicity_adjustment <- function(trial_data, n_simulations) {

  # --- Input Validation ---
  if (!is.list(trial_data) || !all(c("toxicity", "N_patients") %in% names(trial_data))) {
    stop("Input 'trial_data' must be a list containing 'toxicity' and 'N_patients'.")
  }

  # --- Part 1: Calculate Scores (Logic from the first function) ---
  soc_weights <- c(
    "Gastrointestinal disorders" = 1.2,
    "Blood and lymphatic system disorders" = 1.6,
    "General, metabolic, and other disorders" = 1.3,
    "Dermatologic disorders" = 1.1,
    "Infections and infestations" = 1.6,
    "Respiratory, thoracic and mediastinal disorders" = 2.5
  )
  W_prom_1_2 <- 1.5
  W_prom_3_4 <- 6.0

  datos_toxicidad <- trial_data$toxicity
  inc_g1_2_exp <- pmax(0, datos_toxicidad$Incidence_G1_4_Experimental - datos_toxicidad$Incidence_G3_4_Experimental)
  inc_g1_2_con <- pmax(0, datos_toxicidad$Incidence_G1_4_Control - datos_toxicidad$Incidence_G3_4_Control)
  score_base_exp <- (inc_g1_2_exp / 100 * W_prom_1_2) + (datos_toxicidad$Incidence_G3_4_Experimental / 100 * W_prom_3_4)
  score_base_con <- (inc_g1_2_con / 100 * W_prom_1_2) + (datos_toxicidad$Incidence_G3_4_Control / 100 * W_prom_3_4)
  w_soc_applied <- soc_weights[datos_toxicidad$SystemOrganClass]
  w_soc_applied[is.na(w_soc_applied)] <- 1.0
  score_os_exp <- score_base_exp * w_soc_applied
  score_os_con <- score_base_con * w_soc_applied

  summary_scores <- data.frame(
    Score = "WTS-A-OS (Adjusted)",
    Experimental = sum(score_os_exp),
    Control = sum(score_os_con)
  )

  # --- Part 2: Generate Distribution Vector (Logic from the second function) ---
  score_exp <- summary_scores$Experimental
  score_ctrl <- summary_scores$Control
  if (score_ctrl == 0) { score_ctrl <- 1e-9 }

  n_experimental <- trial_data$N_patients['Experimental']
  n_control <- trial_data$N_patients['Control']

  ratio_estimate <- score_exp / score_ctrl
  prob_exp_worse_estimate <- (ratio_estimate ^ 2) / ((ratio_estimate ^ 2) + 1)
  kappa <- n_experimental + n_control
  alpha <- prob_exp_worse_estimate * kappa
  beta <- (1 - prob_exp_worse_estimate) * kappa

  prob_exp_worse_dist <- rbeta(n_simulations, shape1 = alpha, shape2 = beta)
  toxicity_utility_dist <- 1 - prob_exp_worse_dist
  adjustment_vector <- 2 * (toxicity_utility_dist - 0.5)

  # --- Return a list with both results ---
  return(list(
    summary_scores = summary_scores,
    adjustment_vector = adjustment_vector
  ))
}


#' Plot the Toxicity Adjustment Distribution
#'
#' Creates a density ridge plot for the toxicity adjustment factor.
#'
#' @param toxicity_results A list object returned by `calculate_toxicity_adjustment`.
#' @param trial_name A character string for the trial name.
#'
#' @return A ggplot2 object.
#' @importFrom ggplot2 ggplot aes labs theme_minimal theme element_blank element_text annotate after_stat unit geom_vline
#' @importFrom ggplot2 scale_fill_viridis_c
#' @export
#'
plot_toxicity_adjustment <- function(toxicity_results, trial_name = "Trial") {

  # --- Extract the vector from the results object ---
  adjustment_vector <- toxicity_results$adjustment_vector

  # Validate that the vector exists
  if (is.null(adjustment_vector)) {
    stop("The 'toxicity_results' object must contain an element named 'adjustment_vector'.")
  }

  # The rest of the function remains unchanged
  plot_data <- data.frame(
    AdjustmentFactor = adjustment_vector,
    Analysis = "Toxicity"
  )
  mean_val <- mean(adjustment_vector)

  p <- ggplot(plot_data, aes(x = AdjustmentFactor, y = Analysis)) +
    ggridges::geom_density_ridges_gradient(
      aes(fill = after_stat(x)),
      scale = 3, rel_min_height = 0.01, show.legend = FALSE
    ) +
    scale_fill_viridis_c(option = "plasma") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "white", linewidth = 0.7, alpha = 0.8) +
    geom_vline(xintercept = mean_val, linetype = "solid", color = "black", linewidth = 1) +
    annotate(
      "text", x = mean_val, y = 0.15, label = paste("Mean =", format(mean_val, digits = 2)),
      angle = 90, hjust = 0, vjust = 1.5, fontface = "bold", color = "black"
    ) +
    labs(
      title = paste("Toxicity Adjustment Distribution:", trial_name),
      x = "Toxicity Adjustment Factor",
      y = NULL,
      caption = "Factor > 0 indicates experimental arm is less toxic | Factor < 0 indicates more toxic"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.caption = element_text(hjust = 0, size = 10, face = "italic")
    )

  return(p)
}



#' Create an AMIT plot (Adverse Event Dot Plot) from a clinical trial object.
#'
#' This version is flexible and allows specifying the names of the elements
#' within the trial object.
#'
#' @param trial_object The list object containing the trial data.
#' @param grade_type A string: "any_grade" (G1–4) or "severe_grade" (G3–4).
#' @param main_title The main title of the plot.
#' @param data_element The name (as a character string) of the element containing the toxicity data frame.
#' @param n_element The name (as a character string) of the element containing the vector with patient counts.
#' @param ... Additional arguments passed to the HH::AEdotplot function.
#'
#' @importFrom stats na.omit
#' @export
#'
create_amit_plot <- function(trial_object,
                             grade_type = "any_grade",
                             main_title = "",
                             data_element = "datos",   # Default value for backward compatibility
                             n_element = "N_pacientes",# Default value for backward compatibility
                             ...) {

  # 1. Check if the HH package is installed
  if (!requireNamespace("HH", quietly = TRUE)) {
    stop("The 'HH' package is required. Please install it with: install.packages('HH')")
  }

  # 2. Extract the data from the object using the specified element names
  if (!data_element %in% names(trial_object)) stop(paste("The object does not contain an element named '", data_element, "'", sep = ""))
  if (!n_element %in% names(trial_object)) stop(paste("The object does not contain an element named '", n_element, "'", sep = ""))

  datos <- trial_object[[data_element]]
  N_pacientes <- trial_object[[n_element]]

  # 3. Select incidence columns according to grade type
  if (grade_type == "any_grade") {
    inc_exp_col <- "Incidence_G1_4_Experimental"
    inc_ctrl_col <- "Incidence_G1_4_Control"
  } else if (grade_type == "severe_grade") {
    inc_exp_col <- "Incidence_G3_4_Experimental"
    inc_ctrl_col <- "Incidence_G3_4_Control"
  } else {
    stop("The 'grade_type' argument must be 'any_grade' or 'severe_grade'.")
  }

  # 4. Calculate the number of patients with the event (nAE)
  # unname() is used to avoid warnings about row names
  datos$nAE_Exp <- (datos[[inc_exp_col]] / 100) * unname(N_pacientes["Experimental"])
  datos$nAE_Ctrl <- (datos[[inc_ctrl_col]] / 100) * unname(N_pacientes["Control"])

  # 5. Create long-format data frames for each arm
  df_exp <- data.frame(
    TRT = "Experimental",
    AE = datos$EventName,
    nAE = datos$nAE_Exp,
    nTRT = unname(N_pacientes["Experimental"]),
    OrgSys = datos$SystemOrganClass
  )

  df_ctrl <- data.frame(
    TRT = "Control",
    AE = datos$EventName,
    nAE = datos$nAE_Ctrl,
    nTRT = unname(N_pacientes["Control"]),
    OrgSys = datos$SystemOrganClass
  )

  # 6. Combine and prepare for plotting
  long_data <- rbind(df_exp, df_ctrl)
  long_data <- stats::na.omit(long_data) # Using stats:: explicitly for clarity

  if (nrow(long_data) == 0) {
    warning("No data to plot after omitting NAs. Check whether the incidence columns are empty.")
    return(invisible(NULL))
  }

  long_data$TRT <- factor(long_data$TRT, levels = c("Control", "Experimental"))
  long_data$AE <- factor(long_data$AE)
  long_data$OrgSys <- factor(long_data$OrgSys)

  # 7. Generate the plot
  HH::AEdotplot(
    AE ~ nAE / nTRT | OrgSys,
    groups = TRT,
    data = long_data,
    main = main_title,
    ...
  )
}

# --- END OF FUNCTION ---
