

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



#' @title Calculate a complete toxicity analysis
#' @description Calculates weighted toxicity scores (WTS).
#' @param trial_data A list containing the toxicity data.
#' @param n_simulations An integer for the number of simulations.
#' @param soc_weights The default weights for toxicity categories
#' @param unacceptable_rel_increase A number for the unacceptable increase limit.
#' @param k_uncertainty A number for the uncertainty constant.
#' @return A list containing the results.
#' @export
#' @importFrom stats rnorm

calculate_toxicity_analysis <- function(
    trial_data,
    n_simulations,
    unacceptable_rel_increase = 0.5,
    k_uncertainty = 5,
    soc_weights = c(
      "Gastrointestinal disorders" = 1.2,
      "Blood and lymphatic system disorders" = 1.6,
      "General, metabolic, and other disorders" = 1.3,
      "Dermatologic disorders" = 1.1,
      "Infections and infestations" = 1.6,
      "Respiratory, thoracic and mediastinal disorders" = 2.5
    )
) {

  # --- 1. Calculate Weighted Toxicity Scores (WTS) ---
  # The soc_weights are now taken directly from the function's arguments.
  W_prom_1_2 <- 1.5; W_prom_3_4 <- 6.0
  toxicity_data <- trial_data$toxicity
  inc_g1_2_exp <- pmax(0, toxicity_data$Incidence_G1_4_Experimental - toxicity_data$Incidence_G3_4_Experimental)
  inc_g1_2_con <- pmax(0, toxicity_data$Incidence_G1_4_Control - toxicity_data$Incidence_G3_4_Control)
  score_base_exp <- (inc_g1_2_exp / 100 * W_prom_1_2) + (toxicity_data$Incidence_G3_4_Experimental / 100 * W_prom_3_4)
  score_base_con <- (inc_g1_2_con / 100 * W_prom_1_2) + (toxicity_data$Incidence_G3_4_Control / 100 * W_prom_3_4)

  # Match the weights from the argument to the data's SystemOrganClass
  w_soc_applied <- soc_weights[toxicity_data$SystemOrganClass]
  # Any SystemOrganClass not in the soc_weights list gets a default weight of 1.0
  w_soc_applied[is.na(w_soc_applied)] <- 1.0

  wts_scores <- data.frame(
    Experimental = sum(score_base_exp * w_soc_applied),
    Control = sum(score_base_con * w_soc_applied)
  )

  # --- 2. Calculate Parameters for the Normal Distribution ---
  wts_diff_estimate <- wts_scores$Experimental - wts_scores$Control
  dynamic_unacceptable_diff <- pmax(wts_scores$Control, 1) * unacceptable_rel_increase
  mu <- (wts_diff_estimate / dynamic_unacceptable_diff) * -1
  total_n <- sum(trial_data$N_patients)
  sigma <- k_uncertainty / sqrt(total_n)

  # --- 3. Generate Posterior Distribution ---
  toxicity_effect_dist <- rnorm(n_simulations, mean = mu, sd = sigma)

  # --- 4. Return Final Results Object ---
  results <- list(
    wts_scores = wts_scores,
    toxicity_effect_vector = pmax(-1, pmin(1, toxicity_effect_dist))
  )
  return(results)
}



#' @title Plot Toxicity Score Density
#' @description Creates a modern density plot of the toxicity effect score.
#' @param analysis_output A list produced by `calculate_toxicity_analysis`.
#' @return A ggplot object.
#'
plot_toxicity_density <- function(analysis_output) {

  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("ggridges", quietly = TRUE)) {
    stop("Packages 'ggplot2' and 'ggridges' are needed.", call. = FALSE)
  }

  # Extraer el vector necesario de la lista de entrada.
  # Esto resuelve el punto que mencionaste: sabe buscar dentro del objeto.
  toxicity_effect_vector <- analysis_output$toxicity_effect_vector

  plot_data <- data.frame(effect = toxicity_effect_vector, category = "Probability Density")
  mean_effect <- mean(plot_data$effect)

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = effect, y = category)) +
    # AJUSTE: 'scale' reducido a un valor más estándar para evitar picos extremos.
    ggridges::geom_density_ridges_gradient(ggplot2::aes(fill = after_stat(x)), rel_min_height = 0.01, scale = 2) +
    ggplot2::scale_fill_gradient2(low = "#440154FF", mid = "#440154FF", high = "#FDE725FF", midpoint = 0) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
    ggplot2::geom_vline(xintercept = mean_effect, linetype = "dotted", color = "black", size = 1) +
    ggplot2::annotate("text", x = mean_effect, y = Inf, label = paste("Mean =", round(mean_effect, 2)), vjust = 1.5, hjust = if (mean_effect > 0) 1.1 else -0.1, fontface = "bold", size = 4) +
    ggplot2::labs(title = "Density of the Toxicity Effect Score", subtitle = "Distribution of the difference in WTS (Experimental vs. Control)", x = "Toxicity Effect Score", y = "") +
    ggridges::theme_ridges(font_size = 14, grid = TRUE) +
    ggplot2::theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = margin(6, 12, 6, 6))

  return(p)
}
