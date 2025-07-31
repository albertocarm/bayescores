

#' @title Run a Scenario for a Sensitivity Plot
#' @description Generates the data needed for a single sensitivity plot line chart.
#' @param x_var_name The parameter to vary on the x-axis.
#' @param x_range The range for the x-axis variable.
#' @param line_var_name The parameter to vary as different lines.
#' @param line_levels The levels for the line variable.
#' @param fixed_params A list of fixed parameters.
#' @param cal_args A list of calibration arguments.
#' @return A data frame with simulation results.
#' @importFrom stats median
#' @export
run_plot_scenario <- function(x_var_name, x_range, line_var_name, line_levels, fixed_params, cal_args) {

  scenarios <- expand.grid(
    x_val = seq(x_range[1], x_range[2], length.out = 50),
    line_val = line_levels,
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(scenarios)) {
    scenario <- scenarios[i, ]
    sim_params <- fixed_params
    sim_params[[x_var_name]] <- scenario$x_val
    sim_params[[line_var_name]] <- scenario$line_val

    n_samples <- 500
    tr_samples <- rep(sim_params$TR, n_samples)
    cure_samples <- rep(sim_params$Cure, n_samples)
    qol_samples <- rep(sim_params$QoL, n_samples)
    toxicity_samples <- rep(sim_params$Toxicity, n_samples)

    score_list <- get_bayescores(
      efficacy_inputs = list(tr_posterior_samples = tr_samples, cure_posterior_samples = cure_samples),
      qol_scores = qol_samples,
      toxicity_scores = toxicity_samples,
      calibration_args = cal_args
    )
    scenarios[i, "Median_Utility"] <- median(score_list$final_utility_vector)
  }

  colnames(scenarios) <- c(x_var_name, line_var_name, "Median_Utility")
  return(scenarios)
}







#' @title Generate Comprehensive Sensitivity Analysis Dashboard
#' @description Orchestrates all simulations and generates the final 8-panel
#'   publication-ready plot.
#' @param calibration_args A list containing the final calibration parameters
#'   for the utility model (e.g., efficacy targets).
#' @return A composite ggplot object created with the patchwork package.
#' @import ggplot2
#' @import patchwork
#' @import ggsci
#' @importFrom dplyr arrange
#' @export
generate_sensitivity_dashboard <- function(calibration_args) {

  # --- 1. Data generation for the 8 plots ---
  cat("Generating data for all 8 plots...\n")
  plot_a_data <- run_plot_scenario("TR", c(1, 4), "Cure", c(0.0, 0.1, 0.2), fixed_params = list(Toxicity = -0.2, QoL = 0), calibration_args)
  plot_b_data <- run_plot_scenario("Cure", c(0, 0.5), "TR", c(1.0, 1.5, 2.0), fixed_params = list(Toxicity = -0.2, QoL = 0), calibration_args)
  plot_c_data <- run_plot_scenario("Toxicity", c(1, -1), "Cure", c(0.05, 0.1, 0.15), fixed_params = list(TR = 1.3, QoL = 0), calibration_args)
  plot_d_data <- run_plot_scenario("QoL", c(-2, 2), "Cure", c(0.05, 0.1, 0.15), fixed_params = list(TR = 1.3, Toxicity = -0.2), calibration_args)
  plot_e_data <- run_plot_scenario("Toxicity", c(1, -1), "TR", c(1.1, 1.5, 2.0), fixed_params = list(Cure = 0.05, QoL = 0), calibration_args)
  plot_f_data <- run_plot_scenario("QoL", c(-2, 2), "TR", c(1.1, 1.5, 2.0), fixed_params = list(Cure = 0.05, Toxicity = -0.2), calibration_args)
  plot_g_data <- run_plot_scenario("QoL", c(-2, 2), "Toxicity", c(0.2, -0.2, -0.6, -1.0), fixed_params = list(Cure = 0.1, TR = 1.5), calibration_args)
  plot_h_data <- run_plot_scenario("Toxicity", c(1, -1), "QoL", c(-1, 0, 1, 2), fixed_params = list(Cure = 0.1, TR = 1.5), calibration_args)
  cat(" Data generation complete.\n")

  # --- 2. PASO CLAVE: Ordenar los datos internamente ---
  # Esto hace que la función sea robusta y funcione en cualquier entorno.
  plot_a_data <- plot_a_data %>% dplyr::arrange(Cure, TR)
  plot_b_data <- plot_b_data %>% dplyr::arrange(TR, Cure)
  plot_c_data <- plot_c_data %>% dplyr::arrange(Cure, Toxicity)
  plot_d_data <- plot_d_data %>% dplyr::arrange(Cure, QoL)
  plot_e_data <- plot_e_data %>% dplyr::arrange(TR, Toxicity)
  plot_f_data <- plot_f_data %>% dplyr::arrange(TR, QoL)
  plot_g_data <- plot_g_data %>% dplyr::arrange(Toxicity, QoL)
  plot_h_data <- plot_h_data %>% dplyr::arrange(QoL, Toxicity)

  # --- 3. Creation and combination of plots ---
  p_a <- ggplot(plot_a_data, aes(x=TR, y=Median_Utility, color=factor(Cure*100))) + geom_line(linewidth=1) + scale_color_jco() + labs(title="A: Impact of Time Ratio", subtitle="Fixed: Toxicity=-0.2, QoL=0", y="Utility Score", x="Time Ratio (TR)", color="Cure Rate (%)")
  p_b <- ggplot(plot_b_data, aes(x=Cure, y=Median_Utility, color=factor(TR))) + geom_line(linewidth=1) + scale_color_jco() + scale_x_continuous(labels=scales::percent) + labs(title="B: Impact of Cure Rate", subtitle="Fixed: Toxicity=-0.2, QoL=0", y=NULL, x="Cure Rate", color="Time Ratio")
  p_c <- ggplot(plot_c_data, aes(x=Toxicity, y=Median_Utility, color=factor(Cure*100))) + geom_line(linewidth=1) + scale_color_jco() + scale_x_reverse() + labs(title="C: Impact of Toxicity", subtitle="Fixed: TR=1.3, QoL=0", y="Utility Score", x="Toxicity Score", color="Cure Rate (%)")
  p_d <- ggplot(plot_d_data, aes(x=QoL, y=Median_Utility, color=factor(Cure*100))) + geom_line(linewidth=1) + scale_color_jco() + labs(title="D: Impact of Quality of Life", subtitle="Fixed: TR=1.3, Toxicity=-0.2", y=NULL, x="QoL Score", color="Cure Rate (%)")
  p_e <- ggplot(plot_e_data, aes(x=Toxicity, y=Median_Utility, color=factor(TR))) + geom_line(linewidth=1) + scale_color_jco() + scale_x_reverse() + labs(title="E: TR as a Modulator of Toxicity", subtitle="Fixed: Cure=5%, QoL=0", y="Utility Score", x="Toxicity Score", color="Time Ratio")
  p_f <- ggplot(plot_f_data, aes(x=QoL, y=Median_Utility, color=factor(TR))) + geom_line(linewidth=1) + scale_color_jco() + labs(title="F: TR as a Modulator of QoL", subtitle="Fixed: Cure=5%, Toxicity=-0.2", y=NULL, x="QoL Score", color="Time Ratio")
  p_g <- ggplot(plot_g_data, aes(x=QoL, y=Median_Utility, color=factor(Toxicity))) + geom_line(linewidth=1) + scale_color_jco() + labs(title="G: QoL vs. Toxicity Interaction", subtitle="Fixed: TR=1.5, Cure=10%", x="QoL Score", y="Utility Score", color="Toxicity Score")
  p_h <- ggplot(plot_h_data, aes(x=Toxicity, y=Median_Utility, color=factor(QoL))) + geom_line(linewidth=1) + scale_color_jco() + scale_x_reverse() + labs(title="H: Toxicity vs. QoL Interaction", subtitle="Fixed: TR=1.5, Cure=10%", x="Toxicity Score", y=NULL, color="QoL Score")

  dashboard <- (p_a + p_b) / (p_c + p_d) / (p_e + p_f) / (p_g + p_h)

  final_dashboard <- dashboard +
    plot_annotation(
      title = "Bayescores Model: Comprehensive Sensitivity Analysis",
      caption = "All scenarios are simulated using the final, inherently bounded, efficacy-dependent model.",
      theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
    ) &
    theme_bw(base_size = 12) &
    coord_cartesian(ylim = c(0, 100)) &
    theme(plot.title = element_text(size=12), plot.subtitle = element_text(size=9, face="italic"))

  return(final_dashboard)
}





