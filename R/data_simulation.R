# =============================================================================
#
# bayesCure: Package Functions v2.0
#
# This script contains all functions for the bayesCure package.
# Alberto Carmona-Bayonas - Murcia, 2025
#
# =============================================================================


#' Simulate Survival Data from a Cure Model
#'
#' @description
#' Generates a dataset for survival analysis based on a two-arm
#' (Control vs. Experimental) clinical trial scenario. It uses a mixture
#' cure model, where a fraction of subjects are considered "cured". For the
#' "susceptible" subjects, their time-to-event follows a Weibull distribution.
#'
#' @param n_patients The total number of patients to simulate.
#' @param cure_fraction_ctrl The proportion (0-1) of "cured" patients in the Control arm.
#' @param cure_fraction_exp The proportion (0-1) of "cured" patients in the Experimental arm.
#' @param max_follow_up The maximum patient follow-up time (administrative censoring).
#' @param weibull_shape The shape parameter of the Weibull distribution.
#' @param median_survival_ctrl The median survival time for susceptible patients in the Control arm.
#' @param time_ratio_exp The ratio of the experimental arm's median survival time to the control arm's.
#' @param seed An optional integer to set the random seed for reproducibility.
#'
#' @return A data.frame with the class `survival_sim_data`.
#' @importFrom stats rbinom rweibull
#' @export
#' @examples
#' \dontrun{
#'   sim_data <- simulate_weibull_cure_data(n_patients = 100,
#'                                          cure_fraction_ctrl = 0.3,
#'                                          cure_fraction_exp = 0.6,
#'                                          max_follow_up = 60,
#'                                          weibull_shape = 1.5,
#'                                          median_survival_ctrl = 15,
#'                                          time_ratio_exp = 2.0)
#'   head(sim_data)
#' }
simulate_weibull_cure_data <- function(n_patients,
                                       cure_fraction_ctrl,
                                       cure_fraction_exp,
                                       max_follow_up,
                                       weibull_shape,
                                       median_survival_ctrl,
                                       time_ratio_exp,
                                       seed = 42) {

  set.seed(seed)
  median_survival_exp <- median_survival_ctrl * time_ratio_exp
  scale_susc_ctrl <- median_survival_ctrl / (log(2))^(1 / weibull_shape)
  scale_susc_exp <- median_survival_exp / (log(2))^(1 / weibull_shape)
  sim_data <- data.frame(id = 1:n_patients)
  sim_data$arm <- sample(c("Control", "Experimental"),
                         size = n_patients,
                         replace = TRUE,
                         prob = c(0.5, 0.5))
  cure_probability <- ifelse(sim_data$arm == "Control", cure_fraction_ctrl, cure_fraction_exp)
  sim_data$is_cured <- rbinom(n = n_patients, size = 1, prob = cure_probability)
  latent_time <- rep(Inf, n_patients)
  susceptible_ctrl_indices <- which(sim_data$arm == "Control" & sim_data$is_cured == 0)
  susceptible_exp_indices  <- which(sim_data$arm == "Experimental" & sim_data$is_cured == 0)
  latent_time[susceptible_ctrl_indices] <- rweibull(n = length(susceptible_ctrl_indices),
                                                    shape = weibull_shape,
                                                    scale = scale_susc_ctrl)
  latent_time[susceptible_exp_indices]  <- rweibull(n = length(susceptible_exp_indices),
                                                    shape = weibull_shape,
                                                    scale = scale_susc_exp)
  sim_data$time <- pmin(latent_time, max_follow_up)
  sim_data$event <- as.numeric(latent_time <= max_follow_up)
  sim_data$arm <- as.factor(sim_data$arm)
  class(sim_data) <- c("survival_sim_data", "data.frame")
  return(sim_data)
}
