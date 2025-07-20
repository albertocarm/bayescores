

#' @title Fit the Bayesian Cure Model (Engine)
#' @description Prepares data and calls the external Stan cure model.
#'
#' @param data A data frame with time, event, and arm columns.
#' @param time_col,event_col,arm_col Character strings for column names.
#' @param chains,iter,warmup,seed Numeric arguments passed to `rstan::stan`.
#' @param adapt_delta Target acceptance rate for Stan's NUTS algorithm. A higher
#'   value (e.g., 0.99) can help with divergent transitions.
#' @param ... Additional arguments passed to `rstan::stan`.
#'
#' @return A custom S3 object of class `bcm_fit`.
#' @importFrom rstan stan
#' @export
fit_bayesian_cure_model <- function(data,
                                    time_col = "time",
                                    event_col = "event",
                                    arm_col = "arm",
                                    chains = 4,
                                    iter = 2000,
                                    warmup = 1000,
                                    seed = 555,
                                    adapt_delta = 0.99,
                                    ...) {

  # Validate and convert 'arm' to a factor
  if (!is.factor(data[[arm_col]])) {
    message(paste0("Note: Converting column '", arm_col, "' to a factor."))
    data[[arm_col]] <- as.factor(data[[arm_col]])
  }

  stan_model_path <- system.file("stan", "cure_model.stan", package = "bayesCure")
  if (stan_model_path == "") {
    stop("Could not find Stan model file 'cure_model.stan'.")
  }

  stan_data <- list(
    N = nrow(data),
    tiempo = data[[time_col]],
    evento = data[[event_col]],
    arm = as.numeric(data[[arm_col]]) - 1
  )

  # Control list for Stan's algorithm
  control_list <- list(adapt_delta = adapt_delta)

  stan_fit <- rstan::stan(
    file = stan_model_path,
    data = stan_data,
    chains = chains, iter = iter, warmup = warmup, seed = seed,
    control = control_list,
    ...
  )

  result <- list(
    stan_fit = stan_fit,
    original_data = data,
    column_map = list(time_col = time_col, event_col = event_col, arm_col = arm_col)
  )
  class(result) <- "bcm_fit"

  return(result)
}

#' Plot Model Diagnostics
#'
#' @description A generic function to plot model diagnostics.
#' @param x An object, typically a model fit.
#' @param ... Other arguments.
#' @export
model_diagnostics <- function(x, ...) {
  UseMethod("model_diagnostics")
}
