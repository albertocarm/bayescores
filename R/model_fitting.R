
#' @title Fit the Bayesian Cure Model (Engine)
#' @description Prepares data, calls the external Stan cure model, and stores MCMC draws.
#'
#' @param data A data frame with time, event, and arm columns.
#' @param time_col,event_col,arm_col Character strings for column names.
#' @param tail_assumption Character string. Sets the prior belief strategy for the adjuvant
#'    cure effect. One of "neutral", "immature_skeptical", "biologically_null",
#'    "supportive", or "optimistic".
#' @param use_historical_prior Logical. If TRUE, overrides `tail_assumption` and uses
#'    an informative historical prior defined by `historical_prior_params`.
#'    Default is FALSE.
#' @param historical_prior_params Numeric vector of length 2 (Mean, SD).
#'    Used only if `use_historical_prior = TRUE`. Defines the Normal prior
#'    for the treatment effect log-OR. Default is c(0, 1).
#' @param shared_shape Logical. If TRUE, both treatment arms share the same Weibull
#'    shape parameter (proportional hazards). If FALSE (default), allows independent
#'    shapes for each arm.
#' @param chains,iter,warmup,seed Numeric arguments passed to `rstan::stan`.
#' @param adapt_delta Target acceptance rate for Stan's NUTS algorithm.
#' @param ... Additional arguments passed to `rstan::stan`.
#'
#' @return A custom S3 object of class `bcm_fit`, a list with elements:
#'    - `stan_fit`: the original `stanfit` object
#'    - `original_data`: the input data frame
#'    - `column_map`: list mapping `time_col`, `event_col`, `arm_col`
#'    - `posterior_draws`: list of posterior samples for each parameter
#'    - `n_draws`: the total number of post-warmup posterior draws
#'
#' @importFrom rstan stan extract
#' @export
fit_bayesian_cure_model <- function(data,
                                    time_col      = "time",
                                    event_col     = "event",
                                    arm_col       = "arm",
                                    tail_assumption = "neutral",
                                    shared_shape  = FALSE,
                                    chains        = 4,
                                    iter          = 2000,
                                    warmup        = 1000,
                                    seed          = 555,
                                    adapt_delta   = 0.99,
                                    use_historical_prior = FALSE,
                                    historical_prior_params = c(0, 1),
                                    ...) {

  tail_assumption <- match.arg(tail_assumption, c("neutral", "immature_skeptical", "biologically_null", "supportive", "optimistic"))

  # Validate and convert 'arm' to a factor
  if (!is.factor(data[[arm_col]])) {
    message(paste0("Note: converting '", arm_col, "' to factor."))
    data[[arm_col]] <- as.factor(data[[arm_col]])
  }

  # Locate Stan model file
  stan_model_path <- system.file("stan", "cure_model.stan", package = "bayescores")
  if (stan_model_path == "") {
    stop("Cannot find 'cure_model.stan'. Is the package installed correctly?")
  }

  # Updated mapping to match new Stan logic
  prior_int_map <- c("neutral" = 1, "immature_skeptical" = 2, "biologically_null" = 3, "supportive" = 4, "optimistic" = 5)
  cure_prior_type <- prior_int_map[tail_assumption]

  # Prepare data for Stan
  stan_data <- list(
    N      = nrow(data),
    tiempo = data[[time_col]],
    evento = data[[event_col]],
    arm    = as.numeric(data[[arm_col]]) - 1,

    cure_prior_type = as.integer(cure_prior_type),

    # New arguments mapped to Stan
    use_historical_prior = as.integer(use_historical_prior),
    historical_prior_params = historical_prior_params,

    # Argument for shape flexibility
    shared_shape = as.integer(shared_shape)
  )

  # Stan control settings
  control_list <- list(adapt_delta = adapt_delta)

  # Initialization function
  init_fun <- function() {
    if (cure_prior_type >= 4) {
      # Safe positive start for radial priors (4 and 5)
      list(beta_cure_arm_raw = runif(1, 0.05, 0.2))
    } else {
      # Standard random start for others
      list(beta_cure_arm_raw = rnorm(1, 0, 0.1))
    }
  }

  # Fit model
  stan_fit <- rstan::stan(
    file    = stan_model_path,
    data    = stan_data,
    init    = init_fun,
    chains  = chains,
    iter    = iter,
    warmup  = warmup,
    seed    = seed,
    control = control_list,
    ...
  )

  # Extract posterior draws (excluding warmup)
  posterior_draws <- rstan::extract(stan_fit, permuted = TRUE, inc_warmup = FALSE)

  # Calculate the total number of post-warmup draws
  n_draws <- chains * (iter - warmup)

  # Assemble result
  result <- list(
    stan_fit        = stan_fit,
    original_data   = data,
    column_map      = list(
      time_col  = time_col,
      event_col = event_col,
      arm_col   = arm_col
    ),
    posterior_draws = posterior_draws,
    n_draws         = n_draws,
    shared_shape    = shared_shape
  )

  class(result) <- "bcm_fit"

  return(result)
}









#' Plot Model Diagnostics
#'
#' @description A generic function to plot model diagnostics.
#' @param x An object, typically a model fit.
#' @param ... Other arguments.
#'• @export
model_diagnostics <- function(x, ...) {
  UseMethod("model_diagnostics")
}
