
#' @title Fit the Bayesian Cure Model (Engine)
#' @description Prepares data, calls the external Stan cure model, and stores MCMC draws.
#'
#' @param data A data frame with time, event, and arm columns.
#' @param time_col,event_col,arm_col Character strings for column names.
#' @param tail_assumption Character string. Sets the prior belief strategy for the adjuvant
#'    cure effect. One of "neutral", "immature_skeptical", or "biologically_null".
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
#' @param dist Character string. Parametric distribution assumed for the event
#'    times of the non-cured patients. One of "weibull" (default), "exponential",
#'    "lognormal", "loglogistic", "gamma", or "gengamma" (generalized gamma).
#'    All of them are accelerated failure time models, so the time ratio keeps
#'    its interpretation. The Weibull distribution is the pre-specified default;
#'    the alternatives are intended for sensitivity analyses and can optionally
#'    be compared with `compare_models`.
#' @param ... Additional arguments passed to `rstan::stan`.
#'
#' @details
#' For every distribution the location of the non-cured patients on the log-time
#' scale is `beta_surv_intercept + beta_surv_arm * arm`, and `exp(beta_surv_arm)`
#' is the time ratio. The shape parameter `alpha_control` (multiplied by
#' `exp(beta_alpha_arm)` in the experimental arm when `shared_shape = FALSE`)
#' has the following meaning:
#' \itemize{
#'   \item "weibull": Weibull shape; `exp(location)` is the Weibull scale.
#'   \item "exponential": no shape parameter (fixed at 1); `exp(location)` is the mean.
#'   \item "lognormal": inverse of the standard deviation of the log times;
#'     `exp(location)` is the median.
#'   \item "loglogistic": log-logistic shape; `exp(location)` is the median.
#'   \item "gamma": gamma shape; `exp(location)` is the mean.
#'   \item "gengamma": inverse of the scale parameter `sigma` of the generalized
#'     gamma distribution in the parameterization of Prentice (1974), with
#'     `location` as `mu` and the additional shape parameter `q_shape`, common
#'     to both arms, as `Q` (`Q = 1` is the Weibull and `Q = 0` the log-normal
#'     distribution).
#' }
#' The priors of the parameters shared with the Weibull model are the same for
#' every distribution; `q_shape` has a standard normal prior. With
#' "exponential" the argument `shared_shape` has no effect.
#'
#' Priors of the treatment effects. With `tail_assumption = "neutral"` (the
#' default) the log odds ratio of cure (`beta_cure_arm`) and the log time ratio
#' (`beta_surv_arm`) have independent Normal(0, 2.5) priors. A standard
#' deviation of 2.5 is the scale of the weakly informative default of Gelman
#' et al. (2008): a shift of 5 on the logit scale moves a probability from
#' 0.01 to 0.5, so larger effects are implausible but not excluded. The prior
#' places 95\% of its mass on odds ratios and time ratios between 0.007 and
#' 134 and is nearly flat over the range of effects seen in clinical trials,
#' so it stabilizes the computation without pulling the estimates towards the
#' null. With "immature_skeptical" and "biologically_null" the log odds ratio
#' has a Laplace(0, 0.2) and a Laplace(0, 0.1) prior, respectively, and the log
#' time ratio keeps the Normal(0, 2.5) prior. The intercepts have
#' Student-t(4, 0, 2.5) priors and the shape parameter a Gamma(1, 1) prior.
#'
#' @references
#' Gelman A, Jakulin A, Pittau MG, Su YS (2008). A weakly informative default
#' prior distribution for logistic and other regression models. The Annals of
#' Applied Statistics, 2, 1360-1383.
#'
#' Prentice RL (1974). A log gamma model and its maximum likelihood estimation.
#' Biometrika, 61, 539-544.
#'
#' @return A custom S3 object of class `bcm_fit`, a list with elements:
#'    - `stan_fit`: the original `stanfit` object
#'    - `original_data`: the input data frame
#'    - `column_map`: list mapping `time_col`, `event_col`, `arm_col`
#'    - `posterior_draws`: list of posterior samples for each parameter
#'    - `n_draws`: the total number of post-warmup posterior draws
#'    - `dist`: the distribution used for the non-cured patients
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
                                    dist          = "weibull",
                                    ...) {

  tail_assumption <- match.arg(tail_assumption, c("neutral", "immature_skeptical", "biologically_null"))
  dist <- match.arg(dist, c("weibull", "exponential", "lognormal", "loglogistic", "gamma", "gengamma"))
  if (dist == "exponential") {
    shared_shape <- TRUE
  }

  # Validate and convert 'arm' to a factor
  if (!is.factor(data[[arm_col]])) {
    message(paste0("Note: converting '", arm_col, "' to factor."))
    data[[arm_col]] <- as.factor(data[[arm_col]])
  }

  # Locate Stan model file
  stan_model_file <- if (dist == "weibull") "cure_model.stan" else "cure_model_dist.stan"
  stan_model_path <- system.file("stan", stan_model_file, package = "bayescores")
  if (stan_model_path == "") {
    stop(paste0("Cannot find '", stan_model_file, "'. Is the package installed correctly?"))
  }

  prior_int_map <- c("neutral" = 1, "immature_skeptical" = 2, "biologically_null" = 3)
  cure_prior_type <- prior_int_map[tail_assumption]

  # Prepare data for Stan
  stan_data <- list(
    N      = nrow(data),
    tiempo = data[[time_col]],
    evento = data[[event_col]],
    arm    = as.numeric(data[[arm_col]]) - 1,

    cure_prior_type = as.integer(cure_prior_type),

    use_historical_prior = as.integer(use_historical_prior),
    historical_prior_params = historical_prior_params,

    # Argument for shape flexibility
    shared_shape = as.integer(shared_shape)
  )
  if (dist != "weibull") {
    stan_data$dist <- match(dist, c("exponential", "lognormal", "loglogistic", "gamma", "gengamma"))
  }

  # Stan control settings
  control_list <- list(adapt_delta = adapt_delta)

  # Initialization function
  init_fun <- function() {
    list(beta_cure_arm_raw = rnorm(1, 0, 0.1))
  }

  # Fit model
  if (dist == "weibull") {
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
  } else {
    stan_fit <- rstan::sampling(
      .cure_dist_stanmodel(stan_model_path),
      data    = stan_data,
      init    = init_fun,
      chains  = chains,
      iter    = iter,
      warmup  = warmup,
      seed    = seed,
      control = control_list,
      ...
    )
  }

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
    shared_shape    = shared_shape,
    dist            = dist
  )

  class(result) <- "bcm_fit"

  return(result)
}

# Compiled model cache
.cure_model_cache <- new.env(parent = emptyenv())

.cure_dist_stanmodel <- function(stan_model_path) {
  if (is.null(.cure_model_cache$dist_model)) {
    .cure_model_cache$dist_model <- rstan::stan_model(file = stan_model_path, model_name = "anon_model")
  }
  .cure_model_cache$dist_model
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
