# Log density of the standardized generalized gamma
.gengamma_logdens <- function(w, q) {
  q <- as.numeric(q)
  w <- rep_len(as.numeric(w), length(q))
  q2 <- q^2
  u <- q * w
  u_safe <- ifelse(abs(u) < 0.01, 1, u)
  k <- ifelse(
    abs(u) < 0.01,
    0.5 + u * (1 / 6 + u * (1 / 24 + u * (1 / 120 + u * (1 / 720 + u / 5040)))),
    (expm1(u_safe) - u_safe) / u_safe^2
  )
  q2_safe <- pmax(q2, 0.05)
  const <- ifelse(
    q2 < 0.05,
    -0.5 * log(2 * pi) - q2 / 12 + q2^3 / 360 - q2^5 / 1260,
    0.5 * log(q2_safe) - (log(q2_safe) + 1) / q2_safe - lgamma(1 / q2_safe)
  )
  const - w^2 * k
}

# Log survival of the standardized generalized gamma
.gengamma_logsurv <- function(w, q) {
  q <- as.numeric(q)
  w <- rep_len(as.numeric(w), length(q))
  out <- stats::pnorm(w, lower.tail = FALSE, log.p = TRUE)
  pos <- q > 1e-8
  neg <- q < -1e-8
  if (any(pos)) {
    a <- 1 / q[pos]^2
    out[pos] <- stats::pgamma(a * exp(q[pos] * w[pos]), shape = a, lower.tail = FALSE, log.p = TRUE)
  }
  if (any(neg)) {
    a <- 1 / q[neg]^2
    out[neg] <- stats::pgamma(a * exp(q[neg] * w[neg]), shape = a, lower.tail = TRUE, log.p = TRUE)
  }
  out
}

# Log survival of the non-cured patients
.cure_latency_logsurv <- function(t, dist, mu, alpha, q = NULL) {
  mu <- as.numeric(mu)
  alpha <- as.numeric(alpha)
  switch(
    dist,
    weibull = -(t / exp(mu))^alpha,
    exponential = -t / exp(mu),
    lognormal = stats::plnorm(t, meanlog = mu, sdlog = 1 / alpha, lower.tail = FALSE, log.p = TRUE),
    loglogistic = stats::plogis((log(t) - mu) * alpha, lower.tail = FALSE, log.p = TRUE),
    gamma = stats::pgamma(t, shape = alpha, rate = alpha / exp(mu), lower.tail = FALSE, log.p = TRUE),
    gengamma = .gengamma_logsurv((log(t) - mu) * alpha, q),
    stop("Unknown distribution.")
  )
}

# Log density of the non-cured patients
.cure_latency_logdens <- function(t, dist, mu, alpha, q = NULL) {
  mu <- as.numeric(mu)
  alpha <- as.numeric(alpha)
  switch(
    dist,
    weibull = stats::dweibull(t, shape = alpha, scale = exp(mu), log = TRUE),
    exponential = stats::dexp(t, rate = exp(-mu), log = TRUE),
    lognormal = stats::dlnorm(t, meanlog = mu, sdlog = 1 / alpha, log = TRUE),
    loglogistic = log(alpha) - log(t) + stats::dlogis((log(t) - mu) * alpha, log = TRUE),
    gamma = stats::dgamma(t, shape = alpha, rate = alpha / exp(mu), log = TRUE),
    gengamma = log(alpha) - log(t) + .gengamma_logdens((log(t) - mu) * alpha, q),
    stop("Unknown distribution.")
  )
}

# Pointwise log-likelihood (iterations x chains x observations)
.cure_pointwise_loglik <- function(fit) {
  dist <- if (is.null(fit$dist)) "weibull" else fit$dist
  draws <- as.array(fit$stan_fit)
  dims <- dim(draws)
  pars <- dimnames(draws)$parameters
  get_par <- function(p) as.vector(draws[, , p])

  beta_cure_intercept <- get_par("beta_cure_intercept")
  beta_cure_arm <- get_par("beta_cure_arm")
  beta_surv_intercept <- get_par("beta_surv_intercept")
  beta_surv_arm <- get_par("beta_surv_arm")
  alpha_control <- get_par("alpha_control")
  beta_alpha_arm <- if ("beta_alpha_arm" %in% pars) get_par("beta_alpha_arm") else rep(0, length(alpha_control))
  q_shape <- if ("q_shape[1]" %in% pars) get_par("q_shape[1]") else NULL

  data <- fit$original_data
  cols <- fit$column_map
  time <- data[[cols$time_col]]
  event <- data[[cols$event_col]]
  arm <- as.numeric(data[[cols$arm_col]]) - 1
  n_obs <- length(time)

  log_lik <- matrix(NA_real_, nrow = length(alpha_control), ncol = n_obs)
  for (b in unique(arm)) {
    cure_logit <- beta_cure_intercept + beta_cure_arm * b
    log_cured <- stats::plogis(cure_logit, log.p = TRUE)
    log_uncured <- stats::plogis(cure_logit, lower.tail = FALSE, log.p = TRUE)
    mu <- beta_surv_intercept + beta_surv_arm * b
    alpha <- alpha_control * exp(beta_alpha_arm * b)
    for (i in which(arm == b)) {
      if (event[i] == 1) {
        log_lik[, i] <- log_uncured + .cure_latency_logdens(time[i], dist, mu, alpha, q_shape)
      } else {
        log_surv <- log_uncured + .cure_latency_logsurv(time[i], dist, mu, alpha, q_shape)
        m <- pmax(log_cured, log_surv)
        log_lik[, i] <- m + log(exp(log_cured - m) + exp(log_surv - m))
      }
    }
  }
  array(log_lik, dim = c(dims[1], dims[2], n_obs))
}

#' @title Compare Bayesian Cure Models and Recommend a Distribution
#' @description Ranks two or more models fitted with `fit_bayesian_cure_model`
#'   to the same data by their expected out-of-sample predictive accuracy and
#'   reports which distribution for the non-cured patients is preferred. Models
#'   whose predictive accuracy is practically the same as that of the best one
#'   are flagged as practical ties.
#'
#' @param ... Two or more objects of class `bcm_fit`, for example
#'   `compare_models(fit1, fit2, fit3)`, or a single list of them. Models are
#'   labelled with the argument or list names and, failing that, with their
#'   distribution.
#' @param criterion Character string. "loo" (default) uses approximate
#'   leave-one-out cross-validation with Pareto smoothed importance sampling;
#'   "waic" uses the widely applicable information criterion.
#' @param tie_elpd Numeric. A model whose ELPD is lower than that of the best
#'   model by less than this amount (default 4) is a practical tie.
#' @param tie_se Numeric. A model that is not a practical tie is "clearly
#'   worse" when the ELPD difference is at least this many standard errors
#'   (default 2) and "probably worse" otherwise.
#'
#' @details
#' Fitting is the slow step, so the comparison works on models that have
#' already been fitted, and possibly saved with `saveRDS`: nothing is refitted.
#'
#' The pointwise log-likelihood of every model is computed from its posterior
#' draws: an observed event contributes the density of the non-cured patients
#' weighted by the probability of not being cured, and a censored observation
#' contributes the mixture of the cure probability and the survival function of
#' the non-cured patients. The expected log pointwise predictive density (ELPD)
#' is then estimated with the `loo` package (Vehtari, Gelman and Gabry, 2017).
#' Stacking weights (Yao et al., 2018) give the weight of each model in the
#' combination with the best predictive performance. With `criterion = "loo"`,
#' Pareto k values above 0.7 flag models for which the approximation is
#' unreliable.
#'
#' Each model receives a verdict:
#' \itemize{
#'   \item "best": the highest ELPD.
#'   \item "practical tie": the ELPD difference with the best model is below
#'     `tie_elpd`. Differences of this size are small: the models make very
#'     similar predictions and choosing one or the other has little
#'     consequence (Sivula et al., 2020).
#'   \item "probably worse": the difference is at least `tie_elpd` but smaller
#'     than `tie_se` standard errors, so it is not conclusive.
#'   \item "clearly worse": the difference is at least `tie_elpd` and at least
#'     `tie_se` standard errors.
#' }
#' Any threshold is a convention. Among practical ties the choice should rest
#' on other grounds, such as the pre-specified model or parsimony, and the
#' conclusions should be checked to agree across the tied models.
#'
#' @return An object of class `bcm_comparison`, whose print method shows the
#'   report. It is a list with elements:
#'   \itemize{
#'     \item{\code{comparison}: a data frame with one row per model, ordered from
#'       best to worst, with the difference in ELPD with respect to the best
#'       model and its standard error, the ELPD, the effective number of
#'       parameters, the information criterion, the stacking weight, for
#'       "loo" the largest Pareto k value, and the verdict.}
#'     \item{\code{best_model}: the name of the model with the highest ELPD.}
#'     \item{\code{practical_ties}: the names of the models that are practical
#'       ties with the best one.}
#'     \item{\code{recommendation}: the text of the recommendation.}
#'     \item{\code{criterion}, \code{tie_elpd}, \code{tie_se}: the settings used.}
#'     \item{\code{n_obs}: the number of observations.}
#'     \item{\code{estimates}: the list of objects returned by the `loo` package
#'       for each model.}
#'   }
#'
#' @references
#' Vehtari A, Gelman A, Gabry J (2017). Practical Bayesian model evaluation
#' using leave-one-out cross-validation and WAIC. Statistics and Computing,
#' 27, 1413-1432.
#'
#' Yao Y, Vehtari A, Simpson D, Gelman A (2018). Using stacking to average
#' Bayesian predictive distributions. Bayesian Analysis, 13, 917-1007.
#'
#' Sivula T, Magnusson M, Matamoros AA, Vehtari A (2020). Uncertainty in
#' Bayesian leave-one-out cross-validation based model comparison.
#' arXiv:2008.10296.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sim_data <- simulate_weibull_cure_data(
#'   n_patients = 300, cure_fraction_ctrl = 0.15, cure_fraction_exp = 0.40,
#'   max_follow_up = 36, weibull_shape = 1.2, median_survival_ctrl = 12,
#'   time_ratio_exp = 1.25
#' )
#' fit1 <- fit_bayesian_cure_model(sim_data, dist = "weibull")
#' fit2 <- fit_bayesian_cure_model(sim_data, dist = "lognormal")
#' fit3 <- fit_bayesian_cure_model(sim_data, dist = "loglogistic")
#' comparison <- compare_models(fit1, fit2, fit3)
#' comparison
#' comparison$comparison
#' }
compare_models <- function(..., criterion = c("loo", "waic"), tie_elpd = 4, tie_se = 2) {
  criterion <- match.arg(criterion)

  # Validate inputs
  fits <- list(...)
  arg_exprs <- as.list(substitute(list(...)))[-1]
  if (length(fits) == 1 && !inherits(fits[[1]], "bcm_fit") && is.list(fits[[1]])) {
    fits <- fits[[1]]
    arg_exprs <- list()
  }
  if (length(fits) < 2) {
    stop("At least two fitted models are required.")
  }
  if (!all(vapply(fits, inherits, logical(1), what = "bcm_fit"))) {
    stop("All inputs must be objects of class 'bcm_fit'.")
  }

  # Label models
  dists <- vapply(fits, function(f) if (is.null(f$dist)) "weibull" else f$dist, character(1))
  model_names <- names(fits)
  if (is.null(model_names)) {
    model_names <- rep("", length(fits))
  }
  for (i in seq_along(arg_exprs)) {
    if (model_names[i] == "" && is.symbol(arg_exprs[[i]])) {
      model_names[i] <- as.character(arg_exprs[[i]])
    }
  }
  model_names[model_names == ""] <- dists[model_names == ""]
  model_names <- make.unique(model_names)
  names(fits) <- model_names
  names(dists) <- model_names

  # Check that all models share the same data
  model_data <- function(f) {
    cols <- f$column_map
    list(
      as.numeric(f$original_data[[cols$time_col]]),
      as.numeric(f$original_data[[cols$event_col]]),
      as.numeric(f$original_data[[cols$arm_col]])
    )
  }
  reference_data <- model_data(fits[[1]])
  for (f in fits[-1]) {
    if (!isTRUE(all.equal(reference_data, model_data(f)))) {
      stop("All models must be fitted to the same data.")
    }
  }

  # Estimate ELPD
  estimates <- lapply(fits, function(f) {
    log_lik <- .cure_pointwise_loglik(f)
    if (criterion == "loo") {
      r_eff <- loo::relative_eff(exp(log_lik))
      loo::loo(log_lik, r_eff = r_eff)
    } else {
      loo::waic(log_lik)
    }
  })

  # Rank models and compute stacking weights
  ranking <- loo::loo_compare(estimates)
  ordered_names <- rownames(ranking)
  n_obs <- nrow(estimates[[1]]$pointwise)
  pointwise_elpd <- vapply(estimates, function(e) e$pointwise[, 1], numeric(n_obs))
  weights <- as.numeric(loo::stacking_weights(pointwise_elpd))
  names(weights) <- model_names

  comparison <- data.frame(
    model = ordered_names,
    dist = unname(dists[ordered_names]),
    as.data.frame(unclass(ranking))[, c(1, 2, 3, 4, 5, 7, 8), drop = FALSE],
    stacking_weight = unname(weights[ordered_names]),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  if (criterion == "loo") {
    comparison$max_pareto_k <- vapply(
      estimates[ordered_names],
      function(e) max(e$diagnostics$pareto_k),
      numeric(1)
    )
  }

  # Classify models
  gap <- abs(comparison$elpd_diff)
  comparison$verdict <- ifelse(
    gap < tie_elpd,
    "practical tie",
    ifelse(gap < tie_se * comparison$se_diff, "probably worse", "clearly worse")
  )
  comparison$verdict[1] <- "best"

  # Write recommendation
  label <- function(m) ifelse(m == dists[m], m, paste0(m, " (", dists[m], ")"))
  ties <- comparison$model[comparison$verdict == "practical tie"]
  inconclusive <- comparison$model[comparison$verdict == "probably worse"]
  if (length(ties) == 0) {
    recommendation <- paste0("Recommended model: ", label(ordered_names[1]), ".")
  } else {
    recommendation <- paste0(
      "Highest predictive accuracy: ", label(ordered_names[1]), ". Practical ",
      if (length(ties) == 1) "tie: " else "ties: ", paste(label(ties), collapse = ", "),
      ". These models predict the data equally well for practical purposes: ",
      "choose among them on other grounds (pre-specified model, parsimony) ",
      "and check that the conclusions agree."
    )
  }
  if (length(inconclusive) > 0) {
    recommendation <- c(recommendation, paste0(
      "The evidence against ", paste(label(inconclusive), collapse = ", "),
      " is not conclusive (difference within ", tie_se, " standard errors)."
    ))
  }
  if (criterion == "loo" && any(comparison$max_pareto_k > 0.7)) {
    recommendation <- c(recommendation, paste0(
      "Pareto k above 0.7 for ",
      paste(label(comparison$model[comparison$max_pareto_k > 0.7]), collapse = ", "),
      ": the LOO estimate of these models may be unreliable."
    ))
  }
  if (n_obs < 100) {
    recommendation <- c(
      recommendation,
      "Fewer than 100 observations: the standard errors of the differences may be unreliable."
    )
  }

  # Assemble result
  result <- list(
    comparison = comparison,
    best_model = ordered_names[1],
    practical_ties = ties,
    recommendation = recommendation,
    criterion = criterion,
    tie_elpd = tie_elpd,
    tie_se = tie_se,
    n_obs = n_obs,
    estimates = estimates
  )
  class(result) <- "bcm_comparison"
  result
}

#' @rdname compare_models
#' @param x An object of class `bcm_comparison`.
#' @param digits Number of decimal places shown in the report.
#' @export
print.bcm_comparison <- function(x, digits = 1, ...) {
  tab <- x$comparison
  cat(
    "Comparison of Bayesian cure models (",
    if (x$criterion == "loo") "PSIS-LOO" else "WAIC", ")\n",
    x$n_obs, " observations, ", nrow(tab), " models\n\n",
    sep = ""
  )
  report <- data.frame(
    model = tab$model,
    dist = tab$dist,
    elpd_diff = formatC(tab$elpd_diff, digits = digits, format = "f"),
    se_diff = formatC(tab$se_diff, digits = digits, format = "f"),
    verdict = tab$verdict,
    stringsAsFactors = FALSE
  )
  print(report, row.names = FALSE)
  cat("\n")
  rule <- paste0(
    "Rule: practical tie if the ELPD difference is below ", x$tie_elpd,
    "; beyond that, clearly worse if the difference is at least ", x$tie_se,
    " standard errors and probably worse otherwise."
  )
  for (paragraph in x$recommendation) {
    cat(strwrap(paragraph, width = 80), sep = "\n")
  }
  cat("\n")
  cat(strwrap(rule, width = 80), sep = "\n")
  invisible(x)
}
