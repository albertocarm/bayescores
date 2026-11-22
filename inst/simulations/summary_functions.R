# Performance metrics per scenario and parameter, as in Table 1.
summarise_results <- function(res) {
  keys <- unique(res[, "scenario", drop = FALSE])
  out <- lapply(seq_len(nrow(keys)), function(i) {
    d <- res[res$scenario == keys$scenario[i], ]
    tv <- true_values(scenarios[[keys$scenario[i]]])
    do.call(rbind, lapply(pars, function(p) {
      x <- d[d$param == p, ]; truth <- tv[[p]]
      excl0 <- mean(x$ci_low > 0 | x$ci_high < 0)
      data.frame(
        scenario = keys$scenario[i], param = p,
        replicates = nrow(x), true_value = truth,
        mean_estimate = mean(x$estimate),
        bias = mean(x$estimate) - truth,
        coverage = mean(x$ci_low <= truth & x$ci_high >= truth),
        mse = mean((x$estimate - truth)^2),
        empirical_se = stats::sd(x$estimate),
        mean_model_se = mean(x$sd_model),
        excludes_zero = excl0)
    }))
  })
  do.call(rbind, out)
}
