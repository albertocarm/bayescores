# The Stan program of the package is compiled once and sampled with the data
# list and settings of bayescores::fit_bayesian_cure_model() (neutral prior,
# common shape). Both give identical draws for the same data and seed.
compile_model <- function(stan_file) rstan::stan_model(file = stan_file, auto_write = FALSE)

fit_trial <- function(sm, dat, seed, shared_shape = TRUE,
                      chains = 2, iter = 2000, warmup = 1000) {
  set.seed(100000 + seed)
  stan_data <- list(
    N = nrow(dat), tiempo = dat$time, evento = dat$event,
    arm = as.numeric(dat$arm) - 1,
    cure_prior_type = 1L,
    use_historical_prior = 0L,
    historical_prior_params = c(0, 1),
    shared_shape = as.integer(shared_shape))
  init_fun <- function() list(beta_cure_arm_raw = stats::rnorm(1, 0, 0.1))
  t0 <- Sys.time()
  sf <- rstan::sampling(sm, data = stan_data, init = init_fun, chains = chains,
                        iter = iter, warmup = warmup, seed = seed,
                        control = list(adapt_delta = 0.99), refresh = 0, cores = 1)
  sm_ <- rstan::summary(sf, pars = pars)$summary
  sp <- rstan::get_sampler_params(sf, inc_warmup = FALSE)
  dr <- rstan::extract(sf, pars = c("beta_cure_arm", "beta_surv_arm"))
  data.frame(
    param = rownames(sm_), estimate = sm_[, "mean"], sd_model = sm_[, "sd"],
    ci_low = sm_[, "2.5%"], ci_high = sm_[, "97.5%"], rhat = sm_[, "Rhat"],
    n_eff = sm_[, "n_eff"],
    divergent = sum(vapply(sp, function(x) sum(x[, "divergent__"]), numeric(1))),
    seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")),
    # Posterior median, posterior correlation between the two treatment effects
    # (identifiability diagnostic) and transitions at the maximum treedepth.
    median = sm_[, "50%"],
    r_pearson = stats::cor(dr$beta_cure_arm, dr$beta_surv_arm),
    r_spearman = stats::cor(dr$beta_cure_arm, dr$beta_surv_arm, method = "spearman"),
    treedepth_hits = sum(vapply(sp, function(x) sum(x[, "treedepth__"] >= 10), numeric(1))),
    row.names = NULL)
}
