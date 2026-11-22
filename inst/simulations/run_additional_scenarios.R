# Additional simulation scenarios (non-Weibull survival times, delayed
# treatment effect, alternative censoring) fitted with the default Weibull
# mixture cure model of bayescores 1.2.
#
# Usage: Rscript run_additional_scenarios.R [replicates] [output_dir]
# Requires bayescores 1.2 and rstan. Results are saved per scenario, so an
# interrupted run can be resumed.

args <- commandArgs(TRUE)
replicates <- if (length(args) >= 1) as.integer(args[1]) else 500L
out_dir <- if (length(args) >= 2) args[2] else "results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

library(rstan)
source("simulation_functions.R")
source("model_fit.R")
source("summary_functions.R")

sm <- compile_model(system.file("stan", "cure_model.stan", package = "bayescores"))

for (name in names(scenarios)) {
  f <- file.path(out_dir, paste0(name, ".rds"))
  if (file.exists(f)) next
  res <- lapply(seq_len(replicates), function(r) {
    dat <- generate_trial(scenarios[[name]], seed = r)
    fit <- suppressWarnings(fit_trial(sm, dat, seed = r))
    cbind(scenario = name, replicate = r, fit,
          n_events = sum(dat$event), prop_censored = mean(dat$event == 0))
  })
  saveRDS(do.call(rbind, res), f)
  message(name, " done")
}

res <- do.call(rbind, lapply(file.path(out_dir, paste0(names(scenarios), ".rds")), readRDS))
summary_table <- summarise_results(res)
utils::write.csv(res, file.path(out_dir, "per_replicate_results.csv"), row.names = FALSE)
utils::write.csv(summary_table, file.path(out_dir, "summary_by_scenario.csv"), row.names = FALSE)
