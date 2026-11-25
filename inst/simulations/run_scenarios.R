# Simulation scenarios of the BayeScores manuscript (Scenarios 1 to 12) fitted
# with the default Weibull mixture cure model of bayescores.
#
# Usage: Rscript run_scenarios.R [replicates] [output_dir] [scenarios] [worker]
#   replicates  number of replicates per scenario (default 500)
#   output_dir  output folder (default "results")
#   scenarios   "all" (default) or names separated by commas, e.g.
#               sc01_benefit,sc02_null
#   worker      "k/n" runs the replicates with (replicate %% n) == k - 1, so that
#               n processes can share the work (default "1/1")
#
# One file is saved per replicate, so an interrupted run can be resumed. When
# all requested replicates exist, the per-replicate results and the summary are
# written to the output folder. Scenario 5 (5,000 patients) takes about 4
# minutes per replicate; the scenarios with 200 patients about 5 seconds.

args <- commandArgs(TRUE)
replicates <- if (length(args) >= 1) as.integer(args[1]) else 500L
out_dir <- if (length(args) >= 2) args[2] else "results"
which_sc <- if (length(args) >= 3) args[3] else "all"
worker <- if (length(args) >= 4) as.integer(strsplit(args[4], "/")[[1]]) else c(1L, 1L)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

library(rstan)
source("simulation_functions.R")
source("model_fit.R")
source("summary_functions.R")

run <- if (which_sc == "all") names(scenarios) else strsplit(which_sc, ",")[[1]]
stopifnot(all(run %in% names(scenarios)))
rep_file <- function(name, r) file.path(out_dir, name, sprintf("rep_%04d.rds", r))

sm <- compile_model(system.file("stan", "cure_model.stan", package = "bayescores"))

for (name in run) {
  dir.create(file.path(out_dir, name), showWarnings = FALSE)
  for (r in seq_len(replicates)) {
    if (r %% worker[2] != worker[1] - 1 || file.exists(rep_file(name, r))) next
    dat <- generate_trial(scenarios[[name]], seed = r)
    fit <- suppressWarnings(fit_trial(sm, dat, seed = r))
    res <- cbind(scenario = name, replicate = r, fit,
                 n_events = sum(dat$event), prop_censored = mean(dat$event == 0))
    tmp <- paste0(rep_file(name, r), ".tmp")
    saveRDS(res, tmp); file.rename(tmp, rep_file(name, r))
  }
  message(name, " done")
}

files <- unlist(lapply(run, function(name) rep_file(name, seq_len(replicates))))
if (all(file.exists(files))) {
  res <- do.call(rbind, lapply(files, readRDS))
  summary_table <- summarise_results(res)
  utils::write.csv(res, file.path(out_dir, "per_replicate_results.csv"), row.names = FALSE)
  utils::write.csv(summary_table, file.path(out_dir, "summary_by_scenario.csv"), row.names = FALSE)
}
