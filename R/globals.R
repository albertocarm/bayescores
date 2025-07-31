# This tells R CMD check that these variables are used intentionally
# inside dplyr/ggplot2 pipes and are not global variables.
utils::globalVariables(c("x", "y", "value", "ymax_pos", "ymin_pos", "component", "TRT", "AdjustmentFactor", "Analysis", "Cure", "Median_Utility", "QoL", "St", "TR", "Toxicity", "category",
  "curve", "effect", "end", "end_id", "get_bayescores_final", "id",
  "indices", "nrisk", "start", "start_id", "time", "conf.low", "conf.high", "strata", "estimate", "survival", "arm"))

