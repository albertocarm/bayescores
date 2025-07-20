# This tells R CMD check that these variables are used intentionally
# inside dplyr/ggplot2 pipes and are not global variables.
utils::globalVariables(c("x", "y", "value", "ymax_pos", "ymin_pos", "component", "time", "conf.low", "conf.high", "strata", "estimate", "survival", "arm"))
