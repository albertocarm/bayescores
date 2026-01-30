# Package index

## Model Fitting & Diagnostics

Fit Bayesian AFT mixture cure models and assess convergence

- [`fit_bayesian_cure_model()`](https://albertocarm.github.io/bayescores/reference/fit_bayesian_cure_model.md)
  : Fit the Bayesian Cure Model (Engine)
- [`outcomes()`](https://albertocarm.github.io/bayescores/reference/outcomes.md)
  : Generate a Summary Table of Model Outcomes (Final Version 2.1)
- [`diagnose_fit()`](https://albertocarm.github.io/bayescores/reference/diagnose_fit.md)
  : Diagnose Bayesian MCMC fit
- [`model_diagnostics()`](https://albertocarm.github.io/bayescores/reference/model_diagnostics.md)
  : Model Diagnostics
- [`check_false_plateau()`](https://albertocarm.github.io/bayescores/reference/check_false_plateau.md)
  : Check for False Plateau or Structural Instability

## BayeScores Calculation

Compute and summarize Bayesian clinical benefit scores

- [`get_bayescores()`](https://albertocarm.github.io/bayescores/reference/get_bayescores.md)
  : Calculate BayeScores Final Utility Score (Version 4, Corrected
  Toxicity Model)
- [`get_bayescores_draws()`](https://albertocarm.github.io/bayescores/reference/get_bayescores_draws.md)
  : Get Original or Shrunk MCMC Draws for Subsequent Analyses (Updated)
- [`summarize_final_utility()`](https://albertocarm.github.io/bayescores/reference/summarize_final_utility.md)
  : Summarizes the final utility score with advanced dynamic weighting.
- [`summarize_cure_rates()`](https://albertocarm.github.io/bayescores/reference/summarize_cure_rates.md)
  : Summarize Cure Rate Results
- [`time_ratio()`](https://albertocarm.github.io/bayescores/reference/time_ratio.md)
  : Calculate the Time Ratio
- [`calculate_asco_nhb()`](https://albertocarm.github.io/bayescores/reference/calculate_asco_nhb.md)
  : Calculate the ASCO Net Health Benefit (NHB) Score
- [`calculate_esmo_mcbs()`](https://albertocarm.github.io/bayescores/reference/calculate_esmo_mcbs.md)
  : Calculate the ESMO-MCBS v2.0 Score (Guided by Objective)
- [`print(`*`<bayes_summary>`*`)`](https://albertocarm.github.io/bayescores/reference/print.bayes_summary.md)
  : Print method for bayes_summary objects

## Posterior Extraction

Extract MCMC posterior samples

- [`extract_mcmc_time_ratios()`](https://albertocarm.github.io/bayescores/reference/extract_mcmc_time_ratios.md)
  : Extract MCMC Samples for the Time Ratio
- [`extract_mcmc_cure_diffs()`](https://albertocarm.github.io/bayescores/reference/extract_mcmc_cure_diffs.md)
  : Extract MCMC Samples for the Difference in Cure Rates

## Toxicity Analysis

Weighted toxicity scoring and visualization

- [`calculate_toxicity_analysis()`](https://albertocarm.github.io/bayescores/reference/calculate_toxicity_analysis.md)
  : Calculate a complete toxicity analysis
- [`create_amit_plot()`](https://albertocarm.github.io/bayescores/reference/create_amit_plot.md)
  : Create an AMIT plot (Adverse Event Dot Plot) from a clinical trial
  object.

## Quality of Life

Sample and visualize QoL distributions

- [`sample_qol_scores()`](https://albertocarm.github.io/bayescores/reference/sample_qol_scores.md)
  : Sample Quality of Life (QoL) Scores from a Multinomial Distribution
- [`generate_qol_vector()`](https://albertocarm.github.io/bayescores/reference/generate_qol_vector.md)
  : Interactively Generate a Quality of Life (QoL) Probability Vector

## Visualization

Plotting functions for model results and clinical scores

- [`plot_densities()`](https://albertocarm.github.io/bayescores/reference/plot_densities.md)
  : Plot Posterior Densities
- [`plot_correlated_densities()`](https://albertocarm.github.io/bayescores/reference/plot_correlated_densities.md)
  : Plot Correlated Posterior Densities
- [`plot_utility_donut()`](https://albertocarm.github.io/bayescores/reference/plot_utility_donut.md)
  : Plot Utility Score Gauge Visualization
- [`plot_final_utility_density()`](https://albertocarm.github.io/bayescores/reference/plot_final_utility_density.md)
  : Plot Utility Score Density Distribution
- [`plot_km_curves()`](https://albertocarm.github.io/bayescores/reference/plot_km_curves.md)
  : Plot Kaplan-Meier Survival Curves
- [`plot_qol_histogram()`](https://albertocarm.github.io/bayescores/reference/plot_qol_histogram.md)
  : Plot an Elegant Histogram of QoL Outcomes
- [`plot_toxicity_density()`](https://albertocarm.github.io/bayescores/reference/plot_toxicity_density.md)
  : Plot Toxicity Score Density
- [`generate_sensitivity_dashboard()`](https://albertocarm.github.io/bayescores/reference/generate_sensitivity_dashboard.md)
  : Generate Comprehensive Sensitivity Analysis Dashboard
- [`run_plot_scenario()`](https://albertocarm.github.io/bayescores/reference/run_plot_scenario.md)
  : Run a Scenario for a Sensitivity Plot

## Simulation & Data Reconstruction

Generate synthetic trial data and reconstruct IPD from published curves

- [`simulate_weibull_cure_data()`](https://albertocarm.github.io/bayescores/reference/simulate_weibull_cure_data.md)
  : Simulate Survival Data from a Cure Model
- [`simulate_trial_data()`](https://albertocarm.github.io/bayescores/reference/simulate_trial_data.md)
  : Interactively or non-interactively simulate a realistic clinical
  trial dataset
- [`reconstruct_ipd()`](https://albertocarm.github.io/bayescores/reference/reconstruct_ipd.md)
  : Reconstruct Individual Patient Data from Multiple Kaplan-Meier
  Curves

## Included Datasets

Pre-loaded RCT survival and toxicity data

- [`AURA3_1A`](https://albertocarm.github.io/bayescores/reference/AURA3_1A.md)
  : Data from the AURA3_1A Trial
- [`BOLERO3_2`](https://albertocarm.github.io/bayescores/reference/BOLERO3_2.md)
  : Data from the BOLERO3_2 Trial
- [`CONCUR_2A`](https://albertocarm.github.io/bayescores/reference/CONCUR_2A.md)
  : Data from the CONCUR_2A Trial
- [`COUAA202_2`](https://albertocarm.github.io/bayescores/reference/COUAA202_2.md)
  : Data from the COUAA202_2 Trial
- [`DART0105_2A`](https://albertocarm.github.io/bayescores/reference/DART0105_2A.md)
  : Data from the DART0105_2A Trial
- [`GEICAM2003_2A`](https://albertocarm.github.io/bayescores/reference/GEICAM2003_2A.md)
  : Data from the GEICAM2003_2A Trial
- [`GETUG12_2`](https://albertocarm.github.io/bayescores/reference/GETUG12_2.md)
  : Data from the GETUG12_2 Trial
- [`IMELDA_2A`](https://albertocarm.github.io/bayescores/reference/IMELDA_2A.md)
  : Data from the IMELDA_2A Trial
- [`KEYNOTE024_1A`](https://albertocarm.github.io/bayescores/reference/KEYNOTE024_1A.md)
  : Data from the KEYNOTE024_1A Trial
- [`KEYNOTE062`](https://albertocarm.github.io/bayescores/reference/KEYNOTE062.md)
  : Data from the KEYNOTE062 Trial
- [`LUXLung8_2A`](https://albertocarm.github.io/bayescores/reference/LUXLung8_2A.md)
  : Data from the LUXLung8_2A Trial
- [`MA17R_1A`](https://albertocarm.github.io/bayescores/reference/MA17R_1A.md)
  : Data from the MA17R_1A Trial
- [`MONALEESA2_1`](https://albertocarm.github.io/bayescores/reference/MONALEESA2_1.md)
  : Data from the MONALEESA2_1 Trial
- [`NO16968_2A`](https://albertocarm.github.io/bayescores/reference/NO16968_2A.md)
  : Data from the NO16968_2A Trial
- [`NO16968_2C`](https://albertocarm.github.io/bayescores/reference/NO16968_2C.md)
  : Data from the NO16968_2C Trial
- [`NOAH_2A`](https://albertocarm.github.io/bayescores/reference/NOAH_2A.md)
  : Data from the NOAH_2A Trial
- [`NOAH_2B`](https://albertocarm.github.io/bayescores/reference/NOAH_2B.md)
  : Data from the NOAH_2B Trial
- [`NSABPB35_2`](https://albertocarm.github.io/bayescores/reference/NSABPB35_2.md)
  : Data from the NSABPB35_2 Trial
- [`NSABPB35_3`](https://albertocarm.github.io/bayescores/reference/NSABPB35_3.md)
  : Data from the NSABPB35_3 Trial
- [`NSABPB40_3A`](https://albertocarm.github.io/bayescores/reference/NSABPB40_3A.md)
  : Data from the NSABPB40_3A Trial
- [`NSABPB40_3B`](https://albertocarm.github.io/bayescores/reference/NSABPB40_3B.md)
  : Data from the NSABPB40_3B Trial
- [`PALOMA2_1A`](https://albertocarm.github.io/bayescores/reference/PALOMA2_1A.md)
  : Data from the PALOMA2_1A Trial
- [`PROFILE1014_1A`](https://albertocarm.github.io/bayescores/reference/PROFILE1014_1A.md)
  : Data from the PROFILE1014_1A Trial
- [`RECOURSE_1A`](https://albertocarm.github.io/bayescores/reference/RECOURSE_1A.md)
  : Data from the RECOURSE_1A Trial
- [`aura3_toxicity`](https://albertocarm.github.io/bayescores/reference/aura3_toxicity.md)
  : Toxicity and quality of life data from the AURA3 trial
- [`bolero3_toxicity`](https://albertocarm.github.io/bayescores/reference/bolero3_toxicity.md)
  : Toxicity and quality of life data from the BOLERO3 trial
- [`checkmate649_toxicity`](https://albertocarm.github.io/bayescores/reference/checkmate649_toxicity.md)
  : Toxicity and quality of life data from the CheckMate-649 trial
- [`cm649`](https://albertocarm.github.io/bayescores/reference/cm649.md)
  : Data from the CheckMate-649 Trial
- [`cm649_2024_update`](https://albertocarm.github.io/bayescores/reference/cm649_2024_update.md)
  : Data from the CheckMate-649 Trial -2024 update
- [`concur_toxicity`](https://albertocarm.github.io/bayescores/reference/concur_toxicity.md)
  : Toxicity and quality of life data from the CONCUR trial
- [`dart0105_toxicity`](https://albertocarm.github.io/bayescores/reference/dart0105_toxicity.md)
  : Toxicity and quality of life data from the DART0105 trial
- [`esopec`](https://albertocarm.github.io/bayescores/reference/esopec.md)
  : Data from the esopec Trial
- [`geicam2003_toxicity`](https://albertocarm.github.io/bayescores/reference/geicam2003_toxicity.md)
  : Toxicity and quality of life data from the GEICAM 2003 trial
- [`getug12_toxicity`](https://albertocarm.github.io/bayescores/reference/getug12_toxicity.md)
  : Toxicity and quality of life data from the GETUG12 trial
- [`glow`](https://albertocarm.github.io/bayescores/reference/glow.md) :
  Data from the glow Trial
- [`glow_toxicity`](https://albertocarm.github.io/bayescores/reference/glow_toxicity.md)
  : Toxicity and quality of life data from the GLOW trial
- [`imbrave150_toxicity`](https://albertocarm.github.io/bayescores/reference/imbrave150_toxicity.md)
  : Toxicity and quality of life data from the IMBRAVE150 trial
- [`imelda_toxicity`](https://albertocarm.github.io/bayescores/reference/imelda_toxicity.md)
  : Toxicity and quality of life data from the IMELDA trial
- [`keynote024_toxicity`](https://albertocarm.github.io/bayescores/reference/keynote024_toxicity.md)
  : Toxicity and quality of life data from the KEYNOTE024 trial
- [`luxlung3_6_toxicity`](https://albertocarm.github.io/bayescores/reference/luxlung3_6_toxicity.md)
  : Toxicity and quality of life data from the LUX‑Lung 3 and 6 trials
- [`luxlung8_toxicity`](https://albertocarm.github.io/bayescores/reference/luxlung8_toxicity.md)
  : Toxicity and quality of life data from the LUX‑Lung 8 trial
- [`ma17r_toxicity`](https://albertocarm.github.io/bayescores/reference/ma17r_toxicity.md)
  : Toxicity and quality of life data from the MA17R trial
- [`monaleesa2_toxicity`](https://albertocarm.github.io/bayescores/reference/monaleesa2_toxicity.md)
  : Toxicity and quality of life data from the MONALEESA‑2 trial
- [`monarchE`](https://albertocarm.github.io/bayescores/reference/monarchE.md)
  : Data from the monarchE Trial
- [`natalee`](https://albertocarm.github.io/bayescores/reference/natalee.md)
  : Data from the natalee Trial
- [`no16968_toxicity`](https://albertocarm.github.io/bayescores/reference/no16968_toxicity.md)
  : Toxicity and quality of life data from the NO16968 trial
- [`noah_toxicity`](https://albertocarm.github.io/bayescores/reference/noah_toxicity.md)
  : Toxicity and quality of life data from the NOAH trial
- [`nsabpb35_toxicity`](https://albertocarm.github.io/bayescores/reference/nsabpb35_toxicity.md)
  : Toxicity and quality of life data from the NSABP B‑35 trial
- [`nsabpb40_toxicity`](https://albertocarm.github.io/bayescores/reference/nsabpb40_toxicity.md)
  : Toxicity and quality of life data from the NSABP B‑40 trial
- [`paloma2_toxicity`](https://albertocarm.github.io/bayescores/reference/paloma2_toxicity.md)
  : Toxicity and quality of life data from the PALOMA2 trial
- [`profile1014_toxicity`](https://albertocarm.github.io/bayescores/reference/profile1014_toxicity.md)
  : Toxicity and quality of life data from the PROFILE1014 trial
- [`recourse_toxicity`](https://albertocarm.github.io/bayescores/reference/recourse_toxicity.md)
  : Toxicity and quality of life data from the RECOURSE trial
- [`spotlight`](https://albertocarm.github.io/bayescores/reference/spotlight.md)
  : Data from the spotlight Trial
- [`spotlight_toxicity`](https://albertocarm.github.io/bayescores/reference/spotlight_toxicity.md)
  : Toxicity and quality of life data from the SPOTLIGHT trial
- [`toga`](https://albertocarm.github.io/bayescores/reference/toga.md) :
  Data from the toga Trial
- [`toga_toxicity`](https://albertocarm.github.io/bayescores/reference/toga_toxicity.md)
  : Toxicity and quality of life data from the TOGA trial
- [`Keynote062_toxicity`](https://albertocarm.github.io/bayescores/reference/Keynote062_toxicity.md)
  : Toxicity and QoL data from Keynote062_toxicity
- [`esopec_toxicity`](https://albertocarm.github.io/bayescores/reference/esopec_toxicity.md)
  : Toxicity and QoL data from esopec_toxicity
- [`keynote_119_toxicity`](https://albertocarm.github.io/bayescores/reference/keynote_119_toxicity.md)
  : Toxicity and QoL data from KEYNOTE-119
- [`monarchE_toxicity`](https://albertocarm.github.io/bayescores/reference/monarchE_toxicity.md)
  : Toxicity and QoL data from monarchE_toxicity
- [`natalee_toxicity`](https://albertocarm.github.io/bayescores/reference/natalee_toxicity.md)
  : Toxicity and QoL data from natalee_toxicity
- [`kn119`](https://albertocarm.github.io/bayescores/reference/kn119.md)
  : Digitized survival data from KEYNOTE-119
