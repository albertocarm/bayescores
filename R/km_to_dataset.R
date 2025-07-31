
# Internal Processing Function
# -------------------------------------------------------------------------

#' Process a single survival curve and generate IPD
#'
#' This is an internal helper function that takes the digitized points of a single
#' Kaplan-Meier curve, cleans the data, handles intervals with no data points,
#' and uses the `survHE::digitise` function to generate
#' Individual Patient Data (IPD). It is not meant to be called directly by the end-user.
#'
#' @param curve_df A data frame containing the digitized points for a single curve.
#'   Must include 'time' and 'St' (survival probability) columns.
#' @param n_risk A numeric vector representing the number of patients at risk at
#'   the start of each time interval defined by `time_breaks`.
#' @param time_breaks A numeric vector defining the start and end times for the
#'   intervals of the number at risk table.
#' @param arm_name A character string to name the treatment arm (e.g., "control").
#' @param output_directory A character string specifying the path to the directory
#'   where intermediate files (`KMdata` and `IPDdata`) will be saved.
#' @param digitise_function The function to be used for IPD reconstruction.
#'   Defaults to `survHE::digitise`.
#'
#' @importFrom dplyr select filter arrange mutate row_number slice_tail pull bind_rows rowwise ungroup rename
#' @importFrom utils write.table read.delim
#' @importFrom stats na.omit
#'
#' @return A data frame containing the reconstructed Individual Patient Data (IPD)
#'   for the given arm, with 'time', 'event', and 'arm' columns.
#' @noRd
process_curve_and_generate_ipd <- function(
    curve_df,
    n_risk,
    time_breaks,
    arm_name,
    output_directory,
    digitise_function
) {

  cat(paste("\n--- Starting processing for arm:", arm_name, "---\n"))

  # --- STEP A: Cleaning duplicates (ORIGINAL LOGIC RESTORED) ---
  clean_df <- curve_df %>% dplyr::select(-curve)

  # 1. Correct duplicates in 'time'
  clean_df <- clean_df[order(clean_df$time), ]
  for (i in 2:nrow(clean_df)) {
    if (clean_df$time[i] <= clean_df$time[i - 1]) {
      clean_df$time[i] <- clean_df$time[i - 1] + 0.001
    }
  }

  # 2. Correct non-decreasing values in 'St'
  clean_df <- clean_df[order(-clean_df$St), ]
  for (i in 2:nrow(clean_df)) {
    if (clean_df$St[i] >= clean_df$St[i - 1]) {
      clean_df$St[i] <- clean_df$St[i - 1] - 0.001
    }
  }

  # Re-sort by time
  clean_df <- clean_df[order(clean_df$time), ]


  # --- STEP B: Add points in time intervals with no data ---
  intervals <- data.frame(start = time_breaks[-length(time_breaks)], end = time_breaks[-1])
  new_rows <- list()

  for (i in 1:nrow(intervals)) {
    a <- intervals$start[i]
    b <- intervals$end[i]
    if (!any(clean_df$time >= a & clean_df$time < b)) {
      cat(paste("     -> Empty interval detected [", a, ",", b, "). Adding points.\n"))
      previous_st <- clean_df %>% dplyr::filter(time < a) %>% dplyr::slice_tail(n = 1) %>% dplyr::pull(St)
      if (length(previous_st) == 0) previous_st <- 1.0

      new_rows[[length(new_rows) + 1]] <- data.frame(
        time = a + 0.25 * (b - a),
        St = previous_st - 0.001
      )
      new_rows[[length(new_rows) + 1]] <- data.frame(
        time = a + 0.75 * (b - a),
        St = previous_st - 0.002
      )
    }
  }

  if (length(new_rows) > 0) {
    clean_df <- dplyr::bind_rows(clean_df, dplyr::bind_rows(new_rows))
  }

  # --- STEP C: Prepare data for 'digitise' ---
  surv_data <- clean_df %>%
    dplyr::arrange(time) %>%
    dplyr::mutate(id = dplyr::row_number())

  nrisk_data <- intervals %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      indices = list(which(surv_data$time >= start & surv_data$time < end))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      start_id = sapply(indices, function(x) if(length(x) > 0) min(x) else NA_integer_),
      end_id = sapply(indices, function(x) if(length(x) > 0) max(x) else NA_integer_),
      id = dplyr::row_number(),
      time = end,
      nrisk = n_risk
    ) %>%
    dplyr::select(id, time, start_id, end_id, nrisk)

  # --- STEP D: Format and call the 'digitise' function ---
  formatted_surv_data <- surv_data %>% dplyr::rename(Id = id, Time = time, Survival = St)
  formatted_nrisk_data <- nrisk_data %>% dplyr::rename(Interval = id, Time = time, Lower = start_id, Upper = end_id, nrisk = nrisk)

  temp_surv_file <- tempfile(fileext = ".txt")
  temp_nrisk_file <- tempfile(fileext = ".txt")
  utils::write.table(formatted_surv_data, temp_surv_file, sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(formatted_nrisk_data, temp_nrisk_file, sep = "\t", row.names = FALSE, quote = FALSE)

  cat(paste("     -> Calling 'digitise' to generate output files.\n"))
  km_output_path <- file.path(output_directory, paste0("KMdata_", arm_name, "_final.txt"))
  ipd_output_path <- file.path(output_directory, paste0("IPDdata_", arm_name, "_final.txt"))

  digitise_function(
    surv_inp = temp_surv_file,
    nrisk_inp = temp_nrisk_file,
    km_output = km_output_path,
    ipd_output = ipd_output_path
  )

  unlink(c(temp_surv_file, temp_nrisk_file))

  # --- STEP E: Read and return the generated IPD data ---
  final_ipd <- utils::read.delim(ipd_output_path)
  final_ipd$arm <- arm_name

  cat(paste("--- Processing for", arm_name, "completed. ---\n"))
  return(final_ipd)
}


# 3. Main Workflow Function (Final Version)
# ------------------------------------------------

#' Digitize Kaplan-Meier Curves and Reconstruct IPD
#'
#' This function orchestrates the entire workflow for reconstructing Individual
#' Patient Data (IPD) from an image of Kaplan-Meier (KM) curves. It first
#' digitizes the curves from the image using `SurvdigitizeR::survival_digitize`,
#' then processes each curve individually to clean the data and generate IPD
#' using the Guyot et al. algorithm via `survHE::digitise`.
#'
#' @param img_path Character string. The full path to the image file (e.g., PNG, JPG)
#'   containing the KM curves.
#' @param time_breaks Numeric vector. A sequence of time points that define the
#'   intervals for the number-at-risk table. For example, `seq(0, 60, 12)`.
#' @param n_risk_list A named list where each element is a numeric vector
#'   representing the number of patients at risk for one arm. The names of the
#'   list elements will be used as arm names (e.g., `list(control = c(...), treatment = c(...))`).
#'   The order of arms must correspond to the order of curves digitized (curve 1, curve 2).
#' @param output_filename Character string or `NULL`. The name for the output
#'   `.Rda` file that will contain the final combined IPD data frame. If `NULL` or
#'   empty, the data is returned but not saved to a file. The file is saved in the same directory as `img_path`.
#' @param num_curves Integer. The number of KM curves to digitize from the image.
#'   Default is 2.
#' @param x_start,x_end,x_increment Numeric values defining the x-axis (time) scale.
#' @param y_start,y_end,y_increment Numeric values defining the y-axis (survival probability) scale.
#' @param y_text_vertical Logical. Set to `TRUE` if the y-axis tick labels are
#'   written vertically. Default is `TRUE`.
#' @param censoring Logical. Set to `TRUE` if censoring marks are present on the
#'   curves and need to be digitized. Default is `FALSE`.
#'
#' @importFrom SurvdigitizeR survival_digitize
#' @importFrom survHE digitise
#' @importFrom dplyr filter
#' @importFrom tools file_path_sans_ext
#'
#' @return A single data frame containing the combined IPD for all processed arms.
#'   The data frame includes columns for 'time', 'event', and 'arm'.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # This is a conceptual example as it requires an actual image file.
#' # 1. Define parameters
#' image_file <- "path/to/your/km_curve_image.png"
#' output_data_file <- "reconstructed_ipd.Rda"
#'
#' # Time intervals from the x-axis of the publication's graph
#' time_points <- seq(0, 60, by = 12)
#'
#' # Number at risk data from the publication's table
#' risk_data <- list(
#'   Control = c(500, 450, 400, 300, 200, 100),
#'   Treatment = c(502, 480, 460, 410, 350, 250)
#' )
#'
#' # 2. Run the main function
#' final_ipd_dataset <- km_to_dataset(
#'   img_path = image_file,
#'   time_breaks = time_points,
#'   n_risk_list = risk_data,
#'   output_filename = output_data_file,
#'   num_curves = 2,
#'   x_start = 0, x_end = 60, x_increment = 12,
#'   y_start = 0, y_end = 1, y_increment = 0.2
#' )
#'
#' # 3. Analyze the resulting dataset
#' if (exists("final_ipd_dataset")) {
#'   library(survival)
#'   library(survminer)
#'
#'   fit <- survfit(Surv(time, event) ~ arm, data = final_ipd_dataset)
#'   print(summary(fit))
#'
#'   ggsurvplot(
#'      fit,
#'      data = final_ipd_dataset,
#'      risk.table = TRUE,
#'      legend.title = "Group"
#'   )
#' }
#'}
km_to_dataset <- function(
    img_path,
    time_breaks,
    n_risk_list,
    output_filename, # <-- ADDED: Argument for the output filename
    # Arguments for survival_digitize
    num_curves = 2,
    x_start, x_end, x_increment,
    y_start, y_end, y_increment,
    y_text_vertical = TRUE,
    censoring = FALSE
) {

  # Automatically get the output directory from the image path
  output_dir <- dirname(img_path)
  cat(paste("--- Output directory set to:", output_dir, "---\n"))

  # --- STEP 3.1: Digitize the image ONLY ONCE ---
  cat("--- Starting image digitization ---\n")
  digitized_data <- SurvdigitizeR::survival_digitize(
    img_path = img_path, num_curves = num_curves, x_start = x_start, x_end = x_end,
    x_increment = x_increment, y_start = y_start, y_end = y_end, y_increment = y_increment,
    y_text_vertical = y_text_vertical, censoring = censoring
  )
  cat("--- Image digitization complete ---\n")

  # --- STEP 3.2: Process each curve ---
  control_data <- digitized_data %>% dplyr::filter(curve == 1)
  exp_data     <- digitized_data %>% dplyr::filter(curve == 2)

  arm1_name <- names(n_risk_list)[1]
  arm2_name <- names(n_risk_list)[2]
  nrisk_arm1 <- n_risk_list[[1]]
  nrisk_arm2 <- n_risk_list[[2]]

  ipd_arm1 <- process_curve_and_generate_ipd(
    curve_df = control_data, n_risk = nrisk_arm1, time_breaks = time_breaks,
    arm_name = arm1_name, output_directory = output_dir, digitise_function = survHE::digitise
  )

  ipd_arm2 <- process_curve_and_generate_ipd(
    curve_df = exp_data, n_risk = nrisk_arm2, time_breaks = time_breaks,
    arm_name = arm2_name, output_directory = output_dir, digitise_function = survHE::digitise
  )

  # Combine the IPD data from both arms
  total_data <- rbind(ipd_arm1, ipd_arm2)

  # --- MODIFIED: Save file using the provided argument ---
  if (!is.null(output_filename) && nchar(output_filename) > 0) {
    full_output_path <- file.path(output_dir, output_filename)
    save(total_data, file = full_output_path)
    cat(paste("\n--- Full workflow completed. Final dataset saved to:", full_output_path, "---\n"))
  } else {
    cat("\n--- Full workflow completed. No output filename provided, dataset not saved. ---\n")
  }

  return(total_data)
}
