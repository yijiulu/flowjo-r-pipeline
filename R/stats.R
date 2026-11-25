# Replicate summary statistics and error-propagated baseline subtraction.

suppressPackageStartupMessages({
  library(dplyr)
})

#' Compute biological replicate mean, sample SD, and unbiased SEM (sd / sqrt(n)).
#'
#' @param df Joined cytometry data frame.
#' @param group_cols Character vector of grouping columns.
#' @param value_col Numeric column to summarise (default: "Median_GFP").
#' @return Summarised tibble with N, Mean_GFP, SD_GFP, SE_GFP.
summarise_replicates <- function(
  df,
  group_cols = c("Species", "Strain", "Time_h", "Condition"),
  value_col = "Median_GFP"
) {
  missing_cols <- setdiff(c(group_cols, value_col), names(df))
  if (length(missing_cols) > 0) {
    stop(
      sprintf("Missing required column(s): %s", paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  df |>
    group_by(across(all_of(group_cols))) |>
    summarise(
      N        = sum(!is.na(.data[[value_col]])),
      Mean_GFP = mean(.data[[value_col]], na.rm = TRUE),
      SD_GFP   = if (sum(!is.na(.data[[value_col]])) >= 2) stats::sd(.data[[value_col]], na.rm = TRUE) else NA_real_,
      SE_GFP   = if (sum(!is.na(.data[[value_col]])) >= 2) {
        stats::sd(.data[[value_col]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[value_col]])))
      } else {
        NA_real_
      },
      .groups = "drop"
    ) |>
    arrange(across(all_of(group_cols)))
}

#' Subtract a matched baseline condition and propagate standard errors in quadrature.
#'
#' @param summary_df Output of summarise_replicates().
#' @param baseline_condition Value of Condition used as the negative/baseline control.
#' @param key_cols Matching key columns linking each condition to its baseline.
#' @return Data frame augmented with Baseline_GFP, Baseline_SE, Corrected_GFP, Corrected_SE.
subtract_baseline <- function(
  summary_df,
  baseline_condition = "-/-",
  key_cols = c("Species", "Strain", "Time_h")
) {
  baseline <- summary_df |>
    filter(.data$Condition == baseline_condition) |>
    select(
      all_of(key_cols),
      Baseline_GFP = "Mean_GFP",
      Baseline_SE  = "SE_GFP"
    )

  if (nrow(baseline) == 0) {
    stop(
      sprintf("No rows found for baseline condition '%s'.", baseline_condition),
      call. = FALSE
    )
  }

  summary_df |>
    left_join(baseline, by = key_cols, relationship = "many-to-one") |>
    mutate(
      Corrected_GFP = .data$Mean_GFP - .data$Baseline_GFP,
      Corrected_SE  = sqrt(.data$SE_GFP^2 + .data$Baseline_SE^2)
    ) |>
    arrange(across(all_of(c(key_cols, "Condition"))))
}
