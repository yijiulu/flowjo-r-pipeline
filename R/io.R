# FlowJo CSV ingestion and metadata joining utilities.

suppressPackageStartupMessages({
  library(dplyr)
})

resolve_column <- function(df, explicit_name, regex_pattern, label) {
  if (!is.null(explicit_name) && nzchar(explicit_name)) {
    if (!explicit_name %in% names(df)) {
      stop(sprintf("Column '%s' for %s not found in table.", explicit_name, label), call. = FALSE)
    }
    return(explicit_name)
  }
  matches <- grep(regex_pattern, names(df), value = TRUE, ignore.case = TRUE)
  if (length(matches) == 0) {
    stop(
      sprintf("Could not auto-detect %s column matching pattern '%s'.", label, regex_pattern),
      call. = FALSE
    )
  }
  matches[1]
}

#' Read and standardise a FlowJo exported summary CSV table.
#'
#' @param path Path to the FlowJo CSV export.
#' @param file_col Optional column name holding the FCS file path.
#' @param well_col Optional column name holding 96-well coordinates (e.g. A01).
#' @param count_col Optional column name holding gated cell counts.
#' @param fluor_col Optional column name holding median fluorescence intensity.
#' @return A clean data.frame with columns Sample_ID, Run, Well, Cells_count, Median_GFP.
read_flowjo_table <- function(
  path,
  file_col = NULL,
  well_col = NULL,
  count_col = NULL,
  fluor_col = NULL
) {
  if (!file.exists(path)) {
    stop(sprintf("FlowJo table not found: %s", path), call. = FALSE)
  }
  raw <- read.csv(path, header = TRUE, stringsAsFactors = FALSE, check.names = TRUE)

  col_file  <- resolve_column(raw, file_col,  "^(X|File|Sample)$", "file path")
  col_well  <- resolve_column(raw, well_col,  "(WELLID|^Well$)", "well coordinate")
  col_count <- resolve_column(raw, count_col, "(Count|Cells_count)$", "cell count")
  col_fluor <- resolve_column(raw, fluor_col, "(Median.*(GRN|GFP|Fluorescence)|^Median_GFP$)", "median fluorescence")

  cleaned <- data.frame(
    File        = as.character(raw[[col_file]]),
    Well        = trimws(as.character(raw[[col_well]])),
    Cells_count = suppressWarnings(as.numeric(raw[[col_count]])),
    Median_GFP  = suppressWarnings(as.numeric(raw[[col_fluor]])),
    stringsAsFactors = FALSE
  )

  # Keep valid 96-well plate coordinates (A01-H12 or A1-H12) and drop FlowJo Mean/SD footer rows.
  valid_well <- grepl("^[A-H](0?[1-9]|1[0-2])$", cleaned$Well)
  cleaned <- cleaned[valid_well, , drop = FALSE]

  # Support both Windows (\) and POSIX (/) FCS export paths.
  cleaned$Run <- sub("\\.fcs$", "", sub("^.*[/\\\\]", "", cleaned$File), ignore.case = TRUE)
  cleaned$Sample_ID <- paste(cleaned$Run, cleaned$Well, sep = "_")

  if (anyDuplicated(cleaned$Sample_ID) > 0) {
    dups <- unique(cleaned$Sample_ID[duplicated(cleaned$Sample_ID)])
    stop(
      sprintf("Duplicate Sample_ID entries detected: %s", paste(head(dups, 5), collapse = ", ")),
      call. = FALSE
    )
  }

  cleaned[, c("Sample_ID", "Run", "Well", "Cells_count", "Median_GFP"), drop = FALSE]
}

#' Join cleaned FlowJo measurements with plate-layout metadata.
#'
#' @param flow_df Data frame returned by read_flowjo_table().
#' @param metadata_df Data frame containing Run, Well, Species, Strain, Time_h, Condition, etc.
#' @param species Optional species filter applied to metadata_df before joining.
#' @param strict_match If TRUE, error when any FlowJo wells fail to match metadata.
#' @return Joined data.frame.
join_flowjo_metadata <- function(flow_df, metadata_df, species = NULL, strict_match = TRUE) {
  meta <- metadata_df
  if (!is.null(species) && nzchar(species)) {
    if (!"Species" %in% names(meta)) {
      stop("Metadata table must contain a 'Species' column when filtering by species.", call. = FALSE)
    }
    meta <- meta[meta$Species == species, , drop = FALSE]
  }

  if ("Sample_ID" %in% names(meta)) {
    meta <- meta[, setdiff(names(meta), "Sample_ID"), drop = FALSE]
  }

  joined <- flow_df |>
    left_join(meta, by = c("Run", "Well"), relationship = "one-to-one")

  if (strict_match && "Strain" %in% names(joined) && any(is.na(joined$Strain))) {
    n_missing <- sum(is.na(joined$Strain))
    stop(
      sprintf("%d FlowJo well(s) failed to match metadata on (Run, Well).", n_missing),
      call. = FALSE
    )
  }

  joined
}

#' Read an optional two-column CSV (raw, label) for display labels on facets.
read_label_map <- function(path) {
  if (is.null(path) || !nzchar(path)) {
    return(NULL)
  }
  if (!file.exists(path)) {
    stop(sprintf("Label map file not found: %s", path), call. = FALSE)
  }
  df <- read.csv(path, header = TRUE, stringsAsFactors = FALSE)
  if (!all(c("raw", "label") %in% names(df))) {
    stop("Label map CSV must contain columns 'raw' and 'label'.", call. = FALSE)
  }
  stats::setNames(as.character(df$label), as.character(df$raw))
}
