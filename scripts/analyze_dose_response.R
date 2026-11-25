#!/usr/bin/env Rscript
# CLI runner for mediator dose-response FlowJo cytometry assays.

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  getwd()
}

parse_cli_args <- function(args) {
  opts <- list(
    flowjo    = NULL,
    metadata  = NULL,
    outdir    = "outputs/dose_response",
    species   = NULL,
    times     = "72,48",
    count_col = NULL,
    fluor_col = NULL
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    if (key %in% c("-h", "--help")) {
      cat(
        "Usage: Rscript scripts/analyze_dose_response.R --flowjo <csv> --metadata <csv>",
        "[--outdir <dir>] [--species <name>] [--times '72,48']\n"
      )
      quit(status = 0)
    }
    if (i == length(args)) stop(sprintf("Missing value for argument %s", key), call. = FALSE)
    val <- args[i + 1]
    clean_key <- gsub("-", "_", sub("^--", "", key))
    if (!clean_key %in% names(opts)) stop(sprintf("Unknown argument: %s", key), call. = FALSE)
    opts[[clean_key]] <- val
    i <- i + 2
  }
  if (is.null(opts$flowjo) || is.null(opts$metadata)) {
    stop("Both --flowjo and --metadata are required. Pass --help for usage.", call. = FALSE)
  }
  opts
}

repo_root <- normalizePath(file.path(get_script_dir(), ".."))
source(file.path(repo_root, "R", "io.R"))
source(file.path(repo_root, "R", "stats.R"))
source(file.path(repo_root, "R", "plotting.R"))

opts <- parse_cli_args(commandArgs(trailingOnly = TRUE))
dir.create(opts$outdir, recursive = TRUE, showWarnings = FALSE)

times <- as.numeric(strsplit(opts$times, ",", fixed = TRUE)[[1]])

flow_df <- read_flowjo_table(
  opts$flowjo,
  count_col = opts$count_col,
  fluor_col = opts$fluor_col
)
meta_df <- read.csv(opts$metadata, header = TRUE, stringsAsFactors = FALSE)
joined  <- join_flowjo_metadata(flow_df, meta_df, species = opts$species)

write.csv(flow_df, file.path(opts$outdir, "flowjo_clean.csv"), row.names = FALSE)

summary_df <- summarise_replicates(
  joined,
  group_cols = c("Species", "Strain", "Time_h", "Condition", "Concentration")
)
write.csv(summary_df, file.path(opts$outdir, "dose_response_mean_se.csv"), row.names = FALSE)

p <- plot_dose_response(
  summary_df,
  times = times,
  title = "Redox mediator dose-response profile"
)
ggsave(file.path(opts$outdir, "dose_response.png"), p, width = 10, height = 6.5, dpi = 300)

cat(sprintf("[OK] Dose-response outputs written to %s\n", opts$outdir))
