#!/usr/bin/env Rscript
# CLI runner for multi-timepoint FlowJo cytometry analysis and baseline-corrected plots.

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
    outdir    = "outputs/timecourse",
    species   = NULL,
    baseline  = "-/-",
    times     = "96,72,48,24",
    label_map = NULL,
    count_col = NULL,
    fluor_col = NULL
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    if (key %in% c("-h", "--help")) {
      cat(
        "Usage: Rscript scripts/analyze_timecourse.R --flowjo <csv> --metadata <csv>",
        "[--outdir <dir>] [--species <name>] [--baseline '-/-'] [--times '96,72,48,24']",
        "[--label-map <csv>] [--count-col <col>] [--fluor-col <col>]\n"
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
label_map <- read_label_map(opts$label_map)

flow_df <- read_flowjo_table(
  opts$flowjo,
  count_col = opts$count_col,
  fluor_col = opts$fluor_col
)
meta_df <- read.csv(opts$metadata, header = TRUE, stringsAsFactors = FALSE)
joined  <- join_flowjo_metadata(flow_df, meta_df, species = opts$species)

write.csv(flow_df, file.path(opts$outdir, "flowjo_clean.csv"), row.names = FALSE)

summary_df <- summarise_replicates(joined)

has_baseline <- !is.null(opts$baseline) && opts$baseline %in% summary_df$Condition
if (has_baseline) {
  summary_df <- subtract_baseline(summary_df, baseline_condition = opts$baseline)
}

write.csv(summary_df, file.path(opts$outdir, "summary_mean_se.csv"), row.names = FALSE)

cond_levels <- unique(summary_df$Condition)
if (all(c("-/-", "+/-", "+/+") %in% cond_levels)) {
  cond_levels <- c("-/-", "+/-", "+/+")
}

title_prefix <- if (!is.null(opts$species) && nzchar(opts$species)) opts$species else "All species"

p_raw <- plot_layered_timecourse(
  summary_df,
  x_col        = "Condition",
  y_col        = "Mean_GFP",
  se_col       = "SE_GFP",
  times        = times,
  x_levels     = cond_levels,
  facet_col    = "Strain",
  facet_labels = label_map,
  title        = sprintf("%s: Mean GFP across timepoints", title_prefix),
  y_label      = "Mean GFP (a.u.)"
)
ggsave(file.path(opts$outdir, "timecourse_raw_gfp.png"), p_raw, width = 9.2, height = 6.5, dpi = 300)

if (has_baseline) {
  net_df <- summary_df[summary_df$Condition != opts$baseline, , drop = FALSE]
  net_levels <- setdiff(cond_levels, opts$baseline)
  p_net <- plot_layered_timecourse(
    net_df,
    x_col        = "Condition",
    y_col        = "Corrected_GFP",
    se_col       = "Corrected_SE",
    times        = times,
    x_levels     = net_levels,
    facet_col    = "Strain",
    facet_labels = label_map,
    title        = sprintf("%s: Baseline-corrected net GFP", title_prefix),
    y_label      = "Net mean GFP (a.u.)"
  )
  ggsave(file.path(opts$outdir, "timecourse_net_gfp.png"), p_net, width = 9.2, height = 6.5, dpi = 300)
}

cat(sprintf("[OK] Time-course outputs written to %s\n", opts$outdir))
