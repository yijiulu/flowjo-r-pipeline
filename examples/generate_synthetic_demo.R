#!/usr/bin/env Rscript
# Generate synthetic FlowJo export and plate metadata tables, then execute both CLI pipelines.

set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
outdir <- "outputs/demo"
if (length(args) >= 2 && args[1] == "--outdir") {
  outdir <- args[2]
}

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
repo_root <- normalizePath(file.path(script_dir, ".."))

demo_data_dir <- file.path(outdir, "synthetic_data")
dir.create(demo_data_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Synthetic time-course assay (3 constructs x 3 conditions x 4 timepoints x 3 biological replicates)
strains <- c("Construct_A", "Construct_B", "Pos_Control")
conditions <- c("-/-", "+/-", "+/+")
times <- c(24, 48, 72, 96)
reps <- 1:3

grid_tc <- expand.grid(
  Strain = strains,
  Condition = conditions,
  Time_h = times,
  Replicate = reps,
  stringsAsFactors = FALSE
)
grid_tc$Species <- "Microbe_A"
grid_tc$Run <- paste0("Run_T", grid_tc$Time_h)

well_pool <- as.vector(outer(LETTERS[1:8], sprintf("%02d", 1:12), paste0))
grid_tc$Well <- ave(seq_len(nrow(grid_tc)), grid_tc$Run, FUN = function(idx) well_pool[seq_along(idx)])
grid_tc$Sample_ID <- paste(grid_tc$Run, grid_tc$Well, sep = "_")

base_signal <- ifelse(grid_tc$Strain == "Pos_Control", 180, ifelse(grid_tc$Strain == "Construct_B", 110, 55))
cond_boost  <- ifelse(grid_tc$Condition == "+/+", 2.8, ifelse(grid_tc$Condition == "+/-", 1.7, 1.0))
time_factor <- grid_tc$Time_h / 24
gfp_vals    <- pmax(5, rnorm(nrow(grid_tc), mean = base_signal * cond_boost * (0.6 + 0.35 * time_factor), sd = 18))

tc_flowjo <- data.frame(
  X = paste0("C:\\FlowJo\\Demo\\", grid_tc$Run, ".fcs"),
  X.WELLID = grid_tc$Well,
  Microbe.Forward.Scatter..FSC.HLin...Count = sample(8000:15000, nrow(grid_tc), replace = TRUE),
  Microbe.Forward.Scatter..FSC.HLin...Median..Green.B.Fluorescence..GRN.B.HLin.. = round(gfp_vals, 2),
  stringsAsFactors = FALSE
)
# Append FlowJo-style summary footer rows to verify parser filtering
tc_flowjo <- rbind(
  tc_flowjo,
  data.frame(
    X = c("Mean", "SD"),
    X.WELLID = c("", ""),
    Microbe.Forward.Scatter..FSC.HLin...Count = c(11500, 1200),
    Microbe.Forward.Scatter..FSC.HLin...Median..Green.B.Fluorescence..GRN.B.HLin.. = c(200, 50),
    stringsAsFactors = FALSE
  )
)

tc_flowjo_path <- file.path(demo_data_dir, "demo_timecourse_flowjo.csv")
tc_meta_path   <- file.path(demo_data_dir, "demo_timecourse_metadata.csv")
write.csv(tc_flowjo, tc_flowjo_path, row.names = FALSE)
write.csv(grid_tc, tc_meta_path, row.names = FALSE)

# 2. Synthetic dose-response assay (4 mediators x 4 concentrations x 2 timepoints x 3 replicates)
grid_dr <- expand.grid(
  Condition = c("Mediator_1", "Mediator_2", "Mediator_3", "Mediator_4"),
  Concentration = c(0, 5, 25, 100),
  Time_h = c(48, 72),
  Replicate = 1:3,
  stringsAsFactors = FALSE
)
grid_dr$Species <- "Microbe_A"
grid_dr$Strain  <- "Construct_B"
grid_dr$Run     <- paste0("Dose_T", grid_dr$Time_h)
grid_dr$Well    <- ave(seq_len(nrow(grid_dr)), grid_dr$Run, FUN = function(idx) well_pool[seq_along(idx)])
grid_dr$Sample_ID <- paste(grid_dr$Run, grid_dr$Well, sep = "_")

potency <- c(Mediator_1 = 3.2, Mediator_2 = 1.9, Mediator_3 = 2.5, Mediator_4 = 1.1)[grid_dr$Condition]
dr_mean <- 60 + potency * 45 * (grid_dr$Concentration / (15 + grid_dr$Concentration)) * (grid_dr$Time_h / 48)
dr_gfp  <- pmax(5, rnorm(nrow(grid_dr), mean = dr_mean, sd = 12))

dr_flowjo <- data.frame(
  X = paste0("/data/flowjo/demo/", grid_dr$Run, ".fcs"),
  X.WELLID = grid_dr$Well,
  Microbe.Forward.Scatter..FSC.HLin...Count = sample(9000:14000, nrow(grid_dr), replace = TRUE),
  Microbe.Forward.Scatter..FSC.HLin...Median..Green.B.Fluorescence..GRN.B.HLin.. = round(dr_gfp, 2),
  stringsAsFactors = FALSE
)

dr_flowjo_path <- file.path(demo_data_dir, "demo_dose_flowjo.csv")
dr_meta_path   <- file.path(demo_data_dir, "demo_dose_metadata.csv")
write.csv(dr_flowjo, dr_flowjo_path, row.names = FALSE)
write.csv(grid_dr, dr_meta_path, row.names = FALSE)

# Run both CLI scripts on the synthetic data
status1 <- system2(
  "Rscript",
  c(
    file.path(repo_root, "scripts", "analyze_timecourse.R"),
    "--flowjo", tc_flowjo_path,
    "--metadata", tc_meta_path,
    "--species", "Microbe_A",
    "--outdir", file.path(outdir, "timecourse")
  )
)
stopifnot(status1 == 0)

status2 <- system2(
  "Rscript",
  c(
    file.path(repo_root, "scripts", "analyze_dose_response.R"),
    "--flowjo", dr_flowjo_path,
    "--metadata", dr_meta_path,
    "--species", "Microbe_A",
    "--outdir", file.path(outdir, "dose_response")
  )
)
stopifnot(status2 == 0)

cat(sprintf("[DONE] Synthetic demo completed cleanly in %s\n", outdir))
