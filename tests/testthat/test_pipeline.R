library(testthat)

repo_root <- normalizePath(file.path("..", ".."))
source(file.path(repo_root, "R", "io.R"))
source(file.path(repo_root, "R", "stats.R"))
source(file.path(repo_root, "R", "plotting.R"))

test_that("read_flowjo_table handles Windows/POSIX paths and filters footer rows", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(tmp), add = TRUE)

  df_in <- data.frame(
    X = c("C:\\Runs\\plate1.fcs", "/mnt/data/plate2.fcs", "Mean", "SD"),
    X.WELLID = c("A01", "H12", "", "Z99"),
    Gate.Count = c(10000, 12000, 11000, 1000),
    Gate.Median..Green.B.Fluorescence..GRN.B.HLin.. = c(123.4, 456.7, 290.0, 200.0),
    stringsAsFactors = FALSE
  )
  write.csv(df_in, tmp, row.names = FALSE)

  out <- read_flowjo_table(tmp)
  expect_equal(nrow(out), 2L)
  expect_equal(out$Run, c("plate1", "plate2"))
  expect_equal(out$Sample_ID, c("plate1_A01", "plate2_H12"))
  expect_equal(out$Median_GFP, c(123.4, 456.7))
})

test_that("summarise_replicates uses unbiased sample SEM (n - 1 denominator)", {
  df <- data.frame(
    Species = "S1",
    Strain = "C1",
    Time_h = 48,
    Condition = c("-/-", "-/-", "-/-", "+/+", "+/+", "+/+"),
    Median_GFP = c(10, 13, 16, 40, 46, 52),
    stringsAsFactors = FALSE
  )

  res <- summarise_replicates(df)
  base_row <- res[res$Condition == "-/-", ]
  pos_row  <- res[res$Condition == "+/+", ]

  expect_equal(base_row$Mean_GFP, 13)
  expect_equal(base_row$SE_GFP, sd(c(10, 13, 16)) / sqrt(3), tolerance = 1e-12)
  expect_equal(pos_row$SE_GFP, sd(c(40, 46, 52)) / sqrt(3), tolerance = 1e-12)

  net <- subtract_baseline(res, baseline_condition = "-/-")
  pos_net <- net[net$Condition == "+/+", ]
  expect_equal(pos_net$Corrected_GFP, 46 - 13)
  expect_equal(
    pos_net$Corrected_SE,
    sqrt(pos_row$SE_GFP^2 + base_row$SE_GFP^2),
    tolerance = 1e-12
  )
})

test_that("plotting builders construct valid ggplot objects with matched error bars and legends", {
  df <- expand.grid(
    Species = "S1",
    Strain = c("C1", "C2"),
    Time_h = c(48, 96),
    Condition = c("+/-", "+/+"),
    stringsAsFactors = FALSE
  )
  df$Mean_GFP <- seq(100, 240, length.out = nrow(df))
  df$SE_GFP <- 5
  df$Corrected_GFP <- df$Mean_GFP - 40
  df$Corrected_SE <- 7.07

  p_tc <- plot_layered_timecourse(
    df,
    y_col = "Corrected_GFP",
    se_col = "Corrected_SE",
    times = c(96, 48),
    x_levels = c("+/-", "+/+")
  )
  expect_s3_class(p_tc, "ggplot")
  built <- ggplot_build(p_tc)
  # Ensure error bar layer reflects Corrected_SE (7.07), not SE_GFP (5.0)
  err_layer <- built$data[[2]]
  expect_equal(err_layer$ymax[1] - err_layer$ymin[1], 2 * 7.07, tolerance = 1e-6)

  df_dr <- expand.grid(
    Species = "S1",
    Strain = "C1",
    Time_h = c(48, 72),
    Condition = c("M1", "M2"),
    Concentration = c(0, 10, 50),
    stringsAsFactors = FALSE
  )
  df_dr$Mean_GFP <- seq(50, 200, length.out = nrow(df_dr))
  df_dr$SE_GFP <- 4
  p_dr <- plot_dose_response(df_dr, times = c(72, 48))
  expect_s3_class(p_dr, "ggplot")
})
