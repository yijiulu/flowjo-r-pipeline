# Publication-grade layered time-course and dose-response ggplot2 builders.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

TIMEPOINT_PALETTE <- data.frame(
  Time_h = c(96, 72, 48, 24),
  width  = c(0.70, 0.55, 0.40, 0.25),
  fill   = c("#3F1151", "#405187", "#478F8B", "#7AC76E"),
  border = c("#26082F", "#2A3559", "#2F5F5C", "#4A8F41"),
  alpha  = c(0.45, 0.45, 0.45, 0.90),
  label  = c("96h", "72h", "48h", "24h"),
  stringsAsFactors = FALSE
)

resolve_timepoint_style <- function(times) {
  times_sorted <- sort(unique(as.numeric(times)), decreasing = TRUE)
  idx <- match(times_sorted, TIMEPOINT_PALETTE$Time_h)

  if (any(is.na(idx))) {
    n <- length(times_sorted)
    return(data.frame(
      Time_h = times_sorted,
      width  = seq(0.70, 0.25, length.out = n),
      fill   = scales::viridis_pal(option = "viridis", begin = 0.15, end = 0.85)(n),
      border = rep("#1E293B", n),
      alpha  = c(rep(0.45, max(0, n - 1)), 0.90),
      label  = paste0(times_sorted, "h"),
      stringsAsFactors = FALSE
    ))
  }

  TIMEPOINT_PALETTE[idx, , drop = FALSE]
}

add_layered_bars <- function(p, df, style_df, y_col, se_col) {
  for (i in seq_len(nrow(style_df))) {
    t_val <- style_df$Time_h[i]
    dt <- df[df$Time_h == t_val, , drop = FALSE]
    if (nrow(dt) == 0) next

    p <- p +
      geom_col(
        data      = dt,
        aes(y = .data[[y_col]]),
        fill      = style_df$fill[i],
        width     = style_df$width[i],
        alpha     = style_df$alpha[i],
        colour    = style_df$border[i],
        linewidth = 0.4
      ) +
      geom_errorbar(
        data = dt,
        aes(
          ymin = .data[[y_col]] - .data[[se_col]],
          ymax = .data[[y_col]] + .data[[se_col]]
        ),
        width     = style_df$width[i] * 0.25,
        linewidth = 0.4,
        na.rm     = TRUE
      )
  }

  legend_df <- data.frame(
    x    = rep(NA_real_, nrow(style_df)),
    y    = rep(NA_real_, nrow(style_df)),
    Time = factor(style_df$label, levels = rev(style_df$label))
  )

  p +
    geom_point(
      data   = legend_df,
      aes(x = .data$x, y = .data$y, fill = .data$Time, colour = .data$Time),
      shape  = 22,
      size   = 4,
      stroke = 0.5,
      na.rm  = TRUE
    ) +
    scale_fill_manual(
      name   = "Time",
      values = stats::setNames(style_df$fill, style_df$label),
      breaks = rev(style_df$label)
    ) +
    scale_colour_manual(
      name   = "Time",
      values = stats::setNames(style_df$border, style_df$label),
      breaks = rev(style_df$label)
    )
}

#' Build a layered multi-timepoint barplot (raw or baseline-corrected GFP).
plot_layered_timecourse <- function(
  summary_df,
  x_col = "Condition",
  y_col = "Mean_GFP",
  se_col = "SE_GFP",
  times = NULL,
  x_levels = NULL,
  facet_col = "Strain",
  facet_labels = NULL,
  title = NULL,
  y_label = "Mean GFP (a.u.)"
) {
  d <- summary_df
  if (is.null(times)) {
    times <- sort(unique(d$Time_h), decreasing = TRUE)
  }
  d <- d[d$Time_h %in% times, , drop = FALSE]

  if (!is.null(x_levels)) {
    d[[x_col]] <- factor(d[[x_col]], levels = x_levels)
  } else {
    d[[x_col]] <- factor(d[[x_col]])
  }

  style_df <- resolve_timepoint_style(times)

  p <- ggplot(mapping = aes(x = .data[[x_col]]))
  p <- add_layered_bars(p, d, style_df, y_col = y_col, se_col = se_col)

  if (!is.null(facet_col) && facet_col %in% names(d)) {
    p <- p + facet_wrap(
      stats::as.formula(paste("~", facet_col)),
      nrow = 1,
      strip.position = "bottom",
      labeller = if (!is.null(facet_labels)) {
        do.call(labeller, stats::setNames(list(facet_labels), facet_col))
      } else {
        "label_value"
      }
    )
  }

  p +
    labs(title = title, x = NULL, y = y_label) +
    theme_classic(base_size = 11) +
    theme(
      strip.placement  = "outside",
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold"),
      axis.text.x      = element_text(angle = 45, hjust = 1, colour = "black"),
      axis.text.y      = element_text(colour = "black"),
      panel.spacing    = unit(0.5, "cm")
    )
}

#' Build a faceted dose-response barplot across mediators and concentrations.
plot_dose_response <- function(
  summary_df,
  times = c(72, 48),
  condition_levels = NULL,
  y_col = "Mean_GFP",
  se_col = "SE_GFP",
  title = "Mediator dose-response profile",
  y_label = "Mean GFP (a.u.)"
) {
  d <- summary_df[summary_df$Time_h %in% times, , drop = FALSE]
  if (!is.null(condition_levels)) {
    d$Condition <- factor(d$Condition, levels = condition_levels)
  } else {
    d$Condition <- factor(d$Condition)
  }
  d$Concentration_f <- factor(d$Concentration)

  style_df <- resolve_timepoint_style(times)

  p <- ggplot(mapping = aes(x = .data$Concentration_f))
  p <- add_layered_bars(p, d, style_df, y_col = y_col, se_col = se_col)

  p +
    facet_grid(
      cols   = vars(.data$Condition),
      scales = "free_x",
      space  = "free_x",
      switch = "x"
    ) +
    labs(title = title, x = NULL, y = y_label) +
    theme_classic(base_size = 11) +
    theme(
      strip.placement  = "outside",
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold"),
      axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1, size = 8, colour = "black"),
      axis.text.y      = element_text(colour = "black"),
      panel.spacing    = unit(0.15, "cm")
    )
}
