#' Fit and plot a mixed-effects cosinor model for a single group
#'
#' Fits a linear mixed-effects cosinor model (24 h period) to replicate-averaged
#' intensity data for one experimental group, tests rhythmicity against a null
#' (intercept-only) model via likelihood ratio test, and produces a scatter +
#' fitted-curve plot. The fitted curve and significance markers (dashed fit line,
#' dotted 4h gridlines) are only drawn if the cosinor fit is significant (p <= 0.05).
#'
#' @param df Data frame containing raw single-cell/replicate measurements. Must
#'   include columns: `Group`, `Timepoint` (format "T<hours>", e.g. "T12"),
#'   `Replicate`, and `Mean` (the raw intensity measurement to be averaged).
#' @param group_name Character. Value of `Group` to filter and plot (e.g. "WT").
#' @param color_hex Character. Hex color code (e.g. "#1F77B4") used for the
#'   group's summary (mean) points.
#' @param outfile Character. File path to save the plot to (passed to `ggsave`).
#' @param protein_name Character. Protein/marker name, used in the plot title.
#' @param scale Logical. If `TRUE`, min-max scales `Mean_val` (and the fitted
#'   curve) to the [0, 1] range before plotting. Default `FALSE`.
#' @param y_max Numeric or `NULL`. If provided, sets the upper y-axis limit via
#'   `ylim(0, y_max)`. Overrides the default `ylim(0, 1)` set when `scale = TRUE`.
#'   Default `NULL`.
#'
#' @details
#' Processing steps:
#' \enumerate{
#'   \item Filter `df` to `Group == group_name` and parse `Timepoint` into a
#'     numeric `Time_h` column.
#'   \item Average `Mean` within each `Replicate` x `Time_h` combination
#'     (`df_rep_avg`).
#'   \item Fit `Mean_val ~ cos(2*pi*Time_h/24) + sin(2*pi*Time_h/24) + (1 | Replicate)`
#'     via `lme4::lmer`, compared against a random-intercept-only null model
#'     using a likelihood ratio test (`anova`) to obtain `p_val`.
#'   \item Derive amplitude (`sqrt(beta_cos^2 + beta_sin^2)`) and phase (in
#'     hours, wrapped to `[0, 24)`) from the fixed-effect cosine/sine
#'     coefficients.
#'   \item Optionally min-max scale the replicate averages and the fitted
#'     curve to `[0, 1]`.
#'   \item Summarize replicate averages per `Time_h` (mean) for the colored
#'     summary points.
#'   \item Build the plot: gray individual replicate points, colored group
#'     summary points, and (only if `p_val <= 0.05`) a dashed fitted curve
#'     plus dotted vertical gridlines every 4 h.
#'   \item Save the plot to `outfile` via `ggsave` (300 dpi, 9 x 5 in).
#' }
#'
#' @return A `ggplot` object (invisibly also saved to `outfile` as a side effect).
#'
#' @section Required packages:
#' `dplyr`, `ggplot2`, `lme4` (for `lmer`/`fixef`). `anova`, `cos`, `sin`,
#' `atan2`, `sqrt` are base R.
#'
#'
#' @usage
#' \dontrun{
#' plot_cosinor(
#'   df = intensity_data,
#'   group_name = "WT",
#'   color_hex = "#1F77B4",
#'   outfile = "WT_cosinor.png",
#'   protein_name = "BMAL1",
#'   scale = TRUE,
#'   y_max = NULL
#' )
#' }
plot_cosinor <- function(df, group_name, color_hex, outfile, protein_name, scale = FALSE, y_max = NULL) {
  df_group <- df %>%
    filter(Group == group_name) %>%
    mutate(Time_h = as.numeric(gsub("T","",Timepoint)))
  
  df_rep_avg <- df_group %>%
    group_by(Replicate, Time_h) %>%
    summarise(Mean_val = mean(Mean), .groups = "drop")
  
  cosinor_mixed <- lmer(
    Mean_val ~ cos(2*pi*Time_h/24) + sin(2*pi*Time_h/24) + 
      (1 | Replicate),
    data = df_rep_avg
  )
  
  null_model <- lmer(Mean_val ~ 1 + (1 | Replicate), data = df_rep_avg)
  p_val <- anova(null_model, cosinor_mixed)$`Pr(>Chisq)`[2]
  
  coefs <- fixef(cosinor_mixed)
  beta_cos <- coefs["cos(2 * pi * Time_h/24)"]
  beta_sin <- coefs["sin(2 * pi * Time_h/24)"]
  amplitude <- sqrt(beta_cos^2 + beta_sin^2)
  phase_rad <- atan2(-beta_sin, beta_cos)
  phase_h <- (phase_rad * 24)/(2*pi)
  if (phase_h < 0) phase_h <- phase_h + 24
  
  if (scale) {
    global_min <- min(df_rep_avg$Mean_val)
    global_max <- max(df_rep_avg$Mean_val)
    
    df_rep_avg <- df_rep_avg %>%
      mutate(Mean_val = (Mean_val - global_min) / (global_max - global_min))
    
    time_grid <- data.frame(Time_h = seq(0, 36, by = 0.1))
    raw_fit <- coefs["(Intercept)"] +
      coefs["cos(2 * pi * Time_h/24)"] * cos(2*pi*time_grid$Time_h/24) +
      coefs["sin(2 * pi * Time_h/24)"] * sin(2*pi*time_grid$Time_h/24)
    time_grid$Mean_fit <- (raw_fit - global_min) / (global_max - global_min)
    
  } else {
    time_grid <- data.frame(Time_h = seq(0, 36, by = 0.1))
    time_grid$Mean_fit <- coefs["(Intercept)"] +
      coefs["cos(2 * pi * Time_h/24)"] * cos(2*pi*time_grid$Time_h/24) +
      coefs["sin(2 * pi * Time_h/24)"] * sin(2*pi*time_grid$Time_h/24)
  }
  
  df_summary <- df_rep_avg %>%
    group_by(Time_h) %>%
    summarise(
      Mean_val = mean(Mean_val),
      .groups = "drop"
    )
  
  y_label <- ifelse(scale, "Scaled Intensity (0-1)", "Nuclear Intensity (A.U.)")
  
  p <- ggplot() +
    geom_point(
      data = df_rep_avg,
      aes(x = Time_h, y = Mean_val),
      color = "gray70",size = 4, alpha = 0.5
    ) +
    geom_point(
      data = df_summary,
      aes(x = Time_h, y = Mean_val),
      color = color_hex, size = 4, alpha = 0.5
    ) +
    scale_x_continuous(breaks = seq(0, 36, by = 4), limits = c(0, 36)) +
    theme_minimal(base_size = 25) +
    labs(
      title = paste(group_name, "-", protein_name),
      subtitle = sprintf("Cosinor p = %.3g | Amplitude = %.1f | Phase = %.1f h",
                         p_val, amplitude, phase_h),
      y = y_label,
      x = "Time (h)"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 25),
      plot.subtitle = element_text(size = 20, color = "gray30"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  if (scale) p <- p + ylim(0, 1)
  if (!is.null(y_max)) p <- p + ylim(0, y_max)
  
  if (p_val <= 0.05) {
    p <- p +
      geom_line(
        data = time_grid,
        aes(x = Time_h, y = Mean_fit),
        colour = "black", linetype = "dashed", linewidth = 1
      ) +
      geom_vline(
        xintercept = seq(0, 36, by = 4),
        linetype = "dotted", colour = "grey70"
      )
  }
  
  ggsave(outfile, p, dpi = 300, width = 9, height = 5)
  return(p)
}