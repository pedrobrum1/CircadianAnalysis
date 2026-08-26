#' Plot expression rhythm comparison between two groups/genotypes
#'

#'
#' @param gene_name Character. Gene identifier to look up in `echo_df_1`/`echo_df_2`.
#' @param echo_df_1,echo_df_2 Data frames of ECHO/RAIN output for group 1 and
#'   group 2, respectively. Each must contain a gene identifier column matching
#'   `^Gene[ ._]?[Nn]ame$` (e.g. "GeneName", "Gene.Name"), columns matching
#'   `^Original TP` (raw replicate values per timepoint), columns matching
#'   `^Fitted TP` (fitted values per timepoint), a p-value column matching
#'   `^P[ ._-]?Value$` (ECHO p-value), a p-value column matching `^pVal$`
#'   (RAIN p-value, case-insensitive), and an `"Oscillation Type"` column.
#' @param begin,end Numeric. Start and end of the time range in hours.
#'   Defaults `0` and `48`.
#' @param resol Numeric. Time resolution/spacing between timepoints, in hours.
#'   Default `4`.
#' @param num_reps Integer. Number of replicates per timepoint. Default `3`.
#' @param label_1,label_2 Character. Legend labels for group 1 and group 2.
#'   Defaults `"Group 1"`, `"Group 2"`.
#' @param color_1,color_2 Character. Hex color codes for group 1 and group 2.
#'   Defaults `"#663399"`, `"#CC3300"`.
#' @param p_cutoff Numeric. Significance threshold applied to both RAIN and
#'   ECHO p-values for drawing the fitted curve. Default `0.025`.
#'
#'
#' A dashed spline-fitted curve is added for a group only if all of the
#' following hold: RAIN p-value and ECHO p-value are both non-`NA` and below
#' `p_cutoff`; and oscillation type is one of `"Harmonic"`, `"Damped"`, or
#' `"Forced"`.
#'
#' @return A `ggplot` object.
#'
#' @section Required packages:
#' `dplyr` (`group_by`, `summarise`, `bind_rows`), `ggplot2`. `spline` is
#' base R (`stats`).
#'
#' @examples
#' \dontrun{
#' plot_rhyhtms(
#'   gene_name = "Per2",
#'   echo_df_1 = echo_wt, echo_df_2 = echo_ko,
#'   label_1 = "WT", label_2 = "KO"
#' )
#' }
plot_rhyhtms <- function(gene_name, echo_df_1, echo_df_2,
                         begin = 0, end = 48, resol = 4, num_reps = 3,
                         label_1 = "Group 1", label_2 = "Group 2",
                         color_1 = "#663399", color_2 = "#CC3300",
                         p_cutoff = 0.025) {
  
  timepoints <- seq(begin, end, by = resol)
  n_vals <- length(timepoints) * num_reps
  
  # Sort TP columns numerically by the timepoint number in their name
  # (e.g. "Original TP2", "Original TP10" -> ordered 2, 10, not 10, 2)
  sort_by_tp_number <- function(cols) {
    if (length(cols) == 0) return(cols)
    nums <- as.numeric(gsub("\\D", "", cols))
    cols[order(nums)]
  }
  
  get_data <- function(echo_df, label, color) {
    gene_col    <- grep("^Gene[ ._]?[Nn]ame$", names(echo_df), value = TRUE)[1]
    echo_p_col  <- grep("^P[ ._-]?Value$", names(echo_df), value = TRUE)[1]
    rain_p_col  <- grep("^pVal$", names(echo_df), value = TRUE, ignore.case = TRUE)[1]
    fitted_cols <- sort_by_tp_number(grep("^Fitted TP", names(echo_df), value = TRUE))
    raw_cols    <- sort_by_tp_number(grep("^Original TP", names(echo_df), value = TRUE))
    
    if (is.na(gene_col))   stop("Could not find a gene identifier column (expected something like 'GeneName').")
    if (is.na(echo_p_col)) stop("Could not find an ECHO p-value column (expected something like 'P-Value').")
    if (is.na(rain_p_col)) stop("Could not find a RAIN p-value column (expected 'pVal').")
    
    gene_row <- echo_df[echo_df[[gene_col]] == gene_name, , drop = FALSE]
    
    if (nrow(gene_row) > 1) {
      stop(sprintf("Gene '%s' matched %d rows in the ECHO data frame; expected exactly 1.",
                   gene_name, nrow(gene_row)))
    }
    
    if (nrow(gene_row) == 0) {
      return(list(
        pval.ECHO = NA,
        pval.rain = NA,
        osc.type  = "Not detected",
        raw_vals  = rep(0, n_vals),
        fit_vals  = rep(0, length(timepoints)),
        label     = label,
        color     = color
      ))
    }
    
    if (length(raw_cols) != n_vals) {
      stop(sprintf(
        "Gene '%s': found %d 'Original TP*' columns, expected %d (length(timepoints) * num_reps).",
        gene_name, length(raw_cols), n_vals
      ))
    }
    if (length(fitted_cols) != length(timepoints)) {
      stop(sprintf(
        "Gene '%s': found %d 'Fitted TP*' columns, expected %d (length(timepoints)).",
        gene_name, length(fitted_cols), length(timepoints)
      ))
    }
    
    list(
      pval.ECHO = signif(as.numeric(gene_row[[echo_p_col]]), 3),
      pval.rain = signif(as.numeric(gene_row[[rain_p_col]]), 3),
      osc.type  = as.character(gene_row[["Oscillation Type"]]),
      raw_vals  = as.numeric(gene_row[, raw_cols]),
      fit_vals  = as.numeric(gene_row[, fitted_cols]),
      label     = label,
      color     = color
    )
  }
  
  d1 <- get_data(echo_df_1, label_1, color_1)
  d2 <- get_data(echo_df_2, label_2, color_2)
  
  make_df <- function(raw_vals, timepoints, num_reps, label) {
    raw_times <- rep(timepoints, each = num_reps)
    df_raw <- data.frame(Time = raw_times, Expression = raw_vals, Genotype = label)
    df_raw %>%
      group_by(Time, Genotype) %>%
      summarise(mean_expr = mean(Expression), sd_expr = sd(Expression), .groups = "drop")
  }
  
  df_summary <- bind_rows(
    make_df(d1$raw_vals, timepoints, num_reps, d1$label),
    make_df(d2$raw_vals, timepoints, num_reps, d2$label)
  )
  
  fmt_p <- function(p) if (is.na(p)) "NA" else sprintf("%.3g", p)
  
  p <- ggplot(df_summary, aes(x = Time, y = mean_expr, color = Genotype)) +
    geom_point(size = 4) +
    geom_errorbar(aes(ymin = mean_expr - sd_expr, ymax = mean_expr + sd_expr), width = 1) +
    scale_color_manual(values = setNames(c(d1$color, d2$color), c(d1$label, d2$label))) +
    scale_x_continuous(breaks = timepoints) +
    theme_minimal(base_size = 25) +
    labs(
      title = gene_name,
      subtitle = sprintf(
        "%s: RAIN p = %s | ECHO p = %s | Type = %s\n%s: RAIN p = %s | ECHO p = %s | Type = %s",
        d1$label, fmt_p(d1$pval.rain), fmt_p(d1$pval.ECHO), d1$osc.type,
        d2$label, fmt_p(d2$pval.rain), fmt_p(d2$pval.ECHO), d2$osc.type
      ),
      x = "Time (CT, h)", y = "Expression\n(vst-counts)", color = "Genotype"
    ) +
    theme(
      axis.text = element_text(size = 25, face = "bold"),
      plot.title = element_text(size = 30, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 25, color = "gray30"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      plot.margin = margin(15, 15, 15, 15)
    )
  
  for (d in list(d1, d2)) {
    if (!is.na(d$pval.rain) &&
        !is.na(d$pval.ECHO) &&
        d$pval.rain < p_cutoff &&
        d$pval.ECHO < p_cutoff &&
        d$osc.type %in% c("Harmonic", "Damped", "Forced")) {
      
      interp_time <- seq(begin, end, by = 0.1)
      fit_spline  <- spline(x = timepoints, y = d$fit_vals, xout = interp_time)
      df_fit      <- data.frame(Time = fit_spline$x, Fit = fit_spline$y, Genotype = d$label)
      p <- p + geom_line(data = df_fit, aes(x = Time, y = Fit, color = Genotype),
                         linetype = "dashed", linewidth = 1)
    }
  }
  
  return(p)
}
