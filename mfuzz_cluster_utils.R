# mfuzz_cluster_utils.R
#
# Helper functions for extracting MFuzz cluster centroids and visualizing
# their peak phases on a 24-hour polar plot. Used for the analyses
# presented in this study.
#
# Typical usage:
#
#   # 1. Extract centroids across conditions into long format
#   centers_long <- extract_mfuzz_centroids(
#     clusters   = list(cl_wt, cl_ko),
#     conditions = c("WT", "KO")
#   )
#
#   # 2. Plot peak phases for a single clustering object
#   plot_mfuzz_peak_phases(cl_wt, title = "WT active")
#
# Dependencies: dplyr, tidyr, ggplot2, ggrepel

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)


# --- Data loading ------------------------------------------------------------
# Adjust these paths to point to your saved clustering objects (.rds or .RData)

# cl_wt_a   <- readRDS("data/cl_wt_active.rds")
# cl_p2ko_a <- readRDS("data/cl_p2ko_active.rds")
# cl_wt_q   <- readRDS("data/cl_wt_quiescent.rds")
# cl_p2ko_q <- readRDS("data/cl_p2ko_quiescent.rds")


# --- Functions ---------------------------------------------------------------

#' Extract MFuzz cluster centroids in long format
#'
#' Extracts cluster centroids from one or more MFuzz clustering objects,
#' assigns condition labels, orders clusters by the time point of their
#' maximum centroid expression, and returns the result in long format.
#'
#' Assumes colnames(clusters[[i]]$centers) are "CT<number>" strings
#' representing hours (e.g., "CT0", "CT4", ..., "CT20"). Rows of
#' $centers are clusters; columns are timepoints.
#'
#' @param clusters A list of MFuzz clustering objects.
#' @param conditions Character vector giving the condition label for each
#'   clustering object (must be the same length as `clusters`).
#'
#' @return A data frame with columns: cluster (factor, ordered by peak
#'   timepoint), condition, timepoint, and expression.
extract_mfuzz_centroids <- function(clusters, conditions) {

  stopifnot(length(clusters) == length(conditions))

  centers_list <- lapply(seq_along(clusters), function(i) {

    centers           <- as.data.frame(clusters[[i]]$centers)
    centers$cluster   <- paste0("C", seq_len(nrow(centers)))
    centers$condition <- conditions[[i]]

    # Order factor levels by the time point of maximum centroid expression.
    # peak_order[j] = row index whose peak is the j-th earliest, so
    # paste0("C", peak_order) gives level order from earliest to latest peak.
    peak_order        <- order(apply(clusters[[i]]$centers, 1, which.max))
    centers$cluster   <- factor(
      centers$cluster,
      levels = paste0("C", peak_order)
    )

    centers
  })

  bind_rows(centers_list) %>%
    pivot_longer(
      cols      = -c(cluster, condition),
      names_to  = "timepoint",
      values_to = "expression"
    )
}


#' Plot MFuzz cluster peak phases on a 24-hour polar plot
#'
#' Estimates the peak phase of each MFuzz cluster centroid using spline
#' interpolation over the first 24-hour cycle and displays the results
#' on a circular 24-hour plot.
#'
#' Assumes colnames(cl$centers) are "CT<number>" strings (e.g., "CT0",
#' "CT4"). Only timepoints with CT value <= 24 are used for interpolation;
#' if your data spans multiple cycles, values beyond CT24 are silently
#' ignored. peak_hour is returned in hours (0-24).
#'
#' @param cl An MFuzz clustering object.
#' @param nudge_clusters Character vector of cluster names (e.g., "C3")
#'   whose labels should be nudged by +1 hour in the pre-polar Cartesian
#'   space to reduce overlap. Note: nudging behaves uniformly in x-axis
#'   units, not strictly clockwise on the rendered circle.
#' @param title Plot title string. Defaults to NULL (no title).
#'
#' @return A ggplot object.
plot_mfuzz_peak_phases <- function(cl,
                                   nudge_clusters = character(),
                                   title          = NULL) {

  hours   <- as.numeric(gsub("CT", "", colnames(cl$centers)))
  n_clust <- nrow(cl$centers)

  centroid_peaks <- data.frame(
    cluster  = paste0("C", seq_len(n_clust)),
    peak_hour = apply(cl$centers, 1, function(x) {

      # Restrict interpolation to one 24-hour cycle
      keep     <- hours <= 24
      sp       <- splinefun(hours[keep], x[keep])
      opt      <- optimize(sp, interval = c(0, 24), maximum = TRUE)
      opt$maximum
    })
  )

  centroid_peaks$peak_hour <- round(centroid_peaks$peak_hour, 1)
  centroid_peaks$nudge_x   <- ifelse(
    centroid_peaks$cluster %in% nudge_clusters, 1, 0
  )

  ggplot(
    centroid_peaks,
    aes(x = peak_hour, y = 1, color = cluster, label = cluster)
  ) +

    # 24-hour circular reference ring
    geom_path(
      data        = data.frame(x = seq(0, 24, length.out = 300), y = 1),
      aes(x = x, y = y),
      color       = "grey80",
      linewidth   = 0.5,
      inherit.aes = FALSE
    ) +

    geom_point(size = 5) +

    geom_label_repel(
      size          = 5,
      fontface      = "bold",
      show.legend   = FALSE,
      force         = 10,
      max.overlaps  = Inf,
      box.padding   = 0.8,
      point.padding = 0.5,
      nudge_x       = centroid_peaks$nudge_x
    ) +

    coord_polar(theta = "x", start = 0) +

    # breaks stop at 20h to avoid label overlap at the 0/24 boundary
    scale_x_continuous(
      limits = c(0, 24),
      breaks = seq(0, 20, 4),
      labels = paste0(seq(0, 20, 4), "h")
    ) +

    scale_y_continuous(limits = c(0.7, 1)) +

    theme_minimal() +

    theme(
      axis.text.y        = element_blank(),
      axis.title         = element_blank(),
      legend.position    = "none",
      axis.text.x        = element_text(size = 14, face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    ) +

    labs(title = title)
}
