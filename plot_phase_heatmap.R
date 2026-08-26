#' Plot a phase-ordered gene-expression heatmap
#'
#' Averages replicate expression values measured at the same timepoint,
#' smooths each gene's temporal expression profile using a cubic smoothing
#' spline, z-score scales each gene across time, and displays the result as a
#' heatmap. Genes found in `df_results` are ordered by their reported phase
#' (`Hours Shifted`). A second heatmap can optionally be displayed side by side.
#'
#' @param genes Character vector of gene names to include in the heatmap.
#'   Genes are plotted in phase order when they are present in `df_results`;
#'   other requested genes are appended afterward in their original order.
#' @param df_results Data frame containing at least the columns `Gene Name`,
#'   `Hours Shifted`, `P-Value`, and `pVal`. `Gene Name` is used to match
#'   genes, and `Hours Shifted` is used to order genes by phase. Genes
#'   meeting the p-value and oscillation-type criteria are marked with an
#'   asterisk in the heatmap row labels.
#' @param df_counts Numeric matrix or data frame containing expression values,
#'   with genes as row names and timepoints as column names. Replicate columns
#'   measured at the same timepoint must have identical column names.
#' @param title Character string used as the title of the first heatmap.
#'   Defaults to an empty string.
#' @param df_counts2 Optional second numeric expression matrix or data frame
#'   with the same structure as `df_counts`. If supplied, a second heatmap is
#'   plotted next to the first one. Defaults to `NULL`.
#' @param df_results2 Optional second results data frame containing a `Gene Name`
#'   column. It is used to mark genes in the second heatmap. If `NULL`, the
#'   first heatmap's `df_results` is used instead. Defaults to `NULL`.
#' @param title2 Character string used as the title of the second heatmap.
#'   Used only when `df_counts2` is supplied.
#' @param rownames Logical. Whether gene names should be displayed as heatmap
#'   row labels. Defaults to `FALSE`.
#' @param colnames Logical. Whether timepoints should be displayed as heatmap
#'   column labels. Defaults to `TRUE`.
#'
#' @section Timepoint requirements:
#' Column names of `df_counts` and, when supplied, `df_counts2` are interpreted
#' as timepoints and **must be numeric or coercible to numeric**. Valid examples
#' include `"0"`, `"4"`, `"8"`, and `"12"`.
#'
#' Each replicate measured at the same timepoint must use the same column name.
#' For example, three replicates at timepoint 4 should have three columns named
#' `"4"`. These columns are averaged using `rowMeans()`.
#'
#' @section Significance marking:
#' A gene is marked with an asterisk in the row labels only if **both**:
#' \itemize{
#'   \item ECHO `` `P-Value` `` < 0.025, and
#'   \item RAIN `pVal` < 0.025,
#' }
#' and, when an `` `Oscillation Type` `` column is present, the gene's
#' oscillation type is one of `"Damped"`, `"Forced"`, or `"Harmonic"`.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Restricts `df_results` to the requested genes and orders them by
#'     `Hours Shifted`.
#'   \item Retains genes found in the row names of the count matrix.
#'   \item Averages replicate columns sharing the same timepoint.
#'   \item Removes genes whose averaged expression profile contains missing
#'     values.
#'   \item Smooths each gene's temporal profile using
#'     `smooth.spline(..., spar = 0.2)`.
#'   \item Z-score scales each gene across timepoints.
#'   \item Displays the scaled expression values using a fixed blue-white-red
#'     colour scale ranging from -3 to 3.
#' }
#'
#' When `df_counts2` is supplied, the second matrix is processed independently
#' using the same requested gene order and displayed beside the first heatmap.
#'
#' @return The heatmap is drawn to the current graphics device. The function
#'   does not return a data frame or a reusable heatmap object.
#'
#' @section Required packages:
#' `pheatmap`, `grid`, and `gridExtra` (called via namespace, not attached).
#'
#' @note Genes absent from the count matrix, or genes removed because their
#' averaged expression profile contains missing values, are not displayed.
#' The number of genes shown is reported in each heatmap title.
#'
#' @examples
#' \dontrun{
#' plot_phase_heatmap(
#'   genes = c("Per2", "Bmal1", "Cry1"),
#'   df_results = echo_wt,
#'   df_counts = counts_wt,
#'   title = "WT"
#' )
#'
#' plot_phase_heatmap(
#'   genes = c("Per2", "Bmal1", "Cry1"),
#'   df_results = echo_wt,
#'   df_counts = counts_wt,
#'   df_counts2 = counts_ko,
#'   df_results2 = echo_ko,
#'   title = "WT",
#'   title2 = "KO",
#'   rownames = TRUE
#' )
#' }
plot_phase_heatmap <- function(genes, df_results, df_counts, title = "",
                               df_counts2 = NULL, df_results2 = NULL, title2 = "",
                               rownames = FALSE, colnames = TRUE) {
  
  df_results <- df_results[df_results$`Gene Name` %in% genes, ]
  if (nrow(df_results) == 0) stop("No genes from the list found in df_results.")
  df_results <- df_results[order(df_results$`Hours Shifted`), ]
  genes      <- df_results$`Gene Name`
  
  avg_matrix <- function(df_counts_in, genes_in) {
    gene_counts  <- df_counts_in[match(genes_in, rownames(df_counts_in)), ]
    unique_times <- unique(colnames(gene_counts))
    mat <- as.data.frame(sapply(unique_times, function(tp) {
      rowMeans(gene_counts[, colnames(gene_counts) == tp, drop = FALSE], na.rm = TRUE)
    }))
    mat[complete.cases(mat), ]
  }
  
  smooth_mat <- function(mat) {
    t(apply(mat, 1, function(x) smooth.spline(seq_along(x), x, spar = 0.2)$y))
  }
  
  avg1 <- smooth_mat(avg_matrix(df_counts, genes))
  avg1_scaled <- t(scale(t(avg1)))
  
  if (!is.null(df_counts2)) {
    if (!is.null(df_results2)) {
      df_results2 <- df_results2[df_results2$`Gene Name` %in% genes, ]
    } else {
      df_results2 <- df_results
    }
    
    avg2 <- smooth_mat(avg_matrix(df_counts2, genes))
    avg2_scaled <- t(scale(t(avg2)))
  }
  
  make_labels <- function(mat, df_res) {
    sapply(rownames(mat), function(g) {
      row_match <- df_res[df_res$`Gene Name` == g, , drop = FALSE]
      if (nrow(row_match) == 0) return(g)
      
      p_ok <- any(
        is.finite(row_match$`P-Value`) &
          is.finite(row_match$pVal) &
          row_match$`P-Value` < 0.025 &
          row_match$pVal < 0.025
      )
      
      osc_ok <- !("Oscillation Type" %in% colnames(row_match)) ||
        any(row_match$`Oscillation Type` %in% c("Damped", "Forced", "Harmonic"))
      if (p_ok && osc_ok) paste0(g, " *") else g
    })
  }
  
  col_fun <- colorRampPalette(c("#1f77b4", "white", "#CC3300"))
  breaks  <- seq(-3, 3, length.out = 101)
  
  p1 <- pheatmap::pheatmap(
    avg1_scaled,
    cluster_rows  = FALSE,
    cluster_cols  = FALSE,
    color         = col_fun(100),
    breaks        = breaks,
    scale         = "none",
    show_rownames = rownames,
    show_colnames = colnames,
    labels_row    = make_labels(avg1_scaled, df_results),
    main          = paste0(title, " (n = ", nrow(avg1_scaled), " genes)"),
    border_color  = NA,
    use_raster    = FALSE,
    fontsize      = 40,
    silent        = TRUE
  )
  
  if (!is.null(df_counts2)) {
    p2 <- pheatmap::pheatmap(
      avg2_scaled,
      cluster_rows  = FALSE,
      cluster_cols  = FALSE,
      color         = col_fun(100),
      breaks        = breaks,
      scale         = "none",
      show_rownames = rownames,
      show_colnames = colnames,
      labels_row    = make_labels(avg2_scaled, df_results2),
      main          = paste0(title2, " (n = ", nrow(avg2_scaled), " genes)"),
      border_color  = NA,
      use_raster    = FALSE,
      fontsize      = 40,
      silent        = TRUE
    )
    
    gridExtra::grid.arrange(p1$gtable, p2$gtable, ncol = 2)
    
  } else {
    grid::grid.draw(p1$gtable)
  }
}