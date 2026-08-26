# Permutation pipeline — RAIN + ECHO null-distribution filter
# -----------------------------------------------------------------------------
# Runs a timepoint-shuffle permutation test to build a null distribution of
# spurious rhythmic hits. Use apply_null_filter() afterwards to remove genes
# that appear too often by chance.
#
# Input matrices must be formatted for RAIN:
#   rows    = genes  (rownames = gene IDs)
#   columns = samples named  ZT_<timepoint>_<replicate>  (e.g. ZT_0_1, ZT_4_2)


# =============================================================================
## USER CONFIG
# =============================================================================

# Load your matrices here, then add them to perm_inputs.
# Each entry should be a matrix formatted as described above.
#
#   rain_counts_group1 <- readRDS("path/to/group1.rds")
#   rain_counts_group2 <- readRDS("path/to/group2.rds")

perm_inputs <- list(
  group1 = rain_counts_group_1,
  group2 = rain_counts_group_2
)

N_PERM           <- 1000
CHECKPOINT_EVERY <- 5 ## Limits the number of workers!!
SAVE_RDS         <- TRUE
WORKERS          <- max(1, parallel::detectCores() - 1)

PERIODS          <- c(22, 23, 24, 25)
RAIN_P_CUTOFF    <- 0.025
ECHO_P_CUTOFF    <- 0.025
HARM_CUT         <- 0.03
OVER_CUT         <- 0.15
NULL_PROP_CUTOFF <- 0.05

# =============================================================================
## END USER CONFIG
# =============================================================================


library(rain)
library(echo.find)
library(future.apply)


run_permutation <- function(rain_counts, group_label,
                            n_perm           = N_PERM,
                            seed             = 123,
                            periods          = PERIODS,
                            save_rds         = SAVE_RDS,
                            workers          = WORKERS,
                            checkpoint_every = CHECKPOINT_EVERY,
                            rain_p_cutoff    = RAIN_P_CUTOFF,
                            echo_p_cutoff    = ECHO_P_CUTOFF,
                            harm_cut         = HARM_CUT,
                            over_cut         = OVER_CUT) {
  
  gene_names <- rownames(rain_counts)
  
  timepoint_labels  <- gsub("_(\\d+)$", "", colnames(rain_counts))
  unique_timepoints <- unique(timepoint_labels)
  timepoint_groups  <- lapply(unique_timepoints,
                              function(tp) which(timepoint_labels == tp))
  
  prepare_echo_matrix <- function(mat, gnames) {
    colnames(mat) <- gsub("^ZT_(\\d+)_([1-3])$", "CT\\1.\\2", colnames(mat))
    df <- as.data.frame(mat)
    df$`Gene Name` <- gnames
    df <- df[, c("Gene Name", setdiff(names(df), "Gene Name"))]
    rownames(df) <- NULL
    df
  }
  
  rain_multi <- function(mat) {
    all_res <- do.call(rbind, lapply(periods, function(p) {
      res <- rain(t(mat), period = p, deltat = 4, nr.series = 3,
                  method = "independent", adjp.method = "ABH")
      df       <- as.data.frame(res)
      df$Gene  <- rownames(df)
      df
    }))
    ord     <- order(all_res$Gene, all_res$pVal)
    all_res <- all_res[ord, ]
    best    <- all_res[!duplicated(all_res$Gene), ]
    rownames(best) <- best$Gene
    best
  }
  
  one_perm <- function(i) {
    n_tp   <- length(timepoint_groups)
    rm_mat <- rain_counts
    rm_mat[] <- t(apply(rain_counts, 1, function(r) {
      ord <- sample(n_tp)
      r[unlist(timepoint_groups[ord])]
    }))
    colnames(rm_mat) <- colnames(rain_counts)
    
    rain_res  <- rain_multi(rm_mat)
    rain_hits <- rain_res$Gene[!is.na(rain_res$pVal) & rain_res$pVal < rain_p_cutoff]
    
    if (!length(rain_hits)) {
      message(sprintf("[%s] perm %d/%d | RAIN: 0 | ECHO: 0 | overlap: 0",
                      group_label, i, n_perm))
      return(character(0))
    }
    
    sub_mat  <- rm_mat[rain_hits, , drop = FALSE]
    echo_res <- echo_find(
      prepare_echo_matrix(sub_mat, rain_hits),
      begin = 0, end = 48, resol = 4, num_reps = 3,
      low = 22, high = 25,
      run_all_per = FALSE, paired = TRUE, rem_unexpr = FALSE,
      rem_unexpr_amt_below = 0, is_normal = FALSE,
      is_de_linear_trend = FALSE, is_smooth = FALSE,
      run_conf = FALSE, harm_cut = harm_cut, over_cut = over_cut
    )
    
    echo_hits <- echo_res$`Gene Name`[
      !is.na(echo_res$`P-Value`) &
        echo_res$`P-Value` < echo_p_cutoff &
        !echo_res$`Oscillation Type` %in% c("Repressed", "Overexpressed")
    ]
    
    both <- intersect(echo_hits[!is.na(echo_hits)], rain_hits)
    
    message(sprintf("[%s] perm %d/%d | RAIN: %d | ECHO: %d | overlap: %d",
                    group_label, i, n_perm,
                    length(rain_hits), length(echo_hits), length(both)))
    both
  }
  
  overlap_counts <- setNames(integer(length(gene_names)), gene_names)
  
  plan(multisession, workers = workers)
  on.exit(plan(sequential), add = TRUE)
  
  chunks <- split(seq_len(n_perm),
                  ceiling(seq_len(n_perm) / checkpoint_every))
  done <- 0L
  
  for (ch in chunks) {
    hits_list <- future_lapply(ch, one_perm, future.seed = seed + done)
    
    for (hits in hits_list)
      if (length(hits)) overlap_counts[hits] <- overlap_counts[hits] + 1L
    
    done <- done + length(ch)
    gc()
    
    if (save_rds) {
      tmp <- overlap_counts[overlap_counts > 0]
      saveRDS(tmp, sprintf("gene_overlap_counts_%s_checkpoint.rds", group_label))
      message(sprintf("[%s] checkpoint saved (%d/%d done)", group_label, done, n_perm))
    }
  }
  
  overlap_counts <- overlap_counts[overlap_counts > 0]
  if (length(overlap_counts)) {
    ord <- order(overlap_counts, decreasing = TRUE)
    permutation_df <- data.frame(
      Gene       = names(overlap_counts)[ord],
      Frequency  = as.integer(overlap_counts[ord]),
      Proportion = as.numeric(overlap_counts[ord]) / n_perm,
      stringsAsFactors = FALSE
    )
  } else {
    permutation_df <- data.frame(Gene = character(), Frequency = integer(),
                                 Proportion = numeric())
  }
  
  if (save_rds)
    saveRDS(permutation_df, sprintf("permutation_df_%s.rds", group_label))
  
  permutation_df
}


# Drops genes from filtered_df whose null-hit rate exceeds prop_cutoff
apply_null_filter <- function(filtered_df, overlap_df,
                              prop_cutoff = NULL_PROP_CUTOFF,
                              gene_col    = "Gene Name") {
  drop_genes <- overlap_df$Gene[overlap_df$Proportion > prop_cutoff]
  filtered_df[!filtered_df[[gene_col]] %in% drop_genes, , drop = FALSE]
}


# -----------------------------------------------------------------------------
# Run
# -----------------------------------------------------------------------------

verbose <- TRUE

permutation_dfs <- lapply(names(perm_inputs), function(g) {
  if (verbose) message(sprintf("[%s] Starting '%s' (n = %d rows)",
                               format(Sys.time(), "%H:%M:%S"), g, nrow(perm_inputs[[g]])))
  t0  <- Sys.time()
  res <- run_permutation(perm_inputs[[g]], group_label = g)
  if (verbose) message(sprintf("[%s] Done '%s' — %.1f sec, %d genes",
                               format(Sys.time(), "%H:%M:%S"), g,
                               as.numeric(difftime(Sys.time(), t0, units = "secs")),
                               nrow(res)))
  res
})
names(permutation_dfs) <- names(perm_inputs)

if (verbose) message(sprintf("[%s] All %d groups complete.",
                             format(Sys.time(), "%H:%M:%S"), length(permutation_dfs)))

