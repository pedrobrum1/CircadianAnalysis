# Reproducibility files

This repository contains the code and processed data required to reproduce the analyses and figures from:

**PER2- and state-dependent transcriptional programs gate neural stem cell proliferation with niche-specific circadian autonomy**  
DOI: `xxxx.xxxx`

Raw sequencing data can be downloaded from:

ENA accession: `PRJEB122945`

Immunofluorescence datasets can be made available upon request.

For questions, please contact:

- pedro.ozorio@imba.oeaw.ac.at
- pedro.ozorio.brum@univie.ac.at

## Immunofluorescence data _in vitro_

### Script

- `plot_cosinor.R` — fits and plots 24-hour cosinor models for immunofluorescence intensity data.

### Input data

- `BMAL1_if.rds`
- `VENUS_if.rds`

## RNA-seq data

### Input data

Rhythmicity analysis derived the following objects:

#### Active cells

- `active_sgz_wt_combined.rds`
- `active_sgz_wt_filtered.rds`
- `active_sgz_p2ko_combined.rds`
- `active_sgz_p2ko_filtered.rds`
- `active_svz_wt_combined.rds`
- `active_svz_wt_filtered.rds`

#### Quiescent cells

- `quiescent_sgz_wt_combined.rds`
- `quiescent_sgz_wt_filtered.rds`
- `quiescent_sgz_p2ko_combined.rds`
- `quiescent_sgz_p2ko_filtered.rds`
- `quiescent_svz_wt_combined.rds`
- `quiescent_svz_wt_filtered.rds`

### Object names

- `_combined` — rhythmicity-analysis results for all transcripts.
- `_filtered` — results restricted to transcripts with `p < 0.025` in both RAIN and ECHO, and a frequency of less than 5% after applying `permutation_filtering.R`.

### Scripts
- `permutation_filtering.R` - permutation pipeline
- `plot_rhythms.R` — plots results from the rhythmicity analysis of bulk RNA-seq data.

