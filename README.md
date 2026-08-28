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

- `plot_cosinor.R` — fits and plots 24-hour cosinor models for immunofluorescence intensity data (Fig. 2)

### Input data
* Objects contain nuclear intensity values ($Mean, $IndDen, $RawIntDen) objtained through DAPI segmentation utilizing StarDist 2D ([Weigert et al., 2022]https://github.com/stardist/stardist-imagej).
* Includes $Group: SGZ_A, SGZ_Q, SVZ_A, SVZ_Q / $Replicate: BR1, BR2, BR3 / $File: original .vsi (images) and .csv (quantification) file names 

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
* Objects contain RAIN ([Thaben et al., 2014](https://www.bioconductor.org/packages//3.18/bioc/html/rain.html)) and ECHO ([De los Santos et al., 2020], https://github.com/delosh653/ECHO) results, vst-normalized counts, and ECHO-fitted counts.
- `_combined` — rhythmicity-analysis results for all transcripts.
- `_filtered` — results restricted to transcripts with `p < 0.025` in both RAIN and ECHO, and a frequency of less than 5% after applying `permutation_filtering.R`.

### Scripts
- `permutation_filtering.R` - permutation pipeline
- `plot_rhythms.R` — plots results from the rhythmicity analysis of bulk RNA-seq data. (Fig. 2 and Supp. Fig 5 & 7)
- `plot_phase_heatmap.R` - heatmap ordered by phase (Fig. 2)
- `mfuzz_cluster_utils.R` - extracts centroids and plots clusters based on phase (Fig. 2)

