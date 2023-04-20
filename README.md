# DNA Methylation Analysis - Scripts 
### DNA methylation analysis for kidney and liver datasets

## Liver Scripts
- `sort_liver_files.py`: Generate control, 450k and EPIC directories and fill with respective raw liver data.
- `create_sheet-liver.py`: Create a uniform liver sample sheet from the original sample annotation excel file to be used by R. 
- `450k_liver_methylation_analysis.R`: Determines the methylation status for the CpGs in the 450K samples only, normalized via SWAN in minfi. 
- `denconv_boxplots.R`: Generates boxplots when provided the liver sample sheet, containing cell type percentages (deconvolution). The deconvolution is performed in the 450k_liver_methylation_analysis.R script using the EpiDISH package. This outputs the designed sample sheet to be used by this script. 

`Note`: The R code for analysis includes optional functions for clustering via PCA and MDS, generating dendrograms and creating manhattan plots.

## Kidney Scripts
- `sort_kidney_files.py`: Generate control, 450k and EPIC directories and fill with respective raw kidney data.
- `create_sheet-kidney.py`: Create a uniform kidney sample sheet from the original sample annotation excel file to be used by R.
- `EPIC_kidney_methylation_analysis.R`: Determines the methylation status for the CpGs in the EPIC samples only, normalized via SWAN in minfi. 

`Note`: Samples sheets with only paired samples are denoted by a 'paired' prefix for the sample sheet. 

## Other Scripts
- `calc_bio_age.R`: Calculate the biological (DNAm) age given both the beta values and the sample sheet (inputs) using the methylclock package. 
- `methylation_heatmaps.R`: Generates heatmaps to display changes in methylation for paired samples across timepoints given both the beta values and the sample sheet as inputs. 
- `line_graph_DMRs.R`: Generates line graphs to illustrate changes in methylation status at specific CpG sites given the beta values, sample sheet, DMPs and DMRs as inputs. 
- `filter_DMPs.py`: Generates sig_DMPs.csv file, consisting of only CpG sites with an adjusted p-value of 0.05 or less. 

## Supporting Files
- `liver_sample_sheet.csv`: Sample sheet for liver samples with phenotypes.
- `kidney_sample_sheet.csv`: Sample sheet for ALL kidney samples with phenotypes. 
- `kidney_sample_sheet_K3_K12.csv`: Sample sheet for kidney samples with phenotypes. This is specific sample sheet for ALL K3 and K12 samples was generated via a modified version of the create_sheet-kidney.py script.
- `48639-non-specific-probes-Illumina450k.csv`: Identified as non-specific probes, removed during the analysis. 
- `HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt`: Identified as probes with potential hybridization issues, removed during the analysis. 
