# DNA Methylation Analysis - Scripts 
DNA methylation analysis for kidney and liver datasets
Liver Pipeline
- `sort_liver_files.py`: Generate control, 450k and EPIC directories and fill with respective raw data.
- `create_sheet-liver.py`: Create a uniform sample sheet from the original sample annotation excel file to be used by R.
- `liver_sample_sheet.csv`: Sample sheet for liver samples with phenotypes. 
- `liver_methylation_analysis.R`: Determines the methylation status for the CpGs in a combined array set, normalized via ComBat().
- `450k_liver_methylation_analysis.R`: Determines the methylation status for the CpGs in the 450K samples only, normalized via SWAN in minfi. 
- `EPIC_liver_methylation_analysis.R`: Determines the methylation status for the CpGs in the EPIC samples only, normalized via SWAN in minfi. 

`Note`: The R code for analysis includes optional functions for clustering via PCA and MDS, generating dendrograms and creating manhattan plots.
