# dna_methylation_analysis
DNA methylation analysis for kidney and liver datasets
Liver Pipeline
- `sort_liver_files.py`: Generate control, 450k and EPIC directories and fill with respective raw data.
- `create_sheet-liver.py`: Create a uniform sample sheet from the original sample annotation excel file to be used by R.
- `liver_sample_sheet.csv`: Sample sheet for liver samples with phenotypes. 
- `liver_methylation_analysis.R`: Determines the methylation status for a combined array set, normalized via ComBat().
