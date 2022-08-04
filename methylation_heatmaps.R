#Heatmap for DMPs 
rm(list = ls())
dev.off()
library(rio)
library(tidyverse)
library(ggplot2)
library(lumi)
#Liver heatmaps 
liver_betas <- import('/Users/adrianharris/Desktop/large_files_liver_combined/beta_values.csv')
liver_betas <- liver_betas %>% rename(ProbeID = V1)
row.names(liver_betas) <- liver_betas$ProbeID

liver_pheno <- import('/Users/adrianharris/Documents/dna_methylation_analysis/liver_sample_sheet.csv')

#Removal of control samples 
liver_pheno <- subset(liver_pheno, sample_group != 'control')
liver_pheno <- subset(liver_pheno, sample_group != 'Control')

#Removal of one unpaired sample
liver_pheno <- liver_pheno[!(liver_pheno$sample_name == 'V037L1'),]

#Removal of samples not on server 
liver_pheno <- liver_pheno[!(liver_pheno$Basename == '203751390020_R08C01' | liver_pheno$Basename == '203751390020_R01C01'),]

#Drop DCD samples
liver_pheno = liver_pheno[!(liver_pheno$donor_type == 'DCD'),]

#Update sample_group column 
for (row in 1:nrow(liver_pheno)) {
  if (liver_pheno[row, 'sample_group'] == 'High Injury') {
    injury <- 'High'
  } else {
    injury <- 'Low'
  }
  liver_pheno[row, 'sample_group'] <- injury
}

rm(row)
rm(injury)

#Addition and filling of the condition column 
liver_pheno$condition <- 'NA'

#Filling the condition column
for (row in 1:nrow(liver_pheno)) {
  if (liver_pheno[row, 'donor_type'] == 'DD' & liver_pheno[row, 'sample_group'] == 'High' &  liver_pheno[row, 'collection'] == 'L1') {
    liver_pheno[row, 'condition'] <- "DD_HI_L1"
  } else if (liver_pheno[row, 'donor_type'] == 'DD' & liver_pheno[row, 'sample_group'] == 'High' &  liver_pheno[row, 'collection'] == 'L2') {
    liver_pheno[row, 'condition'] <- "DD_HI_L2"
  } else if (liver_pheno[row, 'donor_type'] == 'DD' & liver_pheno[row, 'sample_group'] == 'Low' &  liver_pheno[row, 'collection'] == 'L1') {
    liver_pheno[row, 'condition'] <- "DD_LI_L1"
  } else if (liver_pheno[row, 'donor_type'] == 'DD' & liver_pheno[row, 'sample_group'] == 'Low' &  liver_pheno[row, 'collection'] == 'L2') {
    liver_pheno[row, 'condition'] <- "DD_LI_L2"
  } else if (liver_pheno[row, 'donor_type'] == 'LD' & liver_pheno[row, 'sample_group'] == 'Low' &  liver_pheno[row, 'collection'] == 'L1') {
    liver_pheno[row, 'condition'] <- "LD_LI_L1"
  } else if (liver_pheno[row, 'donor_type'] == 'LD' & liver_pheno[row, 'sample_group'] == 'Low' &  liver_pheno[row, 'collection'] == 'L2') {
    liver_pheno[row, 'condition'] <- "LD_LI_L2"
  }
}

rm(row)

probes1 <- import('/Users/adrianharris/Desktop/Mas_lab/meeting_07292022/liver_csvs_combat0728/DD_HI_L1-DD_HI_L2_DMPs_sig.csv')
probes1 <- probes1$Name

probes2 <- import('/Users/adrianharris/Desktop/Mas_lab/meeting_07292022/liver_csvs_combat0728/DD_HI_L1-DD_LI_L1_DMPs_sig.csv')
probes2 <- probes2$Name

probes3 <- import('/Users/adrianharris/Desktop/Mas_lab/meeting_07292022/liver_csvs_combat0728/DD_HI_L1-LD_LI_L1_DMPs_sig.csv')
probes3 <- probes3$Name

probes4 <- import('/Users/adrianharris/Desktop/Mas_lab/meeting_07292022/liver_csvs_combat0728/DD_HI_L2-LD_LI_L2_DMPs_sig.csv')
probes4 <- probes4$Name

probes5 <- import('/Users/adrianharris/Desktop/Mas_lab/meeting_07292022/liver_csvs_combat0728/DD_LI_L1-DD_LI_L2_DMPs_sig.csv')
probes5 <- probes5$Name

probes6 <- import('/Users/adrianharris/Desktop/Mas_lab/meeting_07292022/liver_csvs_combat0728/DD_LI_L1-LD_LI_L1_DMPs_sig.csv')
probes6 <- probes6$Name

probes7 <- import('/Users/adrianharris/Desktop/Mas_lab/meeting_07292022/liver_csvs_combat0728/DD_LI_L2-LD_LI_L2_DMPs_sig.csv')
probes7 <- probes7$Name

#Filter down datasets 
liver_pheno1 <- liver_pheno[(liver_pheno$condition == 'DD_HI_L1' | liver_pheno$condition == 'DD_HI_L2'),]
liver_pheno1 <- liver_pheno1 %>% arrange(condition)

liver_pheno2 <- liver_pheno[(liver_pheno$condition == 'DD_HI_L1' | liver_pheno$condition == 'DD_LI_L1'),]
liver_pheno2 <- liver_pheno2 %>% arrange(condition)

liver_pheno3 <- liver_pheno[(liver_pheno$condition == 'DD_HI_L1' | liver_pheno$condition == 'LD_LI_L1'),]
liver_pheno3 <- liver_pheno3 %>% arrange(condition)

liver_pheno4 <- liver_pheno[(liver_pheno$condition == 'DD_HI_L2' | liver_pheno$condition == 'LD_LI_L2'),]
liver_pheno4 <- liver_pheno4 %>% arrange(condition)

liver_pheno5 <- liver_pheno[(liver_pheno$condition == 'DD_LI_L1' | liver_pheno$condition == 'DD_LI_L2'),]
liver_pheno5 <- liver_pheno5 %>% arrange(condition)

liver_pheno6 <- liver_pheno[(liver_pheno$condition == 'DD_LI_L1' | liver_pheno$condition == 'LD_LI_L1'),]
liver_pheno6 <- liver_pheno6 %>% arrange(condition)

liver_pheno7 <- liver_pheno[(liver_pheno$condition == 'DD_LI_L2' | liver_pheno$condition == 'LD_LI_L2'),]
liver_pheno7 <- liver_pheno7 %>% arrange(condition)

#Filter down betas given the dataset and match order
liver_betas1 <- liver_betas[,match(liver_pheno1$Basename, colnames(liver_betas))]
colnames(liver_betas1) <- liver_pheno1$sample_name
liver_betas2 <- liver_betas[,match(liver_pheno2$Basename, colnames(liver_betas))]
colnames(liver_betas2) <- liver_pheno2$sample_name
liver_betas3 <- liver_betas[,match(liver_pheno3$Basename, colnames(liver_betas))]
colnames(liver_betas3) <- liver_pheno3$sample_name
liver_betas4 <- liver_betas[,match(liver_pheno4$Basename, colnames(liver_betas))]
colnames(liver_betas4) <- liver_pheno4$sample_name
liver_betas5 <- liver_betas[,match(liver_pheno5$Basename, colnames(liver_betas))]
colnames(liver_betas5) <- liver_pheno5$sample_name
liver_betas6 <- liver_betas[,match(liver_pheno6$Basename, colnames(liver_betas))]
colnames(liver_betas6) <- liver_pheno6$sample_name
liver_betas7 <- liver_betas[,match(liver_pheno7$Basename, colnames(liver_betas))]
colnames(liver_betas7) <- liver_pheno7$sample_name

liver_betas1 <- liver_betas1[rownames(liver_betas1) %in% probes1,]
liver_betas2 <- liver_betas2[rownames(liver_betas2) %in% probes2,]
liver_betas3 <- liver_betas3[rownames(liver_betas3) %in% probes3,]
liver_betas4 <- liver_betas4[rownames(liver_betas4) %in% probes4,]
liver_betas5 <- liver_betas5[rownames(liver_betas5) %in% probes5,]
liver_betas6 <- liver_betas6[rownames(liver_betas6) %in% probes6,]
liver_betas7 <- liver_betas7[rownames(liver_betas7) %in% probes7,]

#install.packages('gplots')
library(gplots)
pdf('/Users/adrianharris/Desktop/liver_heatmaps.pdf')
liver_betas1 <- liver_betas1[sample(nrow(liver_betas1), 1000), ]
heatmap.2(x=as.matrix(liver_betas1), 
          Colv=FALSE, 
          dendrogram="none",
          scale="row",
          col="bluered",
          trace="none",
          labRow=TRUE,
          main="DD_HI_L1-DD_HI_L2",
          ylab="Probes",
          xlab="Samples")

heatmap.2(x=as.matrix(liver_betas2), 
          Colv=FALSE, 
          dendrogram="none",
          scale="row",
          col="bluered",
          trace="none",
          labRow=TRUE,
          main="DD_HI_L1-DD_LI_L1",
          ylab="Probes",
          xlab="Samples")

heatmap.2(x=as.matrix(liver_betas3), 
          Colv=FALSE, 
          dendrogram="none",
          scale="row",
          col="bluered",
          trace="none",
          labRow=TRUE,
          main="DD_HI_L1-LD_LI_L1",
          ylab="Probes",
          xlab="Samples")

liver_betas4 <- liver_betas4[sample(nrow(liver_betas4), 1000), ]
heatmap.2(x=as.matrix(liver_betas4), 
          Colv=FALSE, 
          dendrogram="none",
          scale="row",
          col="bluered",
          trace="none",
          labRow=TRUE,
          main="DD_HI_L2-LD_LI_L2",
          ylab="Probes",
          xlab="Samples")

liver_betas5 <- liver_betas5[sample(nrow(liver_betas5), 1000), ]
heatmap.2(x=as.matrix(liver_betas5), 
          Colv=FALSE, 
          dendrogram="none",
          scale="row",
          col="bluered",
          trace="none",
          labRow=TRUE,
          main="DD_LI_L1-DD_LI_L2",
          ylab="Probes",
          xlab="Samples")

liver_betas6 <- liver_betas6[sample(nrow(liver_betas6), 1000), ]
heatmap.2(x=as.matrix(liver_betas6), 
          Colv=FALSE, 
          dendrogram="none",
          scale="row",
          col="bluered",
          trace="none",
          labRow=TRUE,
          main="DD_LI_L1-LD_LI_L1",
          ylab="Probes",
          xlab="Samples")

heatmap.2(x=as.matrix(liver_betas7), 
          Colv=FALSE, 
          dendrogram="none",
          scale="row",
          col="bluered",
          trace="none",
          labRow=TRUE,
          main="DD_LI_L2-LD_LI_L2",
          ylab="Probes",
          xlab="Samples")

dev.off()

#Kidney heatmaps
rm(list = ls())
kidney_betas <- import('/Users/adrianharris/Desktop/kidney_combined_beta_dir/beta_values.csv')
kidney_betas <- kidney_betas %>% rename(ProbeID = V1)
rownames(kidney_betas) <- kidney_betas$ProbeID

kidney_pheno <- import('/Users/adrianharris/Documents/dna_methylation_analysis/paired_kidney_sample_sheet.csv')

#Must remove outlier sample, 203504430032_R01C01 (and its paired sample 203504430032-R02C01)
kidney_pheno <- kidney_pheno[!(kidney_pheno$Basename == '203504430032_R01C01' | kidney_pheno$Basename == '203504430032_R02C01'),]

#Drop other outlier sample 
kidney_pheno <- kidney_pheno[!(kidney_pheno$sample_id == 'KUT4_K2' | kidney_pheno$sample_id == 'KUT4_K1'),]

#Addition of condition column 
kidney_pheno$condition <- 'NA'

kidney_pheno <- subset(kidney_pheno, select = -c(eGFR_1month, eGFR_24month))
kidney_pheno <- kidney_pheno %>% rename('eGFR' = 'eGFR_12month')

#Filling of condition column
for (row in 1:nrow(kidney_pheno)) {
  if (kidney_pheno[row, 'time'] == 'K1' & kidney_pheno[row, 'eGFR'] == 'High') {
    kidney_pheno[row, 'condition'] <- "K1_High"
  } else if (kidney_pheno[row, 'time'] == 'K1' & kidney_pheno[row, 'eGFR'] == 'Low') {
    kidney_pheno[row, 'condition'] <- "K1_Low"
  } else if (kidney_pheno[row, 'time'] == 'K2' & kidney_pheno[row, 'eGFR'] == 'High') {
    kidney_pheno[row, 'condition'] <- "K2_High"
  } else if (kidney_pheno[row, 'time'] == 'K2' & kidney_pheno[row, 'eGFR'] == 'Low') {
    kidney_pheno[row, 'condition'] <- "K2_Low"
  }
}

#if condition is empty, drop it
kidney_pheno <- kidney_pheno[!(kidney_pheno$condition == 'NA'),]

#Load probes 
rm(row)

probes1 <- import('/Users/adrianharris/Desktop/kidney_csvs_combined0728/eGFR_12month-K1_Low-K1_High_DMPs_sig.csv')
probes1 <- probes1$Name

probes2 <- import('/Users/adrianharris/Desktop/kidney_csvs_combined0728/eGFR_12month-K2_Low-K2_High_DMPs_sig.csv')
probes2 <- probes2$Name

#Filter down datasets 
kidney_pheno1 <- kidney_pheno[(kidney_pheno$condition == 'K1_Low' | kidney_pheno$condition == 'K1_High'),]
kidney_pheno1 <- kidney_pheno1 %>% arrange(condition)

kidney_pheno2 <- kidney_pheno[(kidney_pheno$condition == 'K2_Low' | kidney_pheno$condition == 'K2_High'),]
kidney_pheno2 <- kidney_pheno2 %>% arrange(condition)

#Filter down betas given the dataset and match order
kidney_betas1 <- kidney_betas[,match(kidney_pheno1$Basename, colnames(kidney_betas))]
colnames(kidney_betas1) <- kidney_pheno1$sample_id

kidney_betas2 <- kidney_betas[,match(kidney_pheno2$Basename, colnames(kidney_betas))]
colnames(kidney_betas2) <- kidney_pheno2$sample_id

kidney_betas1 <- kidney_betas1[rownames(kidney_betas1) %in% probes1,]
kidney_betas2 <- kidney_betas2[rownames(kidney_betas2) %in% probes2,]

pdf('/Users/adrianharris/Desktop/kidney_heatmaps.pdf')
heatmap.2(x=as.matrix(kidney_betas1), 
          Colv=FALSE, 
          dendrogram="none",
          scale="row",
          col="bluered",
          trace="none",
          labRow=TRUE,
          main="eGFR_12month-K1_Low-K1_High",
          ylab="Probes",
          xlab="Samples")

heatmap.2(x=as.matrix(kidney_betas2), 
          Colv=FALSE, 
          dendrogram="none",
          scale="row",
          col="bluered",
          trace="none",
          labRow=TRUE,
          main="eGFR_12month-K2_Low-K2_High",
          ylab="Probes",
          xlab="Samples")

dev.off()
