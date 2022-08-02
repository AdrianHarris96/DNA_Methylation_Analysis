#Calculation of biological age 
rm(list = ls())
dev.off()
#BiocManager::install("methylclock")
library(methylclock)
library(rio)
library(tidyverse)
library(ggplot2)

#MethylationData <- get_MethylationDataExample()
#Kidney 
betas <- import('/Users/adrianharris/Desktop/kidney_combined_beta_dir/beta_values.csv')
betas <- betas %>% rename(ProbeID = V1)
age <- DNAmAge(betas)

pheno_df <- import('/Users/adrianharris/Documents/dna_methylation_analysis/paired_kidney_sample_sheet.csv')

age <- age %>% rename(Basename = id)

pheno_df <- merge(pheno_df, age, by = 'Basename')
pheno_df <- subset(pheno_df, time == 'K1')

x1 <- pheno_df$donor_age
y1 <- pheno_df$Horvath

plot <- ggplot(data=pheno_df, mapping = aes(x = x1, y = y1, color=eGFR_1month)) + theme_bw() + geom_point(alpha=0.9) + labs(x='chronological age', y='biological age (DNAm)') + xlim(20, 80) + ylim(20, 80) + ggtitle('kidney - biological v. chronological age')
plot <- ggplot(data=pheno_df, mapping = aes(x = x1, y = y1, color=eGFR_1month)) + theme_bw() + geom_text(aes(label = sample_id), size = 2.5) + labs(x='chronological age', y='biological age (DNAm)') + xlim(20, 80) + ylim(20, 80) + ggtitle('kidney - biological v. chronological age')
plot + geom_abline(slope = 1, intercept=0)
print(plot)

#Liver 
betas <- import('/Users/adrianharris/Desktop/large_files_liver_combined/beta_values.csv')
betas <- betas %>% rename(ProbeID = V1)
age <- DNAmAge(betas)

pheno_df <- import('/Users/adrianharris/Documents/dna_methylation_analysis/liver_sample_sheet.csv')

age <- age %>% rename(Basename = id)

pheno_df <- merge(pheno_df, age, by = 'Basename')
pheno_df <- subset(pheno_df, sample_group != 'control')
pheno_df <- subset(pheno_df, sample_group != 'Control')
pheno_df <- subset(pheno_df, collection == 'L1')

x1 <- pheno_df$donor_age
y1 <- pheno_df$Horvath

plot <- ggplot(data=pheno_df, mapping = aes(x = x1, y = y1, color=sample_group)) + theme_bw() + geom_point(alpha=0.9) + labs(x='chronological age', y='biological age (DNAm)') + xlim(20, 120) + ylim(20, 120) + ggtitle('liver - biological v. chronological age')
plot <- ggplot(data=pheno_df, mapping = aes(x = x1, y = y1, color=sample_group)) + theme_bw() + geom_text(aes(label = sample_name), size = 2.5) + labs(x='chronological age', y='biological age (DNAm)') + xlim(20, 120) + ylim(20, 120) + ggtitle('liver - biological v. chronological age')
plot + geom_abline(slope = 1, intercept=0)
print(plot)
