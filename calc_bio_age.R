#Calculation of biological age 
rm(list = ls())
dev.off()

library(methylclock)
library(rio)
library(tidyverse)
library(ggplot2)

#Kidney section for calculating biological age 
betas_file <- '/Users/adrianharris/Desktop/Mas_lab/week_08262022/kidney_combined_beta_dir/beta_values.csv'
metadata_file <- '/Users/adrianharris/Documents/GitHub/DNA_Methylation_Analysis/paired_kidney_sample_sheet.csv'
betas <- import(betas_file)

#Renaming a column and calculating age 
betas <- betas %>% rename(ProbeID = V1)
age <- DNAmAge(betas)

pheno_df <- import(metadata_file)

#Formatting the age dataframe before merging with the metadata
age <- age %>% rename(Basename = id)
pheno_df <- merge(pheno_df, age, by = 'Basename')

x1 <- pheno_df$donor_age
y1 <- pheno_df$Horvath

#Plotting via points 
plot <- ggplot(data=pheno_df, mapping = aes(x = x1, y = y1, color=eGFR_1month)) + theme_bw() + geom_point(alpha=0.9) + labs(x='chronological age', y='biological age (DNAm)') + xlim(20, 80) + ylim(20, 80) + ggtitle('kidney - biological v. chronological age')
plot + geom_abline(slope = 1, intercept=0)

#Plotting with labels 
plot1 <- ggplot(data=pheno_df, mapping = aes(x = x1, y = y1, color=eGFR_1month)) + theme_bw() + geom_text(aes(label = sample_id), size = 2.5) + labs(x='chronological age', y='biological age (DNAm)') + xlim(20, 80) + ylim(20, 80) + ggtitle('kidney - biological v. chronological age')
plot1 + geom_abline(slope = 1, intercept=0)

#Consistency with sample_name
for (row in 1:nrow(pheno_df)){
  sample <- pheno_df[row, 'sample_id']
  sample <- strsplit(sample, "L")
  sample <- unlist(sample)[1]
  sample <- str_trim(sample, side="both")
  collection <- pheno_df[row, 'collection']
  new_name <- paste(sample, collection, sep = " ")
  pheno_df[row, 'sample_name'] <- new_name
}

#SPlitting the pheno_df
pheno_df1 <- pheno_df[(pheno_df$collection == 'L1'),]
pheno_df2 <- pheno_df[(pheno_df$collection == 'L2'),]

pheno_df1 <- pheno_df1[!(pheno_df1$sample_name == 'V037 L1'),]

#Preparing for the merge by creating pair column
for (row in 1:nrow(pheno_df1)) {
  common <- pheno_df1[row, 'sample_name']
  common <- substr(common,1,nchar(common)-1)
  pheno_df1[row, 'pair'] <- common
}

for (row in 1:nrow(pheno_df2)) {
  common <- pheno_df2[row, 'sample_name']
  common <- substr(common,1,nchar(common)-1)
  pheno_df2[row, 'pair'] <- common
}

#Merging on pair
final_df <- merge(pheno_df1, pheno_df2, by = 'pair')

final_df$diff <- 'na'
for (row in 1:nrow(final_df)) {
  diff <- final_df[row, 'Horvath.y'] - final_df[row, 'Horvath.x']
  final_df[row, 'diff'] <- diff
}

L1 <- final_df$Horvath.x
L2 <- final_df$Horvath.y
group <- final_df$sample_group.x
d <- data.frame(L1 = L1, L2 = L2, group = group)

for (row in 1:nrow(d)){
  if (d[row, 'group'] == "High Injury") {
    d[row, 'group'] <- "red"
  } else {
    d[row, 'group'] <- "blue"
  }
}

x <- factor(d$group)
y <- append(x, x)
library(ggpubr)
ggpaired(d, cond1 = "L1", cond2 = "L2", width = 0.1, point.size = 1.0, line.size = 0.3, line.color = y, ggtheme=theme_bw(), title="Changes in biological age (DNAm) - Injury Status") + labs(x='collection', y='biological age (DNAm)', color="Legend")
line.color = "black"

#Rename columns and write to csv
d <- d %>% rename('L1 Biological Age' = 'L1')
d <- d %>% rename('L2 Biological Age' = 'L2')
write.csv(d, file = '/Users/adrianharris/Desktop/liver_bio_age_changes.csv', row.names = FALSE)

#Liver 
betas <- import('/Users/adrianharris/Desktop/large_files_liver_combined/beta_values.csv')
betas <- betas %>% rename(ProbeID = V1)
age <- DNAmAge(betas)

pheno_df <- import('/Users/adrianharris/Documents/dna_methylation_analysis/liver_sample_sheet.csv')

age <- age %>% rename(Basename = id)

pheno_df <- merge(pheno_df, age, by = 'Basename')
pheno_df <- subset(pheno_df, sample_group != 'control')
pheno_df <- subset(pheno_df, sample_group != 'Control')

# pheno_df <- subset(pheno_df, collection == 'L1')

x1 <- pheno_df$donor_age
y1 <- pheno_df$Horvath

plot <- ggplot(data=pheno_df, mapping = aes(x = x1, y = y1, color=sample_group)) + theme_bw() + geom_point(alpha=0.9) + labs(x='chronological age', y='biological age (DNAm)') + xlim(20, 120) + ylim(20, 120) + ggtitle('liver - biological v. chronological age')
plot <- ggplot(data=pheno_df, mapping = aes(x = x1, y = y1, color=sample_group)) + theme_bw() + geom_text(aes(label = sample_name), size = 2.5) + labs(x='chronological age', y='biological age (DNAm)') + xlim(20, 120) + ylim(20, 120) + ggtitle('liver - biological v. chronological age')
plot + geom_abline(slope = 1, intercept=0)
print(plot)


#Pairs - liver 
#Consistency with sample_name
for (row in 1:nrow(pheno_df)){
  sample <- pheno_df[row, 'sample_name']
  sample <- strsplit(sample, "L")
  sample <- unlist(sample)[1]
  sample <- str_trim(sample, side="both")
  collection <- pheno_df[row, 'collection']
  new_name <- paste(sample, collection, sep = " ")
  pheno_df[row, 'sample_name'] <- new_name
}

#SPlitting the pheno_df
pheno_df1 <- pheno_df[(pheno_df$collection == 'L1'),]
pheno_df2 <- pheno_df[(pheno_df$collection == 'L2'),]

pheno_df1 <- pheno_df1[!(pheno_df1$sample_name == 'V037 L1'),]

#Preparing for the merge by creating pair column
for (row in 1:nrow(pheno_df1)) {
  common <- pheno_df1[row, 'sample_name']
  common <- substr(common,1,nchar(common)-1)
  pheno_df1[row, 'pair'] <- common
}

for (row in 1:nrow(pheno_df2)) {
  common <- pheno_df2[row, 'sample_name']
  common <- substr(common,1,nchar(common)-1)
  pheno_df2[row, 'pair'] <- common
}

#Merging on pair
final_df <- merge(pheno_df1, pheno_df2, by = 'pair')

final_df$diff <- 'na'
for (row in 1:nrow(final_df)) {
  diff <- final_df[row, 'Horvath.y'] - final_df[row, 'Horvath.x']
  final_df[row, 'diff'] <- diff
}

L1 <- final_df$Horvath.x
L2 <- final_df$Horvath.y
group <- final_df$sample_group.x
d <- data.frame(L1 = L1, L2 = L2, group = group)

for (row in 1:nrow(d)){
  if (d[row, 'group'] == "High Injury") {
    d[row, 'group'] <- "red"
  } else {
    d[row, 'group'] <- "blue"
  }
}

x <- factor(d$group)
y <- append(x, x)
library(ggpubr)
ggpaired(d, cond1 = "L1", cond2 = "L2", width = 0.1, point.size = 1.0, line.size = 0.3, line.color = y, ggtheme=theme_bw(), title="Changes in biological age (DNAm) - Injury Status") + labs(x='collection', y='biological age (DNAm)', color="Legend")
line.color = "black"

#Rename columns and write to csv
d <- d %>% rename('L1 Biological Age' = 'L1')
d <- d %>% rename('L2 Biological Age' = 'L2')
write.csv(d, file = '/Users/adrianharris/Desktop/liver_bio_age_changes.csv', row.names = FALSE)


