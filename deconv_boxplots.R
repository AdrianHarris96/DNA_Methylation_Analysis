#Boxplots for the deconvolution 
rm(list = ls())

dev.off()
library(pacman)
p_unload(all)
library(rio)
library(tidyverse)
library(ggplot2)
library(lumi)
library(matrixStats)
library(ggpubr)

#Phenotype for liver samples - Remove reference to local machine
liver_pheno <- import('/Users/adrianharris/Desktop/Mas_lab/epiDISH_liver_450k_output/liver_samples_450k_deconv.csv')

#Update sample_group column 
for (row in 1:nrow(liver_pheno)) {
  if (liver_pheno[row, 'sample_group'] == 'High Injury') {
    injury <- 'High'
  } else {
    injury <- 'Low'
  }
  liver_pheno[row, 'sample_group'] <- injury
}
rm(row, injury)

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

#Iterate through the different comparisons through a function
generate_boxplot <- function(pheno, condition1, condition2) {
  #Filter down according to conditon1
  if (condition1 == "DD_HI_L1") {
    pheno1 <- pheno[(pheno$donor_type == 'DD' & pheno$sample_group == 'High' & pheno$collection == 'L1'),]
  } else if (condition1 == "DD_HI_L2"){
    pheno1 <- pheno[(pheno$donor_type == 'DD' & pheno$sample_group == 'High' & pheno$collection == 'L2'),]
  } else if (condition1 == "DD_LI_L1") {
    pheno1 <- pheno[(pheno$donor_type == 'DD' & pheno$sample_group == 'Low' & pheno$collection == 'L1'),]
  } else if (condition1 == "DD_LI_L2") {
    pheno1 <- pheno[(pheno$donor_type == 'DD' & pheno$sample_group == 'Low' & pheno$collection == 'L2'),]
  } else if (condition1 == "LD_LI_L1") {
    pheno1 <- pheno[(pheno$donor_type == 'LD' & pheno$sample_group == 'Low' & pheno$collection == 'L1'),]
  } else if (condition1 == "LD_LI_L2") {
    pheno1 <- pheno[(pheno$donor_type == 'LD' & pheno$sample_group == 'Low' & pheno$collection == 'L2'),]
  } else {
    cat('Comparison does not exist\n')
  }
  
  #Filter down according to conditon2
  if (condition2 == "DD_HI_L1") {
    pheno2 <- pheno[(pheno$donor_type == 'DD' & pheno$sample_group == 'High' & pheno$collection == 'L1'),]
  } else if (condition2 == "DD_HI_L2"){
    pheno2 <- pheno[(pheno$donor_type == 'DD' & pheno$sample_group == 'High' & pheno$collection == 'L2'),]
  } else if (condition2 == "DD_LI_L1") {
    pheno2 <- pheno[(pheno$donor_type == 'DD' & pheno$sample_group == 'Low' & pheno$collection == 'L1'),]
  } else if (condition2 == "DD_LI_L2") {
    pheno2 <- pheno[(pheno$donor_type == 'DD' & pheno$sample_group == 'Low' & pheno$collection == 'L2'),]
  } else if (condition2 == "LD_LI_L1") {
    pheno2 <- pheno[(pheno$donor_type == 'LD' & pheno$sample_group == 'Low' & pheno$collection == 'L1'),]
  } else if (condition2 == "LD_LI_L2") {
    pheno2 <- pheno[(pheno$donor_type == 'LD' & pheno$sample_group == 'Low' & pheno$collection == 'L2'),]
  } else {
    cat('Comparison does not exist\n')
  }
  
  #Joining pheno1 and pheno2 
  pheno <- rbind(pheno1, pheno2)
  
  pheno$condition <- as.factor(pheno$condition)
  
  #Assumption of equal variance is not exactly met for some comparisons 
  #Consider organizing this portion of the code to avoid redundancy similar to what is done in the RNA section
  t1 <- t.test(pheno$B~pheno$condition)
  p1 <- round(t1$p.value, 6)
  t1 <- round(as.numeric(t1[[1]]), 3)
  t1_title <- paste("t-statistic: ", t1, sep = "")
  p1_title <- paste("p-value: ", p1, sep = "")
  
  t2 <- t.test(pheno$NK~pheno$condition)
  p2 <- round(t2$p.value, 6)
  t2 <- round(as.numeric(t2[[1]]), 3)
  t2_title <- paste("t-statistic: ", t2, sep = "")
  p2_title <- paste("p-value: ", p2, sep = "")
  
  t3 <- t.test(pheno$CD4T~pheno$condition)
  p3 <- round(t3$p.value, 6)
  t3 <- round(as.numeric(t3[[1]]), 3)
  t3_title <- paste("t-statistic: ", t3, sep = "")
  p3_title <- paste("p-value: ", p3, sep = "")
  
  t4 <- t.test(pheno$CD8T~pheno$condition)
  p4 <- round(t4$p.value, 6)
  t4 <- round(as.numeric(t4[[1]]), 3)
  t4_title <- paste("t-statistic: ", t4, sep = "")
  p4_title <- paste("p-value: ", p4, sep = "")
  
  t5 <- t.test(pheno$Mono~pheno$condition)
  p5 <- round(t5$p.value, 6)
  t5 <- round(as.numeric(t5[[1]]), 3)
  t5_title <- paste("t-statistic: ", t5, sep = "")
  p5_title <- paste("p-value: ", p5, sep = "")
  
  t6 <- t.test(pheno$Neutro~pheno$condition)
  p6 <- round(t6$p.value, 6)
  t6 <- round(as.numeric(t6[[1]]), 3)
  t6_title <- paste("t-statistic: ", t6, sep = "")
  p6_title <- paste("p-value: ", p6, sep = "")
  
  t7 <- t.test(pheno$Eosino~pheno$condition)
  p7 <- round(t7$p.value, 6)
  t7 <- round(as.numeric(t7[[1]]), 3)
  t7_title <- paste("t-statistic: ", t7, sep = "")
  p7_title <- paste("p-value: ", p7, sep = "")
  
  #Write to a line in the CSV
  file_suffix <- paste(condition1, condition2, sep="_")
  #return(c(file_suffix, t1, p1, t2, p2, t3, p3, t4, p4, t5, p5, t6, p6, t7, p7)) #Comment out this return to generate the boxplots with code below
  
  #(B, NK, CD4T, CD8T, Mono, Neutro, Eosino)
  p1 <-  ggplot(pheno, aes(x=condition, y=B, fill=condition, alpha = 0.6)) + geom_violin(trim=FALSE) + theme_bw() + geom_boxplot(width=.1) + labs(x='condition', y= "B") + theme(legend.position="none") + ggtitle(paste(t1_title, p1_title, sep = "\n"))
  p2 <-  ggplot(pheno, aes(x=condition, y=NK, fill=condition, alpha = 0.6)) + geom_violin(trim=FALSE) + theme_bw() + geom_boxplot(width=.1) + labs(x='condition', y= "NK") + theme(legend.position="none") + ggtitle(paste(t2_title, p2_title, sep = "\n"))
  p3 <-  ggplot(pheno, aes(x=condition, y=CD4T, fill=condition, alpha = 0.6)) + geom_violin(trim=FALSE) + theme_bw() + geom_boxplot(width=.1) + labs(x='condition', y= "CD4T") + theme(legend.position="none") + ggtitle(paste(t3_title, p3_title, sep = "\n"))
  p4 <-  ggplot(pheno, aes(x=condition, y=CD8T, fill=condition, alpha = 0.6)) + geom_violin(trim=FALSE) + theme_bw() + geom_boxplot(width=.1) + labs(x='condition', y= "CD8T") + theme(legend.position="none") + ggtitle(paste(t4_title, p4_title, sep = "\n"))
  p5 <-  ggplot(pheno, aes(x=condition, y=Mono, fill=condition, alpha = 0.6)) + geom_violin(trim=FALSE) + theme_bw() + geom_boxplot(width=.1) + labs(x='condition', y= "Mono") + theme(legend.position="none") + ggtitle(paste(t5_title, p5_title, sep = "\n"))
  p6 <-  ggplot(pheno, aes(x=condition, y=Neutro, fill=condition, alpha = 0.6)) + geom_violin(trim=FALSE) + theme_bw() + geom_boxplot(width=.1) + labs(x='condition', y= "Neutro") + theme(legend.position="none") + ggtitle(paste(t6_title, p6_title, sep = "\n"))
  p7 <-  ggplot(pheno, aes(x=condition, y=Eosino, fill=condition, alpha = 0.6)) + geom_violin(trim=FALSE) + theme_bw() + geom_boxplot(width=.1) + labs(x='condition', y= "Eosino") + theme(legend.position="none") + ggtitle(paste(t7_title, p7_title, sep = "\n"))
  
  g <- ggarrange(p1, p2, p3, p4, p5, p6, p7, ncol = 4, nrow=2)
  file_suffix <- paste(file_suffix, "_deconv_boxplot.jpeg", sep="")
  ggsave(file = paste("~/Desktop/liver_boxplots/", file_suffix, sep="") , units = c("in"), width=10, height=8, dpi=300, g)
  return('Saved')
}

#Write to CSV 
cols <- c('comparison', 'B_t_val', 'B_p_val', 'NK_t_val', 'NK_p_val', 'CD4T_t_val', 'CD4T_p_val', 'CD8T_t_val', 'CD8T_p_val', 'Mono_t_val', 'Mono_p_val', 'Neutro_t_val', 'Neutro_p_val', 'Eosino_t_val', 'Eosino_p_val')
t_df <- data.frame(matrix(ncol=length(cols), nrow=0))
colnames(t_df) <- cols

#Here, it has been edited to just generate the CSV - REMOVE assignment and return statement to plot
t_df[nrow(t_df)+1,] <- generate_boxplot(liver_pheno, "DD_HI_L1", "DD_HI_L2")
t_df[nrow(t_df)+1,] <- generate_boxplot(liver_pheno, "DD_HI_L1", "DD_LI_L1")
t_df[nrow(t_df)+1,] <- generate_boxplot(liver_pheno, "DD_HI_L1", "LD_LI_L1")
t_df[nrow(t_df)+1,] <- generate_boxplot(liver_pheno, "DD_HI_L2", "DD_LI_L2")
t_df[nrow(t_df)+1,] <- generate_boxplot(liver_pheno, "DD_HI_L2", "LD_LI_L2")
t_df[nrow(t_df)+1,] <- generate_boxplot(liver_pheno, "DD_LI_L1", "DD_LI_L2")
t_df[nrow(t_df)+1,] <- generate_boxplot(liver_pheno, "DD_LI_L1", "LD_LI_L1")
t_df[nrow(t_df)+1,] <- generate_boxplot(liver_pheno, "DD_LI_L2", "LD_LI_L2")
t_df[nrow(t_df)+1,] <- generate_boxplot(liver_pheno, "LD_LI_L1", "LD_LI_L2")

write.csv(t_df, file='~/Desktop/450k_deconv.csv', row.names = FALSE)

#Deconvolution violin plots for RNA data
library(readxl)
rm(list = ls())

deconv_df <- read_excel("~/Desktop/Mas_lab/liver_450k_immunecells.xlsx", sheet = 1)

#Rename first column to sample_name
deconv_df <- deconv_df %>% rename("sample_name" = "...1")

#Phenotype for liver samples
liver_pheno <- import('/Users/adrianharris/Documents/dna_methylation_analysis/liver_sample_sheet.csv')

#String parsing function
substrLeft <- function(x, n) {
  return(substr(x, 1, n))
}

#Fix the sample_name column for deconv_df
for (row in 1:nrow(deconv_df)) {
  old_name <- as.character(deconv_df[row, "sample_name"])
  if (grepl("pre", old_name)) {
    new_name <- substrLeft(old_name, 4)
    new_name <- paste(new_name, "L1", sep="_")
  } else {
    new_name <- substrLeft(old_name, 4)
    new_name <- paste(new_name, "L2", sep="_")
  }
  deconv_df[row, "sample_name"] <- new_name
  rm(old_name, new_name, row)
}

#Fix the sample_name for liver_pheno
for (row in 1:nrow(liver_pheno)) {
  name <- as.character(liver_pheno[row, "sample_name"])
  timepoint <- as.character(liver_pheno[row, "collection"])
  name <- substrLeft(name, 4)
  name <- paste(name, timepoint, sep="_")
  liver_pheno[row, "sample_name"] <- name
  rm(row, name, timepoint)
}

#Merge the two dataframes 
deconv_df <- merge(deconv_df, liver_pheno, by = "sample_name")

#Remove sample V037 
deconv_df <- deconv_df[(deconv_df$sample_name != "V037_L1"),]

#Update sample_group column 
for (row in 1:nrow(deconv_df)) {
  if (deconv_df[row, 'sample_group'] == 'High Injury') {
    injury <- 'High'
  } else {
    injury <- 'Low'
  }
  deconv_df[row, 'sample_group'] <- injury
}
rm(row, injury)

#Addition and filling of the condition column 
deconv_df$condition <- 'NA'

#Filling the condition column
for (row in 1:nrow(deconv_df)) {
  if (deconv_df[row, 'donor_type'] == 'DD' & deconv_df[row, 'sample_group'] == 'High' &  deconv_df[row, 'collection'] == 'L1') {
    deconv_df[row, 'condition'] <- "DD_HI_L1"
  } else if (deconv_df[row, 'donor_type'] == 'DD' & deconv_df[row, 'sample_group'] == 'High' &  deconv_df[row, 'collection'] == 'L2') {
    deconv_df[row, 'condition'] <- "DD_HI_L2"
  } else if (deconv_df[row, 'donor_type'] == 'DD' & deconv_df[row, 'sample_group'] == 'Low' &  deconv_df[row, 'collection'] == 'L1') {
    deconv_df[row, 'condition'] <- "DD_LI_L1"
  } else if (deconv_df[row, 'donor_type'] == 'DD' & deconv_df[row, 'sample_group'] == 'Low' &  deconv_df[row, 'collection'] == 'L2') {
    deconv_df[row, 'condition'] <- "DD_LI_L2"
  } else if (deconv_df[row, 'donor_type'] == 'LD' & deconv_df[row, 'sample_group'] == 'Low' &  deconv_df[row, 'collection'] == 'L1') {
    deconv_df[row, 'condition'] <- "LD_LI_L1"
  } else if (deconv_df[row, 'donor_type'] == 'LD' & deconv_df[row, 'sample_group'] == 'Low' &  deconv_df[row, 'collection'] == 'L2') {
    deconv_df[row, 'condition'] <- "LD_LI_L2"
  }
}
rm(row)

#Plot violin plots - called in the second function
plot_violin <- function(phenotype, cells) {
  y1 <- phenotype[,cells]
  x1 <- phenotype[,"condition"]
  t <- t.test(y1~x1)
  t_val <- round(as.numeric(t[[1]]), 3)
  p_val <- round(t$p.value, 6)
  t_title <- paste("t-statistic: ", t_val, sep = "")
  p_title <- paste("p-value: ", p_val, sep = "")
  plot1 <- ggplot(phenotype, aes(x=condition, y=phenotype[,cells], fill=condition, alpha = 0.6)) + geom_violin(trim=FALSE) + theme_bw() + geom_boxplot(width=0.1) + labs(x='condition', y= as.character(cells)) + theme(legend.position="none") + ggtitle(paste(t_title, p_title, sep = "\n"))
  return(plot1)
}

#Iterate through the different comparisons through a function
generate_boxplot <- function(pheno, condition1, condition2) {
  #Filter down according to conditon1
  if (condition1 == "DD_HI_L1") {
    pheno1 <- pheno[(pheno$donor_type == 'DD' & pheno$sample_group == 'High' & pheno$collection == 'L1'),]
  } else if (condition1 == "DD_HI_L2"){
    pheno1 <- pheno[(pheno$donor_type == 'DD' & pheno$sample_group == 'High' & pheno$collection == 'L2'),]
  } else if (condition1 == "DD_LI_L1") {
    pheno1 <- pheno[(pheno$donor_type == 'DD' & pheno$sample_group == 'Low' & pheno$collection == 'L1'),]
  } else if (condition1 == "DD_LI_L2") {
    pheno1 <- pheno[(pheno$donor_type == 'DD' & pheno$sample_group == 'Low' & pheno$collection == 'L2'),]
  } else if (condition1 == "LD_LI_L1") {
    pheno1 <- pheno[(pheno$donor_type == 'LD' & pheno$sample_group == 'Low' & pheno$collection == 'L1'),]
  } else if (condition1 == "LD_LI_L2") {
    pheno1 <- pheno[(pheno$donor_type == 'LD' & pheno$sample_group == 'Low' & pheno$collection == 'L2'),]
  } else {
    cat('Comparison does not exist\n')
  }
  
  #Filter down according to conditon2
  if (condition2 == "DD_HI_L1") {
    pheno2 <- pheno[(pheno$donor_type == 'DD' & pheno$sample_group == 'High' & pheno$collection == 'L1'),]
  } else if (condition2 == "DD_HI_L2"){
    pheno2 <- pheno[(pheno$donor_type == 'DD' & pheno$sample_group == 'High' & pheno$collection == 'L2'),]
  } else if (condition2 == "DD_LI_L1") {
    pheno2 <- pheno[(pheno$donor_type == 'DD' & pheno$sample_group == 'Low' & pheno$collection == 'L1'),]
  } else if (condition2 == "DD_LI_L2") {
    pheno2 <- pheno[(pheno$donor_type == 'DD' & pheno$sample_group == 'Low' & pheno$collection == 'L2'),]
  } else if (condition2 == "LD_LI_L1") {
    pheno2 <- pheno[(pheno$donor_type == 'LD' & pheno$sample_group == 'Low' & pheno$collection == 'L1'),]
  } else if (condition2 == "LD_LI_L2") {
    pheno2 <- pheno[(pheno$donor_type == 'LD' & pheno$sample_group == 'Low' & pheno$collection == 'L2'),]
  } else {
    cat('Comparison does not exist\n')
  }
  
  #Joining pheno1 and pheno2 
  pheno <- rbind(pheno1, pheno2)
  
  pheno$condition <- as.factor(pheno$condition)
  
  #Assumption of equal variance is not exactly met for some comparisons 
  cell_types <- c("B.cells.naive", 
    "B.cells.memory", 
    "Plasma.cells", 
    "T.cells.CD8", 
    "T.cells.CD4.naive",
    "T.cells.CD4.memory.resting", 
    "T.cells.CD4.memory.activated", 
    "T.cells.follicular.helper", 
    "T.cells.regulatory..Tregs.", 
    "T.cells.gamma.delta", 
    "NK.cells.resting", 
    "NK.cells.activated", 
    "Monocytes", 
    "Macrophages.M0", 
    "Macrophages.M1", 
    "Macrophages.M2", 
    "Dendritic.cells.resting", 
    "Dendritic.cells.activated", 
    "Mast.cells.resting", 
    "Mast.cells.activated", 
    "Eosinophils", 
    "Neutrophils")
  
  #Generate plots for each column
  p1 <- plot_violin(pheno, cell_types[1])
  p2 <- plot_violin(pheno, cell_types[2])
  p3 <- plot_violin(pheno, cell_types[3])
  p4 <- plot_violin(pheno, cell_types[4])
  p5 <- plot_violin(pheno, cell_types[5])
  p6 <- plot_violin(pheno, cell_types[6])
  p7 <- plot_violin(pheno, cell_types[7])
  p8 <- plot_violin(pheno, cell_types[8])
  p9 <- plot_violin(pheno, cell_types[9])
  p10 <- plot_violin(pheno, cell_types[10])
  p11 <- plot_violin(pheno, cell_types[11])
  p12 <- plot_violin(pheno, cell_types[12])
  p13 <- plot_violin(pheno, cell_types[13])
  p14 <- plot_violin(pheno, cell_types[14])
  p15 <- plot_violin(pheno, cell_types[15])
  p16 <- plot_violin(pheno, cell_types[16])
  p17 <- plot_violin(pheno, cell_types[17])
  p18 <- plot_violin(pheno, cell_types[18])
  p19 <- plot_violin(pheno, cell_types[19])
  p20 <- plot_violin(pheno, cell_types[20])
  p21 <- plot_violin(pheno, cell_types[21])
  p22 <- plot_violin(pheno, cell_types[22])
  
  g <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, ncol = 4, nrow=6)
  file_suffix <- paste(condition1, condition2, sep="_")
  file_suffix <- paste(file_suffix, "_deconv_boxplot.jpeg", sep="")
  ggsave(file = paste("~/Desktop/liver_rna_boxplots/", file_suffix, sep="") , units = c("in"), width=10, height=14, dpi=300, g)
  return('Saved')
}

generate_boxplot(deconv_df, "DD_HI_L1", "DD_HI_L2")
generate_boxplot(deconv_df, "DD_HI_L1", "DD_LI_L1")
generate_boxplot(deconv_df, "DD_HI_L1", "LD_LI_L1")
generate_boxplot(deconv_df, "DD_HI_L2", "DD_LI_L2")
generate_boxplot(deconv_df, "DD_HI_L2", "LD_LI_L2")
generate_boxplot(deconv_df, "DD_LI_L1", "DD_LI_L2")
generate_boxplot(deconv_df, "DD_LI_L1", "LD_LI_L1")
generate_boxplot(deconv_df, "DD_LI_L2", "LD_LI_L2")
generate_boxplot(deconv_df, "LD_LI_L1", "LD_LI_L2")


