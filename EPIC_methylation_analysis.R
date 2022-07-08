#Script for DNA Methylation Analysis
rm(list = ls())
library(minfi)
library(rio)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(lumi)
library(limma)

#Example input: Rscript methylation_analysis.R <pheno_file> <base_dir> <output_dir>
args=commandArgs(trailingOnly=TRUE)
pheno_file = args[1]
base_dir = args[2]
git_dir = args[3]
output_dir = args[4]

#Local Machine
# pheno_file = '/Users/adrianharris/Documents/dna_methylation_analysis/paired_kidney_sample_sheet.csv'
# base_dir = '/Users/adrianharris/Desktop/kidney/'
# git_dir = '/Users/adrianharris/Documents/dna_methylation_analysis/'
# output_dir = '/Users/adrianharris/Desktop/epic_kidney/'

if (file.exists(output_dir)) {
  cat("Directory already exists\n")
} else {
  dir.create(output_dir)
}

#Importing manually-curated sample sheet
pheno_df <- import(pheno_file)

#Must remove outlier sample, 203504430032_R01C01 (and its paired sample 203504430032-R02C01)
pheno_df <- pheno_df[!(pheno_df$Basename == '203504430032_R01C01' | pheno_df$Basename == '203504430032_R02C01'),]
pheno_df <- pheno_df[!(pheno_df$sample_id == 'KUT4_K2' | pheno_df$sample_id == 'KUT4_K1'),]

nrow(subset(pheno_df, array_type == 'EPIC'))

#Number of files 
nrow(subset(pheno_df, time == 'K1'))
nrow(subset(pheno_df, time == 'K2'))

#Specify the directories and read in respective IDAT Files 
dirEPIC <- paste(base_dir, "EPIC_array", sep="")
if (file.exists(paste(output_dir, "rgSet.RDS", sep=""))) {
  cat('Loading in rgSet (combined)\n')
  rgSet <- readRDS(paste(output_dir, "rgSet.RDS", sep=""))
} else {
  cat('Generate rgSet\n')
  rgSet <- read.metharray.exp(base=dirEPIC, target=pheno_df, force=TRUE)
  #Saving set as an RDS file 
  saveRDS(rgSet, file = paste(output_dir, "rgSet.RDS", sep=""))
}

#Color Scheme Defined 
pal <- brewer.pal(4,"Dark2")

#Calculate Detection p-values
if (file.exists(paste(output_dir, "beta_values.csv", sep=""))) {
  cat('Skipping p-value step\n')
} else if (file.exists(paste(output_dir, "p-values.csv", sep=""))) {
  cat('P-values.csv is already exists\n')
} else {
  cat('Performing p-value detection\n')
  detP <- detectionP(rgSet)
  detP <- data.frame(detP)
  sample_names <- colnames(detP)
  sample_names <- substr(sample_names, 2, nchar(sample_names))
  detPmeans <- colMeans(detP)
  detP_df <- data.frame(matrix(ncol=2, nrow=ncol(detP)))
  detP_df['X1'] <- sample_names
  detP_df['X2'] <- detPmeans
  colnames(detP_df) <- c('Basename', 'p_value_mean')
  write.csv(detP_df, file = paste(output_dir, "p-values.csv", sep=""), row.names = FALSE)
  rm(sample_names, detPmeans)
  
  #Merging to pheno dataframe for color coding and plot via barplots 
  detP_df <- merge(detP_df, pheno_df, by = 'Basename')
  
  #Plotting EPIC barplot
  barplot((subset(detP_df, array_type == 'EPIC'))$p_value_mean, col=pal[factor(detP_df$time)], names.arg=(subset(detP_df, array_type == 'EPIC'))$sample_name, las=2, cex.names=0.4, cex.axis=0.5, space=0.5, ylab="Mean detection p-values", main='EPIC Array')
  legend("topleft", legend=levels(factor(detP_df$time)), fill=pal,
         cex=0.27, bty = "n", bg="white")
  
  rm(detP_df, detP)
  dev.off()
  par(mfrow=c(1,1))
}

if (file.exists(paste(output_dir, "preprocessQC.jpeg", sep=""))) {
  cat('Preprocessing already performed\n')
} else {
  cat('Performing preprocessing and plotting\n')
  #preprocessing QC and plotting
  mtSet <- preprocessRaw(rgSet)
  qc <- getQC(mtSet)
  #plotQC(qc)
  qc <- data.frame(qc)
  qc['Basename'] <- row.names(qc)
  qc <- merge(qc, pheno_df, by = 'Basename')
  jpeg(paste(output_dir, "preprocessQC.jpeg", sep=""), quality = 90)
  plot <- ggplot(data=qc, mapping = aes(x = mMed, y = uMed, color=time)) + geom_point(aes(shape=array_type), alpha=0.5) + xlim(7, 14) + ylim(7, 14) + theme_bw()+ geom_abline(slope=-1, intercept = 21.25 , color="black", linetype="dashed", size=0.5) + scale_color_manual(values=pal)
  print(plot)
  dev.off()
  jpeg(paste(output_dir, "preprocessDensity.jpeg", sep=""), quality = 90)
  densityPlot(mtSet, sampGroups =qc$array_type)
  dev.off()
  rm(mtSet, qc, plot)
}

#Normalization
if (file.exists(paste(output_dir, "postNormQC.jpeg", sep=""))) {
  cat('Loading normalization\n')
  mtSet <- readRDS(paste(output_dir, "mtSet.RDS", sep=""))
} else {
  cat('Performing normalization and plotting\n')
  #Normalization and plotting 
  mtSet <- preprocessNoob(rgSet)
  saveRDS(mtSet, file = paste(output_dir, "mtSet.RDS", sep=""))
postqc <- getQC(mtSet)
postqc <- data.frame(postqc)
postqc['Basename'] <- row.names(postqc)
postqc <- merge(postqc, pheno_df, by = 'Basename')
jpeg(paste(output_dir, "postNormQC.jpeg", sep=""), quality = 90)
plot2 <- ggplot(data=postqc, mapping = aes(x = mMed, y = uMed, color=time)) + geom_point(aes(shape=array_type), alpha=0.5) + xlim(5, 14) + ylim(5, 14) + theme_bw()+ geom_abline(slope=-1, intercept = 21.25 , color="black", linetype="dashed", size=0.5) + scale_color_manual(values=pal)
print(plot2)
dev.off()
jpeg(paste(output_dir, "postNormDensity.jpeg", sep=""), quality = 90)
densityPlot(mtSet, sampGroups = postqc$array_type)
dev.off()
rm(postqc, plot2)
}

# Map to Genome
if (file.exists(paste(output_dir, "beta_values.csv", sep=""))) {
  cat('Skip genomic methyl set\n')
}else {
  cat('Converting to Genomic Methyl Set\n')
  rSet <- ratioConvert(mtSet, what = "both", keepCN = TRUE)
  gmtSet <- mapToGenome(rSet)
  print(dim(gmtSet)) #Number of probes = 452453
  
  # #Predicted Sex 
  # predictedSex <- getSex(gmtSet, cutoff = -2)
  # gmtSet <- addSex(gmtSet, sex = predictedSex)
  # plotSex(gmtSet, id = row.names(predictedSex))
  # rm(predictedSex)
  
  #Removing probes that include SNPs
  snps <- getSnpInfo(gmtSet)
  gmtSet <- addSnpInfo(gmtSet)
  gr <- granges(gmtSet)
  gmtSet <- dropLociWithSnps(gmtSet, snps=c("SBE","CpG"), maf=0)
  print(dim(gmtSet)) #Number of probes = 436144
  rm(snps, gr)
  
  # Filter unwanted sites 
  # ensure probes are in the same order in the mSetSq and detP objects
  #rgSet <- readRDS(paste(output_dir, "rgSet.RDS", sep=""))
  detP <- detectionP(rgSet)
  detP <- detP[match(featureNames(gmtSet),rownames(detP)),] 
  
  # remove any probes that have failed in >50% of samples
  keep <- detP < 0.01
  keep.probes <- rownames(detP[rowMeans(keep)>=0.5,]) #probes that failed detection in more than half of the samples
  gmtSet <- gmtSet[keep.probes,] 
  print(dim(gmtSet)) #Number of probes = 436128
  rm(keep.probes)
  
  #Remove probes that located on the X or Y chromosome
  annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  
  keep <- !(featureNames(gmtSet) %in% annEPIC$Name[annEPIC$chr %in% c("chrX","chrY")]) #remove probes that are not of the chrom x or y
  table(keep)
  gmtSet <- gmtSet[keep,]
  print(dim(gmtSet)) #Number of probes = 425718
  rm(annEPIC, keep)
  
  #Creation of bad probes character and filter gmtSet
  cross.react <- read.csv(paste(git_dir, '48639-non-specific-probes-Illumina450k.csv', sep=""), head = T, as.is = T)
  cross.react.probes <- as.character(cross.react$TargetID)
  #Probes identified with potential hybridization issues
  multi.map <- read.csv(paste(git_dir, 'HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt', sep=""), head = F, as.is = T)
  multi.map.probes <- as.character(multi.map$V1)
  # determine unique probes between the cross-reactive and multi-map probes
  bad.probes <- unique(c(cross.react.probes, multi.map.probes))
  keep <- !(featureNames(gmtSet) %in% bad.probes) #Removal of these bad probes
  table(keep)
  gmtSet <- gmtSet[keep,]
  print(dim(gmtSet)) #Number of probes = 392871
  rm(cross.react, multi.map, bad.probes, cross.react.probes, multi.map.probes, keep)
  
  #Extract betas and m_values
  beta_values <- getBeta(gmtSet)
  # m_values <- getM(gmtSet)
  
  #Remove probes Hypomethylated in all samples identified by CpG sites having beta < 0.05 in all samples
  beta_values_filtered <- as.data.frame(beta_values) 
  cat('Hypomethylated\n')
  print(dim(filter_all(beta_values_filtered, all_vars(. < 0.05)))) #Number of hypomethylated
  beta_values_filtered <- filter_all(beta_values_filtered, any_vars(. >= 0.05)) 
  #33419 Hypomethylated
  
  #Remove probes Hypermethylated in all samples identified by CpG sites having beta > 0.95 in all samples
  cat("Hypermethylated\n")
  print(dim(filter_all(beta_values_filtered, all_vars(. > 0.95))))
  #3091 hypermethylated
  beta_values_filtered <- filter_all(beta_values_filtered, any_vars(. < 0.95)) 
  dim(beta_values_filtered)
  m_values = beta2m(beta_values_filtered)
}

clustering <- function(pheno, month, comparison, betas) {
  #Drop other columns and rename eGFR_month -> eGFR
  if (month == 'eGFR_1month') {
    pheno <- subset(pheno, select = -c(eGFR_12month, eGFR_24month))
    pheno <- pheno %>% rename('eGFR' = 'eGFR_1month')
  } else if (month == 'eGFR_12month') {
    pheno <- subset(pheno, select = -c(eGFR_1month, eGFR_24month))
    pheno <- pheno %>% rename('eGFR' = 'eGFR_12month')
  } else {
    pheno <- subset(pheno, select = -c(eGFR_1month, eGFR_12month))
    pheno <- pheno %>% rename('eGFR' = 'eGFR_24month')
  }
  
  if (comparison == 'K1_Low_K1_High') {
    pheno <- pheno[((pheno$time == 'K1' & pheno$eGFR == 'Low') | (pheno$time == 'K1' & pheno$eGFR == 'High')),]
  } else if (comparison == 'K2_Low_K2_High') {
    pheno <- pheno[((pheno$time == 'K2' & pheno$eGFR == 'Low') | (pheno$time == 'K2' & pheno$eGFR == 'High')),]
  } else if (comparison == 'K1_High_K2_High') {
    pheno <- pheno[((pheno$time == 'K1' & pheno$eGFR == 'High') | (pheno$time == 'K2' & pheno$eGFR == 'High')),]
  } else if (comparison == 'K1_Low_K2_Low') {
    pheno <- pheno[((pheno$time == 'K1' & pheno$eGFR == 'Low') | (pheno$time == 'K2' & pheno$eGFR == 'Low')),]
  } else if (comparison == 'K1_High_K2_Low') {
    pheno <- pheno[((pheno$time == 'K1' & pheno$eGFR == 'High') | (pheno$time == 'K2' & pheno$eGFR == 'Low')),]
  } else if (comparison == 'K1_Low_K2_High') {
    pheno <- pheno[((pheno$time == 'K1' & pheno$eGFR == 'Low') | (pheno$time == 'K2' & pheno$eGFR == 'High')),]
  } else {
    cat('Comparison request does not exist\n')
    q()
  }
  
  #Filter beta dataframe using column names for the relevant comparison
  beta_values_filtered <- beta_values_filtered[,(colnames(beta_values_filtered) %in% pheno_df$Basename)]
  # dim(beta_values_filtered)
  
  # #EXCLUDE DCD SAMPLES - Kidney dataset
  # beta_values_case <- beta_values_filtered[,!(colnames(beta_values_filtered) %in% c("9296930129_R05C01", "9305651174_R01C01", "9305651174_R03C01", "9305651174_R02C02", "9305651174_R03C02", "9305651191_R02C02", "9305651191_R04C02", "201465900002_R04C01", "202240580106_R03C01", "202240580208_R03C01", "202259340119_R05C01", "202259350016_R04C01", "203496240002_R03C01", "203504430032_R05C01", "204001300109_R07C01", "204001300109_R08C01", "204001350016_R01C01", "202702240079_R06C01"))]
  # #Removing rows based on the sample_name column in phenotype dataframe
  # pheno_df_case <- pheno_df[!(pheno_df$sample_name %in% c("9296930129_R05C01", "9305651174_R01C01", "9305651174_R03C01", "9305651174_R02C02", "9305651174_R03C02", "9305651191_R02C02", "9305651191_R04C02", "201465900002_R04C01", "202240580106_R03C01", "202240580208_R03C01", "202259340119_R05C01", "202259350016_R04C01", "203496240002_R03C01", "203504430032_R05C01", "204001300109_R07C01", "204001300109_R08C01", "204001350016_R01C01", "202702240079_R06C01")),]
  
  #betas <- betas[1:20000, ] #Used for troubleshooting
  #m_values <- beta2m(betas)
  
  cat("PCA plots - Betas\n")
  #Transpose resulting PCA plot - beta values 
  transposed_beta <- t(betas) #t(beta_values_case)
  transposed_beta <- data.frame(transposed_beta)
  pca_general <- prcomp(transposed_beta, center=TRUE)
  var_explained <- pca_general$sdev^2/sum(pca_general$sdev^2)
  scores = as.data.frame(pca_general$x)
  scores['Basename'] <- row.names(scores)
  scores <- merge(scores, pheno, by = 'Basename')
  output <- paste(month, "_", sep="")
  output <- paste(output, comparison, sep="")
  title <- paste(output, "_betas", sep="")
  output_path <- paste(output_dir, title, sep="")
  #jpeg(paste(output_path, "_PCA.jpeg", sep=""), quality = 100)
  plot <- ggplot(data=scores, mapping = aes(x = PC1, y = PC2, color=eGFR)) +theme_bw() + geom_point(aes(shape=time), alpha=0.5, size=2) + labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"), y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) + scale_color_manual(values=pal) + ggtitle(paste(title, "_PCA", sep="")) + geom_text(aes(label = sample_id), size=1.75, colour="black")
  #plot + geom_text(aes(label = sample_name), size=3.5) + xlim(-100, 400)
  print(plot)
  #dev.off()
  #Converting NAs to 0 in donor_age
  vec <- scores$donor_age
  vec[is.na(vec)] <- 0
  scores$donor_age <- vec
  
  if (length(unique(scores$donor_gender)) == 3) {
    plot2 <- ggplot(data=scores, mapping = aes(x = PC1, y = PC2, color=donor_gender)) + geom_point(size = 2, alpha = 0.5) + labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"), y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) + ggtitle(paste(title, " _PCA - donor gender and donor age", sep="")) + geom_text(aes(label=donor_age), size=1.75, colour="black") + scale_color_manual(values=c("grey", pal[1], pal[2])) + theme_bw()
  } else {
    plot2 <- ggplot(data=scores, mapping = aes(x = PC1, y = PC2, color=donor_gender)) + geom_point(size = 2, alpha = 0.5) + labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"), y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) + ggtitle(paste(title, " -PCA - donor gender and donor age", sep="")) + geom_text(aes(label=donor_age), size=1.75, colour="black") + scale_color_manual(values=c(pal[1], pal[2])) + theme_bw()
  }
  print(plot2)
  
  # cat("PCA plots - m_values\n")
  # #Transpose resulting PCA plot - m_values
  # transposed_m <- t(m_values) 
  # transposed_m <- data.frame(transposed_m)
  # pca_general <- prcomp(transposed_m, center=TRUE)
  # var_explained <- pca_general$sdev^2/sum(pca_general$sdev^2)
  # scores = as.data.frame(pca_general$x)
  # scores['Basename'] <- row.names(scores)
  # scores <- merge(scores, pheno, by = 'Basename')
  # title <- paste(output, "_m", sep="")
  # output_path <- paste(output_dir, title, sep="")
  # #jpeg(paste(output_path, "_PCA.jpeg", sep=""), quality = 100)
  # plot <- ggplot(data=scores, mapping = aes(x = PC1, y = PC2, color=eGFR)) +theme_bw() + geom_point(aes(shape=time), alpha=0.5, size=2)+ labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"), y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) + scale_color_manual(values=pal) + ggtitle(paste(title, "_PCA", sep="")) + geom_text(aes(label = sample_id), size=1.75, colour="black")
  # #plot + geom_text(aes(label = sample_name), size=3.5) + xlim(-100, 400)
  # print(plot)
  # #dev.off()
  return('Clustering\n')
}

if (file.exists(paste(output_dir, "beta_values.csv", sep=""))) {
  cat('Skip writing beta to CSV\n')
  beta_values_filtered <- import(paste(output_dir, "beta_values.csv", sep=""))
  row.names(beta_values_filtered) <- beta_values_filtered$V1
  beta_values_filtered <- beta_values_filtered[, 2:ncol(beta_values_filtered)]
} else {
  write.csv(beta_values_filtered, file = paste(output_dir, "beta_values.csv", sep=""), row.names = TRUE)
  write.csv(m_values, file = paste(output_dir, "m_values.csv", sep=""), row.names = TRUE)
} 

generate_dendro <- function(beta, pheno, timepoint){
  beta_t <- data.frame(t(beta))
  row.names(pheno) <- pheno$Basename
  if (timepoint == 'K1') {
    pheno <- pheno[(pheno$time == 'K1'),]
  } else if (timepoint == 'K2') {
    pheno <- pheno[(pheno$time == 'K2'),]
  } else {
    cat("Skip filtering down")
  }
  
  pheno <- pheno %>% select(c('sample_id', 'donor_age', 'eGFR_1month', 'eGFR_12month', 'eGFR_24month'))
  
  final_beta <- merge(pheno, beta_t, by='row.names')
  
  #changing the row.names accordingly 
  final_beta['sample_id_age'] <- 'na'
  final_beta <- final_beta[,c(ncol(final_beta),1:(ncol(final_beta)-1))]
  final_beta['sample_id_eGFR1'] <- 'na'
  final_beta <- final_beta[,c(ncol(final_beta),1:(ncol(final_beta)-1))]
  final_beta['sample_id_eGFR12'] <- 'na'
  final_beta <- final_beta[,c(ncol(final_beta),1:(ncol(final_beta)-1))]
  final_beta['sample_id_eGFR24'] <- 'na'
  final_beta <- final_beta[,c(ncol(final_beta),1:(ncol(final_beta)-1))]
  
  for (row in 1:nrow(final_beta)) {
    age <- paste(final_beta[row, 'sample_id'], final_beta[row, 'donor_age'], sep = " ")
    eGFR1 <- paste(final_beta[row, 'sample_id'], final_beta[row, 'eGFR_1month'], sep = " ")
    eGFR12 <- paste(final_beta[row, 'sample_id'], final_beta[row, 'eGFR_12month'], sep = " ")
    eGFR24 <- paste(final_beta[row, 'sample_id'], final_beta[row, 'eGFR_24month'], sep = " ")
    final_beta[row, 'sample_id_age'] <- age
    final_beta[row, 'sample_id_eGFR1'] <- eGFR1
    final_beta[row, 'sample_id_eGFR12'] <- eGFR12
    final_beta[row, 'sample_id_eGFR24'] <- eGFR24
  }
  
  dendro_out <- paste(timepoint, 'dendrograms.pdf', sep=" ")
  pdf(file = paste(output_dir, dendro_out, sep=""), width = 12, height = 8)
  row.names(final_beta) <- final_beta$sample_id_age
  clusters <- hclust(dist(final_beta[, 11:ncol(final_beta)]))
  plot(clusters, xlab = "Sample ID and Donor Age", main= paste(timepoint, "- age Dendrogram", sep = " "))
  
  row.names(final_beta) <- final_beta$sample_id_eGFR1
  clusters <- hclust(dist(final_beta[, 11:ncol(final_beta)]))
  plot(clusters, xlab = "Sample ID and 1month eGFR Status", main= paste(timepoint, "- eGFR 1month Dendrogram", sep = " "))
  
  row.names(final_beta) <- final_beta$sample_id_eGFR12
  clusters <- hclust(dist(final_beta[, 11:ncol(final_beta)]))
  plot(clusters, xlab = "Sample ID and 12month eGFR Status", main= paste(timepoint, "- eGFR 12month Dendrogram", sep = " "))
  
  row.names(final_beta) <- final_beta$sample_id_eGFR24
  clusters <- hclust(dist(final_beta[, 11:ncol(final_beta)]))
  plot(clusters, xlab = "Sample ID and 24month eGFR Status", main= paste(timepoint, "- eGFR 24month Dendrogram", sep = " "))
  
  dev.off()
}

timeList <- c('K1-K2')
for (time in timeList) {
  generate_dendro(beta_values_filtered, pheno_df, time)
}

q()

eGFR_List <- c('eGFR_1month', 'eGFR_12month', 'eGFR_24month')
#eGFR_List <- c('eGFR_1month')
comp_List <- c('K1_Low_K1_High', 'K2_Low_K2_High', 'K1_High_K2_High', 'K1_Low_K2_Low', 'K1_High_K2_Low', 'K1_Low_K2_High')
#comp_List <- c('K1_Low_K1_High')

pdf(file = paste(output_dir, "out.pdf", sep=""))
for (comp in comp_List) {
  for (outcome in eGFR_List) {
    clustering(pheno_df, outcome, comp, beta_values_filtered)
  }
}
dev.off()

eGFR_List <- c('eGFR_1month', 'eGFR_12month', 'eGFR_24month')
#eGFR_List <- c('eGFR_1month')
comp_List <- c('K1_Low_K1_High', 'K2_Low_K2_High', 'K1_High_K2_High', 'K1_Low_K2_Low', 'K1_High_K2_Low', 'K1_Low_K2_High')
#comp_List <- c('K1_Low_K1_High')

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

for (comp in comp_List) {
  for (outcome in eGFR_List) {
    cat("Identify CpGs\n")
    if (outcome == 'eGFR_1month') {
      pheno <- subset(pheno_df, select = -c(eGFR_12month, eGFR_24month))
      pheno <- pheno %>% rename('eGFR' = 'eGFR_1month')
    } else if (outcome == 'eGFR_12month') {
      pheno <- subset(pheno_df, select = -c(eGFR_1month, eGFR_24month))
      pheno <- pheno %>% rename('eGFR' = 'eGFR_12month')
    } else {
      pheno <- subset(pheno_df, select = -c(eGFR_1month, eGFR_12month))
      pheno <- pheno %>% rename('eGFR' = 'eGFR_24month')
    }
    
    if (comp == 'K1_Low_K1_High') {
      pheno <- pheno[((pheno$time == 'K1' & pheno$eGFR == 'Low') | (pheno$time == 'K1' & pheno$eGFR == 'High')),]
    } else if (comp == 'K2_Low_K2_High') {
      pheno <- pheno[((pheno$time == 'K2' & pheno$eGFR == 'Low') | (pheno$time == 'K2' & pheno$eGFR == 'High')),]
    } else if (comp == 'K1_High_K2_High') {
      pheno <- pheno[((pheno$time == 'K1' & pheno$eGFR == 'High') | (pheno$time == 'K2' & pheno$eGFR == 'High')),]
    } else if (comp == 'K1_Low_K2_Low') {
      pheno <- pheno[((pheno$time == 'K1' & pheno$eGFR == 'Low') | (pheno$time == 'K2' & pheno$eGFR == 'Low')),]
    } else if (comp == 'K1_High_K2_Low') {
      pheno <- pheno[((pheno$time == 'K1' & pheno$eGFR == 'High') | (pheno$time == 'K2' & pheno$eGFR == 'Low')),]
    } else if (comp == 'K1_Low_K2_High') {
      pheno <- pheno[((pheno$time == 'K1' & pheno$eGFR == 'Low') | (pheno$time == 'K2' & pheno$eGFR == 'High')),]
    } else {
      cat('Comparison request does not exist\n')
      q()
    }
    
    if (comp == 'K1_Low_K1_High') {
      condition <- factor(pheno$eGFR)
      cols <- c("High", "Low")
    } else if (comp == 'K2_Low_K2_High') {
      condition <- factor(pheno$eGFR)
      cols <- c("High", "Low")
    } else if (comp == 'K1_High_K2_High') {
      condition <- factor(pheno$time)
      cols <- c("K2", "K1")
    } else if (comp == 'K1_Low_K2_Low') {
      condition <- factor(pheno$time)
      cols <- c("K2", "K1")
    } else if (comp == 'K1_High_K2_Low') {
      condition <- factor(pheno$time)
      cols <- c("K2", "K1")
    } else {
      condition <- factor(pheno$time)
      cols <- c("K2", "K1")
    }
    
    #CpGs - DMPs
    
    beta_values_condition <- beta_values_filtered[,(colnames(beta_values_filtered) %in% pheno$Basename)]
    
    cat("Identify DMPs\n")
    design <- model.matrix(~0+condition, data=pheno)
    print(design)
    colnames(design) <- cols
    fit1 <- lmFit(beta_values_condition, design)
    
    if (comp == 'K1_Low_K1_High') {
      contMatrix <- makeContrasts(Low-High, levels=design)
    } else if (comp == 'K2_Low_K2_High') {
      contMatrix <- makeContrasts(Low-High, levels=design)
    } else if (comp == 'K1_High_K2_High') {
      contMatrix <- makeContrasts(K1-K2, levels=design)
    } else if (comp == 'K1_Low_K2_Low') {
      contMatrix <- makeContrasts(K1-K2, levels=design)
    } else if (comp == 'K1_High_K2_Low') {
      contMatrix <- makeContrasts(K1-K2, levels=design)
    } else {
      contMatrix <- makeContrasts(K1-K2, levels=design)
    }
    
    #contMatrix
    
    fit2 <- contrasts.fit(fit1, contMatrix)
    fit2 <- eBayes(fit2)
    
    #summary(decideTests(fit2))
    annEPICSub <- annEPIC[match(rownames(beta_values_condition),annEPIC$Name), c(1:4,12:19,24:ncol(annEPIC))]
    DMPs <- topTable(fit2, num=Inf, coef=1, genelist=annEPICSub)
    #plotCpg(m_values, cpg=rownames(DMPs)[1:4], pheno=type, ylab = "Beta values") #plots individual probes
    output <- paste(outcome, "_", sep="")
    title <- paste(output, comp, sep="")
    output <- paste(output_dir, title, sep="")
    write.csv(DMPs, file = paste(output, "_DMPs.csv", sep=""), row.names = TRUE)
    
    #Manhattan plot using the DMPs
    cat("Generating manhattan plot from DMPs\n")
    library(qqman)
    library(DMRcate)
    title <- paste(title, " (Adj. P-val)", sep="")
    col=c("black","grey")
    DMPs$chr = str_replace_all(DMPs$chr, 'chr', '')
    DMPs$chr = as.numeric(DMPs$chr)
    DMPs$pos = as.numeric(DMPs$pos)
    jpeg(paste(output, "_manhattan.jpeg", sep=""), quality = 90)
    manhattan(DMPs, chr="chr", bp="pos",, p="adj.P.Val", snp="Islands_Name", col=col, suggestiveline=(-log10(0.05)), main=title)
    dev.off()
    
    # myAnnotation <- cpg.annotate(object = as.matrix(m_values), datatype = "array", what = "M",
    #                              analysis.type = "differential", design = design, 
    #                              contrasts = TRUE, cont.matrix = contMatrix,
    #                              coef = "Low - High", arraytype = "450K")
    # 
    # DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
    # results.ranges <- extractRanges(DMRs, genome = "hg19")
    # write.csv(result.ranges, file=paste(output, "_DMRs.csv", sep=""), row.names = FALSE)
    rm(pheno, beta_values_condition, annEPICSub)
  }
}

q()


