#Script for DNA Methylation Analysis
library(minfi)
library(rio)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(limma)
library(dendextend)

#Example input: Rscript methylation_analysis.R <pheno_file> <paired_pheno_file> <base_dir> <git_dir> <output_dir>
args=commandArgs(trailingOnly=TRUE)
pheno_file = args[1]
pheno_file2 = args[2]
base_dir = args[3]
git_dir = args[4]
output_dir = args[5]

#Local Machine
# pheno_file = '/Users/adrianharris/Documents/dna_methylation_analysis/kidney_sample_sheet.csv'
# pheno_file2 = '/Users/adrianharris/Documents/dna_methylation_analysis/paired_kidney_sample_sheet.csv'
# base_dir = '/Users/adrianharris/Desktop/kidney/'
# git_dir = '/Users/adrianharris/Documents/dna_methylation_analysis/'
# output_dir = '/Users/adrianharris/Desktop/test_kidney/'

#Checking if output directory exists
if (file.exists(output_dir)) {
  cat("Directory already exists\n")
} else {
  dir.create(output_dir)
}

#Importing manually-curated sample sheet
pheno_df <- import(pheno_file)

#Must remove outlier sample, 203504430032_R01C01 (and its paired sample 203504430032-R02C01)
pheno_df <- pheno_df[!(pheno_df$Basename == '203504430032_R01C01' | pheno_df$Basename == '203504430032_R02C01'),]

pheno_df <- pheno_df[(pheno_df$array_type == 'EPIC'),]

#Number of files 
nrow(subset(pheno_df, time == 'K1'))
nrow(subset(pheno_df, time == 'K2'))

#Specify the directories and read in respective IDAT Files 
dirEPIC <- paste(base_dir, "EPIC_array", sep="")
if (file.exists(paste(output_dir, "rgSet.RDS", sep=""))) {
  cat('Loading in rgSet\n')
  rgSet <- readRDS(paste(output_dir, "rgSet.RDS", sep=""))
} else {
  cat('Generate rgSet\n')
  rgSet <- read.metharray.exp(base=dirEPIC, target=pheno_df, force=TRUE)
  #Saving set as an RDS file 
  saveRDS(rgSet, file = paste(output_dir, "rgSet.RDS", sep=""))
}

pheno_df <- data.frame(pData(rgSet))

#Color Scheme Defined 
pal <- brewer.pal(4,"Dark2")

#Complete the detection of p-values if necessary
if (file.exists(paste(output_dir, "beta_values.csv", sep=""))) {
  cat('Skipping p-value step\n')
} else if (file.exists(paste(output_dir, "p-values.csv", sep=""))) {
  cat('P-values.csv already exists\n')
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
  barplot((subset(detP_df, array_type == 'EPIC'))$p_value_mean, col=pal[factor(detP_df$time)], names.arg=(subset(detP_df, array_type == 'EPIC'))$sample_id, las=2, cex.names=0.4, cex.axis=0.5, space=0.5, ylab="Mean detection p-values", main='EPIC Array')
  legend("topleft", legend=levels(factor(detP_df$time)), fill=pal,
         cex=0.27, bty = "n", bg="white")
  
  rm(detP_df, detP)
  dev.off()
  par(mfrow=c(1,1))
}

#Perform preprocessing if necessary
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
  jpeg(paste(output_dir, "preprocessQC.jpeg", sep=""), quality = 100)
  plot <- ggplot(data=qc, mapping = aes(x = mMed, y = uMed, color=time)) + geom_point(aes(shape=array_type), alpha=0.5) + xlim(7, 14) + ylim(7, 14) + theme_bw()+ geom_abline(slope=-1, intercept = 21.25 , color="black", linetype="dashed", size=0.5) + scale_color_manual(values=pal)
  print(plot)
  dev.off()
  jpeg(paste(output_dir, "preprocessDensity.jpeg", sep=""), quality = 100)
  densityPlot(mtSet, sampGroups =qc$array_type)
  dev.off()
  rm(mtSet, qc, plot)
}

#Normalize if necessary
if (file.exists(paste(output_dir, "postNormQC.jpeg", sep=""))) {
  cat('Loading normalization\n')
  mtSet <- readRDS(paste(output_dir, "mtSet.RDS", sep=""))
} else {
  cat('Performing normalization and plotting\n')
  #Normalization and plotting 
  #mtSet <- preprocessNoob(rgSet)
  mtSet <- preprocessSWAN(rgSet)
  saveRDS(mtSet, file = paste(output_dir, "mtSet.RDS", sep=""))
  postqc <- getQC(mtSet)
  postqc <- data.frame(postqc)
  postqc['Basename'] <- row.names(postqc)
  postqc <- merge(postqc, pheno_df, by = 'Basename')
  jpeg(paste(output_dir, "postNormQC.jpeg", sep=""), quality = 100)
  plot2 <- ggplot(data=postqc, mapping = aes(x = mMed, y = uMed, color=time)) + geom_point(aes(shape=array_type), alpha=0.5) + xlim(5, 14) + ylim(5, 14) + theme_bw()+ geom_abline(slope=-1, intercept = 21.25 , color="black", linetype="dashed", size=0.5) + scale_color_manual(values=pal)
  print(plot2)
  dev.off()
  jpeg(paste(output_dir, "postNormDensity.jpeg", sep=""), quality = 100)
  densityPlot(mtSet, sampGroups = postqc$array_type)
  dev.off()
  rm(postqc, plot2)
}

library(lumi) #Load in for beta2m
# Map to genome-generate betas if necessary
if (file.exists(paste(output_dir, "beta_values.csv", sep=""))) {
  cat('Skip genomic methyl set\n')
}else {
  cat('Converting to Genomic Methyl Set\n')
  rSet <- ratioConvert(mtSet, what = "both", keepCN = TRUE)
  gmtSet <- mapToGenome(rSet)
  print(dim(gmtSet)) #Number of probes = 865859
  
  #Removing probes that include SNPs
  snps <- getSnpInfo(gmtSet)
  gmtSet <- addSnpInfo(gmtSet)
  gr <- granges(gmtSet)
  gmtSet <- dropLociWithSnps(gmtSet, snps=c("SBE","CpG"), maf=0)
  print(dim(gmtSet)) #Number of probes = 835424
  rm(snps, gr)
  
  # Filter unwanted sites - ensure probes are in the same order in the gmtSet and detP objects
  detP <- detectionP(rgSet)
  detP <- detP[match(featureNames(gmtSet),rownames(detP)),] 
  
  # remove any probes that have failed in >50% of samples
  keep <- detP < 0.01
  keep.probes <- rownames(detP[rowMeans(keep)>=0.5,]) #exclude probes that failed detection in more than half of the samples
  gmtSet <- gmtSet[keep.probes,] 
  print(dim(gmtSet)) #Number of probes = 835364 
  rm(keep.probes)
  
  #Remove probes that located on the X or Y chromosome
  annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  
  keep <- !(featureNames(gmtSet) %in% annEPIC$Name[annEPIC$chr %in% c("chrX","chrY")]) #remove probes that are of the chrom x or y
  table(keep)
  gmtSet <- gmtSet[keep,]
  print(dim(gmtSet)) #Number of probes = 816068
  rm(annEPIC, keep)
  
  #Creation of bad probes character and filter gmtSet
  cross.react <- read.csv(paste(git_dir, '48639-non-specific-probes-Illumina450k.csv', sep=""), head = T, as.is = T)
  cross.react.probes <- as.character(cross.react$TargetID)
  #Probes identified with potential hybridization issues
  multi.map <- read.csv(paste(git_dir, 'HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt', sep=""), head = F, as.is = T)
  multi.map.probes <- as.character(multi.map$V1)
  # Determine unique probes between the cross-reactive and multi-map probes
  bad.probes <- unique(c(cross.react.probes, multi.map.probes))
  keep <- !(featureNames(gmtSet) %in% bad.probes) #Removal of these bad probes
  table(keep)
  gmtSet <- gmtSet[keep,]
  print(dim(gmtSet)) #Number of probes = 783224
  rm(cross.react, multi.map, bad.probes, cross.react.probes, multi.map.probes, keep)
  
  #Extract betas and m_values
  beta_values <- getBeta(gmtSet)
  # m_values <- getM(gmtSet)
  
  #Remove probes Hypomethylated in all samples identified by CpG sites having beta < 0.05 in all samples
  beta_values_filtered <- as.data.frame(beta_values) 
  cat('Hypomethylated\n')
  print(dim(filter_all(beta_values_filtered, all_vars(. < 0.05)))) #Number of hypomethylated
  beta_values_filtered <- filter_all(beta_values_filtered, any_vars(. >= 0.05)) 
  #20844 Hypomethylated
  
  #Remove probes Hypermethylated in all samples identified by CpG sites having beta > 0.95 in all samples
  cat("Hypermethylated\n")
  print(dim(filter_all(beta_values_filtered, all_vars(. > 0.95))))
  #334 hypermethylated
  beta_values_filtered <- filter_all(beta_values_filtered, any_vars(. < 0.95)) 
  cat("Final Dimensions\n")
  dim(beta_values_filtered)
}

#Writing beta to CSV if necessary
if (file.exists(paste(output_dir, "beta_values.csv", sep=""))) {
  cat('Skip writing beta to CSV\n')
  beta_values_filtered <- import(paste(output_dir, "beta_values.csv", sep=""))
  row.names(beta_values_filtered) <- beta_values_filtered$V1
  beta_values_filtered <- beta_values_filtered[, 2:ncol(beta_values_filtered)]
} else {
  cat('Writing beta to CSV\n')
  write.csv(beta_values_filtered, file = paste(output_dir, "beta_values.csv", sep=""), row.names = TRUE)
} 

#Loading lumi and converting beta2m 
library(lumi)
m_values <- beta2m(beta_values_filtered)

if (file.exists(paste(output_dir, "m_values.csv", sep=""))) {
  cat('Skip writing of m_values\n')
} else {
  write.csv(m_values, file = paste(output_dir, "m_values.csv", sep=""), row.names = TRUE)
}

#Loading in the paired file
paired_pheno <- import(pheno_file2)
paired_pheno <- paired_pheno[!(paired_pheno$sample_id == 'KUT4_K2' | paired_pheno$sample_id == 'KUT4_K1'),]

#Making the phenotype dataframe the paired dataframe 
pheno_df <- pheno_df[(pheno_df$Basename %in% paired_pheno$Basename),]
beta_values_filtered <- beta_values_filtered[,(colnames(beta_values_filtered) %in% pheno_df$Basename)]

#Exclude DCD samples - Kidney data
beta_values_filtered <- beta_values_filtered[,!(colnames(beta_values_filtered) %in% c("9296930129_R05C01", "9305651174_R01C01", "9305651174_R03C01", "9305651174_R02C02", "9305651174_R03C02", "9305651191_R02C02", "9305651191_R04C02", "201465900002_R04C01", "202240580106_R03C01", "202240580208_R03C01", "202259340119_R05C01", "202259350016_R04C01", "203496240002_R03C01", "203504430032_R05C01", "204001300109_R07C01", "204001300109_R08C01", "204001350016_R01C01", "202702240079_R06C01"))]
#Removing rows based on the sample_name column in phenotype dataframe
pheno_df <- pheno_df[!(pheno_df$Basename %in% c("9296930129_R05C01", "9305651174_R01C01", "9305651174_R03C01", "9305651174_R02C02", "9305651174_R03C02", "9305651191_R02C02", "9305651191_R04C02", "201465900002_R04C01", "202240580106_R03C01", "202240580208_R03C01", "202259340119_R05C01", "202259350016_R04C01", "203496240002_R03C01", "203504430032_R05C01", "204001300109_R07C01", "204001300109_R08C01", "204001350016_R01C01", "202702240079_R06C01")),]

m_values <- beta2m(beta_values_filtered)

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
  
  #Filtering the dataframe down
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
  }
  
  #Filter beta dataframe using column names for the relevant comparison
  beta_values_filtered <- beta_values_filtered[,(colnames(beta_values_filtered) %in% pheno$Basename)]
  # dim(beta_values_filtered)
  
  cat("PCA plots - Betas\n")
  #Transpose betas and generate PCA plot
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
  #plot + geom_text(aes(label = sample_id), size=3.5) + xlim(-100, 400)
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
  
  #m_values <- beta2m(betas)
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
  # #plot + geom_text(aes(label = sample_id), size=3.5) + xlim(-100, 400)
  # print(plot)
  # #dev.off()
  return('Clustering\n')
}

#eGFR_List <- c('eGFR_1month', 'eGFR_12month', 'eGFR_24month')
#comp_List <- c('Low_High')
eGFR_List <- c('eGFR_1month')
comp_List <- c('K1_Low_K1_High', 'K2_Low_K2_High', 'K1_High_K2_High', 'K1_Low_K2_Low', 'K1_High_K2_Low', 'K1_Low_K2_High')


# pdf(file = paste(output_dir, "eGFR1month_comparisons.pdf", sep=""))
# for (comp in comp_List) {
#   for (outcome in eGFR_List) {
#     clustering(pheno_df, outcome, comp, beta_values_filtered)
#   }
# }
# 
# dev.off()

#Generate dendrogram 
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
  
  #Creating new labels/corresponding columns and moving toward front of dataframe
  final_beta['sample_id_age'] <- 'na'
  final_beta <- final_beta[,c(ncol(final_beta),1:(ncol(final_beta)-1))]
  final_beta['sample_id_eGFR1'] <- 'na'
  final_beta <- final_beta[,c(ncol(final_beta),1:(ncol(final_beta)-1))]
  final_beta['sample_id_eGFR12'] <- 'na'
  final_beta <- final_beta[,c(ncol(final_beta),1:(ncol(final_beta)-1))]
  final_beta['sample_id_eGFR24'] <- 'na'
  final_beta <- final_beta[,c(ncol(final_beta),1:(ncol(final_beta)-1))]
  
  #Filling these new columns
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
  
  #Generate dendrograms for each label
  dendro_out <- paste(timepoint, 'dendrograms.pdf', sep=" ")
  pdf(file = paste(output_dir, dendro_out, sep=""), width = 12, height = 8)
  row.names(final_beta) <- final_beta$sample_id_age
  clusters <- hclust(dist(final_beta[, 11:ncol(final_beta)]))
  dend <- as.dendrogram(clusters)
  dend <- set(dend, "labels_cex", 0.3)
  plot(dend, xlab = "Sample ID and Donor Age", ylab="Height", main= paste(timepoint, "- age Dendrogram", sep = " "))
  
  row.names(final_beta) <- final_beta$sample_id_eGFR1
  clusters <- hclust(dist(final_beta[, 11:ncol(final_beta)]))
  dend <- as.dendrogram(clusters)
  dend <- set(dend, "labels_cex", 0.3)
  plot(dend, xlab = "Sample ID and 1month eGFR Status", ylab="Height", main= paste(timepoint, "- eGFR 1month Dendrogram", sep = " "))
  
  row.names(final_beta) <- final_beta$sample_id_eGFR12
  clusters <- hclust(dist(final_beta[, 11:ncol(final_beta)]))
  dend <- as.dendrogram(clusters)
  dend <- set(dend, "labels_cex", 0.3)
  plot(dend, xlab = "Sample ID and 12month eGFR Status", ylab="Height",  main= paste(timepoint, "- eGFR 12month Dendrogram", sep = " "))
  
  row.names(final_beta) <- final_beta$sample_id_eGFR24
  clusters <- hclust(dist(final_beta[, 11:ncol(final_beta)]))
  dend <- as.dendrogram(clusters)
  dend <- set(dend, "labels_cex", 0.3)
  plot(dend, xlab = "Sample ID and 24month eGFR Status", ylab="Height", main= paste(timepoint, "- eGFR 24month Dendrogram", sep = " "))
  
  dev.off()
}

# timeList <- c('K1-K2')
# for (time in timeList) {
#   generate_dendro(beta_values_filtered, pheno_df, time)
# }

#Function for generating manhattan plots
library(qqman)
library(DMRcate)
generate_man <- function(DMPs, comp, status) {
  #Manhattan plot using the DMPs
  cat("Generating manhattan plot from DMPs\n")
  title <- paste(status, comp, sep="-")
  title_fig <- paste(title, " (Adj. P-val)", sep="")
  col=c("black","grey")
  DMPs$chr = str_replace_all(DMPs$chr, 'chr', '')
  DMPs$chr = as.numeric(DMPs$chr)
  DMPs$pos = as.numeric(DMPs$pos)
  output <- paste(title, "_manhattan.jpeg", sep="")
  jpeg(paste(output_dir, output, sep=""), quality = 100)
  manhattan(DMPs, chr="chr", bp="pos", p="adj.P.Val", snp="Islands_Name", col=col, suggestiveline=(-log10(0.05)), main=title_fig)
  dev.off()
  title <- paste(status, comp, sep="-")
  title_fig <- paste(title, " (deltaBetas)", sep="")
  output <- paste(title, "_manhattan_deltas.jpeg", sep="")
  jpeg(paste(output_dir, output, sep=""), quality = 100)
  manhattan(DMPs, chr="chr", bp="pos", logp=FALSE, p="deltaBeta", snp="Islands_Name", col=col, suggestiveline=(0.15), main=title_fig)
  dev.off()
  return("Done with manhattan plot")
}

#Copying of the betas dataframe
newBeta_df <- beta_values_filtered

#Changing the column names to sample_id
colnames(newBeta_df) <- pheno_df$sample_id

#Calculating average delta beta per comparison
get_deltaBeta <- function(cond1, cond2, phenotype) {
  pheno_condition1 <- phenotype[(phenotype$condition == cond1),]
  pheno_condition2 <- phenotype[(phenotype$condition == cond2),]
  print(dim(pheno_condition1))
  betas_condition1 <- newBeta_df[,(colnames(newBeta_df) %in% pheno_condition1$sample_id)]
  betas_condition2 <- newBeta_df[,(colnames(newBeta_df) %in% pheno_condition2$sample_id)]
  print(dim(betas_condition1))
  betas_condition1['average'] <- rowSums(betas_condition1[,1:ncol(betas_condition1)])
  betas_condition2['average'] <- rowSums(betas_condition2[,1:ncol(betas_condition2)])
  betas_condition1$average <-as.numeric(as.character(betas_condition1$average)) / (nrow(pheno_condition1))
  betas_condition2$average <-as.numeric(as.character(betas_condition2$average)) / (nrow(pheno_condition2))
  betas_condition1$Name <- row.names(betas_condition1)
  betas_condition2$Name <- row.names(betas_condition2)
  betas_condition1 <- betas_condition1[,c(ncol(betas_condition1), (ncol(betas_condition1)-1))]
  betas_condition2 <- betas_condition2[,c(ncol(betas_condition2), (ncol(betas_condition2)-1))]
  betas_condition <- merge(betas_condition1, betas_condition2, by = "Name")
  betas_condition['deltaBeta'] <- (betas_condition$average.y - betas_condition$average.x)
  betas_condition <- betas_condition[,c(1, ncol(betas_condition))]
  return(betas_condition)
}

#Load library for EPIC 
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

#Vector of eGFR statuses
eGFR_List <- c('eGFR_1month', 'eGFR_12month', 'eGFR_24month')
#eGFR_List <- c('eGFR_1month')
for (outcome in eGFR_List) {
  cat("Identify CpGs\n")
  log_df <- data.frame(comparison = character(), number_of_samples = double(), number_of_sig_DMPs = double())
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
  
  #Addition of condition column 
  pheno$condition <- 'NA'
  
  #Filling of condition column
  for (row in 1:nrow(pheno)) {
    if (pheno[row, 'time'] == 'K1' & pheno[row, 'eGFR'] == 'High') {
      pheno[row, 'condition'] <- "K1_High"
    } else if (pheno[row, 'time'] == 'K1' & pheno[row, 'eGFR'] == 'Low') {
      pheno[row, 'condition'] <- "K1_Low"
    } else if (pheno[row, 'time'] == 'K2' & pheno[row, 'eGFR'] == 'High') {
      pheno[row, 'condition'] <- "K2_High"
    } else if (pheno[row, 'time'] == 'K2' & pheno[row, 'eGFR'] == 'Low') {
      pheno[row, 'condition'] <- "K2_Low"
    }
  }

  # create design matrix
  design <- model.matrix(~0+condition, data=pheno)
  colnames(design) <- c("K1_High", "K1_Low", "K2_High", "K2_Low")
  
  # fit the linear model 
  fit1 <- lmFit(m_values, design)
  # create a contrast matrix for specific comparisons
  contMatrix <- makeContrasts(K1_Low-K1_High,
                              K2_Low-K2_High,
                              K1_High-K2_High,
                              K1_Low-K2_Low,
                              K1_High-K2_Low,
                              K1_Low-K2_High,
                              levels=design)
  
  contMatrix

  # fit the contrasts
  fit2 <- contrasts.fit(fit1, contMatrix)
  fit2 <- eBayes(fit2)
  
  # look at the numbers of DM CpGs at FDR < 0.05
  summary(decideTests(fit2))

  annEPICSub <- annEPIC[match(rownames(m_values),annEPIC$Name),
                        c(1:4,12:19,24:ncol(annEPIC))]
  
  #Identifying and writing output for DMPs
  DMPs1 <- topTable(fit2, num=Inf, coef=1, genelist=annEPICSub)
  DMPs1 <- data.frame(DMPs1)
  deltaBeta_df <- get_deltaBeta("K1_Low", "K1_High", pheno)
  sample_num <- nrow(subset(pheno, (condition == 'K1_Low' | condition == 'K1_High')))
  DMPs1 <- merge(DMPs1, deltaBeta_df, by = 'Name')
  output <- paste(outcome, "K1_Low-K1_High_DMPs.csv", sep="-")
  write.csv(DMPs1, file = paste(output_dir, output, sep=""), row.names = FALSE) 
  DMPs1_sig <- DMPs1[(DMPs1$adj.P.Val < 0.05 & DMPs1$deltaBeta > 0.15),]
  print(dim(DMPs1_sig))
  output <- paste(outcome, "K1_Low-K1_High_DMPs_sig.csv", sep = "-")
  write.csv(DMPs1_sig, file = paste(output_dir, output, sep=""), row.names = FALSE)
  log_df[nrow(log_df) + 1,] <- c("K1_Low-K1_High", sample_num, nrow(DMPs1_sig))
  
  DMPs2 <- topTable(fit2, num=Inf, coef=2, genelist=annEPICSub)
  DMPs2 <- data.frame(DMPs2)
  deltaBeta_df <- get_deltaBeta("K2_Low", "K2_High", pheno)
  sample_num <- nrow(subset(pheno, (condition == 'K2_Low' | condition == 'K2_High')))
  DMPs2 <- merge(DMPs2, deltaBeta_df, by = 'Name')
  output <- paste(outcome, "K2_Low-K2_High_DMPs.csv", sep="-")
  write.csv(DMPs2, file = paste(output_dir, output, sep=""), row.names = FALSE) 
  DMPs2_sig <- DMPs2[(DMPs2$adj.P.Val < 0.05 & DMPs2$deltaBeta > 0.15),]
  print(dim(DMPs2_sig))
  output <- paste(outcome, "K2_Low-K2_High_DMPs_sig.csv", sep="-")
  write.csv(DMPs2_sig, file = paste(output_dir, output, sep=""), row.names = FALSE)
  log_df[nrow(log_df) + 1,] <- c("K2_Low-K2_High", sample_num, nrow(DMPs2_sig))
  
  DMPs3 <- topTable(fit2, num=Inf, coef=3, genelist=annEPICSub)
  DMPs3 <- data.frame(DMPs3)
  deltaBeta_df <- get_deltaBeta("K1_High", "K2_High", pheno)
  sample_num <- nrow(subset(pheno, (condition == 'K1_High' | condition == 'K2_High')))
  DMPs3 <- merge(DMPs3, deltaBeta_df, by = 'Name')
  output <- paste(outcome, "K1_High-K2_High_DMPs.csv", sep="-")
  write.csv(DMPs3, file = paste(output_dir, output, sep=""), row.names = FALSE) 
  DMPs3_sig <- DMPs3[(DMPs3$adj.P.Val < 0.05 & DMPs3$deltaBeta > 0.15),]
  print(dim(DMPs3_sig))
  output <- paste(outcome, "K1_High-K2_High_DMPs_sig.csv", sep="-")
  write.csv(DMPs3_sig, file = paste(output_dir, output, sep=""), row.names = FALSE)
  log_df[nrow(log_df) + 1,] <- c("K1_High-K2_High", sample_num, nrow(DMPs3_sig))
  
  DMPs4 <- topTable(fit2, num=Inf, coef=4, genelist=annEPICSub)
  DMPs4 <- data.frame(DMPs4)
  deltaBeta_df <- get_deltaBeta("K1_Low", "K2_Low", pheno)
  sample_num <- nrow(subset(pheno, (condition == 'K1_Low' | condition == 'K2_Low')))
  DMPs4 <- merge(DMPs4, deltaBeta_df, by = 'Name')
  output <- paste(outcome, "K1_Low-K2_Low_DMPs.csv", sep="-")
  write.csv(DMPs4, file = paste(output_dir, output, sep=""), row.names = FALSE) 
  DMPs4_sig <- DMPs4[(DMPs4$adj.P.Val < 0.05 & DMPs4$deltaBeta > 0.15),]
  print(dim(DMPs4_sig))
  output <- paste(outcome, "K1_Low-K2_Low_DMPs_sig.csv", sep="-")
  write.csv(DMPs4_sig, file = paste(output_dir, output, sep=""), row.names = FALSE)
  log_df[nrow(log_df) + 1,] <- c("K1_Low-K2_Low", sample_num, nrow(DMPs4_sig))
  
  DMPs5 <- topTable(fit2, num=Inf, coef=5, genelist=annEPICSub)
  DMPs5 <- data.frame(DMPs5)
  deltaBeta_df <- get_deltaBeta("K1_High", "K2_Low", pheno)
  sample_num <- nrow(subset(pheno, (condition == 'K1_High' | condition == 'K2_Low')))
  DMPs5 <- merge(DMPs5, deltaBeta_df, by = 'Name')
  output <- paste(outcome, "K1_High-K2_Low_DMPs.csv", sep="-")
  write.csv(DMPs5, file = paste(output_dir, output, sep=""), row.names = FALSE) 
  DMPs5_sig <- DMPs5[(DMPs5$adj.P.Val < 0.05 & DMPs5$deltaBeta > 0.15),]
  print(dim(DMPs5_sig))
  output <- paste(outcome, "K1_High-K2_Low_DMPs_sig.csv", sep="-")
  write.csv(DMPs5_sig, file = paste(output_dir, output, sep=""), row.names = FALSE)
  log_df[nrow(log_df) + 1,] <- c("K1_High-K2_Low", sample_num, nrow(DMPs5_sig))
  
  DMPs6 <- topTable(fit2, num=Inf, coef=6, genelist=annEPICSub)
  DMPs6 <- data.frame(DMPs6)
  deltaBeta_df <- get_deltaBeta("K1_Low", "K2_High", pheno)
  sample_num <- nrow(subset(pheno, (condition == 'K1_Low' | condition == 'K2_High')))
  DMPs6 <- merge(DMPs6, deltaBeta_df, by = 'Name')
  output <- paste(outcome, "K1_Low-K2_High_DMPs.csv", sep="-")
  write.csv(DMPs6, file = paste(output_dir, output, sep=""), row.names = FALSE) 
  DMPs6_sig <- DMPs6[(DMPs6$adj.P.Val < 0.05 & DMPs6$deltaBeta > 0.15),]
  print(dim(DMPs6_sig))
  output <- paste(outcome, "K1_Low-K2_High_DMPs_sig.csv", sep="-")
  write.csv(DMPs6_sig, file = paste(output_dir, output, sep=""), row.names = FALSE)
  log_df[nrow(log_df) + 1,] <- c("K1_Low-K2_High", sample_num, nrow(DMPs6_sig))
  
  final_title <- paste(outcome, "_EPIC_kidney_log.csv", sep="")
  write.csv(log_df, file = paste(output_dir, final_title, sep=""), row.names = FALSE)
  
  generate_man(DMPs1, 'K1_Low-K1_High', outcome)
  generate_man(DMPs2, 'K2_Low-K2_High', outcome)
  generate_man(DMPs3, 'K1_High-K2_High', outcome)
  generate_man(DMPs4, 'K1_Low-K2_Low', outcome)
  generate_man(DMPs5, 'K1_High-K2_Low', outcome)
  generate_man(DMPs6, 'K1_Low-K2_High', outcome)

}

q()

