#Script for DNA Methylation Analysis
library(minfi)
library(rio)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(dendextend)
library(limma)
library(sva)
sessionInfo()

#Example input: Rscript liver_methylation_analysis.R <pheno_file> <base_dir> <git_dir> <output_dir>
args=commandArgs(trailingOnly=TRUE)
pheno_file = args[1]
base_dir = args[2]
git_dir = args[3]
output_dir = args[4]

#Local Machine
# pheno_file = '/Users/adrianharris/Documents/dna_methylation_analysis/liver_sample_sheet.csv'
# base_dir = '/Users/adrianharris/Desktop/liver/'
# git_dir = '/Users/adrianharris/Documents/dna_methylation_analysis/'
# output_dir = '/Users/adrianharris/Desktop/liver/'

#Checking if the output directory exists
if (file.exists(output_dir)) {
  cat("Directory already exists\n")
} else {
  dir.create(output_dir)
}

#Importing manually-curated sample sheet
pheno_df <- import(pheno_file)

#Drop 2 samples that are not in L1 or L2 
pheno_df <- pheno_df[!(pheno_df$Basename == '203751390020_R08C01' | pheno_df$Basename == '203751390020_R01C01'),]

#Update sample_group column 
for (row in 1:nrow(pheno_df)) {
  if (pheno_df[row, 'sample_group'] == 'High Injury') {
    injury <- 'High'
  } else {
    injury <- 'Low'
  }
  pheno_df[row, 'sample_group'] <- injury
}

nrow(subset(pheno_df, array_type == '450K'))
nrow(subset(pheno_df, array_type == 'EPIC'))

#Split into phenotype file 
pheno450k <- pheno_df[!(pheno_df$array_type == 'EPIC'),]
phenoEPIC <- pheno_df[!(pheno_df$array_type == '450K'),]

#Specify the directories and read in respective IDAT Files 
dir450k <- paste(base_dir, "450k_array", sep="")
dirEPIC <- paste(base_dir, "EPIC_array", sep="")

#Load or generate the rgSet if necessary
if (file.exists(paste(output_dir, "rgSet.RDS", sep=""))) {
  cat('Loading in rgSet (combined)\n')
  rgSet <- readRDS(paste(output_dir, "rgSet.RDS", sep=""))
} else {
  cat('Generate rgSet\n')
  rgSet450k <- read.metharray.exp(base=dir450k, target=pheno450k)
  rgSetEPIC <- read.metharray.exp(base=dirEPIC, target=phenoEPIC, force=TRUE)
  rgSet <- combineArrays(rgSet450k, rgSetEPIC)
  ?read.metharray.exp()
  #Saving set as an RDS file 
  saveRDS(rgSet450k, file = paste(output_dir, "rgSet450k.RDS", sep=""))
  saveRDS(rgSetEPIC, file = paste(output_dir, "rgSetEPIC.RDS", sep=""))
  saveRDS(rgSet, file = paste(output_dir, "rgSet.RDS", sep=""))
  rm(rgSet450k)
  rm(rgSetEPIC)
}

#Simply ensure the correct order - If you ever run with a new pheno_df, a new rgSet should be created
pheno_df <- data.frame(pData(rgSet))

#Color scheme defined 
pal <- brewer.pal(4,"Dark2")

#Skip or generate p-values if necessary 
if (file.exists(paste(output_dir, "p-values.csv", sep=""))) {
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
  
  #Plotting 450K barplot 
  jpeg(paste(output_dir, "p_values450k.jpeg", sep=""), quality = 100)
  barplot((subset(detP_df, array_type == '450K'))$p_value_mean, col=pal[factor(detP_df$collection)], names.arg=(subset(detP_df, array_type == '450K'))$sample_name, las=2, cex.names=0.4, cex.axis=0.5, space=0.5, ylab="Mean detection p-values", main='450K Array')
  legend("topleft", legend=levels(factor(detP_df$collection)), fill=pal,
         cex=0.27, bty = "n", bg="white")
  dev.off()
  
  #Plotting EPIC barplot
  jpeg(paste(output_dir, "p_valuesEPIC.jpeg", sep=""), quality = 100)
  barplot((subset(detP_df, array_type == 'EPIC'))$p_value_mean, col=pal[factor(detP_df$collection)], names.arg=(subset(detP_df, array_type == 'EPIC'))$sample_name, las=2, cex.names=0.4, cex.axis=0.5, space=0.5, ylab="Mean detection p-values", main='EPIC Array')
  legend("topleft", legend=levels(factor(detP_df$collection)), fill=pal,
         cex=0.27, bty = "n", bg="white")
  dev.off()
  rm(detP_df, detP)
}

#Skip or perform preprocessing if necessary
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
  plot <- ggplot(data=qc, mapping = aes(x = mMed, y = uMed, color=collection)) + geom_point(aes(shape=array_type), alpha=0.5) + xlim(7, 14) + ylim(7, 14) + theme_bw()+ geom_abline(slope=-1, intercept = 21.25 , color="black", linetype="dashed", size=0.5) + scale_color_manual(values=pal)
  print(plot)
  dev.off()
  jpeg(paste(output_dir, "preprocessDensity.jpeg", sep=""), quality = 100)
  densityPlot(mtSet, sampGroups =qc$array_type)
  dev.off()
  rm(mtSet, qc, plot)
}

#Load or normalize if necessary
if (file.exists(paste(output_dir, "mtSet.RDS", sep=""))) {
  cat('Loading normalization\n')
  mtSet <- readRDS(paste(output_dir, "mtSet.RDS", sep=""))
} else {
  cat('Performing normalization and plotting\n')
  mtSet <- preprocessNoob(rgSet)
  saveRDS(mtSet, file = paste(output_dir, "mtSet.RDS", sep=""))
  postqc <- getQC(mtSet)
  postqc <- data.frame(postqc)
  postqc['Basename'] <- row.names(postqc)
  postqc <- merge(postqc, pheno_df, by = 'Basename')
  jpeg(paste(output_dir, "postNormQC.jpeg", sep=""), quality = 100)
  plot2 <- ggplot(data=postqc, mapping = aes(x = mMed, y = uMed, color=collection)) + geom_point(aes(shape=array_type), alpha=0.5) + xlim(5, 14) + ylim(5, 14) + theme_bw()+ geom_abline(slope=-1, intercept = 21.25 , color="black", linetype="dashed", size=0.5) + scale_color_manual(values=pal)
  print(plot2)
  dev.off()
  jpeg(paste(output_dir, "postNormDensity.jpeg", sep=""), quality = 100)
  densityPlot(mtSet, sampGroups = postqc$array_type)
  dev.off()
  rm(postqc, plot2)
}

#Map to genome-generate betas if necessary 
if (file.exists(paste(output_dir, "beta_values.csv", sep=""))) {
  cat('beta_values exists\n')
} else {
  cat('Converting to Genomic Methyl Set\n')
  rSet <- ratioConvert(mtSet, what = "both", keepCN = TRUE)
  gmtSet <- mapToGenome(rSet)
  print(dim(gmtSet))
  
  #Removing probes that include SNPs
  snps <- getSnpInfo(gmtSet)
  gmtSet <- addSnpInfo(gmtSet)
  gr <- granges(gmtSet)
  gmtSet <- dropLociWithSnps(gmtSet, snps=c("SBE","CpG"), maf=0)
  print(dim(gmtSet))
  rm(snps, gr)
  
  # Filter unwanted sites - ensure probes are in the same order in the gmtSet and detP objects
  detP <- detectionP(rgSet)
  detP <- detP[match(featureNames(gmtSet),rownames(detP)),] 
  
  # remove any probes that have failed in >50% of samples
  keep <- detP < 0.01
  keep.probes <- rownames(detP[rowMeans(keep)>=0.5,]) #exclude probes that failed detection in more than half of the samples
  gmtSet <- gmtSet[keep.probes,] 
  print(dim(gmtSet))
  rm(keep.probes)
  
  #Remove probes that located on the X or Y chromosome
  ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  
  keep <- !(featureNames(gmtSet) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")]) #remove probes that are of the chrom x or y
  table(keep)
  gmtSet <- gmtSet[keep,]
  print(dim(gmtSet))
  rm(ann450k, keep)
  
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
  print(dim(gmtSet))
  rm(cross.react, multi.map, bad.probes, cross.react.probes, multi.map.probes, keep)
  
  #Extract betas and m_values
  beta_values <- getBeta(gmtSet)
  #m_values <- getM(gmtSet)
  
  #Remove probes Hypomethylated in all samples identified by CpG sites having beta < 0.05 in all samples
  beta_values_filtered <- as.data.frame(beta_values) 
  cat('Hypomethylated\n')
  print(dim(filter_all(beta_values_filtered, all_vars(. < 0.05)))) 
  #345 hypomethylated 
  beta_values_filtered <- filter_all(beta_values_filtered, any_vars(. >= 0.05)) 
  
  #Remove probes Hypermethylated in all samples identified by CpG sites having beta > 0.95 in all samples
  cat("Hypermethylated\n")
  print(dim(filter_all(beta_values_filtered, all_vars(. > 0.95))))
  #5 hypermethylated
  beta_values_filtered <- filter_all(beta_values_filtered, any_vars(. < 0.95)) 
  cat("Final Dimensions\n")
  print(dim(beta_values_filtered))
}

#Writing beta to CSV if necessary
if (file.exists(paste(output_dir, "beta_values.csv", sep=""))) {
  cat('Loading beta CSV\n')
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

#Drop 'methylated' and 'unmethylated' sample names
pheno_df <- pheno_df[!(pheno_df$Basename %in% c('200999740023_R05C02', '200999740023_R06C02', '201004900096_R05C02', '201004900096_R06C02', '202702240054_R01C01', '202702240054_R02C01', '202702240079_R07C01', '202702240079_R08C01', '3999442124_R05C02', '3999442124_R06C02')),]

#Drop the one unpaired samples
pheno_df <- pheno_df[!(pheno_df$sample_name == 'V037L1'),]

#Drop two samples Haseeb dropped 
pheno_df <- pheno_df[!(pheno_df$sample_name == 'V024L1' | pheno_df$sample_name == 'V024L2'),]
pheno_df <- pheno_df[!(pheno_df$sample_name == 'V050 L1' | pheno_df$sample_name == 'V050'),]

#Drop DCD samples
pheno_df = pheno_df[!(pheno_df$donor_type == 'DCD'),]

#Exclude control and DCD samples for liver data 
beta_values_filtered <- beta_values_filtered[,(colnames(beta_values_filtered) %in% pheno_df$Basename)]
m_values <- m_values[, colnames(m_values) %in% colnames(beta_values_filtered)]

#Correction with combat()
# pheno_df$Basename 
# colnames(m_values)
batch <- pheno_df$array_type
modCombat <- model.matrix(~1, data=pheno_df)
m_values <- ComBat(dat=m_values, batch=batch, mod=modCombat)
beta_values_filtered <- data.frame(m2beta(m_values))
write.csv(beta_values_filtered, file = paste(output_dir, "beta_values_combat.csv", sep=""), row.names = TRUE)

#Horvath to categorical - Second Combat adjustment (according to old posts on biocond.)
# med <- median(as.numeric(pheno_df$Horvath))
# print(med)
# for (row in 1:nrow(pheno_df)) {
#   if (as.numeric(pheno_df[row, 'Horvath']) < med) {
#     pheno_df[row, 'Horvath'] <- 'Low'
#   } else {
#     pheno_df[row, 'Horvath'] <- 'High'
#   }
# }
# 
# batch <- pheno_df$Horvath
# modCombat <- model.matrix(~1, data=pheno_df)
# m_values <- ComBat(dat=m_values, batch=batch, mod=modCombat)
# beta_values_filtered <- data.frame(m2beta(m_values))
# write.csv(beta_values_filtered, file = paste(output_dir, "beta_values_age.csv", sep=""), row.names = TRUE)

#Function to generate several dendrograms with different labels
generate_dendro <- function(beta, pheno, timepoint){
  cat('Generating dendrograms\n')
  beta_t <- data.frame(t(beta)) #Transform betas
  row.names(pheno) <- pheno$Basename #Basename applied as rownames
  #Filtering down to a single timepoint if necessary
  if (timepoint == 'L1') {
    pheno <- pheno[(pheno$collection == 'L1'),]
  } else if (timepoint == 'L2') {
    pheno <- pheno[(pheno$collection == 'L2'),]
  } else {
    cat("Skip filtering down\n")
  }
  
  pheno <- pheno %>% select(c('sample_name', 'donor_age', 'donor_type', 'sample_group', 'array_type')) #Selecting necessary columns
  
  final_beta <- merge(pheno, beta_t, by='row.names')
  
  #Creating new labels/corresponding columns and moving toward front of dataframe
  final_beta['sample_name_age'] <- 'na'
  final_beta <- final_beta[,c(ncol(final_beta),1:(ncol(final_beta)-1))]
  final_beta['sample_name_donor'] <- 'na'
  final_beta <- final_beta[,c(ncol(final_beta),1:(ncol(final_beta)-1))]
  final_beta['sample_name_group'] <- 'na'
  final_beta <- final_beta[,c(ncol(final_beta),1:(ncol(final_beta)-1))]
  final_beta['sample_name_array'] <- 'na'
  final_beta <- final_beta[,c(ncol(final_beta),1:(ncol(final_beta)-1))]
  
  #Filling these new columns
  for (row in 1:nrow(final_beta)) {
    age <- paste(final_beta[row, 'sample_name'], final_beta[row, 'donor_age'], sep = " ")
    donor <- paste(final_beta[row, 'sample_name'], final_beta[row, 'donor_type'], sep = " ")
    group <- paste(final_beta[row, 'sample_name'], final_beta[row, 'sample_group'], sep = " ")
    array <- paste(final_beta[row, 'sample_name'], final_beta[row, 'array_type'], sep = " ")
    
    final_beta[row, 'sample_name_age'] <- age
    final_beta[row, 'sample_name_donor'] <- donor
    final_beta[row, 'sample_name_group'] <- group
    final_beta[row, 'sample_name_array'] <- array
  }
  
  #Dendrograms with new labels 
  dendro_out <- paste(timepoint, 'dendrograms.pdf', sep=" ")
  pdf(file = paste(output_dir, dendro_out, sep=""), width = 12, height = 8)
  row.names(final_beta) <- final_beta$sample_name_age
  clusters <- hclust(dist(final_beta[, 11:ncol(final_beta)]))
  dend <- as.dendrogram(clusters)
  dend <- set(dend, "labels_cex", 0.4)
  plot(dend, xlab = "Sample ID and Donor Age", ylab="Height", main= paste(timepoint, "- Donor Age Dendrogram", sep = " "))
  
  row.names(final_beta) <- final_beta$sample_name_donor
  clusters <- hclust(dist(final_beta[, 11:ncol(final_beta)]))
  dend <- as.dendrogram(clusters)
  dend <- set(dend, "labels_cex", 0.4)
  plot(dend, xlab = "Sample ID and Donor Status", ylab="Height", main= paste(timepoint, "- Donor Status Dendrogram", sep = " "))
  
  row.names(final_beta) <- final_beta$sample_name_group
  clusters <- hclust(dist(final_beta[, 11:ncol(final_beta)]))
  dend <- as.dendrogram(clusters)
  dend <- set(dend, "labels_cex", 0.4)
  plot(dend, xlab = "Sample ID and Injury Status", ylab="Height", main= paste(timepoint, "- Injury Status Dendrogram", sep = " "))
  
  row.names(final_beta) <- final_beta$sample_name_array
  clusters <- hclust(dist(final_beta[, 11:ncol(final_beta)]))
  dend <- as.dendrogram(clusters)
  dend <- set(dend, "labels_cex", 0.4)
  plot(dend, xlab = "Sample ID and Array Type", ylab="Height", main= paste(timepoint, "- Array Type Dendrogram", sep = " "))
  
  dev.off()
} 

#Loop through parameters to create dendrograms
# timeList <- c('L1', 'L2', 'L1-L2')
# for (time in timeList) {
#   generate_dendro(beta_values_filtered, pheno_df, time)
# }

clustering <- function(pheno, condition1, condition2, beta) {
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
  
  #Filter beta dataframe using column names for the relevant comparison
  beta <- beta[,(colnames(beta) %in% pheno$Basename)]
  
  cat("PCA plots - Betas\n")
  #Transpose beta, run PCA clustering and incorporate phenotype data
  transposed_beta <- t(beta)
  transposed_beta <- data.frame(transposed_beta)
  pca_general <- prcomp(transposed_beta, center=TRUE)
  var_explained <- pca_general$sdev^2/sum(pca_general$sdev^2)
  scores = as.data.frame(pca_general$x)
  scores['Basename'] <- row.names(scores)
  scores <- merge(scores, pheno, by = 'Basename')
  output <- paste(condition1, condition2, sep="-")
  title <- paste(output, "_betas", sep="")
  output_path <- paste(output_dir, title, sep="")
  jpeg(paste(output_path, "_PCA.jpeg", sep=""), quality = 100)
  plot <- ggplot(data=scores, mapping = aes(x = PC1, y = PC2, color=sample_group)) +theme_bw() + geom_point(aes(shape=collection), alpha=0.5, size=2) + labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"), y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) + scale_color_manual(values=pal) + ggtitle(paste(title, "_PCA", sep="")) + geom_text(aes(label = sample_name), size=1.75, colour="black")
  print(plot)
  #dev.off()
  #Converting NAs to 0 in donor_age
  vec <- scores$donor_age
  vec[is.na(vec)] <- 0
  scores$donor_age <- vec
  
  if (length(unique(scores$donor_type)) == 3) {
    plot2 <- ggplot(data=scores, mapping = aes(x = PC1, y = PC2, color=donor_type)) + geom_point(size = 2, alpha = 0.5) + labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"), y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) + ggtitle(paste(title, " _PCA - donor type and donor age", sep="")) + geom_text(aes(label=donor_age), size=1.75, colour="black") + scale_color_manual(values=c("grey", pal[1], pal[2])) + theme_bw()
  } else {
    plot2 <- ggplot(data=scores, mapping = aes(x = PC1, y = PC2, color=donor_type)) + geom_point(size = 2, alpha = 0.5) + labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"), y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) + ggtitle(paste(title, " -PCA - donor type and donor age", sep="")) + geom_text(aes(label=donor_age), size=1.75, colour="black") + scale_color_manual(values=c(pal[1], pal[2])) + theme_bw()
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
  # plot <- ggplot(data=scores, mapping = aes(x = PC1, y = PC2, color=eGFR)) +theme_bw() + geom_point(aes(shape=time), alpha=0.5, size=2)+ labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"), y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) + scale_color_manual(values=pal) + ggtitle(paste(title, "_PCA", sep="")) + geom_text(aes(label = sample_name), size=1.75, colour="black")
  # #plot + geom_text(aes(label = sample_name), size=3.5) + xlim(-100, 400)
  # print(plot)
  dev.off()
  return('Clustering done\n')
}

# clustering(pheno_df, "DD_HI_L1", "DD_HI_L2", beta_values_filtered)
# clustering(pheno_df, "DD_HI_L1", "DD_LI_L1", beta_values_filtered)
# clustering(pheno_df, "DD_HI_L1", "LD_LI_L1", beta_values_filtered)
# clustering(pheno_df, "DD_HI_L2", "DD_LI_L2", beta_values_filtered)
# clustering(pheno_df, "DD_HI_L2", "LD_LI_L2", beta_values_filtered)
# clustering(pheno_df, "DD_LI_L1", "DD_LI_L2", beta_values_filtered)
# clustering(pheno_df, "DD_LI_L1", "LD_LI_L1", beta_values_filtered)
# clustering(pheno_df, "DD_LI_L2", "LD_LI_L2", beta_values_filtered)
# clustering(pheno_df, "LD_LI_L1", "LD_LI_L2", beta_values_filtered)

#Copying of the betas dataframe
newBeta_df <- beta_values_filtered

#Adding the collection to the sample_name to introduce consistency 
for (row in 1:nrow(pheno_df)){
  sample <- pheno_df[row, 'sample_name']
  collection <- pheno_df[row, 'collection']
  new_name <- paste(sample, collection, sep = " ")
  pheno_df[row, 'sample_name'] <- new_name
}

#Changing the column names to sample_name
colnames(newBeta_df) <- pheno_df$sample_name

#Calculating average delta beta per comparison
get_deltaBeta <- function(cond1, cond2) {
  pheno_condition1 <- pheno_df[(pheno_df$condition == cond1),]
  pheno_condition2 <- pheno_df[(pheno_df$condition == cond2),]
  betas_condition1 <- newBeta_df[,(colnames(newBeta_df) %in% pheno_condition1$sample_name)]
  betas_condition2 <- newBeta_df[,(colnames(newBeta_df) %in% pheno_condition2$sample_name)]
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

log_df <- data.frame(comparison = character(), number_of_samples = double(), number_of_sig_DMPs = double())

#Loading in relevant annotation file
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

#Addition of condition column 
pheno_df$condition <- 'NA'

#Filling the condition column
for (row in 1:nrow(pheno_df)) {
  if (pheno_df[row, 'donor_type'] == 'DD' & pheno_df[row, 'sample_group'] == 'High' &  pheno_df[row, 'collection'] == 'L1') {
    pheno_df[row, 'condition'] <- "DD_HI_L1"
  } else if (pheno_df[row, 'donor_type'] == 'DD' & pheno_df[row, 'sample_group'] == 'High' &  pheno_df[row, 'collection'] == 'L2') {
    pheno_df[row, 'condition'] <- "DD_HI_L2"
  } else if (pheno_df[row, 'donor_type'] == 'DD' & pheno_df[row, 'sample_group'] == 'Low' &  pheno_df[row, 'collection'] == 'L1') {
    pheno_df[row, 'condition'] <- "DD_LI_L1"
  } else if (pheno_df[row, 'donor_type'] == 'DD' & pheno_df[row, 'sample_group'] == 'Low' &  pheno_df[row, 'collection'] == 'L2') {
    pheno_df[row, 'condition'] <- "DD_LI_L2"
  } else if (pheno_df[row, 'donor_type'] == 'LD' & pheno_df[row, 'sample_group'] == 'Low' &  pheno_df[row, 'collection'] == 'L1') {
    pheno_df[row, 'condition'] <- "LD_LI_L1"
  } else if (pheno_df[row, 'donor_type'] == 'LD' & pheno_df[row, 'sample_group'] == 'Low' &  pheno_df[row, 'collection'] == 'L2') {
    pheno_df[row, 'condition'] <- "LD_LI_L2"
  }
}

#Preparation for model matrix (multiple regression)
condition <- factor(pheno_df$condition)

print('Begin regression model')
# create design matrix
design <- model.matrix(~0+condition+as.numeric(Horvath), data=pheno_df)
#design <- model.matrix(~0+condition, data=pheno_df)

colnames(design) <- c("DD_HI_L1","DD_HI_L2","DD_LI_L1","DD_LI_L2","LD_LI_L1","LD_LI_L2", "Age")
#colnames(design) <- c("DD_HI_L1","DD_HI_L2","DD_LI_L1","DD_LI_L2","LD_LI_L1","LD_LI_L2")

# fit the linear model 
fit1 <- lmFit(m_values, design)
# create a contrast matrix for specific comparisons
contMatrix <- makeContrasts(DD_HI_L1-DD_HI_L2,
                            DD_HI_L1-DD_LI_L1,
                            DD_HI_L1-LD_LI_L1,
                            DD_HI_L2-DD_LI_L2,
                            DD_HI_L2-LD_LI_L2,
                            DD_LI_L1-DD_LI_L2,
                            DD_LI_L1-LD_LI_L1,
                            DD_LI_L2-LD_LI_L2,
                            LD_LI_L1-LD_LI_L2,
                            levels=design)

#contMatrix

# fit the contrasts
fit2 <- contrasts.fit(fit1, contMatrix)
fit2 <- eBayes(fit2)

#topTable(fit2, adjust="BH")

# look at the numbers of DM CpGs at FDR < 0.05
summary(decideTests(fit2))

ann450kSub <- ann450k[match(rownames(m_values),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]

#Identifying and writing output for DMPs
DMPs1 <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub)
DMPs1 <- data.frame(DMPs1)
deltaBeta_df <- get_deltaBeta("DD_HI_L1", "DD_HI_L2")
sample_num <- nrow(subset(pheno_df, (condition == 'DD_HI_L1' | condition == 'DD_HI_L2')))
DMPs1 <- merge(DMPs1, deltaBeta_df, by = 'Name')
output <- "DD_HI_L1-DD_HI_L2_DMPs.csv"
write.csv(DMPs1, file = paste(output_dir, output, sep=""), row.names = FALSE) 
DMPs1_sig <- DMPs1[(DMPs1$adj.P.Val < 0.05),]
print(dim(DMPs1_sig))
output <- "DD_HI_L1-DD_HI_L2_DMPs_sig.csv"
write.csv(DMPs1_sig, file = paste(output_dir, output, sep=""), row.names = FALSE)
log_df[nrow(log_df) + 1,] <- c("DD_HI_L1-DD_HI_L2", sample_num, nrow(DMPs1_sig))

DMPs2 <- topTable(fit2, num=Inf, coef=2, genelist=ann450kSub)
DMPs2 <- data.frame(DMPs2)
deltaBeta_df <- get_deltaBeta("DD_HI_L1", "DD_LI_L1")
sample_num <- nrow(subset(pheno_df, (condition == "DD_HI_L1" | condition == "DD_LI_L1")))
DMPs2 <- merge(DMPs2, deltaBeta_df, by = 'Name')
output <- "DD_HI_L1-DD_LI_L1_DMPs.csv"
write.csv(DMPs2, file = paste(output_dir, output, sep=""), row.names = FALSE) 
DMPs2_sig <- DMPs2[(DMPs2$adj.P.Val < 0.05),]
print(dim(DMPs2_sig))
output <- "DD_HI_L1-DD_LI_L1_DMPs_sig.csv"
write.csv(DMPs2_sig, file = paste(output_dir, output, sep=""), row.names = FALSE)
log_df[nrow(log_df) + 1,] <- c("DD_HI_L1-DD_LI_L1", sample_num, nrow(DMPs2_sig))

DMPs3 <- topTable(fit2, num=Inf, coef=3, genelist=ann450kSub)
DMPs3 <- data.frame(DMPs3)
deltaBeta_df <- get_deltaBeta("DD_HI_L1", "LD_LI_L1")
sample_num <- nrow(subset(pheno_df, (condition == "DD_HI_L1" | condition == "LD_LI_L1")))
DMPs3 <- merge(DMPs3, deltaBeta_df, by = 'Name')
output <- "DD_HI_L1-LD_LI_L1_DMPs.csv"
write.csv(DMPs3, file = paste(output_dir, output, sep=""), row.names = FALSE) 
DMPs3_sig <- DMPs3[(DMPs3$adj.P.Val < 0.05),]
print(dim(DMPs3_sig))
output <- "DD_HI_L1-LD_LI_L1_DMPs_sig.csv"
write.csv(DMPs3_sig, file = paste(output_dir, output, sep=""), row.names = FALSE)
log_df[nrow(log_df) + 1,] <- c("DD_HI_L1-LD_LI_L1", sample_num, nrow(DMPs3_sig))

DMPs4 <- topTable(fit2, num=Inf, coef=4, genelist=ann450kSub)
DMPs4 <- data.frame(DMPs4)
deltaBeta_df <- get_deltaBeta("DD_HI_L2", "DD_LI_L2")
sample_num <- nrow(subset(pheno_df, (condition == "DD_HI_L2" | condition == "DD_LI_L2")))
DMPs4 <- merge(DMPs4, deltaBeta_df, by = 'Name')
output <- "DD_HI_L2-DD_LI_L2_DMPs.csv"
write.csv(DMPs4, file = paste(output_dir, output, sep=""), row.names = FALSE) 
DMPs4_sig <- DMPs4[(DMPs4$adj.P.Val < 0.05),]
print(dim(DMPs4_sig))
output <- "DD_HI_L2-DD_LI_L2_DMPs_sig.csv"
write.csv(DMPs4_sig, file = paste(output_dir, output, sep=""), row.names = FALSE)
log_df[nrow(log_df) + 1,] <- c("DD_HI_L2-DD_LI_L2", sample_num, nrow(DMPs4_sig))

DMPs5 <- topTable(fit2, num=Inf, coef=5, genelist=ann450kSub)
DMPs5 <- data.frame(DMPs5)
deltaBeta_df <- get_deltaBeta("DD_HI_L2", "LD_LI_L2")
sample_num <- nrow(subset(pheno_df, (condition == "DD_HI_L2" | condition == "LD_LI_L2")))
DMPs5 <- merge(DMPs5, deltaBeta_df, by = 'Name')
output <- "DD_HI_L2-LD_LI_L2_DMPs.csv"
write.csv(DMPs5, file = paste(output_dir, output, sep=""), row.names = FALSE) 
DMPs5_sig <- DMPs5[(DMPs5$adj.P.Val < 0.05),]
print(dim(DMPs5_sig))
output <- "DD_HI_L2-LD_LI_L2_DMPs_sig.csv"
write.csv(DMPs5_sig, file = paste(output_dir, output, sep=""), row.names = FALSE)
log_df[nrow(log_df) + 1,] <- c("DD_HI_L2-LD_LI_L2", sample_num, nrow(DMPs5_sig))

DMPs6 <- topTable(fit2, num=Inf, coef=6, genelist=ann450kSub)
DMPs6 <- data.frame(DMPs6)
deltaBeta_df <- get_deltaBeta("DD_LI_L1", "DD_LI_L2")
sample_num <- nrow(subset(pheno_df, (condition == "DD_LI_L1" | condition == "DD_LI_L2")))
DMPs6 <- merge(DMPs6, deltaBeta_df, by = 'Name')
output <- "DD_LI_L1-DD_LI_L2_DMPs.csv"
write.csv(DMPs6, file = paste(output_dir, output, sep=""), row.names = FALSE) 
DMPs6_sig <- DMPs6[(DMPs6$adj.P.Val < 0.05),]
print(dim(DMPs6_sig))
output <- "DD_LI_L1-DD_LI_L2_DMPs_sig.csv"
write.csv(DMPs6_sig, file = paste(output_dir, output, sep=""), row.names = FALSE)
log_df[nrow(log_df) + 1,] <- c("DD_LI_L1-DD_LI_L2", sample_num, nrow(DMPs6_sig))

DMPs7 <- topTable(fit2, num=Inf, coef=7, genelist=ann450kSub)
DMPs7 <- data.frame(DMPs7)
deltaBeta_df <- get_deltaBeta("DD_LI_L1", "LD_LI_L1")
sample_num <- nrow(subset(pheno_df, (condition == "DD_LI_L1" | condition == "LD_LI_L1")))
DMPs7 <- merge(DMPs7, deltaBeta_df, by = 'Name')
output <- "DD_LI_L1-LD_LI_L1_DMPs.csv"
write.csv(DMPs7, file = paste(output_dir, output, sep=""), row.names = FALSE) 
DMPs7_sig <- DMPs7[(DMPs7$adj.P.Val < 0.05),]
print(dim(DMPs7_sig))
output <- "DD_LI_L1-LD_LI_L1_DMPs_sig.csv"
write.csv(DMPs7_sig, file = paste(output_dir, output, sep=""), row.names = FALSE)
log_df[nrow(log_df) + 1,] <- c("DD_LI_L1-LD_LI_L1", sample_num, nrow(DMPs7_sig))

DMPs8 <- topTable(fit2, num=Inf, coef=8, genelist=ann450kSub)
DMPs8 <- data.frame(DMPs8)
deltaBeta_df <- get_deltaBeta("DD_LI_L2", "LD_LI_L2")
sample_num <- nrow(subset(pheno_df, (condition == "DD_LI_L2" | condition == "LD_LI_L2")))
DMPs8 <- merge(DMPs8, deltaBeta_df, by = 'Name')
output <- "DD_LI_L2-LD_LI_L2_DMPs.csv"
write.csv(DMPs8, file = paste(output_dir, output, sep=""), row.names = FALSE) 
DMPs8_sig <- DMPs8[(DMPs8$adj.P.Val < 0.05),]
print(dim(DMPs8_sig))
output <- "DD_LI_L2-LD_LI_L2_DMPs_sig.csv"
write.csv(DMPs8_sig, file = paste(output_dir, output, sep=""), row.names = FALSE)
log_df[nrow(log_df) + 1,] <- c("DD_LI_L2-LD_LI_L2", sample_num, nrow(DMPs8_sig))

DMPs9<- topTable(fit2, num=Inf, coef=9, genelist=ann450kSub)
DMPs9 <- data.frame(DMPs9)
deltaBeta_df <- get_deltaBeta("LD_LI_L1", "LD_LI_L2")
sample_num <- nrow(subset(pheno_df, (condition == "LD_LI_L1" | condition == "LD_LI_L2")))
DMPs9 <- merge(DMPs9, deltaBeta_df, by = 'Name')
output <- "LD_LI_L1-LD_LI_L2_DMPs.csv"
write.csv(DMPs9, file = paste(output_dir, output, sep=""), row.names = FALSE) 
DMPs9_sig <- DMPs9[(DMPs9$adj.P.Val < 0.05),]
print(dim(DMPs9_sig))
output <- "LD_LI_L1-LD_LI_L2_DMPs_sig.csv"
write.csv(DMPs9_sig, file = paste(output_dir, output, sep=""), row.names = FALSE)
log_df[nrow(log_df) + 1,] <- c("LD_LI_L1-LD_LI_L2", sample_num, nrow(DMPs9_sig))

write.csv(log_df, file = paste(output_dir, "combined_array_liver_log.csv", sep=""), row.names = FALSE)

library(qqman)
library(DMRcate)
generate_man <- function(DMPs, comp) {
  #Manhattan plot using the DMPs
  cat("Generating manhattan plot from DMPs\n")
  title <- paste(comp, " (Adj. P-val)", sep="")
  col=c("black","grey")
  DMPs$chr = str_replace_all(DMPs$chr, 'chr', '')
  DMPs$chr = as.numeric(DMPs$chr)
  DMPs$pos = as.numeric(DMPs$pos)
  output <- paste(comp, "_manhattan.jpeg", sep="")
  jpeg(paste(output_dir, output, sep=""), quality = 100)
  manhattan(DMPs, chr="chr", bp="pos", p="adj.P.Val", snp="Islands_Name", col=col, suggestiveline=(-log10(0.05)), main=title)
  dev.off()
  return("Done with manhattan plot")
}

generate_man(DMPs1, 'DD_HI_L1-DD_HI_L2')
generate_man(DMPs2, 'DD_HI_L1-DD_LI_L1')
generate_man(DMPs3, 'DD_HI_L1-LD_LI_L1')
generate_man(DMPs4, 'DD_HI_L2-DD_LI_L2')
generate_man(DMPs5, 'DD_HI_L2-LD_LI_L2')
generate_man(DMPs6, 'DD_LI_L1-DD_LI_L2')
generate_man(DMPs7, 'DD_LI_L1-LD_LI_L1')
generate_man(DMPs8, 'DD_LI_L2-LD_LI_L2')
generate_man(DMPs9, 'LD_LI_L1-LD_LI_L2')

library(DMRcate)
#Specify groups for DMR.plot 
groups <- pal[1:length(unique(pheno_df$condition))]
names(groups) <- levels(factor(pheno_df$condition))
cols <- groups[as.character(factor(pheno_df$condition))]
samps <- 1:nrow(pheno_df)

#1 DMRs - DD_HI_L1-DD_HI_L2
myAnnotation <- cpg.annotate(datatype = "array", object = as.matrix(m_values), what = "M",
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contMatrix, 
                             coef = "DD_HI_L1 - DD_HI_L2", arraytype = "450K")
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs, genome = "hg19")
write.csv(results.ranges, file=paste(output_dir, "DD_HI_L1-DD_HI_L2_DMRs.csv", sep=""), row.names=FALSE)

# draw the plot for the top DMR
jpeg(file=paste(output_dir, "DD_HI_L1-DD_HI_L2_DMR.jpeg", sep=""), quality = 100)
DMR.plot(ranges = results.ranges, dmr=1, CpGs=as.matrix(m_values), phen.col=cols, 
         what="Beta", arraytype="450K", genome="hg19", cex=0.5, pch=16, toscale=TRUE, plotmedians=TRUE, samps=samps)
dev.off

# #2 DMRs - DD_HI_L1-DD_LI_L1
# myAnnotation <- cpg.annotate(datatype = "array", object = as.matrix(m_values), what = "M",
#                              analysis.type = "differential", design = design, 
#                              contrasts = TRUE, cont.matrix = contMatrix, 
#                              coef = "DD_HI_L1 - DD_LI_L1", arraytype = "450K")
# DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
# results.ranges <- extractRanges(DMRs, genome = "hg19")
# write.csv(results.ranges, file=paste(output_dir, "DD_HI_L1-DD_LI_L1_DMRs.csv", sep=""), row.names=FALSE)
# 
# # draw the plot for the top DMR
# jpeg(file=paste(output_dir, "DD_HI_L1-DD_LI_L1_DMR.jpeg", sep=""), quality = 100)
# DMR.plot(ranges = results.ranges, dmr=1, CpGs=as.matrix(m_values), phen.col=cols, 
#          what="Beta", arraytype="450K", genome="hg19", cex=0.5, pch=16, toscale=TRUE, plotmedians=TRUE, samps=samps)
# dev.off

#3 DMRs - DD_HI_L1-LD_LI_L1
# myAnnotation <- cpg.annotate(datatype = "array", object = as.matrix(m_values), what = "M",
#                              analysis.type = "differential", design = design, 
#                              contrasts = TRUE, cont.matrix = contMatrix, 
#                              coef = "DD_HI_L1 - LD_LI_L1", arraytype = "450K")
# DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
# results.ranges <- extractRanges(DMRs, genome = "hg19")
# write.csv(results.ranges, file=paste(output_dir, "DD_HI_L1-LD_LI_L1_DMRs.csv", sep=""), row.names=FALSE)
# 
# # draw the plot for the top DMR
# jpeg(file=paste(output_dir, "DD_HI_L1-LD_LI_L1_DMR.jpeg", sep=""), quality = 100)
# DMR.plot(ranges = results.ranges, dmr=1, CpGs=as.matrix(m_values), phen.col=cols, 
#          what="Beta", arraytype="450K", genome="hg19", cex=0.5, pch=16, toscale=TRUE, plotmedians=TRUE, samps=samps)
# dev.off

#4 DMRs - DD_HI_L2-DD_LI_L2 - No sig. DMPs
# myAnnotation <- cpg.annotate(datatype = "array", object = as.matrix(m_values), what = "M",
#                              analysis.type = "differential", design = design, 
#                              contrasts = TRUE, cont.matrix = contMatrix, 
#                              coef = "DD_HI_L2 - DD_LI_L2", arraytype = "450K")
# DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
# results.ranges <- extractRanges(DMRs, genome = "hg19")
# write.csv(results.ranges, file=paste(output_dir, "DD_HI_L2-DD_LI_L2_DMRs.csv", sep=""), row.names=FALSE)

#5 DMRs - DD_HI_L2-LD_LI_L2
myAnnotation <- cpg.annotate(datatype = "array", object = as.matrix(m_values), what = "M",
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contMatrix, 
                             coef = "DD_HI_L2 - LD_LI_L2", arraytype = "450K")
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs, genome = "hg19")
write.csv(results.ranges, file=paste(output_dir, "DD_HI_L2-LD_LI_L2_DMRs.csv", sep=""), row.names=FALSE)

# draw the plot for the top DMR
jpeg(file=paste(output_dir, "DD_HI_L2-LD_LI_L2_DMR.jpeg", sep=""), quality = 100)
DMR.plot(ranges = results.ranges, dmr=1, CpGs=as.matrix(m_values), phen.col=cols, 
         what="Beta", arraytype="450K", genome="hg19", cex=0.5, pch=16, toscale=TRUE, plotmedians=TRUE, samps=samps)
dev.off

#6 DMRs - DD_LI_L1-DD_LI_L2
myAnnotation <- cpg.annotate(datatype = "array", object = as.matrix(m_values), what = "M",
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contMatrix, 
                             coef = "DD_LI_L1 - DD_LI_L2", arraytype = "450K")
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs, genome = "hg19")
write.csv(results.ranges, file=paste(output_dir, "DD_LI_L1-DD_LI_L2_DMRs.csv", sep=""), row.names=FALSE)

# draw the plot for the top DMR
jpeg(file=paste(output_dir, "DD_LI_L1-DD_LI_L2_DMR.jpeg", sep=""), quality = 100)
DMR.plot(ranges = results.ranges, dmr=1, CpGs=as.matrix(m_values), phen.col=cols, 
         what="Beta", arraytype="450K", genome="hg19", cex=0.5, pch=16, toscale=TRUE, plotmedians=TRUE, samps=samps)
dev.off

#7 DMRs - DD_LI_L1-LD_LI_L1
myAnnotation <- cpg.annotate(datatype = "array", object = as.matrix(m_values), what = "M",
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contMatrix, 
                             coef = "DD_LI_L1 - LD_LI_L1", arraytype = "450K")
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs, genome = "hg19")
write.csv(results.ranges, file=paste(output_dir, "DD_LI_L1-LD_LI_L1_DMRs.csv", sep=""), row.names=FALSE)

# draw the plot for the top DMR
jpeg(file=paste(output_dir, "DD_LI_L1-LD_LI_L1_DMR.jpeg", sep=""), quality = 100)
DMR.plot(ranges = results.ranges, dmr=1, CpGs=as.matrix(m_values), phen.col=cols, 
         what="Beta", arraytype="450K", genome="hg19", cex=0.5, pch=16, toscale=TRUE, plotmedians=TRUE, samps=samps)
dev.off

#8 DMRs - DD_LI_L2-LD_LI_L2
myAnnotation <- cpg.annotate(datatype = "array", object = as.matrix(m_values), what = "M",
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contMatrix, 
                             coef = "DD_LI_L2 - LD_LI_L2", arraytype = "450K")
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs, genome = "hg19")
write.csv(results.ranges, file=paste(output_dir, "DD_LI_L2-LD_LI_L2_DMRs.csv", sep=""), row.names=FALSE)

# draw the plot for the top DMR
jpeg(file=paste(output_dir, "DD_LI_L2-LD_LI_L2_DMR.jpeg", sep=""), quality = 100)
DMR.plot(ranges = results.ranges, dmr=1, CpGs=as.matrix(m_values), phen.col=cols, 
         what="Beta", arraytype="450K", genome="hg19", cex=0.5, pch=16, toscale=TRUE, plotmedians=TRUE, samps=samps)
dev.off

#9 DMRs - LD_LI_L1-LD_LI_L2 - No sig. DMPs
# myAnnotation <- cpg.annotate(datatype = "array", object = as.matrix(m_values), what = "M",
#                              analysis.type = "differential", design = design, 
#                              contrasts = TRUE, cont.matrix = contMatrix, 
#                              coef = "LD_LI_L1 - LD_LI_L2", arraytype = "450K")
# DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
# results.ranges <- extractRanges(DMRs, genome = "hg19")
# write.csv(results.ranges, file=paste(output_dir, "LD_LI_L1-LD_LI_L2_DMRs.csv", sep=""), row.names=FALSE)

q()

