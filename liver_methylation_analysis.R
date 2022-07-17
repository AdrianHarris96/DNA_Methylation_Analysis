#Script for DNA Methylation Analysis
library(minfi)
library(rio)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(dendextend)
library(limma)

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

#Simply ensure the correct order 
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
  jpeg(paste(output_dir, "p_values450k.jpeg", sep=""), quality = 90)
  barplot((subset(detP_df, array_type == '450K'))$p_value_mean, col=pal[factor(detP_df$collection)], names.arg=(subset(detP_df, array_type == '450K'))$sample_name, las=2, cex.names=0.4, cex.axis=0.5, space=0.5, ylab="Mean detection p-values", main='450K Array')
  legend("topleft", legend=levels(factor(detP_df$collection)), fill=pal,
         cex=0.27, bty = "n", bg="white")
  dev.off()
  
  #Plotting EPIC barplot
  jpeg(paste(output_dir, "p_valuesEPIC.jpeg", sep=""), quality = 90)
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
  jpeg(paste(output_dir, "preprocessQC.jpeg", sep=""), quality = 90)
  plot <- ggplot(data=qc, mapping = aes(x = mMed, y = uMed, color=collection)) + geom_point(aes(shape=array_type), alpha=0.5) + xlim(7, 14) + ylim(7, 14) + theme_bw()+ geom_abline(slope=-1, intercept = 21.25 , color="black", linetype="dashed", size=0.5) + scale_color_manual(values=pal)
  print(plot)
  dev.off()
  jpeg(paste(output_dir, "preprocessDensity.jpeg", sep=""), quality = 90)
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
  jpeg(paste(output_dir, "postNormQC.jpeg", sep=""), quality = 90)
  plot2 <- ggplot(data=postqc, mapping = aes(x = mMed, y = uMed, color=collection)) + geom_point(aes(shape=array_type), alpha=0.5) + xlim(5, 14) + ylim(5, 14) + theme_bw()+ geom_abline(slope=-1, intercept = 21.25 , color="black", linetype="dashed", size=0.5) + scale_color_manual(values=pal)
  print(plot2)
  dev.off()
  jpeg(paste(output_dir, "postNormDensity.jpeg", sep=""), quality = 90)
  densityPlot(mtSet, sampGroups = postqc$array_type)
  dev.off()
  rm(postqc, plot2)
}

#Load betas or map to genome-generate betas
if (file.exists(paste(output_dir, "beta_values.csv", sep=""))) {
  cat('Loading beta_values\n')
  beta_values_filtered <- import(paste(output_dir, "beta_values.csv", sep=""))
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
  
  keep <- !(featureNames(gmtSet) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")]) #remove probes that are not of the chrom x or y
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
  cat('Writing beta to CSV')
  write.csv(beta_values_filtered, file = paste(output_dir, "beta_values.csv", sep=""), row.names = TRUE)
}

#Loading lumi and converting beta2m 
library(lumi)
m_values <- beta2m(beta_values_filtered)

if (file.exists(paste(output_dir, "m_values.csv", sep=""))) {
  cat('Skip writing of m_values\n')
} else {
  write.csv(beta_values_filtered, file = paste(output_dir, "m_values.csv", sep=""), row.names = TRUE)
}

#Drop 'methylated' and 'unmethylated' sample names
pheno_df <- pheno_df[!(pheno_df$Basename %in% c('200999740023_R05C02', '200999740023_R06C02', '201004900096_R05C02', '201004900096_R06C02', '202702240054_R01C01', '202702240054_R02C01', '202702240079_R07C01', '202702240079_R08C01', '3999442124_R05C02', '3999442124_R06C02')),]

#Drop the one unpaired samples
pheno_df <- pheno_df[!(pheno_df$sample_name == 'V037L1'),]

#Drop DCD samples
pheno_df = pheno_df[!(pheno_df$donor_type == 'DCD'),]

#Exclude control and DCD samples for liver data 
beta_values_filtered <- beta_values_filtered[,(colnames(beta_values_filtered) %in% pheno_df$Basename)]
m_values <- m_values[, colnames(m_values) %in% colnames(beta_values_filtered)]

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

clustering <- function(pheno, condition1, condition2, betas) {
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
  plot <- ggplot(data=scores, mapping = aes(x = PC1, y = PC2, color=sample_group)) +theme_bw() + geom_point(aes(shape=collection), alpha=0.5, size=2) + labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"), y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) + scale_color_manual(values=pal) + ggtitle(paste(title, "_PCA", sep="")) + geom_text(aes(label = sample_id), size=1.75, colour="black")
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
  # plot <- ggplot(data=scores, mapping = aes(x = PC1, y = PC2, color=eGFR)) +theme_bw() + geom_point(aes(shape=time), alpha=0.5, size=2)+ labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"), y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) + scale_color_manual(values=pal) + ggtitle(paste(title, "_PCA", sep="")) + geom_text(aes(label = sample_id), size=1.75, colour="black")
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

#Vector of comparisons
comparisons <- c('DD_HI_L1-DD_HI_L2', 'DD_HI_L1-DD_LI_L1', 'DD_HI_L1-LD_LI_L1', 'DD_HI_L2-DD_LI_L2', 'DD_HI_L2-LD_LI_L2', 'DD_LI_L1-DD_LI_L2', 'DD_LI_L1-LD_LI_L1', 'DD_LI_L2-LD_LI_L2', 'LD_LI_L1-LD_LI_L2')

#Loading in relevant annotation file
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

for (comp in comparisons) {
  cat("Identify CpGs\n")
  cond <- unlist(strsplit(comp, "-"))
  condition1 <- cond[1]
  condition2 <- cond[2]
  
  if (condition1 == "DD_HI_L1") {
    pheno1 <- pheno_df[(pheno_df$donor_type == 'DD' & pheno_df$sample_group == 'High' & pheno_df$collection == 'L1'),]
  } else if (condition1 == "DD_HI_L2"){
    pheno1 <- pheno_df[(pheno_df$donor_type == 'DD' & pheno_df$sample_group == 'High' & pheno_df$collection == 'L2'),]
  } else if (condition1 == "DD_LI_L1") {
    pheno1 <- pheno_df[(pheno_df$donor_type == 'DD' & pheno_df$sample_group == 'Low' & pheno_df$collection == 'L1'),]
  } else if (condition1 == "DD_LI_L2") {
    pheno1 <- pheno_df[(pheno_df$donor_type == 'DD' & pheno_df$sample_group == 'Low' & pheno_df$collection == 'L2'),]
  } else if (condition1 == "LD_LI_L1") {
    pheno1 <- pheno_df[(pheno_df$donor_type == 'LD' & pheno_df$sample_group == 'Low' & pheno_df$collection == 'L1'),]
  } else if (condition1 == "LD_LI_L2") {
    pheno1 <- pheno_df[(pheno_df$donor_type == 'LD' & pheno_df$sample_group == 'Low' & pheno_df$collection == 'L2'),]
  } else {
    cat('Comparison does not exist\n')
  }
  
  if (condition2 == "DD_HI_L1") {
    pheno2 <- pheno_df[(pheno_df$donor_type == 'DD' & pheno_df$sample_group == 'High' & pheno_df$collection == 'L1'),]
  } else if (condition2 == "DD_HI_L2"){
    pheno2 <- pheno_df[(pheno_df$donor_type == 'DD' & pheno_df$sample_group == 'High' & pheno_df$collection == 'L2'),]
  } else if (condition2 == "DD_LI_L1") {
    pheno2 <- pheno_df[(pheno_df$donor_type == 'DD' & pheno_df$sample_group == 'Low' & pheno_df$collection == 'L1'),]
  } else if (condition2 == "DD_LI_L2") {
    pheno2 <- pheno_df[(pheno_df$donor_type == 'DD' & pheno_df$sample_group == 'Low' & pheno_df$collection == 'L2'),]
  } else if (condition2 == "LD_LI_L1") {
    pheno2 <- pheno_df[(pheno_df$donor_type == 'LD' & pheno_df$sample_group == 'Low' & pheno_df$collection == 'L1'),]
  } else if (condition2 == "LD_LI_L2") {
    pheno2 <- pheno_df[(pheno_df$donor_type == 'LD' & pheno_df$sample_group == 'Low' & pheno_df$collection == 'L2'),]
  } else {
    cat('Comparison does not exist\n')
  }
  
  pheno <- rbind(pheno1, pheno2)
  #print(dim(pheno))
  
  m_values_condition <- m_values[,(colnames(m_values) %in% pheno$Basename)]
  m_values_condition <- as.matrix(m_values_condition)
  
  if (comp == 'DD_HI_L1-DD_HI_L2') {
    condition <- pheno$collection
  } else if (comp == 'DD_HI_L1-DD_LI_L1') {
    condition <- pheno$sample_group
  } else if (comp == 'DD_HI_L1-LD_LI_L1') {
    condition <- pheno$sample_group
  } else if (comp == 'DD_HI_L2-DD_LI_L2') {
    condition <- pheno$sample_group
  } else if (comp == 'DD_HI_L2-LD_LI_L2') {
    condition <- pheno$sample_group
  } else if (comp == 'DD_LI_L1-DD_LI_L2') {
    condition <- pheno$collection
  } else if (comp == 'DD_LI_L1-LD_LI_L1') {
    condition <- pheno$donor_type
  } else if (comp == 'DD_LI_L2-LD_LI_L2') {
    condition <- pheno$donor_type
  } else if (comp == 'LD_LI_L1-LD_LI_L2') {
    condition <- pheno$collection
  }
  
  DMPs <- dmpFinder(m_values_condition, pheno=condition, type = "categorical", shrinkVar = TRUE)
  DMPs$adj_p <- p.adjust(DMPs$pval, method="BY")
  print(comp)
  DMPs_sig <- DMPs[(DMPs$pval < 0.05),]
  print(dim(DMPs_sig))
  DMPs_sig <- DMPs[(DMPs$adj_p < 0.05),]
  print(dim(DMPs_sig))

  ann450kSub <- ann450k[match(rownames(m_values_condition),ann450k$Name), c(1:4,12:19,24:ncol(ann450k))]
  ann450kSub <- data.frame(ann450kSub)
  
  DMPs <- merge(DMPs, ann450kSub, by = "row.names")
  
  DMPs$probes <- DMPs$Row.names
  DMPs <- DMPs[,c(ncol(DMPs), 2:(ncol(DMPs)-1))]
  
  #output <- paste(comp, "_DMPs.csv", sep="")
  #write.csv(DMPs, file = paste(output_dir, output, sep=""), row.names = FALSE)
  
  #Manhattan plot using the DMPs
  cat("Generating manhattan plot from DMPs\n")
  library(qqman)
  library(DMRcate)
  title <- paste(comp, " (Adj. P-val)", sep="")
  col=c("black","grey")
  DMPs$chr = str_replace_all(DMPs$chr, 'chr', '')
  DMPs$chr = as.numeric(DMPs$chr)
  DMPs$pos = as.numeric(DMPs$pos)
  output <- paste(comp, "_manhattan.jpeg", sep="")
  jpeg(paste(output_dir, output, sep=""), quality = 90)
  manhattan(DMPs, chr="chr", bp="pos",, p="adj_p", snp="Islands_Name", col=col, suggestiveline=(-log10(0.05)), main=title)
  dev.off()
}

q()

