#Script for DNA Methylation Analysis
library(minfi)
library(rio)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

#Example input: Rscript liver_methylation_analysis.R <pheno_file> <base_dir> <output_dir>
args=commandArgs(trailingOnly=TRUE)
pheno_file = args[1]
base_dir = args[2]
git_dir = args[3]
output_dir = args[4]
comparison = args[5]

#Local Machine
# pheno_file = '/Users/adrianharris/Documents/dna_methylation_analysis/liver_sample_sheet.csv'
# base_dir = '/Users/adrianharris/Desktop/liver/'
# git_dir = '/Users/adrianharris/Documents/dna_methylation_analysis/'
# output_dir = '/Users/adrianharris/Desktop/liver/'

if (file.exists(output_dir)) {
  cat("Directory already exists\n")
} else {
  dir.create(output_dir)
}

#Importing manually-curated sample sheet
pheno_df <- import(pheno_file)

#Must remove outlier sample, 203504430032_R01C01 (and its paired sample 203504430032-R02C01)
pheno_df <- pheno_df[!(pheno_df$sample_group == 'control' | pheno_df$sample_group == 'Control'),]
pheno_df <- pheno_df[!(pheno_df$Basename == '203751390020_R08C01' | pheno_df$Basename == '203751390020_R01C01'),]

nrow(subset(pheno_df, array_type == '450K'))
nrow(subset(pheno_df, array_type == 'EPIC'))

#Split into phenotype file 
pheno450k <- pheno_df[!(pheno_df$array_type == 'EPIC'),]
phenoEPIC <- pheno_df[!(pheno_df$array_type == '450K'),]

#Specify the directories and read in respective IDAT Files 
dir450k <- paste(base_dir, "450k_array", sep="")
dirEPIC <- paste(base_dir, "EPIC_array", sep="")

if (file.exists(paste(output_dir, "rgSet.RDS", sep=""))) {
  cat('Loading in rgSet (combined)\n')
  rgSet <- readRDS(paste(output_dir, "rgSet.RDS", sep=""))
} else {
  cat('Generate rgSet\n')
  rgSet450k <- read.metharray.exp(base=dir450k, target=pheno450k)
  rgSetEPIC <- read.metharray.exp(base=dirEPIC, target=phenoEPIC, force=TRUE)
  rgSet <- combineArrays(rgSet450k, rgSetEPIC)
  #Saving set as an RDS file 
  saveRDS(rgSet450k, file = paste(output_dir, "rgSet450k.RDS", sep=""))
  saveRDS(rgSetEPIC, file = paste(output_dir, "rgSetEPIC.RDS", sep=""))
  saveRDS(rgSet, file = paste(output_dir, "rgSet.RDS", sep=""))
  rm(rgSet450k)
  rm(rgSetEPIC)
}

#Color Scheme Defined 
pal <- brewer.pal(4,"Dark2")

#Calculate Detection p-values if 
if (file.exists(paste(output_dir, "p-values.csv", sep=""))) {
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
  
  #Plotting 450K barplot 
  par(mfrow=c(2,1))
  jpeg(paste(output_dir, "p_values.jpeg", sep=""), quality = 90)
  barplot((subset(detP_df, array_type == '450K'))$p_value_mean, col=pal[factor(detP_df$collection)], names.arg=(subset(detP_df, array_type == '450K'))$sample_name, las=2, cex.names=0.4, cex.axis=0.5, space=0.5, ylab="Mean detection p-values", main='450K Array')
  legend("topleft", legend=levels(factor(detP_df$collection)), fill=pal,
         cex=0.27, bty = "n", bg="white")
  
  #Plotting EPIC barplot
  barplot((subset(detP_df, array_type == 'EPIC'))$p_value_mean, col=pal[factor(detP_df$collection)], names.arg=(subset(detP_df, array_type == 'EPIC'))$sample_name, las=2, cex.names=0.4, cex.axis=0.5, space=0.5, ylab="Mean detection p-values", main='EPIC Array')
  legend("topleft", legend=levels(factor(detP_df$collection)), fill=pal,
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
  plot <- ggplot(data=qc, mapping = aes(x = mMed, y = uMed, color=collection)) + geom_point(aes(shape=array_type), alpha=0.5) + xlim(7, 14) + ylim(7, 14) + theme_bw()+ geom_abline(slope=-1, intercept = 21.25 , color="black", linetype="dashed", size=0.5) + scale_color_manual(values=pal)
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
  plot2 <- ggplot(data=postqc, mapping = aes(x = mMed, y = uMed, color=collection)) + geom_point(aes(shape=array_type), alpha=0.5) + xlim(5, 14) + ylim(5, 14) + theme_bw()+ geom_abline(slope=-1, intercept = 21.25 , color="black", linetype="dashed", size=0.5) + scale_color_manual(values=pal)
  print(plot2)
  dev.off()
  jpeg(paste(output_dir, "postNormDensity.jpeg", sep=""), quality = 90)
  densityPlot(mtSet, sampGroups = postqc$array_type)
  dev.off()
  rm(postqc, plot2)
}

# Map to Genome
if (file.exists(paste(output_dir, "beta_values.csv", sep=""))) {
  cat('Loading beta_values')
  beta_values_filtered <- import(paste(output_dir, "beta_values.csv", sep=""))
} else {
  cat('Converting to Genomic Methyl Set\n')
  rSet <- ratioConvert(mtSet, what = "both", keepCN = TRUE)
  gmtSet <- mapToGenome(rSet)
  print(dim(gmtSet))
  
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
  print(dim(gmtSet))
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
  # m_values <- getM(gmtSet)
  
  #Remove probes Hypomethylated in all samples identified by CpG sites having beta < 0.05 in all samples
  beta_values_filtered <- as.data.frame(beta_values) 
  cat('Hypomethylated\n')
  print(dim(filter_all(beta_values_filtered, all_vars(. < 0.05)))) #Number of hypomethylated
  beta_values_filtered <- filter_all(beta_values_filtered, any_vars(. >= 0.05)) 
  
  
  #Remove probes Hypermethylated in all samples identified by CpG sites having beta > 0.95 in all samples
  cat("Hypermethylated\n")
  print(dim(filter_all(beta_values_filtered, all_vars(. > 0.95))))
  #3091 hypermethylated
  beta_values_filtered <- filter_all(beta_values_filtered, any_vars(. < 0.95)) 
  cat("Final Dimensions\n")
  print(dim(beta_values_filtered))
  write.csv(beta_values_filtered, file = paste(output_dir, "beta_values.csv", sep=""), row.names = TRUE)
}

generate_dendro <- function(beta, pheno, timepoint){
  beta_t <- data.frame(t(beta))
  row.names(pheno) <- pheno$Basename
  if (timepoint == 'L1') {
    pheno <- pheno[(pheno$collection == 'L1'),]
  } else if (timepoint == 'L2') {
    pheno <- pheno[(pheno$collection == 'L2'),]
  } else {
    cat("Skip filtering down")
  }
  
  pheno <- pheno %>% select(c('sample_id', 'donor_age', 'donor_type', 'sample_group'))
  
  final_beta <- merge(pheno, beta_t, by='row.names')
  
  #changing the row.names accordingly 
  final_beta['sample_id_age'] <- 'na'
  final_beta <- final_beta[,c(ncol(final_beta),1:(ncol(final_beta)-1))]
  final_beta['sample_id_donor'] <- 'na'
  final_beta <- final_beta[,c(ncol(final_beta),1:(ncol(final_beta)-1))]
  final_beta['sample_id_group'] <- 'na'
  final_beta <- final_beta[,c(ncol(final_beta),1:(ncol(final_beta)-1))]
  
  for (row in 1:nrow(final_beta)) {
    age <- paste(final_beta[row, 'sample_id'], final_beta[row, 'donor_age'], sep = " ")
    donor <- paste(final_beta[row, 'sample_id'], final_beta[row, 'donor_type'], sep = " ")
    if (final_beta[row, 'sample_group'] == 'High Injury') {
      injury <- 'High'
    } else {
      injury <- 'Low'
    }
    group <- paste(final_beta[row, 'sample_id'], injury, sep = " ")
    final_beta[row, 'sample_id_age'] <- age
    final_beta[row, 'sample_id_donor'] <- donor
    final_beta[row, 'sample_id_group'] <- group
  }
  
  dendro_out <- paste(timepoint, 'dendrograms.pdf', sep=" ")
  pdf(file = paste(output_dir, dendro_out, sep=""), width = 12, height = 8)
  row.names(final_beta) <- final_beta$sample_id_age
  clusters <- hclust(dist(final_beta[, 9:ncol(final_beta)]))
  dend <- as.dendrogram(clusters)
  if (timepoint == 'L1-L2'){
    dend <- set(dend, "labels_cex", 0.4)
  }
  plot(dend, xlab = "Sample ID and Donor Age", ylab="Height", main= paste(timepoint, "- age Dendrogram", sep = " "))
  
  row.names(final_beta) <- final_beta$sample_id_donor
  clusters <- hclust(dist(final_beta[, 9:ncol(final_beta)]))
  dend <- as.dendrogram(clusters)
  if (timepoint == 'L1-L2'){
    dend <- set(dend, "labels_cex", 0.4)
  }
  plot(dend, xlab = "Sample ID and Donor Status", ylab="Height", main= paste(timepoint, "- Donor Status Dendrogram", sep = " "))
  
  row.names(final_beta) <- final_beta$sample_id_group
  clusters <- hclust(dist(final_beta[, 9:ncol(final_beta)]))
  dend <- as.dendrogram(clusters)
  if (timepoint == 'L1-L2'){
    dend <- set(dend, "labels_cex", 0.4)
  }
  plot(dend, xlab = "Sample ID and Injury Status", ylab="Height", main= paste(timepoint, "- Injury Status Dendrogram", sep = " "))
  
  dev.off()
}

timeList <- c('L1', 'L2', 'L1-L2')
for (time in timeList) {
  generate_dendro(beta_values_filtered, pheno_df, time)
}

q()

if (comparison == 'K1_Low_K1_High') {
  pheno_df <- pheno_df[((pheno_df$time == 'K1' & pheno_df$eGFR == 'Low') | (pheno_df$time == 'K1' & pheno_df$eGFR == 'High')),]
} else if (comparison == 'K2_Low_K2_High') {
  pheno_df <- pheno_df[((pheno_df$time == 'K2' & pheno_df$eGFR == 'Low') | (pheno_df$time == 'K2' & pheno_df$eGFR == 'High')),]
} else if (comparison == 'K1_High_K2_High') {
  pheno_df <- pheno_df[((pheno_df$time == 'K1' & pheno_df$eGFR == 'High') | (pheno_df$time == 'K2' & pheno_df$eGFR == 'High')),]
} else if (comparison == 'K1_Low_K2_Low') {
  pheno_df <- pheno_df[((pheno_df$time == 'K1' & pheno_df$eGFR == 'Low') | (pheno_df$time == 'K2' & pheno_df$eGFR == 'Low')),]
} else if (comparison == 'K1_High_K2_Low') {
  pheno_df <- pheno_df[((pheno_df$time == 'K1' & pheno_df$eGFR == 'High') | (pheno_df$time == 'K2' & pheno_df$eGFR == 'Low')),]
} else if (comparison == 'K1_Low_K2_High') {
  pheno_df <- pheno_df[((pheno_df$time == 'K1' & pheno_df$eGFR == 'Low') | (pheno_df$time == 'K2' & pheno_df$eGFR == 'High')),]
} else {
  cat('Comparison request does not exist\n')
  q()
}

library(lumi)
m_values = beta2m(beta_values_filtered)
if (file.exists(paste(output_dir, "beta_values.csv", sep="")) & file.exists(paste(output_dir, "m_values.csv", sep=""))) {
  cat('Skip beta and m_value CSV\n')
} else {
  write.csv(beta_values_filtered, file = paste(output_dir, "beta_values.csv", sep=""), row.names = TRUE)
  write.csv(m_values, file = paste(output_dir, "m_values.csv", sep=""), row.names = TRUE)
}

#Filter beta dataframe using column names for the relevant comparison
beta_values_filtered <- beta_values_filtered[,(colnames(beta_values_filtered) %in% pheno_df$Basename)]
dim(beta_values_filtered)

# #EXCLUDE CONTROL SAMPLES AND DCD SAMPLES - Liver dataset
# beta_values_case = beta_values_filtered[,!(colnames(beta_values_filtered) %in% c("200999740023_R05C02","200999740023_R06C02","201004900096_R05C02","201004900096_R06C02","202702240054_R01C01","202702240054_R02C01","202702240079_R07C01","202702240079_R08C01","3999442124_R05C02","3999442124_R06C02","203751390020_R02C01","3999442124_R01C02","201004900096_R02C02","200999740005_R06C02","201004900018_R06C01","203751390017_R07C01","200687170042_R05C02","201004900096_R03C02","200999740023_R01C01","201004900018_R01C02"))]
# pheno_df_case = pheno_df[!(rownames(pheno_df) %in% c("200999740023_R05C02","200999740023_R06C02","201004900096_R05C02","201004900096_R06C02","202702240054_R01C01","202702240054_R02C01","202702240079_R07C01","202702240079_R08C01","3999442124_R05C02","3999442124_R06C02","203751390020_R02C01","3999442124_R01C02","201004900096_R02C02","200999740005_R06C02","201004900018_R06C01","203751390017_R07C01","200687170042_R05C02","201004900096_R03C02","200999740023_R01C01","201004900018_R01C02","NA","NA.1","NA.2","NA.3","NA.4","NA.5","NA.6","NA.7")),]  

# #EXCLUDE DCD SAMPLES - Kidney dataset
# beta_values_case <- beta_values_filtered[,!(colnames(beta_values_filtered) %in% c("9296930129_R05C01", "9305651174_R01C01", "9305651174_R03C01", "9305651174_R02C02", "9305651174_R03C02", "9305651191_R02C02", "9305651191_R04C02", "201465900002_R04C01", "202240580106_R03C01", "202240580208_R03C01", "202259340119_R05C01", "202259350016_R04C01", "203496240002_R03C01", "203504430032_R05C01", "204001300109_R07C01", "204001300109_R08C01", "204001350016_R01C01", "202702240079_R06C01"))]
# #Removing rows based on the sample_name column in phenotype dataframe
# pheno_df_case <- pheno_df[!(pheno_df$sample_name %in% c("9296930129_R05C01", "9305651174_R01C01", "9305651174_R03C01", "9305651174_R02C02", "9305651174_R03C02", "9305651191_R02C02", "9305651191_R04C02", "201465900002_R04C01", "202240580106_R03C01", "202240580208_R03C01", "202259340119_R05C01", "202259350016_R04C01", "203496240002_R03C01", "203504430032_R05C01", "204001300109_R07C01", "204001300109_R08C01", "204001350016_R01C01", "202702240079_R06C01")),]
# 
# #Convert filter beta dataframe to m_value dataframe
m_values = beta2m(beta_values_filtered)

#String parsing for later output 
cat("String parsing\n")
string <- unlist(strsplit(pheno_file, "/"))
string <- rev(string)[1]
string <- substr(string, 1, 5)
if (string == "comp1") {
  string <- "eGFR_1month"
} else if (string == "comp2") {
  string <- "eGFR_12month"
} else if (string == "comp3") {
  string <- "eGFR_24month"
}

cat("PCA plots - Betas\n")
#Transpose resulting PCA plot - beta values 
transposed_beta <- t(beta_values_filtered) #t(beta_values_case)
transposed_beta <- data.frame(transposed_beta)
pca_general <- prcomp(transposed_beta, center=TRUE)
var_explained <- pca_general$sdev^2/sum(pca_general$sdev^2)
scores = as.data.frame(pca_general$x)
scores['Basename'] <- row.names(scores)
scores <- merge(scores, pheno_df, by = 'Basename')
output <- paste(output_dir, string, sep="")
output <- paste(output, "_", sep="")
output <- paste(output, comparison, sep="")
output_path <- paste(output, "_betas", sep="")
jpeg(paste(output_path, "_PCA.jpeg", sep=""), quality = 90)
plot <- ggplot(data=scores, mapping = aes(x = PC1, y = PC2, color=eGFR)) +theme_bw() + geom_point(aes(shape=time), alpha=0.5)+ labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"), y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) + scale_color_manual(values=pal) + ggtitle(paste(output_path, "_PCA", sep=""))
#plot + geom_text(aes(label = sample_name), size=3.5) + xlim(-100, 400)
print(plot)
dev.off()

cat("PCA plots - m_values\n")
#Transpose resulting PCA plot - m_values
transposed_m <- t(m_values) 
transposed_m <- data.frame(transposed_m)
pca_general <- prcomp(transposed_m, center=TRUE)
var_explained <- pca_general$sdev^2/sum(pca_general$sdev^2)
scores = as.data.frame(pca_general$x)
scores['Basename'] <- row.names(scores)
scores <- merge(scores, pheno_df, by = 'Basename')
output <- paste(output_dir, string, sep="")
output <- paste(output, "_", sep="")
output <- paste(output, comparison, sep="")
output_path <- paste(output, "_m", sep="")
jpeg(paste(output_path, "_PCA.jpeg", sep=""), quality = 90)
plot <- ggplot(data=scores, mapping = aes(x = PC1, y = PC2, color=eGFR)) +theme_bw() + geom_point(aes(shape=time), alpha=0.5)+ labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"), y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) + scale_color_manual(values=pal) + ggtitle(paste(output_path, "_PCA", sep=""))
#plot + geom_text(aes(label = sample_name), size=3.5) + xlim(-100, 400)
print(plot)
dev.off()

#CpGs - DMPs
cat("Identify CpGs\n")
if (comparison == 'K1_Low_K1_High') {
  condition <- factor(pheno_df$eGFR)
  cols <- c("High", "Low")
} else if (comparison == 'K2_Low_K2_High') {
  condition <- factor(pheno_df$eGFR)
  cols <- c("High", "Low")
} else if (comparison == 'K1_High_K2_High') {
  condition <- factor(pheno_df$time)
  cols <- c("K2", "K1")
} else if (comparison == 'K1_Low_K2_Low') {
  condition <- factor(pheno_df$time)
  cols <- c("K2", "K1")
} else if (comparison == 'K1_High_K2_Low') {
  condition <- factor(pheno_df$time)
  cols <- c("K2", "K1")
} else {
  condition <- factor(pheno_df$time)
  cols <- c("K2", "K1")
}

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

library(limma)
design <- model.matrix(~0+condition, data=pheno_df)
colnames(design) <- cols
fit1 <- lmFit(beta_values_filtered, design)

if (comparison == 'K1_Low_K1_High') {
  contMatrix <- makeContrasts(Low-High, levels=design)
} else if (comparison == 'K2_Low_K2_High') {
  contMatrix <- makeContrasts(Low-High, levels=design)
} else if (comparison == 'K1_High_K2_High') {
  contMatrix <- makeContrasts(K1-K2, levels=design)
} else if (comparison == 'K1_Low_K2_Low') {
  contMatrix <- makeContrasts(K1-K2, levels=design)
} else if (comparison == 'K1_High_K2_Low') {
  contMatrix <- makeContrasts(K1-K2, levels=design)
} else {
  contMatrix <- makeContrasts(K1-K2, levels=design)
}

#contMatrix

fit2 <- contrasts.fit(fit1, contMatrix)
fit2 <- eBayes(fit2)

#summary(decideTests(fit2))
ann450kSub <- ann450k[match(rownames(m_values),ann450k$Name), c(1:4,12:19,24:ncol(ann450k))]
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub)
DMPs_sig <- DMPs[which(DMPs$adj.P.Val < 0.05),]
write.csv(DMPs, file = paste(output, "_DMPs.csv", sep=""), row.names = TRUE)
write.csv(DMPs, file = paste(output, "_DMPs_sig.csv", sep=""), row.names = TRUE)
#plotCpg(m_values, cpg=rownames(DMPs)[1:4], pheno=type, ylab = "Beta values") #plots individual probes

#Manhattan plot using the DMPs
cat("Generating manhattan plot from DMPs\n")
library(qqman)
library(DMRcate)
title <- paste(output, " (Adj. P-val)", sep="")
col=c("black","grey")
DMPs$chr = str_replace_all(DMPs$chr, 'chr', '')
DMPs$chr = as.numeric(DMPs$chr)
DMPs$pos = as.numeric(DMPs$pos)
jpeg(paste(output, "_manhattan.jpeg", sep=""), quality = 90)
manhattan(DMPs, chr="chr", bp="pos",, p="adj.P.Val", snp="Islands_Name", col=col, suggestiveline=(-log10(0.05)), main=title)
dev.off()

myAnnotation <- cpg.annotate(object = as.matrix(m_values), datatype = "array", what = "M",
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contMatrix,
                             coef = "Low - High", arraytype = "450K")

DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs, genome = "hg19")
write.csv(result.ranges, file=paste(output, "_DMRs.csv", sep=""), row.names = FALSE)

q()

