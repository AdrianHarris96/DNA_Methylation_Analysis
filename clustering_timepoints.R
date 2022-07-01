#Clustering timepoints
rm(list = ls())
library(minfi)
library(rio)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

#Example input: Rscript methylation_analysis.R <pheno_file> <base_dir> <output_dir>
args=commandArgs(trailingOnly=TRUE)
pheno_file = args[1]
base_dir = args[2]
git_dir = args[3]
output_dir = args[4]

#Local Machine
pheno_file = '/Users/adrianharris/Documents/dna_methylation_analysis/eGFR_1month_sample_sheet.csv'
base_dir = '/Users/adrianharris/Desktop/kidney/'
git_dir = '/Users/adrianharris/Documents/dna_methylation_analysis/'
output_dir = '/Users/adrianharris/Desktop/kidney/'

if (file.exists(output_dir)) {
  cat("Directory already exists\n")
} else {
  dir.create(output_dir)
}

#Importing manually-curated sample sheet
pheno_df <- import(pheno_file)

#Must remove outlier sample, 203504430032_R01C01 (and its paired sample 203504430032-R02C01)
pheno_df <- pheno_df[!(pheno_df$Basename == '203504430032_R01C01' | pheno_df$Basename == '203504430032_R02C01'),]

nrow(subset(pheno_df, array_type == '450K'))
nrow(subset(pheno_df, array_type == 'EPIC'))

#Split into phenotype file 
pheno450k <- pheno_df[!(pheno_df$array_type == 'EPIC'),]
phenoEPIC <- pheno_df[!(pheno_df$array_type == '450K'),]

#Specify the directories and read in respective IDAT Files 
dir450k <- paste(base_dir, "450k_array", sep="")
dirEPIC <- paste(base_dir, "EPIC_array", sep="")

#Generation or Loading of rgSet
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

#Generation or Loading of mtSet
if (file.exists(paste(output_dir, "postNormQC.jpeg", sep=""))) {
  cat('Loading normalization\n')
  mtSet <- readRDS(paste(output_dir, "mtSet.RDS", sep=""))
} else {
  cat('Performing normalization and plotting\n')
  #Normalization and plotting 
  mtSet <- preprocessNoob(rgSet)
  saveRDS(mtSet, file = paste(output_dir, "mtSet.RDS", sep=""))
}

# Map to Genome
cat('Converting to Genomic Methyl Set\n')
rSet <- ratioConvert(mtSet, what = "both", keepCN = TRUE)
gmtSet <- mapToGenome(rSet)
dim(gmtSet) #Number of probes = 452453

#Removing probes that include SNPs
snps <- getSnpInfo(gmtSet)
gmtSet <- addSnpInfo(gmtSet)
gr <- granges(gmtSet)
gmtSet <- dropLociWithSnps(gmtSet, snps=c("SBE","CpG"), maf=0)
dim(gmtSet) #Number of probes = 436144
rm(snps, gr)

# Filter unwanted sites 
detP <- detectionP(rgSet)
detP <- detP[match(featureNames(gmtSet),rownames(detP)),]
rm(rgSet)

# remove any probes that have failed in >50% of samples
keep <- detP < 0.01
keep.probes <- rownames(detP[rowMeans(keep)>=0.5,]) #probes that failed detection in more than half of the samples
gmtSet <- gmtSet[keep.probes,] 
dim(gmtSet) #Number of probes = 436128
rm(keep.probes)

#Remove probes that located on the X or Y chromosome
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

keep <- !(featureNames(gmtSet) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")]) #remove probes that are not of the chrom x or y
table(keep)
gmtSet <- gmtSet[keep,]
dim(gmtSet) #Number of probes = 425718
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
dim(gmtSet) #Number of probes = 392871
rm(cross.react, multi.map, bad.probes, cross.react.probes, multi.map.probes, keep)

#Extract betas and m_values
beta_values <- getBeta(gmtSet)

#Remove probes Hypomethylated in all samples identified by CpG sites having beta < 0.05 in all samples
beta_values_filtered <- as.data.frame(beta_values) 
rm(beta_values)
cat('Hypomethylated\n')
dim(filter_all(beta_values_filtered, all_vars(. < 0.05))) #Number of hypomethylated
beta_values_filtered <- filter_all(beta_values_filtered, any_vars(. >= 0.05)) 
#33419 Hypomethylated

#Remove probes Hypermethylated in all samples identified by CpG sites having beta > 0.95 in all samples
cat("Hypermethylated\n")
dim(filter_all(beta_values_filtered, all_vars(. > 0.95)))
#3091 hypermethylated
beta_values_filtered <- filter_all(beta_values_filtered, any_vars(. < 0.95)) 
dim(beta_values_filtered) #Final dimensions = 356361

#Color Scheme Defined 
pal <- brewer.pal(4,"Dark2")

library(lumi)
clustering <- function(pheno, timepoint, betas_df) {
  #Filter beta dataframe using column names for the relevant comparison
  pheno <- pheno[(pheno$time == timepoint),]
  pheno <- pheno[(pheno$eGFR == ""),]
  betas_df <- betas_df[,(colnames(betas_df) %in% pheno$Basename)]
  print(dim(betas_df))
  #m_values = beta2m(betas_df)
  
  #Parsing name
  string <- unlist(strsplit(pheno_file, "/"))
  string <- rev(string)[1]
  print(string)
  string <- unlist(strsplit(string, "_"))
  string <- paste(string[1], string[2], sep="_")
  
  #PCA - Betas
  transposed_beta <- t(betas_df) 
  transposed_beta <- data.frame(transposed_beta)
  pca_general <- prcomp(transposed_beta, center=TRUE)
  var_explained <- pca_general$sdev^2/sum(pca_general$sdev^2)
  scores = as.data.frame(pca_general$x)
  scores['Basename'] <- row.names(scores)
  scores <- merge(scores, pheno, by = 'Basename')
  plot <- ggplot(data=scores, mapping = aes(x = PC1, y = PC2, color=eGFR)) +theme_bw() + geom_point(aes(shape=), alpha=0.5, size=2)+ labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"), y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) + scale_color_manual(values=pal)
  plot + geom_text(aes(label = sample_name), size=3.5)
  print(plot)
}

clustering(pheno_df, 'K1', beta_values_filtered)

#Number of files 
nrow(subset(pheno_df, time == 'K1')) #144
nrow(subset(pheno_df, time == 'K2')) #48
