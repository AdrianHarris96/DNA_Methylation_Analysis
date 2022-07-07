#!/home/amharris/.conda/envs/methylation_env/bin/R Rscript

#Example input: Rscript clustering_timepoints.R -s /home/amharris/dna_methylation_analysis/eGFR_1month_sample_sheet.csv -b /local/projects-t3/XVMAS_Lab/Projects_2022/XVMAS_P09_methylation/01_analysis/kidneyTx_methylation/ -c /home/amharris/dna_methylation_analysis/ -o /local/projects-t3/XVMAS_Lab/Projects_2022/XVMAS_P09_methylation/01_analysis/kidneyTx_methylation/clustering_out/
library(minfi)
library(rio)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(lumi)
library(optparse)
library(gridExtra)

option_list = list(
  make_option(c("-s", "--sample"), type="character", default=NULL, 
              help="sample sheet CSV"),
  make_option(c("-b", "--base_dir"), type="character", default=NULL, 
              help="base directory with IDAT files"),
  make_option(c("-c", "--git_dir"), type="character", default=NULL, 
              help="git directory with script"),
  make_option(c("-o", "--out_dir"), type="character", default=NULL, 
              help="directory to output")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Color Scheme Defined 
pal <- brewer.pal(4,"Dark2")

calculate_betas <- function(pheno_file, base_dir, git_dir, output_dir) {
  #Importing manually-curated sample sheet
  pheno_df <- import(pheno_file)
  
  #Must remove outlier sample, 203504430032_R01C01 (and its paired sample 203504430032-R02C01)
  pheno_df <- pheno_df[!(pheno_df$Basename == '203504430032_R01C01' | pheno_df$Basename == '203504430032_R02C01'),]
  
  nrow(subset(pheno_df, array_type == '450K'))
  nrow(subset(pheno_df, array_type == 'EPIC'))
  
  #Split into phenotype file
  pheno450k <- pheno_df[!(pheno_df$array_type == 'EPIC'),]
  phenoEPIC <- pheno_df[!(pheno_df$array_type == '450K'),]
  
  print(base_dir)
  #Specify the directories and read in respective IDAT Files 
  dir450k <- paste(base_dir, "450k_array", sep="")
  dirEPIC <- paste(base_dir, "EPIC_array", sep="")
  
  #Generation or Loading of rgSet
  if (file.exists(paste(output_dir, "m_values.csv", sep=""))) {
    cat('Skipping rgSet Loading\n')
  } else if (file.exists(paste(output_dir, "rgSet.RDS", sep=""))) {
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
  if (file.exists(paste(output_dir, "m_values.csv", sep=""))) {
    cat('Skipping normalizaion\n')
  }else if (file.exists(paste(output_dir, "mtSet.RDS", sep=""))) {
    cat('Loading normalization\n')
    mtSet <- readRDS(paste(output_dir, "mtSet.RDS", sep=""))
  } else {
    cat('Performing normalization\n')
    #Normalization and plotting 
    mtSet <- preprocessNoob(rgSet)
    saveRDS(mtSet, file = paste(output_dir, "mtSet.RDS", sep=""))
  }
  
  # Map to Genome
  if (file.exists(paste(output_dir, "m_values2.csv", sep=""))) {
    cat('Load m-values\n')
    m_values <- import(paste(output_dir, "m_values.csv", sep=""))
  } else {
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
    rm(keep.probes, detP)
    
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
    library(ChAMP)
    beta_values <- champ.runCombat(beta=beta_values, pd=pData(gmtSet), variablename="array_type", logitTrans=TRUE)
    
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
    m_values <- beta2m(beta_values_filtered)
    write.csv(beta_values_filtered, file = paste(output_dir, "beta_values2.csv", sep=""), row.names = TRUE)
    write.csv(m_values, file = paste(output_dir, "m_values2.csv", sep=""), row.names = TRUE)
    rm(beta_values_filtered)
  }
  
  typeList <- c('450K', 'EPIC', 'Combined')
  for (array in typeList) {
    clustering(array, pheno_df, 'K1', m_values, pheno_file)
    clustering(array, pheno_df, 'K2', m_values, pheno_file)
    clustering(array, pheno_df, 'K1&K2', m_values, pheno_file)
  }
  return(paste("Done with plotting for file: ", pheno_file, sep=""))
}

clustering <- function(type, pheno, timepoint, m_df, pheno_file) {
  #Filter beta dataframe using column names for the relevant comparison
  if (timepoint == 'K1&K2') {
    cat('Skip\n')
  } else {
    pheno <- pheno[(pheno$time == timepoint),]
  }
  
  if (type == 'EPIC') {
    pheno <- pheno[!(pheno$array_type == "450K"),]
  } else if (type == '450K') {
    pheno <- pheno[!(pheno$array_type == "EPIC"),]
  }
  
  m_df <- m_df[,(colnames(m_df) %in% pheno$Basename)]
  #print(dim(m_df))
  
  #Parsing name
  title <- paste(timepoint, type, sep="_")
  
  #PCA - m_values 
  title <- paste(title, "_PCA", sep="_")
  transposed_m <- t(m_df)
  transposed_m <- data.frame(transposed_m)
  pca_general <- prcomp(transposed_m, center=TRUE)
  var_explained <- pca_general$sdev^2/sum(pca_general$sdev^2)
  scores <- as.data.frame(pca_general$x)
  scores['Basename'] <- row.names(scores)
  scores <- merge(scores, pheno, by = 'Basename')
  #Converting NAs to 0 in donor_age
  vec <- scores$donor_age
  vec[is.na(vec)] <- 0
  scores$donor_age <- vec
  plot <- ggplot(data=scores, mapping = aes(x = PC1, y = PC2, color=array_type)) + geom_point(size=3, alpha=0.5) + labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"), y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) + ggtitle(title) + geom_text(aes(label = sample_id), size=1.75, colour="black") + scale_color_manual(values=c(pal[1], pal[2])) + theme_bw()
  if (length(unique(scores$donor_gender)) == 3) {
    plot2 <- ggplot(data=scores, mapping = aes(x = PC1, y = PC2, color=donor_gender)) + geom_point(size = 3, alpha = 0.5) + labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"), y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) + ggtitle(paste(title, " - donor gender and donor age", sep="")) + geom_text(aes(label=donor_age), size=1.75, colour="black") + scale_color_manual(values=c("grey", pal[1], pal[2])) + theme_bw()
  } else {
    plot2 <- ggplot(data=scores, mapping = aes(x = PC1, y = PC2, color=donor_gender)) + geom_point(size = 3, alpha = 0.5) + labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"), y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) + ggtitle(paste(title, " - donor gender and donor age", sep="")) + geom_text(aes(label=donor_age), size=1.75, colour="black") + scale_color_manual(values=c(pal[1], pal[2])) + theme_bw()
  }
  
  if (timepoint == 'K1&K2') {
    plot3 <- ggplot(data=scores, mapping = aes(x = PC1, y = PC2, color=time)) + geom_point(size=3, alpha=0.5) + labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"), y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) + ggtitle(paste(title, " - K1 vs K2", sep="")) + geom_text(aes(label = sample_id), size=1.75, colour="black") + scale_color_manual(values=c(pal[1], pal[2])) + theme_bw()
    grid.arrange(plot, plot2, plot3, ncol=1, nrow = 3)
  } else {
    library(gridExtra)
    grid.arrange(plot, plot2, ncol=1, nrow = 2)
  }
  
  if (length(unique(scores$eGFR_1month)) == 3) {
    month1_col <- c("grey", pal[1], pal[2])
  } else {
    month1_col <- c(pal[1], pal[2])
  }
  
  if (length(unique(scores$eGFR_12month)) == 3) {
    month12_col <- c("grey", pal[1], pal[2])
  } else {
    month12_col <- c(pal[1], pal[2])
  }
  
  if (length(unique(scores$eGFR_24month)) == 3) {
    month24_col <- c("grey", pal[1], pal[2])
  } else {
    month24_col <- c(pal[1], pal[2])
  }
  
  plot4 <- ggplot(data=scores, mapping = aes(x = PC1, y = PC2, color=eGFR_1month)) + geom_point(size=3, alpha=0.5) + labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"), y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) + ggtitle(title) + geom_text(aes(label = sample_id), size=1.75, colour="black") + scale_color_manual(values=month1_col) + theme_bw()
  plot5 <- ggplot(data=scores, mapping = aes(x = PC1, y = PC2, color=eGFR_12month)) + geom_point(size=3, alpha=0.5) + labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"), y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) + ggtitle(title) + geom_text(aes(label = sample_id), size=1.75, colour="black") + scale_color_manual(values=month12_col) + theme_bw()
  plot6 <- ggplot(data=scores, mapping = aes(x = PC1, y = PC2, color=eGFR_24month)) + geom_point(size=3, alpha=0.5) + labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"), y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) + ggtitle(title) + geom_text(aes(label = sample_id), size=1.75, colour="black") + scale_color_manual(values=month24_col) + theme_bw()
  
  grid.arrange(plot4, plot5, plot6, ncol=1)
  # #Classical MDS
  # title <- paste(title, "_MDS", sep="_")
  # transposed_m <- t(m_df) 
  # transposed_m <- data.frame(transposed_m)
  # d <- dist(transposed_m) # euclidean distances between the rows
  # fit <- cmdscale(d,eig=TRUE, k=2)
  # fit <- data.frame(fit$points)
  # fit['Basename'] <- row.names(fit)
  # fit <- merge(fit, pheno, by = 'Basename')
  # plot <- ggplot(data=fit, mapping = aes(x = X1, y = X2, color=eGFR)) + geom_point(size=3, alpha=0.5) + labs(x="Coordinate 1", y="Coordinate 2") + ggtitle(title) + geom_text(aes(label = sample_id), size=1.75, colour="black") + scale_color_manual(values=c(pal[1], pal[2])) + theme_bw()
  # if (length(unique(fit$donor_gender)) == 3) {
  #   plot2 <- ggplot(data=fit, mapping = aes(x = X1, y = X2, color=donor_gender)) + geom_point(size = 3, alpha = 0.5) + labs(x="Coordinate 1", y="Coordinate 2") + ggtitle(paste(title, " - donor gender and donor age", sep="")) + geom_text(aes(label=donor_age), size=1.75, colour="black") + scale_color_manual(values=c("grey", pal[1], pal[2])) + theme_bw()
  # } else {
  #   plot2 <- ggplot(data=fit, mapping = aes(x = X1, y = X2, color=donor_gender)) + geom_point(size = 3, alpha = 0.5) + labs(x="Coordinate 1", y="Coordinate 2") + ggtitle(paste(title, " - donor gender and donor age", sep="")) + geom_text(aes(label=donor_age), size=1.75, colour="black") + scale_color_manual(values=c(pal[1], pal[2])) + theme_bw()
  # }
  # 
  # if (timepoint == 'K1&K2') {
  #   plot3 <- ggplot(data=fit, mapping = aes(x = X1, y = X2, color=time)) + geom_point(size=3, alpha=0.5) + labs(x="Coordinate 1", y="Coordinate 2") + ggtitle(paste(title, " - K1 vs K2", sep="")) + geom_text(aes(label = sample_id), size=1.75, colour="black") + scale_color_manual(values=c(pal[1], pal[2])) + theme_bw()
  #   library(gridExtra)
  #   grid.arrange(plot, plot2, plot3, ncol=1)
  # } else {
  #   library(gridExtra)
  #   grid.arrange(plot, plot2, ncol=1)
  # }
  return('Plotting\n')
}

#Check if output directory exists
if (file.exists(opt$out_dir)) {
  cat("Directory already exists\n")
} else {
  dir.create(opt$out_dir)
}


pdf(file = paste(opt$out_dir, "out_combat.pdf", sep=""))
calculate_betas(pheno_file = opt$sample, 
                base_dir =opt$base_dir, 
                git_dir = opt$git_dir, 
                output_dir = opt$out_dir)

# pdf(file = paste('/Users/adrianharris/Desktop/test_kidney/', "out.pdf", sep=""))
# calculate_betas(pheno_file = '/Users/adrianharris/Documents/dna_methylation_analysis/kidney_sample_sheet.csv', 
#                 base_dir ='/Users/adrianharris/Desktop/kidney/', 
#                 git_dir = '/Users/adrianharris/Documents/dna_methylation_analysis/', 
#                 output_dir = '/Users/adrianharris/Desktop/test_kidney/')

dev.off()

#hierarchical clustering
# rm(list=ls())
# library(rio)
# output_dir = '/Users/adrianharris/Desktop/test_kidney/'
# m_values <- import(paste(output_dir, "m_values.csv", sep=""))
# beta_values <- import(paste(output_dir, "beta_values.csv", sep=""))
# row.names(beta_values) <- beta_values$V1
# beta_values <- beta_values[, 2:ncol(beta_values)]

# pheno_df <- import('/Users/adrianharris/Documents/dna_methylation_analysis/kidney_sample_sheet.csv')
# pheno_df <- pheno_df[!(pheno_df$Basename == '203504430032_R01C01' | pheno_df$Basename == '203504430032_R02C01'),]
# phenoK1 <- pheno_df[(pheno_df$time == 'K1'),]
# phenoK2 <- pheno_df[(pheno_df$time == 'K2'),]
# 
# K1 <- phenoK1$Basename
# K2 <- phenoK2$Basename
# 
# beta_k1 <- beta_values[, colnames(beta_values) %in% K1]
# beta_k2 <- beta_values[, colnames(beta_values) %in% K2]
# 
# beta_k1$mean <- rowMeans(beta_k1)
# beta_k2$mean <- rowMeans(beta_k2)
# 
# df <- data.frame(matrix(ncol=2, nrow=356361))
# colnames(df) <- c('k1', 'k2')
# df$k1 <- beta_k1$mean
# df$k2 <- beta_k2$mean
# 
# library(ggplot2)
# plot <- ggplot(data=df, mapping = aes(x = k1, y = k2)) + geom_point(size=1, alpha = 0.5) + labs(x="average beta for K1", y="average beta for k2")
# print(plot)
# 
# 
# clusters <- hclust(dist(m_values))


