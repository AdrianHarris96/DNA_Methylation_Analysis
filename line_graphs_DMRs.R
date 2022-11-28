#Line graphs for DMRs 
rm(list = ls())

dev.off()
library(pacman)
p_unload(all)
library(rio)
library(tidyverse)
library(ggplot2)
library(lumi)
library(matrixStats)

#High Injury
DMPs_High <- import('/Users/adrianharris/Desktop/Mas_lab/week_07292022/liver_csvs_combat0728/DD_HI_L1-DD_HI_L2_DMPs.csv')
DMRs_High <- import('/Users/adrianharris/Desktop/Mas_lab/week_07292022/liver_csvs_combat0728/DD_HI_L1-DD_HI_L2_DMRs.csv')
DMRs_High <- DMRs_High[(DMRs_High$HMFDR < 0.05),]

#Low injury 
DMPs_Low <- import('/Users/adrianharris/Desktop/Mas_lab/week_07292022/liver_csvs_combat0728/DD_LI_L1-DD_LI_L2_DMPs.csv')
DMRs_Low <- import('/Users/adrianharris/Desktop/Mas_lab/week_07292022/liver_csvs_combat0728/DD_LI_L1-DD_LI_L2_DMRs.csv')
DMRs_Low <- DMRs_Low[(DMRs_Low$HMFDR < 0.05),]
 
#Loading the betas 
liver_betas <- import('/Users/adrianharris/Desktop/Mas_lab/week_08262022/liver_combat_betas/beta_values_combat.csv')
liver_betas <- liver_betas %>% rename(ProbeID = V1)
row.names(liver_betas) <- liver_betas$ProbeID

#Removal of the ProbeID column and adjustment to the colnames
liver_betas <- liver_betas %>% select(-c(ProbeID))
colnames(liver_betas) <- substring(colnames(liver_betas), 2)
# colnames(liver_betas)

#Phenotype for liver samples
liver_pheno <- import('/Users/adrianharris/Documents/dna_methylation_analysis/liver_sample_sheet.csv')

#Removal of control samples 
liver_pheno <- subset(liver_pheno, sample_group != 'control')
liver_pheno <- subset(liver_pheno, sample_group != 'Control')

#Removal of one unpaired sample
liver_pheno <- liver_pheno[!(liver_pheno$sample_name == 'V037L1'),]

#Removal of samples not on server 
liver_pheno <- liver_pheno[!(liver_pheno$Basename == '203751390020_R08C01' | liver_pheno$Basename == '203751390020_R01C01'),]

#Drop DCD samples
liver_pheno = liver_pheno[!(liver_pheno$donor_type == 'DCD'),]

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

#Filter down phenotype 
liver_pheno_High <- liver_pheno[(liver_pheno$condition == 'DD_HI_L1' | liver_pheno$condition == 'DD_HI_L2'),]
liver_pheno_High <- liver_pheno_High %>% arrange(condition)

liver_pheno_Low <- liver_pheno[(liver_pheno$condition == 'DD_LI_L1' | liver_pheno$condition == 'DD_LI_L2'),]
liver_pheno_Low <- liver_pheno_Low %>% arrange(condition)

#Genelist 
genelist <- c('NR4A2', 'JUN', 'FOS', 'CCL3', 'EGR3', 'KLF6', 'EGR2', 'NR4A3', 'DUSP5', 'LOC102724428', 'TRIB1', 'PHLDA1', 'THBS1', 'KLF4', 'PLK2', 'CYR61', 'SOCS3', 'CCL3L3', 'CCL3L1', 'SIK1')
genelist <- c('SLC2A3',
                'ARL14',
                'BAG3',
                'BCL2A1',
                'CXCL3',
                'HBA2',
                'IL1B',
                'PLAUR',
                'PMAIP1',
                'PTX3',
                'SLC2A14',
                'SLC2A3',
                'TNFAIP6')

no_overlap <- c()
genelist <- c('PXDN')

#Plotting loop
for (gene in genelist){
  gene <- paste("^", gene, sep="")
  DMRs_check <- DMRs_High[grep(gene,DMRs_High$overlapping.genes),] #Only thing that requires changing
  if (nrow(DMRs_check) == 0) {
    no_overlap <- c(no_overlap, gene)
    print("not here")
  } else {
    for (row in 1:nrow(DMRs_check)) {
      chrom <- DMRs_check[row, 'seqnames']
      chr_num <- substring(chrom, 4)
      length1 <- DMRs_check[row, 'start']
      length2 <- DMRs_check[row, 'end']
      gene <- DMRs_check[row, 'overlapping.genes']
      gene <- sub(", ", ",", gene)
      DMPs_DMR_High <- DMPs_High[(DMPs_High$chr == chrom & DMPs_High$pos >=	length1 & DMPs_High$pos <= length2),]
      DMPs_DMR_Low <- DMPs_Low[(DMPs_Low$chr == chrom & DMPs_Low$pos >=	length1 & DMPs_Low$pos <= length2),]
      
      #Using the probes of interest to filter the beta dataframe
      liver_betas_High <- liver_betas[(row.names(liver_betas) %in% DMPs_DMR_High$Name),]
      liver_betas_High <- liver_betas_High[,(colnames(liver_betas_High) %in% liver_pheno_High$Basename)]
      liver_betas_Low <- liver_betas[(row.names(liver_betas) %in% DMPs_DMR_Low$Name),]
      liver_betas_Low <- liver_betas_Low[,(colnames(liver_betas_Low) %in% liver_pheno_Low$Basename)]
      
      #Calculate mean/STD
      liver_betas_High$Mean1 <- rowMeans(liver_betas_High[,1:20])
      liver_betas_High$Mean2 <- rowMeans(liver_betas_High[,21:40])
      liver_betas_High$STD1 <- rowSds(as.matrix(liver_betas_High[,1:20]))
      liver_betas_High$STD2 <- rowSds(as.matrix(liver_betas_High[,21:40]))
      
      liver_betas_Low$Mean1 <- rowMeans(liver_betas_Low[,1:15])
      liver_betas_Low$Mean2 <- rowMeans(liver_betas_Low[,16:30])
      liver_betas_Low$STD1 <- rowSds(as.matrix(liver_betas_Low[,1:15]))
      liver_betas_Low$STD2 <- rowSds(as.matrix(liver_betas_Low[,16:30]))
      
      #High/Low
      liver_betas_High$Name <- row.names(liver_betas_High)
      DMPs_DMR_High <- DMPs_DMR_High %>% select(c(Name, pos))
      liver_betas_High <- merge(liver_betas_High, DMPs_DMR_High, by = 'Name')
      final_df1 <- liver_betas_High %>% select(pos, Mean=Mean1, STD=STD1)
      final_df1$group <- 'DD_HI_L1'
      final_df2 <- liver_betas_High %>% select(pos, Mean=Mean2, STD=STD2)
      final_df2$group <- 'DD_HI_L2'
      final_df_High <- rbind(final_df1, final_df2)
      
      liver_betas_Low$Name <- row.names(liver_betas_Low)
      DMPs_DMR_Low <- DMPs_DMR_Low %>% select(c(Name, pos))
      liver_betas_Low <- merge(liver_betas_Low, DMPs_DMR_Low, by = 'Name')
      final_df1 <- liver_betas_Low %>% select(pos, Mean=Mean1, STD=STD1)
      final_df1$group <- 'DD_LI_L1'
      final_df2 <- liver_betas_Low %>% select(pos, Mean=Mean2, STD=STD2)
      final_df2$group <- 'DD_LI_L2'
      final_df_Low <- rbind(final_df1, final_df2)
      
      final_df <- rbind(final_df_High, final_df_Low)
      
      #Plotting 
      filename = paste("~/Desktop/DMR_High_L1-L2graphs/DD_HI_L1-DD_HI_L2-", gene, sep="") #Change name here 
      filename <- paste(filename, "-1.jpeg", sep="")
      
      p <- ggplot(final_df_High, mapping = aes(x=pos, y=Mean, group=group, color=group)) + geom_line() + geom_point()
      p<- p + labs(title=paste("L1 vs L2\nOverlapping genes:", gene, sep=" "), x=paste("Position at Chromosome", chr_num, sep= " "), y = "Methylation Level (β)") + theme_bw() + geom_errorbar(aes(ymin=Mean-STD, ymax=Mean+STD, x=pos), width=((length2-length1)/75))
      substrRight <- function(x, n){
        substr(x, nchar(x)-n+1, nchar(x)-n+1)
      }
      substrLeft <- function(x, n){
        substr(x, 1, nchar(x)-n)
      }
      while (file.exists(filename)) {
        num <- substrRight(filename, 6)
        num <- as.double(num)
        num <- num + 1
        num <- as.character(num)
        filename <- substrLeft(filename, 6)
        filename <- paste(filename, num, sep="")
        filename <- paste(filename, ".jpeg", sep="")
      }
      ggsave(file = filename, units = c("in"), width=5, height=5, dpi=300, p)
    }
  }
}
print(no_overlap)

#Adding genome tracks
library(rtracklayer)
library(Gviz)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
gen <- genome(cpgIslands)
#Ideogram track


head(Hsapiens)
#Ideogram
strack <- SequenceTrack(Hsapiens, chromosome = "chrX")
ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = "chrX")
grtrack <- GeneRegionTrack(Hsapiens,
                           chromosome = "chrX", name = "Gene Model", 
                           transcriptAnnotation = "symbol",
                           background.panel = "#FFFEDB",
                           background.title = "darkblue")
plotTracks(strack, from = 85000000, to = 129000000, showId = FALSE, showBandId = TRUE, cex.bands = 0.4)
data(geneModels)
head(geneModels)
grtrack <- GeneRegionTrack(geneModels, genome = gen,
                           chromosome = chr, name = "Gene Model", 
                           transcriptAnnotation = "symbol",
                           background.panel = "#FFFEDB",
                           background.title = "darkblue")
head(DMPs_Low)

#Clustering 
library(lumi)
liver_betas <- liver_betas[, colnames(liver_betas) %in% liver_pheno$Basename]
m_values <- beta2m(liver_betas)

library(RColorBrewer)
#Color Scheme Defined 
pal <- brewer.pal(4,"Dark2")

clustering <- function(type, pheno, timepoint, m_df) {
  #Filter beta dataframe using column names for the relevant comparison
  if (timepoint == 'L1&L2') {
    cat('Skip\n')
  } else {
    pheno <- pheno[(pheno$collection == timepoint),]
  }
  
  if (type == 'EPIC') {
    pheno <- pheno[!(pheno$array_type == "450K"),]
  } else if (type == '450K') {
    pheno <- pheno[!(pheno$array_type == "EPIC"),]
  }
  
  m_df <- m_df[,(colnames(m_df) %in% pheno$Basename)]
  print(dim(m_df))
  
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
  plot <- ggplot(data=scores, mapping = aes(x = PC1, y = PC2, color=array_type)) + geom_point(size=3, alpha=0.5) + labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"), y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) + ggtitle(title) + geom_text(aes(label = sample_name), size=1.75, colour="black") + scale_color_manual(values=c(pal[1], pal[2])) + theme_bw()
  if (length(unique(scores$donor_gender)) == 3) {
    plot2 <- ggplot(data=scores, mapping = aes(x = PC1, y = PC2, color=donor_gender)) + geom_point(size = 3, alpha = 0.5) + labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"), y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) + ggtitle(paste(title, " - donor gender and donor age", sep="")) + geom_text(aes(label=donor_age), size=1.75, colour="black") + scale_color_manual(values=c("grey", pal[1], pal[2])) + theme_bw()
  } else {
    plot2 <- ggplot(data=scores, mapping = aes(x = PC1, y = PC2, color=donor_gender)) + geom_point(size = 3, alpha = 0.5) + labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"), y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) + ggtitle(paste(title, " - donor gender and donor age", sep="")) + geom_text(aes(label=donor_age), size=1.75, colour="black") + scale_color_manual(values=c(pal[1], pal[2])) + theme_bw()
  }
  
  if (timepoint == 'L1&L2') {
    plot3 <- ggplot(data=scores, mapping = aes(x = PC1, y = PC2, color=collection)) + geom_point(size=3, alpha=0.5) + labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"), y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) + ggtitle(paste(title, " - L1 vs L2", sep="")) + geom_text(aes(label = sample_name), size=1.75, colour="black") + scale_color_manual(values=c(pal[1], pal[2])) + theme_bw()
    grid.arrange(plot, plot2, plot3, ncol=1, nrow = 3)
  } else {
    library(gridExtra)
    grid.arrange(plot, plot2, ncol=1, nrow = 2)
  }
  
 
  # #Classical MDS
  # title <- paste(title, "_MDS", sep="_")
  # transposed_m <- t(m_df) 
  # transposed_m <- data.frame(transposed_m)
  # d <- dist(transposed_m) # euclidean distances between the rows
  # fit <- cmdscale(d,eig=TRUE, k=2)
  # fit <- data.frame(fit$points)
  # fit['Basename'] <- row.names(fit)
  # fit <- merge(fit, pheno, by = 'Basename')
  # plot <- ggplot(data=fit, mapping = aes(x = X1, y = X2, color=eGFR)) + geom_point(size=3, alpha=0.5) + labs(x="Coordinate 1", y="Coordinate 2") + ggtitle(title) + geom_text(aes(label = sample_name), size=1.75, colour="black") + scale_color_manual(values=c(pal[1], pal[2])) + theme_bw()
  # if (length(unique(fit$donor_gender)) == 3) {
  #   plot2 <- ggplot(data=fit, mapping = aes(x = X1, y = X2, color=donor_gender)) + geom_point(size = 3, alpha = 0.5) + labs(x="Coordinate 1", y="Coordinate 2") + ggtitle(paste(title, " - donor gender and donor age", sep="")) + geom_text(aes(label=donor_age), size=1.75, colour="black") + scale_color_manual(values=c("grey", pal[1], pal[2])) + theme_bw()
  # } else {
  #   plot2 <- ggplot(data=fit, mapping = aes(x = X1, y = X2, color=donor_gender)) + geom_point(size = 3, alpha = 0.5) + labs(x="Coordinate 1", y="Coordinate 2") + ggtitle(paste(title, " - donor gender and donor age", sep="")) + geom_text(aes(label=donor_age), size=1.75, colour="black") + scale_color_manual(values=c(pal[1], pal[2])) + theme_bw()
  # }
  # 
  # if (timepoint == 'K1&K2') {
  #   plot3 <- ggplot(data=fit, mapping = aes(x = X1, y = X2, color=collection)) + geom_point(size=3, alpha=0.5) + labs(x="Coordinate 1", y="Coordinate 2") + ggtitle(paste(title, " - K1 vs K2", sep="")) + geom_text(aes(label = sample_name), size=1.75, colour="black") + scale_color_manual(values=c(pal[1], pal[2])) + theme_bw()
  #   library(gridExtra)
  #   grid.arrange(plot, plot2, plot3, ncol=1)
  # } else {
  #   library(gridExtra)
  #   grid.arrange(plot, plot2, ncol=1)
  # }
  return('Plotting\n')
}

pdf(file="~/Desktop/liver_clustering.pdf")
typeList <- c('Combined')
for (array in typeList) {
  clustering(array, liver_pheno, 'L1', m_values)
  clustering(array, liver_pheno, 'L2', m_values)
  clustering(array, liver_pheno, 'L1&L2', m_values)
}

dev.off()

#Kidney
rm(list = ls())
kidney_betas <- import('/Users/adrianharris/Desktop/large_files_kidneyEPIC/beta_values.csv')
kidney_betas <- kidney_betas %>% rename(ProbeID = V1)
rownames(kidney_betas) <- kidney_betas$ProbeID

kidney_pheno <- import('/Users/adrianharris/Documents/dna_methylation_analysis/paired_EPICkidney_sample_sheet.csv')

#Must remove outlier sample, 203504430032_R01C01 (and its paired sample 203504430032-R02C01)
kidney_pheno <- kidney_pheno[!(kidney_pheno$Basename == '203504430032_R01C01' | kidney_pheno$Basename == '203504430032_R02C01'),]

#Drop other outlier sample 
kidney_pheno <- kidney_pheno[!(kidney_pheno$sample_id == 'KUT4_K2' | kidney_pheno$sample_id == 'KUT4_K1'),]



pheno <- subset(kidney_pheno, select = -c(eGFR_24month))
pheno <- pheno %>% rename("eGFR1" = "eGFR_1month")

pheno <- pheno %>% rename("eGFR2" = "eGFR_12month")

#Addition of condition column 
pheno$condition <- 'NA'

#Filling of condition column
for (row in 1:nrow(pheno)) {
  if (pheno[row, 'time'] == 'K1') {
    if (pheno[row, 'eGFR1'] == 'High' & pheno[row, 'eGFR2'] == 'High') {
      pheno[row, 'condition'] <- "K1_High_High"
    } else if (pheno[row, 'eGFR1'] == 'Low' & pheno[row, 'eGFR2'] == 'Low') {
      pheno[row, 'condition'] <- "K1_Low_Low"
    }
  } else {
    if (pheno[row, 'eGFR1'] == 'High' & pheno[row, 'eGFR2'] == 'High') {
      pheno[row, 'condition'] <- "K2_High_High"
    } else if (pheno[row, 'eGFR1'] == 'Low' & pheno[row, 'eGFR2'] == 'Low') {
      pheno[row, 'condition'] <- "K2_Low_Low"
    }
  }
}

#if condition is NA, drop it
pheno <- pheno[!(pheno$condition == 'NA'),]

DMPs <- import('/Users/adrianharris/Desktop/kidney_status_csvs/eGFR_1month-eGFR_12month-K1_Low_Low-K1_High_High_DMPs.csv')


DMPs2 <- import('/Users/adrianharris/Desktop/kidney_status_csvs/eGFR_1month-eGFR_12month-K2_Low_Low-K2_High_High_DMPs.csv')

DMRs <- import('/Users/adrianharris/Desktop/kidney_status_csvs/eGFR_1month-eGFR_12month-K1_Low_Low-K1_High_High_DMRs.csv')
DMRs2 <- import('/Users/adrianharris/Desktop/kidney_status_csvs/eGFR_1month-eGFR_12month-K2_Low_Low-K2_High_High_DMRs.csv')

DMRs <- DMRs[(DMRs$HMFDR < 0.05),]

DMRs2 <- DMRs2[(DMRs2$HMFDR < 0.05),]

#Filter down datasets 
kidney_pheno1 <- pheno[(pheno$condition == 'K1_Low_Low' | pheno$condition == 'K1_High_High'),]
kidney_pheno1 <- kidney_pheno1 %>% arrange(condition)

kidney_pheno2 <- pheno[(pheno$condition == 'K2_Low_Low' | pheno$condition == 'K2_High_High'),]
kidney_pheno2 <- kidney_pheno2 %>% arrange(condition)

#Filter down betas given the dataset and match order
kidney_betas1 <- kidney_betas[,match(kidney_pheno1$Basename, colnames(kidney_betas))]
colnames(kidney_betas1) <- kidney_pheno1$sample_id

kidney_betas2 <- kidney_betas[,match(kidney_pheno2$Basename, colnames(kidney_betas))]
colnames(kidney_betas2) <- kidney_pheno2$sample_id

chrom = 'chr17'
length1 <- 15902694
length2 <- 15903096

#Change this to two
DMPs_DMR <- DMPs[(DMPs$chr == chrom & DMPs$pos >=	length1 & DMPs$pos <= length2),]
betas_DMRs <- kidney_betas1[(row.names(kidney_betas1) %in% DMPs_DMR$Name),]

betas_DMRs$Mean1 <- rowMeans(betas_DMRs[,1:23])
betas_DMRs$Mean2 <- rowMeans(betas_DMRs[,24:27])

betas_DMRs$STD1 <- rowSds(as.matrix(betas_DMRs[,1:23]))
betas_DMRs$STD2 <- rowSds(as.matrix(betas_DMRs[,24:27]))

betas_DMRs$Name <- row.names(betas_DMRs)
DMPs_DMR <- DMPs_DMR %>% select(c(Name, pos))
betas_DMRs <- merge(betas_DMRs, DMPs_DMR, by = 'Name')
final_df1 <- betas_DMRs %>% select(pos, Mean=Mean1, STD=STD1)
final_df1$group <- 'K1_High_High'
final_df2 <- betas_DMRs %>% select(pos, Mean=Mean2, STD=STD2)
final_df2$group <- 'K1_Low_Low'
final_df <- rbind(final_df1, final_df2)

p <- ggplot(final_df, mapping = aes(x=pos, y=Mean, group=group, color=group)) + geom_line() + geom_point() 
p <- p + labs(title="K1_High_High vs K1_Low_Low", x=paste("Position at Chromosome", '17', sep= " "), y = "Methylation Level (β)") + theme_bw() + geom_errorbar(aes(ymin=Mean-STD, ymax=Mean+STD, x=pos), width=((length2-length1)/75))
print(p)
