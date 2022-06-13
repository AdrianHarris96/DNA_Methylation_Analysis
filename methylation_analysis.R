#Script for DNA Methylation Analysis
rm(list = ls())
library(minfi)
library(shinyMethyl)
library(MEAL)
library(rio)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)

#Read data in 
dir450k <- "/Users/adrianharris/Desktop/kidney/450k_array"
dirEPIC <- "/Users/adrianharris/Desktop/kidney/EPIC_array"
rgSet450k <- read.metharray.exp(base=dir450k) #Loading this can be taxing 
rgSetEPIC <- read.metharray.exp(base=dirEPIC, force=TRUE)
rgSet <- combineArrays(rgSet450k, rgSetEPIC)
rm(rgSet450k)
rm(rgSetEPIC)

#Saving set as an RDS file 
#saveRDS(rgSet, file = "/Users/adrianharris/Desktop/kidney/rgSet.RDS")

#Importing manually-curated sample sheet
pheno_df <- import('/Users/adrianharris/Desktop/kidney/kidneyTx_methylation.csv')

#Adding new column - sample name
pheno_df['sample_name'] <- 'NA'
for (row in 1:nrow(pheno_df)) {
  sample_name <- paste(pheno_df[row, 'sentrix_ID'], pheno_df[row, 'sentrix_position'], sep="_")
  pheno_df[row, 'sample_name'] <- sample_name
  rm(sample_name)
}
rm(row)

#Removal of the V1 column, moving "sample_name" column to front, ordering that column
pheno_df <- pheno_df[,c(2:(ncol(pheno_df)))]
pheno_df <- pheno_df[,c(ncol(pheno_df),1:(ncol(pheno_df)-1))]
#pheno_df <- pheno_df[order(pheno_df$sample_name),]
#pheno_df <- pheno_df[!duplicated(pheno_df$sample_name),]

nrow(subset(pheno_df, array_type == '450K'))
nrow(subset(pheno_df, array_type == 'EPIC'))

#Calculate Detection p-values
detP <- detectionP(rgSet)
detP <- data.frame(detP)
sample_names <- colnames(detP)
sample_names <- substr(sample_names, 2, nchar(sample_names))
detPmeans <- colMeans(detP)
detP_df <- data.frame(matrix(ncol=2, nrow=235))
detP_df['X1'] <- sample_names
detP_df['X2'] <- detPmeans
colnames(detP_df) <- c('sample_name', 'p_value_mean')
#write.csv(detP_df,"/Users/adrianharris/Desktop/kidneyTx_p-values.csv", row.names = FALSE)

#Drop the top row before plotting barplot
detP_df<- detP_df[order(-detP_df$p_value_mean),]
detP_df<- detP_df[2:nrow(detP_df),]
detP_df <- merge(detP_df, pheno_df, by = 'sample_name')
pal <- brewer.pal(4,"Dark2")
par(mfrow=c(2,1))
#Plotting 450K barplot 
barplot((subset(detP_df, array_type == '450K'))$p_value_mean, col=pal[factor(detP_df$time)], names.arg=(subset(detP_df, array_type == '450K'))$sample_name, las=2, cex.names=0.4, cex.axis=0.5, space=0.5, ylab="Mean detection p-values", main='450K Array')
legend("topleft", legend=levels(factor(detP_df$time)), fill=pal,
       cex=0.27, bty = "n", bg="white")
#Plotting EPIC barplot
barplot((subset(detP_df, array_type == 'EPIC'))$p_value_mean, col=pal[factor(detP_df$time)], names.arg=(subset(detP_df, array_type == 'EPIC'))$sample_name, las=2, cex.names=0.4, cex.axis=0.5, space=0.5, ylab="Mean detection p-values", main='EPIC Array')
legend("topleft", legend=levels(factor(detP_df$time)), fill=pal,
       cex=0.27, bty = "n", bg="white")

par(mfrow=c(1,1))
#preprocessing QC 
mtSet <- preprocessRaw(rgSet)
qc <- getQC(mtSet)
plotQC(qc)
#?plotQC
qc <- data.frame(qc)
qc['sample_name'] <- row.names(qc)
qc <- merge(qc, pheno_df, by = 'sample_name')
plot <- ggplot(data=qc, mapping = aes(x = mMed, y = uMed, color=time)) + geom_point(aes(shape=array_type), alpha=0.5) + xlim(7, 14) + ylim(7, 14) + theme_bw()+ geom_abline(slope=-1, intercept = 21.25 , color="black", linetype="dashed", size=0.5) + scale_color_manual(values=pal)
print(plot)
densityPlot(mtSet, sampGroups =qc$array_type)
rm(detP_df)
rm(mtSet)
rm(qc)
rm(detP)

#Normalization
mtSet <- preprocessNoob(rgSet)
saveRDS(mtSet, file = "/Users/adrianharris/Desktop/kidney/mtSet.RDS")
postqc <- getQC(mtSet)
postqc <- data.frame(postqc)
postqc['sample_name'] <- row.names(postqc)
postqc <- merge(postqc, pheno_df, by = 'sample_name')
plot2 <- ggplot(data=postqc, mapping = aes(x = mMed, y = uMed, color=time)) + geom_point(aes(shape=array_type), alpha=0.5) + xlim(5, 14) + ylim(5, 14) + theme_bw()+ geom_abline(slope=-1, intercept = 21.25 , color="black", linetype="dashed", size=0.5) + scale_color_manual(values=pal)
print(plot2)
densityPlot(mtSet, sampGroups = postqc$array_type)

# Map to Genome
rSet <- ratioConvert(mtSet, what = "both", keepCN = TRUE)
gmtSet <- mapToGenome(rSet)
#dim(gmtSet)

#Predicted Sex 
predictedSex <- getSex(gmtSet, cutoff = -2)$predictedSex

#Add column to pheno_df
postqc['predicted_sex'] <- predictedSex
predictedsex_file <- postqc %>% select(c(sample_name, predicted_sex, gender_donor, gender_recipient))
write.csv(predictedsex_file,"/Users/adrianharris/Desktop/kidneyTx_methylation_predictedSex.csv", row.names = FALSE)

#Removing probes that include SNPs
snps <- getSnpInfo(gmtSet)
gmtSet <- addSnpInfo(gmtSet)
gr <- granges(gmtSet)
gmtSet <- dropLociWithSnps(gmtSet, snps=c("SBE","CpG"), maf=0)
#dim(gmtSet)

#Grab Betas and transpose resulting PCA plot 
beta <- getBeta(gmtSet)
beta <- data.frame(beta)
beta <- t(beta)
beta <- data.frame(beta)
pca_general <- prcomp(beta, center=TRUE)
var_explained <- pca_general$sdev^2/sum(pca_general$sdev^2)
scores = as.data.frame(pca_general$x)
scores['sample_name'] <- row.names(scores)

scores <- merge(scores, pheno_df, by = 'sample_name')
plot <- ggplot(data=scores, mapping = aes(x = PC1, y = PC2)) +theme_bw() + geom_point()+ labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"), y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))
print(plot)

M <- getM(gmtSet)
CN <- getCN(gmtSet)


