#Script for DNA Methylation Analysis
rm(list = ls())
library(minfi)
library(shinyMethyl)
library(MEAL)
library(rio)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(lumi)

#Importing manually-curated sample sheet
pheno_df <- import('/Users/adrianharris/Desktop/kidney/kidneyTx_methylation.csv')

#Adding new column - sample name
pheno_df['Basename'] <- 'NA'
for (row in 1:nrow(pheno_df)) {
  sample_name <- paste(as.numeric(pheno_df[row, 'sentrix_ID']), pheno_df[row, 'sentrix_position'], sep="_")
  pheno_df[row, 'Basename'] <- sample_name
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

#Split into phenotype file 
pheno450k <- pheno_df[!(pheno_df$array_type == 'EPIC'),]
phenoEPIC <- pheno_df[!(pheno_df$array_type == '450K'),]

#Read data in 
dir450k <- "/Users/adrianharris/Desktop/kidney/450k_array"
dirEPIC <- "/Users/adrianharris/Desktop/kidney/EPIC_array"
rgSet450k <- read.metharray.exp(base=dir450k, target=pheno450k)
rgSetEPIC <- read.metharray.exp(base=dirEPIC, target=phenoEPIC, force=TRUE)
rgSet <- combineArrays(rgSet450k, rgSetEPIC)
#Saving set as an RDS file 
saveRDS(rgSet450k, file = "/Users/adrianharris/Desktop/kidney/rgSet450k.RDS")
saveRDS(rgSetEPIC, file = "/Users/adrianharris/Desktop/kidney/rgSetEPIC.RDS")
saveRDS(rgSet, file = "/Users/adrianharris/Desktop/kidney/rgSet.RDS")
rm(rgSet450k)
rm(rgSetEPIC)

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
saveRDS(gmtSet, file = "/Users/adrianharris/Desktop/kidney/gmtSet.RDS")
dim(gmtSet) #Number of probes = 452453

#Predicted Sex 
predictedSex <- getSex(gmtSet, cutoff = -2)$predictedSex

#Add column to pheno_df
postqc['predicted_sex'] <- predictedSex
predictedsex_file <- postqc %>% select(c(sample_id, sample_name, predicted_sex, gender_donor, gender_recipient))
write.csv(predictedsex_file,"/Users/adrianharris/Desktop/kidneyTx_methylation_predictedSex.csv", row.names = FALSE)

#Removing probes that include SNPs
snps <- getSnpInfo(gmtSet)
gmtSet <- addSnpInfo(gmtSet)
gr <- granges(gmtSet)
gmtSet <- dropLociWithSnps(gmtSet, snps=c("SBE","CpG"), maf=0)
dim(gmtSet) #Number of probes = 436144 

# Filter unwanted sites 
# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(gmtSet),rownames(detP)),] 
detP_df <- data.frame(detP)

# remove any probes that have failed in >50% of samples
keep <- detP < 0.01
keep.probes <- rownames(detP[rowMeans(keep)>=0.5,]) #probes that failed detection in more than half of the samples
gmtSet <- gmtSet[keep.probes,] 
dim(gmtSet) #Number of probes = 436127 

#Remove probes that located on the X or Y chromosome
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

keep <- !(featureNames(gmtSet) %in% ann450k$Name[ann450k$chr %in% 
                                                   c("chrX","chrY")]) #remove probes that are not of the chrom x or y
table(keep)
gmtSet <- gmtSet[keep,]
dim(gmtSet) #Number of probes = 425716

#Creation of bad probes character and filter gmtSet
cross.react <- read.csv('/Users/adrianharris/Downloads/illumina450k_filtering/48639-non-specific-probes-Illumina450k.csv', head = T, as.is = T)
cross.react.probes <- as.character(cross.react$TargetID)
#Probes identified with potential hybridization issues
multi.map <- read.csv('/Users/adrianharris/Downloads/illumina450k_filtering/HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt', head = F, as.is = T)
multi.map.probes <- as.character(multi.map$V1)
# determine unique probes between the cross-reactive and multi-map probes
bad.probes <- unique(c(cross.react.probes, multi.map.probes))
keep <- !(featureNames(gmtSet) %in% bad.probes) #Removal of these bad probes
table(keep)
gmtSet <- gmtSet[keep,]
dim(gmtSet) #Number of probes = 392870

#Grab Betas and m_values
beta_values <- getBeta(gmtSet)
m_values <- getM(gmtSet)

#Remove probes Hypomethylated in all samples identified by CpG sites having beta < 0.05 in all samples�0.05 in all samples
beta_values_filtered <- as.data.frame(beta_values) 
dim(filter_all(beta_values_filtered, all_vars(. < 0.05)))
beta_values_filtered <- filter_all(beta_values_filtered, any_vars(. >= 0.05)) 
#6327 Hypomethylated

#Remove probes Hypermethylated in all samples identified by CpG sites having beta > 0.95 in all samples
dim(filter_all(beta_values_filtered, all_vars(. > 0.95)))
beta_values_filtered <- filter_all(beta_values_filtered, any_vars(. < 0.95)) 
dim(beta_values_filtered)
#204 hypermethylated

# #EXCLUDE CONTROL SAMPLES AND DCD SAMPLES - Liver dataset
# beta_values_case = beta_values_filtered[,!(colnames(beta_values_filtered) %in% c("200999740023_R05C02","200999740023_R06C02","201004900096_R05C02","201004900096_R06C02","202702240054_R01C01","202702240054_R02C01","202702240079_R07C01","202702240079_R08C01","3999442124_R05C02","3999442124_R06C02","203751390020_R02C01","3999442124_R01C02","201004900096_R02C02","200999740005_R06C02","201004900018_R06C01","203751390017_R07C01","200687170042_R05C02","201004900096_R03C02","200999740023_R01C01","201004900018_R01C02"))]
# pheno_df_case = pheno_df[!(rownames(pheno_df) %in% c("200999740023_R05C02","200999740023_R06C02","201004900096_R05C02","201004900096_R06C02","202702240054_R01C01","202702240054_R02C01","202702240079_R07C01","202702240079_R08C01","3999442124_R05C02","3999442124_R06C02","203751390020_R02C01","3999442124_R01C02","201004900096_R02C02","200999740005_R06C02","201004900018_R06C01","203751390017_R07C01","200687170042_R05C02","201004900096_R03C02","200999740023_R01C01","201004900018_R01C02","NA","NA.1","NA.2","NA.3","NA.4","NA.5","NA.6","NA.7")),]  


#EXCLUDE DCD SAMPLES - Kidney dataset
beta_values_case <- beta_values_filtered[,!(colnames(beta_values_filtered) %in% c("9296930129_R05C01", "9305651174_R01C01", "9305651174_R03C01", "9305651174_R02C02", "9305651174_R03C02", "9305651191_R02C02", "9305651191_R04C02", "201465900002_R04C01", "202240580106_R03C01", "202240580208_R03C01", "202259340119_R05C01", "202259350016_R04C01", "203496240002_R03C01", "203504430032_R05C01", "204001300109_R07C01", "204001300109_R08C01", "204001350016_R01C01", "202702240079_R06C01"))]
#Removing rows based on the sample_name column in phenotype dataframe
pheno_df_case <- pheno_df[!(pheno_df$sample_name %in% c("9296930129_R05C01", "9305651174_R01C01", "9305651174_R03C01", "9305651174_R02C02", "9305651174_R03C02", "9305651191_R02C02", "9305651191_R04C02", "201465900002_R04C01", "202240580106_R03C01", "202240580208_R03C01", "202259340119_R05C01", "202259350016_R04C01", "203496240002_R03C01", "203504430032_R05C01", "204001300109_R07C01", "204001300109_R08C01", "204001350016_R01C01", "202702240079_R06C01")),]

#Convert filter beta dataframe to m_value dataframe
m_values_case = beta2m(beta_values_case)


#Transpose resulting PCA plot 
transposed_beta <- t(beta_values_case)
transposed_beta <- data.frame(transposed_beta)
pca_general <- prcomp(transposed_beta, center=TRUE)
var_explained <- pca_general$sdev^2/sum(pca_general$sdev^2)
scores = as.data.frame(pca_general$x)
scores['sample_name'] <- row.names(scores)
scores <- merge(scores, pheno_df, by = 'sample_name')
plot <- ggplot(data=scores, mapping = aes(x = PC1, y = PC2)) +theme_bw() + geom_point()+ labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"), y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))
plot + geom_text(aes(label = sample_name), size=3.5) + xlim(-100, 400)
print(plot)

#Filter out the outlier sample - 203504430032_R01C01
filtered_transposed_beta <- transposed_beta[!((row.names(transposed_beta)) %in% '203504430032_R01C01'),]
pca_general <- prcomp(filtered_transposed_beta, center=TRUE)
var_explained <- pca_general$sdev^2/sum(pca_general$sdev^2)
scores = as.data.frame(pca_general$x)
scores['sample_name'] <- row.names(scores)
scores <- merge(scores, pheno_df, by = 'sample_name')
plot <- ggplot(data=scores, mapping = aes(x = PC1, y = PC2, color=time)) +theme_bw() + geom_point(aes(shape=array_type), alpha=0.5)+ labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"), y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) + scale_color_manual(values=pal)
print(plot)

#Other values you can extract
M <- getM(gmtSet)
CN <- getCN(gmtSet)

#Differential Methylation Positions 
?dmpFinder()
time <- pheno_df$time
dmp <- dmpFinder(beta_values, pheno = time, type = "categorical")
head(dmp)

#Differential Methylation - using MEAL 
rowData(gmtSet) <- getAnnotation(gmtSet)[, -c(1:3)]

#Remove probes measuring SNPs
gmtSet <- dropMethylationLoci(gmtSet)

#Remove probes with SNPs
gmtSet <- dropLociWithSnps(gmtSet)

#Remove probes with NAs for betas
gmtSet <- gmtSet[!apply(getBeta(gmtSet), 1, function(x) any(is.na(x))), ]

#Select a subset of samples
set.seed(0)

?runPipeline()
phenoData <- as.data.frame(pData(gmtSet))
rowData(gmtSet) <- getAnnotation(gmtSet)[, -c(1:3)]
set.seed(0) #reproducible results for simulations 
gmtSet <- gmtSet[sample(nrow(gmtSet), 100000), ]
featureNames(gmtSet)

table(keep)
gmtSet <- gmtSet[keep,]
dim(gmtSet)
sampleNames <- sampleNames(gmtSet)
keep <- !(gmtSet@colData$time == 'K3' | gmtSet@colData$time == 'K12')
gmtSet <- gmtSet@colData[keep,]
as.data.frame(gmtSet)
res <- runPipeline(set = gmtSet, variable_names = "time", analyses = c("DiffMean", "DiffVar"))
design <- model.matrix(~ gmtSet$time)
res
names(res)
?()
gmtSet <- gmtSet[!(gmtSet$time == "K12" | gmtSet$time == "K3")]
gmtSet$time
#Plotting 
targetRange <- GRanges("23:13000000-23000000")
plot(res, rid = "DiffMean", type = "manhattan", highlight = targetRange)
head(getAssociation(res, "DiffMean"))
