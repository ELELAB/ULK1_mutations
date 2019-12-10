### AIM ###
#DEA on RNASEQ data#

### WORKING DIRECTORY ####
#replace with proper location #
setwd("/data/user/shared_projects/ulk1/tcga/rnaseq/BRCA/2-DEA/")

### LOAD LIBRARIES ####

library(TCGAbiolinks)
library(SummarizedExperiment)
library(limma)
library(dplyr)
library(ggplot2)
library(reshape2)

### FUNCTIONS ####

# Get harmonized queries using GDCquery from a specific project, e.g "TCGA-BRCA",
# and different sample types defined in the 'sample_types' list

get_gdcquery <- function(project, sample_types){
  query      <- GDCquery(project = project,
                         data.category = "Transcriptome Profiling",
                         data.type = "Gene Expression Quantification",
                         workflow.type = "HTSeq - Counts",
                         sample.type = sample_types,
                         legacy = FALSE)
  return(query)
}


#### VARIABLES ####

## Cancer type to analyze gene expression
project <- "TCGA-BRCA"

## Which sample types to perform the gene expression analysis on (short letter code)
sample_types_short <- c("TP", "NT")

## Which sample types to perform the gene expression analysis on (full names)
## these must correspond to the sample types in 'sample_types_short'
sample_types_full <- c("Primary solid Tumor", "Solid Tissue Normal")

#recall query

query <- get_gdcquery(project, sample_types_full)

## which samples are solid primary tumor
dataSmTP_p <- TCGAquery_SampleTypes(getResults(query,cols="cases"),"TP")
## which samples are solid tissue normal
dataSmNT_p <- TCGAquery_SampleTypes(getResults(query,cols="cases"),"NT")

#load DataFilt

DataPrep1 <- get(load('../1-download_preprocessing/BRCA_PreprocessedData_all_Prep_11062019.rda'))
DataFilt <- get(load('../1-download_preprocessing/BRCA_PreprocessedData_all_Filt_11062019.rda'))

dataDEGs.TSS <- TCGAanalyze_DEA(mat1 = dataFilt[,dataSmNT_p],
                            mat2 = dataFilt[,dataSmTP_p],
                            pipeline = "limma",
                            batch.factors = "TSS",
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.05,
                            logFC.cut = 0.5,
                            voom = TRUE)

write.table(dataDEGs.TSS, file = "DEA_TSS_BRCA.txt", quote = FALSE)

#9001

dataDEGs.Plate <- TCGAanalyze_DEA(mat1 = dataFilt[,dataSmNT_p],
                                mat2 = dataFilt[,dataSmTP_p],
                                pipeline = "limma",
                                batch.factors = "Plate",
                                Cond1type = "Normal",
                                Cond2type = "Tumor",
                                fdr.cut = 0.05,
                                logFC.cut = 0.5,
                                voom = TRUE)
#9391
write.table(dataDEGs.Plate, file = "DEA_Plate_BRCA.txt", quote = FALSE)

dataDEGs.Center <- TCGAanalyze_DEA(mat1 = dataFilt[,dataSmNT_p],
                                mat2 = dataFilt[,dataSmTP_p],
                                pipeline = "limma",
                                batch.factors = "Center",
                                Cond1type = "Normal",
                                Cond2type = "Tumor",
                                fdr.cut = 0.05,
                                logFC.cut = 0.5,
                                voom = TRUE)
#Not available - same sequencing center

write.table(dataDEGs.Center, file = "DEA_Center_BRCA.txt", quote = FALSE)
#Not available


dataDEGs.TSS.ULK1 <- subset(dataDEGs.TSS, rownames(dataDEGs.TSS) %in% "ENSG00000177169")
dataDEGs.TSS.ULK2 <- subset(dataDEGs.TSS, rownames(dataDEGs.TSS) %in% "ENSG00000083290")
dataDEGs.TSS.ULKs <- rbind(dataDEGs.TSS.ULK1, dataDEGs.TSS.ULK2)

#Not DE

dataDEGs.Plate.ULK1 <- subset(dataDEGs.Plate, rownames(dataDEGs.Plate) %in% "ENSG00000177169")
dataDEGs.Plate.ULK2 <- subset(dataDEGs.Plate, rownames(dataDEGs.Plate) %in% "ENSG00000083290")
dataDEGs.Plate.ULKs <- rbind(dataDEGs.Plate.ULK1, dataDEGs.Plate.ULK2)
#Not DE

dataDEGs.Center.ULK1 <- subset(dataDEGs.Center, rownames(dataDEGs.Center) %in% "ENSG00000177169")
dataDEGs.Center.ULK2 <- subset(dataDEGs.Center, rownames(dataDEGs.Center) %in% "ENSG00000083290")
dataDEGs.Center.ULKs <- rbind(dataDEGs.Center.ULK1, dataDEGs.Center.ULK2)

#Not available 

### to check on specific DE genes with boxplots

## log transformation - we could also use the voom limma transformation instead
## select only the matrix including the genes of interest

ensemblIDs <- c("ENSG00000177169", #ULK1
                "ENSG00000083290" #ULK2
)

dataFilt_sub <- dataFilt[rownames(dataFilt) %in% ensemblIDs,]
## log-transformation. add 1 to avoid infinity values
dataFilt_sub <- log2(dataFilt_sub+1)

# prepare data for ggplot to make a boxplot
data <- melt(dataFilt_sub)

sampleTP <- TCGAquery_SampleTypes(data$Var2,'TP')
sampleNT <- TCGAquery_SampleTypes(data$Var2,'NT')
condition <- c()
condition[which(data$Var2 %in% sampleTP)] <- 'cancer'
condition[which(data$Var2 %in% sampleNT)] <- 'normal'
data <- cbind(data,condition)
colnames(data)[1:2] <- c('gene','barcodes')

pdf(file = "boxplot-DEA-BRCA.pdf")
p<-ggplot(data, aes(x=gene, y=value, color=condition)) +
  geom_boxplot()
p
dev.off()
