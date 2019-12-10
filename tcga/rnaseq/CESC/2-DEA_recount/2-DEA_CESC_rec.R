### AIM ###
#download, aggregation and pre-processing of RNASEQ data#

### WORKING DIRECTORY ####
#replace with proper location #
setwd("/data/user/shared_projects/ulk1/tcga/rnaseq/CESC/2-DEA_recount/")

### LOAD LIBRARIES ####

library(TCGAbiolinks) #2.13.1
library(SummarizedExperiment)
library(devtools)
#library(TCGAutils)
library(recount)
library(biomaRt)

#check version of TCGAbiolinks 2.13.1#

### FUNCTIONS

### FUNCTIONS ####

# Get harmonized queries using GDCquery from a specific project, e.g "TCGA-CESC",
# and different sample types defined in the 'sample_types' list


####Function to convert ENSG to Symbol
convert.ENSG.Symbol<-function(genes){
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
  #merge(DEG.ucs,G_list,by.x="gene",by.y="ensembl_gene_id")
  return(G_list)
}



########Query from Recount2 platform: Uterine Carcinoma#######
CESC.recount.gtex<-TCGAquery_recount2(project="GTEX", tissue="cervix uteri")
CESC.recount.tcga<-TCGAquery_recount2(project="TCGA", tissue="cervix")

#to get the SE object
SE.CESC.recount.gtex <- CESC.recount.gtex$GTEX_cervix_uteri
SE.CESC.recount.tcga <- CESC.recount.tcga$TCGA_cervix

#####Preparing/scaling Recount2 data because it was sequenced using Rail-RNA

eset.gtex<-assays(scale_counts(CESC.recount.gtex$GTEX_cervix_uteri, round = TRUE))$counts
eset.tcga<-assays(scale_counts(CESC.recount.tcga$TCGA_cervix, round = TRUE))$counts

#### Check that the number of reads is less than or equal to 40 million
rse_scaled <- scale_counts(CESC.recount.gtex$GTEX_cervix_uteri, round = TRUE)
summary(colSums(assays(rse_scaled)$counts)) / 1e6

###replacing UUIDs with TCGA barcodes:
colnames(eset.tcga)<-colData(CESC.recount.tcga$TCGA_cervix)$gdc_cases.samples.portions.analytes.aliquots.submitter_id

###RemCESCing version (number after ".")
rownames(eset.gtex) <- gsub("\\..*", "", rownames(eset.gtex))
rownames(eset.tcga) <- gsub("\\..*", "", rownames(eset.tcga))

####Segregate between primary tumors and normal samples
eset.tcga.cancer<-eset.tcga[,which(colData(CESC.recount.tcga$TCGA_cervix)$gdc_cases.samples.sample_type=="Primary Tumor")]
#eset.tcga.normal<-eset.tcga[,which(colData(CESC.recount.tcga$TCGA_uterus)$gdc_cases.samples.sample_type=="Solid Tissue Normal")]
####
dim(eset.gtex)
dim(eset.tcga.cancer)

##merging data by row names
dataPrep2<-merge(as.data.frame(eset.gtex), as.data.frame(eset.tcga.cancer), by=0, all=TRUE)
dataPrep2 <- subset(dataPrep2, !duplicated(subset(dataPrep2, select=c(Row.names))))
rownames(dataPrep2)<-dataPrep2$Row.names
dataPrep2$Row.names<-NULL

dim(dataPrep2)
# [1] 58037 315


### PRE-PROCESSING ###


#save(dataPrep2, file = "CESC_PreprocessedData_all_Prep_11062019.rda")

##Normalization for library size and for gcContent (it can be also done for gene lenght in case needed but they are mutally exclusive)

dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep2,
                                             geneInfo = geneInfoHT,
                                             method = "gcContent")
dim(dataNorm)
#[1] 23258 315

cat ("These ensembl IDs was found in the dataPrep file:", "\n",ensemblIDs[which(ensemblIDs %in% rownames(dataNorm))])

#save(dataNorm, file= "CESC_PreprocessedData_all_Norm_11062019.rda")

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                         method = "quantile",
                                         qnt.cut =  0.20)
dim(dataFilt)

##[1] 18604  315

cat ("These ensembl IDs was found in the dataPrep file:", "\n",ensemblIDs[which(ensemblIDs %in% rownames(dataFilt))])

#save(dataFilt, file = "CESC_PreprocessedData_all_Filt_11062019.rda")


data.DEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,colnames(eset.gtex)],
                             mat2 = dataFilt[,colnames(eset.tcga.cancer)],
                             metadata = FALSE,
                             pipeline = "limma",
                             voom = TRUE,
                             Cond1type = "Normal",
                             Cond2type = "Tumor",
                             fdr.cut = 0.05,
                             logFC.cut = 0.5)

#11890
write.table(data.DEGs, file = "DEA_rec_CESC.txt", quote = FALSE)

dataDEGs.ULK1 <- subset(data.DEGs, rownames(data.DEGs) %in% "ENSG00000177169")
dataDEGs.ULK2 <- subset(data.DEGs, rownames(data.DEGs) %in% "ENSG00000083290")
dataDEGs.ULKs <- rbind(dataDEGs.ULK1, dataDEGs.ULK2)

### to check on specific DE genes with boxplots

## log transformation - we could also use the voom limma transformation instead
## select only the matrix including the genes of interest

dataFilt_normal <- dataFilt[,colnames(eset.gtex)]
dataFilt_normal <- log2(dataFilt_normal+1)
dataFilt_cancer <- dataFilt[,colnames(eset.tcga.cancer)]
dataFilt_cancer <- log2(dataFilt_cancer+1)
dataFilt_normal_ULK1 <- (dataFilt_normal[rownames(dataFilt_normal) %in% "ENSG00000177169",])
dataFilt_cancer_ULK1 <- (dataFilt_cancer[rownames(dataFilt_cancer) %in% "ENSG00000177169",])
dataFilt_normal_ULK2 <- (dataFilt_normal[rownames(dataFilt_normal) %in% "ENSG00000083290",])
dataFilt_cancer_ULK2 <- (dataFilt_cancer[rownames(dataFilt_cancer) %in% "ENSG00000083290",])


pdf(file= "boxplot.CESC.ULK1.pdf")
boxplot(dataFilt_normal_ULK1, dataFilt_cancer_ULK1)
dev.off()

pdf(file= "boxplot.CESC.ULK2.pdf")
boxplot(dataFilt_normal_ULK2, dataFilt_cancer_ULK2)
dev.off()

