### AIM ###
#download, aggregation and pre-processing of RNASEQ data#

### WORKING DIRECTORY ####
#replace with proper location #
setwd("/data/user/shared_projects/ulk1/tcga/rnaseq/DLBC/1-download_preprocessing_recount/")

### LOAD LIBRARIES ####

library(TCGAbiolinks) #2.13.1
library(SummarizedExperiment)
library(devtools)
#library(TCGAutils)
library(recount)
library(biomaRt)

#check version of TCGAbiolinks 2.9.5#

### FUNCTIONS

### FUNCTIONS ####

# Get harmonized queries using GDCquery from a specific project, e.g "TCGA-DLBC",
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

####Function to convert ENSG to Symbol
convert.ENSG.Symbol<-function(genes){
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes,mart= mart)
  #merge(DEG.ucs,G_list,by.x="gene",by.y="ensembl_gene_id")
  return(G_list)
}
#### VARIABLES ####

## Cancer type to analyze gene expression
project <- "TCGA-DLBC"

## Which sample types to perform the gene expression analysis on (short letter code)
sample_types_short <- c("TP", "NT")

## Which sample types to perform the gene expression analysis on (full names)
## these must correspond to the sample types in 'sample_types_short'
sample_types_full <- c("Primary solid Tumor", "Solid Tissue Normal")

########Query from Recount2 platform:#######
DLBC.recount.gtex<-TCGAquery_recount2(project="GTEX", tissue="blood")
DLBC.recount.tcga<-TCGAquery_recount2(project="TCGA", tissue="lymph_nodes")

#to get the SE object
SE.DLBC.recount.gtex <- DLBC.recount.gtex$GTEX_blood
SE.DLBC.recount.tcga <- DLBC.recount.tcga$TCGA_lymph_nodes

###   DOWNLOAD AND AGGREGATE THE DATA ###
## Get GDC queries using the function abDLBCe
## project = "TCGA-DLBC"
## sample_types_full = list of sample types to use for further analysis (full names)

query <- get_gdcquery(project, sample_types_full)

#### The steps below are needed to have the right correspondance beetween barcodes (TCGA) and UUID (recount)
#download data -> they will be saved in a GDCdata subfolder
#GDCdownload(query = query)


samplesDown <- getResults(query,cols=c("cases"))


###tumor samples for uterine cancer
dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                      typesample = "TP")

###to check that there are no NT samples
dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
                                      typesample = "NT")

#####Preparing/scaling Recount2 data because it was sequenced using Rail-RNA

eset.gtex<-assays(scale_counts(DLBC.recount.gtex$GTEX_blood, round = TRUE))$counts
eset.tcga<-assays(scale_counts(DLBC.recount.tcga$TCGA_lymph_nodes, round = TRUE))$counts

#### Check that the number of reads is less than or equal to 40 million
rse_scaled <- scale_counts(DLBC.recount.gtex$GTEX_blood, round = TRUE)
summary(colSums(assays(rse_scaled)$counts)) / 1e6

###replacing UUIDs with TCGA barcodes:
colnames(eset.tcga)<-colData(DLBC.recount.tcga$TCGA_lymph_nodes)$gdc_cases.samples.portions.analytes.aliquots.submitter_id

###RemDLBCing version (number after ".")
rownames(eset.gtex) <- gsub("\\..*", "", rownames(eset.gtex))
rownames(eset.tcga) <- gsub("\\..*", "", rownames(eset.tcga))

####Segregate between primary tumors and normal samples
eset.tcga.cancer<-eset.tcga[,which(colData(DLBC.recount.tcga$TCGA_lymph_nodes)$gdc_cases.samples.sample_type=="Primary Tumor")]
#eset.tcga.normal<-eset.tcga[,which(colData(DLBC.recount.tcga$TCGA_DLBCary)$gdc_cases.samples.sample_type=="Solid Tissue Normal")]
####
dim(eset.gtex)
dim(eset.tcga.cancer)

##merging data by row names
dataPrep2<-merge(as.data.frame(eset.gtex), as.data.frame(eset.tcga.cancer), by=0, all=TRUE)
dataPrep2 <- subset(dataPrep2, !duplicated(subset(dataPrep2, select=c(Row.names))))
rownames(dataPrep2)<-dataPrep2$Row.names
dataPrep2$Row.names<-NULL

dim(dataPrep2)
# [1] 58037  643


### PRE-PROCESSING ###

##if we are interested to a specific group of genes to check at each step if they are there
ensemblIDs <- c("ENSG00000177169", #ULK1
                "ENSG00000083290" #ULK2
)

genes_of_interest <- c("ULK1", "ULK2")

##search outliers - only remDLBCing samples not genes

#dataPrep2a <- TCGAanalyze_Preprocessing(object = dataPrep2, cor.cut = 0.6, filename = "DLBC_dataprep_07062019.png")
#dim(dataPrep2a)
#56734   1201

save(dataPrep2, file = "DLBC_PreprocessedData_all_Prep_11062019.rda")

##Normalization for library size and for gcContent (it can be also done for gene lenght in case needed but they are mutally exclusive)

dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep2,
                                             geneInfo = geneInfoHT,
                                             method = "gcContent")
dim(dataNorm)
#[1] 23258 643

cat ("These ensembl IDs was found in the dataPrep file:", "\n",ensemblIDs[which(ensemblIDs %in% rownames(dataNorm))])

save(dataNorm, file= "DLBC_PreprocessedData_all_Norm_11062019.rda")

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                         method = "quantile",
                                         qnt.cut =  0.20)
dim(dataFilt)

##[1] 18606  643

cat ("These ensembl IDs was found in the dataPrep file:", "\n",ensemblIDs[which(ensemblIDs %in% rownames(dataFilt))])

save(dataFilt, file = "DLBC_PreprocessedData_all_Filt_07062019.rda")


