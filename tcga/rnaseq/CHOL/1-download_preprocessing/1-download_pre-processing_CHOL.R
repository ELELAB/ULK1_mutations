### AIM ###
#download, aggregation and pre-processing of RNASEQ data#

### WORKING DIRECTORY ####
#replace with proper location #
setwd("/data/user/shared_projects/ulk1/tcga/rnaseq/CHOL/1-download_preprocessing/")

### LOAD LIBRARIES ####

library(TCGAbiolinks)
library(SummarizedExperiment)


#check version of TCGAbiolinks 2.9.5#

### FUNCTIONS

### FUNCTIONS ####

# Get harmonized queries using GDCquery from a specific project, e.g "TCGA-CHOL",
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
project <- "TCGA-CHOL"

## Which sample types to perform the gene expression analysis on (short letter code)
sample_types_short <- c("TP", "NT")

## Which sample types to perform the gene expression analysis on (full names)
## these must correspond to the sample types in 'sample_types_short'
sample_types_full <- c("Primary solid Tumor", "Solid Tissue Normal")


###   DOWNLOAD AND AGGREGATE THE DATA ###
## Get GDC queries using the function above
## project = "TCGA-CHOL"
## sample_types_full = list of sample types to use for further analysis (full names)

query <- get_gdcquery(project, sample_types_full)

#download data -> they will be saved in a GDCdata subfolder
GDCdownload(query = query)

###preparation of the whole SE with  TP and NT
dataPrep1 <- GDCprepare(query = query,
                               save = TRUE,
                               save.filename = "TCGA_CHOL_HTSeq_Counts_03062019.rda")
dim(dataPrep1)
#[1] 56602   161

length(which(colData(dataPrep1)$shortLetterCode =="TP"))
#156
length(which(colData(dataPrep1)$shortLetterCode =="NT"))
#5

## which samples are solid primary tumor
dataSmTP_p <- TCGAquery_SampleTypes(getResults(query,cols="cases"),"TP")
## which samples are solid tissue normal
dataSmNT_p <- TCGAquery_SampleTypes(getResults(query,cols="cases"),"NT")

##if you reopen this at a later stage and you need to re-load the dataPrep1_paired object - you can uncomment and start from here
#prep_harmo <- get(load('TCGA_CHOL_HTSeq_Counts_03062019.rda'))

### PRE-PROCESSING ###

##if we are interested to a specific group of genes to check at each step if they are there
ensemblIDs <- c("ENSG00000177169", #ULK1
                "ENSG00000083290" #ULK2
)

genes_of_interest <- c("ULK1", "ULK2")

##search outliers - only removing samples not genes

dataPrep2 <- TCGAanalyze_Preprocessing(object = dataPrep1, cor.cut = 0.6, filename = "CHOL_dataprep_03062019.png")
dim(dataPrep2)
#56602   161

# check genes of interest



save(dataPrep2, file = "CHOL_PreprocessedData_all_Prep_03062019.rda")

##Normalization for library size and for gcContent (it can be also done for gene lenght in case needed but they are mutally exclusive)

dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep2,
                                             geneInfo = geneInfoHT,
                                             method = "gcContent")
dim(dataNorm)
#[1] 23192   161

cat ("These ensembl IDs was found in the dataPrep file:", "\n",ensemblIDs[which(ensemblIDs %in% rownames(dataNorm))])

save(dataNorm, file= "CHOL_PreprocessedData_all_Norm_03062019.rda")

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                         method = "quantile",
                                         qnt.cut =  0.25)
dim(dataFilt)

##[1] 17392   161

cat ("These ensembl IDs was found in the dataPrep file:", "\n",ensemblIDs[which(ensemblIDs %in% rownames(dataFilt))])

save(dataFilt, file = "CHOL_PreprocessedData_all_Filt_03062019.rda")
