### AIM ###
#DEA on RNASEQ data without tcgabiolinks functions#

### WORKING DIRECTORY ####
#replace with proper location #
setwd("/data/user/shared_projects/ulk1/tcga/rnaseq/ESCA/2-DEA/")

### LOAD LIBRARIES ####

library(TCGAbiolinks)
library(SummarizedExperiment)
library(limma)
library(dplyr)
library(ggplot2)
library(reshape2)

### FUNCTIONS ####

source("tcgabiolinks_functions.R")


limma_tss <- function(my_IDs,dataframe,limma_name,up_name,down_name){
  
  condition <- as.factor(my_IDs$condition)
  tss <- as.factor(my_IDs$tss)
  
  #design matrix
  design.matrix <- model.matrix(~0+condition+tss)
  colnames(design.matrix)[c(1,2)] <- c("cancer","normal")
  
  #voom transformation
  dataframe <- voom(dataframe,design.matrix,plot=TRUE)
  
  # Making group contrasts 
  N_C_cotr <- makeContrasts("cancer-normal", levels= design.matrix)
  
  # Filter for significance - set log fold change and fdr cutoff
  N_C_L <- DE_limma(N_C_cotr, dataframe, design.matrix, 0.5, 0.05)
  
  # differentially expressed genes file
  write.csv(N_C_L, limma_name, quote = FALSE)
  
  #number of up-regulated genes cancer vs normal
  print(paste0("the analysis detected ",length(which(N_C_L$direction == "up"))," up-regulated genes"))
  #number of down-regulated genes cancer vs normal
  print(paste0("the analysis detected ",length(which(N_C_L$direction== "down")), " down-regulated genes"))
  
  # up and down regulated genes
  up <- data.frame(rownames(N_C_L[N_C_L$direction == "up", ]))
  down <- data.frame(rownames(N_C_L[N_C_L$direction == "down", ]))
  
  write.table(up, up_name, sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  write.table(down, down_name, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  
}



#### VARIABLES ####


dataframe_ESCA <- get(load("../1-download_preprocessing/ESCA_PreprocessedData_all_Filt_08072019.rda"))
my_IDs <- get_IDs(dataframe_ESCA)
length(which(my_IDs$condition=="cancer"))
length(which(my_IDs$condition=="normal"))
limma_name <- "limma_ESCA_all_tss.csv"
up_name <- "up_limma_ESCA_tss.txt"
down_name <- "down_limma_ESCA_tss.txt"

limma_tss(my_IDs,dataframe_ESCA,limma_name,up_name,down_name)

limma_ESCA <- read.csv("limma_ESCA_all_tss.csv")


#5131 up / 4149 down

limma.TSS.ULK1 <- subset(limma_ESCA, X =="ENSG00000177169")
limma.TSS.ULK2 <- subset(limma_ESCA, X== "ENSG00000083290")
limma.TSS.ULKs <- rbind(limma.TSS.ULK1, limma.TSS.ULK2)

### EDGE-R pipeline - functions ####

#DIFFERENTIAL EXPRESSION ANALYSIS - edgeR

edgeR_tss <- function(dataframe,edgeR_name,up_name,down_name){
  
  my_IDs <- get_IDs(dataframe)
  # edgeR object
  y <- DGEList(counts=dataframe)
  
  # model matrix info
  y$samples$condition <- as.factor(my_IDs$condition)
  y$samples$tss <- as.factor(my_IDs$tss)
  
  # relevel data
  y$samples$condition = relevel(y$samples$condition, ref="normal")
  
  # design matrix
  design.mat <- model.matrix(~condition+tss, data=y$samples)
  
  # estimate dispersion
  y <- estimateDisp(y,design.mat)
  
  # fit overdispersed poisson model
  my.fit <- glmFit(y, design.mat)
  
  # Performing likelihood ratio test
  N_C_coef <- glmLRT(my.fit, coef = 2)
  
  # Filter for significance - set log fold change and fdr cutoff
  N_C_E <- DE_edgeR(N_C_coef, y, 0.5, 0.05)
  
  # differentially expressed genes file
  
  write.csv(N_C_E, edgeR_name, quote = FALSE)
  
  #number of up-regulated genes cancer vs normal
  print(paste0("the analysis detected ",length(which(N_C_E$direction == "up"))," up-regulated genes"))
  #number of down-regulated genes cancer vs normal
  print(paste0("the analysis detected ",length(which(N_C_E$direction == "down"))," down-regulated genes"))
  
  # up and down regulated genes
  up <- data.frame(rownames(N_C_E[N_C_E$direction == "up", ]))
  down <- data.frame(rownames(N_C_E[N_C_E$direction == "down", ]))
  
  write.table(up, up_name, sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
  write.table(down,down_name, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
  
}


### EdgeR pipeline DEA ####

edgeR_name <- "edgeR_ESCA_all.csv"
up_name <- "up_edgeR_ESCA.txt"
down_name <- "down_edgeR_ESCA.txt"

edgeR_tss(dataframe_ESCA,edgeR_name,up_name,down_name)

#5939 up / 3837 down

edgeR.ESCA <- read.csv("edgeR_ESCA_all.csv")

ensemblIDs <- c("ENSG00000177169",  "ENSG00000083290" )
                
edgeR.TSS.ULK1 <- subset(edgeR.ESCA, X =="ENSG00000177169")
edgeR.TSS.ULK2 <- subset(edgeR.ESCA, X== "ENSG00000083290")
edgeR.TSS.ULKs <- rbind(edgeR.TSS.ULK1, edgeR.TSS.ULK2)

                
