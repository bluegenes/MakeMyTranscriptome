#!/usr/bin/Rscript

#####################################################
## Run DESeq2 Analysis on RNASeq datasets 
## Tessa Pierce
## 7.31.2015
#####################################################

# get args from bash
args <- commandArgs(trailingOnly = TRUE)
countTableFile = args[2]
treatmentInfoFile = args[3]
baseDir = args[4]

## read in metadata ## 
treatmentInfo <- read.table(treatmentInfoFile, header=TRUE, row.names = 1, stringsAsFactors =TRUE)
## set model from command line OR from order of factors in the treatmentInfoFile
if (length(args) >= 6){ 
    modelParam = paste("~",paste((args[6:length(args)]), collapse = "+")) 
    baseName = paste(args[5],paste((args[6:length(args)]), collapse = "_"),sep='_')
}else{
    modelParam = paste("~",paste(rev(colnames(treatmentInfo)), collapse = "+")) 
    baseName = paste(args[5],paste(rev(colnames(treatmentInfo)), collapse = "_"),sep='_') # reverse order, so 1st/most impt = last in model
    }  

# set up the R environment
source('autoDESeq2functions.r') #wrapper functions necessary for this script to run
setwd(baseDir)# set output directory 

## read in count data ##  
countTableRaw = read.table(countTableFile, header=TRUE, row.names = 1)
countTable <-as.data.frame(lapply(countTableRaw,as.integer)) # eXpress results are floats ...coerce these to integer -> need ints for DESeq2
row.names(countTable) <- row.names(countTableRaw) #add gene names back to integer dataframe
colnames(countTable) <- gsub(pattern = ".bed", replacement = "", x = colnames(countTable))
colnames(countTable) <- gsub(pattern = ".xprs", replacement = "", x = colnames(countTable))

# subset counts based on filenames, factors provided in the treatmentInfo file
subsetCounts <- countTable[,colnames(countTable)%in%rownames(treatmentInfo)]

## run DESeq2: all pairwise contrasts ## 
dds <- DESeq(DESeqDataSetFromMatrix(countData = subsetCounts, colData = treatmentInfo, design= ~ 1))
design(dds) <- formula(modelParam)
dds <- DESeq(dds)
runAllDESeq2(baseDir,treatmentInfo,baseName, dds,.1)

