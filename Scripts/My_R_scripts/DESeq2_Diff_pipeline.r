#!/usr/bin/Rscript

args <- commandArgs(T)
config.file <- args[1]

#source("http://bioconductor.org/biocLite.R")
#biocLite(c("DESeq2"))

library(DESeq2)

# load the yaml configuration file of the pipeline
library(yaml)
myconfig <- yaml.load_file(config.file) 

rawFastqDir  <-       myconfig[["rawFastqDir"]]
resultsOutputDir <-   myconfig[["resultsDir"]]
bowtie2indexDir <-    myconfig[["BowtieIndexPath"]]
countFeatAlg <-       myconfig[["featureCountMethod"]]
genomeAnnotation <-   myconfig[["AnnotationPath"]]

fastQCthreads <-      myconfig[["fastQCthreads"]]

qualityFastqDir <-    myconfig[["qualityFastqDir"]]
tophat2threads <-     myconfig[["tophat2threads"]]
segmentMismatch <-    myconfig[["segmentMismatch"]]
segmentLength <-      myconfig[["segmentLength"]]
libType <-            myconfig[["libType"]]

paralleljobs <-       myconfig[["paralleljobs"]]

FCfeatureType <-      myconfig[["FCfeatureType"]]
FCmetaType <-         myconfig[["FCmetaType"]]
FCthreads <-          myconfig[["FCthreads"]]
HTSfeatureType <-     myconfig[["HTSfeatureType"]]
HTSmetaType <-        myconfig[["HTSmetaType"]]

countTableOrigin <-   myconfig[["countTableOrigin"]]
CtrlTag <-            myconfig[["CtrlTag"]]
TreatmentTag <-       myconfig[["TreatmentTag"]]

#---------------------------------------------------------------------------------------------------------------------#

## DESeq2 for featCounts counttable
# create path to the featCounts count matrix

if(countTableOrigin == "featureCounts") {
  
  print("FeatureCounts counttable chosen for DESeq2 analysis")
  featC_matrix_path <- paste0(resultsOutputDir,"featc_matrix.txt")
  
  if(file.exists(featC_matrix_path)) {
    
    FCdfr <- read.table(file = featC_matrix_path, header = T)
    
    # function to create DESeq2 count matrix from featureCounts matrix
    createCountMatrix <- function(FCdfr, ctrlTag, treatTag) {
      
      pctrlTag <- paste0("*",ctrlTag,"*")
      ptreatTag <- paste0("*",treatTag,"*")
      
      nr_of_ctrls <- length(grep(pctrlTag, colnames(FCdfr[,-1])))
      nr_of_treat <- length(grep(ptreatTag, colnames(FCdfr[,-1])))
      
      newCtrlnames <- NULL
      for(i in 1:nr_of_ctrls) { newCtrlnames <- c(newCtrlnames, paste0("Ctrl","_rep",i)) }
      
      newTreatnames <- NULL
      for(i in 1:nr_of_treat) {  newTreatnames <- c(newTreatnames, paste0("Treat","_rep",i)) }
      
      colnames(FCdfr)[grepl(pctrlTag, colnames(FCdfr))] <- newCtrlnames
      colnames(FCdfr)[grepl(ptreatTag, colnames(FCdfr))] <- newTreatnames
      
      rownames(FCdfr) <- FCdfr[,1]
      
      FCmat <- FCdfr[,-1]
      FCmat <- as.matrix(FCmat, drop = F)
      
      
      FCmat
      
    }
    
    # create matrix
    myMatr <- createCountMatrix(FCdfr = FCdfr, ctrlTag = CtrlTag, treatTag = TreatmentTag)
    
    # create dataframe with metadata
    condition <- ifelse(grepl("*Ctrl*", colnames(myMatr)), "Control","Treatment")
    myrownames <- colnames(myMatr)
    colFCmat <- data.frame(condition = condition, row.names = myrownames)
    
    
    # create DESeqDataSet based on countmatrix and dataframe with metadata
    library(DESeq2)
    ddsFCmat <- DESeqDataSetFromMatrix(countData = myMatr,
                                       colData = colFCmat,
                                       design = ~ condition)
    
    ## normalization of raw counts (for DEG)
    ddsFCmat <- estimateSizeFactors(ddsFCmat)
    sizeFactors(ddsFCmat)
    normddsFCmat <- counts(ddsFCmat, normalized = TRUE)
    
    # Differentially expressed genes between Control and Treatment #
    ################################################################
    
    dds_tempd <- DESeq(ddsFCmat) 
    res <- results(dds_tempd)
    
    # create outputdir for DESeq2 results
    DESeqDir <- paste0(resultsOutputDir,"/FeatC_DESeq2/")
    
    if(!file.exists(DESeqDir)) {
      dir.create(DESeqDir)
    }
    
    resfile <- paste0(DESeqDir,"FC_DESeq2_results.txt")
    
    # save DESeq2 results
    write.table(as.matrix(res),file = resfile, sep="\t", row.names = TRUE)
    
  }
  else print("featureCounts matrix is not present")
}

 

## DESeq2 for HTSeq counttable
# create path to the HTSeq output folder
if(countTableOrigin == "HTSeq") {
  
  print("HTSeq counttable chosen for DESeq2 analysis")
  HTSeq_directory <- paste0(resultsOutputDir,"OutputHTSeq/")
  
  if(file.exists(HTSeq_directory)) {
    
    print("OutputHTSeq directory exists, start with DESeq2 analysis")
    
    # READ IN and ASSEMBLE HTSeq count dataframe
    myHTseqlist <- list.files(HTSeq_directory)
    
    createHTcountDfr <- function(myHTseqlist) {
      
      newlist <- list()
      
      for(i in 1:length(myHTseqlist)) {
        
        # read HTseq counts from all samples
        tempdfr <- read.table(paste0(HTSeq_directory,"/",myHTseqlist[i]), row.names = 1)
        
        myrownames <- rownames(tempdfr)
        # append dataframes to list
        newlist[i] <- tempdfr
        
      }
      
      # cbind all the dataframes and give appropriate column and row names
      finaldfr <- do.call("cbind",newlist)
      rownames(finaldfr) <- myrownames
      
      finaldfr <- as.data.frame(finaldfr)
      
      return(finaldfr)
    }
    
    HTcountdata <- createHTcountDfr(myHTseqlist = myHTseqlist)
    
    
    pctrlTag <- paste0("*","Ctrl","*")
    ptreatTag <- paste0("*","KO","*")
    
    nr_of_ctrls <- length(grep(pctrlTag, myHTseqlist))
    nr_of_treat <- length(grep(ptreatTag, myHTseqlist))
    
    newCtrlnames <- NULL
    for(i in 1:nr_of_ctrls) { newCtrlnames <- c(newCtrlnames, paste0("Ctrl","_rep",i)) }
    
    newTreatnames <- NULL
    for(i in 1:nr_of_treat) {  newTreatnames <- c(newTreatnames, paste0("Treat","_rep",i)) }
    
    myHTseqlist[grepl(pctrlTag, myHTseqlist)] <- newCtrlnames
    myHTseqlist[grepl(ptreatTag, myHTseqlist)] <- newTreatnames
    
    colnames(HTcountdata) <- myHTseqlist
    
    HTSmat <- as.matrix(HTcountdata, drop = F)
    
    # create dataframe with metadata
    condition <- ifelse(grepl("*Ctrl*", colnames(HTSmat)), "Control","Treatment")
    myrownames <- colnames(HTSmat)
    colHTSmat <- data.frame(condition = condition, row.names = myrownames)
    
    
    # create DESeqDataSet based on countmatrix and dataframe with metadata
    library(DESeq2)
    ddsHTSmat <- DESeqDataSetFromMatrix(countData = HTSmat,
                                        colData = colHTSmat,
                                        design = ~ condition)
    
    ## normalization of raw counts (for DEG)
    ddsHTSmat <- estimateSizeFactors(ddsHTSmat)
    sizeFactors(ddsHTSmat)
    normddsHTSmat <- counts(ddsHTSmat, normalized = TRUE)
    
    # Differentially expressed genes between Control and Treatment #
    ################################################################
    
    dds_tempd <- DESeq(ddsHTSmat) 
    res <- results(dds_tempd)
    
    # create outputdir for DESeq2 results
    DESeqDir <- paste0(resultsOutputDir,"/HTSeq_DESeq2/")
    
    if(!file.exists(DESeqDir)) {
      dir.create(DESeqDir)
    }
    
    resfile <- paste0(DESeqDir,"HTSeq_DESeq2_results.txt")
    
    # save DESeq2 results
    write.table(as.matrix(res),file = resfile, sep="\t", row.names = TRUE)
    
  }  
  else print("OutputHTSeq directory is not present")
}





