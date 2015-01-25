#!/usr/bin/Rscript

args <- commandArgs(T)
config.file <- args[1]

#source("http://bioconductor.org/biocLite.R")
#biocLite(c("DESeq2"))

library(DESeq2)

# load the config yaml to read the config of the pipeline
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


featC_matrix_path <- paste0(resultsOutputDir,"featc_matrix.txt")

if(countTableOrigin == "featureCounts" & file.exists(featC_matrix_path)) {
  
  FCdfr <- read.table(file = featC_matrix_path, header = T)
}
 


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

myMatr <- createCountMatrix(FCdfr = FCdfr, ctrlTag = CtrlTag, treatTag = TreatmentTag)

condition <- ifelse(grepl("*Ctrl*", colnames(myMatr)), "Control","Treatment")
myrownames <- colnames(myMatr)

colFCmat <- data.frame(condition = condition, row.names = myrownames)



# create DESeqDataSet based on countmatrix and dataframe for column information
library(DESeq2)
ddsFCmat <- DESeqDataSetFromMatrix(countData = myMatr,
                                 colData = colFCmat,
                                 design = ~ condition)

## normalization of raw counts (for DEG)
ddsFCmat <- estimateSizeFactors(ddsFCmat)
sizeFactors(ddsFCmat)
normddsFCmat <- counts(ddsFCmat, normalized = TRUE)

####################################
# DEG between different conditions #
####################################
####################################

# annotate with entrezid
dds_tempd <- DESeq(ddsFCmat) # actual DEG test "nBinom test"
res <- results(dds_tempd)
#class(res)

# set working directory and create folder if necessary

DESeqDir <- paste0(resultsOutputDir,"/FeatC_DESeq2/")

if(!file.exists(DESeqDir)) {
  dir.create(DESeqDir)
}

resfile <- paste0(DESeqDir,"FC_DESeq2_results.txt")
print(resfile)

write.table(as.matrix(res),file = resfile, sep="\t", row.names = FALSE)
            
          
