#!/usr/bin/Rscript

args <- commandArgs(TRUE)
config.file <- args[1]

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


inputvec <- c(rawFastqDir, resultsOutputDir, bowtie2indexDir, countFeatAlg, genomeAnnotation, fastQCthreads, 
              qualityFastqDir, tophat2threads, segmentMismatch, segmentLength, libType, paralleljobs,
              FCfeatureType, FCmetaType, FCthreads, HTSfeatureType, HTSmetaType)

commandtest <- paste("do_fastqc.sh", inputvec[1], inputvec[2], inputvec[6])

system(commandtest)

