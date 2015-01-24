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

commandtest <- paste("do_Tophat2_countFeat2.sh", inputvec[1], inputvec[2], inputvec[3], inputvec[4], inputvec[5]
                     ,inputvec[6], inputvec[7], inputvec[8], inputvec[9], inputvec[10], inputvec[11], inputvec[12],
                     inputvec[13], inputvec[14], inputvec[15], inputvec[16], inputvec[17])

system(commandtest)

