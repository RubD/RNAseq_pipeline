#!/usr/bin/env Rscript

# accept xlsx file as first argument and output directory as second
args <- commandArgs(TRUE)
xlsx_file <- args[1]
outputdir <- args[2]

# load necessary libraries
library(xlsx)
library(yaml)

# load xlsx file for preconfigfile
preconfig <- read.xlsx(file = xlsx_file, sheetIndex = 1, header = FALSE,
                       rowIndex = c(2:6,9,12:16,19,22:24,27:28, 31:33),colIndex = c(1,2), stringsAsFactors = F)



# function to create yaml file from xlsx file
createYamlConfig <- function(preconfigfile, outputdirectory) {
  
  rownames(preconfig) <- preconfig[,1]
  newpreconfig <- preconfig[,-1, drop = FALSE]
  tnewpreconfig <- as.data.frame(t(newpreconfig), stringsAsFactors = F)
  yamldfr <- as.yaml(tnewpreconfig)
  
  outputpath <- paste0(outputdirectory,"/RNAseqConfig.yaml")
 
  resultlist <- list(yamldfr, outputpath)
  return(resultlist)
}

# if output directory is not provided on command line, change it to working directory
if(is.na(outputdir)) {
  outputdir <- getwd()
}

# if RNA seq configuration file already exists, print error message and stop
if(file.exists(paste0(outputdir,"/RNAseqConfig.yaml"))) {
  print("RNAseq configuration file already exists (RNAseqConfig.yaml), delete or move this file first")
  stop()
}

myresultlist <- createYamlConfig(preconfig, outputdir)

myfinalconfigfile <- myresultlist[[1]][1]
myfinaloutputdir <- myresultlist[[2]][1]

# write out the result
write.table(myfinalconfigfile,file = myfinaloutputdir, row.names = FALSE, quote = FALSE, col.names = FALSE)


