#!/bin/bash

args=("$@")

# path to the directory of the raw fastq files
rawFastqDir=${args[0]}

# path to the directory where the FastQC output directory must be created
resultsOutputDir=${args[1]} #fastQCoutputDir

# number of threads for fastqc command

fastQCthreads=${args[2]}

# check whether fastqc has already been performed
if [ -d $resultsOutputDir/OutputFastQC ]
then
    echo "FastQC folder is already available, check this folder or delete for new analysis"
    exit
else
    mkdir $resultsOutputDir/OutputFastQC
fi

# set outputfolder
outputfolder=$resultsOutputDir/OutputFastQC



for file in $(ls $rawFastqDir | grep "fastq") #only selects fastq files
do
    case ${file} in
	*.tar.gz)
	    echo "this $file should be first be unpacked and gunzipped"
	    tar -zxvf $rawFastqDir$file -C $rawFastqDir 
	    fastqc $rawFastqDir/*.fastq --threads $fastQCthreads --outdir $outputfolder
	    rm  $rawFastqDir/*.fastq
	    ;;
	*.gz)
	    echo "this file should be gunzipped first"
	    gunzip -c $rawFastqDir$file > $rawFastqDir/$(basename "${file}" .gz)
	    fastqc $rawFastqDir/*.fastq --threads $fastQCthreads --outdir $outputfolder
	    rm  $rawFastqDir/*.fastq
	    ;;
	*.fastq)
	    echo "finally a good file extension"
	    fastqc $rawFastqDir/$file --threads $fastQCthreads --outdir $outputfolder
	    gzip $rawFastqDir/$file
	   
	    ;;
        *)
	    echo "there is another problem"
	    ;;
    esac
done

rm $rawFastqDir/*.fastq

exit
