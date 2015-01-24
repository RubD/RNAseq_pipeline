#!/bin/bash
echo $(date)

args=("$@")

# path to the directory of the raw fastq files
rawFastqDir=${args[0]}

# path to the directory where the FastQC output directory must be created
resultsOutputDir=${args[1]} #tophat2outputDir=$2

# path to bowtie2index genome
bowtie2indexDir=${args[2]}

# which algorithm to use for counting
countFeatMethod=${args[3]} #countFeatAlg

# path to genome annotation
genomeAnnotationDir=${args[4]} #genomeAnnotation

# path to trimmed / quality fastq files
qualityFastqDir=${args[6]}

# number of threads assigned to tophat2 aligment
tophat2threads=${args[7]}

# number of alignment mismatch allowed
segmentMismatch=${args[8]}

# half the length of your reads
segmentLength=${args[9]}

# library type
libType=${args[10]}

# how many jobs in parallel
paralleljobs=${args[11]}


FCfeatureType=${args[12]}

FCmetaType=${args[13]}

FCthreads=${args[14]}

HTSfeatureType=${args[15]}

HTSmetaType=${args[16]}


# check whether tophat2 has already been performed
if [ -d $resultsOutputDir/OutputTophat2 ]; then
    
    if [ -d $resultsOutputDir/OutputHTSeq ] && [ $countFeatMethod == HTSeq ]; then
	echo "tophat2 and HTSeq output already exists, delete for new analysis"
	exit

    elif [ -e $resultsOutputDir/featC_counts.txt ] && [ $countFeatMethod == featureCounts ]; then
        echo "tophat2 and featureCounts output already exists, delete for new analysis"
	exit

    else
	do_countAlg2.sh $countFeatMethod $genomeAnnotationDir $resultsOutputDir $FCfeatureType $FCmetaType $FCthreads $HTSfeatureType $HTSmetaType $paralleljobs
	echo $(date)
	exit
    fi
else
    mkdir $resultsOutputDir/OutputTophat2
fi

# RUN Tophat2
# set outputfolder
outputfolder=$resultsOutputDir/OutputTophat2

cd $qualityFastqDir

for file in $(ls $qualityFastqDir | grep "fastq") #only selects fastq files
do
    case ${file} in
	*.tar.gz)
	    echo "this $file should be first be unpacked and gunzipped"
	    tar -zxvf $qualityFastqDir$file -C $qualityFastqDir
	    ;;
	*.gz)
	    echo "this file should be gunzipped first"
	    gunzip -c $qualityFastqDir$file > $qualityFastqDir/$(basename "${file}" .gz)
	    ;;
	*.fastq)
	    echo "finally a good file extension"
	    ;;
        *)
	    echo "there are other files present"
	    ;;
    esac
done


find *.fastq | parallel -j $paralleljobs "tophat2 -p $tophat2threads --no-coverage-search --segment-mismatches $segmentMismatch --segment-length $segmentLength --library-type $libType --output-dir $outputfolder/{}_tophat $bowtie2indexDir {}" 


cd $outputfolder
# create bam index, necessary for some downstream applications
find */accepted_hits.bam | parallel "samtools index {}"

echo " "
echo "here starts the counting"
echo " "

do_countAlg2.sh $countFeatMethod $genomeAnnotationDir $resultsOutputDir $FCfeatureType $FCmetaType $FCthreads $HTSfeatureType $HTSmetaType $paralleljobs
echo $(date)
