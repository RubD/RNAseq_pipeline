#!/bin/bash

# here the arguments are from within the do_Tophat2_countFeat2.sh script
args=("$@")

countFeatMethod=${args[0]}

genomeAnnotationDir=${args[1]}

resultsOutputDir=${args[2]}

FCfeatureType=${args[3]}

FCmetaType=${args[4]}

FCthreads=${args[5]}

HTSfeatureType=${args[6]}

HTSmetaType=${args[7]}

paralleljobs=${args[8]}

outputfolder=$resultsOutputDir/OutputTophat2

echo $(date)

if [ $countFeatMethod == featureCounts ]; then
    cd $outputfolder
    featureCounts -a $genomeAnnotationDir -o $resultsOutputDir/featC_counts.txt -t $FCfeatureType - g $FCmetaType -T $FCthreads ./*/accepted_hits.bam

elif [ $countFeatMethod == HTSeq ]; then
    cd $outputfolder
    find */accepted_hits.bam | parallel -j $paralleljobs "samtools view {} | htseq-count -t $HTSfeatureType -i $HTSmetaType - $genomeAnnotationDir > {}.HTseqcounts"

    
    # rename the output files of HTSeq to give it sample specific names
    for file in $(ls */*.HTseqcounts)
    do
	newname=$(dirname $file)
        newname=${newname/.fastq_tophat}
	
	mv $file $file$newname
    done

    # move the files to the output directory for HTSeq
    mkdir $resultsOutputDir/OutputHTSeq
    mv */*.HTseqcounts* $resultsOutputDir/OutputHTSeq

    # rename the files for proper naming    
    for file in $(ls $resultsOutputDir/OutputHTSeq)
    do
	cd $resultsOutputDir/OutputHTSeq  # WHY? 'cd' is required for 'mv' to work properly
	finalname=${file##*.}.txt
	mv "$file" "$finalname"
    done
    
else
    echo "no appropriate counting algorithm was chosen: featureCounts or HTSeq"
fi

echo $(date)
