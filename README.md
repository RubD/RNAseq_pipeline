-   Description
    -   Pipeline usage
    -   Pipeline commands
-   Pipeline usage
    -   1. Fill in Excell file
    -   2. Convert Excell file to Yaml file
    -   3. pipeline part I : quality analysis
    -   4. pipeline part II : aligment and counting features
-   Pipeline commands
    -   Quality control of the reads and statistics
        -   1. FastQC tool
        -   2. FastX tool
        -   3. Organism length and coverage
    -   Mapping the reads with Bowtie2 and TopHat
        -   1. Prepare index file (if needed)
        -   2. Download annotated reference genome from iGenomes
        -   3. Mapping with TopHat
    -   Counting mapped reads to features (genes, exons, ...)
        -   1. Counting with featureCounts
        -   2. Counting with HTSeq
    -   Statistical analysis
        -   Differentially expressed genes
        -   Gene Ontology analysis

<br><br>

**Installations**:

-   [R](http://www.r-project.org/) with packages [xlsx](http://cran.r-project.org/web/packages/xlsx/index.html) and [yaml](http://cran.r-project.org/web/packages/yaml/index.html)

-   [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) as standalone tool or from <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.3_devel.zip> (reason see: <http://seqanswers.com/forums/showthread.php?t=4846&page=16>) to use in terminal pipeline

-   [FastX](http://hannonlab.cshl.edu/fastx_toolkit/) to trim bad bases

-   parallel: on MAC OS use `$ brew install parallel`

-   bowtie2

-   [TopHat](http://ccb.jhu.edu/software/tophat/manual.shtml)

-   [subread](http://subread.sourceforge.net/) for [featureCounts](http://bioinf.wehi.edu.au/featureCounts/). Also as R package : [Rsubread](http://bioconductor.org/packages/release/bioc/html/Rsubread.html)

-   [HTSeq](http://www-huber.embl.de/users/anders/HTSeq/doc/install.html#install)

Add to your path if needed

<br><br>

Description
===========

### Pipeline usage

An easy to use general pipeline for RNAseq analysis that can quickly be adapted and run semi-independently. A user first has to specify all the parameters in an excell file for the scripts to run. The excell file also contains some general explanation of the different parameters. Next, this excell file can be converted to a yaml configuration (config) file via an R script. This yaml config file can subsequently be used for the quality analysis of the raw reads. TO BE DONE: if necessary trim raw fastq files.

The same yaml config file can then be used for Tophat2 aligment and creation of feature count tables based on featureCount or HTSeq.

TO BE DONE: the resulting count tables can be imported for Differential expression analysis with DESeq2 package in R.

### Pipeline commands

Explanation of the individual commands / tools that can be - or already are - incorporated in the pipeline.

-   Quality control: FastQC and FastX
-   Mapping and alignment of reads: Bowtie(2) and Tophat2
-   Counting mapped reads to features

* * * * *

<br><br>

Pipeline usage
==============

Put script in your PATH `$ echo $PATH` and make them executable `$ chmod +x script`

1. Fill in Excell file
----------------------

An excell file, that can be used as an example is given. Please do not delete or add additional columns or rows. Only the values of the second column (header = values) can be changed. More information about possible parameters can be found in the 3th column (header = description).

2. Convert Excell file to Yaml file
-----------------------------------

    # creates RNAseqConfig.yaml in working directory
    $ makeYamlConfigFromXlsx.r /path-to-excell-file/preconfigfile.xls

    # creates RNAseqConfig.yaml in designated directory
    $ makeYamlConfigFromXlsx.r /path-to-excell-file/preconfigfile.xls /path-for-outputfile/

3. pipeline part I : quality analysis
-------------------------------------

    $ fastQC_pipeline.r /path-to-dir/RNAseqConfig.yaml

4. pipeline part II : aligment and counting features
----------------------------------------------------

    $ tophat2_count_pipeline.r /path-to-dir/RNAseqConfig.yaml

Pipeline commands
=================

Quality control of the reads and statistics
-------------------------------------------

### 1. FastQC tool

Analyze all the .fastq files with the [FASTQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) standalone tool

In terminal as part of pipeline

    $ cd /folder/with/files.fastq       
    $ mkdir OutputFastQC                                    # create folder to store the fastQC results
    $ fastqc -h                                             # help page with all the available parameters
    $ fastqc *.fastq --threads 4 --outdir ./OutputFastQC    # * for all fastq files, --threads = how many files to analyze at the same time

### 2. FastX tool

If FastQC analysis shows that bases at the end have a systematic low score, we can trim them

    # example for just the first lines and removal of last 5 bases
    $ fastx_timmer -h                                                       # help page
    $ head file.fastq | fastx_trimmer -t 5 -Q33                             # Q33 = use of the sanger/Phred+33 score 
    $ fastx_trimmer -Q33 -t 5 -i file.fastq -o trimmed_file.fastq           # input and output files
    $ find *.fastq | parallel "fastx_trimmer -Q33 -t 5 -i {} -o trimmed_{}" # use parallel to do all files at once (for loop would be alternative)

### 3. Organism length and coverage

Go to [NCBI Genome website](http://www.ncbi.nlm.nih.gov/genome) and identify your genome / organism of interest

Your number of reads for each sample can be found in the FastQC quality report

Rule of thumb: for human genome (3 Gb) aim for at least 10 million reads

<br><br>

* * * * *

Mapping the reads with Bowtie2 and TopHat
-----------------------------------------

### 1. Prepare index file (if needed)

-   open terminal and go to data directory
-   try out bowtie2

        $ bowtie2

-   The reference genome must be in specific format (index) for bowtie to be able to use it.
-   Download genome (fasta file) from [NCBI Genome website](http://www.ncbi.nlm.nih.gov/genome), [Ensembl website](http://www.ensembl.org/info/data/ftp/index.html?redirect=no) or [UCSC website](http://hgdownload.soe.ucsc.edu/downloads.html)
-   Build the index for bowtie

        $ bowtie2-build index_organism_genome.fasta index_organism_genome

### 2. Download annotated reference genome from iGenomes

Several pre-built indexes are available to download on the [bowtie2 webpage](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).

iGenomes provides ready-to-use reference sequences and annotations from Ensembl, NCBI and UCSC. Example for mm9.

-   download annotations and reference genome to current folder

<!-- -->

    $ wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/UCSC/mm9/Mus_musculus_UCSC_mm9.tar.gz 

-   use 'nohup ... &' , executes commands after you exit from the terminal and/or 'nice', which determines the priority of the command -20 to 19, respctl highest to lowest

<!-- -->

    $ nohup nice -n 19 wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/UCSC/mm9/Mus_musculus_UCSC_mm9.tar.gz &

-   example to download mm10 to destination folder; takes very long time. Only necessary once on local or server

<!-- -->

    $ wget -P /Users/ruben/Documents/genomedata/ ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Mus_musculus/UCSC/mm10/Mus_musculus_UCSC_mm10.tar.gz 
    $ tar -zxvf /Users/ruben/Documents/genomedata/Mus_musculus_UCSC_mm10.tar.gz

### 3. Mapping with TopHat

TopHat is a splice junction mapper and can map across exon/exon boundaries.

Best to perform mapping on server with adequate power.

    $tophat2 -p 20 --no-coverage-search --segment-mismatches 1 --segment-length 18 --GTF /mnt/nfs/data/ruben/iGenomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf --library-type fr-unstranded --output-dir s_4_1_I12_1093_01_KO13_b1_d0_tophat2 /mnt/nfs/data/ruben/iGenomes/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome s_4_1_I12_1093_01_KO13_b1_d0.fastq.gz

**Parameters for TopHat(2)**

-   Can also be combined with the use of `nohup` and/or `nice`

-   tophat2 can recognize .fastq.gz files, so no need to convert them first

-   p or --num-threads : number of threads allocated to this command, the more the faster

-   no-coverage search : coverage search makes process really slow, can be switched off if you don't want to discover new splicing or new genes (see [biostars](https://www.biostars.org/p/49224/))

-   --segment-mismatches : defaults to 2, set to 0 or 1 for short reads (36bp) (see [seqanswers](http://seqanswers.com/forums/showthread.php?t=12195))

-   --segment-length : defaults to 25, set to half of read length for short reads (36/2 = 18) (see [seqanswers](http://seqanswers.com/forums/showthread.php?t=12195))

-   --GTF : gene annotation file

-   --library-type : depends on the library preparation kit (see tophat [FAQ](http://ccb.jhu.edu/software/tophat/faq.shtml#library_type))

-   --output-dir : defaults to tophat\_out/ , better and easier to make output folders for file.fastq

-   build genome data

-   fastq or fast.gz file

For Agata's RNAseq data:

    $ cd /folder_with_fastq.gz_files/
    $ nohup find *fastq.gz | parallel -j 10 "tophat2 -p 20 --no-coverage-search --segment-mismatches 1 --segment-length 18 --GTF /mnt/nfs/data/ruben/iGenomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf --library-type fr-unstranded --output-dir {}_tophat2 /mnt/nfs/data/ruben/iGenomes/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome {}" &

Create Index, required for some downstream analyses (like IGV):

    $ find */accepted_hits.bam | parallel "samtools index {}"

<br><br>

* * * * *

Counting mapped reads to features (genes, exons, ...)
-----------------------------------------------------

Several tools/algorithms/software can be used to count the mapped features:

[Here](http://genomespot.blogspot.be/2014/09/read-counting-with-featurecounts.html) you can find a speed/memory comparison.

### 1. Counting with featureCounts

On server:

    $ featureCounts -a /mnt/nfs/data/ruben/iGenomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf -o featC_counts.txt -t exon -g gene_name -T 10 */accepted_hits.bam

### 2. Counting with HTSeq

Most default settings are ok.

On server:

    $ samtools view accepted_hits.bam | htseq-count - /mnt/nfs/data/ruben/iGenomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf > sample_accepted_hits.htcounts

    # for all accepted_hits.bam, make sure you are in the right folder
    $ nohup find */accepted_hits.bam | parallel -j 10 "samtools view {} | htseq-count - /mnt/nfs/data/ruben/iGenomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf > {}.HTseqcounts" &

Or one by one with for loop: TO BE TESTED:

    $ cd /folder_with_tophatOutput_folders/

    # create script
    $ nano HTseqForloop.sh

    --------------------------
    #!/usr/bin/bash

    date
    for i in `find */accepted_hits.bam`;
    do
      bam=$i
      count="${bam/.bam}/_HTseqcount"
      echo -e "Counting"$bam
      samtools view $bam | htseq-count - /mnt/nfs/data/ruben/iGenomes/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf > $count
    done
    date
    --------------------------

    $ chmod +x HTseqForloop.sh
    $ ./HTseqForloop.sh

<br><br>

* * * * *

Statistical analysis
--------------------

### Differentially expressed genes

-   DESeq2

-   EdgeR

### Gene Ontology analysis

-   topGO

-   Gage & Pathview

<br><br>

* * * * *

**RESOURCES USED**

-   complete and well explained tutorial: <http://bioconnector.github.io/workshops/lessons/rnaseq-1day/>

-   <http://homer.salk.edu/homer/basicTutorial/index.html>

-   iGenomes reference genomes: <http://support.illumina.com/sequencing/sequencing_software/igenome.html>

-   RNAseq with iGenomes: <http://www.illumina.com/documents/products/technotes/RNASeqAnalysisTopHat.pdf>

-   <https://biobeat.wordpress.com/2013/05/16/bowtie2-htseq-count-reference-downloaded-from-ncbi/>

-   complete pipeline in R: <http://www.nature.com/nprot/journal/v8/n9/full/nprot.2013.099.html>

-   super blogspot: <http://genomespot.blogspot.be/>

-   interesting blogspot: <http://crazyhottommy.blogspot.be/>

-   complete tutorial: <http://web.science.uu.nl/pmi/publications/PDF/2013/TiPS-VanVerk-Hickman-2013-Hands_on_Tutorial.pdf>
