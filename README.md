# _tPipe version 2_ 


This software was especially designed and constructed to analyse small RNA-Seq data and detect tRNA and tsRNA molecules. There are 6 main steps of preprocessing and analysing the data that include:
	* Quality Control
    * Preprocessing
    * Aligning reads against the reference tLibraries
    * Classification
    * Expression Matrix
    * Basic Statistics

<br />

## tPipe In Details

####1. Quality Control
FastQC will perform a series of preliminary quality assessments on the raw data and provide detailed quality control reports for each imputed sample. Such analysis may highlight issues and artefacts in the raw data. 

####2. Preprocessing
Preprocessing steps are critical for handling technical induced errors. Preprocessing includes Low Quality Trimming, Adapter Trimming,Ambiguous bases removal and Minimum length filtering.

####3. Aligning reads against the reference tLibraries
Mapping against the reference tLibraries will help us identifying the reads that belong to the tRNA and tsRNA class. Aligning includes a filtering step that will check all reads that include 1 and 2 mismatches. 

####4. Classification
The detected molecules will be classified and summarised. This step is entirely based on the secondary structure of the detected molecules and their origin (secondary structure of tgenes)

####5. Expression Matrix
The classified molecules will be summarised in several expression matrices. (different versions)

####6. Basic Statistics
Several secondary files will be generated that may help the user to continue with further analysing the data. This step includes a statistical analysis and quantification of the detected mismatches along with the frequency of occurrence.

<br />

## __Installation__ and __Info__

The software is designed for Linux (x64) operating system. 

To install:

Enter repository ("scr" folder) and run:
	
    $ python setup.py install

<br />

### __Help__

For further information about each compartment of the pipeline you can run:

	$ tpipe --help

<br />

### __Dependencies__

tPipe is dependent on xlsxwriter, xlrd, natsort, multiqc, cutadapt, bowtie1, fastx_toolkit and fastqc

<br />

### __Important Notes__
The pipeline needs 2 folders to run.
Folder1 - "RefLibrary" that hosts the indexes of the reference genome (_hg19_) - important for the filtering step
Folder2 - "tLibraries" that hosts the indexes of the reference tLibraries (_p_ and _m_) along with the files _reflibs_mu.fa_ and _sslib_mu.fa_

Folder1 (path) _along with the prefix of the reference genome_ should be mentioned on the configuration file (_Ref_Genome_)
Folder2 (only path) should be mentioned on the configuration file (_tLibs_)

For further details check configuration file or contact the author.