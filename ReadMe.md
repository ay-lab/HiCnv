# HiCnv version 2

HiCnv is a pipeline to call CNVs from Hi-C data.
"scripts" folder contains all the scripts to call CNVs from Hi-C data.

# Steps: Process the HiC data and generate coverage, GC content, mappability and fragment length information file

1) Process your Hi-C fastq files with HiCPro pipline (https://github.com/nservant/HiC-Pro)


2) Processed HiCPro output directory will have a folder like "<path>/hicresult/bowtie_results/bwt2/data". Under this data folder there will be two *.bwt2merged.bam file (*_1_*.bwt2merged.bam and *_2_*.bwt2merged.bam). 


3) "Read_coverage_generation/run_1DReadCoverage.pl" require these two bam files for CNV finding. Change the "$hic_bwt2_folder" variable of "run_1DReadCoverage.pl" script to "<path>/hicresult/bowtie_results/bwt2/data" and type "perl run_1DReadCoverage.pl" inside "Read_coverage_generation" directory. This will create a "1DReadCoverage.*.sh" file.


4) "1DReadCoverage.*.sh" will require a restriction fragment specific *.fragments.F_GC_MAP.bed file(Fragment length, GC content and Mappability information file). To generate this file as per your experiment, please go to "scripts/F_GC_MAP_Files/" directory and type "./F_GC_MAP.file.sh". Change the variables of the script to match your requirement. 


5) Running "1DReadCoverage.*.sh" will create the *.perREfragStats and *.F_GC_MAP.bed files from Hi-C data for downstream processing.

6) Run hicnv_v2.R, check the usage for more details.

```bash
Usage: hicnv_v2.R [options]


Options:
        --refeature=REFEATURE
                A five column restriction enzyme cutsite bed file with GC content, mappability, fragment length infomation

                <chr>   <start> <end>   <GC frequency>  <mappability>   <fragment length>


        --coverage=COVERAGE
                A bedGraph file with read coverage signal (5th column).

                Alternatively, the .perREfragStats file


        --gccutoff=GCCUTOFF
                GC content cutoff. Anything below <gccutoff> will be removed [Default is 0.2].


        --mapcutoff=MAPCUTOFF
                Mappability cutoff. Anything below <mapcutoff> will be removed [Default is 0.5].


        --fragcutoff=FRAGCUTOFF
                Fragment length cutoff. Anything below <fragcutoff> will be removed [Default is 150].

                For Hi-C experiments with 4bp cut enzyme, this value is 150, for 6bp enzymes this value should be 1000.


        --refchrom=REFCHROM
                Name of the reference chromosome for CNV detection.

                If no name is provided then HiCnv proceed to estimate proportion of interaction count or PIC to decide a reference chromosome.


        --bandwidth=BANDWIDTH
                Genomic distance bandwidth with which Kernel smoothing will be performed [Default 1Mb].


        --hmmstate=HMMSTATE
                Number of HMM states to be searched in each chromosome [Default 10].


        --threshold=THRESHOLD
                Threshold value to define amplification and deletion with respect to normal region [Default is 0.2 i.e. deviation of 20% from mean normal value will be labeled as CNV].


        --prefix=PREFIX
                All the files and folders will be created with this name.


        --cpu=CPU
                Number of threads to be used [Default 1].


        -h, --help
                Show this help message and exit
```

# Note:

CNV calling requires GC content, mappability and fragment length information of every RE fragments. The *.F_GC_MAP.bed file contains all these information.
The file can be created using F_GC_MAP.file.sh script available under "scripts/F_GC_MAP_Files/" folder. To create the file, please run the script inside the "scripts/F_GC_MAP_Files/" folder. Also, change the variables as per the Hi-C experiment. 

For more details check "Rscript hicnv_v2.R --help".

# Covert an aligned Hi-C sam file into HiCnv usable format:

The script "samToHiCProFormat.pl" under "scripts/" folder takes an aligned HiC file in sam format. This can be a merged.sam (both forward and reverse reads merged/paired format) file or in single format file where forward and reverse reads are mapped into separate files (e.g forward.sam and reverse.sam).

The script an be run like the following 

perl samToHiCProFormat.pl -format paired -sam_file merged.sam -strand 0 -chr 3 -pos 4 -mapq 5 -read_len 6 -read_pos 10 -out_file test

or 

perl samToHiCProFormat.pl -format single -sam_file forward.sam,reverse.sam -strand 2 -chr 3 -pos 4 -mapq 5 -read_len 6 -read_pos 10 -out_file test

Check for the full details by "perl samToHiCProFormat.pl -help"

Script will create test.mate1.sam and test.mate2.sam file. Then add the header lines using samtools

samtools view -bT hg19.fa test.mate1.sam > mate1.bam

samtools view -bT hg19.fa test.mate2.sam > mate2.bam

Then in the run_1DReadCoverage.pl script of HiCnv, change the following variable to 

$hic_bwt2_folder_FWD = "mate1.bam";

$hic_bwt2_folder_REV = "mate2.bam";

This will enable HiCnv to read the bam files.

The "samToHiCProFormat_Example" folder contains example data to convert Hi-C sam files into HiCnv format.

Note: In forward.sam and reverse,sam files the total number of reads should be equal and assumed that each row of the two files should represent the same read.

# Double Minute (DM) and Homogeneously Staining Regions (HSR) scanning

Run ./scripts/dm_hsr.r to scan DM and HSR regions. 

For more details check "Rscript ./scripts/dm_hsr.r --help"

# Contact

abhijit@lji.org (Abhijit Chakraborty)
