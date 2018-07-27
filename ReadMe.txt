HiCnv is a pipeline to call CNVs from Hi-C data.
"scripts" folder contains all the scripts to call CNVs from Hi-C data.
__________________________________________________________________

Step 1:

A. Process your Hi-C fastq files with HiCPro pipline (https://github.com/nservant/HiC-Pro)

B. Processed HiCPro output directory will have a folder like "<path>/hicresult/bowtie_results/bwt2/data". Under this data folder there will be two *.bwt2merged.bam file (*_1_*.bwt2merged.bam and *_2_*.bwt2merged.bam).

C. "Read_coverage_generation/run_1DReadCoverage.pl" require these two bam files for CNV finding. Change the "$hic_bwt2_folder" variable of "run_1DReadCoverage.pl" script to "<path>/hicresult/bowtie_results/bwt2/data" and type "perl run_1DReadCoverage.pl" inside "Read_coverage_generation" directory. This will create a "1DReadCoverage.*.sh" file.

D. "1DReadCoverage.*.sh" will require a restriction fragment specific *.fragments.F_GC_MAP.bed file(Fragment length, GC content and Mappability information file). To generate this file as per your experiment, please go to "scripts/F_GC_MAP_Files/" directory and type "perl create_F_GC_MAP_file.pl". Change the variables of the script to match your requirement. 

E. Running "1DReadCoverage.*.sh" will create the *.perREfragStats and *.F_GC_MAP.bed files from Hi-C data for downstream processing.
__________________________________________________________________

Step2:

A. Copy the *.perREfragStats and *.F_GC_MAP.bed files from "Read_coverage_generation" directory to "CNV_calling" directory.

B. Run "CNV_calling/run_HiCnv.pl" script to call CNVs from Hi-C data.

C. "T47D.example.input.zip" under "CNV_calling" has the example files.
__________________________________________________________________

Note:
CNV calling requires GC content, mappability and fragment length information of every RE fragments. The *.F_GC_MAP.bed file contains all these information.
The file can be created using create_F_GC_MAP_file.pl script available under "scripts/F_GC_MAP_Files/" folder. To create the file, please run the create_F_GC_MAP_file.pl file inside the "scripts/F_GC_MAP_Files/" folder. 

_____________________________________________________________________

Covert an aligned Hi-C sam file into HiCnv usable format:

The script "samToHiCProFormat.pl" under "scripts/" folder takes an aligned HiC file in sam format. This can be a merged.sam (both forward and reverse reads merged/paired format) file or in single format file where forward and reverse reads are mapped into separate files (e.g forward.sam and reverse.sam).

The script an be run like the following 

perl samToHiCProFormat.pl -format paired -sam_file merged.sam -strand 0 -chr 3 -pos 4 -mapq 5 -read_len 6 -read_pos 10 -out_file test

or 

perl samToHiCProFormat.pl -format single -sam_file forward.sam,reverse.sam -strand 2 -chr 3 -pos 4 -mapq 5 -read_len 6 -read_pos 10 -out_file test

perl samToHiCProFormat.pl -help

-format = either paired or single
-sam_file = Name of the paired sam file. When format is single provide comma separated forward and reverse mapped read files e.g. forward.sam,reverse.sam.
-strand = sam file field that denotes the strand information for mate pair. If unknown make it 0 and everything will be one the positive strand mapping.
-chr = sam file field that denotes the chromosome information for pair.
-pos = sam file field that denotes the chromosome position information for pair.
-mapq = sam file field that denotes the quality score.
-read_len = sam file field that denotes the read length.
-read_pos = sam file field that denotes the read position.

Script will create test.mate1.sam and test.mate2.sam file. Then add the header lines using samtools

samtools view -bT hg19.fa test.mate1.sam > mate1.bam
samtools view -bT hg19.fa test.mate2.sam > mate2.bam

Then in the run_1DReadCoverage.pl script of HiCnv, change the following variable to 

$hic_bwt2_folder_FWD = "mate1.bam";
$hic_bwt2_folder_REV = "mate2.bam";

This is enable HiCnv to read the bam files.

The "samToHiCProFormat_Example" folder contains example data to convert Hi-C sam files into HiCnv format.

Note: In forward.sam and reverse,sam files the total number of reads should be equal and assumed that each row of the two files should represent the same read.
_____________________________________________________________________

Contact

abhijit@lji.org (Abhijit Chakraborty)
