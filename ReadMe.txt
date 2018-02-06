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
