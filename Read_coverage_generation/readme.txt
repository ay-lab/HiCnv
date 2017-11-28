"run_1DReadCoverage.pl" script will create the *.perREfragStats and *.F_GC_MAP.bed files from Hi-C data.
This script will requre a *.F_GC_MAP.bed file, which is restriction fragment and genome specific file. To create this file go to "../scripts/F_GC_MAP_Files/" directory and type "perl create_F_GC_MAP_file.pl". Change the required variables of the script as per your requirement.

Note: 
If you have multiple replicates and generated multiple *.F_GC_MAP.bed and *.perREfragStats files, then please run the "perl combine_multiple_replicates.pl" program.
This script will combine all the *.F_GC_MAP.bed and *.perREfragStats files and generate a single "Sample.combined.F_GC_MAP.bed" and "Sample.combined.perREfragStats".
Use these combined files to run the next scripts under CNV_calling folder.
