Requirement:

A. R/perl environment

B. R packages required

	1. ineq : install.packages("https://cran.r-project.org/src/contrib/ineq_0.2-13.tar.gz")
	2. KernSmooth : install.packages("https://cran.r-project.org/src/contrib/KernSmooth_2.23-15.tar.gz")
	3. MASS : install.packages("https://cran.r-project.org/src/contrib/MASS_7.3-47.tar.gz")
	4. RHmm : install.packages("https://cran.r-project.org/src/contrib/Archive/RHmm/RHmm_2.0.3.tar.gz")
	5. GenVisR : https://bioconductor.org/packages/release/bioc/html/GenVisR.html

First copy the *.perREfragStats and *.F_GC_MAP.bed file from "Read_coverage_generation" folder. Type "perl run_HiCnv.pl", and this will create a bash script file with all the required commands to call CNV from Hi-C data. Please change the run_HiCnv.pl file variables as per your requirement. The example files T47D.R1.F_GC_MAP.bed and T47D.R2.F_GC_MAP.bed represent two replicates from T47D Hi-C experiment. These files were generated from "run_1DReadCoverage.pl" script.  

Note:
If you have multiple replicates and generated multiple *.F_GC_MAP.bed and *.perREfragStats files, then please run the "perl combine_multiple_replicates.pl" program in "Read_coverage_generation" folder.
The script will combine all the *.F_GC_MAP.bed and *.perREfragStats files and generate a single "Sample.combined.F_GC_MAP.bed" and "Sample.combined.perREfragStats".
Use these combined files to run the next scripts under CNV_calling folder. Change the "$sample_name" variable inside "run_HiCnv.pl" script to "Sample.combined" for downstream analysis.
