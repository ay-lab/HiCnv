$stamp = `date|awk '{print \$3_\$2_\$6}'`;
$job   = $$;
chomp $stamp;

################ User specific parameters ###############
$prefix = "WM2664"; #For replicate 1 and change the name while doing for replicate 2

#$hic_bwt2_folder = "path_to/hicresult/bowtie_results/bwt2/data"; #HiC-Pro output file path of replicate 1

#$hic_bwt2_folder_FWD = `ls $hic_bwt2_folder/*1_*.bwt2merged.bam`;

$hic_bwt2_folder_FWD = "/mnt/BioAdHoc/Groups/vd-ay/abhijit/overflow/proj_overflow/KadirMelanomaHiC_2017/WM2664_102617/hicresult_rep3/bowtie_results/bwt2/data/WM2664_Rep3_R1_hg19.bwt2merged.bam";

#$hic_bwt2_folder_REV = `ls $hic_bwt2_folder/*2_*.bwt2merged.bam`;

$hic_bwt2_folder_REV = "/mnt/BioAdHoc/Groups/vd-ay/abhijit/overflow/proj_overflow/KadirMelanomaHiC_2017/WM2664_102617/hicresult_rep3/bowtie_results/bwt2/data/WM2664_Rep3_R2_hg19.bwt2merged.bam";

$fragment_file = "HindIII_resfrag_hg19.bed"; # Change this as per your experiment

$F_GC_MAP_file = "HindIII.hg19.fragments.F_GC_MAP.bed"; #Creat this file as per your restriction fragment. Scripts to create this file are under ../scripts/F_GC_MAP_Files/ directory.
##########################################################

chomp ($hic_bwt2_folder_FWD,$hic_bwt2_folder_REV);
open (out, ">1DReadCoverage.$stamp.$job.sh");
print out "python ../scripts/mergeSAM-singletons.py -f $hic_bwt2_folder_FWD -r $hic_bwt2_folder_REV -o $prefix.bwt2pairs.withSingles.mapq30.bam -v --single -q 30\n";
print out "mkdir $prefix\_allMap2FragmentsOutput\n";
print out "python ../scripts/mapped_2hic_fragments.py -f ../scripts/$fragment_file -s 100 -l 800 -d 1000 -r $prefix.bwt2pairs.withSingles.mapq30.bam -o $prefix\_allMap2FragmentsOutput --all -v\n";
print out "cat $prefix\_allMap2FragmentsOutput/$prefix.bwt2pairs.withSingles.mapq30.perREfragStats | sort -k1,1 -k2,2n > $prefix.perREfragStats\n";
print out "perl ../scripts/mapCount.pl $prefix ../scripts/$F_GC_MAP_file";
close out;
`chmod 755 1DReadCoverage.$stamp.$job.sh`;
