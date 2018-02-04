## Example to create perREfragStats and F_GC_MAP.bed files from Hi-C pro output. R1 is the forward read and R2 is the reverse read mapped file ## 
mkdir allMap2FragmentsOutput
python ../../scripts/mergeSAM-singletons.py -f MCF7.R1_hg19.bwt2merged.bam -r MCF7.R2_hg19.bwt2merged.bam -o MCF7.bwt2pairs.withSingles.mapq30.bam -v --single -q 30
python ../../scripts/mapped_2hic_fragments.py -f ../../scripts/HindIII_resfrag_hg19.bed -s 100 -l 800 -d 1000 -r MCF7.bwt2pairs.withSingles.mapq30.bam -o allMap2FragmentsOutput --all -v
cat allMap2FragmentsOutput/MCF7.bwt2pairs.withSingles.mapq30.perREfragStats | sort -k1,1 -k2,2n > MCF7.perREfragStats

## This script maps and combines the overall mapped reads with GC, Mappability and Fragment length information ##
perl ../../scripts/mapCount.pl MCF7 ../../scripts/HindIII.hg19.fragments.F_GC_MAP.bed
