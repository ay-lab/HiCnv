#!/bin/bash

## Check the https://www.pmgenomics.ca/hoffmanlab/proj/bismap/trackhub/ to download genome specific mappability files
wget https://www.pmgenomics.ca/hoffmanlab/proj/bismap/trackhub/hg38/k50.Umap.MultiTrackMappability.bw

## Process the bigWig file and create bedGraph file
bigWigToBedGraph k50.Umap.MultiTrackMappability.bw k50.Umap.MultiTrackMappability.bedGraph
sort -k 1,1 -k2,2n k50.Umap.MultiTrackMappability.bedGraph > k50.Umap.MultiTrackMappability.sorted.bedGraph

## Download the hg38 fasta files. Change the hg38 to hg19/mm9/mm10 as per the experiment. The chr.txt file contains the chromsome names
cat hg38_chr.txt |perl -ne 'chomp $_; system("curl http://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/$_.fa.gz|zcat > hg38.$_.fa");'
cat hg38.chr*.fa > hg38.fa

## Create extended 500bp restriction fragment file
awk '{print $1"\t"$2"\t"$2+250"\t"$4"\t"$2"\t"$3"\t"$3-$2"\n"$1"\t"$3-250"\t"$3"\t"$4"\t"$2"\t"$3"\t"$3-$2}' ../HindIII_resfrag_hg38.bed|awk '{if($2 >= 0){print}}'|sortBed > HindIII_hg38.500.50.bed

## Find GC percentage of 500bp regions
bedtools nuc -fi hg38.fa -bed HindIII_hg38.500.50.bed > HindIII_hg38.500.50.GC.bed

## Map the mappability over GC content file
bedtools map -a HindIII_hg38.500.50.GC.bed -b k50.Umap.MultiTrackMappability.sorted.bedGraph -c 4 -o mean > HindIII_hg38.500.50.GC_Map.bed

## Create F_GC_MAP file
perl gc_map_per_fragment.pl HindIII_hg38.500.50.GC_Map.bed ../HindIII_resfrag_hg38.bed > ../HindIII_hg38.500.50.F_GC_MAP.bed
