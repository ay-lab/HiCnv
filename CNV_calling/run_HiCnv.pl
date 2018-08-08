$stamp = `date|awk '{print \$3_\$2_\$6}'`; 
$job   = $$;   
############ User specific parameters ###########
$sample_name = "T47D";	#Sample file name. All the output will aslo have this prefix

$gc_limit    = 0.2;	#20% or above GC content

$map_limit   = 0.5;	#50% or above mappability

$frag_limit  = 1000;	#1000bp or greater than fragment length for 6bp cutter 
#Change to 150 for 4bp cutter
#$frag_limit = 150;

$ref_chrom   = "chr2";	#Reference chromosome to calculate PIC and labelling the CNVs. 
#This is a critical choise, please refer PMID 29048467, for reference choice.

@chr_list    = `cat ../scripts/chr.list`;
################################################
chomp $stamp;

open (out,">HiCnv_$stamp.$job.sh");
print out "#HiCnv script to call CNVs from Hi-C data\n";
print out "#Job performed on $stamp; Job id: $job\n\n";
print out "#Normalizing the raw counts\n";
print out "Rscript ../scripts/normalize.r $sample_name $gc_limit $map_limit $frag_limit\n\n";
print out "#PIC calculation and copy number determination based on reference chromosome\n";
print out "Rscript ../scripts/pic.r $sample_name 0 $ref_chrom $gc_limit $map_limit $frag_limit\n\n";

print out "#KDE estimation and HMM segmentation\n";
$i = 0;
while ($i <= $#chr_list){
	chomp $chr_list[$i];
	$k = $i+1;	
	print out "echo \"chr\tstart\tend\tgc\tmap\tfrag\tcount\tnorm.count\" > $sample_name.$chr_list[$i].F_GC_MAP.NormCount.$gc_limit\_$map_limit\_$frag_limit.Index.bed\n";
	print out "grep -v \"start\" $sample_name.F_GC_MAP.NormCount.$gc_limit\_$map_limit\_$frag_limit.bed|grep \"$chr_list[$i]\\b\"|awk '{if(\$4 >= $gc_limit && \$5 >= $map_limit && \$6 >= $frag_limit){print}}' |awk '{c++;print \$1\"\\t\"c\"\\t\"(c+1)\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7\"\\t\"\$8}' >> $sample_name.$chr_list[$i].F_GC_MAP.NormCount.$gc_limit\_$map_limit\_$frag_limit.Index.bed\n";
	print out "grep -v \"start\" $sample_name.F_GC_MAP.NormCount.$gc_limit\_$map_limit\_$frag_limit.bed|grep \"$chr_list[$i]\\b\"|awk '{if(\$4 >= $gc_limit && \$5 >= $map_limit && \$6 >= $frag_limit){print}}' |awk '{c++;print \$1\"\\t\"c\"\\t\"(c+1)\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7\"\\t\"\$8\"\\t\"\$2\"\\t\"\$3}' > $sample_name.$chr_list[$i].Index.file\n";
	print out "Rscript --max-ppsize=500000 ../scripts/kde_hmm_segmentation.r $gc_limit $map_limit $frag_limit $sample_name 0 $k\n";
	$i++;
}
print out "\n";
print out "#Figure generation and CNV lalbelling at 45 and 40% threshold of the reference chromosome mean\n";
print out "Rscript ../scripts/createFigure.r $sample_name\n";
print out "Rscript ../scripts/cnv_labelling.r $sample_name 0.45\n";
print out "Rscript ../scripts/cnv_labelling.r $sample_name 0.40\n";
print out "perl ../scripts/bp_calc.pl\n\n";
print out "Rscript ../scripts/genome_wide_copy_number.r $sample_name.chrWise.txt\n\n";
print out "#Scanning for double minutes and homogeneously staining regions\n";
print out "Rscript ../scripts/dm_hsr.r $sample_name ../scripts/HindIII.hg19.fragments.F_GC_MAP.bed $sample_name $gc_limit $map_limit $frag_limit\n\n";
print out "#Organizing result\n";
print out "mkdir $sample_name\_$stamp\_$job\_result\n";
print out "mkdir $sample_name\_$stamp\_$job\_result/raw_files\n";
print out "mkdir $sample_name\_$stamp\_$job\_result/figures\n";
print out "mkdir $sample_name\_$stamp\_$job\_result/breakpoints\n";
print out "mv $sample_name.*.optimal.* $sample_name\_$stamp\_$job\_result/raw_files/\n";
print out "mv $sample_name.*.count.bed $sample_name\_$stamp\_$job\_result/raw_files/\n";
print out "mv $sample_name.*.Index.* $sample_name\_$stamp\_$job\_result/raw_files/\n";
print out "mv $sample_name*.F_GC_MAP.bed $sample_name\_$stamp\_$job\_result/raw_files/\n";
print out "mv $sample_name.*.NormCount.* $sample_name\_$stamp\_$job\_result/raw_files/\n";
print out "mv $sample_name.chrWise.txt $sample_name\_$stamp\_$job\_result/raw_files/\n";
print out "mv $sample_name.merged_*.WithOutFilt.bed $sample_name\_$stamp\_$job\_result/raw_files/\n";
print out "mv $sample_name.copyNumber.txt $sample_name\_$stamp\_$job\_result/raw_files/\n";
print out "mv $sample_name.$gc_limit\_$map_limit\_$frag_limit.pdf $sample_name\_$stamp\_$job\_result/raw_files/\n";
print out "mv $sample_name.*.png $sample_name\_$stamp\_$job\_result/figures\n";
print out "mv $sample_name.*.WithOutFilt.bed.bp_assigned.txt $sample_name\_$stamp\_$job\_result/breakpoints\n";
print out "mv *.copy_number.pdf $sample_name\_$stamp\_$job\_result/figures/\n";
print out "mv *_copynumber.txt $sample_name\_$stamp\_$job\_result/breakpoints/\n";
print out "mv dm.$sample_name.pval.txt $sample_name\_$stamp\_$job\_result/breakpoints\n";
print out "mv hsr.$sample_name.pval.txt $sample_name\_$stamp\_$job\_result/breakpoints\n";
close out;
`chmod 755 HiCnv_$stamp.$job.sh`;
