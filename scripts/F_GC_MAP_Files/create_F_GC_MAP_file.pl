#Download the file wget ftp://hgdownload.cse.ucsc.edu/gbdb/hg19/bbi/wgEncodeCrgMapabilityAlign50mer.bw
#bigWigToBedGraph wgEncodeCrgMapabilityAlign50mer.bw hg19.MapabilityAlign50mer.bedGraph
#
############ User specific parameters ############
$frag = "HindIII_resfrag_hg19.bed"; #Change this file as per your experiment. Example HindIII restriction fragment files are provided in the "HindIII_resfrag_files.zip" folder.

$read_length = 50; #This is the read length of your fastq file.

$frag_file_name = "HindIII_hg19.$frag_size.$read_length";

$fasta_file = "hg19.fa"; #Needs hg19.fa or equivalent genome sequence file

$map_bedgraph = "hg19.MapabilityAlign50mer.bedGraph"; #Need to download the mappability file separately similar to fastq read length
#####################################################
$frag_size = 500;

$frag_half = $frag_size/2;

open (out,">F_GC_MAP.file.sh");
print out "awk '{print \$1\"\\t\"\$2\"\\t\"\$2+$frag_half\"\\t\"\$4\"\\t\"\$2\"\\t\"\$3\"\\t\"\$3-\$2\"\\n\"\$1\"\\t\"\$3-$frag_half\"\\t\"\$3\"\\t\"\$4\"\\t\"\$2\"\\t\"\$3\"\\t\"\$3-\$2}' $frag|awk '{if(\$2 >= 0){print}}'|sortBed > $frag_file_name.bed\n";
print out "bedtools nuc -fi $fasta_file -bed $frag_file_name.bed > $frag_file_name.GC.bed\n";
print out "bedtools map -a $frag_file_name.GC.bed -b ../$map_bedgraph -c 4 -o mean > $frag_file_name.GC_Map.bed\n";
print out "perl gc_map_per_fragment.pl $frag_file_name.GC_Map.bed $frag > ../$frag_file_name.F_GC_MAP.bed";
close out;
`chomp 755 F_GC_MAP.file.sh`;
