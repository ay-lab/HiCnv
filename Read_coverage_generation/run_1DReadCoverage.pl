$stamp = `date|awk '{print \$3_\$2_\$6}'`;
chomp $stamp;

sub GenerateMasterJobsFiles {

  ################ User specific parameters ###############
  my $prefix = $_[0]; #For replicate 1 and change the name while doing for replicate 2
  my $replicate_prefix = $_[1];
  my $hicpro_runpath   = $_[2];
  my $hic_bwt2_folder_FWD = "$hicpro_runpath/bowtie_results/bwt2/data/$replicate_prefix\_1_*.bwt2merged.bam";
  my $hic_bwt2_folder_REV = "$hicpro_runpath/bowtie_results/bwt2/data/$replicate_prefix\_2_*.bwt2merged.bam";
  my $fragment_file = "HindIII_resfrag_hg38.bed"; # Change this as per your experiment
  my $F_GC_MAP_file = "HindIII_hg38.500.50.F_GC_MAP.bed"; #Creat this file as per your restriction fragment. Scripts to create this file are under ../scripts/F_GC_MAP_Files/ directory.
  ##########################################################


  $pwd = `pwd`;
  chomp $pwd;

my @job_param ="#!/bin/bash -ex
#PBS -l nodes=1:ppn=2
#PBS -l mem=50GB
#PBS -l walltime=48:00:00
#PBS -m ae
#PBS -j eo
#PBS -V
#hostname
#TMPDIR=/scratch
cd $pwd\n";


  chomp ($hic_bwt2_folder_FWD,$hic_bwt2_folder_REV);
  open (out, ">$prefix.1DReadCoverage.$stamp.sh");
  print out @job_param;
  print out "python ../../scripts/mergeSAM-singletons.py -f $hic_bwt2_folder_FWD -r $hic_bwt2_folder_REV -o $prefix.bwt2pairs.withSingles.mapq30.bam -v --single -q 30\n";
  print out "mkdir $prefix\_allMap2FragmentsOutput\n";
  print out "python ../../scripts/mapped_2hic_fragments.py -f ../../scripts/$fragment_file -s 100 -l 800 -d 1000 -r $prefix.bwt2pairs.withSingles.mapq30.bam -o $prefix\_allMap2FragmentsOutput --all -v\n";
  print out "cat $prefix\_allMap2FragmentsOutput/$prefix.bwt2pairs.withSingles.mapq30.perREfragStats | sort -k1,1 -k2,2n > $prefix.perREfragStats\n";
  close out;
  `chmod 755 $prefix.1DReadCoverage.$stamp.sh`;
  `qsub $prefix.1DReadCoverage.$stamp.sh`;
}

my @Files = qw(A549);
my $HiCpro_run = "/mnt/BioAdHoc/Groups/vd-ay/abhijit/overflow/proj_overflow/HiCProcessingPipeLine/Cancer_Cell_Lines/HiC-Pro_Run"; 

foreach my $f (@Files) {
  chomp $f;
  if (! -d "$f\_Read_Coverage") {
    system("mkdir $f\_Read_Coverage");
  }
  my @repFiles = `ls $HiCpro_run/$f\_data_HiCPro/bowtie_results/bwt2/data/*_1_*.bwt2merged.bam |sed -e 's/\\// /g' |awk '{print \$NF}' |cut -d_ -f 1`;
  chdir "$f\_Read_Coverage";
  my $i = 0;
  while ($i <= $#repFiles) {
    chomp $repFiles[$i]; 
    print "$f\_$i\t$repFiles[$i]\n";
    GenerateMasterJobsFiles("$f\_Rep$i", "$repFiles[$i]", "$HiCpro_run/$f\_data_HiCPro");
    $i++;
  }
  chdir "../";
}
undef @Files;
