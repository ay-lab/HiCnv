use Getopt::Long;

GetOptions(
        'format=s'      => \my $format,
        'sam_file=s'    => \my $sam_file,
	'strand=s'   	=> \my $strand,
	'chr=s'      	=> \my $chr,
	'pos=s'      	=> \my $pos,
	'mapq=s'	=> \my $mapq,
	'read_len=s'    => \my $read_len,
	'read_pos=s'    => \my $read_pos,
	'out_file=s'	=> \my $out,
	'help'		=> \my $help
);

## -format = either paired or single
## -sam_file = Name of the paired sam file. When format is single provide comma separated forward and reverse mapped read files e.g. forward.sam,reverse.sam
## -strand = sam file field that denotes the strand information for mate1 pair. If unknown make it 0 and everything will be one the positive strand mapping. 
## -chr = sam file field that denotes the chromosome information for pair.
## -pos = sam file field that denotes the chromosome position information for pair.
## -mapq = sam file field that denotes the quality score.
## -read_len = sam file field that denotes the read length.
## -read_pos = sam file field that denotes the read position.

## Paired sam file obtained from HiC-Pro output example:
## 19435780	67	chr14	72419306	42	50M	chr3	64908008	0	TGTAAAACCCCCGAGAACCTATTCAATTGTTATTTCTTTGATCTATATTG	BCBFFFFFGHHGGGJGIJJJJJJJJJJJJHIJJJJIJJJIGIJJJJIJJI	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:50	YT:Z:UU	RG:Z:BMG
## 19435780	131	chr3	64908008	42	50M	chr14	72419306	0	CTTTTAAAATAGACCCGGAAGCACAGTCAACTAGCCTTAACTCTTTTTCA	JJJJJJIHIHIJJJIJJIHHFFFFFEEDEEEDDDDDDDDDDDDDDDCDDD	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:50	YT:Z:UU	RG:Z:BMG
## 19435782	99	chrX	46444057	42	50M	=	46444136	0	CTGAAATCCAAATGCCTAGAGCGAGTGCTTGTTTGATTTTAATATCACTG	CCCDFFFFHGHHHJJJJJJIJJJJJFHHIIIHIJJIJJJJJJJJJJJIIJ	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:50	YT:Z:UU	RG:Z:BMG
## 19435782	147	chrX	46444136	42	50M	=	46444057	0	GGAAGCCAACCTAACTGAAGGGCTGAAGTGTGTAATTCTGCCCTGAAATA	CDDDDDDDEEEEFEEFFFFFHHHHHJIJJJJJJJJJJJJJIJJJJJJJJI	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:50	YT:Z:UU	RG:Z:BMG
##
## 19435780 and 19435782 are the read ids
## No strand id avaiable for HiC-Pro output like file. Strand ID is desirable. -strand option
## chr14, chr3, chrX and chrX represents the chromosome information field, i.e. -chr option
## 72419306, 64908008, 46444057, 46444136 shows the chromosome position. -pos option
## 42s' are phred quality score position, i.e. -mapq option
## 50Ms' are length of the aligned read, i.e. -read_len option
## "TGTAAAACCCCCGAGAACCTATTCAATTGTTATTTCTTTGATCTATATTG" and respective reads represent the -read_pos position. Here it is 10th column.
## 

chomp (
	$format, $sam_file, $strand, $chr,
	$pos, $mapq, $read_len, $read_pos,
	$out, $help
);

if ($help ne ""){
	print "
-format = either paired or single\n
-sam_file = Name of the paired sam file. When format is single provide comma separated forward and reverse mapped read files e.g. forward.sam,reverse.sam.\n
-strand = sam file field that denotes the strand information for mate1 pair. If unknown make it 0 and everything will be one the positive strand mapping.\n
-chr = sam file field that denotes the chromosome information for pair.\n
-pos = sam file field that denotes the chromosome position information for pair.\n
-mapq = sam file field that denotes the quality score.\n
-read_len = sam file field that denotes the read length.\n
-read_pos = sam file field that denotes the read position.\n

## Paired sam file obtained from HiC-Pro output example:
19435780     67      chr14   72419306        42      50M     chr3    64908008        0       TGTAAAACCCCCGAGAACCTATTCAATTGTTATTTCTTTGATCTATATTG      BCBFFFFFGHHGGGJGIJJJJJJJJJJJJHIJJJJIJJJIGIJJJJIJJI      AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:50 YT:Z:UU RG:Z:BMG
19435780     131     chr3    64908008        42      50M     chr14   72419306        0       CTTTTAAAATAGACCCGGAAGCACAGTCAACTAGCCTTAACTCTTTTTCA      JJJJJJIHIHIJJJIJJIHHFFFFFEEDEEEDDDDDDDDDDDDDDDCDDD      AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:50 YT:Z:UU RG:Z:BMG
19435782     99      chrX    46444057        42      50M     =       46444136        0       CTGAAATCCAAATGCCTAGAGCGAGTGCTTGTTTGATTTTAATATCACTG      CCCDFFFFHGHHHJJJJJJIJJJJJFHHIIIHIJJIJJJJJJJJJJJIIJ      AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:50 YT:Z:UU RG:Z:BMG
19435782     147     chrX    46444136        42      50M     =       46444057        0       GGAAGCCAACCTAACTGAAGGGCTGAAGTGTGTAATTCTGCCCTGAAATA      CDDDDDDDEEEEFEEFFFFFHHHHHJIJJJJJJJJJJJJJIJJJJJJJJI      AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:50 YT:Z:UU RG:Z:BMG

## Description ###
## 19435780 and 19435782 are the read ids
## No strand id avaiable for HiC-Pro output like file. Strand ID is desirable. -strand option
## chr14, chr3, chrX and chrX represents the chromosome information field, i.e. -chr option
## 72419306, 64908008, 46444057, 46444136 shows the chromosome position. -pos option
## 42s' are phred quality score position, i.e. -mapq option
## 50Ms' are length of the aligned read, i.e. -read_len option
## \"TGTAAAACCCCCGAGAACCTATTCAATTGTTATTTCTTTGATCTATATTG\" and respective reads represent the -read_pos position. Here it is 10th column.
## END ##\n
\n"
}

if ($format eq "paired"){
	
	open (mateA,">$out.mate1.sam");
	open (mateB,">$out.mate2.sam");
	open (sam,"$sam_file") or die "No such file $sam_file\n";
	while (<sam>){
		$c++;
		chomp $_;
		if ($c == 2){
			$id++;
			print "$id\n";
			$mateB = $_;
			@mateA = split(/\s+/,$mateA);
			@mateB = split(/\s+/,$mateB);
			
			## mate1
			if ($strand == 0){$strand_1 = 0;}
			else {$strand_1 = @mateA[$strand-1];}
			$chrom_1 = @mateA[$chr-1];
			$pos_1   = @mateA[$pos-1];
			$mapq_1  = @mateA[$mapq-1];
			$read_len_1 = @mateA[$read_len-1];
			$read_pos_1 = @mateA[$read_pos-1];
			if ($chrom_1 ne "*"){
				print mateA $id,"\t",$strand_1,"\t",$chrom_1,"\t",$pos_1,"\t",$mapq_1,"\t",$read_len_1,"\t*\t0\t0";
				$i = $read_pos-1;
				while ($i <= $#mateA){
					chomp $mateA[$i];
					print mateA "\t$mateA[$i]";
					$i++;
				}
				print mateA "\n";
			} else {
				print mateA $id,"\t4\t*\t0\t0\t*\t*\t0\t0";
				$i = $read_pos-1;
                                while ($i <= $#mateA){
                                        chomp $mateA[$i];
                                        print mateA "\t$mateA[$i]";
                                        $i++;
                                }
                                print mateA "\n";	
			}
			
			## mate2
			if ($strand == 0){$strand_2 = 0;}
                        else {$strand_2 = @mateB[$strand-1];}
                        $chrom_2 = @mateB[$chr-1];
                        $pos_2   = @mateB[$pos-1];
                        $mapq_2  = @mateB[$mapq-1];
                        $read_len_2 = @mateB[$read_len-1];
                        $read_pos_2 = @mateB[$read_pos-1];
                        if ($chrom_2 ne "*"){
                                print mateB $id,"\t",$strand_2,"\t",$chrom_2,"\t",$pos_2,"\t",$mapq_2,"\t",$read_len_2,"\t*\t0\t0";
                                $i = $read_pos-1;
                                while ($i <= $#mateB){
                                        chomp $mateB[$i];
                                        print mateB "\t$mateB[$i]";
                                        $i++;
                                }
                                print mateB "\n";
                        } else {
                                print mateB $id,"\t4\t*\t0\t0\t*\t*\t0\t0";
                                $i = $read_pos-1;
                                while ($i <= $#mateB){
                                        chomp $mateB[$i];
                                        print mateB "\t$mateB[$i]";
                                        $i++;
                                }
                                print mateB "\n";
                        }	
			undef @mateA;
			undef @mateB;
			$c = 0;
		}
		else {
			$mateA = $_;
		}
	}	
	close sam;
	close mateA;
	close mateB;
}

elsif ($format eq "single"){
	@samfile = split(/,/,$sam_file);
	open(f1, "@samfile[0]") or die "No such file @samfile[0]\n";
	open(f2, "@samfile[1]") or die "No such file @samfile[1]\n";
	open (mateA,">$out.mate1.sam");
        open (mateB,">$out.mate2.sam");

	while(($mateA = <f1>) and ($mateB = <f2>)){
		$id++;
	        chomp $mateA;
	        chomp $mateB;
		@mateA = split(/\s+/,$mateA);
                @mateB = split(/\s+/,$mateB);
		
		#Forward read processing
		if ($strand == 0){$strand_1 = 0;}
            	else {$strand_1 = @mateA[$strand-1];}
               	$chrom_1 = @mateA[$chr-1];
               	$pos_1   = @mateA[$pos-1];
            	$mapq_1  = @mateA[$mapq-1];
                $read_len_1 = @mateA[$read_len-1];
                $read_pos_1 = @mateA[$read_pos-1];
                if ($chrom_1 ne "*"){
                	print mateA $id,"\t",$strand_1,"\t",$chrom_1,"\t",$pos_1,"\t",$mapq_1,"\t",$read_len_1,"\t*\t0\t0";
                        $i = $read_pos-1;
                        while ($i <= $#mateA){
                        	chomp $mateA[$i];
                                print mateA "\t$mateA[$i]";
                                $i++;
           		}
                        print mateA "\n";
            	} else {
                	print mateA $id,"\t4\t*\t0\t0\t*\t*\t0\t0";
                        $i = $read_pos-1;
                        while ($i <= $#mateA){
                        	chomp $mateA[$i];
                                print mateA "\t$mateA[$i]";
                         	$i++;
                     	}
                        print mateA "\n";
         	}
		
		#Reverse read processing
		if ($strand == 0){$strand_2 = 0;}
                else {$strand_2 = @mateB[$strand-1];}
		$chrom_2 = @mateB[$chr-1];
                $pos_2   = @mateB[$pos-1];
                $mapq_2  = @mateB[$mapq-1];
                $read_len_2 = @mateB[$read_len-1];
                $read_pos_2 = @mateB[$read_pos-1];
                if ($chrom_2 ne "*"){
                	print mateB $id,"\t",$strand_2,"\t",$chrom_2,"\t",$pos_2,"\t",$mapq_2,"\t",$read_len_2,"\t*\t0\t0";
                        $i = $read_pos-1;
                        while ($i <= $#mateB){
                        	chomp $mateB[$i];
                                print mateB "\t$mateB[$i]";
                                $i++;
                      	}
                        print mateB "\n";
       		} else {
             		print mateB $id,"\t4\t*\t0\t0\t*\t*\t0\t0";
                        $i = $read_pos-1;
                        while ($i <= $#mateB){
                        	chomp $mateB[$i];
                                print mateB "\t$mateB[$i]";
                                $i++;
                    	}
                        print mateB "\n";
           	}

	        undef @mateA;
	        undef @mateB;
	}
	close f1;
	close f2;
	close mate1;
	close mate2;
	undef @samfile;
}
