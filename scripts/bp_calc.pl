@bp_file = `ls *.WithOutFilt.bed`;
$i = 0;
while ($i <= $#bp_file){
	chomp $bp_file[$i];
	open (out,">$bp_file[$i].bp_assigned.txt");
	@bpData = split(/\./,$bp_file[$i]);
	@file = `cat $bp_file[$i]`;
	$j = 0;
	while ($j <= $#file){
		chomp $file[$j];
		@fileData = split(/\s+/,$file[$j]);
		if ($file[$j] =~ /Amplification/){
			$total = (`cat @bpData[0].@fileData[0].Index.file|wc -l`)-20;
			$start = `awk '{if(\$9==@fileData[1] || \$10==@fileData[1]){print}}' @bpData[0].@fileData[0].Index.file|head -1|awk '{print \$2}'`;
			$end   = `awk '{if(\$9==@fileData[2] || \$10==@fileData[2]){print}}' @bpData[0].@fileData[0].Index.file|head -1|awk '{print \$2}'`;
			chomp ($start,$end,$total);
			if ($start < 20 && $end > $total){
				print out "$file[$j]\n";
			}
			elsif ($start < 20 && $end < $total){
				$start_mod = $start+20;
				@bp = `Rscript ../scripts/bp_assignment.r @bpData[0].@fileData[0].xy.optimal.max amp $start_mod $end|tail -1|awk '{print \$3}'`;
				chomp @bp[0];
				$bp_end = `awk '{if(\$2==@bp[0] || \$3==@bp[0]){print}}' @bpData[0].@fileData[0].Index.file|head -1|awk '{print \$9}'`;
				chomp $bp_end;
				print out "@fileData[0]\t@fileData[1]\t$bp_end\t@fileData[3]\t@fileData[4]\n";
				undef @bp;
			}
			elsif ($start > 20 && $end > $total){
				$end_mod = $end-20;
				@bp = `Rscript ../scripts/bp_assignment.r @bpData[0].@fileData[0].xy.optimal.max amp $start $end_mod|tail -2|head -1|awk '{print \$3}'`;
				chomp @bp[0];
				$bp_start = `awk '{if(\$2==@bp[0] || \$3==@bp[0]){print}}' @bpData[0].@fileData[0].Index.file|head -1|awk '{print \$9}'`;
				chomp $bp_start;
				print out "@fileData[0]\t$bp_start\t@fileData[2]\t@fileData[3]\t@fileData[4]\n";
				undef @bp;
			}
			elsif ($start > 20 && $end < $total){
				@bp = `Rscript ../scripts/bp_assignment.r @bpData[0].@fileData[0].xy.optimal.max amp $start $end|tail -2|awk '{print \$3}'`;
				chomp (@bp[0],@bp[1]);
				$bp_start = `awk '{if(\$2==@bp[0] || \$3==@bp[0]){print}}' @bpData[0].@fileData[0].Index.file|head -1|awk '{print \$9}'`;
				$bp_end   = `awk '{if(\$2==@bp[1] || \$3==@bp[1]){print}}' @bpData[0].@fileData[0].Index.file|head -1|awk '{print \$9}'`;
				chomp ($bp_start,$bp_end);
				print out "@fileData[0]\t$bp_start\t$bp_end\t@fileData[3]\t@fileData[4]\n";
				undef @bp;
			}
		}
		elsif ($file[$j] =~ /Deletion/){
                        $total = (`cat @bpData[0].@fileData[0].Index.file|wc -l`)-20;
                        $start = `awk '{if(\$9==@fileData[1] || \$10==@fileData[1]){print}}' @bpData[0].@fileData[0].Index.file|head -1|awk '{print \$2}'`;
                        $end   = `awk '{if(\$9==@fileData[2] || \$10==@fileData[2]){print}}' @bpData[0].@fileData[0].Index.file|head -1|awk '{print \$2}'`;
                        chomp ($start,$end,$total);
                        if ($start < 20 && $end > $total){
                                print out "$file[$j]\n";
                        }
                        elsif ($start < 20 && $end < $total){
                                $start_mod = $start+20;
                                @bp = `Rscript ../scripts/bp_assignment.r @bpData[0].@fileData[0].xy.optimal.max del $start_mod $end|tail -1|awk '{print \$3}'`;
                                chomp @bp[0];
                                $bp_end = `awk '{if(\$2==@bp[0] || \$3==@bp[0]){print}}' @bpData[0].@fileData[0].Index.file|head -1|awk '{print \$9}'`;
                                chomp $bp_end;
                                print out "@fileData[0]\t@fileData[1]\t$bp_end\t@fileData[3]\t@fileData[4]\n";
                                undef @bp;
                        }
                        elsif ($start > 20 && $end > $total){
                                $end_mod = $end-20;
                                @bp = `Rscript ../scripts/bp_assignment.r @bpData[0].@fileData[0].xy.optimal.max del $start $end_mod|tail -2|head -1|awk '{print \$3}'`;
                                chomp @bp[0];
                                $bp_start = `awk '{if(\$2==@bp[0] || \$3==@bp[0]){print}}' @bpData[0].@fileData[0].Index.file|head -1|awk '{print \$9}'`;
                                chomp $bp_start;
                                print out "@fileData[0]\t$bp_start\t@fileData[2]\t@fileData[3]\t@fileData[4]\n";
                                undef @bp;
                        }
                        elsif ($start > 20 && $end < $total){
                                @bp = `Rscript ../scripts/bp_assignment.r @bpData[0].@fileData[0].xy.optimal.max del $start $end|tail -2|awk '{print \$3}'`;
                                chomp (@bp[0],@bp[1]);
                                $bp_start = `awk '{if(\$2==@bp[0] || \$3==@bp[0]){print}}' @bpData[0].@fileData[0].Index.file|head -1|awk '{print \$9}'`;
                                $bp_end   = `awk '{if(\$2==@bp[1] || \$3==@bp[1]){print}}' @bpData[0].@fileData[0].Index.file|head -1|awk '{print \$9}'`;
                                chomp ($bp_start,$bp_end);
                                print out "@fileData[0]\t$bp_start\t$bp_end\t@fileData[3]\t@fileData[4]\n";
                                undef @bp;
                        }
                }
		undef @fileData;
		$j++;
	}
	close out;
	undef @bpData;
	undef @file;
	$i++;
}
