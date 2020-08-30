#This script will combine all the replicates to generate a signle combined replicate
# 

my $path   = $ARGV[0];
my $prefix = $ARGV[1];
chomp ($path, $prefix);

@map_list  = `ls $path/*.F_GC_MAP.bed`;
@frag_list = `ls $path/*.perREfragStats`;

$i = 0;
while ($i <= $#map_list){
	chomp $map_list[$i];
	open (in,"$map_list[$i]");
	while (<in>){
		chomp $_;
		@data = split(/\s+/,$_);
		if ($map_list_check{@data[0]}{@data[1]}{@data[2]} eq ""){
			$map_list_check{@data[0]}{@data[1]}{@data[2]} = "yes";
			$map_list_gc{@data[0]}{@data[1]}{@data[2]}  = @data[3];
			$map_list_map{@data[0]}{@data[1]}{@data[2]} = @data[4];
			$map_list_len{@data[0]}{@data[1]}{@data[2]} = @data[5];
			push @map_list_re,"@data[0]\t@data[1]\t@data[2]\n";
		}
		$map_list_count{@data[0]}{@data[1]}{@data[2]} += @data[6];
		undef @data;
	}
	close in;
	$i++;
}

open (out,">$path/$prefix.combined.F_GC_MAP.bed");
$i = 0;
while ($i <= $#map_list_re){
	chomp $map_list_re[$i];
	@map_list_re_data = split(/\s+/,$map_list_re[$i]);
	$chr   = @map_list_re_data[0];
	$start = @map_list_re_data[1];
	$end   = @map_list_re_data[2];
	$gc    = $map_list_gc{$chr}{$start}{$end};
	$map   = $map_list_map{$chr}{$start}{$end};
	$len   = $map_list_len{$chr}{$start}{$end};
	$count = int($map_list_count{$chr}{$start}{$end}/($#map_list+1));
        if ($chr ne "chrM" && $chr ne "chrY") {
  	  print out "$chr\t$start\t$end\t$gc\t$map\t$len\t$count\n";
        }
	undef @map_list_re_data;
	$i++;
}
close out;
undef %map_list_check;
undef %map_list_gc;
undef %map_list_map;
undef %map_list_len;
undef %map_list_count;
undef @map_list_re;
undef @map_list;

$i = 0;
while ($i <= $#frag_list){
        chomp $frag_list[$i];
        open (in,"$frag_list[$i]");
        while (<in>){
                chomp $_;
                @data = split(/\s+/,$_);
                if ($frag_list_check{@data[0]}{@data[1]}{@data[2]} eq ""){
                        $frag_list_check{@data[0]}{@data[1]}{@data[2]} = "yes";
                        push @frag_list_re,"@data[0]\t@data[1]\t@data[2]\n";
                }
		$j = 3;
		while ($j <= 11){
			chomp $data[$j];
			$frag_list_count{$j}{@data[0]}{@data[1]}{@data[2]} += $data[$j];
			$j++;
		}
                undef @data;
        }
        close in;
        $i++;
}

open (out,">$path/$prefix.combined.perREfragStats");
$i = 0;
while ($i <= $#frag_list_re){
        chomp $frag_list_re[$i];
        @frag_list_re_data = split(/\s+/,$frag_list_re[$i]);
        $chr   = @frag_list_re_data[0];
        $start = @frag_list_re_data[1];
        $end   = @frag_list_re_data[2];
        if ($chr ne "chrM" && $chr ne "chrY") {
	  print out "$chr\t$start\t$end\t.";
	  $j = 3;
	  while ($j <= 11){
	        $count = int($frag_list_count{$j}{$chr}{$start}{$end}/($#frag_list+1));
		print out "\t$count";
		$j++;
	  } 
          print out "\n";
        }
        undef @frag_list_re_data;
        $i++;
}
close out;
undef %frag_list_check;
undef %frag_list_count;
undef @frag_list_re;
undef @frag_list;
