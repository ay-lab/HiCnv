$frag_file_in  = @ARGV[0];
$frag_file_out = @ARGV[1]; 
open (inA,"$frag_file_in");
while (<inA>){
	chomp $_;
	@data = split(/\s+/,$_);
	if (@data[16] ne "."){
		$gc_sum{@data[0]}{@data[4]}{@data[5]} += @data[8];
		$gc_num{@data[0]}{@data[4]}{@data[5]} += 1;
		$map_sum{@data[0]}{@data[4]}{@data[5]} += @data[16];
		$map_num{@data[0]}{@data[4]}{@data[5]} += 1;
	}
	undef @data;
}
close inA;

open (inB,"$frag_file_out");
while (<inB>){
	chomp $_;
	@data = split(/\s+/,$_);
	if ($gc_num{@data[0]}{@data[1]}{@data[2]} == 0){$gc_num{@data[0]}{@data[1]}{@data[2]} = 1;}
	if ($map_num{@data[0]}{@data[1]}{@data[2]} == 0){$map_num{@data[0]}{@data[1]}{@data[2]} = 1;}
	$gc_mean = sprintf("%0.3f",$gc_sum{@data[0]}{@data[1]}{@data[2]}/$gc_num{@data[0]}{@data[1]}{@data[2]});
	$map_mean = sprintf("%0.3f",$map_sum{@data[0]}{@data[1]}{@data[2]}/$map_num{@data[0]}{@data[1]}{@data[2]});
	$frag = @data[2]-@data[1];
	print "@data[0]\t@data[1]\t@data[2]\t$gc_mean\t$map_mean\t$frag\n";
	undef @data;
}
close inB;
