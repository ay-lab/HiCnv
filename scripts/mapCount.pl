#open (inA,"T470.L1_hg19.perREfragStats");
open (inA,"@ARGV[0].perREfragStats");
while (<inA>){
	chomp $_;
	@data = split(/\s+/,$_);
	$count{@data[0]}{@data[1]}{@data[2]} = @data[3]; 
	undef @data;
}
close inA;
open (out,">@ARGV[0].F_GC_MAP.bed");
open (inB,"@ARGV[1]");
while (<inB>){
	chomp $_;
	@data = split(/\s+/,$_);
	if ($count{@data[0]}{@data[1]}{@data[2]} ne ""){
		print out "$_\t$count{@data[0]}{@data[1]}{@data[2]}\n";
	}
	undef @data;
}
close inB;
close out;
undef %count;
