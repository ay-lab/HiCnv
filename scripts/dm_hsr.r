## Find DM from perREfragStats file

library(optparse)
option_list = list(

  make_option("--refeature", type="character", help="A five column restriction enzyme cutsite bed file with GC content, mappability, fragment length infomation\n
  \t\t<chr>\t<start>\t<end>\t<GC frequency>\t<mappability>\t<fragment length>\n"),
  make_option("--coverage", type="character", help="A .perREfragStats file\n"),
  make_option("--chrlist", type="character", help="A single column chromosome name file.\n"),
  make_option("--gccutoff", type="numeric", default=0.2, help="GC content cutoff. Anything below <gccutoff> will be removed [Default is 0.2].\n"),
  make_option("--mapcutoff", type="numeric", default=0.5, help="Mappability cutoff. Anything below <mapcutoff> will be removed [Default is 0.5].\n"),
  make_option("--fragcutoff", type="integer", default=150, help="Fragment length cutoff. Anything below <fragcutoff> will be removed [Default is 150].\n
  \t\tFor Hi-C experiments with 4bp cut enzyme, this value is 150, for 6bp enzymes this value should be 1000 [Optional].\n"),
  make_option("--prefix", type="character", help="All the files and folders will be created with this name.\n")
)

opt <- parse_args(OptionParser(option_list=option_list))

cmd1 <- paste0("sort -k 1,1 -k2,2n ",opt$refeature," > ",opt$refeature,".sorted.bed")
cmd2 <- paste0("sort -k 1,1 -k2,2n ",opt$coverage," > ",opt$coverage,".sorted.bedGraph")
cmd3 <- paste0("bedtools map -a ",opt$refeature,".sorted.bed -b ",opt$coverage,".sorted.bedGraph -c 5,6,7 -o sum -null 0 > ",opt$coverage,".signal.bedGraph")
print (cmd1)
system(cmd1, wait=T)
print (cmd2)
system(cmd2, wait=T)
print (cmd3)
system(cmd3, wait=T)

lane1 <- read.table(paste0(opt$coverage,".signal.bedGraph"),h=F,as.is=T)
chr <- read.table(opt$chrlist, h=F, as.is=T)

chrom <- lane1$V1
start <- lane1$V2
end   <- lane1$V3

mapped.count <- as.numeric(lane1$V7)
intra  <- as.numeric(lane1$V8)
inter  <- as.numeric(lane1$V9)

gc   <- lane1$V4
map  <- lane1$V5
frag <- lane1$V6

gc_limit   <- opt$gccutoff
map_limit  <- opt$mapcutoff
frag_limit <- opt$fragcutoff

r <- c()
data <- data.frame(chrom,start,end,mapped.count,intra,inter,gc,map,frag)
head(data)

i = 1
while (i <= length(chr$V1)){
	filt = subset(data,chrom==as.character(chr$V1[i]) & gc >= gc_limit & map >= map_limit & frag >= frag_limit)
	filt$mapped.count[filt$mapped.count==0] = median(filt$mapped.count)
	filt$intra[filt$intra==0] = median(filt$intra)
	filt$inter[filt$inter==0] = median(filt$inter)
	mapped_norm = filt$mapped.count/median(filt$mapped.count)
	intra_norm = filt$intra/median(filt$intra)
	inter_norm = filt$inter/median(filt$inter)
	filt = cbind(filt,mapped_norm,inter_norm)
	r = rbind(r,filt)
	i = i+1
}
mapped_norm_pval <- p.adjust(1-ppois(r$mapped_norm,mean(r$mapped_norm)))
inter_norm_pval  <- p.adjust(1-ppois(r$inter_norm,mean(r$inter_norm)))
mapped_norm_pval[is.na(mapped_norm_pval)] <- 1
inter_norm_pval[is.na(inter_norm_pval)]   <- 1
r <- cbind(r,mapped_norm_pval,inter_norm_pval)
hsr <- r[r$inter_norm_pval > 0.05 & r$mapped_norm_pval < 0.05,]
dm  <- r[r$inter_norm_pval < 0.05 & r$mapped_norm_pval < 0.05,]
dm  <- dm[table(dm$chrom) > 3,]
hsr <- hsr[table(hsr$chrom) > 3,]
write.table(dm,file=paste0("dm.",opt$prefix,".pval.txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
write.table(hsr,file=paste0("hsr.",opt$prefix,".pval.txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
