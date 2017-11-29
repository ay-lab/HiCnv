args <- commandArgs(trailingOnly = TRUE)
lane1 = read.table(paste(as.character(args[1]),".perREfragStats",sep=""),header=F)
#lane2 = read.table(paste(as.character(args[2]),".perREfragStats",sep=""),header=F)
hindIII = read.table(paste(as.character(args[2]),sep=""),header=F)
chr = read.table("../scripts/chr.list",header=F)
chrom = lane1$V1
start = lane1$V2
end   = lane1$V3
#mapped.count = (lane1$V4+lane2$V4)/2
#intra  = (lane1$V5+lane2$V5)/2
#inter  = (lane1$V6+lane2$V6)/2
mapped.count = lane1$V4
intra  = lane1$V5
inter  = lane1$V6
gc = hindIII$V4
map = hindIII$V5
frag = hindIII$V6
r = c()
gc_limit = as.numeric(args[4])
map_limit = as.numeric(args[5])
frag_limit = as.numeric(args[6])
data = data.frame(chrom,start,end,mapped.count,intra,inter,gc,map,frag)
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
mapped_norm_pval = p.adjust(1-ppois(r$mapped_norm,mean(r$mapped_norm)))
inter_norm_pval  = p.adjust(1-ppois(r$inter_norm,mean(r$inter_norm)))
mapped_norm_pval[is.na(mapped_norm_pval)] = 1
inter_norm_pval[is.na(inter_norm_pval)] = 1
r = cbind(r,mapped_norm_pval,inter_norm_pval)
hsr = r[r$inter_norm_pval > 0.05 & r$mapped_norm_pval < 0.05,]
dm  = r[r$inter_norm_pval < 0.05 & r$mapped_norm_pval < 0.05,]
dm  = dm[table(dm$chrom) > 3,]
hsr = hsr[table(hsr$chrom) > 3,]
write.table(dm,file=paste("dm.",as.character(args[3]),".pval.txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
write.table(hsr,file=paste("hsr.",as.character(args[3]),".pval.txt",sep=""),row.names=F,col.names=T,quote=F,sep="\t")
