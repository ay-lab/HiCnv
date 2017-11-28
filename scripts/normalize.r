#sample = read.table("sample.list",header=F)
chr = read.table("../scripts/chr.list",header=F)
args <- commandArgs(trailingOnly = TRUE) # GC, MAP and Then Frag length
sample = data.frame(V1=as.character(args[1]))
gc_limit = as.numeric(args[2])
map_limit = as.numeric(args[3])
frag_limit = as.numeric(args[4])
pdf(paste(as.character(args[1]),".",gc_limit,"_",map_limit,"_",frag_limit,".pdf",sep=""))
i = 1
while (i <= length(sample$V1)){
	#sample.Rep1 = read.table(as.character(paste(sample$V1[i],".R1.F_GC_MAP.bed",sep="")),header=F)
	#sample.Rep2 = read.table(as.character(paste(sample$V1[i],".R2.F_GC_MAP.bed",sep="")),header=F)
	sample.Rep1 = read.table(as.character(paste(sample$V1[i],".F_GC_MAP.bed",sep="")),header=F)
	chrom.data  = c()
	bic.data = c()
	j = 1
	while (j <= length(chr$V1)){
		#sample.R1 = subset(sample.Rep1,V1==as.character(chr$V1[j]) & V4 >= gc_limit & V5 >= map_limit & V6 >= frag_limit)
		#sample.R2 = subset(sample.Rep2,V1==as.character(chr$V1[j]) & V4 >= gc_limit & V5 >= map_limit & V6 >= frag_limit)
		#count = as.integer((sample.R1$V7+sample.R2$V7)/2)
		sample.R1 = subset(sample.Rep1,V1==as.character(chr$V1[j]) & V4 >= gc_limit & V5 >= map_limit & V6 >= frag_limit)
		count = sample.R1$V7
		gcc_vec = sample.R1$V4
		map_vec = sample.R1$V5
		gcc_vec[gcc_vec==0] = 0.001
		map_vec[map_vec==0] = 0.001
		len_vec = log(sample.R1$V6)
	        gcc_vec = log(gcc_vec)
	        map_vec = log(map_vec)
		len_vec = (len_vec-mean(c(len_vec)))/sd(c(len_vec))
		gcc_vec = (gcc_vec-mean(c(gcc_vec)))/sd(c(gcc_vec))
		map_vec = (map_vec-mean(c(map_vec)))/sd(c(map_vec))
		fit = glm(count ~ len_vec+gcc_vec+map_vec,family="poisson")
		coeff = round(fit$coeff,4)
		res = round(count/exp(coeff[1]+coeff[2]*len_vec+coeff[3]*gcc_vec+coeff[4]*map_vec), 2)*100
		data = data.frame(chr=sample.R1$V1,start=sample.R1$V2,end=sample.R1$V3,gc=sample.R1$V4,map=sample.R1$V5,frag=sample.R1$V6,count=count,norm.count=res)
		filt.out = subset(sample.Rep1,V1==as.character(chr$V1[j]) & (V4 < gc_limit | V5 < map_limit | V6 < frag_limit))
		filt.out.count = rep(0,length(filt.out$V1))
		filt.out = cbind(filt.out,filt.out.count)
		final.data = data.frame(chr=c(as.character(data$chr),as.character(filt.out$V1)),start=c(data$start,filt.out$V2),end=c(data$end,filt.out$V3),gc=c(data$gc,filt.out$V4),map=c(data$map,filt.out$V5),frag=c(data$frag,filt.out$V6),count=c(data$count,filt.out$V7),norm.count=c(data$norm.count,filt.out$filt.out.count))
		final.data = final.data[order(final.data$start),]
		chrom.data = rbind(chrom.data,final.data)
		sample.name = as.character(sample$V1[i])
		chr.name = as.character(chr$V1[j])
		gc.limit = gc_limit
		map.limit = map_limit
		frag.limit = frag_limit
		bic.val = fit$aic
		bic.data.tmp = cbind(sample.name,chr.name,gc.limit,map.limit,frag.limit,bic.val)
		bic.data = rbind(bic.data,bic.data.tmp)
		j = j+1
	}
	write.table(chrom.data,file=as.character(paste(sample$V1[i],".F_GC_MAP.NormCount.",gc_limit,"_",map_limit,"_",frag_limit,".bed",sep="")),row.names=F,quote=F,sep="\t")	
	#rite.table(bic.data,file=as.character(paste(sample$V1[i],".F_GC_MAP.AIC.",gc_limit,"_",map_limit,"_",frag_limit,".bed",sep="")),row.names=F,quote=F,sep="\t")
	count.log = log10(chrom.data$count/mean(chrom.data$count))
	norm.count.log = log10(chrom.data$norm.count/mean(chrom.data$norm.count))
	par(mfrow=c(3,2))
	smoothScatter(chrom.data$gc,count.log,colramp = colorRampPalette(c("white", "azure", "azure2", "blue", "green", "yellow", "red")), nrpoints = 0, nbins=1000, ylim=c(-1,1), ylab="Raw Count distribution", xlab="GC%", main=as.character(paste(sample$V1[i])))
	lines(chrom.data$gc,rep(0,length(chrom.data$gc)))
	smoothScatter(chrom.data$gc,norm.count.log,colramp = colorRampPalette(c("white", "azure", "azure2", "blue", "green", "yellow", "red")), nrpoints = 0, nbins=1000, ylim=c(-1,1), ylab="Normalized Count distribution", xlab="GC%", main=as.character(paste(sample$V1[i])))
	lines(chrom.data$gc,rep(0,length(chrom.data$gc)))
	smoothScatter(chrom.data$map,count.log,colramp = colorRampPalette(c("white", "azure", "azure2", "blue", "green", "yellow", "red")), nrpoints = 0, nbins=1000, ylim=c(-1,1), ylab="Raw Count distribution", xlab="Mappability", main=as.character(paste(sample$V1[i])))
	lines(chrom.data$map,rep(0,length(chrom.data$map)))
	smoothScatter(chrom.data$map,norm.count.log,colramp = colorRampPalette(c("white", "azure", "azure2", "blue", "green", "yellow", "red")), nrpoints = 0, nbins=1000, ylim=c(-1,1), ylab="Normalized Count distribution", xlab="Mappability", main=as.character(paste(sample$V1[i])))
	lines(chrom.data$map,rep(0,length(chrom.data$map)))
	smoothScatter(log10(chrom.data$frag),count.log,nrpoints = 0,colramp = colorRampPalette(c("white", "azure", "azure2", "blue", "green", "yellow", "red")), nbins=1000, ylim=c(-1,1), ylab="Raw Count distribution", xlab="Log10 Fragment size", main=as.character(paste(sample$V1[i])))
	lines(log10(chrom.data$frag),rep(0,length(chrom.data$frag)))
	smoothScatter(log10(chrom.data$frag),norm.count.log,nrpoints = 0,colramp = colorRampPalette(c("white", "azure", "azure2", "blue", "green", "yellow", "red")), nbins=1000, ylim=c(-1,1), ylab="Normalized Count distribution", xlab="Log10 Fragment size", main=as.character(paste(sample$V1[i])))
	lines(log10(chrom.data$frag),rep(0,length(chrom.data$frag)))
	i = i+1	
}
dev.off()
