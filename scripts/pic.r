library(ineq)
args <- commandArgs(trailingOnly = TRUE)
all.chr = read.table("../scripts/chr.list",header=F)
colnames(all.chr) <- c("chrom")
gc_limit = as.numeric(args[4])
map_limit = as.numeric(args[5])
frag_limit = as.numeric(args[6])
copyEstimate_with_ref <- function(chr.in,drop=TRUE) {
	chr.name <- chr.in[[1]]
	count <- chr.in[[2]]
	ref.chrom <- chr.in[[3]]
	g_name = c()
	m_name = c()
        c_name = c()
	len = c()
	rpm = c()
	data = c()
	i = 1
	while (i <= length(chr.name$chrom)){
		d.chr = subset(count,chr==as.character(chr.name$chrom[i]) & gc >= gc_limit & map >= map_limit & frag >= frag_limit)
		g_name = c(g_name,ineq(d.chr$count,type="Gini"))
		m_name = c(m_name,sum(d.chr$count))
		x = length(d.chr$frag)*1000
		rpm = c(rpm,(sum(d.chr$count)/(x/1e6)))
		c_name = c(c_name,as.character(chr.name$chrom[i]))
		i = i+1
	}
	tpm = round(rpm/sum(rpm),4)
	data = data.frame(c_name,g_name,tpm)
	gain = round(tpm/subset(data,c_name==ref.chrom)$tpm,1)
	gain[gain==0] = 1.0
	copy.num = 2+(2*(gain-1))
	data = cbind(data,gain,copy.num,m_name,rpm)
	return(data)
}

sample = as.character(args[1])
copyNumber = c()
d = read.table(as.character(paste(sample,".F_GC_MAP.NormCount.",gc_limit,"_",map_limit,"_",frag_limit,".bed",sep="")),header=T)
if (as.numeric(args[2]) == 0){
	l = list()
	ref.chr = as.character(args[3])
	l[[1]] = all.chr
	l[[2]] = d
	l[[3]] = ref.chr
	ref.chrom = rep(ref.chr,length(all.chr$chrom))
	large.copyEstimate = copyEstimate_with_ref(l)
	large.copyEstimate = cbind(large.copyEstimate,ref.chrom)
	copyNumber = rbind(large.copyEstimate)
	sample.name = rep(sample,length(copyNumber$c_name))
	copyNumber  = cbind(sample.name,copyNumber)
	copyNumber = copyNumber[order(-1*copyNumber$copy.num),]
	write.table(copyNumber,file=as.character(paste(sample,".copyNumber.txt",sep="")),row.names=F,quote=F,sep="\t")
}else if (as.numeric(args[2]) == 1){
	l = list()
	ref.chr = as.character("chr1")
   	l[[1]] = all.chr
   	l[[2]] = d
     	l[[3]] = ref.chr
	large.copyEstimate = copyEstimate_with_ref(l)
	m = matrix(large.copyEstimate$tpm,nrow=length(large.copyEstimate$tpm), ncol=length(large.copyEstimate$tpm), byrow=T)
	g = ((m/large.copyEstimate$tpm)*2)-2
	k = apply(g[1:23,],1,sd)
	p = k/large.copyEstimate$g_name
	Group = rep("large",length(large.copyEstimate$c_name))
	ratio = large.copyEstimate$tpm/(1/23)
	large.copyEstimate = cbind(large.copyEstimate,k,p,Group,ratio)
	if (as.character(args[3]) == "normal"){
		ref.chr = large.copyEstimate[order(large.copyEstimate$k),]$c_name
	}else if (as.character(args[3]) == "cancer"){
        	#ref.chr = large.copyEstimate[order(-1*large.copyEstimate$p),]$c_name
		ref.chr = large.copyEstimate[order(-1*large.copyEstimate$k),]$c_name
		large.copyEstimate = large.copyEstimate[order(-1*large.copyEstimate$k),]
      	}
	if (ref.chr[1]=="chrX"){
		ref.chr = ref.chr[2]
	}else {
		ref.chr = ref.chr[1]	
	}
	copyNumber = rbind(large.copyEstimate)
	copyNumber = subset(copyNumber, selec = -c(gain,copy.num))
 	sample.name = rep(sample,length(copyNumber$c_name))
	write.table(copyNumber,file=as.character(paste(sample,".PreEstimate.txt",sep="")),row.names=F,quote=F,sep="\t")
	print (paste("Please Check",sample,".PreEstimate.txt file to select the reference",sep=""))
}
