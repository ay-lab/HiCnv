args <- commandArgs(trailingOnly = TRUE)
sample = as.character(args[1])
print (sample)
bw = c("optimal")
chr = read.table(as.character(paste(sample,".copyNumber.txt",sep="")),header=T)
zscore = c()
chr.name = c()
start = c()
end = c()
score = c()
data = c()
i = 1 
while (i <= length(chr$c_name)){
	print (chr$c_name[i])
	x = read.table(as.character(paste(sample,".",chr$c_name[i],".x.",bw,".axis",sep="")),header=F)
	y = read.table(as.character(paste(sample,".",chr$c_name[i],".y.",bw,".axis",sep="")),header=F)
	z = as.matrix(read.table(as.character(paste(sample,".",chr$c_name[i],".xy.",bw,".matrix",sep="")),header=F))
	b = read.table(as.character(paste(sample,".",chr$c_name[i],".Index.file",sep="")),header=F)
	sub = read.table(as.character(paste(sample,".",chr$c_name[i],".xy.",bw,".max",sep="")),header=F)
	ref = read.table(as.character(paste(sample,".",chr$ref.chrom[i],".xy.",bw,".max",sep="")),header=F)
        ref.mean = mean(ref$V6)
	zscore = (sub$V6-ref.mean)/100
	chr.name = as.character(sub$V1)
	upper = ref.mean+(ref.mean*0.45)
        lower = ref.mean-(ref.mean*0.45)
	start = b$V9
	end   = b$V10
	score = sub$V6
	d = cbind(chr.name,start,end,zscore,score,ref.mean)
	data = rbind(data,d)
	png(filename=as.character(paste(sample,".",chr$c_name[i],".chrom.density.png",sep="")),width=8000,height=7000,res=1000)
        pos = start/1e6
        filled.contour(pos, y$V1, z, ylim=range(0,quantile(y$V1,0.99)), 
	plot.axes={
	#axis(1);axis(2);points(pos,sub$V4,type="p",cex=0.5);
        axis(1);axis(2);points(pos,score,type="l",lwd=4,col="black");
        axis(1);axis(2);points(pos,rep(ref.mean,length(sub$V1)),type="l",lwd=2,col="red")
        axis(1);axis(2);points(pos,rep(upper,length(sub$V1)),type="l",lwd=2,lty=2,col="red")
        axis(1);axis(2);points(pos,rep(lower,length(sub$V1)),type="l",lwd=2,lty=2,col="red")
	},color.palette=topo.colors)
        dev.off()
	i = i+1
}
write.table(data,file=as.character(paste(sample,".chrWise.txt",sep="")),row.names=F,col.names=F,quote=F,sep="\t")
