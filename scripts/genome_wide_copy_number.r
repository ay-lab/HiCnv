library(ggplot2)
library(GenVisR)
args <- commandArgs(trailingOnly = TRUE)
chr = read.table("../scripts/chr.list",header=F)
cp_perc = read.table(as.character(args[1]),header=F)
i = 1
while (i <= length(chr$V1)){
	cp_perc_chr = subset(cp_perc,V1==as.character(chr$V1[i]))
	cn = 2+(2*cp_perc_chr$V4)
	coordinate = as.integer(cp_perc_chr$V2+((cp_perc_chr$V3-cp_perc_chr$V2)/2))
	data = data.frame(chromosome=as.character(cp_perc_chr$V1), coordinate, cn)
	pdf(paste(as.character(chr$V1[i]),".copy_number.pdf",sep=""))
	plot_theme <- ylim(0, 5)
        cnView(data, chr = as.character(chr$V1[i]), genome = "hg19", ideogram_txtSize = 2, plotLayer = plot_theme)
	write.table(data,file=paste(as.character(chr$V1[i]),"_copynumber.txt",sep=""), row.name=F, quote=F, sep="\t")
	dev.off()
	i = i+1
}
