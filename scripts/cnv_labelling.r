args <- commandArgs(trailingOnly = TRUE)
sample = as.character(args[1])
chr = read.table("../scripts/chr.list",header=F)
r = read.table(paste(sample,".chrWise.txt",sep=""),header=F)
perc = as.numeric(args[2])
bed.final = c()
bed.full  = c()
k = 1
while (k <= length(chr$V1)){
	r.filt = subset(r,V1==as.character(chr$V1[k]))

	chrom.name = c()
	chrom.start = c()
	chrom.end = c()
	chrom.perc = c()
	chrom.val = c()
	chrom.val.mean = c()
	chrom.cnv = c()

	nrm = 0
	nrm_chk = 0
	del = 0
	del_chk = 0
	amp = 0
	amp_chk = 0
	cnv = c()
	i = 1
	while (i <= length(r.filt$V1)){
		if (r.filt$V4[i] > -1*perc & r.filt$V4[i] < perc){
			if (nrm_chk == 0){
				nrm = nrm+1
				nrm_chk = 1		
				del_chk = 0
				amp_chk = 0
			}
			cnv = c(cnv,paste("Normal_",nrm,sep=""))						
		}else if (r.filt$V4[i] <= -1*perc){
			if (del_chk == 0){
				del = del+1
				del_chk = 1
				nrm_chk = 0
				amp_chk = 0
			}
			cnv = c(cnv,paste("Deletion_",del,sep=""))
		}else if (r.filt$V4[i] >= perc){
			if (amp_chk == 0){
				amp = amp+1
				amp_chk = 1
				nrm_chk = 0
				del_chk = 0
			}
			cnv = c(cnv,paste("Amplification_",amp,sep=""))
		}
		i = i+1
	}
	r.filt = cbind(r.filt,cnv)
	pos.start = c()
	pos.end   = c()
	m.g = c()
	cnv.type = c()
	segments = c()
	i = 1
	while (i <= length(levels(r.filt$cnv))){
		v = subset(r.filt,cnv==levels(r.filt$cnv)[[i]])
		m.g = c(m.g,mean(v$V4))
		pos.start = c(pos.start,v$V2[1])
		pos.end = c(pos.end,v$V3[length(v$cnv)])
		cnv.type = c(cnv.type,levels(r.filt$cnv)[[i]])
		i = i+1
	}
	i = 1
        while (i <= length(levels(r.filt$cnv))){
                v = subset(r.filt,cnv==levels(r.filt$cnv)[[i]])
                chrom.name = c(chrom.name,as.character(v$V1))
                chrom.start = c(chrom.start,v$V2)
                chrom.end = c(chrom.end,v$V3)
		chrom.perc = c(chrom.perc,v$V4)
		chrom.val = c(chrom.val,v$V5)
		v.tmp = v[which(diff(v$V5) >= sd(r.filt$V5)),]
		if (length(v.tmp$V1) == 0){
			x = mean(v$V5)
			chrom.val.mean = c(chrom.val.mean,rep(x,length(v$V5)))
		}else {
			s = v$V2[1]		
			e = 1
			while (e <= length(v.tmp$V1)){
				x.tmp = mean(subset(v,V2 >= s & V2 < v.tmp$V2[e])$V5)
				chrom.val.mean = c(chrom.val.mean,rep(x.tmp,length(subset(v,V2 >= s & V2 < v.tmp$V2[e])$V5)))
				s = v.tmp$V2[e]
				e = e+1
			}	
			x.tmp = mean(subset(v,V2 >= s & V3 <= max(v$V3))$V5)
			chrom.val.mean = c(chrom.val.mean,rep(x.tmp,length(subset(v,V2 >= s & V3 <= max(v$V3))$V5)))
		}
		chrom.cnv = c(chrom.cnv,as.character(v$cnv))
                i = i+1
        }

	chromosome = rep(as.character(r.filt$V1[1]),length(levels(r.filt$cnv)))
	bed.a = data.frame(chromosome,pos.start,pos.end,m.g,cnv.type)
	bed.a = bed.a[order(bed.a$pos.start),]
	bed.final = rbind(bed.final,bed.a)
	
	tmp = data.frame(chrom.name,chrom.start,chrom.end,chrom.perc,chrom.val,chrom.val.mean,chrom.cnv)
	tmp = tmp[order(tmp$chrom.start),]
	bed.full = rbind(bed.full,tmp)

	k = k+1
}
write.table(bed.final,file=as.character(paste(sample,".merged_",perc,".WithOutFilt.bed",sep="")),row.names=F,col.names=F,quote=F,sep="\t")
