library(KernSmooth)
library(MASS)
library(RHmm)

#sample = read.table("list",header=F)
chr.list = read.table("../scripts/chr.list",header=F)
args <- commandArgs(trailingOnly = TRUE) # GC, MAP and Then Frag length
gc_limit = as.numeric(args[1])
map_limit = as.numeric(args[2])
frag_limit = as.numeric(args[3])
sample = data.frame(V1=as.character(args[4]))
i = 1
while (i <= length(sample$V1)){
        #sample.combined = read.table(as.character(paste(sample$V1[i],".F_GC_MAP.NormCount.",gc_limit,"_",map_limit,"_",frag_limit,".Index.bed",sep="")),header=T)
	sample.tpm = read.table(as.character(paste(sample$V1[i],".copyNumber.txt",sep="")),header=T)
        j = as.numeric(args[6])
        while (j <= as.numeric(args[6])){
		sample.combined = read.table(as.character(paste(sample$V1[i],".",chr.list$V1[j],".F_GC_MAP.NormCount.",gc_limit,"_",map_limit,"_",frag_limit,".Index.bed",sep="")),header=T)
                sample.filt = subset(sample.combined,chr==as.character(chr.list$V1[j]) & gc >= gc_limit & map >= map_limit & frag >= frag_limit)
		sample.chr.tpm = subset(sample.tpm,c_name==as.character(chr.list$V1[j]))
                count = sample.filt$norm.count*sample.chr.tpm$gain
                chr.name = sample.filt$chr
                start = as.integer(sample.filt$start)
                end   = as.integer(sample.filt$end)
                data  = data.frame(chr.name,start,end,count)
		x.dir = c(min(end),max(end))
		y.dir = c(min(count),max(count))
		range.x = list(x.dir,y.dir)
		if (as.numeric(args[5]) != 0){
			bw = as.numeric(args[5])  #dpik(end); 1e5; 2e6
		}else if (as.numeric(args[5]) == 0){
			bw = dpik(start)
		}
		l = bkde2D(data[,3:4], gridsize=c(as.integer(length(end)),length(unique(count))), bandwidth=c(bw,dpik(count)), range.x = range.x, truncate=TRUE)	
		pos.max = apply(l$fhat, 1, which.max)
		y = l$x2[pos.max]
		hmm.list = list()
		vitpath.list = list()
		bic.values = c()
		iter.num = c()
		na.num = 0
		
		## April 26, 2018 ##
		##########################################################		
		#k = 2
		#while (k <= 10){
		#	print (k)
		#	error  = try(HMMFit(y, nStates=k),silent=TRUE)
		#	if (!inherits(error, "try-error")){
	        #		ResFit = HMMFit(y, nStates=k)
        	#		VitPath = viterbi(ResFit, y)
        	#		hmm.list[[k]] = ResFit
        	#		vitpath.list[[k]] = VitPath
        	#		if (hmm.list[[k]]$convergence){
        	#       		bic.values = c(bic.values,ResFit$BIC)
        	#        		iter.num = c(iter.num,k)
        	#		}else {
		#			na.num = na.num+1
		#		}
		#		print (paste(k,ResFit$BIC,sep=""))
		#	}
       		#	k=k+1
		#}
		
		k = 2
                while (k <= 10){
                        ResFit = try(HMMFit(y, nStates=k),silent=TRUE)
                        if (!inherits(ResFit, "try-error")){
                                VitPath = viterbi(ResFit, y)
                                hmm.list[[k]] = ResFit
                                VitPath = viterbi(ResFit, y)
                                vitpath.list[[k]] = VitPath
                                if (hmm.list[[k]]$convergence){
                                        bic.values = c(bic.values,ResFit$BIC)
                                        iter.num = c(iter.num,k)
                                }else {
                                        na.num = na.num+1
                                }
                                print (paste(k,ResFit$BIC,sep=""))
                        } else {
                                na.num = na.num+1
                        }
                        k=k+1
                }
		##########################################################
		
		num.states = iter.num[which.min(bic.values)]
		if (na.num == 9 | length(bic.values) == 0){
			k = 2
			ResFit = HMMFit(y, nStates=k)
                      	VitPath = viterbi(ResFit, y)
                     	hmm.list[[k]] = ResFit
               		vitpath.list[[k]] = VitPath
			num.states = 2
		}
		state.mean = as.integer(vitpath.list[[num.states]]$states)
		m = 1
		while (m <= num.states){
        		state.mean[state.mean == m] = as.integer(hmm.list[[num.states]]$HMM$distribution$mean[m])
        		m=m+1
		}
		start.seg = as.integer(l$x1)
	        width.seg = diff(start.seg)
	        width.seg[length(start.seg)] = width.seg[1]
	        end.seg   = as.integer(start.seg+width.seg)
	        chr.name = rep(as.character(chr.list$V1[j]),length(start.seg))
	        seg.score = l$x2[pos.max]
		hmm_state = vitpath.list[[num.states]]$states
		hmm.state.new = c()
		state.mean.new = c()
		number_of_sites_up = c()	
		number_of_sites_dw = c()
		state.mean.adjusted = c()
		re_theoretic = (bw/2)/mean(diff(end))
		state.sd = c()
	        d = data.frame(chr.name,start.seg,end.seg,seg.score,hmm_state,state.mean)
		k = 1
		while (k <= length(d$start)){
			number_of_sites.tmp = sample.filt[sample.filt$start <= (sample.filt$start[k]+(bw/2)) & sample.filt$start >= (sample.filt$start[k]-(bw/2)),]$start
			number_of_sites.tmp = number_of_sites.tmp-sample.filt$start[k]
			number_of_sites_up = c(number_of_sites_up,length(number_of_sites.tmp[number_of_sites.tmp < 0]))
			number_of_sites_dw = c(number_of_sites_dw,length(number_of_sites.tmp[number_of_sites.tmp > 0]))

			k = k+1
		}
		up_diff = re_theoretic-number_of_sites_up
		dw_diff = re_theoretic-number_of_sites_dw
		d = cbind(d,number_of_sites_up,number_of_sites_dw,count)
		if (as.numeric(args[5]) == 0){
			bw = c("optimal")
		}
		bic = data.frame(iter.num,bic.values)
		#If there is no breakpoint, use mean normalized contact count to represent the chromosome.
		#if(is.na(d$state.mean)){d$state.mean = mean(d$count); d$seg.score = mean(d$count);}
		d$state.mean[is.na(d$state.mean)]= mean(d$count)
                d$seg.score[is.na(d$state.mean)] = mean(d$count)
		write.table(d,file=as.character(paste(sample$V1[i],".",chr.list$V1[j],".xy.",bw,".max",sep="")),row.names=F,col.names=F,quote=F,sep="\t")
		write.table(l$x1,file=as.character(paste(sample$V1[i],".",chr.list$V1[j],".x.",bw,".axis",sep="")),row.names=F,col.names=F,quote=F,sep="\t")
		write.table(l$x2,file=as.character(paste(sample$V1[i],".",chr.list$V1[j],".y.",bw,".axis",sep="")),row.names=F,col.names=F,quote=F,sep="\t")
		write.table(l$fhat,file=as.character(paste(sample$V1[i],".",chr.list$V1[j],".xy.",bw,".matrix",sep="")),row.names=F,col.names=F,quote=F,sep="\t")
                write.table(data,file=as.character(paste(sample$V1[i],".",chr.list$V1[j],".count.bed",sep="")),row.names=F,col.names=F,quote=F,sep="\t")
		j = j+1
        }
        i = i+1
}
