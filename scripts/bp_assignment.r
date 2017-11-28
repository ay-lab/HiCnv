args <- commandArgs(trailingOnly = TRUE)
sample = read.table(as.character(args[1]),header=F)
type = as.character(args[2])
start_boundary = c()
end_boundary = c()
override = max(sample$V8)+50
if (type == "del"){
	start_pos = as.numeric(args[3])
	sample_start = subset(sample,V2 >= start_pos-11 & V2 <= start_pos+(sample[sample$V2==start_pos,]$V8+override)+11)
	i = 11
	while (i <= length(sample_start$V1)-11){
		count_up = sample_start[(i-10):(i-1),]$V9
		count_dw = sample_start[(i+1):(i+10),]$V9
		m_up = median(count_up)
		m_dw = median(count_dw)
		r = m_up/m_dw
		x = data.frame(chrom=as.character(sample_start$V1[i]),start=sample_start$V2[i],end=sample_start$V3[i],count=sample_start$V9[i],ratio=r)
		start_boundary = rbind(start_boundary,x)
		i = i+1
	}
	
	end_pos = as.numeric(args[4])
	sample_end = subset(sample,V2 <= end_pos+11 & V2 >= end_pos-(sample[sample$V2==end_pos,]$V7+override)+11)
	i = 11
	while (i <= length(sample_end$V1)-11){
	        count_up = sample_end[(i-10):(i-1),]$V9
	        count_dw = sample_end[(i+1):(i+10),]$V9
	        m_up = median(count_up)
	        m_dw = median(count_dw)
	        r = m_dw/m_up
	        x = data.frame(chrom=as.character(sample_end$V1[i]),start=sample_end$V2[i],end=sample_end$V3[i],count=sample_end$V9[i],ratio=r)
	        end_boundary = rbind(end_boundary,x)
	        i = i+1
	}
}else if (type == "amp"){
	start_pos = as.numeric(args[3])
	sample_start = subset(sample,V2 <= start_pos+11 & V2 >= start_pos-(sample[sample$V2==start_pos,]$V7+override)+11)
	i = 11
        while (i <= length(sample_start$V1)-11){
                count_up = sample_start[(i-10):(i-1),]$V9
                count_dw = sample_start[(i+1):(i+10),]$V9
                m_up = median(count_up)
                m_dw = median(count_dw)
                r = m_dw/m_up
                x = data.frame(chrom=as.character(sample_start$V1[i]),start=sample_start$V2[i],end=sample_start$V3[i],count=sample_start$V9[i],ratio=r)
                start_boundary = rbind(start_boundary,x)
                i = i+1
        }
	
	end_pos = as.numeric(args[4])
        sample_end = subset(sample,V2 >= end_pos-11 & V2 <= end_pos+(sample[sample$V2==end_pos,]$V8+override)+11)
        i = 11
        while (i <= length(sample_end$V1)-11){
                count_up = sample_end[(i-10):(i-1),]$V9
                count_dw = sample_end[(i+1):(i+10),]$V9
                m_up = median(count_up)
                m_dw = median(count_dw)
                r = m_up/m_dw
                x = data.frame(chrom=as.character(sample_end$V1[i]),start=sample_end$V2[i],end=sample_end$V3[i],count=sample_end$V9[i],ratio=r)
                end_boundary = rbind(end_boundary,x)
                i = i+1
        }
}

start_boundary = start_boundary[order(-start_boundary$ratio),]
end_boundary = end_boundary[order(-end_boundary$ratio),]
start_boundary = start_boundary[1,]
end_boundary = end_boundary[1,]
d = rbind(start_boundary,end_boundary)
print (d)
