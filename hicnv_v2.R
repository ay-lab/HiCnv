################################################################################
## HiCnv scans Hi-C data to find copy number variations 
## Authors : Abhijit Chakraborty (abhijit@lji.org) & Ferhat Ay (ferhatay@lji.org)
## DOI: 10.1093/bioinformatics/btx664
################################################################################

#### Get the options ###
library(optparse)
library(data.table)
library(parallel)
library(EbayesThresh)

option_list = list(
  
  make_option("--refeature", type="character", help="A five column restriction enzyme cutsite bed file with GC content, mappability, fragment length infomation\n
  \t\t<chr>\t<start>\t<end>\t<GC frequency>\t<mappability>\t<fragment length>\n"),
  make_option("--coverage", type="character", help="A bedGraph file with read coverage signal (5th column).\n
  \t\tAlternatively, the .perREfragStats file\n"),
  make_option("--gccutoff", type="numeric", default=0.2, help="GC content cutoff. Anything below <gccutoff> will be removed [Default is 0.2].\n"),
  make_option("--mapcutoff", type="numeric", default=0.5, help="Mappability cutoff. Anything below <mapcutoff> will be removed [Default is 0.5].\n"),
  make_option("--fragcutoff", type="integer", default=150, help="Fragment length cutoff. Anything below <fragcutoff> will be removed [Default is 150].\n
  \t\tFor Hi-C experiments with 4bp cut enzyme, this value is 150, for 6bp enzymes this value should be 1000.\n"),
  make_option("--refchrom", type="character", default=NA, help="Name of the reference chromosome for CNV detection.\n
  \t\tIf no name is provided then HiCnv proceed to estimate proportion of interaction count or PIC to decide a reference chromosome.\n"),
  make_option("--bandwidth", type="integer", default=1e6, help="Genomic distance bandwidth with which Kernel smoothing will be performed [Default 1Mb].\n"),
  make_option("--hmmstate", type="integer", default=10, help="Number of HMM states to be searched in each chromosome [Default 10].\n"),
  make_option("--threshold", type="numeric", default=0.20, help="Threshold value to define amplification and deletion with respect to normal region [Default is 0.2 i.e. deviation of 20% from mean normal value will be labeled as CNV].\n"),
  make_option("--prefix", type="character", help="All the files and folders will be created with this name.\n"),
  make_option("--cpu", type="integer", default=1, help="Number of threads to be used [Default 1].\n")
)

opt <- parse_args(OptionParser(option_list=option_list))
refeature <- normalizePath(opt$refeature)
coverage  <- normalizePath(opt$coverage)
prefix    <- opt$prefix
gc_cutoff <- opt$gccutoff
map_coutoff <- opt$mapcutoff
frag_cutoff <- opt$fragcutoff
refchrom  <- opt$refchrom
genome_width <- opt$bandwidth
cnv_threshold <- opt$threshold
hmm_state <- opt$hmmstate

if (opt$cpu > 1) {
  cl <- makeCluster(opt$cpu)
} else { 
  cl <- NA
}

## Parameters provided ##
cat ("________________________________________\n")
cat ("Parameters provided\n")
cat ("________________________________________\n")
cat ("Restriction enzyme bed file:\t",refeature,"\n")
cat ("Read coverage file:\t",coverage,"\n")
cat ("Prefix:\t",prefix,"\n")
cat ("GC cutoff:\t",gc_cutoff,"\n")
cat ("Mappability cutoff:\t",map_coutoff,"\n")
cat ("Fragment length cutoff:\t",frag_cutoff,"\n")
cat ("Reference chromosome:\t",refchrom,"\n")
cat ("Genomic bandwidth:\t",genome_width,"\n")
cat ("CNV threshold:\t",cnv_threshold*100,"% +/- Normal signal\n")
cat ("HMM state:\t",hmm_state,"\n")
cat ("CPU:\t",opt$cpu,"\n")
cat ("________________________________________\n\n")

## Create a folder <prefix>_hicnv and perform everything under that ##
folder <- paste0(prefix,"_hicnv")
if (!file.exists(folder)) {
  dir.create(folder)
  setwd(folder)

} else {
  cat (folder," already exists. Writing everything inside it\n")
  setwd(folder)

}

## Find the required programs ##
bedtools <- Sys.which("bedtools")


## Map the coverage values over the refeature file ##
CoverageMapping <- function(refeature, coverage, bedtools) {
  
  sort_cmd <- paste0(bedtools," sort -i ",refeature," > gc_map_frag.bed")
  cat ("Running ",sort_cmd,"\n")
  system(sort_cmd, wait=T)
 
  sort_cmd <- paste0(bedtools," sort -i ",coverage," > coverage.bedGraph")
  cat ("Running ",sort_cmd,"\n")
  system(sort_cmd, wait=T) 

  map_coverage <- paste0(bedtools," map -a gc_map_frag.bed -b coverage.bedGraph -c 5 -o sum -null 0 > coverage.gc_map_frag.bedGraph")
  cat ("Running ",map_coverage,"\n")
  system(map_coverage, wait=T)
}


## 1D regression function implemneted here ##
#ChromWiseRegresssion <- function(i, chrom, bdg, g, m, f) {
ChromWiseRegresssion <- function(bdg, g, m, f, chrom) {
  
  #cat ("Running 1D regression for ",chrom[i],"\n")
  #bdg_keep <- bdg[bdg$chr==chrom[i] & (bdg$gc >= g & bdg$map >= m & bdg$frag >= f),]
  #bdg_bad  <- bdg[bdg$chr==chrom[i] & (bdg$gc < g | bdg$map < m | bdg$frag < f),]

  bdg_keep <- bdg[(bdg$gc >= g & bdg$map >= m & bdg$frag >= f),]
  bdg_bad  <- bdg[(bdg$gc < g | bdg$map < m | bdg$frag < f),]

  count   <- bdg_keep$count
  gcc_vec <- bdg_keep$gc
  map_vec <- bdg_keep$map
  gcc_vec[gcc_vec==0] <- 0.001
  map_vec[map_vec==0] <- 0.001
  len_vec <- log(bdg_keep$frag)
  gcc_vec <- log(gcc_vec)
  map_vec <- log(map_vec)
  len_vec <- (len_vec-mean(c(len_vec)))/sd(c(len_vec))
  gcc_vec <- (gcc_vec-mean(c(gcc_vec)))/sd(c(gcc_vec))
  map_vec <- (map_vec-mean(c(map_vec)))/sd(c(map_vec))

  fit   <- glm(count ~ len_vec+gcc_vec+map_vec,family="poisson")
  coeff <- round(fit$coeff,4)
  res   <- round(count/exp(coeff[1]+coeff[2]*len_vec+coeff[3]*gcc_vec+coeff[4]*map_vec), 2)*100

  bdg_keep[,"norm.count"] <- res 
  bdg_bad[,"norm.count"]  <- 0
  bdg_combined <- as.data.frame(rbind(bdg_keep,bdg_bad))
  bdg_combined_chromwise <- list() 
  for(i in 1:length(chrom)) {
    bdg_combined_chromwise[[i]] <- bdg_combined[bdg_combined$chr == as.character(chrom[i]),] 
    bdg_combined_chromwise[[i]] <- bdg_combined_chromwise[[i]][order(bdg_combined_chromwise[[i]]$start),]
  }
  #bdg_combined <- bdg_combined[order(bdg_combined$start),]
  bdg_combined <- as.data.frame(rbindlist(bdg_combined_chromwise))
  return(bdg_combined)
}

OneDReg <- function(gc_limit, map_limit, frag_limit, prefix) {
  
  if (!file.exists("normalized_data")) {
    dir.create("normalized_data")
  }

  cov_bdg <- as.data.frame(fread("coverage.gc_map_frag.bedGraph",h=F))
  colnames(cov_bdg) <- c("chr","start","end","gc","map","frag","count")
   
  chr_list <- unique(as.character(cov_bdg$chr))
  #print (chr_list)
  #norm_data <- lapply(c(1:length(chr_list)), ChromWiseRegresssion, chr_list, cov_bdg, gc_limit, map_limit, frag_limit)
  #norm_data <- rbindlist(norm_data)
  norm_data <- ChromWiseRegresssion(cov_bdg, gc_limit, map_limit, frag_limit, chr_list) 
  fwrite(norm_data, file=paste0("normalized_data/",prefix,".F_GC_MAP.NormCount.",gc_limit,"_",map_limit,"_",frag_limit,".bedGraph"),row.names=F,quote=F,sep="\t")
  cat ("Wrote ",paste0("normalized_data/",prefix,".F_GC_MAP.NormCount.",gc_limit,"_",map_limit,"_",frag_limit,".bedGraph file"),"\n")

}


## PIC estimate of chromosome wise copy number variation ##
CopyEstimate_func <- function(chr.in, gc_limit, map_limit, frag_limit, ref, drop=TRUE) {

  chr.name  <- chr.in[[1]]
  count     <- chr.in[[2]]
  ref.chrom <- chr.in[[3]]

  segment_density <- c()
  cov_sum <- c()
  chr  <- c()
  len  <- c()
  rpm  <- c()

  i <- 1
  while (i <= length(chr.name)){
    if (is.na(ref)) {
      segment <- suppressWarnings(changepoint::cpt.mean(count[count$chr==chr.name[i],]$norm.count, method="PELT", minseglen=10))
      segment <- length(segment@cpts) + 1
      segment_density <- c(segment_density,segment/nrow(count[count$chr==chr.name[i],]))
    }
    d.chr   <- count[count$chr==chr.name[i] & count$gc >= gc_limit & count$map >= map_limit & count$frag >= frag_limit,]
    cov_sum <- c(cov_sum,sum(d.chr$count))
    len <- length(d.chr$frag)*1000
    rpm <- c(rpm,(sum(d.chr$count)/(len/1e6)))
    chr <- c(chr,as.character(chr.name[i]))
    i <- i+1
  }
  tpm  <- round(rpm/sum(rpm),4)
  data <- data.frame(chr,tpm)
  if (is.na(ref)) {
    data <- data.frame(data, segment_density)
  }
  data[,"ratio"] <- data$tpm/(1/length(chr.name))
  if (!is.na(ref)) {
    #gain <- round(tpm/subset(data,chr==ref.chrom)$tpm,1)
    #gain[gain==0] = 1.0
    #copy_num = 2+(2*(gain-1))
    #data = data.frame(data,gain,copy_num,cov_sum)
    ref.normfactor <- 1/data[data$chr == ref.chrom,]$ratio
    data[,"gain"]  <- data$ratio * ref.normfactor
    data$gain[data$gain==0]  <- 0.25
    data[,"copy_num"] <- round(c(2 * data$gain),1)
    data[,"ref_chrom"] <- ref.chrom
    return(data)

  } else {
    #data[,"ratio"] <- data$tpm/(1/length(chr.name))
    data[,"maha"]  <- mahalanobis(data[,3:4],c(mean(data[,3]),1),cov(data[,3:4]), df=1)
    data <- data[order(data$maha),]
    ref.chrom <- data$chr[1]
    #gain <- round(tpm/subset(data,chr==ref.chrom)$tpm,1)
    #gain[gain==0] = 1.0
    #copy_num = 2+(2*(gain-1))
    #data = data.frame(data,gain,copy_num,cov_sum)
    ref.normfactor <- 1/data[data$chr == ref.chrom,]$ratio
    data[,"gain"]  <- data$ratio * ref.normfactor
    data$gain[data$gain==0]  <- 0.25
    data[,"copy_num"] <- round(c(2 * data$gain),1)
    data[,"ref_chrom"] <- ref.chrom   
    return(data)
  }
}

PICEstimate <- function(gc_limit, map_limit, frag_limit, prefix, ref) {

  cat ("Estimating copy number information for each chromosomes\n")
  norm_bdg <- as.data.frame(fread(paste0("normalized_data/",prefix,".F_GC_MAP.NormCount.",gc_limit,"_",map_limit,"_",frag_limit,".bedGraph"),h=T))
  chr_list <- unique(as.character(norm_bdg$chr))
  
  if (!is.na(ref)) { 
    l <- list()
    l[[1]] <- chr_list
    l[[2]] <- norm_bdg
    l[[3]] <- ref
    ref.chrom <- rep(ref,length(chr_list))
    copyEstimate <- CopyEstimate_func(l, gc_limit, map_limit, frag_limit, ref)
    copyEstimate[,"sample"] <- prefix
    write.table(copyEstimate,file=paste0(prefix,".copyNumber.txt"),row.names=F,quote=F,sep="\t") 
  } else {
    l <- list()
    l[[1]] <- chr_list
    l[[2]] <- norm_bdg
    l[[3]] <- chr_list[1]
    ref.chrom <- rep(ref,length(chr_list))
    copyEstimate <- CopyEstimate_func(l, gc_limit, map_limit, frag_limit, ref)
    copyEstimate[,"sample"] <- prefix
    
    write.table(copyEstimate,file=paste0(prefix,".copyNumber.txt"),row.names=F,quote=F,sep="\t") 
  }
  cat ("Wrote ",paste0(prefix,".copyNumber.txt")," file\n")
}

## Kernal smoothing followed by HMM segmentation is performed ##
Hmm_func <- function(n,chr_obj) {

  h <- depmixS4::depmix(smooth_count~1,family=gaussian(),nstates=n,data=chr_obj)
  result <- tryCatch({
    f <- depmixS4::fit(h)
    b <- BIC(f)      
    return(list(f,b))
  }, warning = function(e) {
    cat ("HMM did not converged for ",n,"\n")
    f <- NA
    b <- NA
  }, error = function(e) {
    cat ("HMM did not converged for ",n,"\n")
    f <- NA
    b <- NA
  })
}

Min_bic <- function(hmm_obj,nstates,chr_obj) {

  d <- list(bic=c(),nstate=c())
  i <- 1
  while (i < nstates) {
    result <- tryCatch({
      d$bic[[i]] <- hmm_obj[[i]][[2]]
      d$nstate[[i]] <- i 
    }, warning = function(e) {
      cat ("No hmm segments for nstate ",i,"\n")
    }, error = function(e) {
      cat ("No hmm segments for nstate ",i,"\n")
    })
    i <- i + 1
  }
  d <- na.omit(as.data.frame(do.call(cbind,d)))
  chr_obj[,"state"] <- hmm_obj[[d[which.min(d$bic),]$nstate]][[1]]@posterior$state
  return(chr_obj)
}

CptsInBedFormat <- function(cpt.obj) {

  start.index <- list()
  end.index   <- list()
  norm.count  <- list()
  start.index[[1]] <- 1
  end.index[[1]]   <- cpt.obj@cpts[1]
  norm.count[[1]]  <- cpt.obj@param.est$mean[1]
  i <- 2
  while (i <= length(cpt.obj@cpts)) {
    start.index[[i]] <- cpt.obj@cpts[i-1]
    end.index[[i]]   <- cpt.obj@cpts[i]
    norm.count[[i]]  <- cpt.obj@param.est$mean[i]
    i <- i + 1
  }
  start.index <- unlist(start.index)
  end.index   <- unlist(end.index)
  norm.count  <- unlist(norm.count)
  return(data.frame(start.index,end.index,norm.count))
}

GetRegionWiseMean <- function(j, seg, bdg) {

  gc    <- mean(bdg[bdg$start.index %in% c(seg$start.index[j]:seg$end.index[j]),]$gc)
  map   <- mean(bdg[bdg$start.index %in% c(seg$start.index[j]:seg$end.index[j]),]$map)
  count <- mean(bdg[bdg$start.index %in% c(seg$start.index[j]:seg$end.index[j]),]$count)
  return(data.frame(gc,map,count))
}

ChromWise_KDE_HMM_Seg <- function(i, chr, bdg, tpm, g, m, f, pfx, ref, hmm_state, win=1e6, seg="yes", y_grid=1.0) {
 
  cat ("Running KDE and HMM segmentation for ",as.vector(chr[i]),"\n")
  bdg_chr <- bdg[bdg$chr == chr[i],]
  bdg_chr$norm.count <- bdg_chr$norm.count * tpm[tpm$chr==chr[i],]$gain
  if (seg == "yes") {
    segment <- suppressWarnings(changepoint::cpt.mean(bdg_chr$norm.count, method="PELT", minseglen=10))
    segment <- CptsInBedFormat(segment)
    bdg_chr[,"start.index"] <- 1:nrow(bdg_chr)
    bdg_chr[,"end.index"]   <- (1:nrow(bdg_chr)) + 1
    segment[,"chr"]   <- chr[i]
    segment[,"start"] <- bdg_chr[bdg_chr$start.index %in% segment$start.index,]$start
    segment[,"end"]   <- bdg_chr[bdg_chr$end.index %in% segment$end.index,]$end
    #feature_mean <- parLapply(cl,c(1:nrow(segment)), GetRegionWiseMean, segment, bdg_chr)
    #feature_mean <- as.data.frame(rbindlist(feature_mean))
    #segment <- data.frame(segment, feature_mean)
    segment[,"frag"]  <- segment$end - segment$start
    #segment <- segment[,c("chr","start.index","end.index","norm.count","start","end","gc","map","frag","count")]
    segment <- segment[,c("chr","start.index","end.index","norm.count","start","end","frag")]
    segment$start.index <- 1:nrow(segment)
    segment$end.index   <- segment$start.index + 1
    bdg_chr <- segment
  } else {
    bdg_chr[,"start.index"] <- 1:nrow(bdg_chr)
    bdg_chr[,"end.index"]   <- (1:nrow(bdg_chr)) + 1
    bdg_chr <- bdg_chr[,c("chr","start.index","end.index","norm.count","start","end","gc","map","frag","count")]
  }

  ## Calculate bandwidth and gridsize here
  frag_mean <- mean(bdg_chr$frag)
  bw <- win/frag_mean
  
  ## Change the y_grid parameter to reduce the grid size and thus reduces time for smoothing
  cat ("Performing Kernel smoothing on ",as.vector(chr[i]),"\n")
  gridsize <- c(nrow(bdg_chr),round(max(bdg_chr$norm.count) * y_grid))
  bandwidth <- c(bw, 2)
  x.dir <- c(min(bdg_chr$start.index),max(bdg_chr$start.index))
  y.dir <- c(min(bdg_chr$norm.count),max(bdg_chr$norm.count))

  ## Kernel smoothing 
  l <- KernSmooth::bkde2D(bdg_chr[,c(2,4)], gridsize=gridsize, bandwidth=bandwidth, range.x = list(x.dir,y.dir), truncate=TRUE)
  pos.max <- apply(l$fhat, 1, which.max)
  y <- l$x2[pos.max]  
 
  ## Correct the edge effect 
  smc.left  <- caTools::runmean(y,k=round(bw),endrule="keep",align="left")
  smc.right <- caTools::runmean(y,k=round(bw),endrule="keep",align="right")
  bdg_chr[,"smooth_count"] <- round((smc.left + smc.right)/2,3)

  ## Segment with HMM 
  cat ("Performing HMM segmenatation on ",as.vector(chr[i]),"\n")
  if (!is.na(cl)) {
    chr_hmm <- parLapply(cl, c(2:hmm_state), Hmm_func, bdg_chr)
  } else {
    chr_hmm <- lapply(c(2:hmm_state), Hmm_func, bdg_chr)
  }
  bdg_chr <- Min_bic(chr_hmm, hmm_state, bdg_chr) 
  hmm_param <- aggregate(smooth_count ~ state, mean, data=bdg_chr)
  hmm_param[,"sd"] <- aggregate(smooth_count ~ state, sd, data=bdg_chr)$smooth_count 

  write.table(bdg_chr, file=paste0("Kernel_Smoothing/",prefix,".", chr[i],".counts.txt"),row.names=F,sep="\t",quote=F)
  write.table(hmm_param, file=paste0("Kernel_Smoothing/",prefix,".", chr[i],".param.txt"),row.names=F,sep="\t",quote=F)
  write.table(l$x1,file=paste0("Kernel_Smoothing/",prefix,".", chr[i],".kde2d_x.txt"),row.names=F,col.names=F,sep="\t",quote=F)
  write.table(l$x2,file=paste0("Kernel_Smoothing/",prefix,".", chr[i],".kde2d_y.txt"),row.names=F,col.names=F,sep="\t",quote=F)
  write.table(l$fhat,file=paste0("Kernel_Smoothing/",prefix,".", chr[i],".kde2d_z.txt"),row.names=F,col.names=F,sep="\t",quote=F)

}

KDE_HMM_Segmentation <- function(gc_limit, map_limit, frag_limit, prefix, ref, genome_width, hmm_state) {

  if (!is.na(cl)) {
    clusterExport(cl, c("Hmm_func","GetRegionWiseMean"))
  }
  if (!file.exists(paste0("Kernel_Smoothing"))) {
    dir.create(paste0("Kernel_Smoothing"))
  }
  norm_bdg   <- as.data.frame(fread(paste0("normalized_data/",prefix,".F_GC_MAP.NormCount.",gc_limit,"_",map_limit,"_",frag_limit,".bedGraph"),h=T))
  sample_tpm <- read.table(paste0(prefix,".copyNumber.txt",sep=""),h=T)
  chr_list   <- unique(as.character(norm_bdg$chr))
  lapply(c(1:length(chr_list)), ChromWise_KDE_HMM_Seg, chr_list, norm_bdg, sample_tpm, gc_limit, map_limit, frag_limit, prefix, ref, hmm_state, genome_width)
}

## CNV labelling against the reference chromosome
CNV <- function(i, chr, m, pfx, cnv_threshold=0.2) {

  cat ("Running for ",as.vector(chr[i]),"\n")
  cnv <- list()
  count_file <- read.table(paste0("Kernel_Smoothing/",pfx,".",chr[i],".counts.txt"), h=T, as.is=T)
  for(i in 1:nrow(count_file)) {
    cnv[[i]] <- round((((count_file$smooth_count[i]/m) - 1) * 2), 3) + 2
  }
  cnv_threshold <- cnv_threshold * 2
  cnv <- unlist(cnv)
  count_file[,"cnv"] <- cnv
  count_file[,"cnv_label"] <- "N"
  count_file[count_file$cnv < (2 - cnv_threshold),"cnv_label"] <- "D"
  count_file[count_file$cnv > (2 + cnv_threshold),"cnv_label"] <- "A"
  return(count_file)
}

CNV_LABEL <- function(prefix, cnv_threshold) {

  if (!file.exists("CNV_Estimation")) {
    dir.create("CNV_Estimation")
  }
  sample <- read.table(paste0(prefix,".copyNumber.txt",sep=""),h=T,as.is=T)
  ref_chrom <- sample$ref_chr[1]
  ref_count <- read.table(paste0("Kernel_Smoothing/",prefix,".",ref_chrom,".counts.txt"), h=T, as.is=T)
  ref_param <- read.table(paste0("Kernel_Smoothing/",prefix,".",ref_chrom,".param.txt"), h=T, as.is=T)
  ref_param$smooth_count <- postmean(x=ref_param$smooth_count, s=ref_param$sd)
  ref_seg   <- as.integer(names(which.max(table(ref_count$state))))
  print (ref_seg)
  ref_mean  <- ref_param[ref_param$state==ref_seg,]$smooth_count
  print (ref_mean)
  print (ref_param)
  cnv <- lapply(c(1:nrow(sample)), CNV, sample$chr, ref_mean, prefix, cnv_threshold)
  cnv <- rbindlist(cnv)
  write.table(cnv, file=paste0("CNV_Estimation/",prefix,".cnv.txt"), row.names=F, sep="\t", quote=F)
  write.table(cnv[,c("chr","start","end","smooth_count","cnv","cnv_label")], file=paste0("CNV_Estimation/",prefix,".cnv.bedGraph"), row.names=F, col.names=F, sep="\t", quote=F)
  
  bedtools <- Sys.which("bedtools")
  if (nchar(bedtools) > 0) {
    cat ("bedtools found! Grouping the cnv.bedGGraph file based on cnv label\n")
    for(i in 1:nrow(sample)) {
      chr <- sample$chr[i]
      groupby_cmd <- paste0("grep -w ",chr," CNV_Estimation/",prefix,".cnv.bedGraph |",bedtools," groupby -g 6 -c 1,2,3,4,5 -o distinct,min,max,mean,mean |awk -v OFS='\t' '{print $2,$3,$4,$5,$6,$1}' |",bedtools," sort > CNV_Estimation/",prefix,".",chr,".cnv.bedGraph")
      print (groupby_cmd)
      system(groupby_cmd, wait=T)
    }
  }
}

## Call all the functions from here ##
## Create the coverage, GC, Mappability and fragment length file
cat ("Calculating coverage values\n")
CoverageMapping(refeature, coverage, bedtools)

## Normalize the counts using 1D regression 
cat ("Performing 1D Regression\n")
OneDReg(gc_cutoff, map_coutoff, frag_cutoff, prefix)

## Estimate an intial whole chromosome copy number change ##
## Estimating whole chromosome copy numbers\n")
cat ("Running PIC calculation\n")
PICEstimate(gc_cutoff, map_coutoff, frag_cutoff, prefix, refchrom)

## Kernal smoothing and HMM segmentation ##
cat ("Performing Kernel smoothing and HMM segmentation\n")
KDE_HMM_Segmentation(gc_cutoff, map_coutoff, frag_cutoff, prefix, refchrom, genome_width, hmm_state)

## CNV labelling
cat ("Performing CVN labelling\n")
CNV_LABEL(prefix, cnv_threshold)
