require(HiTC)
require(rtracklayer)
require(BSgenome.Hsapiens.UCSC.hg38)

human_chr <- seqlevels(BSgenome.Hsapiens.UCSC.hg38)[1:24]
resFrag <- getRestrictionFragmentsPerChromosome(resSite="AGATCT", chromosomes=human_chr, overhangs5=1, genomePack="BSgenome.Hsapiens.UCSC.hg38")
allRF <- do.call("c",resFrag)
names(allRF) <- unlist(sapply(resFrag, function(x){paste0("HIC_", seqlevels(x), "_", 1:length(x))}))
export(allRF, format="bed", con="Bg1II_resfrag_hg38.bed")
