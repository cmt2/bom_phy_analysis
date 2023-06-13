# remove bad seqs from alignments
# setwd("my/dir")

# dir for new alignments
dir.create("loci_highcoverage_good")

# read in bad sequences
bad_seqs <- read.table(file = "seqs_to_remove.txt",
                       sep = "\t")
colnames(bad_seqs) <- c("locus", "seq")

# list of loci
loci <- list.files("loci_highcoverage")
loci <- loci[grep(".fasta",loci)]
loci <- loci[-grep(".bionj", loci)]
loci <- loci[-grep("gz", loci)]
loci <- loci[-grep(".contree", loci)]
loci <- loci[-grep(".iqtree", loci)]
loci <- loci[-grep(".log", loci)]
loci <- loci[-grep(".mldist", loci)]
loci <- loci[-grep(".splits.nex", loci)]
loci <- loci[-grep(".uniqueseq.phy", loci)]

# output information on bad seqs 
#bad_info <- list()

# for each locus
for (i in 41:length(loci)) {
  
  # read in alignment
  aln <- ape::read.FASTA(paste0("loci_highcoverage/", loci[i]))
  
  # pull out bad sequences
  bds <- bad_seqs[which(bad_seqs$locus == 
                          unlist(strsplit(loci[i], split = "\\."))[1][1]),2]
  bds <- paste0(bds, "\t")
  #bad_aln <- aln[attributes(aln)$names %in% bds]
  
  # create and save alignment without bads 
  good_aln <- aln[!attributes(aln)$names %in% bds]
  ape::write.FASTA(good_aln, 
                   file = paste0("loci_highcoverage_good/", 
                                 "good_",
                                 loci[i]))

  ## write bad info out to file 
  #bad_info[[i]] <- data.frame(matrix(nrow = length(bds),
  #                                   ncol = 3))
  #bad_info[[i]][,1] <- unlist(strsplit(loci[i], ".fasta"))
  #bad_info[[i]][,2] <- unlist(strsplit(bds, "\t"))
  # for (j in 1:length(bds)) {
  #  bad_info
  #}
  
}
  
  