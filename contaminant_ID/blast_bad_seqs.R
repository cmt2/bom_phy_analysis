# blast bad seqs and save info on top result
setwd("~/Desktop/carrie")

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
bad_info <- list()
Sys.setenv(BLASTDB = "/Users/zenillab/Documents/nt_db_new/")

# for each locus
for (i in 47:length(loci)) {
  
  # read in alignment
  aln <- ape::read.FASTA(paste0("loci_highcoverage/", loci[i]))
  
  # pull out bad sequences
  bds <- bad_seqs[which(bad_seqs$locus == 
                          unlist(strsplit(loci[i], split = "\\."))[1][1]),2]
  if (length(bds) > 0) {
    bds <- paste0(bds, "\t")
    bad_aln <- aln[attributes(aln)$names %in% bds]
    ape::write.FASTA(ape::del.gaps(bad_aln), file = "seqs.fasta")
    
    ## write bad info 
    bad_info[[i]] <- data.frame(matrix(nrow 
                                       = length(bds),
                                       ncol = 2))
    bad_info[[i]][,1] <- unlist(strsplit(loci[i], ".fasta"))
    bad_info[[i]][,2] <- unlist(strsplit(names(bad_aln), "\t"))
    colnames(bad_info[[i]]) <- c("locus","query_name")
    
    # blast! 
    system( "blastn -db nt -query seqs.fasta -num_threads 8 -max_target_seqs 1 -outfmt '6 evalue sscinames qseqid pident' -out results.out ")
    
    # read in results
    t <- try(read.table( "results.out",  sep = "\t"), silent = T)
    if (class(t) == "try-error") {
      
      bad_info[[i]]$evalue <- NA
      bad_info[[i]]$res_name <- NA
      bad_info[[i]]$perc_id <- NA
      
    } else {
    
      results <- read.table( "results.out",  sep = "\t")
      colnames(results) <- c("evalue","res_name","query_name","perc_id")
    
      # add to bad_info
      bad_info[[i]] <- dplyr::full_join(bad_info[[i]], results) 
    }
    
    print(paste0("finished processing locus: ", loci[i], ", locus number ", i))
    
  } else {
    bad_info[[i]] <- NA
  }
}

#combined output and remove NA rows
bad_info_df <- do.call(rbind, bad_info)
bad_info_df <- bad_info_df[complete.cases(bad_info_df[,1:2]),]

# remove the outgroup contaminants 
# macusani, sclerophylla, and huanuco

mac <- grep("macusani", bad_info_df$query_name)
scler <- grep("sclerophylla", bad_info_df$query_name)
huan <- grep("huanuco", bad_info_df$query_name)

write.csv(bad_info_df, file = "contamination_blast.csv", row.names = F)
