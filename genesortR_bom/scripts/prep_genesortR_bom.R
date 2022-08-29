# concatenate high coverage loci for genesortR 
setwd("~/Desktop/genesortR_bom/")
library(ape)

### custom functions
gap_perc <- function(x){
  vec <- strsplit(x, split = "")[[1]]
  return(length(which(vec == "-"))/ length(vec))
}

`%notin%` <- Negate(`%in%`)

loci <- list.files("alns_all/")

loci_to_concat <- list() 
for (i in 1:length(loci)) {
  loci_to_concat[[i]] <- read.dna(paste0("alns_all/",loci[i]), format = "fasta", as.matrix = T)
}

### subset by number of sequences per locus

num_seqs <- unlist(lapply(loci_to_concat, nrow))

# only include genes with > 75% of accessions
loci_to_concat_subset <- loci_to_concat[num_seqs > 90]

### create partitioning file 
lengths_of_loci <- unlist(lapply(loci_to_concat_subset, ncol))
names_of_loci <- unlist(strsplit(loci[num_seqs > 90], ".fasta"))
names_of_loci <- gsub("good_", "", names_of_loci)
lines <- character() 
for (i in 1:length(lengths_of_loci)) {
  if (i == 1) {
    first <- 1
  } else {
    first <- second + 1
  }
  second <- first + lengths_of_loci[i] -1
  lines <- append(lines, paste0(names_of_loci[i],  " = ", first, "-", second))
}

fileConn <- file("~/Desktop/genesortR_bom/partitions.txt")
writeLines(lines, fileConn)
close(fileConn)

### concatenate loci
names(loci_to_concat_subset) <- NULL
loci_to_concat_subset[[length(lengths_of_loci) + 1]] <- TRUE
names(loci_to_concat_subset)[length(lengths_of_loci) + 1] <- "check.names"
loci_to_concat_subset[[length(lengths_of_loci) + 2]] <- TRUE
names(loci_to_concat_subset)[length(lengths_of_loci) + 2] <- "fill.with.gaps"
loci_to_concat_subset <- loci_to_concat_subset[!is.na(loci_to_concat_subset)]
concat <- do.call(cbind.DNAbin, loci_to_concat_subset)

write.FASTA(concat, file = "concat.fasta")

### get corresponding gene trees into one newick string 
names_of_gt <- paste0("good_", names_of_loci, "_rooted.tre")

gt <- list()
for (i in 1:length(names_of_gt)) {
  gt[[i]] <- ape::read.tree(paste0("rooted_gt/", names_of_gt[i]))
}

class(gt) <- "multiPhylo"

write.tree(gt, file = "all_gt.tre")
