setwd("~/Desktop/relaxed_rev_dating/data/genes_selected_15/")

library(ape)

# read in all loci
loci <- list()
files <- list.files()
for (i in 1:length(files)) {
  loci[[i]] <- read.dna(file = files[i], format = "fasta", as.matrix = T)
}

# split into codon positions 
partitioned_loci <- list()
n = 1
for (i in 1:length(loci)) {
  partitioned_loci[[n]] <- loci[[i]][,rep(c(T,F,F), length.out = ncol(loci[[i]]))]
  n = n + 1
  partitioned_loci[[n]] <- loci[[i]][,rep(c(F,T,F), length.out = ncol(loci[[i]]))]
  n = n + 1
  partitioned_loci[[n]] <- loci[[i]][,rep(c(F,F,T), length.out = ncol(loci[[i]]))]
  n = n + 1
}
names <- unlist(strsplit(files, ".fasta"))
names <- paste0(rep(names, each = 3), rep(c("_c1","_c2","_c3"), times = 15))

# save as fasta files

for (i in 1:length(partitioned_loci)) {
  write.FASTA(partitioned_loci[[i]], paste0("../", names[i], ".fasta"))
}