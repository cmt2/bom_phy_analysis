# concatenate loci
setwd("~/Desktop/bom_phy/cleaned_13April/target/loci_highcoverage_good/modified_renamed/")

library(ape)

gap_perc <- function(x){
  vec <- strsplit(x, split = "")[[1]]
  return(length(which(vec == "-"))/ length(vec))
}

`%notin%` <- Negate(`%in%`)


loci <- list.files()
#loci <- loci[101:150]
loci_to_concat <- list() 
for (i in 1:length(loci)) {
  loci_to_concat[[i]] <- read.FASTA(loci[i])
  loci_to_concat[[i]] <- ips::mafft(loci_to_concat[[i]])
  
  # remove symbol if present
  rownames(loci_to_concat[[i]]) <- gsub("\t","", rownames(loci_to_concat[[i]]), fixed = T)
}

###### write out partitioning file for partition finder ######
lengths_of_loci <- unlist(lapply(loci_to_concat, ncol))
names_of_loci <- unlist(strsplit(loci, ".fasta"))
names_of_loci <- gsub("good_", "", names_of_loci)
lines <- character() 
for (i in 1:length(lengths_of_loci)) {
  if (i == 1) {
    first <- 1
  } else {
      first <- second + 1
    }
  second <- first + lengths_of_loci[i] -1

  lines <- append(lines, paste0(names_of_loci[i], "_codon1", " = ", first, "-", second,"\\3;"))
  lines <- append(lines, paste0(names_of_loci[i], "_codon2", " = ", first+1, "-", second,"\\3;"))
  lines <- append(lines, paste0(names_of_loci[i], "_codon3", " = ", first+2, "-", second,"\\3;"))
}

fileConn <- file("~/Desktop/partitions.txt")
writeLines(lines, fileConn)
close(fileConn)

###### concatenate good loci ######

names(loci_to_concat) <- NULL
loci_to_concat[[length(loci) + 1]] <- TRUE
names(loci_to_concat)[length(loci) + 1] <- "check.names"
loci_to_concat[[length(loci) + 2]] <- TRUE
names(loci_to_concat)[length(loci) + 2] <- "fill.with.gaps"
loci_to_concat <- loci_to_concat[!is.na(loci_to_concat)]
concat <- do.call(cbind.DNAbin, loci_to_concat)


###### remove accessions with > .90% gaps ######

concat_trimmed <- as.alignment(concat)
r <- unlist(lapply(concat_trimmed$seq, gap_perc))
taxa_keep <- which(r < 0.90)

###### also remove contaminant accessions ######

contams <- c(grep("macusani", concat_trimmed$nam),
             grep("huanuco", concat_trimmed$nam),
             grep("sclerophylla", concat_trimmed$nam))
taxa_keep <- taxa_keep[taxa_keep %notin% contams]

concat_trimmed$seq <- concat_trimmed$seq[taxa_keep]
concat_trimmed$nam <- concat_trimmed$nam[taxa_keep]
concat_trimmed$nb <- length(taxa_keep)

###### concatenate accessions ######
concat_trimmed_DNA <- as.DNAbin(concat_trimmed)
write.dna(concat_trimmed_DNA, file = "concat.phy",format = "sequential")
write.FASTA(concat_trimmed_DNA, file = "concat.fasta")

