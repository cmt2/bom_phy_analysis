setwd("~/Documents/bom_phy/loci/cleaned_13April/target/loci_highcoverage_good/modified_renamed/")

t_iq <- RevGadgets::readTrees("~/Documents/bom_phy/analyses/iqtree_part/cipres_rapid1000/output.tre")
tips_to_keep <- t_iq[[1]][[1]]@phylo$tip.label

files <- list.files()

get_fastas <- function(f) {return(length(unlist(strsplit(f, split = "\\."))) == 2)}
r <- unlist(lapply(files, get_fastas))
files_fasta <- files[r]

for (i in 1:length(files_fasta)) {
  aln <- ape::read.dna(paste0(files_fasta[i]), format="fasta", as.matrix=TRUE)
  rownames(aln) <- gsub("\t","", rownames(aln), fixed = T)
  rownames(aln) <- gsub("*", "_", rownames(aln), fixed = T)
  aln_subbed <- aln[rownames(aln) %in% tips_to_keep,]
  aln_subbed <- ips::mafft(aln_subbed)
  ape::write.FASTA(aln_subbed, file = paste0("edited_alns/", files_fasta[i]))
}
