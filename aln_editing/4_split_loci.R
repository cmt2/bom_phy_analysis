setwd("~/Desktop/bom_phy/cleaned_13April/target/loci_highcoverage_good/")
# First, split  alignments as below: 

rm_num <- function(x) {
  s <- strsplit(x, split = "_")[[1]]
  paste0(s[1:length(s)-1], collapse = "_")
}

# fill in locus name and names of taxa to split below 
locus <- "good_L396.fasta"

names_a <- c(

)

names_b <- c(
  
  
)

# add star symbol back into names for superba sample...sigh
names_a[grep("superba", names_a)] <- gsub("52_A", "52*A", names_a[grep("superba", names_a)])
names_b[grep("superba", names_b)] <- gsub("52_A", "52*A", names_b[grep("superba", names_b)])

# add tilde symbol back into names for edulis colombia sample...sigh
names_a[grep("Colombia_Madri", names_a)] <- gsub("ri_n", "ri~n", names_a[grep("Colombia_Madri", names_a)])
names_b[grep("Colombia_Madri", names_b)] <- gsub("ri_n", "ri~n", names_b[grep("Colombia_Madri", names_b)])

# read in locus
l <- ape::read.FASTA(file = paste0("loci_with_dups/", locus))

# pull out a group
l_a <- l[paste0(names_a, "\t")]
n_a <- unlist(strsplit(names(l_a), "\t"))
names(l_a) <- paste0(unlist(lapply(n_a, rm_num)), "\t")
if (length(names(l_a)) != length(unique(names(l_a)))) warning("more dups!")

# realign 
l_a <- ips::mafft(l_a)
# save 
filename_a <- paste0(unlist(strsplit(locus, split = ".fasta")), 
                   "_a.fasta")
ape::write.FASTA(l_a, file = filename_a)

# pull out b group 
l_b <- l[paste0(names_b, "\t")]
n_b <- unlist(strsplit(names(l_b), "\t"))
names(l_b) <- paste0(unlist(lapply(n_b, rm_num)), "\t")
if (length(names(l_b)) != length(unique(names(l_b)))) warning("more dups!")
# realign 
l_b <- ips::mafft(l_b)
# save 
filename_b <- paste0(unlist(strsplit(locus, split = ".fasta")), 
                     "_b.fasta")
ape::write.FASTA(l_b, file = filename_b)

########## 
setwd("~/Desktop/bom_phy/cleaned_13April/target/loci_highcoverage_good/")
# rename and move good alignments 

rm_num <- function(x) {
  s <- strsplit(x, split = "_")[[1]]
  paste0(s[1:length(s)-1], collapse = "_")
}

files <- list.files()
files <- files[grep(".fasta", files)]

for (i in 1:length(files)) {
  l <- ape::read.FASTA(files[i])
  n_a <- unlist(strsplit(names(l), "\t"))
  names(l) <- paste0(unlist(lapply(n_a, rm_num)), "\t")
  ape::write.FASTA(l, file = paste0("modified_renamed/",files[i]))
  
}
