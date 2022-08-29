# read in n loci and write as phylip 
setwd("~/Documents/bom_phy/")
library(ape)
n <- 20

# delete previous version of the phylip file
file.remove("analyses/bpp_bomarea/bomarea_no_acut.txt")

# get list of files 
files <- list.files("loci/bom_edu_subset/")
files <- files[grep(".nexus", files)]
files <- files[1:n]

# get list of all tips 
tips <- read.nexus(file = "analyses/bpp_bomarea/guide_tree.tre")$tip.label
tips <- tips[grep("acutifolia", tips, invert = T)]
alns <- list()
for (i in 1:length(files)) {
  
  this_file <- paste0("loci/bom_edu_subset/", files[i])
  
  # read the nexus file
  alns[[i]] <- try(read.nexus.data(this_file))
  if (class(alns[[i]]) == "try-error") {
    warning("Could not read alignment ", phased_alignments[i], ". Make sure it is a nexus file")
    next
  }
  
  # drop excluded characters
  aln_lines <- readLines(this_file)
  exclude_line <- grep("EXSET", aln_lines, value = TRUE)
  exclude_line <- strsplit(exclude_line, "=")[[1]][2]
  exclude_line <- strsplit(exclude_line, " ")[[1]]
  exclude_line <- gsub(";", "", exclude_line)
  exclude_line <- exclude_line[exclude_line != ""]
  if ( length(exclude_line) > 0 ) {
    excludes <- unlist(lapply(exclude_line, function(x) {
      tmp <- as.numeric(strsplit(x,"-")[[1]])
      tmp <- seq(head(tmp, 1), tail(tmp, 1))
      return(tmp)
    }))
    
    for(j in 1:length(alns[[i]])) {
      alns[[i]][[j]][excludes] <- "-"
    }
  }
  
  aln_matrix <- matrix(unlist(alns[[i]]),
                       ncol = length(alns[[i]][[1]]),
                       byrow = TRUE)
  rownames(aln_matrix) <- names(alns[[i]])
  aln <- as.DNAbin(aln_matrix, as.matrix = TRUE)
  aln <- del.colgapsonly(aln)
  
  # remove acutifolia
  aln <- aln[grep("acutifolia", rownames(aln), invert = T),]
  
  missing_taxa <- tips[!tips %in% rownames(aln)]
  
  if (length(missing_taxa) > 0 ) {
    missing_rows <- matrix(data = "-",
                           ncol = ncol(aln),
                           nrow = length(missing_taxa))
    rownames(missing_rows) <- missing_taxa
    
    aln_complete <- rbind(aln, as.DNAbin(missing_rows))
  } else {
    aln_complete <- aln
  }
  
  attributes(aln_complete)$dimnames[[1]] <- 
    paste0("^", 
           attributes(aln_complete)$dimnames[[1]])
  
  write.dna(aln_complete, 
            file = "analyses/bpp_bomarea/bomarea_no_acut.txt", 
            format = "sequential",
            append = TRUE)
}

