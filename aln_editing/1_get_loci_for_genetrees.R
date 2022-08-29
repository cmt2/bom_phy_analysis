# Identify high-coverage loci for gene tree analysis 

##### set up #####
# packages
library(ape)
library(tidyverse)

# custom functions
get_extr_num <- function(names) {
  n1 <- matrix(unlist(strsplit(names, split = "\\.")), ncol = 3, byrow = T)[,2]
  gsub("_sp", "", x = n1)
}

match_names <- function(n) {
  e <- get_extr_num(n)
  copy <- unlist(strsplit(n, split = "\\."))[3]
  paste0(name_dict$new_name[which(name_dict$extr_num == e)], "_", copy, "\t")
}

# directories
setwd("~/Desktop/Tribble_13April_Target/")
file_path <- "all_loci/"

##### read in data #####
alns_all <- list()
files <- list.files(path = file_path)
files_all <- files[grep("afterMerge.mafft.oneline.fasta", files)]

for (i in 1:length(files_all)) {
  a <- read.dna(file = paste0(file_path, files_all[i]), format = "fasta")
  if (is.null(a) == FALSE) {
    alns_all[[i]] <- as.alignment(a)
    names(alns_all)[i] <- strsplit(files_all[i],split = "\\.")[[1]][[1]]
  } else {
    alns_all[[i]] <- NA
    names(alns_all)[i] <- strsplit(files_all[i],split = "\\.")[[1]][[1]]
  }
}

##### identify high-coverage loci #####
num_seqs <- integer()
for (i in 1:length(alns_all)) {
  if (is.na(alns_all[[i]])) {
    num_seqs[i] <- NA
  } else {
    names <- unique(get_extr_num(alns_all[[i]]$nam))
    num_seqs[i] <- length(names)  
  }
}
# keep loci with at least a third of samples in alignment 
loci_highcoverage <- names(alns_all)[which(num_seqs >= 172/3)]
alns_subset <- alns_all[loci_highcoverage]

##### rename selected loci #####
name_dict <- read.csv("~/Desktop/bom_phy/metadata/sample_info.csv")
reorder <- name_dict$extr_num == "re-order"
name_dict$extr_num[reorder] <- 
  paste0(name_dict$genus[reorder],  sep = "_", name_dict$species[reorder])
name_dict$new_name <- NA
for (i in 1:nrow(name_dict)) {
  if (strsplit(name_dict$extr_num[i], split = "_")[[1]][1] == "CMT" | 
      strsplit(name_dict$extr_num[i], split = "_")[[1]][1] == "KMW") {
    if (is.na(name_dict$region[i])) {
      name_dict$new_name[i] <- paste0(name_dict$genus[i], "_",
                                      name_dict$species[i], "_",
                                      name_dict$collection_num[i])
    } else {
      name_dict$new_name[i] <- paste0(name_dict$genus[i], "_",
                                      name_dict$species[i], "_",
                                      name_dict$region[i], "_",
                                      name_dict$collection_num[i])
    }
    name_dict$new_name[i] <- iconv(name_dict$new_name[i], from = 'UTF-8', to = 'ASCII//TRANSLIT')
    name_dict$new_name[i] <- gsub("\\(", "_", name_dict$new_name[i])
    name_dict$new_name[i] <- gsub("\\)", "", name_dict$new_name[i])
    name_dict$new_name[i] <- gsub("\\.", "_", name_dict$new_name[i])
    name_dict$new_name[i] <- gsub("'", "_", name_dict$new_name[i])
    name_dict$new_name[i] <- gsub("_NA", "", name_dict$new_name[i])
    name_dict$new_name[i] <- gsub("-", "", name_dict$new_name[i])
    name_dict$new_name[i] <- gsub("__", "_", name_dict$new_name[i])
    name_dict$new_name[i] <- gsub(" ", "", name_dict$new_name[i])
    name_dict$new_name[i] <- gsub(" ", "", name_dict$new_name[i])
  } else {
    name_dict$new_name[i] <- paste0(name_dict$extr_num[i],"_KEW")
  }
}

for (i in 1:length(alns_subset)) {
  alns_subset[[i]]$nam <- unlist(lapply(alns_subset[[i]]$nam, match_names))
}

##### move selected loci to new folder #####
if (!dir.exists("loci_highcoverage")) {
  dir.create("loci_highcoverage")
}

for (i in 1:length(alns_subset)) {
  write.FASTA(as.DNAbin(alns_subset[[i]]), 
              file = paste0("loci_highcoverage/", 
                            names(alns_subset)[i], 
                            ".fasta"))
}
