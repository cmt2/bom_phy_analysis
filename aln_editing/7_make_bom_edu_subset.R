setwd("~/Documents/bom_phy/")
library(ape)
library(ips)
library(dplyr)

# get list of alignments with no dups  
files <- list.files("loci/cleaned_13April/target/loci_highcoverage_good/")
files <- files[grep(".fasta", files)]
  
aln_length <- integer()
num_seqs <- integer()
pis_abs <- integer()
pis_frac <- integer()

for (i in 1:length(files)) {
  # read in the final version of each good alignment alignment
  aln <- read.dna(paste0("loci/cleaned_13April/target/loci_highcoverage_good/modified_renamed/edited_alns/",
                         files[i]), format = "fasta", as.matrix = TRUE)
  
  # trim to just edulis
  aln <- aln[grep("edulis", rownames(aln)),]
  
  # get metrics
  aln_length[i] <- dim(aln)[2]
  num_seqs[i]   <- dim(aln)[1]
  pis_abs[i]    <- pis(aln, what = "absolute")
  pis_frac[i]   <- pis(aln, what = "frac")
}

# combine into data frame 
data <- data.frame(aln_length = aln_length,
                   num_seqs   = num_seqs,
                   pis_abs    = pis_abs,
                   pis_frac   = pis_frac,
                   file       = files)

# get each files rank for each metric 
data$rank_aln_length <- rank(data$aln_length, ties.method = "average")
data$rank_num_seqs   <- rank(data$num_seqs,   ties.method = "average")
data$rank_pis_abs    <- rank(data$pis_abs,    ties.method = "average")
data$rank_pis_frac   <- rank(data$pis_frac,   ties.method = "average")

# get the average rank across metrics 
data$ave_rank <- rowMeans(data[,c("rank_aln_length",
                                   "rank_num_seqs",
                                   "rank_pis_abs",
                                   "rank_pis_frac")])

# select the 20 files with highest average rank 
data <- data[order(data$ave_rank, decreasing = TRUE),]
good_files <- head(data$file, n = 20)
good_files <- gsub("good_", "", good_files)
good_files_base <- gsub(".fasta", "", good_files)
good_new_files <- paste0(good_files_base, ".afterMerge.mafft.oneline.fasta")

# create a new directory with those files 
dir.create("loci/bom_edu_subset")

# for each alignment, we want to only keep the tips that were present in the 
# original alignments

# first, create the name dictionary
name_dict <- read.csv("~/Documents/bom_phy/loci/original/metadata/sample_info.csv")
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

#  some custom function
rm_num <- function(x) {
  s <- strsplit(x, split = "\\.")[[1]]
  paste0(s[1:length(s)-1], collapse = ".")
}
get_extr_num <- function(names) {
  n1 <- matrix(unlist(strsplit(names, split = "\\.")), ncol = 2, byrow = T)[,2]
  gsub("_sp", "", x = n1)
}

match_names <- function(n) {
  e <- get_extr_num(n)
  name_dict$new_name[which(name_dict$extr_num == e)]
}

edulis_tips <- c(
  "Bomarea_edulis_Brazil_CaxiasdoSul_Kegler89",
  "Bomarea_dolichocarpa_Peru_Barbour5069",
  "Bomarea_edulis_Queretaro_Tribble701",
  "Bomarea_edulis_Haiti_Proctor10821",
  "Bomarea_edulis_MorelosN_Tribble325",
  "Bomarea_edulis_Suriname_Lindeman455",
  "Bomarea_edulis_Venezuela_AymardC_1314",
  "Bomarea_edulis_DR_Jimenez3903",
  "Bomarea_edulis_Veracruz_Tribble552",
  "Bomarea_edulis_MorelosS_Tribble377",
  "Bomarea_edulis_CR_J_F_Morales6635",
  "Bomarea_edulis_Nicaragua_Stevens21839",
  "Bomarea_edulis_Panama_Duke13712",
  "Bomarea_edulis_Argentina_M_ulguradeRomero1876",
  "Bomarea_edulis_Bolivia_Solomon7616",
  "Bomarea_edulis_Brazil_MatoGrosso_Eiten9877",
  "Bomarea_edulis_Argentina_Deginani1587",
  "Bomarea_edulis_Brazil_Bahia_Thomas13631",
  "Bomarea_edulis_Honduras_Harmon3784",
  "Bomarea_edulis_Cuba_Smith3264",
  "Bomarea_edulis_Guatemala_Breckon2075",
  "Bomarea_edulis_Chiapas_Tribble90",
  "Bomarea_edulis_Nayarit_Tribble113",
  "Bomarea_edulis_Oaxaca_Tribble76",
  "Bomarea_acutifolia_Mexico_Tribble79",
  "Bomarea_edulis_Veracruz_Tribble65",
  "Bomarea_edulis_Argentina_Hunziker12327",
  "Bomarea_edulis_Argentina_Novara2378"
)

# For each alignment: 
for (i in 1:length(good_new_files)) {
  
  # read it in 
  aln <- read.dna(file = paste0("loci/cleaned_13April/full/", good_new_files[i]),
                  format = "fasta", 
                  as.matrix = TRUE)
  
  # remove number from seqs 
  rownames(aln) <- unlist(lapply(rownames(aln), rm_num))
  
  # rename 
  rownames(aln) <- unlist(lapply(rownames(aln), match_names))
  
  # only keep edulis 
  aln <- aln[rownames(aln) %in% edulis_tips,]
  
  # realign with mafft
  aln <- ips::mafft(x = aln, method = "localpair", maxiterate = 1000)
  
  n <- strsplit(good_new_files[i], "\\.")[[1]][1]
  write.FASTA(aln, 
              file = paste0("loci/bom_edu_subset/", n, ".fasta"))
  
  # check for dups 
  if (any(duplicated(rownames(new_aln)))) stop("duplicated seqs")
}

