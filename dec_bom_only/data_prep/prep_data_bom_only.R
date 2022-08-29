setwd("~/Documents/bom_phy/analyses/dec_bom_only/data_prep/")

tree <- RevGadgets::readTrees("relaxed_clock_partition_45_2.0_map.tree")[[1]][[1]]

tips <- tree@phylo$tip.label
tips <- tips[grep("Bomarea",tips)]

distr <- readLines("~/Documents/liliales_paper_big_files/old_version/liliales_checklists/checklists/Alstroemeriaceae.txt")

# start at Bomarea 
distr_bom <- distr[grep("Bomarea Mirb.", distr):length(distr)]
# end at unplaced species
distr_bom <- distr_bom[1:grep("Unplaced Names:", distr_bom)]
# delete empty rows
distr_bom <- distr_bom[distr_bom != ""]
# start after number of species 
distr_bom <- distr_bom[(grep(" Species", distr_bom) + 1 ): 
                         (grep("Synonyms:", distr_bom) - 1)]
# get rid of synonyms 
distr_bom <- distr_bom[grep("===", distr_bom, invert = T, fixed = T)]
distr_bom <- distr_bom[grep("* ", distr_bom, invert = T, fixed = T)]

# get rid of subsp and their corresponding info
subsp <- grep("subsp.", distr_bom, fixed = T)
subsp <- c(subsp, subsp + 1)
distr_bom <- distr_bom[-subsp]

# deal with species 
get_gen_sp <- function(n) {
  s <- strsplit(n, " ")[[1]]
  return(paste0(s[1:2], collapse = " "))
}
species <- distr_bom[seq(from = 1, to = (length(distr_bom) - 1), by = 2)]
species <- unlist(lapply(species, get_gen_sp))

# deal with ranges 
get_ranges <- function(n) {
  r <- strsplit(n, " ")[[1]]
  r <- gsub(".", "", r, fixed = TRUE)
  ranges <- r[grep("^[A-Z]{3}$", r)]
  return(paste0(ranges, collapse = ", "))
}
range <- distr_bom[seq(from = 2, to = length(distr_bom), by = 2)]
range_parsed <- unlist(lapply(range, get_ranges))

distr_refined <- data.frame(species, range_parsed)

get_genus <- function(n) {unlist(strsplit(n, "_"))[[1]]}
get_species <- function(n) {unlist(strsplit(n, "_"))[[2]]}
tips_df <- data.frame(full_name = tips)
tips_df$genus <- unlist(lapply(tips_df$full_name, get_genus))
tips_df$sp <- unlist(lapply(tips_df$full_name, get_species))
tips_df$sp[which(tips_df$sp=="cf")] <- "anceps"
tips_df$species <- paste0(tips_df$genus, " ", tips_df$sp)

combined <- dplyr::left_join(tips_df, distr_refined)
combined$range <- combined$range_parsed
all_ranges <- unique(strsplit(paste0(combined$range_parsed, collapse = ", "), ", ")[[1]])

# recode to broader ranges
# Mexico, CenAm, and Carib is CenAm
# Northern Andes (Colombia, Ecuador, Venezuela) is NAnd
# Central Andes (Peru, Bolivia) is CAnd
# Southern Andes (Western Argentina, Chile) is SAnd
# Eastern S.A. (Brazil, etc) is EaAm

# "BOL" 
combined$range <- gsub("BOL", "CAnd", combined$range)
# "PER" 
combined$range <- gsub("PER", "CAnd", combined$range)
# "CLM" 
combined$range <- gsub("CLM", "NAnd", combined$range)
# "VEN" 
combined$range <- gsub("VEN", "NAnd", combined$range)
# "MXC" 
combined$range <- gsub("MXC", "CenAm", combined$range)
# "MXE" 
combined$range <- gsub("MXE", "CenAm", combined$range)
# "MXG" 
combined$range <- gsub("MXG", "CenAm", combined$range)
# "MXS" 
combined$range <- gsub("MXS", "CenAm", combined$range)
# "MXT" 
combined$range <- gsub("MXT", "CenAm", combined$range)
# "BLZ" 
combined$range <- gsub("BLZ", "CenAm", combined$range)
# "COS"
combined$range <- gsub("COS", "CenAm", combined$range)
# "ELS" 
combined$range <- gsub("ELS", "CenAm", combined$range)
# "GUA" 
combined$range <- gsub("GUA", "CenAm", combined$range)
# "HON" 
combined$range <- gsub("HON", "CenAm", combined$range)
# "NIC"
combined$range <- gsub("NIC", "CenAm", combined$range)
# "PAN" 
combined$range <- gsub("PAN", "CenAm", combined$range)
# "CUB" 
combined$range <- gsub("CUB", "CenAm", combined$range)
# "DOM" 
combined$range <- gsub("DOM", "CenAm", combined$range)
# "HAI" 
combined$range <- gsub("HAI", "CenAm", combined$range)
# "JAM" 
combined$range <- gsub("JAM", "CenAm", combined$range)
# "PUE" 
combined$range <- gsub("PUE", "CenAm", combined$range)
# "TRT" 
combined$range <- gsub("TRT", "CenAm", combined$range)
# "FRG"
combined$range <- gsub("FRG", "NAnd", combined$range)
# "GUY" 
combined$range <- gsub("GUY", "NAnd", combined$range)
# "SUR" 
combined$range <- gsub("SUR", "NAnd", combined$range)
# "BZC" 
combined$range <- gsub("BZC", "EaAm", combined$range)
# "BZE" 
combined$range <- gsub("BZE", "EaAm", combined$range)
# "BZL" 
combined$range <- gsub("BZL", "EaAm", combined$range)
# "BZS" 
combined$range <- gsub("BZS", "EaAm", combined$range)
# "AGE" 
combined$range <- gsub("AGE", "EaAm", combined$range)
# "AGW" 
combined$range <- gsub("AGW", "SAnd", combined$range)
# "PAR" 
combined$range <- gsub("PAR", "EaAm", combined$range)
# "ECU" 
combined$range <- gsub("ECU", "NAnd", combined$range)
# "CLN" 
combined$range <- gsub("CLN", "SAnd", combined$range)
# "CLC"
combined$range <- gsub("CLC", "SAnd", combined$range)
# "CLS"
combined$range <- gsub("CLS", "SAnd", combined$range)

# simpify: 
get_unique_ranges <- function(x) {
  paste0(unique(strsplit(x, ", ")[[1]]), collapse = ", ")
}
combined$range <- unlist(lapply(combined$range, get_unique_ranges))


### check out NAs 
nas <- combined[combined$range == "NA",]

#if sp., it's from Peru so code as NAnd
combined[combined$sp == "sp", "range"] <- "CAnd"
# Kew lists B. distichophylla as a synonym of B. distichifolia, so Andes
combined[combined$sp == "distichophylla", "range"] <- "CAnd, NAnd"
# Bomarea alstroemerioides, spelled wrong in dataset, in C. Andes
combined[combined$sp == "alstroemeriodes", "range"] <- "CAnd"
# Bomarea chimboracensis, spelled wrong in dataset, in N Andes
combined[combined$sp == "chimborazensis", "range"] <- "NAnd"
# Bomarea tribrachiata, spelled wrong in dataset, in N and C Andes
combined[combined$sp == "tribachiata", "range"] <- "CAnd, NAnd"
# Bomarea lehmanii is a synonym of B. straminea according to Kew, 
# which lists it in Andes + Suriname, but according to gbif/ Alzate,
# B. lehmannii and B. straminea are both synonyms of B. pauciflora, 
# which is only in the Andes 
combined[combined$sp == "straminea", "range"] <- "NAnd"
combined[combined$sp == "lehmannii", "range"] <- "NAnd"
combined[combined$sp == "pauciflora", "range"] <- "NAnd"
# Bomarea angustipetala is also a synonym of B. pauciflora according to kew
combined[combined$sp == "angustipetala", "range"] <- "NAnd"
# Bomarea vestita is a synonym of B. multiflora according to kew
combined[combined$sp == "vestita", "range"] <- "NAnd"
# Bomarea foliosa is EITHER a synonym of B. pataconensis or B. multiflora.
# B. foliosa Sodiro is B. pataconensis, B. foliosa Kraenzl. is B. multiflora.
# Either way, Andes! 
combined[combined$sp == "foliosa", "range"] <- "NAnd"
# Bomarea bredemeyerana is spelled different in the checklist, but
# I'm pretty sure my spelling is correct. N Andes 
combined[combined$sp == "bredemeyerana", "range"] <- "NAnd"
# Bomarea acuminata is tricky. Kew lists it as a synonym of B. andreana,
# which is found in CenAm + Andes. But there's one record of B. polyantha
# which has also been synonimized with B. andreana in Venezuela. I'm sticking
# with CenAm + Andes but might need to revisit.
combined[combined$sp == "acuminata", "range"] <- "NAnd, CenAm"
# Bomarea killipii is supposedly a synonym of B. dolichocarpa, but B. dolichocarpa
# is falling with B. edulis in my tree and B. killipi is distinct. B. killipii is 
# found in Peru, so going with Andes
combined[combined$sp == "killipii", "range"] <- "CAnd"
#setacea should be in northern as well as central 
combined[combined$sp == "setacea", "range"] <- "NAnd, CAnd"

# fix B edulis tips 
combined[which(combined$full_name == "Bomarea_edulis_MorelosN_Tribble325"),  "range"] <- "CenAm, CAnd, NAnd, EaAm"
combined[which(combined$full_name == "Bomarea_edulis_Argentina_Novara2378"),  "range"] <- "SAnd"

# fix missID edulis tips
combined[combined$full_name == "Bomarea_edulis_Venezuela_Bunting4817",  "range"] <- "NAnd"
combined[combined$full_name == "Bomarea_edulis_Brazil_Campbell8900",  "range"] <- "EaAm"

### convert to nexus formatted data! 
lines <- character()
lines[1] <- "#NEXUS"
lines[2] <- "Begin data;"
lines[3] <- paste0("Dimensions ntax=", length(combined$species), " nchar=5;")
lines[4] <- "Format datatype=Standard missing=? gap=- labels='01';"
lines[5] <- "Matrix"
lines[6] <- "                         [EANSC: EaAm CAnd NAnd SAnd CenAm]"

for (i in 1:nrow(combined)){
  l <- i + 6
  
  unique_ranges <- unique(unlist(strsplit(combined[i,"range"], split = ", ")))
  binary <- c(0,0,0,0,0)
  if (grepl("EaAm", combined[i,"range"]))    binary[1] <- 1
  if (grepl("CAnd", combined[i,"range"]))   binary[2] <- 1
  if (grepl("NAnd", combined[i,"range"]))    binary[3] <- 1
  if (grepl("SAnd", combined[i,"range"]))      binary[4] <- 1
  if (grepl("CenAm", combined[i,"range"]))     binary[5] <- 1  
 
  
  lines[l] <- paste0("     ", 
                     combined$full_name[i], 
                     "          ",
                     paste0(binary, collapse = "")
                     )
  
}
lines <- append(lines, ";")
lines <- append(lines, "End;")

fileConn<-file("../data/bomarea_only_ranges.nex")
writeLines(lines, fileConn)
close(fileConn)


# write out trimmed tree file too
phylo <- tree@phylo
phylo <- ape::keep.tip(phylo, combined$full_name)
ape::write.tree(phylo, "../data/bom_only_MAP.tre")

# drop tips from subsetted posterior 
trees <- RevGadgets::readTrees("relaxed_clock_partition_45_2.0_subsampled.trees",
                               tree_name = "phylogeny")[[1]]
trees_bom_only <- list()
for (i in 1:length(trees)) {
  print(i)
  trees_bom_only[[i]] <- ape::keep.tip(trees[[i]]@phylo, combined$full_name)
}
ape::write.tree(trees_bom_only, file = "../data/bom_only_posterior.tre")
