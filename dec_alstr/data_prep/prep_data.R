setwd("~/Documents/bom_phy/analyses/dec/data_prep/")

tree <- RevGadgets::readTrees("relaxed_clock_partition_45_2.0_map.tree")[[1]][[1]]

tips <- tree@phylo$tip.label

distr <- read.csv("~/Documents/liliales_paper_big_files/old_version/liliales_checklists/liliales_traits_edited.csv")

distr <- distr[c(grep("Alstroemeria", distr$species), 
                 grep("Bomarea", distr$species),
                 grep("Luzuriaga", distr$species),
                 grep("Drymophila", distr$species)),
               c("species", "range")]


get_genus <- function(n) {unlist(strsplit(n, "_"))[[1]]}
get_species <- function(n) {unlist(strsplit(n, "_"))[[2]]}
tips_df <- data.frame(full_name = tips)
tips_df$genus <- unlist(lapply(tips_df$full_name, get_genus))
tips_df$sp <- unlist(lapply(tips_df$full_name, get_species))
tips_df$sp[which(tips_df$sp=="cf")] <- "anceps"
tips_df$species <- paste0(tips_df$genus, " ", tips_df$sp)

combined <- dplyr::left_join(tips_df, distr)

# recode to broader ranges

# Australia 
combined$range_names <- gsub("50", "Aust", combined$range)
# Mexico -> CenAm+Carib
combined$range_names <- gsub("79", "CenAm+Carib", combined$range_names)
# Central America -> CenAm+Carib
combined$range_names <- gsub("80", "CenAm+Carib", combined$range_names)
# Caribbean -> CenAm+Carib
combined$range_names <- gsub("81", "CenAm+Carib", combined$range_names)
# Suriname -> EastSA
combined$range_names <- gsub("82", "EastSA", combined$range_names)
# Western Southern America -> Andes
combined$range_names <- gsub("83", "Andes", combined$range_names)
# Suriname, etc. + Brazil
combined$range_names <- gsub("84", "EastSA", combined$range_names)
# Southern South America -> SSA
combined$range_names <- gsub("85", "SSA", combined$range_names)
# Subantarctic Islands -> SSA
combined$range_names <- gsub("90", "SSA", combined$range_names)

### check out NAs 

nas <- combined[is.na(combined$range_names),]

#if sp., it's from Peru so code as Andes
combined[combined$sp == "sp", "range_names"] <- "Andes"
# Kew lists B. distichophylla as a synonym of B. distichofolia, so Andes
combined[combined$sp == "distichophylla", "range_names"] <- "Andes"
# Bomarea alstroemerioides, spelled wrong in dataset, in Andes
combined[combined$sp == "alstroemeriodes", "range_names"] <- "Andes"
# Bomarea chimboracensis, spelled wrong in dataset, in Andes
combined[combined$sp == "chimborazensis", "range_names"] <- "Andes"
# Bomarea tribrachiata, spelled wrong in dataset, in Andes
combined[combined$sp == "tribachiata", "range_names"] <- "Andes"
# Bomarea lehmanii is a synonym of B. straminea according to Kew, 
# which lists it in Andes + Suriname, but according to gbif/ Alzate,
# B. lehmannii and B. straminea are both synonyms of B. pauciflora, 
# which is only in the Andes 
combined[combined$sp == "straminea", "range_names"] <- "Andes"
combined[combined$sp == "lehmannii", "range_names"] <- "Andes"
combined[combined$sp == "pauciflora", "range_names"] <- "Andes"
# Bomarea angustipetala is also a synonym of B. pauciflora according to kew
combined[combined$sp == "angustipetala", "range_names"] <- "Andes"
# Bomarea vestita is a synonym of B. multiflora according to kew
combined[combined$sp == "vestita", "range_names"] <- "Andes"
# Bomarea foliosa is EITHER a synonym of B. pataconensis or B. multiflora.
# B. foliosa Sodiro is B. pataconensis, B. foliosa Kraenzl. is B. multiflora.
# Either way, Andes! 
combined[combined$sp == "foliosa", "range_names"] <- "Andes"
# Bomarea bredemeyerana is spelled different in the checklist, but
# I'm pretty sure my spelling is correct. Andes
combined[combined$sp == "bredemeyerana", "range_names"] <- "Andes"
# Bomarea acuminata is tricky. Kew lists it as a synonym of B. andreana,
# which is found in CenAm + Andes. But there's one record of B. polyantha
# which has also been synonimized with B. andreana in Venezuela. I'm sticking
# with CenAm + Andes but might need to revisit.
combined[combined$sp == "acuminata", "range_names"] <- "Andes CenAm+Carib"
# Bomarea killipii is supposedly a synonym of B. dolichocarpa, but B. dolichocarpa
# is falling with B. edulis in my tree and B. killipi is distinct. B. killipii is 
# found in Peru, so going with Andes
combined[combined$sp == "killipii", "range_names"] <- "Andes"
# Alstroemeria isabelleana spelled wrong in dataset, SSA + EastSA
combined[combined$sp == "isabellana", "range_names"] <- "SSA EastSA"
# B. amilcariana should be Andes rather than East SA 
combined[combined$sp == "amilcariana", "range_names"] <- "Andes"

# fix B edulis tips 
combined[which(combined$full_name == "Bomarea_edulis_MorelosN_Tribble325"),  "range_names"] <- "CenAm+Carib Andes EastSA"
combined[which(combined$full_name == "Bomarea_edulis_Argentina_Novara2378"),  "range_names"] <- "SSA"


# fix missID edulis tips
combined[combined$full_name == "Bomarea_edulis_Venezuela_Bunting4817",  "range_names"] <- "Andes"
combined[combined$full_name == "Bomarea_edulis_Brazil_Campbell8900",  "range_names"] <- "EastSA"

nchar <- length(unique(unlist(strsplit(paste0(combined$range_names, collapse = " "), " "))))
### convert to nexus formatted data! 
lines <- character()
lines[1] <- "#NEXUS"
lines[2] <- "Begin data;"
lines[3] <- paste0("Dimensions ntax=", length(combined$species), " nchar=", nchar, ";")
lines[4] <- "Format datatype=Standard missing=? gap=- labels='01';"
lines[5] <- "Matrix"
lines[6] <- "                         [CENSA: CenAm+Carib,EastSA,Andes,SSA,Aust]"

for (i in 1:nrow(combined)){
  l <- i + 6
  
  unique_ranges <- unique(unlist(strsplit(combined[i,"range_names"], split = " ")))
  binary <- c(0,0,0,0,0)
  if (grepl("CenAm", combined[i,"range_names"]))    binary[1] <- 1
  if (grepl("EastSA", combined[i,"range_names"]))   binary[2] <- 1
  if (grepl("Andes", combined[i,"range_names"]))    binary[3] <- 1
  if (grepl("SSA", combined[i,"range_names"]))      binary[4] <- 1
  if (grepl("Aust", combined[i,"range_names"]))     binary[5] <- 1  
  
  lines[l] <- paste0("     ", 
                     combined$full_name[i], 
                     "          ",
                     paste0(binary, collapse = "")
  )
  
}
lines <- append(lines, ";")
lines <- append(lines, "End;")

fileConn<-file("../data/bomarea_ranges.nex")
writeLines(lines, fileConn)
close(fileConn)

# subset tree posterior 
#trees <- RevGadgets::readTrees("relaxed_clock_partition_45_2.0_combined.trees",
#                               tree_name = "phylogeny")[[1]]

lines <- readLines(con = "relaxed_clock_partition_45_2.0_combined.trees")
trees_subsample <- sample(lines[2:length(lines)], 1000)
lines_subsample <- c(lines[1], trees_subsample)
writeLines(lines_subsample, 
           con = "../data/relaxed_clock_partition_45_2.0_subsampled.trees")
