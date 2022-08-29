t <- "good_L85_b_rooted.tre"
t_new <- gsub("_rooted",".fasta", t)
tre <- read.tree(paste0("~/Desktop/FePhyFoFum/trees/", t_new))

outgroup <- c("Luzuriaga_marginata_Bonifacino550",
              "Luzuriaga_polyphylla_CultivatedinCA_UCBG90_2401",
              "Luzuriaga_radicans_Layton84",
              "Drymophila_moorei_Australia_Copeland4560",       
              "Drymophila_cyanocarpa_Australia_Messina989")

has_outgroup <- any(outgroup %in% tre$tip.label) == T
has_alstr <- any(grepl("Alstroemeria",tre$tip.label))

if (has_outgroup) {
  outgroup_custom <- outgroup[outgroup %in% tre$tip.label]
  num <- ape::getMRCA(tre, outgroup_custom)
  midpoint <- 0.5 * tre$edge.length[which(tre$edge[, 2] == num)]
  t_rooted <- phytools::reroot(tre, num, position = midpoint)
} else if (has_alstr) {
  alstr <- tre$tip.label[grep("Alstroemeria",tre$tip.label)]
  num <- ape::getMRCA(tre, alstr)
  midpoint <- 0.5 * tre$edge.length[which(tre$edge[, 2] == num)]
  t_rooted <- phytools::reroot(tre, num, position = midpoint)
} else {
  "root manually"
  tip <- "Bomarea_obovata_Ecuador_Clark4985"
  num <- which(tre$tip.label == tip)
  midpoint <- 0.5 * tre$edge.length[which(tre$edge[, 2] == num)]
  t_rooted <- phytools::reroot(tre, num, position = midpoint)
}


if ("t_rooted" %in% ls()) {
  write.tree(t_rooted, file = paste0("~/Desktop/genesortR_bom/rooted_gt/", t))
}

