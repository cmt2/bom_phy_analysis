library(RevGadgets)

astral_tree <- readTrees("~/Documents/bom_phy/analyses/astral_multilocus2/output/processed_output/Bomarea_multilocus2_BS10_best.tre")
ml_tree <- readTrees("~/Documents/bom_phy/analyses/iqtree_part/cipres_rapid1000/output.tre")

# reroot 
astral_tree <- rerootPhylo(astral_tree, outgroup = c("Luzuriaga_radicans_Layton84",
                                                     "Drymophila_moorei_Australia_Copeland4560"))
ml_tree <- rerootPhylo(ml_tree, c("Luzuriaga_radicans_Layton84", 
                                  "Drymophila_moorei_Australia_Copeland4560"))


# then move over the node label data 

astral_tree <- matchNodeLabels(tree1 = astral_tree)
ml_tree <- matchNodeLabels(tree1 = ml_tree)


trees <- list(astral_tree,
              svdq_tree,
              ml_tree)
names(trees) <- c("astral", "ml")

m = 2
for (t in 1:length(trees)) {
  pdf(paste0("~/Desktop/treeplot_", names(trees)[t], ".pdf"), height = 11*m, width = 8.5*m)
  p <- plotTree(trees[[t]], 
                node_labels = "posterior", 
                node_labels_color = "red")
  print(p)
  dev.off() 
}

