# plot ASTRAL tree
library(RevGadgets)
library(ggtree)
library(ggplot2)

setwd("~/Documents/bom_phy/analyses/astral_multilocus2/")
t_astral <- readTrees("output/processed_output/Bomarea_multilocus2_BS10_best.tre")


tree_astral <- t_astral[[1]][[1]]
b_edulis_mrca <- ape::getMRCA(tree_astral@phylo,
                              tip = c("Bomarea_edulis_Panama_Duke13712",
                                      "Bomarea_edulis_Argentina_Hunziker12327"))
boms <- tree_astral@phylo$tip.label[grep("Bomarea",
                                         tree_astral@phylo$tip.label)]
bomarea_mrca <- ape::getMRCA(tree_astral@phylo, tip = boms)

alst <- tree_astral@phylo$tip.label[grep("Alstroemeria",
                                         tree_astral@phylo$tip.label)]
alstroemeria_mrca <- ape::getMRCA(tree_astral@phylo, tip = alst)

luz <- tree_astral@phylo$tip.label[grep("Luzuriaga",
                                        tree_astral@phylo$tip.label)]
dry <- tree_astral@phylo$tip.label[grep("Drymophila",
                                        tree_astral@phylo$tip.label)]
luzdry_mrca <- ape::getMRCA(tree_astral@phylo, tip = c(luz,dry))

tree_astral <- list(list(tree_astral))

# get bs values
bs_values <- tree_astral[[1]][[1]]@phylo$node.label
bs_values <- c(rep(NA, times = length(tree_astral[[1]][[1]]@phylo$tip.label)), bs_values)
bs_values <- as.numeric(bs_values)
bs_values[bs_values < .6] <- NA
# bs_values[bs_values == 0.61] <- 0.60
tree_astral[[1]][[1]]@data$posterior <- bs_values

colors <- colFun(4)

pdf("astral.pdf", height = 7, width = 6)
O <- 0
OT <- 0.5
p1 <- plotTree(tree = tree_astral,
               line_width = .5,
               tip_labels = F,
               tree_layout = "cladogram",
               tip_labels_offset =  0.02) +
  # ggplot2::scale_x_continuous(limits = c( -8,2)) +
  ggplot2::scale_x_continuous(limits = c(0, 27)) +
  ggplot2::scale_y_continuous(limits = c(-1, 122)) +
  #ggplot2::ggtitle("ASTRAL") +
  geom_hilight(node = b_edulis_mrca, fill=colors[1], alpha=.6) +
  ggplot2::annotate(geom = "text",
                    label = "B. edulis",
                    color = colors[1],
                    #x = 1, y = 30, 
                    x = 24, y = 30) +
  geom_cladelabel(node=bomarea_mrca,
                  label = "Bomarea",
                  offset.text = OT,
                  color = colors[2],
                  align=T,
                  barsize = 1,
                  offset = O) +
  geom_cladelabel(node=alstroemeria_mrca,
                  label = "Alstroemeria",
                  offset.text = OT,
                  color = colors[3],
                  align=T,
                  barsize = 1,
                  offset = O) +
  geom_cladelabel(node=luzdry_mrca,
                  label = "Luzuriageae",
                  offset.text = OT,
                  color = colors[4],
                  align=T,
                  barsize = 1,
                  offset = O) +
  geom_nodepoint( aes(fill = posterior, size = posterior), 
                  shape = 21, color = "black") +
  guides(fill = guide_legend("Local posterior"), size=guide_legend("Local posterior")) +
  scale_fill_gradient(low = "white", high = "black",
                      breaks = c(1, .90, .80, .70, .61),
                      labels = c(1, .90, .80, .70, .60)) +
  scale_size_continuous(range = c(.5,4),
                        breaks = c(1, .90, .80, .70, .61),
                        labels = c(1, .90, .80, .70, .60)) +
  theme(legend.position = c(.1, .85)) 

# switch order of layers so blue box is below tree
p2 <- p1 
p2$layers <- p2$layers[c(3,1,2,4:11)]

print(p2)

dev.off()
