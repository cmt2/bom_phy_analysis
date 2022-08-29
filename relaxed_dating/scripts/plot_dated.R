##### set up #####
library(RevGadgets)
library(treeio)
library(ggplot2)
library(ggtree)
library(deeptime)

setwd("~/Documents/bom_phy/analyses/relaxed_dating/")

##### source a tip labels document #####
source("~/Documents/bom_phy/scripts/plotting/make_tip_labels.R")

##### read in data #####

# dated results
maptree <- readTrees("output/relaxed_clock_partition_45_2.0_map.tree")[[1]][[1]]

# original topology
astraltree <- readTrees("data/rev_starting.tre")[[1]][[1]]

##### get Astral posterior values into dated tree #####

# set root manually
astraltree@phylo$node.label[which(astraltree@phylo$node.label == "")] <- 1 

# get the node labels as numeric posteriors 
astraltree_tib <- as_tibble(astraltree)
astraltree_tib$posterior <- as.numeric(astraltree_tib$label)
astraltree <- as.treedata(astraltree_tib)

# check that the node PP values are correctly mapped onto the topology 
# plotTree(list(list(astraltree)), 
#          node_pp = TRUE, 
#          node_pp_size = 3,
#          node_pp_color = colFun(2))
# they are! 

# remove the false posteriors from maptree 
maptree@data$posterior <- NULL

# add the correct posteriors from ASTRAL, linking with strings of descending tips

astraltree@data$tips <- NA
for (i in 1:nrow(astraltree@data)) {
  
  tips <- geiger::tips(astraltree@phylo, 
                       as.integer(astraltree@data$node[i]))
  tips <- paste0(sort(tips), collapse = "; ")
  astraltree@data$tips[i] <- tips
}

maptree@data$tips <- NA
for (i in 1:nrow(maptree@data)) {
  
  tips <- geiger::tips(maptree@phylo, 
                       as.integer(maptree@data$node[i]))
  tips <- paste0(sort(tips), collapse = "; ")
  maptree@data$tips[i] <- tips
}

maptree@data <- dplyr::full_join(maptree@data, astraltree@data[ ,c("tips","posterior")])
maptree <- list(list(maptree))

##### identify nodes for node symbols #####
# bom_mrca <- ape::getMRCA(maptree[[1]][[1]]@phylo, tip = c("Bomarea_ovallei_KEW", "Bomarea_sp__oso_Peru_Graham12613"))
# alstoid_mrca <- ape::getMRCA(maptree[[1]][[1]]@phylo, tip = c("Alstroemeria_ligtu_Ornduff9112", "Bomarea_sp__oso_Peru_Graham12613"))
# als_mrca <- ape::getMRCA(maptree[[1]][[1]]@phylo, tip = c("Alstroemeria_ligtu_Ornduff9112", "Alstroemeria_apertiflora_Hatschbach17552"))
# luz_mrca <- ape::getMRCA(maptree[[1]][[1]]@phylo, tip = c("Luzuriaga_marginata_Bonifacino550", "Drymophila_cyanocarpa_Australia_Messina989"))
# all_mrca <- ape::getMRCA(maptree[[1]][[1]]@phylo, tip = c("Luzuriaga_marginata_Bonifacino550", "Bomarea_sp__oso_Peru_Graham12613"))
# mrcas <- c(bom_mrca, alstoid_mrca, als_mrca, luz_mrca, all_mrca)
# 
# maptree[[1]][[1]]@data$target_node <- maptree[[1]][[1]]@data$node %in% mrcas


##### fix tip labels #####
maptree[[1]][[1]]@phylo$tip.label <- make_tip_labels(maptree[[1]][[1]])

##### plot #####

pdf("figures/bomChrono.pdf", height = 11, width = 8.5)
plotTree(maptree, 
         timeline = TRUE, 
         geo_units = "epochs",
         node_age_bars = TRUE, 
         tip_labels_formatted = TRUE,
         age_bars_colored_by = "posterior") + 
  coord_geo(dat = list("epochs", "eras"), 
            size = list(4, 4),
            pos = list("bottom", "bottom"),
            xlim = c(-115, 50),
            ylim = c(-4, 92),
            height = grid::unit(2, "line"),
            abbrv = TRUE,
            rot = 90,
            center_end_labels = TRUE,
            bord = c("right", "top", "bottom"),
            neg  = TRUE) +
  guides(colour = guide_colorbar(title.position = "top",
                                 title.vjust = 0,
                                 title.hjust = 0.5,
                                 barwidth = 8,
                                 barheight = 1.5)) +
  # geom_point2(aes(subset = target_node), 
  #             color = "firebrick3",
  #             pch = 8,
  #             size = 5) +
  ggplot2::scale_color_gradient(
    low = RevGadgets::colFun(2)[2],
    high = RevGadgets::colFun(2)[1],
    name = "Posterior",
    breaks = c(.21, .4, .6, .8, 1),
    labels = c(.2, .4, .6, 8., 1)
  ) +
  theme(legend.position = c(.85, -0.03), #
        legend.direction = "horizontal", #
        legend.background = element_blank(),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 14))
dev.off()

