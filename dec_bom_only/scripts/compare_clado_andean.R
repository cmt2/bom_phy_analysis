setwd("~/Documents/bom_phy/analyses/")
source("dec_bom_only/scripts/new_tip_labels_formatted.R")
source("dec_bom_only/scripts/makeStateLabels.R")

library(RevGadgets)
library(tidyverse)

setwd("~/Documents/bom_phy/analyses/dec_bom_only/")

analysis <- "multitrees"
file <- paste0("output/dec_bom_only_", analysis, "_combined.ase.tre")
labs <- makeStateLabels("data/bomarea_only_ranges.nex", 
                        paste0("output/dec_bom_only_", analysis, "_state_labels.txt"))
dec_example  <- processAncStates(file , state_labels = labs)

dat <- dec_example@data[, c("end_state_1", "start_state_1", "index", "node")]


# look a node -- check its start state and that node's ancester's end states. 
# if diff, then clado 

nodes <- data.frame(ancestor = as.character(dec_example@phylo$edge[,1]), 
                    node = as.character(dec_example@phylo$edge[,2]))

test <- full_join(dat, nodes)
type <- character() 
for (i in 1:nrow(test)) {
  node <- test[i,"node"]
  node_start_state <- test[i, "start_state_1"]
  ancestor <- test[i, "ancestor"]
  ancestor_end_state <- test[which(test$node == as.integer(ancestor)), "end_state_1"]
  
  if (!is.na(ancestor)) {
    if (node_start_state == ancestor_end_state) {
      if (grepl("N", ancestor_end_state) | grepl("A", ancestor_end_state)) {
        type[i] <- "not clado andean"
      } else { type[i] <- "not clado not andean" }
    } else {
      if (grepl("N", ancestor_end_state) | grepl("A", ancestor_end_state)) {
        type[i] <- "clado andean"
      } else { type[i] <- "clado not andean" }
    }
  }

}

res <- table(type)

res_df <- data.frame(andean = c(19, 131), 
                     non_andean = c(0, 2),
                     row.names = c("clado", "non_clado"))
fisher.test(res_df, simulate.p.value = T)
