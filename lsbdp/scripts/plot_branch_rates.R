# get Andean vs. non-Andean rates 

setwd("~/Documents/bom_phy/analyses/")
source("dec_bom_only/scripts/new_tip_labels_formatted.R")
source("dec_bom_only/scripts/makeStateLabels.R")

library(RevGadgets)
library(tidyverse)
library(cowplot)

# load bomarea mapped rates 
rates <- readTrace(paths = c("lsbdp/output/bomarea_only_BDS_tp_rates_1.log", 
                             "lsbdp/output/bomarea_only_BDS_tp_rates_2.log",
                             "lsbdp/output/bomarea_only_BDS_tp_rates_3.log"),
                   burnin = 0)
rates <- combineTraces(rates)
tree <- readTrees("lsbdp/data/bom_only_MAP.tre")
rate_dat <- processBranchData(tree, rates, burnin = 0, net_div = TRUE)

# load rate anc_states 
rates_anc_states <- readTrace(paths = c("lsbdp/output/bomarea_only_BDS_tp_ancstates_1.log", 
                                        "lsbdp/output/bomarea_only_BDS_tp_ancstates_2.log",
                                        "lsbdp/output/bomarea_only_BDS_tp_ancstates_3.log"),
                               burnin = 0)
rates_anc_states <- combineTraces(rates_anc_states)

# load model log for rates
model <- readTrace(paths = c("lsbdp/output/bomarea_only_tp_BDS_1.log", 
                             "lsbdp/output/bomarea_only_tp_BDS_2.log",
                             "lsbdp/output/bomarea_only_tp_BDS_3.log"),
                   burnin = 0)
model <- combineTraces(model)

# load bomarea area reconstructions 

analysis <- "multitrees"
file <- paste0("dec_bom_only/output/dec_bom_only_", analysis, "_combined.ase.tre")
labs <- makeStateLabels("dec_bom_only/data/bomarea_only_ranges.nex", 
                         paste0("dec_bom_only/output/dec_bom_only_", analysis, "_state_labels.txt"))
dec_example  <- processAncStates(file , state_labels = labs)
dec_example@data$node <- as.integer(dec_example@data$node)

geo_states <- readTrace("dec_bom_only/output/dec_bom_only_multitrees_combined.states.log")


# calculate rates for core bomarea and not-core bomarea 
all_tips <- rate_dat[[1]][[1]]@phylo$tip.label
non_core_tips <- c("Bomarea_edulis_Argentina_Novara2378", 
                   "Bomarea_edulis_MorelosN_Tribble325",
                   "Bomarea_obovata_Ecuador_Clark4985", 
                   "Bomarea_ovallei_KEW",
                   "Bomarea_salsilla_Chile_Ackerman545")
core_tips <- all_tips[!all_tips %in% non_core_tips]

core_bom <- RevGadgets::dropTip(rate_dat, non_core_tips)
mean(core_bom[[1]][[1]]@data$net_div)

non_core_bom <- RevGadgets::dropTip(rate_dat, core_tips)
mean(non_core_bom[[1]][[1]]@data$net_div, na.rm = T)


geo_states_letters <- geo_states[[1]]
# recode to area letters:
for (i in 1:nrow(geo_states[[1]])) {
  for (j in 2:ncol(geo_states[[1]])) {
    
    geo_states_letters[i,j] <- labs[ which(names(labs) == geo_states[[1]][i,j]) ]
    
  }
}

# vectors to store final rates
mean_andean <- vector("numeric", length = nrow(geo_states_letters))
mean_non_andean <- vector("numeric", length = nrow(geo_states_letters))

# grab a generation of the biogeo model 
for (i in 1:nrow(geo_states_letters)) {
  
  # grab the end states only 
  end_states <- geo_states_letters[i, grep("end", colnames(geo_states_letters)) ]
  
  # split in andean vs. non-Andean 
  andean <- colnames(end_states)[end_states == "A" | 
                                 end_states == "N" |
                                 end_states == "AN"]
  
  non_andean <- colnames(end_states)[!colnames(end_states) %in% andean]
  
  # choose a random generation from the rates model 
  n <- sample(1:nrow(rates_anc_states[[1]]), 1)

  # get those rate categories for the andean and non-andean nodes 
  andean_rate_cats <- rates_anc_states[[1]][n, andean] + 1 # correct for 0 indexing 
  non_andean_rate_cats <- rates_anc_states[[1]][n, non_andean] + 1 # correct for 0 indexing 
  
  # get the corresponding speciation rates
  net_div_rates <- model[[1]][n, grep("speciation\\[", colnames(model[[1]]))] -
                   model[[1]][n, grep("extinction\\[", colnames(model[[1]]))]
  
  andean_rates <- net_div_rates[as.integer(andean_rate_cats)]
  mean_andean[i] <- mean(as.numeric(andean_rates))
  non_andean_rates <- net_div_rates[as.integer(non_andean_rate_cats)]
  mean_non_andean[i] <- mean(as.numeric(non_andean_rates))

}

# get diff 
diff <- mean_andean - mean_non_andean


# calculate 'bayes factor' 
p <- sum(diff > 0)/ length(diff)
bf <- 2 * log(p / (1 - p))

# plot tree with rates and if Andean 
t_combined <- list(list(dec_example))
t_combined[[1]][[1]]@data <- full_join(t_combined[[1]][[1]]@data, 
                                       rate_dat[[1]][[1]]@data)

t_combined[[1]][[1]]@data$andean <- t_combined[[1]][[1]]@data$end_state_1 == "A" |
                                    t_combined[[1]][[1]]@data$end_state_1 == "N" |
                                    t_combined[[1]][[1]]@data$end_state_1 == "AN"

pdf("lsbdp/figures/andean_rates.pdf", height = 11, width = 8.5)
p <- plotTree(t_combined, 
              color_branch_by = "net_div", 
              branch_color = c("royalblue4", "red3"),
              line_width = 1.5,
              tip_labels = F) +
       ggtree::geom_nodepoint(aes(subset = andean), 
                              color = "white", 
                              shape = 17,
                              size = 3) +
       ggtree::geom_nodepoint(aes(subset = andean), 
                              color = "black", 
                              shape = 2,
                              size = 3) +
       ggtree::geom_nodepoint(aes(subset = !andean), 
                              color = "white",
                              shape = 16,
                              size = 3) +
       ggtree::geom_nodepoint(aes(subset = !andean), 
                              color = "black",
                              shape = 1,
                              size = 3) +
       ggtree::geom_tippoint( aes(subset = andean), 
                              color = "white", 
                              shape = 17,
                              size = 3) +
       ggtree::geom_tippoint( aes(subset = andean), 
                              color = "black", 
                              shape = 2,
                              size = 3) +
       ggtree::geom_tippoint( aes(subset = !andean), 
                              color = "white",
                              shape = 16,
                              size = 3) +
       ggtree::geom_tippoint( aes(subset = !andean), 
                              color = "black",
                              shape = 1,
                              size = 3) +
       guides(color = guide_colorbar(title.position = "left",
                                     #title = "(ELMyrs)")) +
                                     title = "Net-Diversification\nRate (ELMyrs)")) +
       theme(legend.position = c(0.22, 0.35), 
             legend.key.height = unit(30, "native"),
             legend.key.width = unit(30, "native"),
             legend.text = element_text(size = 16),
             legend.title = element_text(size = 16, 
                                         hjust = 0.5,
                                         vjust = 0.5,
                                         angle = 90)) 
p
dev.off()

# also plot the density of the diff to add in 

dens <- density(diff)
df <- data.frame(x = dens$x,
                 y = dens$y,
                 cond = dens$x > 0 )

x_axis <- round(seq(from = min(df$x),
                    to = max(df$x),
                    length.out = 10), 
                digits = 1)

pdf("lsbdp/figures/diff_density.pdf", height = 3, width = 4)
ggplot(df, aes(x, y)) + 
  geom_line() +
  geom_area(data = filter(df, cond == TRUE), fill = 'aquamarine3') +
  geom_line() +
  geom_vline(xintercept = 0, lty = "dashed") +
  xlab("Difference in average net-diversification rate (ELMyrs)") +
  ylab("Density") +
  scale_x_continuous(breaks = x_axis, labels = x_axis) +
  annotate("text", 
           x = .275, y = 4, 
           label = paste0("2lnBF: ", round(bf, digits = 2)), 
           size = 5) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
dev.off()
