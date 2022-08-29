library(RevGadgets)
library(ggplot2)
setwd("~/Documents/bom_phy/analyses/dec_bom_only/")

# analysis <- "single"
analysis <- "multitrees"
run <- "combined"

## plot tree 
file <- paste0("output/dec_bom_only_", analysis,"_", run, ".ase.tre")

## source new tip labels
source("~/Documents/bom_phy/scripts/plotting/make_tip_labels.R")

source("scripts/makeStateLabels.R")
# labels that correspond to each region/ possible combination of regions
labs <- makeStateLabels("data/bomarea_only_ranges.nex", 
                        paste0("output/dec_bom_only_", analysis, "_state_labels.txt"))

# labs[1:5] <- c("E: Eastern S. America", "A: Central Andes", "N: Northern Andes", 
#                "S: Southern Andes", "C: Central America")
dec_example  <- processAncStates(file , state_labels = labs)

ncol <- length(dec_example@state_labels)


# A (central Andes), C (Cen Am), E (Eastern SA), N (northern andes), S (southern andes)
colors_main <- c("#D81B60", "#1E88E5", "#004D40", "#8214a0", "#FFC107")
colors_combined <- colorRampPalette(c("#94C08B","#fa78fa","#fae6be", "darkblue", "#a0fa82", "red4","darkorange", "turquoise3",
                                    "slategray1", "thistle1"))(length(dec_example@state_labels)-5)
colors <- c(colors_main, colors_combined)           

state_labels <- dec_example@state_labels
state_labels <- state_labels[c(1, 8, 9, 15, 17,
                               2, 3, 6, 10, 14, 16,
                               4, 5, 7, 11,
                               12, 13)]
names(colors) <- state_labels



# Andean diversification scale: 
andean_timescale <- data.frame(name =c("Final", "CAN ends", 
                                       "NAN east", "S. marine inc.", 
                                       "CAN east", "W. Patagonia ends", 
                                       "NAN starts", "Marginal arc"),
                               abbr = c("Final", "CAN ends", 
                                        "NAN east", "S. marine inc.", 
                                        "CAN east", "W. Patagonia ends", 
                                        "NAN starts", "Marginal arc"),
                               max_age = c(5,10,15,23,34,55,66,80),
                               min_age = c(0,5,10,15,23,34,55,66),
                               color = rev(terrain.colors(8)))

timeseries <- data.frame(name = c(""), 
                         max_age = 10, 
                         min_age = 0, 
                         color = "#FFFFFF00")

# fix tip labels 
dec_example@phylo$tip.label <- make_tip_labels(dec_example)

pdf(paste0("figures/plot_pies_", analysis, "_combined.pdf"),height = 11, width = 8.5)
p <- plotAncStatesPie(t = dec_example, 
                      pie_colors = colors, 
                      tip_pie_nudge_x = 0.45, 
                      tip_labels_states_offset = 0.2,
                      tip_labels_size = 3, 
                      tip_labels_formatted = TRUE,
                      node_pie_size = 0.5,
                      tip_pie_size = 0.4, 
                      tip_labels_states = TRUE,
                      cladogenetic = TRUE, 
                      tip_labels_offset = 0.55, 
                      state_transparency = 1,
                      geo_units = list("eras", "epochs"),
                      time_bars = TRUE,
                      timeline = T) +
  
  # change state labels in legend 
  scale_color_discrete(labels = c ("A: Central Andes",
                                   "C: Central America",
                                   "E: Eastern S. America",
                                   "N: Northern Andes",
                                   "S: Southern Andes",
                                   "AC", "AN", "AS", "EA",                
                                   "EN", "NC", "ANC", "ANS", 
                                   "ASC", "EAN", "EANC", "EANS", "other"),
                       breaks = c(names(colors), "other"),
                       type = c(colors, "grey50")) +

  # isthmus 
  ggplot2::geom_vline(xintercept = -2.8, linetype = "longdash") +
  ggplot2::annotate(geom = "text", 
                    x = -4.6, y = 77, 
                    label = "Isthmus Rises",
                    hjust = 0,
                    size = 5) +
  
  # add geo scale manually
  deeptime::coord_geo(dat = list(timeseries,"epochs", andean_timescale), 
                      size = list(2,3,3),
                      pos = list("bottom", "bottom", "bottom"),
                      xlim = c(- max(ape::branching.times(dec_example@phylo)) * 1.01, 
                               3.75),
                      ylim = c(-3.5, 78),
                      height = list(grid::unit(5, "line"),
                                    grid::unit(2, "line"),
                                    grid::unit(2, "line")),
                      abbrv = FALSE,
                      skip = "Holocene",
                      center_end_labels = TRUE,
                      bord = c("right", "top", "bottom"),
                      neg  = TRUE) +
  
  # theme adjustments
  ggplot2::theme(legend.position = c(0.15, 0.7),
                 legend.background = 
                   ggplot2::element_rect(fill="white", color = "black"),
                 legend.text = element_text(size = 10),
                 legend.title = element_blank(),
                 legend.key = element_blank(),
                 text = element_text(size = 10)) 

print(p)
dev.off()
