library(RevGadgets)
library(ggplot2)
setwd("~/Documents/bom_phy/analyses/dec_alstr/")

#analysis <- "single"
analysis <- "multitrees"

## plot tree 
file <- paste0("output/dec_", analysis, "_combined.ase.tre")

## source new tip labels
source("~/Documents/bom_phy/scripts/plotting/make_tip_labels.R")

source("scripts/makeStateLabels.R")
# labels that correspond to each region/ possible combination of regions
labs <- makeStateLabels("data/bomarea_ranges.nex", 
                        paste0("output/dec_", analysis, "_state_labels.txt"))
# labs[1:5] <- c("C: Central America", "E: Eastern South America", "N: Andes", 
#                "S: Southern South America", "A: Australasia")
dec_example  <- processAncStates(file , state_labels = labs)


# let's drop a low probability state (EA)
for (i in 1:nrow(dec_example@data)) {
  # check end states 
  for (j in c(1, 3, 5)) {
    if (is.na(dec_example@data[i,j]) == FALSE && dec_example@data[i,j] == "EA") {
      print(paste0("EA in end state. i = ", i, ", j = ", j))
      dec_example@data[i,j] <- NA # change state to EA 
      dec_example@data[i,7] <- as.character( as.numeric(dec_example@data[i,7]) + 
                                                as.numeric(dec_example@data[i,j+1]) )  # change other pp to include it      dec_example@data[i,j+1] <- "0"
    }
  }
  # check start states 
  for (j in c(10, 12, 14)) {
    if (is.na(dec_example@data[i,j]) == FALSE && dec_example@data[i,j] == "EA") {
      print(paste0("EA in start state. i = ", i, ", j = ", j))
      dec_example@data[i,j] <- NA # change state to EA 
      dec_example@data[i,16] <- as.character( as.numeric(dec_example@data[i,16]) + 
                                              as.numeric(dec_example@data[i,j+1]) )  # change other pp to include it
      dec_example@data[i,j+1] <- "0"
    }
  }
}


ncol <- length(dec_example@state_labels) -1


# A, C, E, N, S
colors <- c("#94C08B", "#1E88E5", "#004D40", "#D81B60", "#FFC107",
            "#fa78fa", "#8214a0", "#fae6be", "darkblue", "#a0fa82", "red4","darkorange", 
            "turquoise3",#"slategray1", 
            "thistle1")
# ::show_col(colors)

state_labels <- dec_example@state_labels
state_labels <- state_labels[c(1,2,6,11,14,
                               4,8,10,12,13,15, # get rid of 7 (EA)
                               3,5,9)]
attributes(dec_example)$state_labels <- state_labels
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

# fix tip labels
dec_example@phylo$tip.label <- make_tip_labels(dec_example)

pdf(paste0("figures/plot_pies_", analysis, "_alst.pdf"),height = 11, width = 8.5)
p <- plotAncStatesPie(t = dec_example, 
                      pie_colors = colors, 
                      tip_pie_nudge_x = 4.5, 
                      tip_labels_states_offset = 1.8,
                      tip_labels_size = 3, 
                      node_pie_size = 0.5,
                      tip_pie_size = 0.5, 
                      tip_labels_states = TRUE,
                      cladogenetic = TRUE, 
                      tip_labels_offset = 5.5, 
                      timeline = T,
                      tip_labels_formatted = TRUE,
                      time_bars = FALSE,
                      state_transparency = 1,
                      geo_units = list("eras","epochs")) +
  
  # add geo scale manually
  deeptime::coord_geo(dat = list(andean_timescale, "epochs"), 
                      size = list(2, 3),
                      pos = list("bottom", "bottom"),
                      xlim = c(-75, 40),
                      ylim = c(-4, 92),
                      height = grid::unit(2, "line"),
                      abbrv = TRUE,
                      center_end_labels = TRUE,
                      bord = c("right", "top", "bottom"),
                      neg  = TRUE) +

  # isthmus 
  ggplot2::geom_vline(xintercept = -2.8, linetype = "longdash") +
  ggplot2::annotate(geom = "text", 
                    x = -5, y = 90, 
                    label = "Isthmus Rises",
                    hjust = 1,
                    size = 5) +
  # change state labels in legend 
  scale_color_discrete(labels = c ("A: Australasia",
                                   "C: Central America", 
                                   "E: Eastern South America", 
                                   "N: Andes", 
                                   "S: Southern South America", 
                                   "CN", "EN", "ES", "NA", "NS",
                                   "SA", "CEN", "CNS", "ENS", "other"),
                       breaks = c(names(colors), "other"),
                       type = c(colors, "grey50")) + 
  # theme adjustments
  ggplot2::theme(legend.position = c(0.32, 0.35),
                 legend.background = 
                   ggplot2::element_rect(fill="white", color = "black"),
                 legend.text = element_text(size = 11),
                 legend.title = element_blank(),
                 legend.key = element_blank(),
                 text = element_text(size = 12),
                 legend.spacing.x = unit(15, 'native')) +
  # adjust legend rows 
  guides(colour = guide_legend("State",
                               override.aes = list(size = 4, alpha = 1.0),
                               ncol = 3,
                               byrow = FALSE)) 

# add grey bars manually
x_pos <- -rev(c(0,andean_timescale$max_age))

for (k in 2:(length(x_pos))) {
  box_col <- "white"
  if (k %% 2 == 1)
    box_col <- "gray92"
  box <-
    ggplot2::geom_rect(
      xmin = x_pos[k - 1],
      xmax = x_pos[k],
      ymin = -5 * 5,
      ymax = length(dec_example@phylo$tip.label),
      fill = box_col
    )
p <- gginnards::append_layers(p, box, position = "bottom")
}

print(p)

dev.off()`
