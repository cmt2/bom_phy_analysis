# plot Bomarea ranges 
library(maps)
library(tidyverse)
setwd("~/Documents/bom_phy/analyses/dec_bom_only//")

data("wrld_simpl")

# A (central Andes), C (Cen Am), E (Eastern SA), N (northern andes), S (southern andes)
cols <- c("#D81B60", "#1E88E5", "#004D40", "#8214a0", "#FFC107")

png("figures/dory_map_colored_bom_regions.png", width = 12, height = 8, units = "in", res = 350)
maps::map('world', wrap = c(0,360))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "aliceblue")
maps::map('world', 
          wrap = c(0,360), 
          fill = TRUE, 
          add = T, 
          col = "ivory")
maps::map('world', 
          regions = c("Peru", "Bolivia"),
          wrap = c(0,360), 
          fill = TRUE, 
          add = T, 
          col = cols[1])
maps::map('world', 
          regions = c("Mexico", "Guatemala", 
                      "Puerto Rico", "Honduras",
                      "Nicaragua", "Belize", "Cuba",
                      "Haiti", "Dominican Republic",
                      "Panama", "Costa Rica", "El Salvador"),
          wrap = c(0,360), 
          fill = TRUE, 
          add = T, 
          col = cols[2])
maps::map('world', 
          regions = c("Brazil", "Guyana", 
                      "Suriname", "French Guiana",
                      "Paraguay"),
          wrap = c(0,360), 
          fill = TRUE, 
          add = T, 
          col = cols[3])
maps::map('world', 
          regions = c("Venezuela", "Colombia", 
                      "Ecuador"),
          wrap = c(0,360), 
          fill = TRUE, 
          add = T, 
          col = cols[4])
maps::map('world', 
          regions = c("Chile", "Argentina", 
                      "Uruguay"),
          wrap = c(0,360), 
          fill = TRUE, 
          add = T, 
          col = cols[5])
dev.off()