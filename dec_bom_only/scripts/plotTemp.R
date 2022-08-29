# plot temperature data from Hansen et al 2013 via 
# data from Tierney et al 2020 GitHub 
# https://github.com/jesstierney/PastClimates
setwd("~/Documents/bom_phy/analyses/dec_bom_only/")

temp2 <- read.csv("data/tierney2020_temp.csv")

y <- diff(temp2$surfaceTemperature)^2
x <- temp2$AgeMa[-1]
temp2$tempDiffs <- c(NA, y)

window <- 10
temp2$var <- NA

for (i in 1:nrow(temp2)) {
  start <- i 
  end <- start + (window - 1)
  temp2$var[i] <- var(temp2$surfaceTemperature[start:end])
}

temp2_trimmed <- temp2[temp2$AgeMa < 10, ]

p1 <- ggplot(temp2_trimmed, aes(AgeMa, surfaceTemperature)) +
  geom_line(color = "firebrick4") +
  xlim(c(max(ape::branching.times(dec_example@phylo)) * 1.01, 
         0)) +
  scale_x_reverse(expand = c(0, 0)) + 
  scale_y_continuous(position = "right",
                     name = "GAT (°C)",
                     expand = c(0, 0)) +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    axis.line.y = element_line(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    )

png("figures/cen_temp.png", height = 1.5, 
    width = 8.5, units = "in", res = 350,
    bg = "transparent")
p1
dev.off()
