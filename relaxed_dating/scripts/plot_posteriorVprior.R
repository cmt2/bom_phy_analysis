library(RevGadgets)
library(treeio)
library(ggplot2)
library(deeptime)

setwd("~/Documents/bom_phy/analyses/relaxed_dating/")

#### read in traces ####

# posterior trace (empirical)
posterior_trace <- readTrace(c("output/relaxed_clock_partition_45_2.0.log",
                               "output/relaxed_clock_partition_45_2.0_run2.log"))
# combine runs 
posterior_trace_combined <- do.call("rbind", posterior_trace)

# trace from under prior 
prior_trace <- readTrace(c("output/relaxed_clock_partition_45_2.0_underPrior.log",
                           "output/relaxed_clock_partition_45_2.0_underPrior_2.log"))
prior_trace_combined <- combineTraces(prior_trace)[[1]]

# trace from just fossil (no molecular data)
fossil_trace <- readTrace(c("output/relaxed_clock_partition_45_2.0_onlyfossil.log",
                            "output/relaxed_clock_partition_45_2.0_onlyfossil_2.log"))
fossil_trace_combined <- combineTraces(fossil_trace)[[1]]


fossil_trace_long <- rbind(fossil_trace_combined, 
                           dplyr::slice_sample(fossil_trace_combined,
                                               n = (nrow(posterior_trace_combined) - 
                                                      nrow(fossil_trace_combined)),
                                               replace = TRUE)
                           )


new_trace <- data.frame(posterior_Bom = posterior_trace_combined$tmrca_Bomarea,
                        prior_Bom = prior_trace_combined$tmrca_Bomarea, 
                        fossil_Bom = fossil_trace_long$tmrca_Bomarea,
                        
                        posterior_Als = posterior_trace_combined$tmrca_Alstroemeria,
                        prior_Als = prior_trace_combined$tmrca_Alstroemeria, 
                        fossil_Als = fossil_trace_long$tmrca_Alstroemeria,
                        
                        posterior_Alstoid = posterior_trace_combined$tmrca_Alstroemerioideae,
                        prior_Alstoid = prior_trace_combined$tmrca_Alstroemerioideae, 
                        fossil_Alstoid = fossil_trace_long$tmrca_Alstroemerioideae,
                        
                        posterior_Luz = posterior_trace_combined$tmrca_Luzuriaga,
                        prior_Luz = prior_trace_combined$tmrca_Luzuriaga, 
                        fossil_Luz = fossil_trace_long$tmrca_Luzuriaga,
                        
                        posterior_root = posterior_trace_combined$root_time,
                        prior_root = prior_trace_combined$root_time,
                        fossil_root = fossil_trace_long$root_time)

new_trace <- list(trace1 = new_trace)

t_als <- plotTrace(new_trace, 
                   vars = c("prior_Als", "fossil_Als", "posterior_Als"),
                   color = c("prior_Als" = "burlywood3", 
                             "fossil_Als" = "cadetblue3", 
                             "posterior_Als" = "darkblue"))[[1]] +
  scale_color_manual(values = c("burlywood3", "cadetblue3", "darkblue"),
                     labels = c("Tree prior", "Calibrated prior", "Posterior")) +
  scale_x_reverse(limits = c(125,0)) +
  ggtitle("Alstroemeria crown age") +
  theme(legend.title = element_blank(),
        legend.position = c(.3 ,.75),
        axis.title.x = element_blank(),
        plot.background = element_rect(colour = "black", fill="white", size=1.5))

t_alstoid <- plotTrace(new_trace, 
                       vars = c("prior_Alstoid","fossil_Alstoid", "posterior_Alstoid"),
                       color = c("prior_Alstoid" = "burlywood3", 
                                 "fossil_Alstoid" = "cadetblue3",
                                 "posterior_Alstoid" = "darkblue"))[[1]] +
  scale_color_manual(values = c("burlywood3", "cadetblue3", "darkblue"),
                     labels = c("Tree prior", "Calibrated prior", "Posterior")) +
  scale_x_reverse(limits = c(125,0)) +
  ggtitle("Alstroemerioideae crown age") +
  theme(legend.title = element_blank(),
        legend.position = c(.30,.75),
        axis.title.x = element_blank(),
        plot.background = element_rect(colour = "black", fill="white", size=1.5))

t_luz <- plotTrace(new_trace, 
                   vars = c("prior_Luz","fossil_Luz", "posterior_Luz"),
                   color = c("prior_Luz" = "burlywood3", 
                             "fossil_Luz" = "cadetblue3",
                             "posterior_Luz" = "darkblue"))[[1]] +
  scale_color_manual(values = c("burlywood3", "cadetblue3", "darkblue"),
                     labels = c("Tree prior", "Calibrated prior", "Posterior")) +
  scale_x_reverse(limits = c(125,0)) +
  ggtitle("Luzuriageae crown age") +
  theme(legend.title = element_blank(),
        legend.position = c(.3,.75),
        axis.title.x = element_blank(),
        plot.background = element_rect(colour = "black", fill="white", size=1.5))

t_root <- plotTrace(new_trace, 
                   vars = c("prior_root","fossil_root", "posterior_root"),
                   color = c("prior_root" = "burlywood3",
                             "fossil_root" = "cadetblue3", 
                             "posterior_root" = "darkblue"))[[1]] +
  scale_color_manual(values = c("burlywood3","cadetblue3", "darkblue"),
                     labels = c("Tree prior", "Calibrated prior", "Posterior")) +
  scale_x_reverse(limits = c(125,0)) +
  ggtitle("Root age") +
  theme(legend.title = element_blank(),
        legend.position = c(.3,.75),
        plot.background = element_rect(colour = "black", fill="white", size=1.5))

t_bom <- plotTrace(new_trace, 
                   vars = c("prior_Bom","fossil_Bom", "posterior_Bom"),
                   color = c("prior_Bom" = "burlywood3",
                             "fossil_Bom" = "cadetblue3", 
                             "posterior_Bom" = "darkblue"))[[1]] +
  scale_color_manual(values = c("burlywood3","cadetblue3", "darkblue"),
                     labels = c("Tree prior", "Calibrated prior", "Posterior")) +
  scale_x_reverse(limits = c(125,0)) +
  ggtitle("Bomarea crown age") +
  theme(legend.title = element_blank(),
        legend.position = c(.3,.75),
        axis.title.x = element_blank(),
        plot.background = element_rect(colour = "black", fill="white", size=1.5))

png("figures/priorvsposterior_bom.png", height = 3, width = 6, units = "in", res = 350)
t_bom
dev.off()

png("figures/priorvsposterior_als.png", height = 3, width = 6, units = "in", res = 350)
t_als
dev.off()

png("figures/priorvsposterior_alstoid.png", height = 3, width = 6, units = "in", res = 350)
t_alstoid
dev.off()

png("figures/priorvsposterior_luz.png", height = 3, width = 6, units = "in", res = 350)
t_luz
dev.off()

png("figures/priorvsposterior_root.png", height = 3, width = 6, units = "in", res = 350)
t_root
dev.off()



