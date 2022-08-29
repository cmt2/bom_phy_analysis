# test for convergence in relaxed clock dating of Bomarea
library(ape)
library(coda)
tres_text1 <- read.table("~/Desktop/relaxed_dating/output/relaxed_clock_partition_45_2.0.trees",
                        header = T)
tres_text2 <- read.table("~/Desktop/relaxed_dating/output/relaxed_clock_partition_45_2.0_run2.trees",
                         header = T)
trees1 <- read.tree(text = tres_text1$phylogeny)
trees2 <- read.tree(text = tres_text2$phylogeny)

# test if indexed the same
get_mat <- function(x) {x$edge}
get_tips <- function(x) {x$tip.label}
trees1_mat <- lapply(trees1, get_mat)
trees1_tips <- lapply(trees1, get_tips)
trees2_mat <- lapply(trees2, get_mat)
trees2_tips <- lapply(trees2, get_tips)
length(unique(c(trees1_mat, trees2_mat)))
length(unique(c(trees1_tips, trees2_tips)))

# get branch lengths from trees
get_lengths <- function(x) {x$edge.length}
trees1_lengths <- lapply(trees1, get_lengths)
trees2_lengths <- lapply(trees2, get_lengths)

length1_mat <- matrix(unlist(trees1_lengths), ncol = length(trees1_lengths[[1]]), byrow = T)
length2_mat <- matrix(unlist(trees2_lengths), ncol = length(trees2_lengths[[1]]), byrow = T)

mcmc_lengths1 <- as.mcmc(length1_mat)
mcmc_lengths2 <- as.mcmc(length2_mat)
mcmc_lengths_combined <- as.mcmc(rbind(length1_mat, length2_mat))
mcmc_lengths_list <- as.mcmc.list(list(mcmc_lengths1, mcmc_lengths2))

ess_lengths1 <- effectiveSize(mcmc_lengths1)
ess_lengths2 <- effectiveSize(mcmc_lengths2)
ess_lengths_combined <- effectiveSize(mcmc_lengths_combined)
range(ess_lengths1)
range(ess_lengths2)
range(ess_lengths_combined)

# gelman
gelman.diag(mcmc_lengths_list, multivariate = FALSE, autoburnin = FALSE)

# get ages from trees
trees_heights1 <- lapply(trees1, branching.times)
trees_heights2 <- lapply(trees2, branching.times)

heights_mat1 <- matrix(unlist(trees_heights1), ncol = length(trees_heights1[[1]]), byrow = T)
heights_mat2 <- matrix(unlist(trees_heights2), ncol = length(trees_heights2[[1]]), byrow = T)


mcmc_heights1 <- as.mcmc(heights_mat1)
mcmc_heights2 <- as.mcmc(heights_mat2)
mcmc_heights_combined <- as.mcmc(rbind(heights_mat1, heights_mat2))
mcmc_heights_list <- as.mcmc.list(list(mcmc_heights1, mcmc_heights2))


ess_heights1 <- effectiveSize(mcmc_heights1)
ess_heights2 <- effectiveSize(mcmc_heights2)
ess_heights_combined <- effectiveSize(mcmc_heights_combined)
range(ess_heights1)
range(ess_heights2)
range(ess_heights_combined)

# gelman
gelman.diag(mcmc_heights_list, multivariate = FALSE, autoburnin = FALSE)
