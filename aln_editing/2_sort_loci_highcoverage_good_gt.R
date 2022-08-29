setwd("~/Desktop/bom_phy/cleaned_13April/target/loci_highcoverage_good/")

# read in each gene tree, see if there are dups,
# if so sort into another dir

# dir for gene trees of loci with duplicates
dir.create("gt_output/gt_with_dups")
# dir for loci with duplicates 
dir.create("loci_with_dups")

# list of files
gfiles <- list.files("gt_output")
gfiles <- gfiles[grep(".tre", gfiles)]

# function to remove number suffixes from names
rm_num <- function(x) {
  s <- strsplit(x, split = "_")[[1]]
  paste0(s[1:length(s)-1], collapse = "_")
}

for (i in 1:length(gfiles)) {
  tre <- ape::read.tree(file = paste0("gt_output/",gfiles[i]))
  tips <- tre$tip.label
  tips_no_num <- unlist(lapply(tips, rm_num))
  # if there are dups,
  # (1) move gene trees into new folder
  # (2) move loci into new folder
  if (length(unique(tips_no_num)) != length(tips_no_num)) {
    # genetree 
    file.copy(from = paste0("gt_output/",gfiles[i]),
              to = paste0("gt_output/gt_with_dups/",gfiles[i]))
    file.remove(paste0("gt_output/",gfiles[i]))
    
    # loci 
    loci_name <- gsub(".tre", "", gfiles[i])
    file.copy(from = loci_name,
              to = paste0("loci_with_dups/",loci_name))
    file.remove(loci_name)
  } 
}
