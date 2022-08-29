# make specimen table
library(knitr)
library(kableExtra)
library(tidyverse)

# gather data 
astral_tree <- ape::read.tree(file = "~/Documents/bom_phy/analyses/astral_multilocus2/output/processed_output/Bomarea_multilocus2_BS10_best.tre")
dated_tree <- RevGadgets::readTrees("~/Documents/bom_phy/analyses/relaxed_dating/output/relaxed_clock_partition_45_2.0_map.tree")[[1]][[1]]
data <- read.csv("~/Documents/bom_phy/specific_metadata/sample_info.csv")
loci_per_sample <- read.table("~/Documents/bom_phy/loci/original/metadata/LociPerSample.txt",
                              header = FALSE, sep = "\t", col.names = c("name", "num_loci"))

# prep data 
data$scientific_name <- paste(data$genus, data$species)

# prep loci per sample 
loci_per_sample$extr_num <- matrix(unlist(strsplit(loci_per_sample$name, "\\.")),
                                   ncol = 2, byrow = TRUE)[,2]
loci_per_sample$extr_num[loci_per_sample$extr_num == "Alstroemeria_pallida"] <- "KEW_pal"
loci_per_sample$extr_num[loci_per_sample$extr_num == "Bomarea_ovallei"] <- "KEW_ova"
loci_per_sample$extr_num[loci_per_sample$extr_num == "Alstroemeria_presliana"] <- "KEW_pres"

# combine data 
data$sequenced <- data$extr_num %in% loci_per_sample$extr_num
data$in_chrono <- data$matcher %in% dated_tree@phylo$tip.label

data$sequenced <- ifelse(data$sequenced == TRUE, "Y", "")
data$in_astral <- ifelse(is.na(data$matcher) == FALSE, "Y", "")
data$in_chrono <- ifelse(data$in_chrono == TRUE, "Y", "")

data$year_edited <- paste0("herb. (", data$year, ")")
data$year_edited <- ifelse(data$material_type == "herbarium", data$year_edited, data$material_type)


data[ , c("scientific_name",
          "collection_num", 
          "region", 
          "year_edited", 
          "insitituion_of_origin",
          "sequenced",
          "in_astral",
          "in_chrono")] %>% 
  arrange(scientific_name, region) %>%
  kable("latex", 
        longtable = T, 
        booktabs = T, 
        linesep = "",
        align=c(rep('l',5), rep('c', 3)),
        col.names = c("Scientific Name",
                      "Number",
                      "Region",
                      "Material",
                      "Institution",
                      "Sequenced",
                      "In phylogram",
                      "In chronogram"),
        caption = "Voucher information for specimens used in phylogenetic reconstruction.") %>%
  column_spec(1, italic = T) %>%
  add_header_above(c(" " = 1, "Collection Information" = 4, "Included" = 3)) %>%
  kable_styling(latex_options = c("striped",
                                  "repeat_header") ) %>%
  row_spec(0, bold=TRUE)

