# make tip labels

make_tip_labels <- function(tree_data, formatted = TRUE) {
  
  tip_info <- read.csv("~/Documents/bom_phy/specific_metadata/sample_info.csv")
  
  current_tips <- tree_data@phylo$tip.label
  
  matches <- unlist(lapply(current_tips, name_match_fn, tp = tip_info))
  
  if (formatted == TRUE) {
    new_tip_labels <- paste0("italic(`",
                                       tip_info$genus, "_",
                                       tip_info$species, "`)~",
                                       tip_info$region)
    for (i in 1:length(new_tip_labels)) {
      tiplab <- new_tip_labels[i]
      if (grepl("sp.", tiplab, fixed = T)) {
        tiplab <- gsub("`)", "", tiplab, fixed = T)
        tiplab <- gsub(" ", "~", tiplab, fixed = T)
        tiplab <- gsub("_", "`)~", tiplab, fixed = T)
        
        # get rid of ~ between quotes 
        tiplab_split <- unlist(strsplit(tiplab, "'"))
        tiplab_split[2] <- gsub("~", " ", tiplab_split[2])
        tiplab <- paste0(tiplab_split, collapse = "'")
      }
      new_tip_labels[i] <- tiplab
    }
  } else {
    new_tip_labels <- paste0(tip_info$genus, "_",
                             tip_info$species, "_",
                             tip_info$region)
  }
  
  new_names <- new_tip_labels[matches]
  return(new_names)
}


name_match_fn <- function(x, tp) {
  
  x <- gsub("__", "_", x, fixed = TRUE)
  x <- unlist(strsplit(x, "_"))
  
  # check by matching genus + species 
  gen_sp_new <- paste0(tp$genus, "_", tp$species)
  x_gen_sp <- paste0(x[1:2], collapse = "_")
  match_gen_sp <- which(gen_sp_new == x_gen_sp)
  
  if (length(match_gen_sp)== 1) {
    
    match <- match_gen_sp
    
  } else {
    
    # check by matching collection number
    col_num_new <- tp$collection_num
    col_num_new <- gsub(" ", "", col_num_new, fixed = T)
    x_col_num <- x[length(x)]
    match_col_num <- which(col_num_new == x_col_num)
    
    if (length(match_col_num) == 1) {
      
      match <- match_col_num
      
    } else if (x[length(x)] == "8913"){ 
      
      match <- grep("8913", tp$collection_num)
      
    } else if (x[length(x)] == "Tribble325") {
      
      match <- grep("32-5", tp$collection_num )
      
    } else if (x[2] == "multiflora" & x[3] == "AlzateS") {
      
      match <- which(tp$species == "multiflora" & grepl("Alzate", tp$collection_num)) 
      
    } else {match <- NA}
  } 
  
  return(match)
}



