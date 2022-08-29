# match node labels and data in treedata object 

# tree1 is a list of list of treedata objects
# tree2 is an optional second list of of list of treedata objects
# if tree2 is not provided, the function will map node labels in 
# tree1@phylo with the data object 
# if tree2 is provided, it will use node labels from tree2 
# recursive lengths of tree1 and tree2 must match 

matchNodeLabels <- function(tree1, tree2 = NULL) {
  
  # if no tree2 is provided, use the labels from tree1
  if (is.null(tree2)) {
    tree2 <- tree1
  }
  
  # set up dummy output data structure 
  output_tree <- tree1

  
  for (i in 1:length(tree1)) {
    for (j in 1:length(tree1[[i]])) {
      
      input_tree <- tree1[[i]][[j]]
      
      label_tree <- tree2[[i]][[j]]
     
      # change label_tree to tibble
      label_tib <- treeio::as_tibble(label_tree)
      # get the labels and make them into posteriors
      label_tib$posterior <- as.numeric(label_tib$label)
      # convert back to tree 
      label_new <- treeio::as.treedata(label_tib)

      # remove posteriors if present from input 
      if ("posterior" %in% colnames(input_tree@data)) {
        input_tree@data$posterior <- NULL
      }
      
      # add the correct posteriors from label_tree, 
      # linking with strings of descending tips
      
      label_new@data$tips <- NA
      # this loop makes a new column with strings of tip labels as unique node IDs
      for (n in 1:nrow(label_new@data)) {
        
        tips <- geiger::tips(label_new@phylo, 
                             as.integer(label_new@data$node[n]))
        tips <- paste0(sort(tips), collapse = "; ")
        label_new@data$tips[n] <- tips
      }
      
      # same as above for the new tree 
      input_tree@data$tips <- NA
      for (n in 1:nrow(input_tree@data)) {
        
        tips <- geiger::tips(input_tree@phylo, 
                             as.integer(input_tree@data$node[n]))
        tips <- paste0(sort(tips), collapse = "; ")
        input_tree@data$tips[n] <- tips
      }
      
      # move posteriors over to input tree, matching with the tips column 
      input_tree@data <- dplyr::full_join(input_tree@data, label_new@data[ ,c("tips","posterior")])
      
      # move into the output data structure 
      output_tree[[i]][[j]] <- input_tree
      
    }
  }
  return(output_tree)
}
