# make labels vector using input nexus and intermediate text file

#nexus_path <- "data/bomarea_ranges_full.nex"
#labels_path <- "output/simple_full.state_labels.txt"


makeStateLabels <- function(nexus_path, labels_path) {
  if (file.exists(nexus_path) == F) stop("nexus_path file doesn't exist. Check your working directory and path.")
  if (file.exists(labels_path) == F) stop("nexus_path file doesn't exist. Check your working directory and path.")
  
  nexus_text <- readLines(nexus_path)
  labels_text <- read.csv(labels_path, colClasses = c("character", "character"))
  
  nexus_key <- nexus_text[which(nexus_text == "Matrix")+1]
  
  #### this will try to parse out the key automatically, but if you've 
  # modified the format you'll have to do this manually 
  
  # remove spaces
  nexus_key <- gsub(" ", "", nexus_key)
  # remove brackes
  nexus_key <- gsub("\\[", "", nexus_key)
  nexus_key <- gsub("\\]", "", nexus_key)
  # remove anything after a colon
  nexus_key <- unlist(strsplit(nexus_key, split = ":"))[[1]]
  
  # did we do this right? 
  print(paste0("We've parsed ", "'", nexus_key, "'", " as the state codes from your nexus file."))
  
  good <- readline(prompt="Is this correct? Enter Y for yes or N for no.")
  if (good == "Y") {
    nexus_key <- unlist(strsplit(nexus_key, split = ""))
    
    convert_int_to_state_label <- function(n) {
      paste0(nexus_key[as.logical(as.integer(unlist(strsplit(n, ""))))], collapse = "")
    }
    
    labels <- labels_text
    # remove 0 state
    labels <- labels[-which(labels$state == "0"),] 
    # get state labels
    labels$state_labels <- unlist(lapply(labels$range,convert_int_to_state_label))
    
    # produce named vector
    state_labels <- labels$state_labels
    names(state_labels) <- labels$state
    return(state_labels)
  } else if (good == "N") {
    improve <- readline(prompt="Do you want to enter the state codes manually? Enter Y for yes or N for no.")
    if (improve == "Y") {
      
      nexus_key <- readline(prompt="Enter your state codes as a single character string in the order of the nexus file binaries:")
      if (is.character(nexus_key) == F) stop("Your state codes must be a single character string. Try again.")
      if (length(nexus_key) != 1) stop("Your state codes must be a single character string. Try again.")
      nexus_key <- unlist(strsplit(nexus_key, split = ""))
      convert_int_to_state_label <- function(n) {
        paste0(nexus_key[as.logical(as.integer(unlist(strsplit(n, ""))))], collapse = "")
      }
      
      labels <- labels_text
      # remove 0 state
      labels <- labels[-which(labels$state == "0"),] 
      # get state labels
      labels$state_labels <- unlist(lapply(labels$range,convert_int_to_state_label))
      
      # produce named vector
      state_labels <- labels$state_labels
      names(state_labels) <- labels$state
      return(state_labels)
      
    } else if (improve == "N") {
      stop("Sorry, you'll have to create the state labels vector manually")
      } else stop("Try again.")
  } 
}

