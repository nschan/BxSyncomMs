### read_fasta from the metacoder package.


read_fasta <- function(file_path) {
  # Read raw string
  raw_data <- readr::read_file(file_path)
  
  # Return an empty vector an a warning if no sequences are found
  if (raw_data == "") {
    warning(paste0("No sequences found in the file: ", file_path))
    return(character(0))
  }
  
  # Find location of every header start 
  split_data <- stringr::str_split(raw_data, pattern = "\n>", simplify = TRUE)
  
  # Split the data for each sequence into lines
  split_data <- stringr::str_split(split_data, pattern = "\n")
  
  # The first lines are headers, so remvove those
  headers <- vapply(split_data, FUN = `[`, FUN.VALUE = character(1), 1)
  split_data <- lapply(split_data, FUN = `[`, -1)
  
  # Remove the > from the first sequence. The others were removed by the split
  headers[1] <- sub(headers[1], pattern = "^>", replacement = "")
  
  # Combine multiple lines into single sequences
  seqs <- vapply(split_data, FUN = paste0, FUN.VALUE = character(1), collapse = "")
  
  # Combine and return results 
  return(stats::setNames(seqs, headers))
}