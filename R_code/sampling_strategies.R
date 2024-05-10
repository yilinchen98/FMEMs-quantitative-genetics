diverse_sampling <- function(A, sample_size){
  #'This function maximises genetic diversity in a sub-sample.
  #'
  #'@param A: additive genetic relationship matrix
  #'@param sample_size: number of individuals in a sub-sample.
  #'
  #'@return return selected individuals' IDs
  
  N <- dim(A)[1] # number of total individuals
  indices <- 1:N 
  selected <- sample(indices, 1)  # Start by selecting one individual randomly
  remaining <- setdiff(indices, selected)  # Remaining individuals to choose from
  
  while (length(selected) < sample_size) {
    candidate_selected <- NULL
    optimal_relation_score <- Inf
    
    # Iterate through remaining candidates
    for (candidate in remaining) {
      candidate_relation_score <- mean(A[candidate, selected, drop = FALSE])
      if (candidate_relation_score < optimal_relation_score) {
        candidate_selected <- candidate
        optimal_relation_score <- candidate_relation_score
      }
    }
    
    # Update the selected and remaining lists
    selected <- c(selected, candidate_selected)
    remaining <- setdiff(remaining, candidate_selected)
  }
  
  return(selected)
}

related_sampling <- function(A, sample_size){
  #'This function minimises genetic diversity in a sub-sample.
  #'
  #'@param A: additive genetic relationship matrix
  #'@param sample_size: number of individuals in a sub-sample.
  #'
  #'@return return selected individuals' IDs
  
  N <- dim(A)[1] # number of total individuals
  indices <- 1:N 
  selected <- sample(indices, 1)  # Start by selecting one individual randomly
  remaining <- setdiff(indices, selected)  # Remaining individuals to choose from
  
  while (length(selected) < sample_size) {
    candidate_selected <- NULL
    optimal_relation_score <- -Inf
    
    # Iterate through remaining candidates
    for (candidate in remaining) {
      candidate_relation_score <- mean(A[candidate, selected, drop = FALSE])
      if (candidate_relation_score > optimal_relation_score) {
        candidate_selected <- candidate
        optimal_relation_score <- candidate_relation_score
      }
    }
    
    # Update the selected and remaining lists
    selected <- c(selected, candidate_selected)
    remaining <- setdiff(remaining, candidate_selected)
  }
  
  return(selected)
}