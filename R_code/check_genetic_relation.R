check_genetic_relation <- function(id1, id2){
  #'This function examines how two individuals are genetic related. 
  #'Two individuals are full siblings if they have both sire and dam in common.
  #'Two individuals are half siblings if they either have a sire or dam in common.
  #'Two individuals are not related if they have different sires and dams.
  #'
  #'@id1,id2 numeric vector. First element is the offspring's id. 
  #'Second element is the sire id.
  #'Third element is the dam id.
  #'
  #'@return logic values: TRUE means full/half siblings; FALSE means not related.
  
  # Full/Half siblings
  if (any(id1[2:3] %in% id2[2:3])){
    return(TRUE)
  }
  else{
  return(FALSE)
  }
}