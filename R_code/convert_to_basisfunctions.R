convert_to_basisfunctions <- function(t, eigenvecs, tout) 
{
  #'This function converts an eigenvector to an eigenfunction using linear interpolation.
  #'
  #'@param t a numeric vector containing the fine time grid points used to compute the principal components.
  #'@param eigenvecs vector or matrix of eigenvectors obtained from FPCA.
  #'@param tout an optional set of numeric values specifying where interpolation is to take place.
  #'@return matrix where each column represents an eigenfunction of time.
  
  # Check if eigenvecs is a vector or matrix
  if (is.vector(eigenvecs)) {
    eigen_functions <- approx(x = t, y = eigenvecs, xout = tout)$y
  } else {
    eigen_functions <- matrix(0, nrow = length(tout), ncol = ncol(eigenvecs))
    for (i in 1:ncol(eigenvecs)) {
      eigen_functions[,i] <- approx(x = t, y = eigenvecs[,i], xout = tout)$y
    }
  }
  
  return(eigen_functions)
}

