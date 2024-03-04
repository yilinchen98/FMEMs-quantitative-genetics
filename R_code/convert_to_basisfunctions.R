convert_to_basisfunctions <- function(t, eigenvecs, tout, method = "linear") 
  {
  #'This function converts an eigenvector to an eigenfunction using interpolation.
  #'
  #'@param t a numeric vector containing the fine time grid points used to compute the principal components.
  #'@param eigenvecs vector or matrix of eigenvectors obtained from FPCA.
  #'@param tout an optional set of numeric values specifying where interpolation is to take place.
  #'@param method specifies the interpolation method to be used. Choices are "linear" or "spline".
  #'@return matrix where each column represents an eigenfunction of time.
  
  
  # Initialize an empty matrix to store eigenfunctions
  
  if (is.vector(eigenvecs) == TRUE){
   
    if (method =="linear"){

      for (i in 1: length(eigenvecs)){
        eigen_functions <- approx(x = t, y = eigenvecs, xout = tout)$y
      }
    }
    else if (method == "spline"){
      for (i in 1: length(eigenvecs)){
        eigen_functions <- spline(x = t, y = eigenvecs, xout = tout, method = "natural")$y
      }
    }
    else{
      stop("Invalid method. Supported methods are 'linear' and 'spline'.")
    }
  }
  
  else{
    eigen_functions <- matrix(0, nrow = length(tout), ncol = ncol(eigenvecs))
    
    for (i in 1:ncol(eigenvecs)) {
      
      if (method == "linear"){
        eigen_functions[,i] <- approx(x = t, y = eigenvecs[,i], xout = tout)$y
      }
      else if (method == "spline"){
        eigen_functions[,i] <- spline(x = t, y = eigenvecs[,i], xout = tout, method="natural")$y
      }
      else {
        stop("Invalid method. Supported methods are 'linear' and 'spline'.")
      }
    }
  }
  
  # Interpolate eigenvectors to the original time points
    
  return(eigen_functions)

}
