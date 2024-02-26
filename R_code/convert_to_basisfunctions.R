convert_to_basisfunctions <- function(t, eigenvecs, tout, method = "linear") 
  {
  #'This function converts an eigenvector to an eigenfunction using interpolation.
  #'
  #'@param t a numeric vector containing the fine time grid points used to compute the principal components.
  #'@param eigenvecs matrix of eigenvectors obtained from FPCA.
  #'@param tout an optional set of numeric values specifying where interpolation is to take place.
  #'@param method specifies the interpolation method to be used. Choices are "linear" or "spline".
  #'@return matrix where each column represents an eigenfunction of time.
  
  
  # Initialize an empty matrix to store eigenfunctions
  eigen_functions <- matrix(0, nrow = length(tout), ncol = ncol(eigenvecs))
  
  # Interpolate eigenvectors to the original time points
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
  return(eigen_functions)
}
