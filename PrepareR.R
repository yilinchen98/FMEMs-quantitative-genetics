library(Matrix)
library(MASS)
library(fda)
library(lme4)
library(fdasrvf)
library(fdapace)
library(pedigreemm)
library(mvnfast)
library(ggplot2)
library(plotly)
library(fBasics)
library(latex2exp)
library(abind)
library(factoextra)
library(colorspace)
library(viridis)
library(patchwork)
library(tidyr)
library(dplyr)
library(pracma)


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

fit_genetic_fmm <- function(formula, data, A, nbasis, control = lmerControl()) {
  #' This function uses the lme4 package to fit a linear mixed-effect model to genetic data, 
  #' with a specified additive genetic relationship matrix A.
  #' In this particular format, we fit both fixed effects and random effects using the same
  #' basis functions (principal components obtained from running FPCA).
  #'
  #' @param formula a two-sided linear formula object describing both the fixed-effects
  #' and random-effects of the model (in the same form as used in lmer).
  #' @param data a data frame containing the variables named in the formula.
  #' @param A a sparse matrix: an additive genetic relationship matrix 
  #' which models the genetic relationship in the dataset.
  #' @param nbasis number/vector: number of basis used to fit the random effects 
  #'        [basis for genetic, no. of basis for environment]
  #' @param control lmerControl object: control parameters for the optimiser.
  #' @return returns a fitted mixed-effect model
  
  # Load required packages
  library(lme4)
  library(Matrix)
  
  # Cholesky decomposition of A
  LA <- as(t(chol(A)), "sparseMatrix")
  if (length(nbasis) == 1){
    I_p <- as(diag(nbasis), "sparseMatrix")
  }
  else{I_p <- as(diag(nbasis[1]), "sparseMatrix")}
  MA <- kronecker(LA, I_p) # used to update the genetic design matrix Z_G = ZM
  
  # Define the mixed-model formula
  fmmParsedForm <- lFormula(formula=formula, data = data, control = control)
  
  # Compute the random-effect matrix
  Z <- t(fmmParsedForm$reTrms$Zt)
  Z[,1:dim(MA)[1]] <- Z[, 1:dim(MA)[1]] %*% MA # update the genetic-random effect matrix
  
  # Modularisation
  fmmParsedForm$reTrms$Zt <- t(Z) # Update Z in the reTrms term
  fmmDevFun <- do.call(mkLmerDevfun, fmmParsedForm) # update the objective function
  
  fmmOptimize <- optimizeLmer(devfun = fmmDevFun, control = control) # update the optimisation module
  
  # Return the mixed-effect model
  fmm <- mkMerMod(rho = environment(fmmDevFun), opt = fmmOptimize, reTrms = fmmParsedForm$reTrms, fr = fmmParsedForm$fr)
  
  return(fmm)
}

