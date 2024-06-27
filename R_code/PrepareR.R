library(Matrix)
library(MASS)
library(fda)
library(lme4)
library(pedigreemm)
library(mvnfast)
library(ggplot2)
library(plotly)
library(fBasics)
library(roahd)
library(robustlmm)
library(lme4pureR)

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

fit_genetic_fmm <- function(formula, data, A, phi, control = lmerControl()) {
  #' This function uses the lme4 package to fit a linear mixed-effect model to genetic data, 
  #' with a specified additive genetic relationship matrix A.
  #' In this particular format, we fit both fixed effects and random effects using the same
  #' basis functions (principal components obtained from running FPCA).
  #'
  #' @param formula a two-sided linear formula object describing both the fixed-effects
  #' and random-effects of the model (as the same form used in lmer).
  #' @param data a data frame containing the variables named in formula.
  #' @param A a sparse matrix: an additive genetic relationship matrix 
  #' which models the genetic relationship in the dataset.
  #' @param phi functional basis: a matrix where each column represents a basis element.
  #' @param control lmerControl object: control parameters for the optimizer.
  #' @return returns a fitted mixed-effect model
  
  # Load required packages
  library(lme4)
  library(Matrix)
  
  # Cholesky decomposition of A
  LA <- as(t(chol(A)), "sparseMatrix")
  p <- dim(phi)[2] # number of elements of the functional basis
  I_p <- as(diag(p), "sparseMatrix")
  M <- kronecker(LA, I_p) # used to update the genetic design matrix Z_G = ZM
  
  # Define the mixed-model formula
  fmmParsedForm <- lFormula(formula, data = data, control = control)
  
  # Compute the random-effect matrix
  Z_pre <- t(fmmParsedForm$reTrms$Zt)
  ZE <- Z_pre[, 1:dim(M)[1]] # environmental random-effect matrix
  ZG <- Z_pre[, 1:dim(M)[1]] %*% M # update the genetic-random effect matrix
  Z <- cbind(ZG, ZE) # the updated random effect design matrix
  
  # Modularisation
  fmmParsedForm$reTrms$Zt <- t(Z) # Update Z in the reTrms term
  fmmDevFun <- do.call(mkLmerDevfun, fmmParsedForm) # update the objective function
  
  fmmOptimize <- optimizeLmer(devfun = fmmDevFun, control = control) # update the optimisation module
  
  # Return the mixed-effect model
  fmm <- mkMerMod(rho = environment(fmmDevFun), opt = fmmOptimize, reTrms = fmmParsedForm$reTrms, fr = fmmParsedForm$fr)
  
  return(fmm)
}
