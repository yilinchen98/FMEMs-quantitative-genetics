fit_genetic_fmm <- function(formula, data, A, phi)
  {
  #'This function uses the lme4 package to fit a linear mixed-effect model to genetic data, 
  #'with a specified additive genetic relationship matrix A.
  #'In this particular format, we fit both fixed effects and random effects using the same
  #'basis functions (principal components obtained from running FPCA).
  #'
  #'@param fromula a two-sided linear formula object describing both the fixed-effects
  #'and random-effects of the model (as the same form used in lmer).
  #'@param data an data frame containing the variables named in formula.
  #'@param A a sparse matrix: an additive genetic relationship matrix 
  #'which models the genetic relationship in the dataset.
  #'@param phi functional basis: a matrix where each column represents a basis element.
  #'@return returns a fitted mixed-effect model
  
  # Random effect parameterisation
  
  Lt <- as(chol(A), "dgCMatrix") # cholesky decomposition of A
  p <- dim(phi)[2] # number of elements of the functional basis
  I_p <- as(diag(p), "dgCMatrix")
  M <- kronecker(Lt, I_p) # used to update the genetic design matrix Z_E = ZM
  
  # Fit mixed-effect model
  
  ## define the mixed-model formula
  fmmParsedForm <- lFormula(formula, data=data)
  
  ### Compute the random-effect matrix
  Z_pre <- t(fmmParsedForm$reTrms$Zt)
  ZE <- Z_pre[,1:dim(M)[1]] %*% M # update the genetic-random effect matrix
  ZG <- Z_pre[,1:dim(M)[1]] # environmental random-effect matrix
  Z <- cbind(ZG, ZE) # the updated random effect design matrix
  
  ### Modularisation
  fmmParsedForm$reTrms$Zt <- t(Z) # Update Z in the reTrms term
  fmmDevFun <- do.call(mkLmerDevfun,fmmParsedForm) # update the objective function
  fmmOpitimize <- optimizeLmer(devfun=fmmDevFun)# update the optimisation module
  
  ### returns the mixed-effect model
  fmm <- mkMerMod(rho=environment(fmmDevFun),opt=fmmOpitimize, reTrms=fmmParsedForm$reTrms, fr=fmmParsedForm$fr)

  return(fmm)
}