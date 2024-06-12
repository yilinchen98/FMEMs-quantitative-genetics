library(Matrix)
library(MASS)
library(fda)
library(fdasrvf)
library(lme4)
library(pedigreemm)
library(mvnfast)
library(ggplot2)
library(plotly)
library(fBasics)
library(latex2exp)
library(roahd)

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
  require(lme4)
  require(Matrix)
  
  L <- as(t(chol(A)), "dgCMatrix") # cholesky decomposition of A
  p <- dim(phi)[2] # number of elements of the functional basis
  I_p <- as(diag(p), "dgCMatrix")
  M <- kronecker(L, I_p) # used to update the genetic design matrix Z_E = ZM
  
  # Fit mixed-effect model
  
  ## define the mixed-model formula
  fmmParsedForm <- lFormula(formula, data=data)
  
  ### Compute the random-effect matrix
  Z_pre <- t(fmmParsedForm$reTrms$Zt)
  ZE <- Z_pre[,1:dim(M)[1]] # environmental random-effect matrix
  ZG <- Z_pre[,1:dim(M)[1]] %*% M # update the genetic-random effect matrix
  Z <- cbind(ZG, ZE) # the updated random effect design matrix
  
  ### Modularisation
  fmmParsedForm$reTrms$Zt <- t(Z) # Update Z in the reTrms term
  fmmDevFun <- do.call(mkLmerDevfun,fmmParsedForm) # update the objective function
  fmmOpitimize <- optimizeLmer(devfun=fmmDevFun)# update the optimisation module
  
  ### returns the mixed-effect model
  fmm <- mkMerMod(rho=environment(fmmDevFun),opt=fmmOpitimize, reTrms=fmmParsedForm$reTrms, fr=fmmParsedForm$fr)
  
  return(fmm)
}

convert_to_basisfunctions <- function(t, eigenvecs, tout) 
{
  #'This function converts an eigenvector to an eigenfunction using interpolation.
  #'
  #'@param t a numeric vector containing the fine time grid points used to compute the principal components.
  #'@param eigenvecs vector or matrix of eigenvectors obtained from FPCA.
  #'@param tout an optional set of numeric values specifying where interpolation is to take place.
  #'@return matrix where each column represents an eigenfunction of time.
  
  
  # Initialize an empty matrix to store eigenfunctions
  
  if (is.vector(eigenvecs) == TRUE){
    for (i in 1: length(eigenvecs)){
      eigen_functions <- approx(x = t, y = eigenvecs, xout = tout)$y
    }
  }
  
  else{
    eigen_functions <- matrix(0, nrow = length(tout), ncol = ncol(eigenvecs))
    
    for (i in 1:ncol(eigenvecs)) {
      eigen_functions[,i] <- approx(x = t, y = eigenvecs[,i], xout = tout)$y
    }
  }
  # Interpolate eigenvectors to the original time points
  
  return(eigen_functions)
  
}

compute_scaling_para <- function(warping_data, m = 3){
  #'
  #'@param warping_data fdawarp object from time_warping of aligned data
  #'@param m number of principal components 
  #'
  #'@return C positive scaling parameter which joins the curve data and warping functions
  
  require(fdasrvf)
  
  f0 <- warping_data$f0 # original curve
  fn <- warping_data$fn # aligned curves
  time <- warping_data$time # time grid which curve data evaluated on
  gam <- warping_data$warping_functions # warping functions
  Tgam <- SqrtMean(gam) 
  vec <- Tgam$vec # shooting vectors
  
  jointFPCAobj <- function(fn, vec, C=1, m=m){
    n_col <- ncol(fn)
    n_row <- nrow(fn)
    time_new <- seq(0,2, length.out = n_row + nrow(vec))
    ## joint curves and shooting vector
    g <- rbind(fn, C*vec)
    ## run FPCA on g
    fpcaobj <- prcomp(x=t(g), retx = TRUE, center = TRUE, rank. = m)
    eigenvecs <- fpcaobj$rotation # principal components
    scores <- fpcaobj$x # component scores
    
    mu_fn <- rowMeans(fn) # mean of aligned curves
    fn_hat <- matrix(0, n_row, n_col) # estimated aligned curves
    for (i in 1:n_col ){
      fn_hat[,i] <- mu_fn + scores[i,] %*% t(eigenvecs[1:length(time),])
    }
    
    vec_hat <- matrix(0, n_row, n_col) # estimated shooting vectors
    for (i in 1:n_col){
      vec_hat[,i] <- (scores[i,]/C) %*% t(eigenvecs[length(time)+1:length(time),])
    }
    
    gam_hat <- v_to_gam(vec_hat) # estimated warping functions
    
    return(list(f0 = f0, fn_hat = fn_hat, vec_hat = vec_hat, gam_hat = gam_hat, eigenvecs = eigenvecs, time = time))
  }
  
  compute_C <- function(C, warping_data, m) {
    
    require(pracma)
    
    out.pca <- jointFPCAobj(fn, vec, C = C, m = m)
    f0 <- out.pca$f0
    fn_hat <- out.pca$fn_hat
    gam_hat <- out.pca$gam_hat
    n_col <- ncol(fn_hat)
    time <- out.pca$time
    d <- rep(0, n_col)
    for (i in 1:n_col) {
      tmp <- warp_f_gamma(fn_hat[,i], time, invertGamma(gam_hat[,i]))
      d[i] <- sum(trapz(time, (tmp - f0[, i])^2))
    }
    return(sum(d^2)/n_col)
  }
  
  C <- stats::optimize(compute_C, c(0, 100), warping_data = warping_data, m = m)$minimum
  
  return(C)
}

compute_C_2 <- function(warping_data, m = 3){
  #'
  #'@param warping_data fdawarp object from time_warping of aligned data
  #'@param m number of principal components 
  #'
  #'@return C positive scaling parameter which joins the curve data and warping functions
  
  require(fdasrvf)
  
  f0 <- warping_data$f0 # original curve
  fn <- warping_data$fn # aligned curves
  time <- warping_data$time # time grid which curve data evaluated on
  gam <- warping_data$warping_functions # warping functions
  Tgam <- SqrtMean(gam) 
  vec <- Tgam$vec # shooting vectors
  
  jointFPCAobj <- function(fn, vec, C=1, m=m){
    n_col <- ncol(fn)
    n_row <- nrow(fn)
    time_new <- seq(0,2, length.out = n_row + nrow(vec))
    ## joint curves and shooting vector
    g <- rbind(fn, C*vec)
    ## run FPCA on g
    fpcaobj <- prcomp(x=t(g), retx = TRUE, center = TRUE, rank. = m)
    eigenvecs <- fpcaobj$rotation # principal components
    scores <- fpcaobj$x # component scores
    
    mu_fn <- rowMeans(fn) # mean of aligned curves
    fn_hat <- matrix(0, n_row, n_col) # estimated aligned curves
    for (i in 1:n_col ){
      fn_hat[,i] <- mu_fn + scores[i,] %*% t(eigenvecs[1:length(time),])
    }
    
    vec_hat <- matrix(0, n_row, n_col) # estimated shooting vectors
    for (i in 1:n_col){
      vec_hat[,i] <- (scores[i,]/C) %*% t(eigenvecs[length(time)+1:length(time),])
    }
    
    gam_hat <- v_to_gam(vec_hat) # estimated warping functions
    
    return(list(f0 = f0, fn = fn, fn_hat = fn_hat, vec_hat = vec_hat, gam_hat = gam_hat, eigenvecs = eigenvecs, time = time))
  }
  
  compute_C <- function(C, warping_data, m) {
    
    require(pracma)
    
    out.pca <- jointFPCAobj(fn, vec, C = C, m = m)
    f0 <- out.pca$f0
    fn_hat <- out.pca$fn_hat
    gam_hat <- out.pca$gam_hat
    n_col <- ncol(fn_hat)
    time <- out.pca$time
    d <- rep(0, n_col)
    for (i in 1:n_col) {
      tmp <- warp_f_gamma(fn_hat[,i], time, gam_hat[,i])
      d[i] <- sum(trapz(time, (tmp - fn[, i])^2))
    }
    return(sum(d^2)/n_col)
  }
  
  C <- stats::optimize(compute_C, c(0, 1000), warping_data = warping_data, m = m)$minimum
  
  return(C)
}

## load data
setwd("D:/KCL_2023-2027_PhD/Genetics_FMEM_Project/FMEMs-quantitative-genetics/R_code")
TRFUN25PUP4 = read.delim("TRFUN25PUP4.DAT",header = FALSE)
names(TRFUN25PUP4)<-c("id","sire","dam","trait","x")
df <- data.frame(TRFUN25PUP4)

## calculate genetic relationship matrix
FirstUniqueIdPos <- which(duplicated(df$id) == FALSE)

pos = df$id[FirstUniqueIdPos] # extract ids for all subjects
sire_id = df$sire[FirstUniqueIdPos] # extract ids for sire
dam_id = df$dam[FirstUniqueIdPos] # extract ids for dam

pede <- editPed(sire = sire_id, dam = dam_id, label = pos)
ped<- with(pede, pedigree(label=label, sire=sire, dam=dam))
A <- getA(ped)[163:1035,163:1035]
