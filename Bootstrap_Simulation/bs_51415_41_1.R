.libPaths("/users/k2368651/R/4.3")
## load packages 
library(Matrix)
library(MASS)
library(fda)
library(lme4)
library(pedigreemm)
library(mvnfast)

## load functions
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
  #' and random-effects of the model (as the same form used in lmer).
  #' @param data a data frame containing the variables named in formula.
  #' @param A a sparse matrix: an additive genetic relationship matrix 
  #' which models the genetic relationship in the dataset.
  #' @param nbasis numeber/vector: number of basis used to fit the random effects 
  #'        [basis for genetic, no. of basis for environment]
  #' @param control lmerControl object: control parameters for the optimizer.
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

load("GeneticRelationMatrix.Rdata")

N <- 873 # total number of individuals
n <- 14 # number of measurements per individual
nbasis <- 5 # number of basis
npc <- 6

time_rang <- seq(0,1,length=n) # time points t_j

basisObj <- create.bspline.basis(c(0,1), nbasis = nbasis, norder = 4)
basis <- eval.basis(time_rang, basisObj)

### true covariance matrix
C5_true <- matrix(c(750, 10 ,130, 80, 250,
                    10, 800, 30, 15, 40,
                    130, 30, 700, 50, 130,
                    80, 15, 50, 420, 50,
                    250, 40, 130, 50, 330), nrow = 5, byrow = T)

### residual variance
sigma2 <- 15

## fit group 41
load("CG_hat_41.Rdata")
load("CE_hat_41.Rdata")
load("sig2_hat_41.Rdata")
load("pcs41.Rdata")
load("fef_hat_41.Rdata")

set.seed(1)
## Use estimated covariance matrices CE and CG and residual variance to generate bootstrap samples

### model form
gpf <- gl(N,n) # grouping factor
Ji <- t(as(gpf, Class = "sparseMatrix")) # Indicator matrix of grouping factor indices

Xi <- as(cbind(rep(pcs41[,1], times = N),rep(pcs41[,2], times = N),rep(pcs41[,3], times = N),rep(pcs41[,4], times = N),rep(pcs41[,5], times = N)), Class = "sparseMatrix" )# raw random effect matrix

Zi <- t(KhatriRao(t(Ji), t(Xi))) # random effect design matrix

### reparameterise covariance structure
I_N <- as(diag(N), "sparseMatrix")
I_n <- as(diag(N*n), "sparseMatrix")

CG_para <- as(kronecker(A, CG_hat_41), "sparseMatrix") # reparameterised genetic covariance 
CE_para <- as(kronecker(I_N, CE_hat_41), "sparseMatrix") # reparameterised environmental covariance 

mu_hat <- rep(0, dim(CG_para)[1])
mu_res_hat <- rep(0, N*n)
res_cov_hat <- as(sig2_hat_41 * I_n, "sparseMatrix") 

### generate bootstrap samples
n_iter <- 40 # number of iterations (parallel bootstrap jobs)
uniqueIds <- seq(1,N, length=N)
id <- rep(uniqueIds, each=n)

Y_est41_1 <- matrix(0, n*N,n_iter)

lam_sap41_1 <- matrix(0, N, n_iter) # store the smoothing parameter
pcs_sap41_1 <- array(data = 0, c(n,npc,n_iter)) # store the PCs
fef_sap41_1 <- matrix(0,n,n_iter) # store estimated fixed-effect
CG_fun_sap41_1 <- array(data = 0, c(n,n,n_iter)) # 3d matrix store genetic covariance function of 100 bootstrap replicates
CE_fun_sap41_1 <- array(data = 0, c(n,n,n_iter)) # 3d matrix store environmental covariance function of 100 bootstrap replicates
res_sap41_1 <- rep(0,n_iter) # store estimated residual
sfit_sap41_1 <- rep(NA, n_iter)

for (k in 1:n_iter){
  
  alpha_hat <- rmvn(n=1, mu = mu_hat, sigma = CG_para) # genetic random effect
  gamma_hat <- rmvn(n=1, mu = mu_hat, sigma = CE_para) # environmental random effect
  res_hat <- rmvn(n=1, mu=mu_res_hat, sigma = res_cov_hat)# error vector
  
  ### simulated data
  Yk <-  rep(fef_hat_41, times = N) + Zi %*% t(alpha_hat) + Zi %*% t(gamma_hat) + t(res_hat) 
  Y_est41_1[,k] <- Yk@x
  save(Y_est41_1, file = "bs_sample_51415_41_1.Rdata")
  
  ## data smoothing
  yk_hat <- matrix(0,n,N) 
  lamk <- rep(0, N)
  
  Yk_list <- split(Yk, id)
  for (i in 1:N) {
    ss <- smooth.spline(time_rang, Yk_list[[i]], cv = FALSE)
    lamk[i] <- ss$lambda
    yk_hat[,i] <- ss$y
  }
  lam_sap41_1[,k] <- lamk
  save(lam_sap41_1, file = "lam_sap41_1.Rdata")
  
  ## FPCA
  fpcaobjk <- prcomp(t(yk_hat), center = TRUE, retx = TRUE, rank. = npc)
  pcs <- fpcaobjk$rotation
  pcs_sap41_1[,,k] <- pcs
  save(pcs_sap41_1, file = "pcs_sap41_1.Rdata")
  
  ## Combine the smoothed data and basis functions to a dataframe
  bs1 <- rep(pcs[,1], times = N) 
  bs2 <- rep(pcs[,2], times = N) 
  bs3 <- rep(pcs[,3], times = N) 
  bs4 <- rep(pcs[,4], times = N) 
  bs5 <- rep(pcs[,5], times = N) 
  bs6 <- rep(pcs[,6], times = N) 
  
  df_sm <- data.frame(id = id, y= Yk@x, 
                          phi1 = bs1, 
                          phi2 = bs2, 
                          phi3 = bs3, 
                          phi4 = bs4,
                          phi5 = bs5,
                          phi6 = bs6) 
  
  fform <- y ~ -1 + df_sm$phi1 + df_sm$phi2 + df_sm$phi3 + df_sm$phi4 + df_sm$phi5 + df_sm$phi6 +
    (-1 + df_sm$phi1 + df_sm$phi2 + df_sm$phi3  + df_sm$phi4 + df_sm$phi5 | df_sm$id) + 
    (-1 + df_sm$phi1 + df_sm$phi2 + df_sm$phi3  + df_sm$phi4 + df_sm$phi5 | df_sm$id)
  
  ft <-  fit_genetic_fmm(fform, df_sm, A, nbasis)
  
  ### Extract fixed-effect
  beta <- fixef(ft)
  fef_sap41_1[,k] <- pcs %*% beta
  
  ### Extract covariance
  vc <- VarCorr(ft)
  CG <- vc[["df_sm.id"]] # genetic covariance
  CE <- vc[["df_sm.id.1"]] # environmental covariance
  res_sap41_1[k] <- sigma(ft)^2
  
  CG_fun_sap41_1[,,k] <- pcs[,1:5] %*% CG %*% t(pcs[,1:5]) # estimated gen cov function
  CE_fun_sap41_1[,,k] <- pcs[,1:5] %*% CE %*% t(pcs[,1:5]) # estimated env cov function
  
  sfit_sap41_1[k] <- isSingular(ft)
  
  save(fef_sap41_1, file = "fef_sap41_1.Rdata")
  save(CG_fun_sap41_1, file = "CG_fun_sap41_1.Rdata")
  save(CE_fun_sap41_1, file = "CE_fun_sap41_1.Rdata")
  save(sfit_sap41_1, file = "sfit_sap41_1.Rdata")
  save(res_sap41_1, file = "res_sap41_1.Rdata")

}
