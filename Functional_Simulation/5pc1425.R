## load packages 
library(Matrix)
library(MASS)
library(fda)
library(lme4)
library(pedigreemm)
library(mvnfast)
library(ggplot2)

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
  #' @param nbasis numeber: number of basis used to fit the random effects
  #' @param control lmerControl object: control parameters for the optimizer.
  #' @return returns a fitted mixed-effect model
  
  # Load required packages
  library(lme4)
  library(Matrix)
  
  # Cholesky decomposition of A
  LA <- as(t(chol(A)), "sparseMatrix")
  I_p <- as(diag(nbasis), "sparseMatrix")
  MA <- kronecker(LA, I_p) # used to update the genetic design matrix Z_G = ZM
  
  # Define the mixed-model formula
  fmmParsedForm <- lFormula(formula=formula, data = data, control = control)
  
  # Compute the random-effect matrix
  Z_pre <- t(fmmParsedForm$reTrms$Zt)
  ZE <- Z_pre[, 1:dim(MA)[1]] # environmental random-effect matrix
  ZG <- Z_pre[, 1:dim(MA)[1]] %*% MA # update the genetic-random effect matrix
  Z <- cbind(ZG, ZE) # the updated random effect design matrix
  
  # Modularisation
  fmmParsedForm$reTrms$Zt <- t(Z) # Update Z in the reTrms term
  fmmDevFun <- do.call(mkLmerDevfun, fmmParsedForm) # update the objective function
  
  fmmOptimize <- optimizeLmer(devfun = fmmDevFun, control = control) # update the optimisation module
  
  # Return the mixed-effect model
  fmm <- mkMerMod(rho = environment(fmmDevFun), opt = fmmOptimize, reTrms = fmmParsedForm$reTrms, fr = fmmParsedForm$fr)
  
  return(fmm)
}


load("GeneticRelationMatrix.Rdata")
load("TC_fixed_effect.Rdata")

N <- 873 # total number of individuals
n <- 14 # number of measurements per individual
nbasis <- 5 # number of basis
npc <- 6
ngroups <- 50

time_rang <- seq(0,1,length=n) # time points t_j
timefine <- seq(0,1, length = 100)

basisObj <- create.bspline.basis(c(0,1), nbasis = nbasis, norder = 4)
basis <- eval.basis(time_rang, basisObj)

### true covariance matrix
C5_true <- matrix(c(750, 10 ,130, 80, 250,
                    10, 800, 30, 15, 40,
                    130, 30, 700, 50, 130,
                    80, 15, 50, 420, 50,
                    250, 40, 130, 50, 330), nrow = 5, byrow = T)

### residual variance
sigma2 <- 25

set.seed(123)

## generate data

C_gen_para <- as(kronecker(A, C5_true), "sparseMatrix") 

### reparameterised environmental covariance 
I_N <- as(diag(N), "dgCMatrix")
C_env_para <- as(kronecker(I_N, C5_true), "sparseMatrix")

### distribution of genetic effect
mu <- rep(0, dim(C_gen_para)[1])
alpha <- rmvn(n=ngroups, mu = mu, sigma = C_gen_para)

### distribution of environmental effect
gamma <- rmvn(n=ngroups, mu=mu, sigma = C_env_para)

### distribution of the error vector
mu_res <- rep(0, N*n)
I_n <- as(diag(N*n), "dgCMatrix")
res_cov <- as(sigma2 * I_n, "dgCMatrix") 

res <- rmvn(n=ngroups, mu=mu_res, sigma = res_cov)# error vector

## generating curve data and fit data to mixed-effect model
gpf <- gl(N,n) # grouping factor
Ji <- t(as(gpf, Class = "sparseMatrix")) # Indicator matrix of grouping factor indices

Xi <- as(cbind(rep(basis[,1], times = N),rep(basis[,2], times = N),rep(basis[,3], times = N),rep(basis[,4], times = N),rep(basis[,5], times = N)), Class = "sparseMatrix" )# raw random effect matrix

Zi <- t(KhatriRao(t(Ji), t(Xi))) # random effect design matrix

Y_rand<- matrix(0, N*n, ngroups)
for (k in 1:ngroups){
  Yk <- Zi %*% alpha[k,] + Zi %*% gamma[k,] + res[k,]
  Y_rand[,k] <- Yk@x
}

f_true <- convert_to_basisfunctions(timefine, fixed_effect, time_rang)
Y_51425 <- Y_rand + rep(f_true, times=N) # simulated data (fixed + random)
save(Y_51425, file = "Y_51425.Rdata")

lam_5pc1425 <- matrix(0, N, ngroups) # store smoothing parameter
fef_5pc1425 <- matrix(0,n,ngroups) # store estimated fixed-effect
pcs_5pc1425 <- array(data = 0, c(n,npc,ngroups)) # store principal components 
CG_fun_5pc1425 <- array(data = 0, c(n,n,ngroups)) # store genetic covariance function
CE_fun_5pc1425 <- array(data = 0, c(n,n,ngroups)) # store environmental covariance function
res_5pc1425 <- rep(0,ngroups) # store estimated residual
sfit_5pc1425 <- rep(NA, ngroups) # whether the cov is singular

## Model fitting

uniqueIds <- seq(1,N, length=N)
id <- rep(uniqueIds, each=n)

for (k in 1:ngroups){
  yk_hat <- matrix(0,n,N) # directly smooth data
  lamk <- rep(0, N)
  
  Yk <- split(Y_51425[,k],id)
  for (i in 1:N) {
    ss <- smooth.spline(time_rang, Yk[[i]], cv = FALSE)
    lamk[i] <- ss$lambda
    yk_hat[,i] <- ss$y
  }
  lam_5pc1425[,k] <- lamk
  
  ## FPCA
  fpcaobj <- prcomp(t(yk_hat), center = TRUE, retx = TRUE, rank. = npc)
  pcs <- fpcaobj$rotation
  pcs_5pc1425[,,k] <- pcs
  
  b1 <- rep(pcs[,1], times = N) 
  b2 <- rep(pcs[,2], times = N) 
  b3 <- rep(pcs[,3], times = N) 
  b4 <- rep(pcs[,4], times = N) 
  b5 <- rep(pcs[,5], times = N) 
  b6 <- rep(pcs[,6], times = N) 
  
  df_sm <- data.frame(id = id, y= Y_51425[,k], 
                      phi1 = b1, 
                      phi2 = b2, 
                      phi3 = b3, 
                      phi4 = b4,
                      phi5 = b5,
                      phi6 = b6)  ## fit raw data
  
  fform <- y ~ -1 + df_sm$phi1 + df_sm$phi2 + df_sm$phi3 + df_sm$phi4 + df_sm$phi5 + df_sm$phi6 +
    (-1 + df_sm$phi1 + df_sm$phi2 + df_sm$phi3  + df_sm$phi4 + df_sm$phi5 | df_sm$id) + 
    (-1 + df_sm$phi1 + df_sm$phi2 + df_sm$phi3  + df_sm$phi4 + df_sm$phi5 | df_sm$id)
  
  ft <-  fit_genetic_fmm(fform, df_sm, A, nbasis)
  
  ## Extract fixed-effect
  beta <- fixef(ft)
  fef_5pc1425[,k] <- pcs %*% beta
  
  vc <- VarCorr(ft)
  CG <- vc[["df_sm.id"]] # genetic covariance
  CE <- vc[["df_sm.id.1"]] # environmental covariance
  sig2 <- sigma(ft)^2
  
  CG_fun <- pcs[,1:5] %*% CG %*% t(pcs[,1:5]) # estimated gen cov function
  CE_fun <- pcs[,1:5] %*% CE %*% t(pcs[,1:5]) # estimated env cov function
  
  CG_fun_5pc1425[,,k] <- CG_fun# estimated gen cov function
  CE_fun_5pc1425[,,k] <- CE_fun# estimated env cov function
  
  ### Extract estimated residual variance
  res_5pc1425[k] <- sigma(ft)^2
  
  ### Check singular fit
  sfit_5pc1425[k] <- isSingular(ft)
  
  save(lam_5pc1425, file = "lam_5pc1425.Rdata")
  save(pcs_5pc1425, file = "pcs_5pc1425.Rdata")
  save(fef_5pc1425, file = "fef_5pc1425.Rdata")
  save(CG_fun_5pc1425, file = "CG_fun_5pc1425.Rdata")
  save(CE_fun_5pc1425, file = "CE_fun_5pc1425.Rdata")
  save(res_5pc1425, file = "res_5pc1425.Rdata")
  save(sfit_5pc1425, file = "sfit_5pc1425.Rdata")
}



