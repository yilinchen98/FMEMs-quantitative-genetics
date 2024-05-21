# Functional Bootstrapping Simulation

## load packages 
library(Matrix)
library(MASS)
library(fda)
library(lme4)
library(pedigreemm)
library(mvnfast)
library(ggplot2)
library(plotly)
library(fBasics)
library(latex2exp)

# Extract the estimated genetic, environmental covariance functions and residual variance
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

## generate functional data (50 sets of response as used in functional simulation study)

set.seed(123)

N <- 873 # total number of individuals
n <- 10 # number of measurements per individual
nbasis <- 5 # number of basis

### use b-spline basis
basisObj <- create.bspline.basis(c(0,1), nbasis = nbasis, norder = 4)

### genetic covariance matrix
C_gen <- matrix(c(750, 10 ,130, 80, 250,
                  10, 800, 30, 15, 40,
                  130, 30, 700, 50, 130,
                  80, 15, 50, 420, 50,
                  250, 40, 130, 50, 330), nrow = 5, byrow = T)
### environmental covariance matrix
C_env <- matrix(c(750, 10 ,130, 80, 250,
                  10, 800, 30, 15, 40,
                  130, 30, 700, 50, 130,
                  80, 15, 50, 420, 50,
                  250, 40, 130, 50, 330), nrow = 5, byrow = T)
### residual variance
sigma2 <- 50

### Calculate the genetic relationship matrix A (use the TC dataset)
setwd("D:/KCL_2023-2027_PhD/Genetics_FMEM_Project/FMEMs-quantitative-genetics/R_code")
TRFUN25PUP4 = read.delim("TRFUN25PUP4.DAT",header = FALSE)
names(TRFUN25PUP4)<-c("id","sire","dam","trait","x")
df <- data.frame(TRFUN25PUP4)

FirstUniqueIdPos <- which(duplicated(df$id) == FALSE)

pos = df$id[FirstUniqueIdPos] # extract ids for all subjects
sire_id = df$sire[FirstUniqueIdPos] # extract ids for sire
dam_id = df$dam[FirstUniqueIdPos] # extract ids for dam

pede <- editPed(sire = sire_id, dam = dam_id, label = pos)
ped<- with(pede, pedigree(label=label, sire=sire, dam=dam))
A <- getA(ped)[163:1035,163:1035]

### basis for our functional space
time_rang <- seq(0,1,length=n) # time points t_j
basis <- eval.basis(time_rang, basisObj) # evaluate basis at time points

### true genetic covariance function
CG_fun_true <- basis %*% C_gen %*% t(basis)
#save(CG_fun_true, file="true_G.Rdata")
### true environmental covariance function
CE_fun_true <- basis %*% C_env %*% t(basis) 
#save(CE_fun_true, file = "true_E.Rdata")

### reparameterised genetic covariance 
C_gen_para <- as(kronecker(A, C_gen), "dgCMatrix") 

### reparameterised environmental covariance 
I_N <- as(diag(N), "dgCMatrix")
C_env_para <- as(kronecker(I_N, C_env), "dgCMatrix")

n_groups <- 50
### distribution of genetic effect
mu <- rep(0, dim(C_gen_para)[1])
alpha <- rmvn(n=n_groups, mu = mu, sigma = C_gen_para)

### distribution of environmental effect
gamma <- rmvn(n=n_groups, mu=mu, sigma = C_env_para)

### distribution of the error vector
mu_res <- rep(0, N*n)
I_n <- as(diag(N*n), "dgCMatrix")
res_cov <- as(sigma2 * I_n, "dgCMatrix") 

res <- rmvn(n=n_groups, mu=mu_res, sigma = res_cov)# error vector

### generating curve data
f <- gl(N,n) # grouping factor
Ji <- t(as(f, Class = "sparseMatrix")) # Indicator matrix of grouping factor indices

Xi <- as(cbind(rep(basis[,1], times = N),rep(basis[,2], times = N),rep(basis[,3], times = N),rep(basis[,4], times = N),rep(basis[,5], times = N)), Class = "sparseMatrix" )# raw random effect matrix

Zi <- t(KhatriRao(t(Ji), t(Xi))) # random effect design matrix

Y <- matrix(0, N*n, n_groups) # store simulated data (each column stores one set of simulated response of 873 individuals)
for (j in 1:n_groups){
  Yj <- Zi %*% alpha[j,] + Zi %*% gamma[j,] + res[j,]
  Y[,j] <- Yj@x 
}

## fit 41 set of generated data
time_rang <- seq(0,1,length=n)
timefine <- seq(0,1,length=100) 
uniqueIds <- seq(1,N, length=N)
id <- rep(uniqueIds, each=n)
y_list <- split(Y[,41], id)

### Data smoothing
timefine <- seq(0,1,length=100) # dense time grid
y_hat <- matrix(0, 100, N)

for (i in 1:N){
  ss <- smooth.spline(time_rang, y_list[[i]], cv = FALSE)
  y_hat[,i] <- predict(ss, timefine)$y
}

### FPCA
### fit smoothed data at original time measurement points
trait <- convert_to_basisfunctions(t = timefine, eigenvecs = y_hat,
                                   tout = time_rang)

fpcaobj <- prcomp(x=t(trait), retx = TRUE, center = TRUE, rank. = 3)
pcs <- fpcaobj$rotation # eigen vectors as basis functions for model-fitting

### Combine the smoothed data and basis functions to a data frame
response_test <- c(trait)
basis1 <- rep(pcs[,1], times = N)
basis2 <- rep(pcs[,2], times = N)
basis3 <- rep(pcs[,3], times = N)

df_test <- data.frame(id = id, y= response_test, phi1 = basis1, phi2 = basis2, phi3 = basis3) 

### Fit FMEM
fform <- y ~ 1 + df_test$phi1 + df_test$phi2 + df_test$phi3 + 
  (-1 + df_test$phi1 + df_test$phi2 + df_test$phi3 | df_test$id) + 
  (-1 + df_test$phi1 + df_test$phi2 + df_test$phi3 | df_test$id)

system.time(
  ft <-  fit_genetic_fmm(fform, df_test, A, pcs)
) # user   system   elapsed

summary(ft)

### Extract covariance function
vc <- VarCorr(ft)
CG <- vc[["df_test.id"]] ## estimated genetic covariance
CE <- vc[["df_test.id.1"]] ## estimated environmental covariance
sig2 <- sigma(ft)^2 ## estimated residual variance

### Convert to genetic covariance function
CG_fun_hat <- pcs %*% CG %*% t(pcs)
save(CG_fun_hat, file = "GenFun_hat_31050_41.Rdata")
### environmental covariance function
CE_fun_hat <- pcs%*% CE %*% t(pcs)
save(CE_fun_hat, file = "EnvFun_hat_31050_41.Rdata")

### Trace of genetic covariance
tr_gen_hat <- tr(CG_fun_hat)
save(tr_gen_hat, file = "TrGen_hat_31050_1.Rdata")
### Trace of environmental covariance
tr_env_hat <- tr(CE_fun_hat)
save(tr_env_hat, file = "TrEnv_hat_31050_1.Rdata")

## Use estimated covariance matrices CE and CG and residual variance to generate bootstrap samples

### calculating the design matrix Z
f <- gl(N,n) # grouping factor
J <- t(as(f, Class = "sparseMatrix")) # Indicator matrix of grouping factor indices

X <- as(cbind(rep(pcs[,1], times = N),rep(pcs[,2], times = N),rep(pcs[,3], times = N)), Class = "sparseMatrix" )# raw random effect matrix

Z <- t(KhatriRao(t(J), t(X))) # random effect design matrix

## reparameterise covariance structure
CG_para <- as(kronecker(A, CG), "dgCMatrix") # reparameterised genetic covariance 
CE_para <- as(kronecker(I_N, CE), "dgCMatrix") # reparameterised environmental covariance 

mu_hat <- rep(0, dim(CG_para)[1])
mu_res_hat <- rep(0, N*n)
res_cov_hat <- as(sig2 * I_n, "dgCMatrix") 

## generate bootstrap samples
time_rang <- seq(0,1,length=n)
timefine <- seq(0,1,length=100) # dense time grid
uniqueIds <- seq(1,N, length=N)
id <- rep(uniqueIds, each=n)
y_hat <- matrix(0, 100, N)

n_iter <- 500 # number of iterations

gen_cov_fun <- array(data = 0, c(n,n,n_iter)) # 3d matrix store genetic covariance function of 500 iterations
env_cov_fun <- array(data = 0, c(n,n,n_iter)) # 3d matrix store environmental covariance function of 500 iterations

tr_gen <- rep(0, n_iter) # store the trace of genetic covariance function
tr_env <- rep(0, n_iter) # store the trace of environmental covariance function

for (k in 1:n_iter){
  alpha_hat <- rmvn(n=1, mu = mu_hat, sigma = CG_para) # genetic random effect
  gamma_hat <- rmvn(n=1, mu = mu_hat, sigma = CE_para) # environmental random effect
  res_hat <- rmvn(n=1, mu=mu_res_hat, sigma = res_cov_hat)# error vector
  
  ### simulated data
  Y <- Z %*% t(alpha_hat) + Z %*% t(gamma_hat) + t(res_hat) 
  
  ### data smoothing
  y_list <- split(Y, id)
  for (i in 1:N){
    ss <- smooth.spline(time_rang, y_list[[i]], cv = FALSE)
    y_hat[,i] <- predict(ss, timefine)$y
  }
  
  ### fit smoothed data at original time measurement points
  y <- convert_to_basisfunctions(t = timefine, eigenvecs = y_hat,
                                 tout = time_rang)
  
  ### FPCA
  fpcaobj <- prcomp(x=t(y), retx = TRUE, center = TRUE, rank. = 3)
  pcs <- fpcaobj$rotation
  
  ## Combine the smoothed data and basis functions to a dataframe
  response_test <- c(y)
  basis1 <- rep(pcs[,1], times = N)
  basis2 <- rep(pcs[,2], times = N)
  basis3 <- rep(pcs[,3], times = N)
  
  df_test <- data.frame(id = id, y= response_test, phi1 = basis1, 
                        phi2 = basis2, phi3 = basis3) 
  
  ### Fit FMEM
  fform <- y ~ 1 + df_test$phi1 + df_test$phi2 + df_test$phi3 + 
    (-1 + df_test$phi1 + df_test$phi2 + df_test$phi3 | df_test$id) + 
    (-1 + df_test$phi1 + df_test$phi2 + df_test$phi3 | df_test$id)
  
  ft <-  fit_genetic_fmm(fform, df_test, A, pcs)
  
  ### Estimated covariance
  vc <- VarCorr(ft)
  CG <- vc[["df_test.id"]] # genetic covariance
  CE <- vc[["df_test.id.1"]] # environmental covariance
  
  CG_fun <- pcs %*% CG %*% t(pcs) # estimated gen cov function
  CE_fun <- pcs %*% CE %*% t(pcs) # estimated env cov function
  
  gen_cov_fun[,,k] <- CG_fun
  env_cov_fun[,,k] <- CE_fun
  
  tr_gen[k] <- tr(CG_fun)
  tr_env[k] <- tr(CE_fun)
  
  save(gen_cov_fun, file = "gen_31050_41.Rdata")
  save(env_cov_fun, file = "env_31050_41.Rdata")
  save(tr_gen, file = "tr_gen_31050_41.Rdata")
  save(tr_env, file = "tr_env_31050_41.Rdata")
  
}



