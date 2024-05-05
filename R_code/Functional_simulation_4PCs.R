# Import library
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

# Load functions
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

# Calculate the genetic relationship matrix A (use the TC dataset)
setwd("D:/KCL_2023-2027_PhD/Year 1/Genetics/FMEMs-quantitative-genetics/R_code")
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

# Simulation

## Set true covariance matrices
set.seed(123)

N <- 873 # total number of individuals
n <- 20 # number of measurements per individual
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

## Set distribution of random effect vectors
time_rang <- seq(0,1,length=n) # time points t_j
basis <- eval.basis(time_rang, basisObj) # evaluate basis at time points
C_gen_fun <- basis %*% C_gen %*% t(basis) # true genetic covariance function
C_env_fun <- basis %*% C_env %*% t(basis) # true environmental covariance function

### reparameterised genetic covariance 
C_gen_para <- as(kronecker(A, C_gen), "dgCMatrix") 

### reparameterised environmental covariance 
I_N <- as(diag(N), "dgCMatrix")
C_env_para <- as(kronecker(I_N, C_env), "dgCMatrix")

### distribution of genetic effect
mu <- rep(0, dim(C_gen_para)[1])
alpha <- rmvn(n=50, mu = mu, sigma = C_gen_para)

### distribution of environmental effect
gamma <- rmvn(n=50, mu=mu, sigma = C_env_para)

### distribution of the error vector
mu_res <- rep(0, N*n)
I_n <- as(diag(N*n), "dgCMatrix")
res_cov <- as(sigma2 * I_n, "dgCMatrix") 

res <- rmvn(n=50, mu=mu_res, sigma = res_cov)# error vector

## generating curve data and fit data to mixed-effect model
f <- gl(N,n) # grouping factor
Ji <- t(as(f, Class = "sparseMatrix")) # Indicator matrix of grouping factor indices

Xi <- as(cbind(rep(basis[,1], times = N),rep(basis[,2], times = N),rep(basis[,3], times = N),rep(basis[,4], times = N),rep(basis[,5], times = N)), Class = "sparseMatrix" )# raw random effect matrix

Zi <- t(KhatriRao(t(Ji), t(Xi))) # random effect design matrix

time_rang <- seq(0,1,length=n)
timefine <- seq(0,1,length=100) # dense time grid
uniqueIds <- seq(1,N, length=N)
id <- rep(uniqueIds, each=n)
y_hat <- matrix(0, 100, N)
gen_err4_20 <- rep(0, 50) ## store the error on genetic covariance in Frobenius norm
env_err4_20 <- rep(0, 50) ## store the error on environmental covariance in Frobenius norm
tr_gen4_20 <- rep(0,50) ## store the trace of genetic covariance function
tr_env4_20 <- rep(0,50) ## store the trace of environmental covariance function

system.time(
for (j in 31:50){
  ### generating curve data
  Y <- Zi %*% alpha[j,] + Zi %*% gamma[j,] + res[j,] 
  ### data smoothing
  Y_list <- split(Y, id)
  for (i in 1:N){
    ss <- smooth.spline(time_rang, Y_list[[i]], cv = FALSE)
    y_hat[,i] <- predict(ss, timefine)$y
  }
  
  ### fit smoothed data at original time measurement points
  y <- convert_to_basisfunctions(t = timefine, eigenvecs = y_hat,
                                 tout = time_rang)
  
  ### FPCA
  fpcaobj <- prcomp(x=t(y), retx = TRUE, center = TRUE, rank. = 4)
  pcs <- fpcaobj$rotation
  
  ## Combine the smoothed data and basis functions to a dataframe
  response_test <- c(y)
  basis1 <- rep(pcs[,1], times = N)
  basis2 <- rep(pcs[,2], times = N)
  basis3 <- rep(pcs[,3], times = N)
  basis4 <- rep(pcs[,4], times = N)
  
  df_test <- data.frame(id = id, y= response_test, phi1 = basis1, 
                        phi2 = basis2, phi3 = basis3, phi4 = basis4) 
  
  ### Fit FMEM
  fform <- y ~ 1 + df_test$phi1 + df_test$phi2 + df_test$phi3 + df_test$phi4 +
    (-1 + df_test$phi1 + df_test$phi2 + df_test$phi3 + df_test$phi4 | df_test$id) + 
    (-1 + df_test$phi1 + df_test$phi2 + df_test$phi3 + df_test$phi4 | df_test$id)
  
  ft <-  fit_genetic_fmm(fform, df_test, A, pcs)
  
  ### Estimated covariance
  vc <- VarCorr(ft)
  CG <- vc[["df_test.id"]] # genetic covariance
  CE <- vc[["df_test.id.1"]] # environmental covariance
  
  CG_fun <- pcs %*% CG %*% t(pcs) # estimated gen cov function
  CE_fun <- pcs %*% CE %*% t(pcs) # estimated env cov function
  
  tr_gen4_20[j] <- tr(CG_fun)
  tr_env4_20[j] <- tr(CE_fun)
  
  ### Error on genetic and environmental covariance
  gen_err4_20[j] <- norm((C_gen_fun - CG_fun), type = "F")
  env_err4_20[j] <- norm((C_env_fun - CE_fun), type = "F")
}
)

save(tr_gen4_20,file = "tr_gen4_20.Rdata")
save(tr_env4_20,file = "tr_env4_20.Rdata")
save(gen_err4_20,file = "gen_err4_20.Rdata")
save(env_err4_20,file = "env_err4_20.Rdata")