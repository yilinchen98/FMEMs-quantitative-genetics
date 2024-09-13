# Bootstrap TC data
## Use log scale

.libPaths("/scratch/users/k2368651/software/R/4.3")
## load packages 
library(Matrix)
library(MASS)
library(fda)
library(lme4)
library(fdasrvf)
library(pedigreemm)
library(mvnfast)
library(ggplot2)
library(fBasics)

## load functions
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
  #' @param nbasis numerber: number of basis used to fit the model
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

### load data
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

## Rescale time interval to [0,1]
N = length(FirstUniqueIdPos) # N = 873 subjects
n = length(df$id) # n = 6860 observations
age_list <- split(df$x,df$id)
trait_list <- split(df$trait,df$id)

### x = (x -min(x))/(max(x) - min(x))
age_list_new <- list()
for (i in 1:N){
  age_list_new[[i]] = (age_list[[i]]-min(age_list[[i]]))/(max(age_list[[i]])-min(age_list[[i]]))
}

df$x_rescaled <- unsplit(age_list_new,df$id)
df$logx <- log10(df$trait)

### Data smoothing
timefine <- seq(0,1,length=100) # dense time grid
mass_smoothed <-list()
pred_mass_fine <- matrix(0,100,N) # store the smoothed mass predicted on the dense grid
pred_logmass_fine <- matrix(0,100,N) # store the smoothed logmass
lam <- rep(0,N) # smoothing parameter used for each log growth curve

for (i in 1:N) {
  ss_logmass <- smooth.spline(age_list_new[[i]], log10(trait_list[[i]]), cv=FALSE, all.knots=TRUE)
  
  # Check if lambda is greater than 1e-4
  if (ss_logmass$lambda > 1e-4) {
    # Redo smoothing with lambda set to 1e-4
    ss_logmass <- smooth.spline(age_list_new[[i]], log10(trait_list[[i]]), lambda=1e-4, all.knots=TRUE)
  }
  
  mass_smoothed[[i]] <- 10^(ss_logmass$y)
  pred_logmass_fine[,i] <- predict(ss_logmass, timefine)$y
  pred_mass_fine[,i] <- 10^(predict(ss_logmass, timefine)$y) ## predict smoothed data on a regular dense grid
  lam[i] <- ss_logmass$lambda
}

### curve alignment
aligned_logmass_process <- time_warping(pred_logmass_fine, timefine)
aligned_logmass_curve <- aligned_logmass_process$fn
aligned_logmass_mean <- aligned_logmass_process$fmean
warping_logmass_funs <- aligned_logmass_process$warping_functions

### FPCA
fpcaobj_logmass <- prcomp(x=t(aligned_logmass_curve), retx = TRUE, center = TRUE, rank. = 3)
pcs_logmass <- fpcaobj_logmass$rotation # eigen vectors

### Model Fitting (Log scale)

### Align raw data
aligned_logtrait <- list()
gamma_logmass <- list()
for (i in 1:N){
  gamma_logmass_inter <- convert_to_basisfunctions(timefine, warping_logmass_funs[,i], age_list_new[[i]])
  gamma_logmass[[i]] <- gamma_logmass_inter
  aligned_logtrait[[i]] <- warp_f_gamma(log10(trait_list[[i]]), age_list_new[[i]], gamma_logmass_inter)
}

### basis functions
phi_logmass_list <- list() 
# create an empty list which stores eigenfunctions for 873 subjects
# evaluated at the original time points.

for (i in 1:N){
  phi_logmass <- convert_to_basisfunctions(t = timefine, eigenvecs = pcs_logmass,
                                           tout = age_list_new[[i]])
  phi_logmass_list[[i]] <- phi_logmass
}

phi_logmass <- do.call(rbind,phi_logmass_list)
colnames(phi_logmass) <- c("phi1", "phi2", "phi3")

### Reform dataframe
df_logmass <- data.frame(id = df$id, trait = unsplit(aligned_logtrait,df$id), phi_logmass)

fmmFormL <- trait ~ -1 + df_logmass$phi1 + df_logmass$phi2 + df_logmass$phi3 +
  (-1 + df_logmass$phi1 + df_logmass$phi2 | df_logmass$id) + 
  (-1 + df_logmass$phi1 + df_logmass$phi2 | df_logmass$id) 

ffL <- fit_genetic_fmm(fmmFormL, df_logmass, A, 2)
isSingular(ffL)

### Extract model results
betaL_hat <- fixef(ffL)
fefL_hat <- pcs_logmass %*% betaL_hat # fixed-effect
save(betaL_hat, file = "betaL_hat.Rdata")
save(fefL_hat, file = "fefL_hat.Rdata")

vcL <- VarCorr(ffL)
CGL <- vcL[["df_logmass.id"]] # genetic covariance
CEL <- vcL[["df_logmass.id.1"]] # environmental covariance
sig2 <- sigma(ffL)^2 # residual variance

CG_funL_hat <- pcs_logmass[,1:2] %*% CGL %*% t(pcs_logmass[,1:2]) # estimated gen cov function
CE_funL_hat <- pcs_logmass[,1:2] %*% CEL %*% t(pcs_logmass[,1:2]) # estimated env cov function
save(CG_funL_hat, file = "CG_fun_hatL.Rdata")
save(CE_funL_hat, file = "CE_fun_hatL.Rdata")

## Bootstrapping

### model form
id <- df$id
gpf <- factor(id) # grouping factor
Ji <- t(as(gpf, Class = "sparseMatrix")) # Indicator matrix of grouping factor indices

Xi <- as(phi_logmass[,1:2], Class = "sparseMatrix") #raw random effect design matrix
Zi <- t(KhatriRao(t(Ji), t(Xi))) # random effect design matrix

### reparameterise covariance structure
I_N <- as(diag(N), "sparseMatrix")
I_n <- as(diag(n), "sparseMatrix")

CG_para <- as(kronecker(A, CGL), "sparseMatrix") # reparameterised genetic covariance 
CE_para <- as(kronecker(I_N, CEL), "sparseMatrix") # reparameterised environmental covariance 

mu_hat <- rep(0, dim(CG_para)[1])
mu_res_hat <- rep(0, n)
res_cov_hat <- as(sig2 * I_n, "sparseMatrix") 

set.seed(123)
### generate bootstrap samples
n_iter <- 800 # number of iterations
npc <- 3

Y_estL <- matrix(0, n, n_iter) # store bootstrap sample
lam_sapL <- matrix(0, N, n_iter) # store the smoothing parameter
alignedYL <- array(data = 0, c(100,N,n_iter))
pcs_sapL <- array(data = 0, c(100,npc,n_iter)) # store the PCs
fef_sapL <- matrix(0,100,n_iter) # store estimated fixed-effect
CG_fun_sapL <- array(data = 0, c(100,100,n_iter)) # 3d matrix store genetic covariance function of 500 iterations
CE_fun_sapL <- array(data = 0, c(100,100,n_iter)) # 3d matrix store environmental covariance function of 500 iterations
res_sapL <- rep(0,n_iter) # store estimated residual
sfit_sapL <- rep(NA, n_iter)


for (k in 1:n_iter){
  alpha_hat <- rmvn(n=1, mu = mu_hat, sigma = CG_para) # genetic random effect
  gamma_hat <- rmvn(n=1, mu = mu_hat, sigma = CE_para) # environmental random effect
  res_hat <- rmvn(n=1, mu=mu_res_hat, sigma = res_cov_hat)# error vector
  
  Yk <- phi_logmass %*% betaL_hat + Zi %*% t(alpha_hat) + Zi %*% t(gamma_hat) + t(res_hat) # generate bootstrap sample (aligned)
  Y_estL[,k] <- Yk@x
  save(Y_estL, file = "bs_TC_log.Rdata")
  
  Yk_list <- split(Yk, id)
  
  ## data smoothing
  pred_logmassk_fine <- matrix(0,100,N) # store the smoothed logmass
  lamk <- rep(0,N) # smoothing parameter used for each log growth curve
  
  for (i in 1:N) {
    ssk_logmass <- smooth.spline(age_list_new[[i]], Yk_list[[i]], cv=FALSE, all.knots=TRUE)
    
    # Check if lambda is greater than 1e-4
    if (ssk_logmass$lambda > 1e-4) {
      # Redo smoothing with lambda set to 1e-4
      ssk_logmass <- smooth.spline(age_list_new[[i]], Yk_list[[i]], lambda=1e-4, all.knots=TRUE)
    }
    
    pred_logmassk_fine[,i] <- predict(ssk_logmass, timefine)$y
    lamk[i] <- ssk_logmass$lambda
  }
  lam_sapL[,k] <- lamk
  save(lam_sapL, file = "lam_sapL.Rdata")
  
  ## FPCA
  fpcaobj_logmassk <- prcomp(x=t(pred_logmassk_fine), retx = TRUE, center = TRUE, rank. = 3)
  pcs_logmassk <- fpcaobj_logmassk$rotation # eigen vectors
  
  pcs_sapL[,,k] <- pcs_logmassk
  save(pcs_sapL, file = "pcs_sapL.Rdata")
  
  ### basis functions
  phi_logmassk_list <- list() 
  # create an empty list which stores eigenfunctions for 873 subjects
  # evaluated at the original time points.
  
  for (i in 1:N){
    phi_logmassk <- convert_to_basisfunctions(t = timefine, eigenvecs = pcs_logmassk,
                                             tout = age_list_new[[i]])
    phi_logmassk_list[[i]] <- phi_logmassk
  }
  
  phi_logmassk <- do.call(rbind,phi_logmassk_list)
  colnames(phi_logmassk) <- c("phi1", "phi2", "phi3")
  
  ### Reform dataframe
  df_logmassk <- data.frame(id = id, trait = Yk@x, phi_logmassk)
  
  fmmFormLk <- trait ~ -1 + df_logmassk$phi1 + df_logmassk$phi2 + df_logmassk$phi3 + 
    (-1 + df_logmassk$phi1 + df_logmassk$phi2 | df_logmassk$id) + 
    (-1 + df_logmassk$phi1 + df_logmassk$phi2 | df_logmassk$id) 
  
  ffLk <- fit_genetic_fmm(fmmFormLk, df_logmassk, A, 2)
  
  sfit_sapL[k] <- isSingular(ffLk)
  save(sfit_sapL, file = "sfit_sapL.Rdata")
  
  ### Extract model results
  betaLk <- fixef(ffLk)
  fef_sapL[,k] <- pcs_logmassk%*% betaLk # fixed-effect
  save(fef_sapL, file = "fef_sapL.Rdata")
  
  vcLk <- VarCorr(ffLk)
  CGLk <- vcLk[["df_logmassk.id"]] # genetic covariance
  CELk <- vcLk[["df_logmassk.id.1"]] # environmental covariance
  res_sapL[k] <- sigma(ffLk)^2 # residual variance
  save(res_sapL, file = "res_sapL.Rdata")
  
  CG_funLk <- pcs_logmassk[,1:2] %*% CGLk %*% t(pcs_logmassk[,1:2]) # estimated gen cov function
  CE_funLk <- pcs_logmassk[,1:2] %*% CELk %*% t(pcs_logmassk[,1:2]) # estimated env cov function
  CG_fun_sapL[,,k] <- CG_funLk
  CE_fun_sapL[,,k] <- CE_funLk
  save(CG_fun_sapL, file = "CG_fun_sapL.Rdata")
  save(CE_fun_sapL, file = "CE_fun_sapL.Rdata")
}
