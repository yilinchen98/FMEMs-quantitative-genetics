# Bootstrap TC data
## Use log scale

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
aligned_mass_process <- time_warping(f=pred_mass_fine, time=timefine)
aligned_mass_curve <- aligned_mass_process$fn
aligned_mean <- aligned_mass_process$fmean
warping_funs <- aligned_mass_process$warping_functions

### FPCA
fpcaobj_mass <- prcomp(x=t(aligned_mass_curve), retx = TRUE, center = TRUE, rank. = 4)
pcs_mass <- fpcaobj_mass$rotation # eigen vectors

### Model Fitting (Log scale)

### Align raw data
aligned_trait <- list()
gamma_mass <- list()
for (i in 1:N){
  gamma_mass_inter <- convert_to_basisfunctions(timefine, warping_funs[,i], age_list_new[[i]])
  gamma_mass[[i]] <- gamma_mass_inter
  aligned_trait[[i]] <- warp_f_gamma(trait_list[[i]], age_list_new[[i]], gamma_mass_inter)
}

### basis functions
phi_list <- list() 
# create an empty list which stores eigenfunctions for 873 subjects
# evaluated at the original time points.

for (i in 1:N){
  phi <- convert_to_basisfunctions(t = timefine, eigenvecs = pcs_mass,
                                   tout = age_list_new[[i]])
  phi_list[[i]] <- phi
}

phi <- do.call(rbind,phi_list)
colnames(phi) <- c("phi1", "phi2", "phi3", "phi4")

## Reform dataframe
df_mass <- data.frame(id = df$id, trait = unsplit(aligned_trait,df$id), phi)

fmmForm <- trait ~ -1 + df_mass$phi1 + df_mass$phi2 + 
  (-1 + df_mass$phi1 + df_mass$phi2 | df_mass$id) + 
  (-1 + df_mass$phi1 + df_mass$phi2 | df_mass$id) 

ff <- fit_genetic_fmm(fmmForm, df_mass, A, 2)
isSingular(ff)

### Extract model results
## fixed effect
beta_hat <- fixef(ff)
fef_hat <- pcs_mass[,1:2] %*% beta_hat
save(beta_hat, file = "beta_hat.Rdata")
save(fef_hat, file = "fef_hat.Rdata")

## Extract covariance
vc <- VarCorr(ff)
CG <- vc[["df_mass.id"]] # genetic covariance
CE <- vc[["df_mass.id.1"]] # environmental covariance
sig2 <- sigma(ff)^2 # residual variance

CG_fun_hat <- pcs_mass[,1:2] %*% CG %*% t(pcs_mass[,1:2]) # estimated gen cov function
CE_fun_hat <- pcs_mass[,1:2] %*% CE %*% t(pcs_mass[,1:2]) # estimated env cov function
save(CG_fun_hat, file = "CG_fun_hat.Rdata")
save(CE_fun_hat, file = "CE_fun_hat.Rdata")

## Bootstrapping

### model form
id <- df$id
gpf <- factor(id) # grouping factor
Ji <- t(as(gpf, Class = "sparseMatrix")) # Indicator matrix of grouping factor indices

Xi <- as(phi[,1:2], Class = "sparseMatrix") #raw random effect design matrix
Zi <- t(KhatriRao(t(Ji), t(Xi))) # random effect design matrix

### reparameterise covariance structure
I_N <- as(diag(N), "sparseMatrix")
I_n <- as(diag(n), "sparseMatrix")

CG_para <- as(kronecker(A, CG), "sparseMatrix") # reparameterised genetic covariance 
CE_para <- as(kronecker(I_N, CE), "sparseMatrix") # reparameterised environmental covariance 

mu_hat <- rep(0, dim(CG_para)[1])
mu_res_hat <- rep(0, n)
res_cov_hat <- as(sig2 * I_n, "sparseMatrix") 

set.seed(123)
### generate bootstrap samples
n_iter <- 800 # number of iterations
npc <- 4

Y_est <- matrix(0, n, n_iter) # store bootstrap sample
lam_sap <- matrix(0, N, n_iter) # store the smoothing parameter
pcs_sap <- array(data = 0, c(100,npc,n_iter)) # store the PCs
fef_sap <- matrix(0,100,n_iter) # store estimated fixed-effect
CG_fun_sap <- array(data = 0, c(100,100,n_iter)) # 3d matrix store genetic covariance function of 500 iterations
CE_fun_sap <- array(data = 0, c(100,100,n_iter)) # 3d matrix store environmental covariance function of 500 iterations
res_sap <- rep(0,n_iter) # store estimated residual
sfit_sap <- rep(NA, n_iter)


for (k in 1:n_iter){
  alpha_hat <- rmvn(n=1, mu = mu_hat, sigma = CG_para) # genetic random effect
  gamma_hat <- rmvn(n=1, mu = mu_hat, sigma = CE_para) # environmental random effect
  res_hat <- rmvn(n=1, mu=mu_res_hat, sigma = res_cov_hat)# error vector
  
  Yk <- phi[,1:2] %*% beta_hat + Zi %*% t(alpha_hat) + Zi %*% t(gamma_hat) + t(res_hat) # generate bootstrap sample (aligned)
  Y_est[,k] <- Yk@x
  save(Y_est, file = "bs_TC.Rdata")
  
  Yk_list <- split(Yk, id)
  
  ## data smoothing
  pred_massk_fine <- matrix(0,100,N) # store the smoothed logmass
  lamk <- rep(0,N) # smoothing parameter used for each log growth curve
  
  for (i in 1:N) {
    ssk_mass <- smooth.spline(age_list_new[[i]], Yk_list[[i]], cv=FALSE, all.knots=TRUE)
    
    # Check if lambda is greater than 1e-4
    if (ssk_mass$lambda > 1e-4) {
      # Redo smoothing with lambda set to 1e-4
      ssk_mass <- smooth.spline(age_list_new[[i]], Yk_list[[i]], lambda=1e-4, all.knots=TRUE)
    }
    
    pred_massk_fine[,i] <- predict(ssk_mass, timefine)$y
    lamk[i] <- ssk_mass$lambda
  }
  lam_sap[,k] <- lamk
  save(lam_sap, file = "lam_sap.Rdata")
  
  ## FPCA
  fpcaobj_massk <- prcomp(x=t(pred_massk_fine), retx = TRUE, center = TRUE, rank. = 4)
  pcs_massk <- fpcaobj_massk$rotation # eigen vectors
  
  pcs_sap[,,k] <- pcs_massk
  save(pcs_sap, file = "pcs_sap.Rdata")
  
  ### basis functions
  phi_massk_list <- list() 
  # create an empty list which stores eigenfunctions for 873 subjects
  # evaluated at the original time points.
  
  for (i in 1:N){
    phi_massk <- convert_to_basisfunctions(t = timefine, eigenvecs = pcs_massk,
                                             tout = age_list_new[[i]])
    phi_massk_list[[i]] <- phi_massk
  }
  
  phi_massk <- do.call(rbind,phi_massk_list)
  colnames(phi_massk) <- c("phi1", "phi2", "phi3", "phi4")
  
  ### Reform dataframe
  df_massk <- data.frame(id = id, trait = Yk@x, phi_massk)
  
  fmmFormk <- trait ~ -1 + df_massk$phi1 + df_massk$phi2 + 
    (-1 + df_massk$phi1 + df_massk$phi2 | df_massk$id) + 
    (-1 + df_massk$phi1 + df_massk$phi2 | df_massk$id) 
  
  ffk <- fit_genetic_fmm(fmmFormk, df_massk, A, 2)
  
  sfit_sap[k] <- isSingular(ffk)
  save(sfit_sap, file = "sfit_sap.Rdata")
  
  ### Extract model results
  betak <- fixef(ffk)
  fef_sap[,k] <- pcs_massk[,1:2] %*% betak # fixed-effect
  save(fef_sap, file = "fef_sap.Rdata")
  
  vck <- VarCorr(ffk)
  CGk <- vck[["df_massk.id"]] # genetic covariance
  CEk <- vck[["df_massk.id.1"]] # environmental covariance
  res_sap[k] <- sigma(ffk)^2 # residual variance
  save(res_sap, file = "res_sap.Rdata")
  
  CG_funk <- pcs_massk[,1:2] %*% CGk %*% t(pcs_massk[,1:2]) # estimated gen cov function
  CE_funk <- pcs_massk[,1:2] %*% CEk %*% t(pcs_massk[,1:2]) # estimated env cov function
  CG_fun_sap[,,k] <- CG_funk
  CE_fun_sap[,,k] <- CE_funk
  save(CG_fun_sap, file = "CG_fun_sap.Rdata")
  save(CE_fun_sap, file = "CE_fun_sap.Rdata")
}


