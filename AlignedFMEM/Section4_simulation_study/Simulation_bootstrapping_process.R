source("PrepareR.R")
# This is an example of fit aligned FMEM to one of the simulated dataset and 
# then carry on bootstrapping.

group_id <- 1

N <- 873 # total number of individuals
n <- 14 # number of measurements per individual
nbasis <- 5 # number of basis
npc <- 6
time_rang <- seq(0,1,length=n)

A <- readRDS("sim_results/A.rds")
Y <- readRDS("sim_results/sim_data.rds")

## Model fitting
uniqueIds <- seq(1,N, length=N)
id <- rep(uniqueIds, each=n)

yk_hat <- matrix(0,n,N) # directly smooth data
lam <- rep(0, N)

Yk <- split(Y[,group_id],id)
for (i in 1:N) {
  ss <- smooth.spline(time_rang, Yk[[i]], cv = FALSE)
  lam[i] <- ss$lambda
  yk_hat[,i] <- ss$y
}

## FPCA
fpcaobj <- prcomp(t(yk_hat), center = TRUE, retx = TRUE, rank. = npc)
fpcs <- fpcaobj$rotation

b1 <- rep(fpcs[,1], times = N) 
b2 <- rep(fpcs[,2], times = N) 
b3 <- rep(fpcs[,3], times = N) 
b4 <- rep(fpcs[,4], times = N) 
b5 <- rep(fpcs[,5], times = N) 
b6 <- rep(fpcs[,6], times = N) 

df_sm <- data.frame(id = id, y= Y[,group_id], 
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
fef_hat <- fpcs %*% beta

vc <- VarCorr(ft)
CG_hat <- vc[["df_sm.id"]] # genetic covariance
CE_hat <- vc[["df_sm.id.1"]] # environmental covariance
sig2_hat <- sigma(ft)^2 # residual variance

CG_fun <- fpcs[,1:5] %*% CG_hat %*% t(fpcs[,1:5]) # estimated gen cov function
CE_fun <- fpcs[,1:5] %*% CE_hat %*% t(fpcs[,1:5]) # estimated env cov function

### Check singular fit
sfit <- isSingular(ft)

## Generate bootstrap samples
set.seed(1)
### design matrices
gpf <- gl(N, n)
Ji  <- t(as(gpf, Class = "sparseMatrix"))
Xi  <- as(cbind(b1,b2,b3,b4,b5), Class = "sparseMatrix")
Zi  <- t(KhatriRao(t(Ji), t(Xi)))

### covariance
IN <- as(diag(N), "sparseMatrix")
CG_para <- as(kronecker(A,  CG_hat), "sparseMatrix")
CE_para <- as(kronecker(IN, CE_hat), "sparseMatrix")

## Bootstrapping
n_iter <- 300 # total nubmer of bootstrap samples generated
for(j in 1:n_iter){
  alpha_hat <- rmvn(n = 1, mu = rep(0, nrow(CG_para)), sigma = CG_para)
  gamma_hat <- rmvn(n = 1, mu = rep(0, nrow(CE_para)), sigma = CE_para)
  res_hat   <- t(rmvn(n = N, mu = rep(0, n), sigma = sig2_hat * diag(n)))
  
  yk <- rep(fef_hat, N) + Zi %*% t(alpha_hat) + Zi %*% t(gamma_hat) + c(res_hat)
  y_bs <- matrix(yk@x, nrow = n, ncol = N)
  
  ## fit the model
  ### 1.smoothing
  yk_hat <- matrix(0, n, N)
  lam_sap <- rep(0, N)
  for (i in 1:N) {
    ss <- smooth.spline(time_rang, y_bs[, i], cv = FALSE)
    lam_sap[i] <- ss$lambda
    yk_hat[, i] <- ss$y
  }
  
  ### 2.FPCA
  fpcaobjk <- prcomp(t(yk_hat), center = TRUE, retx = TRUE, rank. = npc)
  pcs <- fpcaobjk$rotation
  
  ### 3.reform model framework
  bs1 <- rep(pcs[,1], times = N)
  bs2 <- rep(pcs[,2], times = N)
  bs3 <- rep(pcs[,3], times = N)
  bs4 <- rep(pcs[,4], times = N)
  bs5 <- rep(pcs[,5], times = N)
  bs6 <- rep(pcs[,6], times = N)
  
  df_sm <- data.frame(
    id = id, y = yk@x,
    phi1 = bs1, phi2 = bs2, phi3 = bs3,
    phi4 = bs4, phi5 = bs5, phi6 = bs6
  )
  
  fform <- y ~ -1 + df_sm$phi1 + df_sm$phi2 + df_sm$phi3 + df_sm$phi4 + df_sm$phi5 + df_sm$phi6 +
    (-1 + df_sm$phi1 + df_sm$phi2 + df_sm$phi3 + df_sm$phi4 + df_sm$phi5 | df_sm$id) +
    (-1 + df_sm$phi1 + df_sm$phi2 + df_sm$phi3 + df_sm$phi4 + df_sm$phi5 | df_sm$id)
  
  ### 4.fit the model and extract results
  ft <- fit_genetic_fmm(fform, df_sm, A, nbasis)
  
  beta_sap <- fixef(ft)
  fef_sap <- pcs %*% beta_sap
  
  vc <- VarCorr(ft)
  CG_sap <- vc[["df_sm.id"]] 
  CE_sap <- vc[["df_sm.id.1"]]
  sig2_sap <- sigma(ft)^2
  
  CG_fun_sap <- pcs[,1:5] %*% CG_sap %*% t(pcs[,1:5]) # estimated gen cov function
  CE_fun_sap <- pcs[,1:5] %*% CE_sap %*% t(pcs[,1:5]) # estimated env cov function
  
  sfit_sap <- isSingular(ft)
  
  bs_res <- list(y_bs = y_bs,
                 lam_sap = lam_sap,
                 fef_sap = fef_sap,
                 CG_sap = CG_sap,
                 CE_sap = CE_sap,
                 sig2_sap = sig2_sap,
                 CG_fun_sap = CG_fun_sap,
                 CE_fun_sap = CE_fun_sap,
                 sfit_sap = sfit_sap)
}
