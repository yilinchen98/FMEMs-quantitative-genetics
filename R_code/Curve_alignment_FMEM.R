# Combined analysis of time and phase variations in FMEMs

## load packages 
library(Matrix)
library(MASS)
library(fda)
library(fdasrvf)
library(lme4)
library(pedigreemm)
library(mvnfast)
library(ggplot2)
library(plotly)

## Import data and rescale to unit interval
setwd("D:/KCL_2023-2027_PhD/Genetics_FMEM_Project/FMEMs-quantitative-genetics/R_code")
TRFUN25PUP4 = read.delim("TRFUN25PUP4.DAT",header = FALSE)
names(TRFUN25PUP4)<-c("id","sire","dam","trait","x")
df <- data.frame(TRFUN25PUP4)

FirstUniqueIdPos <- which(duplicated(df$id) == FALSE)
N = length(FirstUniqueIdPos) # N = 873 subjects
n = length(df$id) # n = 6860 observations
age_list <- split(df$x,df$id)
trait_list <- split(df$trait,df$id)

age_list_new <- list()
for (i in 1:N){
  age_list_new[[i]] = (age_list[[i]]-min(age_list[[i]]))/(max(age_list[[i]])-min(age_list[[i]]))
}

df$x_rescaled <- unsplit(age_list_new,df$id)

## Data smoothing 
### Use positive smoothing and adjusting the smoothing parameters lambda <= 10^-4
agefine <- seq(0,1,length=100) # dense time grid
logmass <- matrix(0,100,N) # store smooth log growth curves
pred_mass <- matrix(0,100,N) # store the smoothed mass recovered by taking exp
lam <- rep(0,N) # smoothing parameter used for each log growth curve

for (i in 1:N){
  ss_logmass <- smooth.spline(age_list_new[[i]], log(trait_list[[i]]), cv=FALSE,
                              all.knots=TRUE)
  logmass[,i] <- predict(ss_logmass, agefine)$y
  pred_mass[,i] <- exp(predict(ss_logmass,agefine)$y)
  lam[i] <- ss_logmass$lambda
}

large_lam <- which(lam > 1e-4)
mass_adjusted <- matrix(0,100, length(large_lam))
for ( i in 1: length(large_lam)){
  ss_new <- smooth.spline(age_list_new[[large_lam[i]]], log(trait_list[[large_lam[i]]]), 
                          all.knots = TRUE, lambda = 1e-4) # restrict lambda <=1e-4
  mass_adjusted[,i] <- exp(predict(ss_new, agefine)$y)
}
colnames(mass_adjusted) <- as.character(large_lam)

## Let reform the smoothed growth curves into a matrix
trait_hat <- pred_mass 
for (i in 1: length(large_lam)){
  trait_hat[,large_lam[i]] <- mass_adjusted[,i]
}

## Curve alignment
### use fdasrvf to compute the aligned functions and warping functions
aligned_mass_process <- time_warping(f=trait_hat, time=agefine)
aligned_mass_curve <- aligned_mass_process$fn # aligned curves (amplitude function)
qn <- aligned_mass_process$qn # aligned curves in SRVF space
aligned_mean <- aligned_mass_process$fmean
warping_funs <- aligned_mass_process$warping_functions # warping function

## Compute phase function
jointFPCA_process <- jointFPCA(aligned_mass_process, no = 3)
C <- jointFPCA_process$C # scaling coefficient
psi <- SqrtMean(warping_funs)$psi # srvf of warping functions
v <- SqrtMean(warping_funs)$vec #phase function = shooting vector

# The eigen vector extracted from the function joinFPCA
# it joins the qn (srvf of original curves) and the warping functions

eigen_vec <- jointFPCA_process$U # eigen vector

par(mfrow = c(3,1))
for (i in 1:3) {
  plot(seq(0,1,length=201), eigen_vec[, i], type = "l", 
       xlab = "time", ylab = "",
       main = paste(" Eigen Vector", i))
}

g2 <- rbind(qn, C*v)

comb_pca <-  prcomp(x=t(g2), retx = TRUE, center = TRUE, rank. = 3)$rotation
par(mfrow=c(3,1))
for (i in 1:3) {
  plot(seq(0,1,length=200), comb_pca[, i], type = "l", 
       xlab = "time", ylab = "",
       main = paste(" Combined Principal Component", i))
}


### joint FPCA
#jointFPCA_process <- jointFPCA(aligned_mass_process, no = 3)
#C <- jointFPCA_process$C # scaling coefficient
g <- rbind(aligned_mass_curve, C*v) # joint aligned curves (amplitude functions) and corresponding phase function

Joint_fpcaobj <- prcomp(x=t(g), retx = TRUE, center = TRUE, rank. = 3)
joint_eig_vecs <- Joint_fpcaobj$rotation # eigen vectors
par(mfrow=c(3,1))
for (i in 1:3) {
  plot(seq(0,1,length=200), joint_eig_vecs[, i], type = "l", 
       xlab = "time", ylab = "",
       main = paste(" Joint Principal Component", i))
}

joint_eig_vecs[,1] %*% joint_eig_vecs[,2]
joint_eig_vecs[,1] %*% joint_eig_vecs[,3]
joint_eig_vecs[,3] %*% joint_eig_vecs[,2]

fpcaobj_mass <- prcomp(x=t(aligned_mass_curve), retx = TRUE, center = TRUE) 
eigen_mass <- fpcaobj_mass$rotation # eigen vectors

## Visualise the first three eigenfunctions
par(mfrow=c(3,1))
for (i in 1:3) {
  matplot(agefine, eigen_mass[, i], type = "l", 
          xlab = "time", ylab = "",
          main = paste("Growth Principal Component", i))
}

## Fit mixed-effect model

### Compute the genetic relationship matrix A
pos = df$id[FirstUniqueIdPos] # extract ids for all subjects
sire_id = df$sire[FirstUniqueIdPos] # extract ids for sire
dam_id = df$dam[FirstUniqueIdPos] # extract ids for dam

pede <- editPed(sire = sire_id, dam = dam_id, label = pos)
ped<- with(pede, pedigree(label=label, sire=sire, dam=dam))
A <- getA(ped)[163:1035,163:1035]

### interpolate smoothed data, corresponding warping functions and the joint PCs at original time points
trait_hat_list <- list() # smoothed data
for (i in 1:N){
  inter_trait <- convert_to_basisfunctions(t = agefine, eigenvecs = trait_hat[,i],
                                           tout = age_list_new[[i]])
  trait_hat_list[[i]] <- inter_trait
}

phasefun_list <- list() # shooting vectors (srvf transformed warping functions)
for(i in 1:N){
  inter_phase <- convert_to_basisfunctions(t=agefine, eigenvecs = C*v[,i], 
                                           tout = age_list_new[[i]])
  phasefun_list[[i]] <- inter_phase
}

basisfun_list <- list() # joint PCs
for (i in 1:N){
  PCs_curve <- convert_to_basisfunctions(t = agefine, eigenvecs = joint_eig_vecs[1:100,],
                                   tout = age_list_new[[i]])
  PCs_warping <- convert_to_basisfunctions(t = agefine, eigenvecs = joint_eig_vecs[101:200,],
                                           tout = age_list_new[[i]])
  basisfun_list[[i]] <- rbind(PCs_curve,PCs_warping)
}

### re-arrange the smoothed data, warping funtions and joint PCs into a dataframe
smoothed_curve <- list()
for (i in 1:N){
  curve_data <- c(trait_hat_list[[i]], phasefun_list[[i]]) # combined smoothed curve data and warping functions (shooting vector)
  smoothed_curve[[i]] <- curve_data
}

id_list <- list()
for (i in 1:N){
  id <- rep(pos[i], length(smoothed_curve[[i]])) # reform the id 
  id_list[[i]] <- id
}

subject_id <- unlist(id_list) 
trait <- unsplit(smoothed_curve, subject_id)
basis <- do.call(rbind,basisfun_list)
colnames(basis) <- c("phi1", "phi2", "phi3")

df_model <- data.frame(id = subject_id, trait = trait, basis) # this is new dataframe





