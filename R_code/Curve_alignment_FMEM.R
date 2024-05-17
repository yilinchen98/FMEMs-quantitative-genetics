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
library(GeodRegr)

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
aligned_mass_curve <- aligned_mass_process$fn # aligned function (amplitude function)
aligned_mean <- aligned_mass_process$fmean
warping_funs <- aligned_mass_process$warping_functions # warping function

psi <- SqrtMean(warping_funs)$psi # square root of warping functions
gam_mu <- SqrtMean(warping_funs)$gam_mu # center of warping functions Karcher mean 

mu <- SqrtMean(warping_funs)$mu # Karcher mean of psi
plot(agefine, gam_mu, type = "l") # identity map
grid(nx = NULL, ny = NULL,lty = 2,col = "gray",lwd = 2)

## Compute phase function
x <- gam_to_v(warping_funs) #map warping function to tangent space at identity

### jointFPCA

jointFPCA_process <- jointFPCA(aligned_mass_process, no = 3)
eigen_vec <- jointFPCA_process$U
C <- jointFPCA_process$C

par(mfrow = c(3,1))
plot(seq(0,1, length=201), eigen_vec[,1], type = "l")
plot(seq(0,1, length=201), eigen_vec[,2], type = "l")
plot(seq(0,1, length=201), eigen_vec[,3], type = "l")

eigen_vec[,2] %*% eigen_vec[,3]

############################
g <- rbind(aligned_mass_curve, C*x)

fpcaobj_comb <- prcomp(x=t(g), retx = TRUE, center = TRUE, rank. = 3)
eigen_comb <- fpcaobj_comb$rotation # eigen vectors

par(mfrow=c(3,1))
for (i in 1:3) {
  plot(seq(0,1,length=200), eigen_comb[, i], type = "l", 
       xlab = "time", ylab = "",
       main = paste(" Combined Principal Component", i))
}

eigen_comb[,1] %*% eigen_comb[,2]

## Fit mixed-effect model

pos = df$id[FirstUniqueIdPos] # extract ids for all subjects
sire_id = df$sire[FirstUniqueIdPos] # extract ids for sire
dam_id = df$dam[FirstUniqueIdPos] # extract ids for dam

pede <- editPed(sire = sire_id, dam = dam_id, label = pos)
ped<- with(pede, pedigree(label=label, sire=sire, dam=dam))
A <- getA(ped)[163:1035,163:1035]

## fit smoothed data at original time points
trait_hat_list <- list()

for (i in 1:N){
  inter_trait <- convert_to_basisfunctions(t = agefine, eigenvecs = trait_hat[,i],
                                           tout = age_list_new[[i]])
  trait_hat_list[[i]] <- inter_trait
}

phasefun_list <- list()

for(i in 1:N){
  inter_phase <- convert_to_basisfunctions(t=agefine, eigenvecs = C*x[,i], 
                                           tout = age_list_new[[i]])
  phasefun_list[[i]] <- inter_phase
}




