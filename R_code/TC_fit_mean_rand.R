library(fdasrvf)
library(Matrix)
library(pedigreemm)
library(lme4)

setwd("D:/KCL_2023-2027_PhD/Year 1/Genetics/FMEMs-quantitative-genetics/R_code")
# Import data
TRFUN25PUP4 <- read.delim("TRFUN25PUP4.DAT",header = FALSE)
names(TRFUN25PUP4)<-c("id","sire","dam","trait","x")
df_center <- data.frame(TRFUN25PUP4) 

FirstUniqueIdPos <- which(duplicated(df_center$id) == FALSE)
N = length(FirstUniqueIdPos) # N = 873 subjects
n = length(df_center$id) # n = 6860 observations
age_list <- split(df_center$x,df_center$id)
trait_list <- split(df_center$trait,df_center$id)

# Rescale time interval to [0,1]  
## x = (x -min(x))/(max(x) - min(x))
age_list_new <- list()
for (i in 1:N){
  age_list_new[[i]] = (age_list[[i]]-min(age_list[[i]]))/(max(age_list[[i]])-min(age_list[[i]]))
}

df_center$x <- unsplit(age_list_new,df_center$id)

## use smooth.spline to smooth the data
fitted_mass_list <- list() #Store fitted values for each curve in a list
lambda_i <- rep(0,N) # store smoothing parameters for each curve
for (i in 1:N){
  fit_mass <- smooth.spline(age_list_new[[i]], trait_list[[i]], cv=FALSE)
  # use GCV to choose smoothing parameters
  fitted_mass_list[[i]] <- fit_mass$y
  # update smoothing parameters
  lambda_i[i] <- fit_mass$lambda
}

## Now let's smooth the growth curve on log scale
fitted_logmass_list <- list()
lambda_logmass <- rep(0,N)
for (i in 1:N){
  fit_logmass <- smooth.spline(age_list_new[[i]],log10(trait_list[[i]]),cv=FALSE)
  fitted_logmass_list[[i]] <- fit_logmass$y
  lambda_logmass[i] <- fit_logmass$lambda
}

# Curve alignment
## Use align_fPCA from the package fdasrvf
### Transform the fitted mass to a fine time grid
agefine = seq(0,1,length=100)

### for the growth curve
mass_hat <- matrix(0,100,N) 
for (i in 1:N){
  fit_mass <- smooth.spline(age_list_new[[i]], trait_list[[i]], cv=FALSE)
  fit_mass_hat <- predict(fit_mass,agefine)
  mass_hat[,i] <- fit_mass_hat$y
}

### for log growth curve
logmass_hat <- matrix(0,100,N) 
for (i in 1:N){
  fit_logmass <- smooth.spline(age_list_new[[i]], log10(trait_list[[i]]), cv=FALSE)
  fit_logmass_hat <- predict(fit_logmass,agefine)
  logmass_hat[,i] <- fit_logmass_hat$y
}
# Based these plot, the prediction of mass/log(mass) based on smoothing spline works well.

## Use align_fPCA to align all curves
register_growth_curve <- align_fPCA(f = mass_hat, t = agefine, showplot = FALSE)
mass_aligned <- register_growth_curve$fn
mass_warp <- register_growth_curve$gam

par(mfrow=c(1,2))
matplot(agefine,mass_warp, type = "l", col=1:N, 
        xlab = "time", ylab="warping function",
        main = "Warping Functions")
matplot(agefine,mass_aligned, type = "l", col=1:N, 
        xlab = "time", ylab = "mass",
        main="Aligned Growth Curves")

## Compute the mean curve
mean_mass <- apply(mass_aligned, 1, mean)

## Plot the mean curve
lines(agefine, mean_mass, col="red", lwd=2)
legend(x="bottomright", legend="mean", lwd=2.0, col="red")

## Repeat the process for log(mass)
register_logrowth_curve <- align_fPCA(f = logmass_hat, t = agefine, showplot = FALSE)
logmass_aligned <- register_logrowth_curve$fn
logmass_warp <- register_logrowth_curve$gam

par(mfrow=c(1,2))
matplot(agefine,logmass_warp, type = "l", col=1:N, xlim = c(0,1),
        xlab = "time", ylab="warping function",
        main = "Warping Functions")
matplot(agefine,logmass_aligned, type = "l", col=1:N, xlim = c(0,1),
        xlab = "time", ylab = "mass",
        main="Aligned Log Growth Curves")

mean_logmass <- apply(logmass_aligned, 1, mean)
lines(agefine, mean_logmass, col="red", lwd=2)
legend(x="bottomright", legend="mean", lwd=2.0, col="red")

# Run FPCA on aligned data
## For the growth curve with zero mean
fpcaobj_mass <- prcomp(x=t(mass_aligned), retx = TRUE, center = TRUE) 
sdeval_mass <- fpcaobj_mass$sdev # the singular value of X
eigen_mass <- fpcaobj_mass$rotation # eigen vectors
pcscores_mass <- fpcaobj_mass$x # principal component scores

## For the log growth curve
fpcaobj_logmass <- prcomp(x=t(logmass_aligned), retx = TRUE, center = TRUE) 
sdeval_logmass <- fpcaobj_logmass$sdev # the singular value of X
eigen_logmass <- fpcaobj_logmass$rotation # eigen vectors
pcscores_logmass <- fpcaobj_logmass$x # principal component scores

# Fit FMEMs
## Centered the data
mean_mass_list <- list()
for (i in 1:N){
  mean_mass_curve <- convert_to_basisfunctions(agefine, mean_mass, age_list_new[[i]])
  mean_mass_list[[i]] <- mean_mass_curve
}
mean_mass_curve <- unsplit(mean_mass_list, df_center$id)

## convert basis vectors to basis functions
p =3 # number of basis

### define the basis functions (eigenfunctions)
phi_list <- list() 
# create an empty list which stores eigenfunctions for 873 subjects
# evaluated at the original time points.

for (i in 1:N){
  phi <- convert_to_basisfunctions(t = agefine, eigenvecs = eigen_mass[,1:p],
                                   tout = age_list_new[[i]], method = "linear")
  phi_list[[i]] <- phi
}

phi <- do.call(rbind,phi_list)
colnames(phi) <- c("phi1", "phi2", "phi3")


## Update the data frame 
df_center <- cbind(df_center,phi)
df_center$trait <- df_center$trait - mean_mass_curve

## calculate the pedigree and relationship matrix
pos = df_center$id[FirstUniqueIdPos] # extract ids for all subjects
sire_id = df_center$sire[FirstUniqueIdPos] # extract ids for sire
dam_id = df_center$dam[FirstUniqueIdPos] # extract ids for dam

pede <- editPed(sire = sire_id, dam = dam_id, label = pos)
ped<- with(pede, pedigree(label=label, sire=sire, dam=dam))
A <- getA(ped)[163:1035,163:1035]

## Fit mixed-effect model
mrformula <- trait ~ -1 + (-1 + df_center$phi1 + df_center$phi2 + df_center$phi3 | df_center$id) + 
  (-1 + df_center$phi1 + df_center$phi2 + df_center$phi3 | df_center$id)
system.time(
fit_mean_rand <- fit_genetic_fmm(mrformula, df_center, A, phi)
) # user 57.75 system 3.10 elapsed 71.97
