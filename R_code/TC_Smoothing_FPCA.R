library(fdasrvf)
library(pedigreemm)
library(lme4)
library(Matrix)

# Import data
setwd("D:/KCL_2023-2027_PhD/Year 1/Genetics/FMEMs-quantitative-genetics/R_code")
TRFUN25PUP4 = read.delim("TRFUN25PUP4.DAT",header = FALSE)
names(TRFUN25PUP4)<-c("id","sire","dam","trait","x")
df <- data.frame(TRFUN25PUP4) # n = 6860 observations

# Plot raw data
FirstUniqueIdPos <- which(duplicated(df$id) == FALSE)
N = length(FirstUniqueIdPos) # N = 873 subjects
n = length(df$id) # n = 6860 observations
age_list <- split(df$x,df$id)
trait_list <- split(df$trait,df$id)

# check missing data
which(is.na(df$x))
which(is.na(df$trait)) # no missing data

## Plot the Growth curve on raw data
par(mfrow=c(2,1))
plot(age_list[[1]],trait_list[[1]],type = 'l',xlab = 'Days', 
     ylab = 'Mass', xlim = c(1,25), ylim = c(0,400),
     main = 'Growth of Tribolium Castaneum')
for (i in 2:N) {
  lines(age_list[[i]], trait_list[[i]], col = i)
}

## Plot the Growth curve on log scale
plot(age_list[[1]],log10(trait_list[[1]]),type = 'l',xlab = 'Days', 
     ylab = 'log(Mass)', xlim = c(1,25), ylim = c(0,3),
     main = 'Log Growth Curve of Tribolium Castaneum')
for (i in 2:N) {
  lines(age_list[[i]], log10(trait_list[[i]]), col = i)
}

# Rescale time interval to [0,1]
## x = (x -min(x))/(max(x) - min(x))
age_list_new <- list()
for (i in 1:N){
  age_list_new[[i]] = (age_list[[i]]-min(age_list[[i]]))/(max(age_list[[i]])-min(age_list[[i]]))
}
## check the first subject
print(age_list_new[[1]])

## use smooth.spline to smooth the data
fitted_mass_list <- list() #Store fitted values for each curve in a list
lambda_i <- rep(0,N) # store smoothing parameters for each curve
### Plot the smoothed curves
plot(c(0,1), c(0,400), type="n", xlab = "Days", ylab="Mass",
     main="Smoothed Growth Plot of Tribolium Castaneum")
for (i in 1:N){
  fit_mass <- smooth.spline(age_list_new[[i]], trait_list[[i]], cv=FALSE)
  # use GCV to choose smoothing parameters
  lines(fit_mass, col=i)
  # update the fitted values
  fitted_mass_list[[i]] <- fit_mass$y
  # update smoothing parameters
  lambda_i[i] <- fit_mass$lambda
}

## Now let's smooth the growth curve on log scale
fitted_logmass_list <- list()
lambda_logmass <- rep(0,N)
plot(c(0,1),c(0,3),type="n",xlab="Days",ylab="log(Mass)",
     main="Smoothed Log Growth of Tribolium Castaneum")
for (i in 1:N){
  fit_logmass <- smooth.spline(age_list_new[[i]],log10(trait_list[[i]]),cv=FALSE)
  lines(fit_logmass,col=i)
  fitted_logmass_list[[i]] <- fit_logmass$y
  lambda_logmass[i] <- fit_logmass$lambda
}

# Curve alignment
## Use align_fPCA from the package fdasrvf
### Transform the fitted mass to a fine time grid
agefine = seq(0,1,length=100)

### for the growth curve
mass_hat <- matrix(0,100,N) 
plot(c(0,1),c(0,400), type="n", xlab="Days",ylab="Mass",
     main = "Growth plot on a fine grid")
for (i in 1:N){
  fit_mass <- smooth.spline(age_list_new[[i]], trait_list[[i]], cv=FALSE)
  fit_mass_hat <- predict(fit_mass,agefine)
  lines( fit_mass_hat,col=i)
  mass_hat[,i] <- fit_mass_hat$y
}

### for log growth curve
logmass_hat <- matrix(0,100,N) 
plot(c(0,1),c(0,3), type="n", xlab="Days",ylab="log(Mass)",
     main = "Log growth plot on a fine grid")
for (i in 1:N){
  fit_logmass <- smooth.spline(age_list_new[[i]], log10(trait_list[[i]]), cv=FALSE)
  fit_logmass_hat <- predict(fit_logmass,agefine)
  lines( fit_logmass_hat,col=i)
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

# Run FPCA on the aligned curves

## For the growth curve with zero mean
fpcaobj_mass <- prcomp(x=t(mass_aligned), retx = TRUE, center = TRUE) 
sdeval_mass <- fpcaobj_mass$sdev # the singular value of X
eigen_mass <- fpcaobj_mass$rotation # eigen vectors
pcscores_mass <- fpcaobj_mass$x # principal component scores
summary(fpcaobj_mass)
# Cumulative Proportion:  0.9475  0.99017  0.99808 (up to PC3)

## Visualise the first three eigenfunctions
par(mfrow=c(3,1))
for (i in 1:3) {
  matplot(agefine, eigen_mass[, i], type = "l", 
          xlab = "time", ylab = "",
          main = paste("Growth: Principal Component", i))
}

## Test for orthogonality
eigen_mass[,1] %*% eigen_mass[,2] # -1.63064e-16
eigen_mass[,1] %*% eigen_mass[,3] # 1.249001e-16
eigen_mass[,2] %*% eigen_mass[,3] # -8.326673e-17

## For the log growth curve
fpcaobj_logmass <- prcomp(x=t(logmass_aligned), retx = TRUE, center = TRUE) 
sdeval_logmass <- fpcaobj_logmass$sdev # the singular value of X
eigen_logmass <- fpcaobj_logmass$rotation # eigen vectors
pcscores_logmass <- fpcaobj_logmass$x # principal component scores
summary(fpcaobj_logmass)
# Cumulative Proportion  0.9348 0.99781 0.99877 (up to PC3)

par(mfrow=c(3,1))
for (i in 1:3) {
  matplot(agefine, eigen_logmass[, i], type = "l", 
          xlab = "time", ylab = "",
          main = paste("Log Growth: Principal Component", i))
}

## Test for orthogonality
eigen_logmass[,1] %*% eigen_logmass[,2] # -7.914676e-17
eigen_logmass[,1] %*% eigen_logmass[,3] # 5.20417e-17
eigen_logmass[,2] %*% eigen_logmass[,3] # -2.081668e-17
