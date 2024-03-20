library(Matrix)
library(fdasrvf)
library(lme4)
library(pedigreemm)
library(npreg)
library(ggplot2)
library(plotly)

# Import data
setwd("D:/KCL_2023-2027_PhD/Year 1/Genetics/FMEMs-quantitative-genetics/R_code")
TRFUN25PUP4 = read.delim("TRFUN25PUP4.DAT",header = FALSE)
names(TRFUN25PUP4)<-c("id","sire","dam","trait","x")
df <- data.frame(TRFUN25PUP4)

FirstUniqueIdPos <- which(duplicated(df$id) == FALSE)
N = length(FirstUniqueIdPos) # N = 873 subjects
n = length(df$id) # n = 6860 observations
age_list <- split(df$x,df$id)
trait_list <- split(df$trait,df$id)

## Rescale time interval to [0,1]
## x = (x -min(x))/(max(x) - min(x))
age_list_new <- list()
for (i in 1:N){
  age_list_new[[i]] = (age_list[[i]]-min(age_list[[i]]))/(max(age_list[[i]])-min(age_list[[i]]))
}

df$x_rescaled <- unsplit(age_list_new,df$id)

# Plot the raw data
par(mfrow=c(1,2))
plot(c(0,25), c(0,400), xlab="Time", ylab="Mass", main="Growth Plot", type="n")
for (i in 1:N){
  lines(age_list[[i]], trait_list[[i]],type="l", col=i)
}

plot(c(0,25), c(0,6), type = 'n',xlab = 'Time', 
     ylab = 'Mass', main = 'Log Growth Plot ')
for (i in 1:N) {
  lines(age_list[[i]], log(trait_list[[i]]), col = i)
}

# Data smoothing
## We get negative values when smoothing on original scale.
## Smooth growth curves on log scale and then take exp to recover body mass.
## This imposes positive smoothing (same methodology as smooth.pos() in fda package).

agefine <- seq(0,1,length=100)
mass_hat <- matrix(0,100,N) # store smooth growth curves on original scale
logmass_hat <- matrix(0,100,N) # store smooth log growth curves
pred_mass <- matrix(0,100,N) # store the smoothed mass recovered by taking exp
lambda_mass <- rep(0,N) # smoothing parameter used for each growth curve
lambda_logmass <- rep(0,N) # smoothing parameter used for each log growth curve
for (i in 1:N){
  ss_mass <- smooth.spline(age_list_new[[i]], trait_list[[i]], cv=FALSE, 
                           all.knots=TRUE)
  ss_logmass <- smooth.spline(age_list_new[[i]], log(trait_list[[i]]), cv=FALSE,
                              all.knots=TRUE)
  mass_hat[,i] <- predict(ss_mass,agefine)$y
  logmass_hat[,i] <- predict(ss_logmass, agefine)$y
  pred_mass[,i] <- exp(predict(ss_logmass,agefine)$y)
  lambda_mass[i] <- ss_mass$lambda
  lambda_logmass[i] <- ss_logmass$lambda
}

## Some smoothed curves have negative values near t=0, let us select 
## the corresponding curves and plot the original measurements.
neg_cols <- unique(which(mass_hat < 0, arr.ind = TRUE)[,2])
length(neg_cols) 
neg_mass<- mass_hat[,neg_cols] # matrix each column represents wrongly smoothed growth curves
colnames(neg_mass) <- as.character(neg_cols)

## There are 77 questionable smoothed growth curves on the original scale.
## On each plot, black points represent the data measurements;
## blue curve represents smoothing on the original scale;
## red curve represents smoothing by imposing positive constraint. 
par(mfrow=c(3,3))
for (i in 1:length(neg_cols)){
  plot(age_list_new[[neg_cols[i]]], trait_list[[neg_cols[i]]], type="p", xlab="Time", ylab="Mass", main=paste("Subject",neg_cols[i]))
  lines(agefine,mass_hat[,neg_cols[i]],type="l", col="blue")
  lines(agefine, pred_mass[,neg_cols[i]],type="l", col="red")
}

## Check the largest and smallest number of measurements per subject
##Initialize variables to hold the maximum and minimum lengths
max_length <- -Inf
min_length <- Inf

## Iterate through each list in trait_list
for (i in 1:length(trait_list)) {
  current_length <- length(trait_list[[i]])
  if (current_length > max_length) {
    max_length <- current_length
  }
  if (current_length < min_length) {
    min_length <- current_length
  }
}

print(paste("Largest number of elements in a list:", max_length))
print(paste("Smallest number of elements in a list:", min_length))

# Plot the smoothed growth curve with positive constraint
par(mfrow=c(3,3))
for (i in 1:N){
  plot(age_list_new[[i]], trait_list[[i]], type = "p", xlab="Time", ylab="Mass", main=paste("Subject", i))
  lines(agefine, pred_mass[,i], type="l", col="red")
  legend("bottomright", legend = paste("lambda=", lambda_logmass[i]))
}

## Plot graph with large lambda
par(mfrow=c(3,3))
large_lam <- which(lambda_logmass > 1e-4)
for (i in 1:length(large_lam)){
  plot(age_list_new[[large_lam[i]]], trait_list[[large_lam[i]]], type="p", xlab="Time", ylab="Mass", main=paste("Subject",large_lam[i]))
  lines(agefine, pred_mass[,large_lam[i]],type="l", col="red")
}

## restrict lambda <=1e-4
mass_adjusted <- matrix(0,100, length(large_lam))
for ( i in 1: length(large_lam)){
  ss_new <- smooth.spline(age_list_new[[large_lam[i]]], log(trait_list[[large_lam[i]]]), all.knots = TRUE, lambda = 1e-4)
  mass_adjusted[,i] <- exp(predict(ss_new, agefine)$y)
}
colnames(mass_adjusted) <- as.character(large_lam)

par(mfrow=c(3,3))
for (i in 1:length(large_lam)){
  plot(age_list_new[[large_lam[i]]], trait_list[[large_lam[i]]], type="p", xlab="Time", ylab="Mass", main=paste("Subject",large_lam[i]))
  lines(agefine, mass_adjusted[,i],type="l", col="red")
}

## Let reform the smoothed growth curves into a matrix
trait_hat <- pred_mass 
for (i in 1: length(large_lam)){
  trait_hat[,large_lam[i]] <- mass_adjusted[,i]
}

## Plot the smoothed growth curves
par(mfrow=c(1,2))

plot(c(0,1), c(0,400), xlab="Time", ylab="Mass", main="Growth Plot (rescaled time)", type="n")
for (i in 1:N){
  lines(age_list_new[[i]], trait_list[[i]],type="l", col=i)
}

matplot(agefine, trait_hat, col=1:N, type = "l", xlab="Time", ylab="Mass", main="Smoothed Growth Plot")


# Curve alignment
aligned_mass_process <- time_warping(f=trait_hat, time=agefine)
aligned_mass_curve <- aligned_mass_process$fn
aligned_mean <- aligned_mass_process$fmean
warping_funs <- aligned_mass_process$warping_functions


par(mfrow=c(1,3))
plot(c(0,1), c(0,1), type = 'n',xlab = 'Time', 
     ylab = 'warping functions')
for (i in 1:N){
  lines(agefine, warping_funs[,i], type="l" , col=i)
}


plot(c(0,1), c(0,400), type = 'n',xlab = 'Time', 
     ylab = 'Mass', main = 'Aligned Growth Plot')
for (i in 1:N){
  lines(agefine, aligned_mass_curve[,i], type="l",col=i)
}
lines(agefine, aligned_mean, lwd = 3.0, col="red")
legend("bottomright", legend="mean", lwd=3.0, col="red")

matplot(agefine, trait_hat, col=1:N, type = "l", xlab="Time", ylab="Mass", main="Smoothed Growth Plot")

# FPCA
fpcaobj_mass <- prcomp(x=t(aligned_mass_curve), retx = TRUE, center = TRUE, rank. = 3)
eigen_mass <- fpcaobj_mass$rotation # eigen vectors

par(mfrow=c(3,1))
for (i in 1:3) {
  plot(agefine, eigen_mass[, i], type = "l", 
       xlab = "time", ylab = "",
       main = paste("Growth Curve Principal Component", i))
}

## Test for orthogonality
eigen_mass[,1] %*% eigen_mass[,2] 
eigen_mass[,1] %*% eigen_mass[,3] 
eigen_mass[,2] %*% eigen_mass[,3] 

## Run FPCA on smoothed but not aligned curves
fpcaobj_smass <- prcomp(x=t(trait_hat), retx = TRUE, center = TRUE, rank. = 3)
pcs <- fpcaobj_smass$rotation

par(mfrow=c(3,1))
for (i in 1:3) {
  plot(agefine, pcs[, i], type = "l", 
       xlab = "time", ylab = "",
       main = paste("Growth Curve (not aligned) Principal Component", i))
}

# Reform the aligned data to a dataframe for model fitting
subjectID <- rep(unique(df$id), each=100)
trait_pred <- c(aligned_mass_curve)
basis1 <- rep(eigen_mass[,1], times = 873)
basis2 <- rep(eigen_mass[,2], times = 873)
basis3 <- rep(eigen_mass[,3], times = 873)
new_df <- data.frame(subjectID, trait_pred, basis1, basis2, basis3)
names(new_df) <- c("subjectID", "trait_hat", "basis1", "basis2", "basis3")

# Pedigree and genetic relationship matrix A
pos = df$id[FirstUniqueIdPos] # extract ids for all subjects
sire_id = df$sire[FirstUniqueIdPos] # extract ids for sire
dam_id = df$dam[FirstUniqueIdPos] # extract ids for dam

pede <- editPed(sire = sire_id, dam = dam_id, label = pos)
ped<- with(pede, pedigree(label=label, sire=sire, dam=dam))
A <- getA(ped)[163:1035,163:1035]

# Fit FMEMs
## Fit fixed and random effects with the same basis
fmeFormula <- trait_hat ~ new_df$basis1 + new_df$basis2 + new_df$basis3 + 
  (-1 + new_df$basis1 + new_df$basis2 + new_df$basis3 | new_df$subjectID) + 
  (-1 + new_df$basis1 + new_df$basis2 + new_df$basis3 | new_df$subjectID)
system.time(
  fit_sameBasis <- fit_genetic_fmm(formula= fmeFormula, data=new_df, A = A, phi = eigen_mass)
) # user  system  elapsed 
# 1094.82 92.27   1361.03 

summary(fit_sameBasis)

# Visualise result
## Extract the fixed effect in fit_sameBasis and compare with the population mean calculated from the aligned data
fe_coefs <- fixef(fit_sameBasis) # extract the fixed-effect coefs
fixef <- eigen_mass%*%fe_coefs[2:4] + fe_coefs[1]

par(mfrow=c(1,1))
plot(agefine, fixef, type = "l", xlab="time", ylab="mass")
lines(agefine, aligned_mean, type="l", col="red")
legend("bottomright", legend=c("fixed-effect", "population mean"), col=c("black", "red"), lwd=1.0)

## Plot the fitted curves of RR and compare with the original data
fitted_list <- split(fitted(fit_sameBasis), new_df$subjectID) # fitted value 
par(mfrow=c(1,2))
plot(c(0,1), c(0,400), type="n", xlab="time", ylab="mass", main="Fitted Value of RR")
for (i in 1:N){
  lines(agefine, fitted_list[[i]], type="l", col=i)
}
plot(c(0,1), c(0,400), type="n", xlab="time", ylab="mass", main="Smoothed Growth Plot")
for (i in 1:N){
  lines(agefine, aligned_mass_curve[,i], type="l", col=i)
}

## Extract the genetic and environment covariance matrices
VC <- as.matrix(as.data.frame(VarCorr(fit_sameBasis))["vcov"])

CG <- matrix(0, 3,3) # genetic covariance matrix
CG[1:3, 1:3] <- VC[1:3]
CG[1,2] <- VC[4]
CG[1,3] <- VC[5]
CG[2,3] <- VC[6]
CG[2,1] <- CG[1,2]
CG[3,2] <- CG[2,3]
CG[3,1] <- CG[1,3]

CE <- matrix(0, 3,3) # environment covariance matrix
CE[1:3, 1:3] <- VC[7:9]
CE[1,2] <- VC[10]
CE[1,3] <- VC[11]
CE[2,3] <- VC[12]
CE[2,1] <- CE[1,2]
CE[3,2] <- CE[2,3]
CE[3,1] <- CE[1,3]
 
### Convert to genetic covariance function
CG_fun <- eigen_mass %*% CG %*% t(eigen_mass)
### environmental covariance function
CE_fun <- eigen_mass %*% CE %*% t(eigen_mass)
### Phenotypic covariance function
P_fun <- CG_fun + CE_fun

# Plot the genetic covariance function
fig_RR1 <- plot_ly(x = agefine, y = agefine, z = ~CG_fun, 
                   type = "surface", showscale = FALSE) %>% 
  layout(title = "Genetic Variance Function", 
         scene = list(zaxis = list(title=""),
                      xaxis = list(title = "time"),
                      yaxis = list(title = "time")))

# Plot the environmental covariance function
fig_RR2 <- plot_ly(x = agefine, y = agefine, z = ~CE_fun, 
                   type = "surface", showscale = FALSE) %>% 
  layout(title = "Environment Variance Function", 
         scene = list(zaxis = list(title=""),
                      xaxis = list(title = "time"),
                      yaxis = list(title = "time")))

# Plot the phenotypic covariance function
fig_RR3 <- plot_ly(x = agefine, y = agefine, z = ~P_fun, 
                   type = "surface", showscale = FALSE) %>% 
  layout(title = "Phenotypic Variance Function", 
         scene = list(zaxis = list(title=""),
                      xaxis = list(title = "time"),
                      yaxis = list(title = "time")))

# Display each plot separately

fig_RR1
fig_RR2
fig_RR3
