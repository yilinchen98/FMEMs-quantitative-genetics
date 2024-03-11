library(fdasrvf)
library(lme4)
library(Matrix)
library(pedigreemm)
library(plotly)
library(ggplot2)

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

df$x <- unsplit(age_list_new,df$id)

## calculate the log mass
df$logtrait <- log10(df$trait)

# Data visualisation
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

#Data smoothing
## Use smooth.spline as basis functions and transform the fitted mass to a fine time grid
agefine = seq(0,1,length=100)

## for the growth curve
mass_hat <- matrix(0,100,N) 
lambda_massi <- rep(0,N) 
plot(c(0,1),c(0,400), type="n", xlab="Days",ylab="Mass",
     main = "Growth plot on a fine grid")
for (i in 1:N){
  fit_mass <- smooth.spline(age_list_new[[i]], trait_list[[i]], cv=FALSE)
  fit_mass_hat <- predict(fit_mass,agefine)
  lines( fit_mass_hat,col=i)
  mass_hat[,i] <- fit_mass_hat$y
  lambda_massi[i] <- fit_mass$lambda
} ## The aligned growth plot have negative values near t=0. 
## Use the log(mass) for further analysis

## for log growth curve
logmass_hat <- matrix(0,100,N)
lambda_logmassi <- rep(0,N)
plot(c(0,1),c(0,3), type="n", xlab="Days",ylab="log(Mass)",
     main = "Log growth plot on a fine grid")
for (i in 1:N){
  fit_logmass <- smooth.spline(age_list_new[[i]], log10(trait_list[[i]]), cv=FALSE)
  fit_logmass_hat <- predict(fit_logmass,agefine)
  lines( fit_logmass_hat,col=i)
  logmass_hat[,i] <- fit_logmass_hat$y
  lambda_logmassi[i] <- fit_logmass$lambda
}

# Curve registration
register_growth_curve <- align_fPCA(f = mass_hat, t = agefine, showplot = FALSE)
mass_aligned <- register_growth_curve$fn
mass_warp <- register_growth_curve$gam

register_logrowth_curve <- align_fPCA(f = logmass_hat, t = agefine, showplot = FALSE)
logmass_aligned <- register_logrowth_curve$fn
logmass_warp <- register_logrowth_curve$gam

## Compute the mean curve
mean_mass <- apply(mass_aligned, 1, mean)
mean_logmass <- apply(logmass_aligned, 1, mean)

## Plot the aligned (log)mass curve and the population mean
par(mfrow=c(1,2))
matplot(agefine,mass_warp, type = "l", col=1:N, 
        xlab = "time", ylab="warping function",
        main = "Warping Functions")
matplot(agefine,mass_aligned, type = "l", col=1:N, 
        xlab = "time", ylab = "mass",
        main="Aligned Growth Curves")
lines(agefine, mean_mass, col="red", lwd=3)
legend(x="bottomright", legend="mean", lwd=3.0, col="red")

par(mfrow=c(1,2))
matplot(agefine,logmass_warp, type = "l", col=1:N, xlim = c(0,1),
        xlab = "time", ylab="warping function",
        main = "Warping Functions")
matplot(agefine,logmass_aligned, type = "l", col=1:N, xlim = c(0,1),
        xlab = "time", ylab = "mass",
        main="Aligned Log Growth Curves")
lines(agefine, mean_logmass, col="red", lwd=3)
legend(x="bottomright", legend="mean", lwd=3.0, col="red")

# Run FPCA on the aligned curve
fpcaobj_mass <- prcomp(x=t(mass_aligned), retx = TRUE, center = TRUE, rank. = 3)
eigen_mass <- fpcaobj_mass$rotation # eigen vectors

fpcaobj_logmass <- prcomp(x=t(logmass_aligned), retx = TRUE, center = TRUE, rank. = 3) 
eigen_logmass <- fpcaobj_logmass$rotation # eigen vectors

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

# Interpolate eigenvectors into eigenfunctions
## Here use data on log scale
basis_list <- list()
for (i in 1:N){
  basis_list[[i]] <- convert_to_basisfunctions(t=agefine, eigenvecs=eigen_logmass, tout=age_list_new[[i]])
}
basis_funcs <- do.call(rbind, basis_list)
colnames(basis_funcs) <- c("basis1", "basis2", "basis3")

# Centered the data
sys_logmean_list <- list()
for (i in 1:N){
  sys_logmean_list[[i]] <- convert_to_basisfunctions(t=agefine, eigenvecs = mean_logmass, tout=age_list_new[[i]])
}
sys_logmean_curve <- unsplit(sys_logmean_list,df$id)
df$centered_logtrait <- df$logtrait - sys_logmean_curve

sys_mean_list <- list()
for (i in 1:N){
  sys_mean_list[[i]] <- convert_to_basisfunctions(t=agefine, eigenvecs = mean_mass, tout=age_list_new[[i]])
}
sys_mean_curve <- unsplit(sys_mean_list,df$id)
df$centered_trait <- df$trait - sys_mean_curve


# update the dataframe to include basis functions
df <- cbind(df, basis_funcs)

# Extract the relationship matrix A from the data

pos = df$id[FirstUniqueIdPos] # extract ids for all subjects
sire_id = df$sire[FirstUniqueIdPos] # extract ids for sire
dam_id = df$dam[FirstUniqueIdPos] # extract ids for dam

pede <- editPed(sire = sire_id, dam = dam_id, label = pos)
ped<- with(pede, pedigree(label=label, sire=sire, dam=dam))
A <- getA(ped)[163:1035,163:1035]

# Fit mixed-effect model
## First fit both fixed effect and random effects using the same basis
meFormula <- logtrait ~ df$basis1 + df$basis2 + df$basis3 + (-1 + df$basis1 + df$basis2 + df$basis3 | df$id) + 
  (-1 + df$basis1 + df$basis2 + df$basis3 | df$id) # mixed-effect model formula
system.time(
  fit_sameBasis <- fit_genetic_fmm(formula= meFormula, data=df, A = A, phi = basis_funcs)
) # user  system  elapsed 
# 108.44   7.45    129.18

summary(fit_sameBasis)

## fit the centered data
meFormula1 <- centered_logtrait ~ -1 + (-1 + df$basis1 + df$basis2 + df$basis3 | df$id) +
  (-1 + df$basis1 + df$basis2 + df$basis3 | df$id) 
system.time(
  fit_centered <- fit_genetic_fmm(formula= meFormula1, data=df, A = A, phi = basis_funcs)
)# user  system  elapsed 
# 97.86   6.58    133.02 

summary(fit_centered)

# Visualise result
## Extract the fixed effect in fit_sameBasis and compare with the population mean calculated from the aligned data
fe_coefs <- fixef(fit_sameBasis) # extract the fixed-effect coefs
fixef <- eigen_logmass%*%fe_coefs[2:4] + fe_coefs[1]

par(mfrow=c(1,1))
plot(agefine, fixef, type = "l", xlab="time", ylab="mass")
lines(agefine, mean_logmass, type="l", col="red")
legend("bottomright", legend=c("fixed-effect", "population mean"), col=c("black", "red"), lwd=1.0)

## Plot the fitted curves of RR and compare with the original data

fitted_list <- split(fitted(fit_sameBasis), df$id) # fitted value 
par(mfrow=c(1,2))
plot(c(0,1), c(0,3), type="n", xlab="time", ylab="mass", main="Fitted Value of RR")
for (i in 1:N){
  lines(age_list_new[[i]], fitted_list[[i]], type="l", col=i)
}
plot(c(0,1), c(0,3), type="n", xlab="time", ylab="mass", main="Smoothed Log Growth Plot")
for (i in 1:N){
  lines(agefine, logmass_aligned[,i], type="l", col=i)
}

fitted_list1 <- split(fitted(fit_centered), df$id) # fitted value 
par(mfrow=c(1,2))
plot(c(0,1), c(-1,1), type="n", xlab="time", ylab="mass", main="Fitted Value of RR")
for (i in 1:N){
  lines(age_list_new[[i]], fitted_list1[[i]], type="l", col=i)
}

plot(c(0,1), c(-1,1), type="n", xlab="time", ylab="mass", main="Smoothed Log Growth Plot")
for (i in 1:N){
  lines(agefine, logmass_aligned[,i]- mean_logmass, type="l", col=i)
}

## Convert covariance matrix to covariance matrix

VC1 <- as.matrix(as.data.frame(VarCorr(fit_sameBasis))["vcov"])

CG1 <- matrix(0, 3,3)
CG1[1:3, 1:3] <- VC1[1:3]
CG1[1,2] <- VC1[4]
CG1[1,3] <- VC1[5]
CG1[2,3] <- VC1[6]
CG1[2,1] <- CG1[1,2]
CG1[3,2] <- CG1[2,3]
CG1[3,1] <- CG1[1,3]


CE1 <- matrix(0, 3,3)
CE1[1:3, 1:3] <- VC1[7:9]
CE1[1,2] <- VC1[10]
CE1[1,3] <- VC1[11]
CE1[2,3] <- VC1[12]
CE1[2,1] <- CG1[1,2]
CE1[3,2] <- CG1[2,3]
CE1[3,1] <- CG1[1,3]

### Convert to genetic covariance function
CG_fun1 <- eigen_logmass %*% CG1 %*% t(eigen_logmass)
### Plot the genetic covariance function
fig_RR <- plot_ly(x = agefine, y = agefine, z = ~CG_fun1, 
                  type = "surface", showscale = FALSE)

fig_RR %>% layout(title = "Genetic Variance Function")

### environmental covariance function
CE_fun1 <- eigen_logmass %*% CE1 %*% t(eigen_logmass)
### Plot the genetic covariance function
fig_RR1 <- plot_ly(x = agefine, y = agefine, z = ~CE_fun1, 
                   type = "surface", showscale = FALSE)

fig_RR1 %>% layout(title = "Environmental Variance Function")

### Phenotypic covariance function
P_fun1 <- CG_fun1 + CE_fun1
### Plot the phenotypic covariance function
fig_RR2 <- plot_ly(x = agefine, y = agefine, z = ~P_fun1, 
                   type = "surface", showscale = FALSE)
fig_RR2 %>% layout(title = "Phenotypic Variance Function")

##################################################################

VC2 <- as.matrix(as.data.frame(VarCorr(fit_centered))["vcov"])

CG2 <- matrix(0, 3,3)
CG2[1:3, 1:3] <- VC2[1:3]
CG2[1,2] <- VC2[4]
CG2[1,3] <- VC2[5]
CG2[2,3] <- VC2[6]
CG2[2,1] <- CG2[1,2]
CG2[3,2] <- CG2[2,3]
CG2[3,1] <- CG2[1,3]


CE2 <- matrix(0, 3,3)
CE2[1:3, 1:3] <- VC2[7:9]
CE2[1,2] <- VC2[10]
CE2[1,3] <- VC2[11]
CE2[2,3] <- VC2[12]
CE2[2,1] <- CG2[1,2]
CE2[3,2] <- CG2[2,3]
CE2[3,1] <- CG2[1,3]

### Convert to genetic covariance function
CG_fun2 <- eigen_logmass %*% CG2 %*% t(eigen_logmass)
### Plot the genetic covariance function
fig <- plot_ly(x = agefine, y = agefine, z = ~CG_fun2, 
                  type = "surface", showscale = FALSE)

fig %>% layout(title = "Genetic Variance Function")

### environmental covariance function
CE_fun2 <- eigen_logmass %*% CE2 %*% t(eigen_logmass)
### Plot the genetic covariance function
fig1 <- plot_ly(x = agefine, y = agefine, z = ~CE_fun2, 
                   type = "surface", showscale = FALSE)

fig1 %>% layout(title = "Environmental Variance Function")
 
### Phenotypic covariance function
P_fun2 <- CG_fun2 + CE_fun2
### Plot the phenotypic covariance function
fig3 <- plot_ly(x = agefine, y = agefine, z = ~P_fun2, 
                   type = "surface", showscale = FALSE)
fig3 %>% layout(title = "Phenotypic Variance Function")


