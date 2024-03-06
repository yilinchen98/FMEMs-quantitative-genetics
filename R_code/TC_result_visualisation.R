library(fdasrvf)
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

# Rescale time interval to [0,1]
## x = (x -min(x))/(max(x) - min(x))
age_list_new <- list()
for (i in 1:N){
  age_list_new[[i]] = (age_list[[i]]-min(age_list[[i]]))/(max(age_list[[i]])-min(age_list[[i]]))
}

df$x <- unsplit(age_list_new,df$id) 

# Curve alignment Use align_fPCA from the package fdasrvf

## Use smooth.spline as basis functions and transform the fitted mass to a fine time grid
agefine = seq(0,1,length=100)

### for the growth curve
mass_hat <- matrix(0,100,N) 

for (i in 1:N){
  fit_mass <- smooth.spline(age_list_new[[i]], trait_list[[i]], cv=FALSE)
  fit_mass_hat <- predict(fit_mass,agefine)
  mass_hat[,i] <- fit_mass_hat$y
}

## Use align_fPCA to align all curves
register_growth_curve <- align_fPCA(f = mass_hat, t = agefine, showplot = FALSE)
mass_aligned <- register_growth_curve$fn
mass_warp <- register_growth_curve$gam

## Compute the mean curve
mean_mass <- apply(mass_aligned, 1, mean)

## Plot the aligned mass curve and the population mean
par(mfrow=c(1,2))
matplot(agefine,mass_warp, type = "l", col=1:N, 
        xlab = "time", ylab="warping function",
        main = "Warping Functions")
matplot(agefine,mass_aligned, type = "l", col=1:N, 
        xlab = "time", ylab = "mass",
        main="Aligned Growth Curves")

lines(agefine, mean_mass, col="red", lwd=3)
legend(x="bottomright", legend="mean", lwd=3.0, col="red")

# Run FPCA on the aligned curves
## For the growth curve with zero mean
fpcaobj_mass <- prcomp(x=t(mass_aligned), retx = TRUE, center = TRUE) 
sdeval_mass <- fpcaobj_mass$sdev # the singular value of X
eigen_mass <- fpcaobj_mass$rotation # eigen vectors
pcscores_mass <- fpcaobj_mass$x # principal component scores

# Convert discrete eigenvectors to eigenfunctions 

## create an empty list which stores eigenfunctions for 873 subjects
## evaluated at the original time points.
phi_list <- list() 

#$ Interpolate eigenvectors for all subjects
for (i in 1:N){
  phi <- convert_to_basisfunctions(t = agefine, eigenvecs = eigen_mass[,1:3],
                                   tout = age_list_new[[i]], method = "linear")
  phi_list[[i]] <- phi
}

phi <- do.call(rbind,phi_list)
colnames(phi) <- c("phi1", "phi2", "phi3")

## update the dataframe to include basis functions
df <- cbind(df,phi)

# Extract the relationship matrix A from the data

pos = df$id[FirstUniqueIdPos] # extract ids for all subjects
sire_id = df$sire[FirstUniqueIdPos] # extract ids for sire
dam_id = df$dam[FirstUniqueIdPos] # extract ids for dam

pede <- editPed(sire = sire_id, dam = dam_id, label = pos)
ped<- with(pede, pedigree(label=label, sire=sire, dam=dam))
A <- getA(ped)[163:1035,163:1035]


#Fit mixed-effect models
rrFormula <- trait ~ df$phi1 + df$phi2 + df$phi3 + (-1 + df$phi1 + df$phi2 + df$phi3 | df$id) + 
  (-1 + df$phi1 + df$phi2 + df$phi3 | df$id) # mixed-effect model formula
system.time(
  fit_RR <- fit_genetic_fmm(formula= rrFormula, data=df, A = A, phi = phi)
) # user  system  elapsed 
  # 61.79  3.45    78.64 


# Let's fit the genetic mixed-effect model with centered data

## Centered the data
mean_mass_list <- list()
for (i in 1:N){
  mean_interp <- convert_to_basisfunctions(agefine, mean_mass, age_list_new[[i]])
  mean_mass_list[[i]] <- mean_interp
}
mean_mass_curve <- unsplit(mean_mass_list,df$id)
df$centered_trait <- df$trait - mean_mass_curve

## Fit mixed-effect model
mrformula <- centered_trait ~ -1 + (-1 + df$phi1 + df$phi2 + df$phi3 | df$id) + 
  (-1 + df$phi1 + df$phi2 + df$phi3 | df$id)
system.time(
  fit_mean_rand <- fit_genetic_fmm(mrformula, df, A, phi)
) # user  system  elapsed 
  # 51.56  3.40    72.53      

#####################################################################################################

# Compare the two models
## Extract the fixed effect in fit_RR and compare with the population mean calculated from the aligned data
summary(fit_RR)
summary(fit_mean_rand)
fixed_effect_coefs <- fixef(fit_RR) # extract the fixed-effect coefs
X<- getME(fit_RR, name="X") # fixed-effect design matrix
fixed_effect <- X %*% fixed_effect_coefs
fixef_list <- split(fixed_effect, df$id)

par(mfrow=c(1,1))

plot(c(0,1),c(0,300), type="n", xlab="time", ylab="mass")
for (i in 1:873){
  lines(age_list_new[[i]], fixef_list[[i]], col=i)
}

lines(agefine, mean_mass, type="l", lwd=3.0, col="red")


mean_estimated <- matrix(0, 100, 873)
for (i in 1:873){
  mean_estimated[,i] <- approx(age_list_new[[i]], fixef_list[[i]], agefine)$y
}

plot(agefine, apply(mean_estimated, 1, mean), type = "l", xlab="time", ylab="mass")
lines(agefine, mean_mass, type="l", col="red")
legend("bottomright", legend=c("fixed-effect", "population mean"), col=c("black", "red"), lwd=1.0)

## Plot the fitted curves of RR and compare with the original data

traitRR_list <- split(fitted(fit_RR), df$id) # fitted value unsmoothed 
par(mfrow=c(1,3))
plot(c(0,1), c(-30,400), type="n", xlab="time", ylab="mass", main="Fitted Value of RR")
for (i in 1:N){
  lines(age_list_new[[i]], traitRR_list[[i]], type="l", col=i)
}

traitRR_smooth <- list() # use spline smoothing fitted data
plot(c(0,1), c(-30,400), type="n", xlab="time", ylab="mass", main="Fitted Value (RR) on a fine grid")
for (i in 1:N){
  smofitTraitRR <- smooth.spline(age_list_new[[i]], traitRR_list[[i]], cv=FALSE)
  smofitTraitRR_hat <- predict(smofitTraitRR,agefine)
  lines(smofitTraitRR_hat,col=i)
  traitRR_smooth[[i]] <- smofitTraitRR_hat$y
}

plot(c(0,1), c(-30,400), type="n", xlab="time", ylab="mass", main="Smoothed Growth Plot")
for (i in 1:N){
lines(agefine, mass_aligned[,i], type="l", col=i)
}

## Plot the fitted curves of MR and compare with the original data

traitMR_list <- split(fitted(fit_mean_rand), df$id) # fitted value unsmoothed 
par(mfrow=c(1,3))
plot(c(0,1), c(-150,150), type="n", xlab="time", ylab="mass", main="Fitted Value of MR")
for (i in 1:N){
  lines(age_list_new[[i]], traitMR_list[[i]], type="l", col=i)
}

traitMR_smooth <- list() # use spline smoothing fitted data
plot(c(0,1), c(-150,150), type="n", xlab="time", ylab="mass", main="Fitted Value (MR) on a fine grid")
for (i in 1:N){
  smofitTraitMR <- smooth.spline(age_list_new[[i]], traitMR_list[[i]], cv=FALSE)
  smofitTraitMR_hat <- predict(smofitTraitMR,agefine)
  lines(smofitTraitMR_hat,col=i)
  traitMR_smooth[[i]] <- smofitTraitMR_hat$y
}

plot(c(0,1), c(-150,150), type="n", xlab="time", ylab="mass", main="Smoothed Centered Growth Plot")
for (i in 1:N){
  lines(agefine, mass_aligned[,i]- mean_mass, type="l", col=i)
}

## Extract the genetic and environmental covariance matrix
## variance-covariance matrices from RR model

VC_RR <- as.matrix(as.data.frame(VarCorr(fit_RR))["vcov"])

CG_RR <- matrix(0, 3,3)
CG_RR[1:3, 1:3] <- VC_RR[1:3]
CG_RR[1,2] <- VC_RR[4]
CG_RR[1,3] <- VC_RR[5]
CG_RR[2,3] <- VC_RR[6]
CG_RR[2,1] <- CG_RR[1,2]
CG_RR[3,2] <- CG_RR[2,3]
CG_RR[3,1] <- CG_RR[1,3]


CE_RR <- matrix(0, 3,3)
CE_RR[1:3, 1:3] <- VC_RR[7:9]
CE_RR[1,2] <- VC_RR[10]
CE_RR[1,3] <- VC_RR[11]
CE_RR[2,3] <- VC_RR[12]
CE_RR[2,1] <- CE_RR[1,2]
CE_RR[3,2] <- CE_RR[2,3]
CE_RR[3,1] <- CE_RR[1,3]

### Convert to genetic covariance function
CG_RR_fun <- estimate_geneticVar_funcs(agefine, CG_RR, eigen_mass[,1:3])
### Plot the genetic covariance function
fig_RR <- plot_ly(x = agefine, y = agefine, z = ~CG_RR_fun, 
                  type = "surface", showscale = FALSE)

fig_RR %>% layout(title = "Genetic Variance Function")

### environmental covariance function
CE_RR_fun <- estimate_geneticVar_funcs(agefine, CE_RR, eigen_mass[,1:3])
### Plot the genetic covariance function
fig_RR1 <- plot_ly(x = agefine, y = agefine, z = ~CE_RR_fun, 
                  type = "surface", showscale = FALSE)

fig_RR1 %>% layout(title = "Encironmental Variance Function")

## genetic covariance matrices from MR model

VC_MR<- as.matrix(as.data.frame(VarCorr(fit_mean_rand))["vcov"])

CG_MR <- matrix(0, 3,3)
CG_MR[1:3, 1:3] <- VC_MR[1:3]
CG_MR[1,2] <- VC_MR[4]
CG_MR[1,3] <- VC_MR[5]
CG_MR[2,3] <- VC_MR[6]
CG_MR[2,1] <- CG_MR[1,2]
CG_MR[3,2] <- CG_MR[2,3]
CG_MR[3,1] <- CG_MR[1,3]


CE_MR <- matrix(0, 3,3)
CE_MR[1:3, 1:3] <- VC_MR[7:9]
CE_MR[1,2] <- VC_MR[10]
CE_MR[1,3] <- VC_MR[11]
CE_MR[2,3] <- VC_MR[12]
CE_MR[2,1] <- CE_MR[1,2]
CE_MR[3,2] <- CE_MR[2,3]
CE_MR[3,1] <- CE_MR[1,3]

CG_MR_fun <- estimate_geneticVar_funcs(agefine, CG_MR, eigen_mass[,1:3])
fig_MR <- plot_ly(x = agefine, y = agefine, z = ~CG_MR_fun, type = "surface", showscale = FALSE)

fig_MR %>% layout(title = "Genetic Variance Function")

CE_MR_fun <- estimate_geneticVar_funcs(agefine, CE_MR, eigen_mass[,1:3])
### Plot the genetic covariance function
fig_MR1 <- plot_ly(x = agefine, y = agefine, z = ~CE_MR_fun, 
                   type = "surface", showscale = FALSE)

fig_MR1 %>% layout(title = "Encironmental Variance Function")

