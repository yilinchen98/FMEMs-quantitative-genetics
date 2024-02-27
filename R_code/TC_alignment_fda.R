library(fda)
library(Matrix)

# Import data
setwd("D:/KCL_2023-2027_PhD/Year 1/Genetics/FMEMs-quantitative-genetics/R_code")
TRFUN25PUP4 = read.delim("TRFUN25PUP4.DAT",header = FALSE)
names(TRFUN25PUP4)<-c("id","sire","dam","trait","x")
df <- data.frame(TRFUN25PUP4) # n = 6860 observations

# Plot raw data
indices <- df$id
trait <- df$trait
time <- df$x
FirstUniqueIdPos <- which(duplicated(indices) == FALSE)
N = length(FirstUniqueIdPos) # N = 873 subjects
age_list <- split(time,indices)
trait_list <- split(trait,indices)

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

# Data Smoothing
## Rescale time interval to [0,1]
## x = (x -min(x))/(max(x) - min(x))
age_list_new <- list()
for (i in 1:N){
  age_list_new[[i]] = (age_list[[i]]-min(age_list[[i]]))/(max(age_list[[i]])-min(age_list[[i]]))
}
## check the first subject
print(age_list_new[[1]])

basis <- create.bspline.basis(c(0,1), 10, 4)  # create cubic bspline basis
## Choose smoothing parameter lambda
loglambda <- seq(-6, 0, by=0.25)
n <- length(loglambda)
gcv <- matrix(0,n,5)
for (j in 1:5){
  for (i in 1:n ){
    lambda <- 10^loglambda[i]
    fdPar <- fdPar(basis, 2, lambda)
    gcvi <- smooth.basis(age_list_new[[j]],trait_list[[j]],fdPar)$gcv
    gcv[i,j] <- sum(gcvi)
  }
}
par(mfrow = c(1,1))
matplot(loglambda,gcv,type="l", col = 1:5,xlab="log(lambda)")

## smoothing by roughness penalty
fdPar <- fdPar(basis, 2, 1e-4) # with penalty level = 10^-4
mass_coefs <- matrix(0,10,N) 
logmass_coefs <- matrix(0,10,N) 
for (i in 1:N){
  mass_rp_i <- smooth.basis(age_list_new[[i]], trait_list[[i]],fdPar)$fd
  mass_coefs[,i] <- mass_rp_i$coefs
  logmass_rp_i <- smooth.basis(age_list_new[[i]], log10(trait_list[[i]]),fdPar)$fd
  logmass_coefs[,i] <- logmass_rp_i$coefs
}

## calculating the mean
mass.fd <- fd(mass_coefs, basis)
logmass.fd <- fd(logmass_coefs,basis)
mean_mass <- mean.fd(mass.fd)
mean_logmass <- mean.fd(logmass.fd)

## Visualise the smoothed curves
par(mfrow=c(2,1))
plot.fd(mass.fd, ylab="mass", main = "Smoothed Mass Plot")
lines(mean_mass, type = "l", lwd = 2.0, col="red")
legend(x="topleft", legend="Mean", col="red", lwd = 2)
plot.fd(logmass.fd, ylab = "log(mass)", main="Smoothed Log Mass Plot")
lines(mean_logmass, type = "l", lwd = 2.0, col = "red")
legend(x="bottomright", legend="Mean", col="red",lwd = 2)

# Curve registration using register.fd
## all curves are registered to the mean
wbasisCR = create.bspline.basis(c(0,1), 10, 4)
Wfd0CR = fd(matrix(0,10,N),wbasisCR)
WfdParCR = fdPar(Wfd0CR, 2, 10^-0.5)

regListMass = register.fd(mean_mass,mass.fd, WfdParCR) # register mass curve
massfdCR = regListMass$regfd # aligned mass curves
mass_warpfdCR = regListMass$warpfd # warping functions
massWfdCR = regListMass$Wfd

regListLogmass = register.fd(mean_logmass,logmass.fd, WfdParCR) # register log mass curve
logmassfdCR = regListLogmass$regfd # aligned log mass curves
logmass_warpfdCR = regListLogmass$warpfd # warping function
logmassWfdCR = regListLogmass$Wfd

par(mfrow=c(1,2))
plot.fd(mass_warpfdCR, main="warping functions")
plot.fd(massfdCR, main="aligned mass curves")

par(mfrow=c(1,2))
plot.fd(logmass_warpfdCR, main="warping functions")
plot.fd(logmassfdCR, main="aligned log mass curves")
