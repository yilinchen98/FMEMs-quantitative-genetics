library(fdasrvf)
library(pedigreemm)
library(lme4)
library(Matrix)

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
# Curve alignment
## Use align_fPCA from the package fdasrvf
## Use smooth.spline as basis functions
## Transform the fitted mass to a fine time grid
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

## Use align_fPCA to align all curves
register_growth_curve <- align_fPCA(f = mass_hat, t = agefine, showplot = FALSE)
mass_aligned <- register_growth_curve$fn
mass_warp <- register_growth_curve$gam

## Repeat the process for log(mass)
register_logrowth_curve <- align_fPCA(f = logmass_hat, t = agefine, showplot = FALSE)
logmass_aligned <- register_logrowth_curve$fn
logmass_warp <- register_logrowth_curve$gam

# Run FPCA on the aligned curves

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

# Fit FMMs
## Extract the relationship matrix A from the data

pos = df$id[FirstUniqueIdPos] # extract ids for all subjects
sire_id = df$sire[FirstUniqueIdPos] # extract ids for sire
dam_id = df$dam[FirstUniqueIdPos] # extract ids for dam

pede <- editPed(sire = sire_id, dam = dam_id, label = pos)
ped<- with(pede, pedigree(label=label, sire=sire, dam=dam))
A <- getA(ped)[163:1035,163:1035]

## Reparameterise the design matrix
Lt <- chol(A) # upper triangular matrix of Cholesky decomposition
p = 3 # dim of genetic random effect vector
I <- as(diag(p), "dgCMatrix")

M <- kronecker(Lt, I) # use in Z* = ZM

## Define the LME formula

### define the basis functions (eigenfunctions)
phi_list <- list() 
# create an empty list which stores eigenfunctions for 873 subjects
# evaluated at the original time points.

for (i in 1:N){
  phi <- convert_to_basisfunctions(t = agefine, eigenvecs = eigen_mass[,1:p],
                                       tout = age_list_new[[i]])
  phi_list[[i]] <- phi
}

phi <- do.call(rbind,phi_list)
colnames(phi) <- c("phi1", "phi2", "phi3")

### update the dataframe to include basis functions
df <- cbind(df,phi)

### now only focus on the environmental random effect 
Zformula <- trait ~ -1 + (-1 + df$phi1 + df$phi2 + df$phi3 | df$id)
ZparsedForm <- lFormula(formula = Zformula, data = df)
ZE <- t(ZparsedForm$reTrms$Zt) # transpose of the environmental random effect matrix

### compute the genetic random effect
Z_prep <- ZE # duplicate the environmental random effect matrix
ZG <- Z_prep %*% M

### combine two random effect matrices
Z <- cbind(ZG, ZE)

######################################################################
#Test on a smaller set, just containing 3 subjects
test_df <- df[1:25,]

### Calculate the relation matrix A
test_pede <- editPed(sire = as.character(c(1,1,1)), dam = as.character(c(101,101,101)), 
                     label = as.character(c(10001,10002,10003)))
test_ped <- with(test_pede, pedigree(label=label,sire = sire, dam=dam))
test_A <- as.matrix(getA(test_ped)[3:5,3:5])
test_M <- as.matrix(kronecker(chol(test_A),I))
all(test_M == as.matrix(M)[1:9,1:9]) # True

### Fit only random effect
test_formula <- trait ~ -1 + (-1 + test_df$phi1 + test_df$phi2 + test_df$phi3 | test_df$id)
test_parsedFormula <- lFormula(formula = test_formula, data = test_df)
testZEt <- test_parsedFormula$reTrms$Zt 
testZEt_standard <- t(as.matrix(testZEt)) # environmental random effect matrix

testZG_prep <- t(testZEt) # duplicate ZE
testZG <- testZG_prep %*% test_M # update the genetic random effect design matrix
testZG_standard <- as.matrix(testZG_prep)
testZ <- cbind(testZG_standard,testZEt_standard) # genetic random effect matrix
