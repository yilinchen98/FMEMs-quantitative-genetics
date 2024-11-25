# Here we modify the algorithm for genetic FMEM which allows different numbers 
# of basis functions for different random effects

fit_genetic_fmm2 <- function(formula, data, A, nbasis, control = lmerControl()) {
  #' This function uses the lme4 package to fit a linear mixed-effect model to genetic data, 
  #' with a specified additive genetic relationship matrix A.
  #' In this particular format, we fit both fixed effects and random effects using the same
  #' basis functions (principal components obtained from running FPCA).
  #'
  #' @param formula a two-sided linear formula object describing both the fixed-effects
  #' and random-effects of the model (as the same form used in lmer).
  #' @param data a data frame containing the variables named in formula.
  #' @param A a sparse matrix: an additive genetic relationship matrix 
  #' which models the genetic relationship in the dataset.
  #' @param nbasis numeber/vector: number of basis used to fit the random effects 
  #'        [basis for genetic, no. of basis for environment]
  #' @param control lmerControl object: control parameters for the optimizer.
  #' @return returns a fitted mixed-effect model
  
  # Load required packages
  library(lme4)
  library(Matrix)
  
  # Cholesky decomposition of A
  LA <- as(t(chol(A)), "sparseMatrix")
  if (length(nbasis) == 1){
    I_p <- as(diag(nbasis), "sparseMatrix")
  }
  else{I_p <- as(diag(nbasis[1]), "sparseMatrix")}
  MA <- kronecker(LA, I_p) # used to update the genetic design matrix Z_G = ZM
  
  # Define the mixed-model formula
  fmmParsedForm <- lFormula(formula=formula, data = data, control = control)
  
  # Compute the random-effect matrix
  Z <- t(fmmParsedForm$reTrms$Zt)
  Z[,1:dim(MA)[1]] <- Z[, 1:dim(MA)[1]] %*% MA # update the genetic-random effect matrix
  
  # Modularisation
  fmmParsedForm$reTrms$Zt <- t(Z) # Update Z in the reTrms term
  fmmDevFun <- do.call(mkLmerDevfun, fmmParsedForm) # update the objective function
  
  fmmOptimize <- optimizeLmer(devfun = fmmDevFun, control = control) # update the optimisation module
  
  # Return the mixed-effect model
  fmm <- mkMerMod(rho = environment(fmmDevFun), opt = fmmOptimize, reTrms = fmmParsedForm$reTrms, fr = fmmParsedForm$fr)
  
  return(fmm)
}

#########################################################################################################

## Test on the change of Z matrix

setwd("D:/KCL_2023-2027_PhD/FMEM_QuantitativeGenetics_Project/PhD_Project_Contents/R_code/TC_data")
TRFUN25PUP4 = read.delim("TRFUN25PUP4.DAT",header = FALSE)
names(TRFUN25PUP4)<-c("id","sire","dam","trait","x")
df <- data.frame(TRFUN25PUP4)

## calculate genetic relationship matrix
id <- df$id
FirstUniqueIdPos <- which(duplicated(id) == FALSE)
pos = id[FirstUniqueIdPos] # extract ids for all subjects
sire_id = df$sire[FirstUniqueIdPos] # extract ids for sire
dam_id = df$dam[FirstUniqueIdPos] # extract ids for dam

pede <- editPed(sire = sire_id, dam = dam_id, label = pos)
ped<- with(pede, pedigree(label=label, sire=sire, dam=dam))
A <- getA(ped)[163:1035,163:1035]

N = length(FirstUniqueIdPos) # N = 873 subjects
n = length(id) # n = 6860 observations
age_list <- split(df$x,id)
trait_list <- split(df$trait,id)

## Rescale time interval to [-1,1]
df$x_new <-  df$x*(1/13)-1
x_new <- split(df$x_new, id)
df$logtrait <- log10(df$trait)
log_trait_list <- split(df$logtrait,id)

## Test on a smaller dataset, containing the first 3 subjects
test_df <- df[1:25,]

## construct basis function for the first three subjects
age_test_list <- split(test_df$x_new, test_df$id)

## Legendre polynomial basis

library(sommer)

legbasis3_test <- list()
for(i in 1:3){
  legbasis3_test[[i]] <- leg(age_test_list[[i]], n = 3)
}

test_basis <- do.call(rbind, legbasis3_test)
test_df <- cbind(test_df, test_basis)

## Calculate the relation matrix A
test_pede <- editPed(sire = as.character(c(1,1,1)), dam = as.character(c(101,101,101)), 
                     label = as.character(c(10001,10002,10003)))
test_ped <- with(test_pede, pedigree(label=label,sire = sire, dam=dam))
test_A <- as.matrix(getA(test_ped)[3:5,3:5])

# lmm formula
test_fform <- logtrait ~ -1 + test_df$leg0 + test_df$leg1 + test_df$leg2 + test_df$leg3 +
  (-1 + test_df$leg0 + test_df$leg1 + test_df$leg2 + test_df$leg3 | test_df$id) + # genetic random effect, use 3rd order
  (-1 + test_df$leg0 + test_df$leg1 + test_df$leg2 | test_df$id) # environmental random effect, 2nd order
test_parsedFormula <- lFormula(formula = test_fform, data = test_df)

## manually update matrix ZG
I4 <- as(diag(4), "sparseMatrix")
test_LA <- as(t(chol(test_A)), "sparseMatrix")
test_M <- kronecker(test_LA, I4)

test_Z <- t(test_parsedFormula$reTrms$Zt)
test_Z_standard <- as.matrix(test_Z)
test_Z[,1:dim(test_M)[1]] <- test_Z[,1:dim(test_M)[1]] %*% test_M
ZG <- as.matrix(test_Z[,1:dim(test_M)[1]])
ZE <- as.matrix((test_Z[,13:21]))

######################################################################################

## use genetic FMEM

test_ff <- fit_genetic_fmm2(formula = test_fform, data = test_df, A = test_A, nbasis = c(4,3))
Zt <- test_ff@pp$`.->Zt`
Z <- t(Zt)
ZG_model <- as.matrix(Z[,1:12])
ZE_model <- as.matrix(Z[,13:21])

all(ZG == ZG_model) # TRUE
all(ZE == ZE_model) # TRUE