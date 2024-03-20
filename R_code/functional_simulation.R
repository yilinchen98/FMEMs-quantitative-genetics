library(Matrix)
library(MASS)
library(fda)
library(lme4)
library(pedigreemm)
library(mvnfast)
library(ggplot2)
library(plotly)
#library(clusterGeneration)


# Functional Simulation Study
## Y = X * beta + Z^G * alpha + Z^E * gamma + epsilon
## Step 1: Fix a functional basis 
## Here we choose B-spline for our simulation
basisObj <- create.bspline.basis(c(0,1), nbasis = 5, norder = 4)

## Step 2: Set the distribution of the random effect vectors
### Calculate the genetic relationship matrix A (use the TC dataset)is
setwd("D:/KCL_2023-2027_PhD/Year 1/Genetics/FMEMs-quantitative-genetics/R_code")
TRFUN25PUP4 = read.delim("TRFUN25PUP4.DAT",header = FALSE)
names(TRFUN25PUP4)<-c("id","sire","dam","trait","x")
df <- data.frame(TRFUN25PUP4)

FirstUniqueIdPos <- which(duplicated(df$id) == FALSE)
N = length(FirstUniqueIdPos) # N = 873 subjects
n = 8730 # each individual have 10 observations, in total 8730

pos = df$id[FirstUniqueIdPos] # extract ids for all subjects
sire_id = df$sire[FirstUniqueIdPos] # extract ids for sire
dam_id = df$dam[FirstUniqueIdPos] # extract ids for dam

pede <- editPed(sire = sire_id, dam = dam_id, label = pos)
ped<- with(pede, pedigree(label=label, sire=sire, dam=dam))
A <- getA(ped)[163:1035,163:1035]

### genetic covariance matrix
C_gen <- matrix(c(750, 10 ,130, 80, 250,
              10, 800, 30, 15, 40,
              130, 30, 700, 50, 130,
              80, 15, 50, 420, 50,
              250, 40, 130, 50, 330), nrow = 5, byrow = T)
#C_gen <- 100*as.matrix(genPositiveDefMat(dim = 5, covMethod = "eigen", eigenvalue = c(10, 8, 6, 4, 2),  rangeVar = c(1, 10))$Sigma)

### environmental covariance matrix
C_env <- C_gen

### reparameterised genetic covariance 
C_gen_para <- as(kronecker(A, C_gen), "dgCMatrix") 

### reparameterised environmental covariance 
I_N <- as(diag(N), "dgCMatrix")
C_env_para <- as(kronecker(I_N, C_env), "dgCMatrix") 

## Plot the true covariances functions
time_rang <- seq(0,1,length=10)
basis <- eval.basis(time_rang, basisObj)

C_gen_fun <- basis %*% C_gen %*% t(basis)
C_env_fun <- basis %*% C_env %*% t(basis)
P_true <- C_gen_fun + C_env_fun

fit4 <- plot_ly(x = time_rang, y = time_rang, z = ~C_gen_fun,
                type = 'surface', showscale = FALSE) %>% 
  layout(title = "TRUE Gen Variance", 
         scene = list(zaxis = list(title=""),
                      xaxis = list(title = "time"),
                      yaxis = list(title = "time")))
fit5 <- plot_ly(x = time_rang, y = time_rang, z = ~C_env_fun,
                type = 'surface', showscale = FALSE) %>% 
  layout(title = "TRUE Env Variance", 
         scene = list(zaxis = list(title=""),
                      xaxis = list(title = "time"),
                      yaxis = list(title = "time")))
fit6 <- plot_ly(x = time_rang, y = time_rang, z = ~P_true,
                type = 'surface', showscale = FALSE) %>% 
  layout(title = "TRUE Phen Variance", 
         scene = list(zaxis = list(title=""),
                      xaxis = list(title = "time"),
                      yaxis = list(title = "time")))

fit4
fit5
fit6

### residual variance
sigma2 <- 50

### distribution of the genetic random vector
mu <- rep(0, dim(C_gen_para)[1])
alpha <- rmvn(n=50, mu = mu, sigma = C_gen_para) #ith row is the ith simulated vector
gamma <- rmvn(n=50, mu = mu, sigma = C_env_para)

### distribution of the error vector
mu_res <- rep(0, n)
I_n <- as(diag(n), "dgCMatrix")
res_cov <- as(sigma2 * I_n, "dgCMatrix") 
  
res <- rmvn(n=50, mu=mu_res, sigma = res_cov) # error vector

## Step 3: Calculate the random effect design matrices
uniqueIds <- seq(1,873, length=873)
id <- rep(uniqueIds, each=10)
time_rang <- seq(0,1,length=10)
basis <- eval.basis(time_rang, basisObj)
b1 <- rep(basis[,1], times = 873)
b2 <- rep(basis[,2], times = 873)
b3 <- rep(basis[,3], times = 873)
b4 <- rep(basis[,4], times = 873)
b5 <- rep(basis[,5], times = 873)

y_pre <- rep(0, 8730) # dummy response values, will be updated 
df_pre <- data.frame(id = id, y_pre = y_pre, b1 = b1, b2 = b2, b3 = b3, b4 = b4, b5 = b5)

parsedFormula <- y_pre ~ (-1 + df_pre$b1 + df_pre$b2 + df_pre$b3 + df_pre$b4 + df_pre$b5 | df_pre$id) + 
  (-1 + df_pre$b1 + df_pre$b2 + df_pre$b3 + df_pre$b4 + df_pre$b5 | df_pre$id)

Z_pre <- t(lFormula(formula=parsedFormula, data=df_pre)$reTrms$Zt)

### update the genetic covariance matrix
L <- as(t(chol(A)), "dgCMatrix")
I3 <- as(diag(5), "dgCMatrix")
M <- kronecker(L, I3)
ZE <- Z_pre[,1:4365] 
ZG <- Z_pre[,1:4365] %*% M 
df_pre$y_pre <- as.vector(ZG %*% alpha[1,] + ZE %*% gamma[1,] + res[1,]) # random effects + residual

## Step 4: Calculate fixed effect
f <- lm(y_pre ~ -1 + df_pre$b1 + df_pre$b2 + df_pre$b3 + df_pre$b4 + df_pre$b5, data = df_pre)
fix_ef <- f$coefficients[1] * df_pre$b1 + f$coefficients[2] * df_pre$b2 + f$coefficients[3] * df_pre$b3 +
  f$coefficients[4] * df_pre$b4 + f$coefficients[5] + df_pre$b5

## Step 5: Generate response 
response <- as.vector(fix_ef + df_pre$y_pre) 
df_simu <- data.frame(id = id, response = response)

y_list <- split(response, df_simu$id)
plot(x=c(0,1), y=c(-200,200), type="n", xlab="time", ylab = "response")
for (i in 1:873){
  lines(time_rang, y_list[[i]], type = "l", col = i)
}

##################################################################################################################################

# Let's fit mixed-effect model using the simulated response
## Data smoothing
timefine <- seq(0,1,length=25)
y_hat <- matrix(0, 25, 873)
lambda <- rep(0, 873)
for (i in 1:873){
  ss <- smooth.spline(time_rang, y_list[[i]], cv = FALSE)
  y_hat[,i] <- predict(ss, timefine)$y
  lambda[i] <- ss$lambda
}

matplot(timefine, y_hat, col = 1:873, type= "l")

## FPCA
fpcaobj <- prcomp(x=t(y_hat), retx = TRUE, center = TRUE, rank. = 3)
pcs <- fpcaobj$rotation # eigen vectors as basis functions for moedel-fitting
summary(pcs)

par(mfrow=c(3,1))
for (i in 1:3) {
  plot(timefine, pcs[, i], type = "l", 
       xlab = "time", ylab = "",
       main = paste("Growth Curve Principal Component", i))
}

## Evaluate on a slightly dense time grid and reform the dataframe
subjectID <- rep(unique(df_simu$id), each=25)
response_test <- c(y_hat)
basis1 <- rep(pcs[,1], times = 873)
basis2 <- rep(pcs[,2], times = 873)
basis3 <- rep(pcs[,3], times = 873)
basis4 <- rep(pcs[,4], times = 873)
df_test <- data.frame(id = subjectID, y= response_test, phi1 = basis1, phi2 = basis2, phi3 = basis3)

## Fit FMEM
fform <- y ~ -1 + df_test$phi1 + df_test$phi2 + df_test$phi3 +  
  (-1 + df_test$phi1 + df_test$phi2 + df_test$phi3 | df_test$id) + 
  (-1 + df_test$phi1 + df_test$phi2 + df_test$phi3 | df_test$id)

system.time(
 ft <-  fit_genetic_fmm(fform, df_test, A, pcs)
) # user   system   elapsed
# 109.42   5.59     218.17 
summary(ft)

# Extract covariance function
## Extract the genetic and environment covariance matrices
vc <- VarCorr(ft)
CG <- vc[["df_test.id"]]
CE <- vc[["df_test.id.1"]]


### Convert to genetic covariance function
CG_fun <- pcs %*% CG %*% t(pcs)
### environmental covariance function
CE_fun <- pcs %*% CE %*% t(pcs)
### Phenotypic covariance function
P_fun <- CG_fun + CE_fun

# Plot the genetic covariance function
fig_RR1 <- plot_ly(x = timefine, y = timefine, z = ~CG_fun, 
                   type = "surface", showscale = FALSE) %>% 
  layout(title = "Genetic Variance Function", 
         scene = list(zaxis = list(title=""),
                      xaxis = list(title = "time"),
                      yaxis = list(title = "time")))

# Plot the environmental covariance function
fig_RR2 <- plot_ly(x = timefine, y = timefine, z = ~CE_fun, 
                   type = "surface", showscale = FALSE) %>% 
  layout(title = "Environment Variance Function", 
         scene = list(zaxis = list(title=""),
                      xaxis = list(title = "time"),
                      yaxis = list(title = "time")))

# Plot the phenotypic covariance function
fig_RR3 <- plot_ly(x = timefine, y = timefine, z = ~P_fun, 
                   type = "surface", showscale = FALSE) %>% 
  layout(title = "Phenotypic Variance Function", 
         scene = list(zaxis = list(title=""),
                      xaxis = list(title = "time"),
                      yaxis = list(title = "time")))

# Display each plot separately

fig_RR1
fig_RR2
fig_RR3

###############################################################################

## Let's use 4 PCs as basis functions to fit model
fpcaobj1 <- prcomp(x=t(y_hat), retx = TRUE, center = TRUE, rank. = 4)
pcs1 <- fpcaobj1$rotation # eigen vectors as basis functions for moedel-fitting
summary(pcs1)

par(mfrow=c(4,1))
for (i in 1:4) {
  plot(timefine, pcs1[, i], type = "l", 
       xlab = "time", ylab = "",
       main = paste("Growth Curve Principal Component", i))
}

bas1 <- rep(pcs1[,1], times = 873)
bas2 <- rep(pcs1[,2], times = 873)
bas3 <- rep(pcs1[,3], times = 873)
bas4 <- rep(pcs1[,4], times = 873)
df_test1 <- data.frame(id = subjectID, y= response_test, phi1 = bas1, phi2 = bas2, phi3 = bas3, phi4 = bas4)

fform1 <- y ~ -1 + df_test1$phi1 + df_test1$phi2 + df_test1$phi3 + df_test1$phi4 + 
  (-1 + df_test1$phi1 + df_test1$phi2 + df_test1$phi3 + df_test1$phi4 | df_test1$id) + 
  (-1 + df_test1$phi1 + df_test1$phi2 + df_test1$phi3 + df_test1$phi4 | df_test1$id)

system.time(
  ft1 <-  fit_genetic_fmm(fform1, df_test1, A, pcs1)
) # user   system   elapsed
# 109.42   5.59     218.17 
summary(ft1)

# Extract covariance function
## Extract the genetic and environment covariance matrices
vc1 <- VarCorr(ft1)
CG1 <- vc1[["df_test1.id"]]
CE1 <- vc1[["df_test1.id.1"]]

### Convert to genetic covariance function
CG_fun1 <- pcs1 %*% CG1 %*% t(pcs1)
### environmental covariance function
CE_fun1 <- pcs1 %*% CE1 %*% t(pcs1)
### Phenotypic covariance function
P_fun1 <- CG_fun1 + CE_fun1

# Plot the genetic covariance function
fig_RR4 <- plot_ly(x = timefine, y = timefine, z = ~CG_fun1, 
                   type = "surface", showscale = FALSE) %>% 
  layout(title = "Genetic Variance Function", 
         scene = list(zaxis = list(title=""),
                      xaxis = list(title = "time"),
                      yaxis = list(title = "time")))

# Plot the environmental covariance function
fig_RR5 <- plot_ly(x = timefine, y = timefine, z = ~CE_fun1, 
                   type = "surface", showscale = FALSE) %>% 
  layout(title = "Environment Variance Function", 
         scene = list(zaxis = list(title=""),
                      xaxis = list(title = "time"),
                      yaxis = list(title = "time")))

# Plot the phenotypic covariance function
fig_RR6 <- plot_ly(x = timefine, y = timefine, z = ~P_fun1, 
                   type = "surface", showscale = FALSE) %>% 
  layout(title = "Phenotypic Variance Function", 
         scene = list(zaxis = list(title=""),
                      xaxis = list(title = "time"),
                      yaxis = list(title = "time")))

# Display each plot separately

fig_RR4
fig_RR5
fig_RR6

## Let's use 2 PCs as basis functions to fit model
fpcaobj2 <- prcomp(x=t(y_hat), retx = TRUE, center = TRUE, rank. = 2)
pcs2 <- fpcaobj2$rotation # eigen vectors as basis functions for moedel-fitting
summary(pcs2)

bs1 <- rep(pcs2[,1], times = 873)
bs2 <- rep(pcs2[,2], times = 873)

df_test2 <- data.frame(id = subjectID, y= response_test, phi1 = bs1, phi2 = bs2)

fform2 <- y ~ -1 + df_test2$phi1 + df_test2$phi2 + 
  (-1 + df_test2$phi1 + df_test2$phi2 | df_test2$id) + 
  (-1 + df_test2$phi1 + df_test2$phi2 | df_test2$id)

system.time(
  ft2 <-  fit_genetic_fmm(fform2, df_test2, A, pcs2)
) 

summary(ft2)

# Extract covariance function
## Extract the genetic and environment covariance matrices
vc2 <- VarCorr(ft2)
CG2 <- vc2[["df_test2.id"]]
CE2 <- vc2[["df_test2.id.1"]]

### Convert to genetic covariance function
CG_fun2 <- pcs2 %*% CG2 %*% t(pcs2)
### environmental covariance function
CE_fun2 <- pcs2 %*% CE2 %*% t(pcs2)
### Phenotypic covariance function
P_fun2 <- CG_fun2 + CE_fun2

# Plot the genetic covariance function
fig_RR7 <- plot_ly(x = timefine, y = timefine, z = ~CG_fun2, 
                   type = "surface", showscale = FALSE) %>% 
  layout(title = "Genetic Variance Function", 
         scene = list(zaxis = list(title=""),
                      xaxis = list(title = "time"),
                      yaxis = list(title = "time")))

# Plot the environmental covariance function
fig_RR8 <- plot_ly(x = timefine, y = timefine, z = ~CE_fun2, 
                   type = "surface", showscale = FALSE) %>% 
  layout(title = "Environment Variance Function", 
         scene = list(zaxis = list(title=""),
                      xaxis = list(title = "time"),
                      yaxis = list(title = "time")))

# Plot the phenotypic covariance function
fig_RR9 <- plot_ly(x = timefine, y = timefine, z = ~P_fun2, 
                   type = "surface", showscale = FALSE) %>% 
  layout(title = "Phenotypic Variance Function", 
         scene = list(zaxis = list(title=""),
                      xaxis = list(title = "time"),
                      yaxis = list(title = "time")))

# Display each plot separately

fig_RR7
fig_RR8
fig_RR9

