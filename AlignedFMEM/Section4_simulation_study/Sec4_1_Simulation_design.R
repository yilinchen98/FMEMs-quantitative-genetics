source("PrepareR.R")
A <- readRDS("A.rds") # relationship matrix
fixed_effect <- readRDS("fixed_effect.rds")

N <- 873 # total number of individuals
n <- 14 # number of measurements per individual
nbasis <- 5 # number of basis
npc <- 6
ngroups <- 50 # generate 50 groups of aligned curves

time_rang <- seq(0,1,length=n) # time points t_j
timefine <- seq(0,1, length = 100)

basisObj <- create.bspline.basis(c(0,1), nbasis = nbasis, norder = 4)
basis <- eval.basis(time_rang, basisObj)

### true covariance matrix
C_true <- matrix(c(750, 10 ,130, 80, 250,
                    10, 800, 30, 15, 40,
                    130, 30, 700, 50, 130,
                    80, 15, 50, 420, 50,
                    250, 40, 130, 50, 330), nrow = 5, byrow = T)

### residual variance
sigma2 <- 15

set.seed(123)

## generate data
C_gen_para <- as(kronecker(A, C_true), "sparseMatrix") 

### reparameterised environmental covariance 
I_N <- as(diag(N), "sparseMatrix")
C_env_para <- as(kronecker(I_N, C_true), "sparseMatrix")

### distribution of genetic effect
mu <- rep(0, dim(C_gen_para)[1])
alpha <- rmvn(n=ngroups, mu = mu, sigma = C_gen_para)

### distribution of environmental effect
gamma <- rmvn(n=ngroups, mu=mu, sigma = C_env_para)

### distribution of the error vector
mu_res <- rep(0, N*n)
I_n <- as(diag(N*n), "sparseMatrix")
res_cov <- as(sigma2 * I_n, "sparseMatrix") 

res <- rmvn(n=ngroups, mu=mu_res, sigma = res_cov)# error vector

## generating curve data and fit data to mixed-effect model
gpf <- gl(N,n) # grouping factor
Ji <- t(as(gpf, Class = "sparseMatrix")) # Indicator matrix of grouping factor indices

Xi <- as(cbind(rep(basis[,1], times = N),rep(basis[,2], times = N),rep(basis[,3], times = N),rep(basis[,4], times = N),rep(basis[,5], times = N)), Class = "sparseMatrix" )# raw random effect matrix

Zi <- t(KhatriRao(t(Ji), t(Xi))) # random effect design matrix

Y_rand<- matrix(0, N*n, ngroups)
for (k in 1:ngroups){
  Yk <- Zi %*% alpha[k,] + Zi %*% gamma[k,] + res[k,]
  Y_rand[,k] <- Yk@x
}

f_true <- convert_to_basisfunctions(timefine, fixed_effect, time_rang)
Y <- Y_rand + rep(f_true, times=N) # simulated data (fixed + random)
saveRDS(Y, file = "sim_data.rds")
