# Functional simulation result

## Simulation setting: n = 7, 10, 14 measurement points for each individual 
## N = 873 subjects
## res var = 15

## set up true covariance
N <- 873 # total number of individualsl
nbasis <- 5 # number of basis
ngroups <- 50

### genetic and environmental covariance matrix
C_true <- matrix(c(750, 10 ,130, 80, 250,
                   10, 800, 30, 15, 40,
                   130, 30, 700, 50, 130,
                   80, 15, 50, 420, 50,
                   250, 40, 130, 50, 330), nrow = 5, byrow = T)

### model reform
t7 <- seq(0,1,length=7) # time points t_j
t14 <- seq(0,1,length=14) # time points t_j
timefine <- seq(0,1, length =100)

basisObj <- create.bspline.basis(c(0,1), nbasis = nbasis, norder = 4)
basis7 <- eval.basis(t7,basisObj)
basis14 <- eval.basis(t14,basisObj)

C_fun_true7 <- basis7 %*% C_true %*% t(basis7)
C_fun_true14 <- basis14 %*% C_true %*% t(basis14)

# load result
setwd("D:/KCL_2023-2027_PhD/Genetics_FMEM_Project/FMEMs-quantitative-genetics/R_code/Functional_simulation/Simulation_PCs/PCs_FitRawData/simulation_results")
files <- list.files(pattern = "\\.Rdata$")

for (file in files) {
  load(file)
}

load("D:/KCL_2023-2027_PhD/Genetics_FMEM_Project/FMEMs-quantitative-genetics/R_code/Functional_simulation/Simulation_PCs/PCs_FitRawData/simulated_data/TC_fixed_effect.Rdata")
f_true7 <- convert_to_basisfunctions(timefine, fixed_effect, t7)
f_true14 <- convert_to_basisfunctions(timefine, fixed_effect, t14)

## Errors on FE

ferr_5705 <- rep(0, ngroups)
ferr_5715 <- rep(0, ngroups)
ferr_5725 <- rep(0, ngroups)

ferr_51405 <- rep(0, ngroups)
ferr_51415 <- rep(0, ngroups)
ferr_51425 <- rep(0, ngroups)

for (k in 1:ngroups){
  ferr_5705[k] <- sqrt(sum((f_true7 - fef_5pc705[,k])^2)/7) # root-mean-square errors
  ferr_5715[k] <- sqrt(sum((f_true7 - fef_5pc715[,k])^2)/7) # root-mean-square errors
  ferr_5725[k] <- sqrt(sum((f_true7 - fef_5pc725[,k])^2)/7) # root-mean-square errors
  
  ferr_51405[k] <- sqrt(sum((f_true14 - fef_5pc1405[,k])^2)/14)
  ferr_51415[k] <- sqrt(sum((f_true14 - fef_5pc1415[,k])^2)/14)
  ferr_51425[k] <- sqrt(sum((f_true14 - fef_5pc1425[,k])^2)/14)
}
par(mfrow = c(1, 1), bty = "l")
FefErr_list7 <- list(ferr_5705, ferr_5715, ferr_5725)
# Generate the boxplot
boxplot(FefErr_list7,
        xlab = "",
        ylab = "RMSE",
        ylim = c(1,8),
        names = c(TeX("$\\sigma^2 = 5$"), TeX("$\\sigma^2 = 15$"), TeX("$\\sigma^2 = 25$")))
mtext("Sampling points per subject = 7", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

FefErr_list14 <- list(ferr_51405, ferr_51415, ferr_51425)
# Generate the boxplot
boxplot(FefErr_list14,
        xlab = "",
        ylab = "RMSE",
        ylim = c(1,8),
        names = c(TeX("$\\sigma^2 = 5$"), TeX("$\\sigma^2 = 15$"), TeX("$\\sigma^2 = 25$")))
mtext("Sampling points per subject = 14", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

###################################################

## Error on the REs

GenErr_5705 <- rep(0, ngroups)
GenErr_5715 <- rep(0, ngroups)
GenErr_5725 <- rep(0, ngroups)

EnvErr_5705 <- rep(0, ngroups)
EnvErr_5715 <- rep(0, ngroups)
EnvErr_5725 <- rep(0, ngroups)

GenErr_51405 <- rep(0, ngroups)
GenErr_51415 <- rep(0, ngroups)
GenErr_51425 <- rep(0, ngroups)

EnvErr_51405 <- rep(0, ngroups)
EnvErr_51415 <- rep(0, ngroups)
EnvErr_51425 <- rep(0, ngroups)

for(k in 1:ngroups){
  GenErr_5705[k] <-norm((C_fun_true7 - CG_fun_5pc705[,,k]), type = "M") # max absolute error
  GenErr_5715[k] <-norm((C_fun_true7 - CG_fun_5pc715[,,k]), type = "M")
  GenErr_5725[k] <-norm((C_fun_true7 - CG_fun_5pc725[,,k]), type = "M")
  
  EnvErr_5705[k] <-norm((C_fun_true7 - CE_fun_5pc705[,,k]), type = "M")
  EnvErr_5715[k] <-norm((C_fun_true7 - CE_fun_5pc715[,,k]), type = "M")
  EnvErr_5725[k] <-norm((C_fun_true7 - CE_fun_5pc725[,,k]), type = "M")
  
  GenErr_51405[k] <-norm((C_fun_true14 - CG_fun_5pc1405[,,k]), type = "M") # max absolute error
  GenErr_51415[k] <-norm((C_fun_true14 - CG_fun_5pc1415[,,k]), type = "M")
  GenErr_51425[k] <-norm((C_fun_true14 - CG_fun_5pc1425[,,k]), type = "M")
  
  EnvErr_51405[k] <-norm((C_fun_true14 - CE_fun_5pc1405[,,k]), type = "M")
  EnvErr_51415[k] <-norm((C_fun_true14 - CE_fun_5pc1415[,,k]), type = "M")
  EnvErr_51425[k] <-norm((C_fun_true14 - CE_fun_5pc1425[,,k]), type = "M")
}

CovErr_list7 <- list(GenErr_5705, EnvErr_5705, GenErr_5715, EnvErr_5715, GenErr_5725, EnvErr_5725)

# Positions for the boxplots
positions <- c(1, 2, 3.5, 4.5, 6, 7)

# Group labels (midpoints)
group_labels <- c(mean(c(1, 2)), mean(c(3.5, 4.5)), mean(c(6, 7)))

# Names for x-axis labels
names <- c(TeX("$\\sigma^2 = 5$"), "", TeX("$\\sigma^2 = 15$"), "", TeX("$\\sigma^2 = 25$"), "")

# Generate the boxplot
boxplot(CovErr_list7,
        at = positions,
        names = names,
        xlab = "",
        ylab = "Error (sup metric)",
        ylim = c(0,1100),
        col = rep(c("salmon", "lightblue"), 3),
        xaxt = 'n')  # Disable default x-axis

# Add custom x-axis labels
axis(1, at = group_labels, labels = c(TeX("$\\sigma^2 = 5$"), TeX("$\\sigma^2 = 15$"), TeX("$\\sigma^2 = 25$")))

# Add legend
legend("topright", legend = c("Genetic", "Environmental"), fill = c("salmon", "lightblue"), bty = "n")
mtext("Sampling points per subject = 7", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

#############################################################################################################
CovErr_list14 <- list(GenErr_51405, EnvErr_51405, GenErr_51415, EnvErr_51415, GenErr_51425, EnvErr_51425)

# Positions for the boxplots
positions <- c(1, 2, 3.5, 4.5, 6, 7)

# Group labels (midpoints)
group_labels <- c(mean(c(1, 2)), mean(c(3.5, 4.5)), mean(c(6, 7)))

# Names for x-axis labels
names <- c(TeX("$\\sigma^2 = 5$"), "", TeX("$\\sigma^2 = 15$"), "", TeX("$\\sigma^2 = 25$"), "")

# Generate the boxplot
boxplot(CovErr_list14,
        at = positions,
        names = names,
        xlab = "",
        ylab = "Error (sup metric)",
        ylim = c(0,1100),
        col = rep(c("salmon", "lightblue"), 3),
        xaxt = 'n')  # Disable default x-axis

# Add custom x-axis labels
axis(1, at = group_labels, labels = c(TeX("$\\sigma^2 = 5$"), TeX("$\\sigma^2 = 15$"), TeX("$\\sigma^2 = 25$")))

# Add legend
legend("topright", legend = c("Genetic", "Environmental"), fill = c("salmon", "lightblue"), bty = "n")
mtext("Sampling points per subject = 14", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")


#######################################################################################
# Functional Boxplot

par(mfrow = c(3,3))
### FE
FE5705 <- fda::fbplot(fef_5pc705, x = t7, xlab = "Time", ylab = "", xlim = c(0,1),
                      ylim = c(-5,300), color = "grey", barcol = "lightblue")
lines(t7,f_true7, lty = "solid", col = "red")

FE51005 <- fda::fbplot(fef_5pc1005, x = t10, xlab = "Time", ylab = "", xlim = c(0,1),
                      ylim = c(-5,300), color = "grey", barcol = "lightblue")
lines(t10,f_true10, lty = "solid", col = "red")

FE51405 <- fda::fbplot(fef_5pc1405, x = t14, xlab = "Time", ylab = "", xlim = c(0,1),
                       ylim = c(-5,300), color = "grey", barcol = "lightblue")
lines(t14,f_true14, lty = "solid", col = "red")

FE5715 <- fda::fbplot(fef_5pc715, x = t7, xlab = "Time", ylab = "", xlim = c(0,1),
                      ylim = c(-5,300), color = "grey", barcol = "lightblue")
lines(t7,f_true7, lty = "solid", col = "red")

FE51015 <- fda::fbplot(fef_5pc1015, x = t10, xlab = "Time", ylab = "", xlim = c(0,1),
                       ylim = c(-5,300), color = "grey", barcol = "lightblue")
lines(t10,f_true10, lty = "solid", col = "red")

FE51415 <- fda::fbplot(fef_5pc1415, x = t14, xlab = "Time", ylab = "", xlim = c(0,1),
                      ylim = c(-5,300), color = "grey", barcol = "lightblue")
lines(t14,f_true14, lty = "solid", col = "red")

FE5725 <- fda::fbplot(fef_5pc725, x = t7, xlab = "Time", ylab = "", xlim = c(0,1),
                      ylim = c(-5,300), color = "grey", barcol = "lightblue")
lines(t7,f_true7, lty = "solid", col = "red")

FE51025 <- fda::fbplot(fef_5pc1025, x = t10, xlab = "Time", ylab = "", xlim = c(0,1),
                       ylim = c(-5,300), color = "grey", barcol = "lightblue")
lines(t10,f_true10, lty = "solid", col = "red")

FE51425 <- fda::fbplot(fef_5pc1425, x = t14, xlab = "Time", ylab = "", xlim = c(0,1),
                       ylim = c(-5,300), color = "grey", barcol = "lightblue")
lines(t14,f_true14, lty = "solid", col = "red")

## Genetic and Environemental Eigenfunction
eigen_true7 <- eigen(C_fun_true7)$vectors[,1]
eigen_true10 <- eigen(C_fun_true10)$vectors[,1]
eigen_true14 <- eigen(C_fun_true14)$vectors[,1]

Geigen5705 <- matrix(0,7,ngroups)
Geigen5715 <- matrix(0,7,ngroups)
Geigen5725 <- matrix(0,7,ngroups)

Eeigen5705 <- matrix(0,7,ngroups)
Eeigen5715 <- matrix(0,7,ngroups)
Eeigen5725 <- matrix(0,7,ngroups)

Geigen51005 <- matrix(0,10,ngroups)
Geigen51015 <- matrix(0,10,ngroups)
Geigen51025 <- matrix(0,10,ngroups)

Eeigen51005 <- matrix(0,10,ngroups)
Eeigen51015 <- matrix(0,10,ngroups)
Eeigen51025 <- matrix(0,10,ngroups)

Geigen51405 <- matrix(0,14,ngroups)
Geigen51415 <- matrix(0,14,ngroups)
Geigen51425 <- matrix(0,14,ngroups)

Eeigen51405 <- matrix(0,14,ngroups)
Eeigen51415 <- matrix(0,14,ngroups)
Eeigen51425 <- matrix(0,14,ngroups)

for (i in 1:ngroups){
  Geigen5705[,i] <- eigen(CG_fun_5pc705[,,i])$vectors[,1]
  Geigen5715[,i] <- eigen(CG_fun_5pc715[,,i])$vectors[,1]
  Geigen5725[,i] <- eigen(CG_fun_5pc725[,,i])$vectors[,1]
  
  Eeigen5705[,i] <- eigen(CE_fun_5pc705[,,i])$vectors[,1]
  Eeigen5715[,i] <- eigen(CE_fun_5pc715[,,i])$vectors[,1]
  Eeigen5725[,i] <- eigen(CE_fun_5pc725[,,i])$vectors[,1]
  
  Geigen51005[,i] <- eigen(CG_fun_5pc1005[,,i])$vectors[,1]
  Geigen51015[,i] <- eigen(CG_fun_5pc1015[,,i])$vectors[,1]
  Geigen51025[,i] <- eigen(CG_fun_5pc1025[,,i])$vectors[,1]
  
  Eeigen51005[,i] <- eigen(CE_fun_5pc1005[,,i])$vectors[,1]
  Eeigen51015[,i] <- eigen(CE_fun_5pc1015[,,i])$vectors[,1]
  Eeigen51025[,i] <- eigen(CE_fun_5pc1025[,,i])$vectors[,1]
  
  Geigen51405[,i] <- eigen(CG_fun_5pc1405[,,i])$vectors[,1]
  Geigen51415[,i] <- eigen(CG_fun_5pc1415[,,i])$vectors[,1]
  Geigen51425[,i] <- eigen(CG_fun_5pc1425[,,i])$vectors[,1]
  
  Eeigen51405[,i] <- eigen(CE_fun_5pc1405[,,i])$vectors[,1]
  Eeigen51415[,i] <- eigen(CE_fun_5pc1415[,,i])$vectors[,1]
  Eeigen51425[,i] <- eigen(CE_fun_5pc1425[,,i])$vectors[,1]
}

## Genetic eigenfunction
GE5705 <- fda::fbplot(Geigen5705, x = t7, xlab = "Time", ylab = "", 
                      xlim = c(0,1), ylim = c(-1,1),
                      color = "grey", barcol = "lightblue",outliercol="orange")
lines(t7,eigen_true7, lty = "solid", col = "red")

GE5715 <- fda::fbplot(Geigen5715, x = t7, xlab = "Time", ylab = "", 
                      xlim = c(0,1), ylim = c(-1,1),
                      color = "grey", barcol = "lightblue",outliercol="orange")
lines(t7,eigen_true7, lty = "solid", col = "red")

GE5725 <- fda::fbplot(Geigen5725, x = t7, xlab = "Time", ylab = "", 
                      xlim = c(0,1), ylim = c(-1,1),
                      color = "grey", barcol = "lightblue",outliercol="orange")
lines(t7,eigen_true7, lty = "solid", col = "red")

GE51005 <- fda::fbplot(Geigen51005, x = t10, xlab = "Time", ylab = "", 
                      xlim = c(0,1), ylim = c(-1,1),
                      color = "grey", barcol = "lightblue",outliercol="orange")
lines(t10,eigen_true10, lty = "solid", col = "red")

GE51015 <- fda::fbplot(Geigen51015, x = t10, xlab = "Time", ylab = "", 
                      xlim = c(0,1), ylim = c(-1,1),
                      color = "grey", barcol = "lightblue",outliercol="orange")
lines(t10,eigen_true10, lty = "solid", col = "red")

GE51025 <- fda::fbplot(Geigen51025, x = t10, xlab = "Time", ylab = "", 
                      xlim = c(0,1), ylim = c(-1,1),
                      color = "grey", barcol = "lightblue",outliercol="orange")
lines(t10,eigen_true10, lty = "solid", col = "red")

GE51405 <- fda::fbplot(Geigen51405, x = t14, xlab = "Time", ylab = "", 
                       xlim = c(0,1), ylim = c(-1,1),
                       color = "grey", barcol = "lightblue",outliercol="orange")
lines(t14,eigen_true14, lty = "solid", col = "red")

GE51415 <- fda::fbplot(Geigen51415, x = t14, xlab = "Time", ylab = "", 
                       xlim = c(0,1), ylim = c(-1,1),
                       color = "grey", barcol = "lightblue",outliercol="orange")
lines(t14,eigen_true14, lty = "solid", col = "red")

GE51425 <- fda::fbplot(Geigen51425, x = t14, xlab = "Time", ylab = "", 
                       xlim = c(0,1), ylim = c(-1,1),
                       color = "grey", barcol = "lightblue",outliercol="orange")
lines(t14,eigen_true14, lty = "solid", col = "red")

###############################################################
## Environmental eigenfunction
EE5705 <- fda::fbplot(Eeigen5705, x = t7, xlab = "Time", ylab = "", 
                      xlim = c(0,1), ylim = c(-1,1),
                      color = "grey", barcol = "lightblue",outliercol="orange")
lines(t7,eigen_true7, lty = "solid", col = "red")

EE5715 <- fda::fbplot(Eeigen5715, x = t7, xlab = "Time", ylab = "", 
                      xlim = c(0,1), ylim = c(-1,1),
                      color = "grey", barcol = "lightblue",outliercol="orange")
lines(t7,eigen_true7, lty = "solid", col = "red")

EE5725 <- fda::fbplot(Eeigen5725, x = t7, xlab = "Time", ylab = "", 
                      xlim = c(0,1), ylim = c(-1,1),
                      color = "grey", barcol = "lightblue",outliercol="orange")
lines(t7,eigen_true7, lty = "solid", col = "red")

EE51005 <- fda::fbplot(Eeigen51005, x = t10, xlab = "Time", ylab = "", 
                       xlim = c(0,1), ylim = c(-1,1),
                       color = "grey", barcol = "lightblue",outliercol="orange")
lines(t10,eigen_true10, lty = "solid", col = "red")

EE51015 <- fda::fbplot(Eeigen51015, x = t10, xlab = "Time", ylab = "", 
                       xlim = c(0,1), ylim = c(-1,1),
                       color = "grey", barcol = "lightblue",outliercol="orange")
lines(t10,eigen_true10, lty = "solid", col = "red")

EE51025 <- fda::fbplot(Eeigen51025, x = t10, xlab = "Time", ylab = "", 
                       xlim = c(0,1), ylim = c(-1,1),
                       color = "grey", barcol = "lightblue",outliercol="orange")
lines(t10,eigen_true10, lty = "solid", col = "red")

EE51405 <- fda::fbplot(Eeigen51405, x = t14, xlab = "Time", ylab = "", 
                       xlim = c(0,1), ylim = c(-1,1),
                       color = "grey", barcol = "lightblue",outliercol="orange")
lines(t14,eigen_true14, lty = "solid", col = "red")

EE51415 <- fda::fbplot(Eeigen51415, x = t14, xlab = "Time", ylab = "", 
                       xlim = c(0,1), ylim = c(-1,1),
                       color = "grey", barcol = "lightblue",outliercol="orange")
lines(t14,eigen_true14, lty = "solid", col = "red")

GE51425 <- fda::fbplot(Geigen51425, x = t14, xlab = "Time", ylab = "", 
                       xlim = c(0,1), ylim = c(-1,1),
                       color = "grey", barcol = "lightblue",outliercol="orange")
lines(t14,eigen_true14, lty = "solid", col = "red")
#################################################################################

## Simultaneous confidence band
## FE
mean_FE5705 <- rowMeans(fef_5pc705)
mean_FE5715 <- rowMeans(fef_5pc715)
mean_FE5725 <- rowMeans(fef_5pc725)

mean_FE51005 <- rowMeans(fef_5pc1005)
mean_FE51015 <- rowMeans(fef_5pc1015)
mean_FE51025 <- rowMeans(fef_5pc1025)

mean_FE51405 <- rowMeans(fef_5pc1405)
mean_FE51415 <- rowMeans(fef_5pc1415)
mean_FE51425 <- rowMeans(fef_5pc1425)

FE_5705_diff <- rep(0,ngroups)
FE_5715_diff <- rep(0,ngroups)
FE_5725_diff <- rep(0,ngroups)

FE_51005_diff <- rep(0,ngroups)
FE_51015_diff <- rep(0,ngroups)
FE_51025_diff <- rep(0,ngroups)

FE_51405_diff <- rep(0,ngroups)
FE_51415_diff <- rep(0,ngroups)
FE_51425_diff <- rep(0,ngroups)

for (i in 1:ngroups){
  FE_5705_diff[i] <- max(abs(fef_5pc705[,i] - mean_FE5705))
  FE_5715_diff[i] <- max(abs(fef_5pc715[,i] - mean_FE5715))
  FE_5725_diff[i] <- max(abs(fef_5pc725[,i] - mean_FE5725))
  
  FE_51005_diff[i] <- max(abs(fef_5pc1005[,i] - mean_FE51005))
  FE_51015_diff[i] <- max(abs(fef_5pc1015[,i] - mean_FE51015))
  FE_51025_diff[i] <- max(abs(fef_5pc1025[,i] - mean_FE51025))
  
  FE_51405_diff[i] <- max(abs(fef_5pc1405[,i] - mean_FE51405))
  FE_51415_diff[i] <- max(abs(fef_5pc1415[,i] - mean_FE51415))
  FE_51425_diff[i] <- max(abs(fef_5pc1425[,i] - mean_FE51425))
}

FE_5705_sup <- quantile(FE_5705_diff, probs = 0.95) 
FE_5715_sup <- quantile(FE_5715_diff, probs = 0.95)
FE_5725_sup <- quantile(FE_5725_diff, probs = 0.95)

FE_51005_sup <- quantile(FE_51005_diff, probs = 0.95) 
FE_51015_sup <- quantile(FE_51015_diff, probs = 0.95)
FE_51025_sup <- quantile(FE_51025_diff, probs = 0.95)

FE_51405_sup <- quantile(FE_51405_diff, probs = 0.95) 
FE_51415_sup <- quantile(FE_51415_diff, probs = 0.95)
FE_51425_sup <- quantile(FE_51425_diff, probs = 0.95)

###############################################################

## Random effect
mean_GE5705 <- rowMeans(Geigen5705)
mean_GE5715 <- rowMeans(Geigen5715)
mean_GE5725 <- rowMeans(Geigen5725)

mean_EE5705 <- rowMeans(Eeigen5705)
mean_EE5715 <- rowMeans(Eeigen5715)
mean_EE5725 <- rowMeans(Eeigen5725)

mean_GE51405 <- rowMeans(Geigen51405)
mean_GE51415 <- rowMeans(Geigen51415)
mean_GE51425 <- rowMeans(Geigen51425)

mean_EE51405 <- rowMeans(Eeigen51405)
mean_EE51415 <- rowMeans(Eeigen51415)
mean_EE51425 <- rowMeans(Eeigen51425)

### Genetic Eigen function

Geigen_diff_5705 <- rep(0,ngroups)
Geigen_diff_5715 <- rep(0,ngroups)
Geigen_diff_5725 <- rep(0,ngroups)

Eeigen_diff_5705 <- rep(0,ngroups)
Eeigen_diff_5715 <- rep(0,ngroups)
Eeigen_diff_5725 <- rep(0,ngroups)

Geigen_diff_51405 <- rep(0,ngroups)
Geigen_diff_51415 <- rep(0,ngroups)
Geigen_diff_51425 <- rep(0,ngroups)

Eeigen_diff_51405 <- rep(0,ngroups)
Eeigen_diff_51415 <- rep(0,ngroups)
Eeigen_diff_51425 <- rep(0,ngroups)

for (i in 1:ngroups){
  Geigen_diff_5705[i] <-  max(abs(Geigen5705[,i] - mean_GE5705))
  Geigen_diff_5715[i] <-  max(abs(Geigen5715[,i] - mean_GE5715))
  Geigen_diff_5725[i] <-  max(abs(Geigen5725[,i] - mean_GE5725))
  
  Eeigen_diff_5705[i] <-  max(abs(Eeigen5705[,i] - mean_EE5705))
  Eeigen_diff_5715[i] <-  max(abs(Eeigen5715[,i] - mean_EE5715))
  Eeigen_diff_5725[i] <-  max(abs(Eeigen5725[,i] - mean_EE5725))
  
  Geigen_diff_51405[i] <-  max(abs(Geigen51405[,i] - mean_GE51405))
  Geigen_diff_51415[i] <-  max(abs(Geigen51415[,i] - mean_GE51415))
  Geigen_diff_51425[i] <-  max(abs(Geigen51425[,i] - mean_GE51425))
  
  Eeigen_diff_51405[i] <-  max(abs(Eeigen51405[,i] - mean_EE51405))
  Eeigen_diff_51415[i] <-  max(abs(Eeigen51415[,i] - mean_EE51415))
  Eeigen_diff_51425[i] <-  max(abs(Eeigen51425[,i] - mean_EE51425))
}

Geigen_5705_sup <- quantile(Geigen_diff_5705, probs = 0.95) 
Geigen_5715_sup <- quantile(Geigen_diff_5715, probs = 0.95) 
Geigen_5725_sup <- quantile(Geigen_diff_5725, probs = 0.95)

Eeigen_5705_sup <- quantile(Eeigen_diff_5705, probs = 0.95) 
Eeigen_5715_sup <- quantile(Eeigen_diff_5715, probs = 0.95) 
Eeigen_5725_sup <- quantile(Eeigen_diff_5725, probs = 0.95)

Geigen_51405_sup <- quantile(Geigen_diff_51405, probs = 0.95) 
Geigen_51415_sup <- quantile(Geigen_diff_51415, probs = 0.95) 
Geigen_51425_sup <- quantile(Geigen_diff_51425, probs = 0.95)

Eeigen_51405_sup <- quantile(Eeigen_diff_51405, probs = 0.95) 
Eeigen_51415_sup <- quantile(Eeigen_diff_51415, probs = 0.95) 
Eeigen_51425_sup <- quantile(Eeigen_diff_51425, probs = 0.95)

###############################################################

mean_CG5705 <- apply(CG_fun_5pc705, c(1, 2), mean)
mean_CG5715 <- apply(CG_fun_5pc715, c(1, 2), mean)
mean_CG5725 <- apply(CG_fun_5pc725, c(1, 2), mean)

mean_CE5705 <- apply(CE_fun_5pc705, c(1, 2), mean)
mean_CE5715 <- apply(CE_fun_5pc715, c(1, 2), mean)
mean_CE5725 <- apply(CE_fun_5pc725, c(1, 2), mean)

mean_CG51405 <- apply(CG_fun_5pc1405, c(1, 2), mean)
mean_CG51415 <- apply(CG_fun_5pc1415, c(1, 2), mean)
mean_CG51425 <- apply(CG_fun_5pc1425, c(1, 2), mean)

mean_CE51405 <- apply(CE_fun_5pc1405, c(1, 2), mean)
mean_CE51415 <- apply(CE_fun_5pc1415, c(1, 2), mean)
mean_CE51425 <- apply(CE_fun_5pc1425, c(1, 2), mean)


CG_diff_5705 <- rep(0,ngroups)
CG_diff_5715 <- rep(0,ngroups)
CG_diff_5725 <- rep(0,ngroups)

CE_diff_5705 <- rep(0,ngroups)
CE_diff_5715 <- rep(0,ngroups)
CE_diff_5725 <- rep(0,ngroups)

CG_diff_51405 <- rep(0,ngroups)
CG_diff_51415 <- rep(0,ngroups)
CG_diff_51425 <- rep(0,ngroups)

CE_diff_51405 <- rep(0,ngroups)
CE_diff_51415 <- rep(0,ngroups)
CE_diff_51425 <- rep(0,ngroups)

for (i in 1:ngroups){
  CG_diff_5705[i] <- max(abs(CG_fun_5pc705[,,i] - mean_CG5705))
  CG_diff_5715[i] <- max(abs(CG_fun_5pc715[,,i] - mean_CG5715))
  CG_diff_5725[i] <- max(abs(CG_fun_5pc725[,,i] - mean_CG5725))
  
  CE_diff_5705[i] <- max(abs(CE_fun_5pc705[,,i] - mean_CE5705))
  CE_diff_5715[i] <- max(abs(CE_fun_5pc715[,,i] - mean_CE5715))
  CE_diff_5725[i] <- max(abs(CE_fun_5pc725[,,i] - mean_CE5725))
  
  CG_diff_51405[i] <- max(abs(CG_fun_5pc1405[,,i] - mean_CG51405))
  CG_diff_51415[i] <- max(abs(CG_fun_5pc1415[,,i] - mean_CG51415))
  CG_diff_51425[i] <- max(abs(CG_fun_5pc1425[,,i] - mean_CG51425))
  
  CE_diff_51405[i] <- max(abs(CE_fun_5pc1405[,,i] - mean_CE51405))
  CE_diff_51415[i] <- max(abs(CE_fun_5pc1415[,,i] - mean_CE51415))
  CE_diff_51425[i] <- max(abs(CE_fun_5pc1425[,,i] - mean_CE51425))
}

CG_5705_sup <- quantile(CG_diff_5705, probs = 0.95)
CG_5715_sup <- quantile(CG_diff_5715, probs = 0.95)
CG_5725_sup <- quantile(CG_diff_5725, probs = 0.95)

CE_5705_sup <- quantile(CE_diff_5705, probs = 0.95)
CE_5715_sup <- quantile(CE_diff_5715, probs = 0.95)
CE_5725_sup <- quantile(CE_diff_5725, probs = 0.95)

CG_51405_sup <- quantile(CG_diff_51405, probs = 0.95)
CG_51415_sup <- quantile(CG_diff_51415, probs = 0.95)
CG_51425_sup <- quantile(CG_diff_51425, probs = 0.95)

CE_51405_sup <- quantile(CE_diff_51405, probs = 0.95)
CE_51415_sup <- quantile(CE_diff_51415, probs = 0.95)
CE_51425_sup <- quantile(CE_diff_51425, probs = 0.95)
