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
files <- list.files(pattern = "\\.Rdata$")

for (file in files) {
  load(file)
}
load("TC_fixed_effects.Rdata")
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

### FE
FE5705 <- fda::fbplot(fef_5pc705, x = t7, xlab = "Time", ylab = "", xlim = c(0,1),
                      ylim = c(-5,300), color = "grey", barcol = "lightblue")
lines(t7,f_true7, lty = "solid", col = "red")

FE51405 <- fda::fbplot(fef_5pc1405, x = t14, xlab = "Time", ylab = "", xlim = c(0,1),
                       ylim = c(-5,300), color = "grey", barcol = "lightblue")
lines(t14,f_true14, lty = "solid", col = "red")

FE5715 <- fda::fbplot(fef_5pc715, x = t7, xlab = "Time", ylab = "", xlim = c(0,1),
                      ylim = c(-5,300), color = "grey", barcol = "lightblue")
lines(t7,f_true7, lty = "solid", col = "red")
mtext("Sampling points per subject = 7", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

FE51415 <- fda::fbplot(fef_5pc1415, x = t14, xlab = "Time", ylab = "", xlim = c(0,1),
                      ylim = c(-5,300), color = "grey", barcol = "lightblue")
lines(t14,f_true14, lty = "solid", col = "red")
mtext("Sampling points per subject = 14", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

FE5725 <- fda::fbplot(fef_5pc725, x = t7, xlab = "Time", ylab = "", xlim = c(0,1),
                      ylim = c(-5,300), color = "grey", barcol = "lightblue")
lines(t7,f_true7, lty = "solid", col = "red")

FE51425 <- fda::fbplot(fef_5pc1425, x = t14, xlab = "Time", ylab = "", xlim = c(0,1),
                       ylim = c(-5,300), color = "grey", barcol = "lightblue")
lines(t14,f_true14, lty = "solid", col = "red")

## Genetic and Environemental Eigenfunction
eigen_true7 <- eigen(C_fun_true7)$vectors[,1]
eigen_true14 <- eigen(C_fun_true14)$vectors[,1]

Geigen5705 <- matrix(0,7,ngroups)
Geigen5715 <- matrix(0,7,ngroups)
Geigen5725 <- matrix(0,7,ngroups)

Eeigen5705 <- matrix(0,7,ngroups)
Eeigen5715 <- matrix(0,7,ngroups)
Eeigen5725 <- matrix(0,7,ngroups)

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
                      xlim = c(0,1), ylim = c(-0.8,0),
                      color = "grey", barcol = "lightblue",outliercol="orange")
lines(t7,eigen_true7, lty = "solid", col = "red")

GE5725 <- fda::fbplot(Geigen5725, x = t7, xlab = "Time", ylab = "", 
                      xlim = c(0,1), ylim = c(-1,1),
                      color = "grey", barcol = "lightblue",outliercol="orange")
lines(t7,eigen_true7, lty = "solid", col = "red")

GE51405 <- fda::fbplot(Geigen51405, x = t14, xlab = "Time", ylab = "", 
                       xlim = c(0,1), ylim = c(-1,1),
                       color = "grey", barcol = "lightblue",outliercol="orange")
lines(t14,eigen_true14, lty = "solid", col = "red")

GE51415 <- fda::fbplot(Geigen51415, x = t14, xlab = "Time", ylab = "", 
                       xlim = c(0,1), ylim = c(-0.8,0),
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
                      xlim = c(0,1), ylim = c(-0.8,1),
                      color = "grey", barcol = "lightblue",outliercol="orange")
lines(t7,eigen_true7, lty = "solid", col = "red")

EE5725 <- fda::fbplot(Eeigen5725, x = t7, xlab = "Time", ylab = "", 
                      xlim = c(0,1), ylim = c(-1,1),
                      color = "grey", barcol = "lightblue",outliercol="orange")
lines(t7,eigen_true7, lty = "solid", col = "red")

EE51405 <- fda::fbplot(Eeigen51405, x = t14, xlab = "Time", ylab = "", 
                       xlim = c(0,1), ylim = c(-1,1),
                       color = "grey", barcol = "lightblue",outliercol="orange")
lines(t14,eigen_true14, lty = "solid", col = "red")

EE51415 <- fda::fbplot(Eeigen51415, x = t14, xlab = "Time", ylab = "", 
                       xlim = c(0,1), ylim = c(-0.8,0.4),
                       color = "grey", barcol = "lightblue",outliercol="orange")
lines(t14,eigen_true14, lty = "solid", col = "red")

EE51425 <- fda::fbplot(Eeigen51425, x = t14, xlab = "Time", ylab = "", 
                       xlim = c(0,1), ylim = c(-1,1),
                       color = "grey", barcol = "lightblue",outliercol="orange")
lines(t14,eigen_true14, lty = "solid", col = "red")
#################################################################################

## Simultaneous confidence band
## FE
mean_FE5705 <- rowMeans(fef_5pc705)
mean_FE5715 <- rowMeans(fef_5pc715)
mean_FE5725 <- rowMeans(fef_5pc725)

mean_FE51405 <- rowMeans(fef_5pc1405)
mean_FE51415 <- rowMeans(fef_5pc1415)
mean_FE51425 <- rowMeans(fef_5pc1425)

FE_5705_diff <- rep(0,ngroups)
FE_5715_diff <- rep(0,ngroups)
FE_5725_diff <- rep(0,ngroups)

FE_51405_diff <- rep(0,ngroups)
FE_51415_diff <- rep(0,ngroups)
FE_51425_diff <- rep(0,ngroups)

for (i in 1:ngroups){
  FE_5705_diff[i] <- max(abs(fef_5pc705[,i] - mean_FE5705))
  FE_5715_diff[i] <- max(abs(fef_5pc715[,i] - mean_FE5715))
  FE_5725_diff[i] <- max(abs(fef_5pc725[,i] - mean_FE5725))
  
  FE_51405_diff[i] <- max(abs(fef_5pc1405[,i] - mean_FE51405))
  FE_51415_diff[i] <- max(abs(fef_5pc1415[,i] - mean_FE51415))
  FE_51425_diff[i] <- max(abs(fef_5pc1425[,i] - mean_FE51425))
}

FE_5705_sup <- quantile(FE_5705_diff, probs = 0.95) 
FE_5715_sup <- quantile(FE_5715_diff, probs = 0.95)
FE_5725_sup <- quantile(FE_5725_diff, probs = 0.95)

FE_51405_sup <- quantile(FE_51405_diff, probs = 0.95) 
FE_51415_sup <- quantile(FE_51415_diff, probs = 0.95)
FE_51425_sup <- quantile(FE_51425_diff, probs = 0.95)

plot(c(0,1), c(-5,300), type = "n", xlab = "Time", ylab = "")
lines(t7, mean_FE5715, lty = "solid")
lines(t7, mean_FE5715 + FE_5715_sup, lty = "dashed")
lines(t7, mean_FE5715 - FE_5715_sup, lty = "dashed")
lines(t7, f_true7, lty = "dotdash", col = "red")
mtext("Samping points per subject = 7", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

par(mfrow = c(1, 1), bty = "l")

plot(c(0,1), c(-5,300), type = "n", xlab = "Time", ylab = "")
lines(t14, mean_FE51415, lty = "solid")
lines(t14, mean_FE51415 + FE_51415_sup, lty = "dashed")
lines(t14, mean_FE51415 - FE_51415_sup, lty = "dashed")
lines(t14, f_true14, lty = "solid", col = "red")
mtext("Samping points per subject = 14", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

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

par(mfrow = c(1, 2), bty = "l")

plot(c(0,1), c(-1.5,1.5), type = "n", xlab = "Time", ylab = "")
lines(t7, mean_GE5715, lty = "solid")
lines(t7, mean_GE5715 + Geigen_5715_sup, lty = "dashed")
lines(t7, mean_GE5715 - Geigen_5715_sup, lty = "dashed")
lines(t7, eigen_true7, lty = "solid", col = "red")
mtext("Genetic Eigenfunction", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

plot(c(0,1), c(-1.5,1.5), type = "n", xlab = "Time", ylab = "")
lines(t7, mean_EE5715, lty = "solid")
lines(t7, mean_EE5715 + Eeigen_5715_sup, lty = "dashed")
lines(t7, mean_EE5715 - Eeigen_5715_sup, lty = "dashed")
lines(t7, eigen_true7, lty = "solid", col = "red")
mtext("Environmental Eigenfunction", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

########################################################################
plot(c(0,1), c(-0.8,0.4), type = "n", xlab = "Time", ylab = "")
lines(t14, mean_GE51415, lty = "solid")
lines(t14, mean_GE51415 + Geigen_51415_sup, lty = "dashed")
lines(t14, mean_GE51415 - Geigen_51415_sup, lty = "dashed")
lines(t14, eigen_true14, lty = "solid", col = "red")
mtext("Genetic Eigenfunction", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

plot(c(0,1), c(-0.8,0.4), type = "n", xlab = "Time", ylab = "")
lines(t14, mean_EE51415, lty = "solid")
lines(t14, mean_EE51415 + Eeigen_51415_sup, lty = "dashed")
lines(t14, mean_EE51415 - Eeigen_51415_sup, lty = "dashed")
lines(t14, eigen_true14, lty = "solid", col = "red")
mtext("Environmental Eigenfunction", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

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

CG_5715_upper <- mean_CG5715 + CG_5715_sup
CG_5715_lower <- mean_CG5715 - CG_5715_sup

CE_5715_upper <- mean_CE5715 + CE_5715_sup
CE_5715_lower <- mean_CE5715 - CE_5715_sup

fig <- plot_ly(showscale = FALSE)
fig <- fig %>% add_surface(z = ~C_fun_true7, x = ~t7, y = ~t7, colorscale = list(c(0, 'rgb(255, 0, 0)'), c(1, 'rgb(255, 150, 150)')))
fig <- fig %>% add_surface(z = ~CG_5715_upper, x = ~t7, y = ~t7, 
                           opacity = 0.98,colorscale = list(c(0, 'rgb(169, 169, 169)'), c(1, 'rgb(211, 211, 211)')))
fig <- fig %>% add_surface(z = ~CG_5715_lower, x = ~t7, y = ~t7, opacity = 0.98,colorscale = list(c(0, 'rgb(169, 169, 169)'), c(1, 'rgb(211, 211, 211)')))
fig <- fig %>% layout(scene = list(
  xaxis = list(title = 'Time'),
  yaxis = list(title = 'Time'),
  zaxis = list(range=c(-300,1100),title = 'Genetic Covariance Functon')
))

fig

fig <- plot_ly(showscale = FALSE)
fig <- fig %>% add_surface(z = ~C_fun_true7, x = ~t7, y = ~t7, colorscale = list(c(0, 'rgb(255, 0, 0)'), c(1, 'rgb(255, 150, 150)')))
fig <- fig %>% add_surface(z = ~CE_5715_upper, x = ~t7, y = ~t7, 
                           opacity = 0.98,colorscale = list(c(0, 'rgb(169, 169, 169)'), c(1, 'rgb(211, 211, 211)')))
fig <- fig %>% add_surface(z = ~CE_5715_lower, x = ~t7, y = ~t7, opacity = 0.98,colorscale = list(c(0, 'rgb(169, 169, 169)'), c(1, 'rgb(211, 211, 211)')))
fig <- fig %>% layout(scene = list(
  xaxis = list(title = 'Time'),
  yaxis = list(title = 'Time'),
  zaxis = list(range=c(-300,1100),title = 'Environmental Covariance Functon')
))

fig


######################################################################

CG_51415_upper <- mean_CG51415 + CG_51415_sup
CG_51415_lower <- mean_CG51415 - CG_51415_sup

CE_51415_upper <- mean_CE51415 + CE_51415_sup
CE_51415_lower <- mean_CE51415 - CE_51415_sup

fig <- plot_ly(showscale = FALSE)
fig <- fig %>% add_surface(z = ~C_fun_true14, x = ~t14, y = ~t14, colorscale = list(c(0, 'rgb(255, 0, 0)'), c(1, 'rgb(255, 150, 150)')))
fig <- fig %>% add_surface(z = ~CG_51415_upper, x = ~t14, y = ~t14, 
                           opacity = 0.98,colorscale = list(c(0, 'rgb(169, 169, 169)'), c(1, 'rgb(211, 211, 211)')))
fig <- fig %>% add_surface(z = ~CG_51415_lower, x = ~t14, y = ~t14, opacity = 0.98,colorscale = list(c(0, 'rgb(169, 169, 169)'), c(1, 'rgb(211, 211, 211)')))
fig <- fig %>% layout(scene = list(
  xaxis = list(title = 'Time'),
  yaxis = list(title = 'Time'),
  zaxis = list(range=c(-400,1200),title = 'Genetic Covariance Functon')
))

fig

fig <- plot_ly(showscale = FALSE)
fig <- fig %>% add_surface(z = ~C_fun_true14, x = ~t14, y = ~t14, colorscale = list(c(0, 'rgb(255, 0, 0)'), c(1, 'rgb(255, 150, 150)')))
fig <- fig %>% add_surface(z = ~CE_51415_upper, x = ~t14, y = ~t14, 
                           opacity = 0.98,colorscale = list(c(0, 'rgb(169, 169, 169)'), c(1, 'rgb(211, 211, 211)')))
fig <- fig %>% add_surface(z = ~CE_51415_lower, x = ~t14, y = ~t14, opacity = 0.98,colorscale = list(c(0, 'rgb(169, 169, 169)'), c(1, 'rgb(211, 211, 211)')))
fig <- fig %>% layout(scene = list(
  xaxis = list(title = 'Time'),
  yaxis = list(title = 'Time'),
  zaxis = list(range=c(-400,1200),title = 'Environmental Covariance Functon')
))

fig

