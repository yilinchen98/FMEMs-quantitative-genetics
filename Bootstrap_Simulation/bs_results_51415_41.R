# Bootstrap Simulations

## set up true covariance
N <- 873 # total number of individualsl
nbasis <- 5 # number of basis
n_iter <- 300

### genetic and environmental covariance matrix
C_true <- matrix(c(750, 10 ,130, 80, 250,
                   10, 800, 30, 15, 40,
                   130, 30, 700, 50, 130,
                   80, 15, 50, 420, 50,
                   250, 40, 130, 50, 330), nrow = 5, byrow = T)

### model reform
t14 <- seq(0,1,length=14) # time points t_j
timefine <- seq(0,1, length =100)

basisObj <- create.bspline.basis(c(0,1), nbasis = nbasis, norder = 4)
basis14 <- eval.basis(t14,basisObj)

C_fun_true14 <- basis14 %*% C_true %*% t(basis14)


## load data

load("TC_fixed_effect.Rdata")
f_true14 <- convert_to_basisfunctions(timefine, fixed_effect, t14)

## Functional Boxplot

### FE
par(mfrow = c(1,1), bty = "l")
FE51415_41 <- fda::fbplot(fef_sap41_300, x = t14, xlab = "Time", ylab = "", xlim = c(0,1),
                      ylim = c(-5,300), color = "grey", barcol = "lightblue")
lines(t14,f_true14, lty = "solid", col = "red")
mtext("ID 41", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))

### Eigenfunctions
eigen_true14 <- eigen(C_fun_true14)$vectors[,1]

Geigen51415_41 <- matrix(0,14,n_iter)
Eeigen51415_41 <- matrix(0,14,n_iter)

for (i in 1:n_iter){
  Geigen51415_41[,i] <- eigen(CG_fun_sap41_300[,,i])$vectors[,1]
  Eeigen51415_41[,i] <- eigen(CE_fun_sap41_300[,,i])$vectors[,1]
}
par(mfrow = c(1,2), bty = "l")
GE571415_41 <- fda::fbplot(Geigen51415_41, x = t14, xlab = "Time", ylab = "", 
                      xlim = c(0,1), ylim = c(-1,0.5),
                      color = "grey", barcol = "lightblue",outliercol="orange")
lines(t14,eigen_true14, lty = "solid", col = "red")
mtext("Genetic Eigenfunction", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))

EE51415_41 <- fda::fbplot(Eeigen51415_41, x = t14, xlab = "Time", ylab = "", 
                         xlim = c(0,1),ylim = c(-1,0.5),
                         color = "grey", barcol = "lightblue",outliercol="orange")
lines(t14,eigen_true14, lty = "solid", col = "red")
mtext("Environmental Eigenfunction", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))

Outlier <- sort(unique(c(EE51415_41$outpoint, GE571415_41$outpoint)))
#################################################################################
## SCB

## FE
FE51415_mean41 <- rowMeans(fef_sap41_300) ## mean curve
FE51415_diff41<- rep(0,n_iter)
for (i in 1:n_iter){
  FE51415_diff41[i] <- max(abs(fef_sap41_300[,i] - FE51415_mean41))
}

FE51415_bound_41 <- quantile(FE51415_diff41, probs = 0.95)

plot(c(0,1), c(-5,300), type = "n", xlab = "Time", ylab = "")
lines(t14, FE51415_mean41, lty = "solid")
lines(t14, FE51415_mean41 + FE51415_bound_41, lty = "dashed")
lines(t14, FE51415_mean41 - FE51415_bound_41, lty = "dashed")
lines(t14, f_true14, lty = "solid", col = "red")
mtext("ID 41", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))

## Eigenfuncitons
mean_GE51415_41 <- rowMeans(Geigen51415_41)
mean_EE51415_41 <- rowMeans(Eeigen51415_41)

Geigen_diff_51415_41 <- rep(0, n_iter)
Eeigen_diff_51415_41 <- rep(0, n_iter)
for (i in 1:n_iter){
  Geigen_diff_51415_41[i] <-  max(abs(Geigen51415_41[,i] - mean_GE51415_41))
  Eeigen_diff_51415_41[i] <-  max(abs(Eeigen51415_41[,i] - mean_EE51415_41))
}

Geigen_5141541_sup <- quantile(Geigen_diff_51415_41, probs = 0.95) 
Eeigen_5141541_sup <- quantile(Eeigen_diff_51415_41, probs = 0.95) 

plot(c(0,1), c(-1,0.5), type = "n", xlab = "Time", ylab = "")
lines(t14, mean_GE51415_41, lty = "solid")
lines(t14, mean_GE51415_41 + Geigen_5141541_sup, lty = "dashed")
lines(t14, mean_GE51415_41 - Geigen_5141541_sup, lty = "dashed")
lines(t14, eigen_true14, lty = "solid", col = "red")
mtext("Genetic Eigenfunction", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))

plot(c(0,1), c(-1,0.5), type = "n", xlab = "Time", ylab = "")
lines(t14, mean_EE51415_41, lty = "solid")
lines(t14, mean_EE51415_41 + Eeigen_5141541_sup, lty = "dashed")
lines(t14, mean_EE51415_41 - Eeigen_5141541_sup, lty = "dashed")
lines(t14, eigen_true14, lty = "solid", col = "red")
mtext("Environmental Eigenfunction", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))

## Covariance functions
mean_CG51415_41 <- apply(CG_fun_sap41_300, c(1, 2), mean)
mean_CE51415_41 <- apply(CE_fun_sap41_300, c(1, 2), mean)

CG_diff_51415_41 <- rep(0,n_iter)
CE_diff_51415_41 <- rep(0,n_iter)

for (i in 1:n_iter){
  CG_diff_51415_41[i] <- max(abs(CG_fun_sap41_300[,,i] - mean_CG51415_41))
  CE_diff_51415_41[i] <- max(abs(CE_fun_sap41_300[,,i] - mean_CE51415_41))
}
CG_5141541_sup <- quantile(CG_diff_51415_41, probs = 0.95)
CE_5141541_sup <- quantile(CE_diff_51415_41, probs = 0.95)

CG_5141541_upper <- mean_CG51415_41 + CG_5141541_sup
CG_5141541_lower <- mean_CG51415_41 - CG_5141541_sup

CE_5141541_upper <- mean_CE51415_41 + CE_5141541_sup
CE_5141541_lower <- mean_CE51415_41 - CE_5141541_sup

fig <- plot_ly(showscale = FALSE)
fig <- fig %>% add_surface(z = ~C_fun_true14, x = ~t14, y = ~t14)
fig <- fig %>% add_surface(z = ~CG_5141541_upper, x = ~t14, y = ~t14, 
                           opacity = 0.5,colorscale = list(c(0, 1),c("grey", "lightgrey")))
fig <- fig %>% add_surface(z = ~CG_5141541_lower, x = ~t14, y = ~t14, 
                           opacity = 0.5,colorscale = list(c(0, 1),c("grey", "lightgrey")))
fig <- fig %>% layout(scene = list(
  xaxis = list(title = 'Time'),
  yaxis = list(title = 'Time'),
  zaxis = list(range=c(-300,1100),title = 'Genetic Covariance Functon')
))

fig

fig <- plot_ly(showscale = FALSE)
fig <- fig %>% add_surface(z = ~C_fun_true14, x = ~t14, y = ~t14)
fig <- fig %>% add_surface(z = ~CE_5141541_upper, x = ~t14, y = ~t14, 
                           opacity = 0.5,colorscale = list(c(0, 1),c("grey", "lightgrey")))
fig <- fig %>% add_surface(z = ~CE_5141541_lower, x = ~t14, y = ~t14, 
                           opacity = 0.5,colorscale = list(c(0, 1),c("grey", "lightgrey")))
fig <- fig %>% layout(scene = list(
  xaxis = list(title = 'Time'),
  yaxis = list(title = 'Time'),
  zaxis = list(range=c(-300,1100),title = 'Environmental Covariance Functon')
))

fig