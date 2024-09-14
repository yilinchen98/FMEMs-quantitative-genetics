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
t7 <- seq(0,1,length=7) # time points t_j
t14 <- seq(0,1,length=14) # time points t_j
timefine <- seq(0,1, length =100)

basisObj <- create.bspline.basis(c(0,1), nbasis = nbasis, norder = 4)
basis7 <- eval.basis(t7,basisObj)
basis14 <- eval.basis(t14,basisObj)

C_fun_true7 <- basis7 %*% C_true %*% t(basis7)
C_fun_true14 <- basis14 %*% C_true %*% t(basis14)


## load data
setwd("D:/KCL_2023-2027_PhD/Genetics_FMEM_Project/FMEMs-quantitative-genetics/R_code/Bootstrapping_simulation/bs_PCs/5715_43")

files <- list.files(pattern = "\\.Rdata$")

for (file in files) {
  load(file)
}

load("D:/KCL_2023-2027_PhD/Genetics_FMEM_Project/FMEMs-quantitative-genetics/R_code/Functional_simulation/Simulation_PCs/PCs_FitRawData/simulated_data/TC_fixed_effect.Rdata")
f_true7 <- convert_to_basisfunctions(timefine, fixed_effect, t7)
f_true14 <- convert_to_basisfunctions(timefine, fixed_effect, t14)

## Functional Boxplot

### FE
FE5715_43 <- fda::fbplot(fef_sap43_300, x = t7, xlab = "Time", ylab = "", xlim = c(0,1),
                      ylim = c(-5,300), color = "grey", barcol = "lightblue")
lines(t7,f_true7, lty = "solid", col = "red")
mtext("ID 43", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

### Eigenfunctions
eigen_true7 <- eigen(C_fun_true7)$vectors[,1]
eigen_true14 <- eigen(C_fun_true14)$vectors[,1]

Geigen5715_43 <- matrix(0,7,n_iter)
Eeigen5715_43 <- matrix(0,7,n_iter)

for (i in 1:n_iter){
  Geigen5715_43[,i] <- eigen(CG_fun_sap43_300[,,i])$vectors[,1]
  Eeigen5715_43[,i] <- eigen(CE_fun_sap43_300[,,i])$vectors[,1]
}

GE5715_43 <- fda::fbplot(Geigen5715_43, x = t7, xlab = "Time", ylab = "", 
                      xlim = c(0,1), ylim = c(-1,0.5),
                      color = "grey", barcol = "lightblue",outliercol="orange")
lines(t7,eigen_true7, lty = "solid", col = "red")
mtext("Genetic Eigenfunction", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

EE5715_43 <- fda::fbplot(Eeigen5715_43, x = t7, xlab = "Time", ylab = "", 
                         xlim = c(0,1),ylim = c(-1,0.5),
                         color = "grey", barcol = "lightblue",outliercol="orange")
lines(t7,eigen_true7, lty = "solid", col = "red")
mtext("Environmental Eigenfunction", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

#################################################################################

## SCB

## FE
FE5715_mean43 <- rowMeans(fef_sap43_300) ## mean curve
FE5715_diff43<- rep(0,n_iter)
for (i in 1:n_iter){
  FE5715_diff43[i] <- max(abs(fef_sap43_300[,i] - FE5715_mean43))
}

FE5715_bound_43 <- quantile(FE5715_diff43, probs = 0.95)

plot(c(0,1), c(-5,300), type = "n", xlab = "Time", ylab = "")
lines(t7, FE5715_mean43, lty = "solid")
lines(t7, FE5715_mean43 + FE5715_bound_43, lty = "dashed")
lines(t7, FE5715_mean43 - FE5715_bound_43, lty = "dashed")
lines(t7, f_true7, lty = "solid", col = "red")
mtext("ID 43", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

## Eigenfuncitons
mean_GE5715_43 <- rowMeans(Geigen5715_43)
mean_EE5715_43 <- rowMeans(Eeigen5715_43)

Geigen_diff_5715_43 <- rep(0, n_iter)
Eeigen_diff_5715_43 <- rep(0, n_iter)
for (i in 1:n_iter){
  Geigen_diff_5715_43[i] <-  max(abs(Geigen5715_43[,i] - mean_GE5715_43))
  Eeigen_diff_5715_43[i] <-  max(abs(Eeigen5715_43[,i] - mean_EE5715_43))
}

Geigen_571543_sup <- quantile(Geigen_diff_5715_43, probs = 0.95) 
Eeigen_571543_sup <- quantile(Eeigen_diff_5715_43, probs = 0.95) 

plot(c(0,1), c(-1,0.5), type = "n", xlab = "Time", ylab = "")
lines(t7, mean_GE5715_43, lty = "solid")
lines(t7, mean_GE5715_43 + Geigen_571543_sup, lty = "dashed")
lines(t7, mean_GE5715_43 - Geigen_571543_sup, lty = "dashed")
lines(t7, eigen_true7, lty = "solid", col = "red")
mtext("Genetic Eigenfunction", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

plot(c(0,1), c(-1,0.5), type = "n", xlab = "Time", ylab = "")
lines(t7, mean_EE5715_43, lty = "solid")
lines(t7, mean_EE5715_43 + Eeigen_571543_sup, lty = "dashed")
lines(t7, mean_EE5715_43 - Eeigen_571543_sup, lty = "dashed")
lines(t7, eigen_true7, lty = "solid", col = "red")
mtext("Environmental Eigenfunction", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

## Covariance functions
mean_CG5715_43 <- apply(CG_fun_sap43_300, c(1, 2), mean)
mean_CE5715_43 <- apply(CE_fun_sap43_300, c(1, 2), mean)

CG_diff_5715_43 <- rep(0,n_iter)
CE_diff_5715_43 <- rep(0,n_iter)

for (i in 1:n_iter){
  CG_diff_5715_43[i] <- max(abs(CG_fun_sap43_300[,,i] - mean_CG5715_43))
  CE_diff_5715_43[i] <- max(abs(CE_fun_sap43_300[,,i] - mean_CE5715_43))
}
CG_571543_sup <- quantile(CG_diff_5715_43, probs = 0.95)
CE_571543_sup <- quantile(CE_diff_5715_43, probs = 0.95)

CG_571543_upper <- mean_CG5715_43 + CG_571543_sup
CG_571543_lower <- mean_CG5715_43 - CG_571543_sup

CE_571543_upper <- mean_CE5715_43 + CE_571543_sup
CE_571543_lower <- mean_CE5715_43 - CE_571543_sup

fig <- plot_ly(showscale = FALSE)
fig <- fig %>% add_surface(z = ~C_fun_true7, x = ~t7, y = ~t7, colorscale = list(c(0, 'rgb(255, 0, 0)'), c(1, 'rgb(255, 150, 150)')))
fig <- fig %>% add_surface(z = ~CG_571543_upper, x = ~t7, y = ~t7, 
                           opacity = 0.98,colorscale = list(c(0, 'rgb(169, 169, 169)'), c(1, 'rgb(211, 211, 211)')))
fig <- fig %>% add_surface(z = ~CG_571543_lower, x = ~t7, y = ~t7, opacity = 0.98,colorscale = list(c(0, 'rgb(169, 169, 169)'), c(1, 'rgb(211, 211, 211)')))
fig <- fig %>% layout(scene = list(
  xaxis = list(title = 'Time'),
  yaxis = list(title = 'Time'),
  zaxis = list(range=c(-300,1200),title = 'Genetic Covariance Functon')
))

fig

fig <- plot_ly(showscale = FALSE)
fig <- fig %>% add_surface(z = ~C_fun_true7, x = ~t7, y = ~t7, colorscale = list(c(0, 'rgb(255, 0, 0)'), c(1, 'rgb(255, 150, 150)')))
fig <- fig %>% add_surface(z = ~CE_571543_upper, x = ~t7, y = ~t7, 
                           opacity = 0.98,colorscale = list(c(0, 'rgb(169, 169, 169)'), c(1, 'rgb(211, 211, 211)')))
fig <- fig %>% add_surface(z = ~CE_571543_lower, x = ~t7, y = ~t7, opacity = 0.98,colorscale = list(c(0, 'rgb(169, 169, 169)'), c(1, 'rgb(211, 211, 211)')))
fig <- fig %>% layout(scene = list(
  xaxis = list(title = 'Time'),
  yaxis = list(title = 'Time'),
  zaxis = list(range=c(-300,1200),title = 'Environmental Covariance Functon')
))

fig