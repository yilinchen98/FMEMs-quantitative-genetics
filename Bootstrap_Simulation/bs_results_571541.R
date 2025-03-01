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
setwd("D:/KCL_2023-2027_PhD/FMEM_QuantiativeGenetics_Project/PhD_Project_Contents/R_code/Bootstrapping_simulation/5715_41")
files <- list.files(pattern = "\\.Rdata$")

for (file in files) {
  load(file)
}

load("D:/KCL_2023-2027_PhD/FMEM_QuantiativeGenetics_Project/PhD_Project_Contents/R_code/Functional_simulation/simulation_data/TC_fixed_effect.Rdata")
f_true7 <- convert_to_basisfunctions(timefine, fixed_effect, t7)
f_true14 <- convert_to_basisfunctions(timefine, fixed_effect, t14)

## Functional Boxplot

### FE
FE5715_41 <- fda::fbplot(fef_sap41_final[,1:300], x = t7, xlab = "Time", ylab = "", xlim = c(0,1),
                      ylim = c(-5,300), color = "grey", barcol = "lightblue")
lines(t7,f_true7, lty = "solid", col = "red")
mtext("ID 41", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))

### Eigenfunctions
eigen_true7 <- eigen(C_fun_true7)$vectors[,1]
eigen_true14 <- eigen(C_fun_true14)$vectors[,1]

Geigen5715_41 <- matrix(0,7,n_iter)
Eeigen5715_41 <- matrix(0,7,n_iter)

for (i in 1:n_iter){
  Geigen5715_41[,i] <- eigen(CG_fun_sap41_final[,,i])$vectors[,1]
  Eeigen5715_41[,i] <- eigen(CE_fun_sap41_final[,,i])$vectors[,1]
}

GE5715_41 <- fda::fbplot(Geigen5715_41, x = t7, xlab = "Time", ylab = "", 
                      xlim = c(0,1), ylim = c(-1,0.5),
                      color = "grey", barcol = "lightblue",outliercol="orange")
lines(t7,eigen_true7, lty = "solid", col = "red")
mtext("Genetic Eigenfunction", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))

EE5715_41 <- fda::fbplot(Eeigen5715_41, x = t7, xlab = "Time", ylab = "", 
                         xlim = c(0,1),ylim = c(-1,0.5),
                         color = "grey", barcol = "lightblue",outliercol="orange")
lines(t7,eigen_true7, lty = "solid", col = "red")
mtext("Environmental Eigenfunction", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))

#################################################################################

## SCB

## FE
FE5715_mean41 <- rowMeans(fef_sap41_300) ## mean curve
FE5715_diff41<- rep(0,n_iter)
for (i in 1:n_iter){
  FE5715_diff41[i] <- max(abs(fef_sap41_300[,i] - FE5715_mean41))
}

FE5715_bound_41 <- quantile(FE5715_diff41, probs = 0.95)

plot(c(0,1), c(-5,300), type = "n", xlab = "Time", ylab = "")
lines(t7, FE5715_mean41, lty = "solid")
lines(t7, FE5715_mean41 + FE5715_bound_41, lty = "dashed")
lines(t7, FE5715_mean41 - FE5715_bound_41, lty = "dashed")
lines(t7, f_true7, lty = "solid", col = "red")
mtext("ID 41", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))

## Eigenfuncitons
mean_GE5715_41 <- rowMeans(Geigen5715_41)
mean_EE5715_41 <- rowMeans(Eeigen5715_41)

Geigen_diff_5715_41 <- rep(0, n_iter)
Eeigen_diff_5715_41 <- rep(0, n_iter)
for (i in 1:n_iter){
  Geigen_diff_5715_41[i] <-  max(abs(Geigen5715_41[,i] - mean_GE5715_41))
  Eeigen_diff_5715_41[i] <-  max(abs(Eeigen5715_41[,i] - mean_EE5715_41))
}

Geigen_571541_sup <- quantile(Geigen_diff_5715_41, probs = 0.95) 
Eeigen_571541_sup <- quantile(Eeigen_diff_5715_41, probs = 0.95) 

plot(c(0,1), c(-1,0.5), type = "n", xlab = "Time", ylab = "")
lines(t7, mean_GE5715_41, lty = "solid")
lines(t7, mean_GE5715_41 + Geigen_571541_sup, lty = "dashed")
lines(t7, mean_GE5715_41 - Geigen_571541_sup, lty = "dashed")
lines(t7, eigen_true7, lty = "solid", col = "red")
mtext("Genetic Eigenfunction", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))

plot(c(0,1), c(-1,0.5), type = "n", xlab = "Time", ylab = "")
lines(t7, mean_EE5715_41, lty = "solid")
lines(t7, mean_EE5715_41 + Eeigen_571541_sup, lty = "dashed")
lines(t7, mean_EE5715_41 - Eeigen_571541_sup, lty = "dashed")
lines(t7, eigen_true7, lty = "solid", col = "red")
mtext("Environmental Eigenfunction", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))

## Covariance functions
mean_CG5715_41 <- apply(CG_fun_sap41_300, c(1, 2), mean)
mean_CE5715_41 <- apply(CE_fun_sap41_300, c(1, 2), mean)

CG_diff_5715_41 <- rep(0,n_iter)
CE_diff_5715_41 <- rep(0,n_iter)

for (i in 1:n_iter){
  CG_diff_5715_41[i] <- max(abs(CG_fun_sap41_300[,,i] - mean_CG5715_41))
  CE_diff_5715_41[i] <- max(abs(CE_fun_sap41_300[,,i] - mean_CE5715_41))
}
CG_571541_sup <- quantile(CG_diff_5715_41, probs = 0.95)
CE_571541_sup <- quantile(CE_diff_5715_41, probs = 0.95)

CG_571541_upper <- mean_CG5715_41 + CG_571541_sup
CG_571541_lower <- mean_CG5715_41 - CG_571541_sup

CE_571541_upper <- mean_CE5715_41 + CE_571541_sup
CE_571541_lower <- mean_CE5715_41 - CE_571541_sup

fig <- plot_ly(showscale = FALSE)
fig <- fig %>% add_surface(z = ~C_fun_true7, x = ~t7, y = ~t7, colorscale = list(c(0, 'rgb(255, 0, 0)'), c(1, 'rgb(255, 150, 150)')))
fig <- fig %>% add_surface(z = ~CG_571541_upper, x = ~t7, y = ~t7, 
                           opacity = 0.98,colorscale = list(c(0, 'rgb(169, 169, 169)'), c(1, 'rgb(211, 211, 211)')))
fig <- fig %>% add_surface(z = ~CG_571541_lower, x = ~t7, y = ~t7, opacity = 0.98,colorscale = list(c(0, 'rgb(169, 169, 169)'), c(1, 'rgb(211, 211, 211)')))
fig <- fig %>% layout(scene = list(
  xaxis = list(title = 'Time'),
  yaxis = list(title = 'Time'),
  zaxis = list(range=c(-300,1100),title = 'Genetic Covariance Functon')
))

fig

fig <- plot_ly(showscale = FALSE)
fig <- fig %>% add_surface(z = ~C_fun_true7, x = ~t7, y = ~t7, colorscale = list(c(0, 'rgb(255, 0, 0)'), c(1, 'rgb(255, 150, 150)')))
fig <- fig %>% add_surface(z = ~CE_571541_upper, x = ~t7, y = ~t7, 
                           opacity = 0.98,colorscale = list(c(0, 'rgb(169, 169, 169)'), c(1, 'rgb(211, 211, 211)')))
fig <- fig %>% add_surface(z = ~CE_571541_lower, x = ~t7, y = ~t7, opacity = 0.98,colorscale = list(c(0, 'rgb(169, 169, 169)'), c(1, 'rgb(211, 211, 211)')))
fig <- fig %>% layout(scene = list(
  xaxis = list(title = 'Time'),
  yaxis = list(title = 'Time'),
  zaxis = list(range=c(-300,1100),title = 'Environmental Covariance Functon')
))

fig
