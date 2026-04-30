# Functional simulation result
source("PrepareR.R")
## Simulation setting: 873 individuals, 14 measurement points for each individual, res var = 15 

N <- 873 # total number of individuals
n <- 14
nbasis <- 5 # number of basis
ngroups <- 50
#ngroups <- 20

### genetic and environmental covariance matrix
C_true <- matrix(c(750, 10 ,130, 80, 250,
                   10, 800, 30, 15, 40,
                   130, 30, 700, 50, 130,
                   80, 15, 50, 420, 50,
                   250, 40, 130, 50, 330), nrow = 5, byrow = T)

t14 <- seq(0,1,length=14) # time points t_j
timefine <- seq(0,1, length =100)

basisObj <- create.bspline.basis(c(0,1), nbasis = nbasis, norder = 4)
basis14 <- eval.basis(t14,basisObj)
C_fun_true14 <- basis14 %*% C_true %*% t(basis14)

### fixed-effect
fixed_effect <- readRDS("fixed_effect.rds")

## load estimation results (for 50 groups)
CE_list <- readRDS("sim_results/CE_list.rds") # estimated env cov matrix
CG_list <- readRDS("sim_results/CG_list.rds") # estimated gen cov matrix
fef_list <- readRDS("sim_results/fef_list.rds") # estimated fixed effect
fpcs_list <- readRDS("sim_results/fpcs_list.rds") # FPC basis
sfit <- readRDS("sim_results/sfit.rds")

f_true14 <- convert_to_basisfunctions(timefine, fixed_effect, t14)

## Plot of the true fixed effect
par(mfrow = c(1,1), bty = "l")
plot(c(0,1), c(0, 280), type = "n", 
     xlab = "Time", 
     ylab = "Fixed Effect",
     xlim = c(0, 1), ylim = c(0,280), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
lines(t14, f_true14, type = "l", col = "black", lwd = 1.5)
mtext("True Fixed Effect", side = 3, adj = 0, line = 1, font = 2)
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0, lwd = 2) 
axis(side = 2, at = seq(0, 280, by = 40), pos = 0, lwd = 2) 

## Plot of the true covariance function
plt_CG <- plot_ly(z = ~C_fun_true14, x = ~t14, y = ~t14, showscale = FALSE) %>%
  add_surface() %>%
  layout(
    scene = list(
      xaxis = list(title = "Time"),
      yaxis = list(title = "Time"),
      zaxis = list(title = "Covariance Function")))

plt_CG
################################################################################
## Compute the genetic and environmental covariance functions
CG_fun_5pc1415 <- array(data = 0, dim = c(n,n,ngroups))
CE_fun_5pc1415 <- array(data = 0, dim = c(n,n,ngroups))

for(k in 1:ngroups){
  CG_fun_5pc1415[,,k] <- fpcs_list[[k]][,1:nbasis] %*% CG_list[[k]] %*% t(fpcs_list[[k]][,1:nbasis])
  CE_fun_5pc1415[,,k] <- fpcs_list[[k]][,1:nbasis] %*% CE_list[[k]] %*% t(fpcs_list[[k]][,1:nbasis])
}
################################################################################
## Functional boxplot 
### Fixed effect
fef_5pc1415 <- do.call(cbind,fef_list)
FE51415 <- fda::fbplot(fef_5pc1415, x = t14, xlab = "Time", ylab = "Fixed Effect", xlim = c(0,1),
                       ylim = c(-5,300),color = "grey",
                       outliercol = "#E87722",
                       barcol = "#005BAA")
lines(t14,f_true14, lty = "dotted", col = "red", lwd = 2)
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
mtext("Fixed Effect", side = 3, adj = 0, line = 1, font = 2)

### Genetic and environmental eigenfunction 
eigen_true14 <- eigen(C_fun_true14)$vectors[,1]
Geigen51415 <- matrix(0,14,ngroups)
Eeigen51415 <- matrix(0,14,ngroups)

for(k in 1:ngroups){
  Geigen51415[,k] <- eigen(CG_fun_5pc1415[,,k])$vectors[,1]
  Eeigen51415[,k] <- eigen(CE_fun_5pc1415[,,k])$vectors[,1]
}

par(mfrow = c(1,2))
GE51415 <- fda::fbplot(Geigen51415, x = t14, xlab = "Time", ylab = "Genetic Eigenfunction", 
                       xlim = c(0,1), ylim = c(-0.8,0.4),
                       color = "grey", barcol = "#005BAA",outliercol="#E87722")
lines(t14,eigen_true14, lty = "dotted", col = "red", lwd = 2)
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
mtext("Genetic Eigenfunction", side = 3, adj = 0, line = 1, font = 2)

EE51415 <- fda::fbplot(Eeigen51415, x = t14, xlab = "Time", ylab = "Environmental Eigenfunction", 
                       xlim = c(0,1), ylim = c(-0.8,0.4),
                       color = "grey", barcol = "#005BAA",outliercol="#E87722",
                       cex.axis = 0.9)
lines(t14,eigen_true14, lty = "dotted", col = "red", lwd = 2)
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
mtext("Environmental Eigenfunction", side = 3, adj = 0, line = 1, font = 2)

################################################################################
## Check singular results and outliers
OutID <- unique(sort(c(GE51415$outpoint, EE51415$outpoint)))
sfitID <- which(sfit)

## Functional distribution band
### Genetic and environmental covariance function
mean_CG51415 <- apply(CG_fun_5pc1415, c(1, 2), mean)
mean_CE51415 <- apply(CE_fun_5pc1415, c(1, 2), mean)
CG_diff_51415 <- rep(0,ngroups)
CE_diff_51415 <- rep(0,ngroups)

for(k in 1:ngroups){
  CG_diff_51415[k] <- max(abs(CG_fun_5pc1415[,,k] - mean_CG51415))
  CE_diff_51415[k] <- max(abs(CE_fun_5pc1415[,,k] - mean_CE51415))
}
CG_51415_sup <- quantile(CG_diff_51415, probs = 0.95)
CE_51415_sup <- quantile(CE_diff_51415, probs = 0.95)

CG_51415_upper <- mean_CG51415 + CG_51415_sup
CG_51415_lower <- mean_CG51415 - CG_51415_sup

CE_51415_upper <- mean_CE51415 + CE_51415_sup
CE_51415_lower <- mean_CE51415 - CE_51415_sup

plt_CG <- plot_ly(showscale = FALSE, scene='scene1') 
plt_CG <- plt_CG %>% add_surface(z = ~C_fun_true14, x = ~t14, y = ~t14)
plt_CG <- plt_CG %>% add_surface(z = ~CG_51415_upper, x = ~t14, y = ~t14, 
                                 opacity = 0.5,colorscale =list(c(0, 1), c("grey", "lightgrey")))
plt_CG <- plt_CG %>% add_surface(z = ~CG_51415_lower, x = ~t14, y = ~t14, 
                                 opacity = 0.5,colorscale =list(c(0, 1), c("grey", "lightgrey")))

plt_CE <- plot_ly(showscale = FALSE, scene='scene2') 
plt_CE <- plt_CE %>% add_surface(z = ~C_fun_true14, x = ~t14, y = ~t14)
plt_CE <- plt_CE %>% add_surface(z = ~CE_51415_upper, x = ~t14, y = ~t14, 
                                 opacity = 0.5,colorscale =list(c(0, 1), c("grey", "lightgrey")))
plt_CE <- plt_CE %>% add_surface(z = ~CE_51415_lower, x = ~t14, y = ~t14, 
                                 opacity = 0.5,colorscale =list(c(0, 1), c("grey", "lightgrey")))

plt <- subplot(plt_CG, plt_CE)
plt <- plt %>% layout(scene = list(domain=list(x=c(0,0.5),y=c(0,1)),
                                   xaxis=list(title = "Time"),
                                   yaxis =list(title = "Time") , 
                                   zaxis=list(range=c(-400,1200),title = "Genetic Covariance Function"),
                                   aspectmode='cube'),
                      scene2 = list(domain=list(x=c(0.5,1),y=c(0,1)),
                                    xaxis=list(title = "Time"),
                                    yaxis =list(title = "Time"),
                                    zaxis=list(range=c(-400,1200),title = "Envrionmental Covariance Function"),
                                    aspectmode='cube'),
                      title = list(
                        text = "Functional Distribution Bands",
                        x = 0.5,
                        xanchor = "center"
                      ))
plt
################################################################################