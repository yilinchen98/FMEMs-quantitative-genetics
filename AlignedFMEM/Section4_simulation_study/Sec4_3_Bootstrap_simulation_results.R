# Bootstrap Simulation Result (Best estimation of the additive genetic covariance function)
source("PrepareR.R")
## model set up 
N <- 873 # total number of individuals
nbasis <- 5 # number of basis
n_iter <- 300

### true genetic and environmental covariance matrix
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
eigen_true14 <- eigen(C_fun_true14)$vectors[,1]
### true fixed effect
fixed_effect <- readRDS("sim_results/fixed_effect.rds")
f_true14 <- convert_to_basisfunctions(timefine, fixed_effect, t14)

### Functional estimation
fef_list <- readRDS("sim_results/fef_list.rds")
CG_fun_mat <- readRDS("sim_results/CG_fun_5pc1415.rds")
CE_fun_mat <- readRDS("sim_results/CE_fun_5pc1415.rds")

Geigen_mat <- matrix(0,14,50)
Eeigen_mat <- matrix(0,14,50)
for(i in 1:50){
  Geigen_est <- eigen(CG_fun_mat[,,i])$vectors[,1]
  if(trapz(t14,eigen_true14 * Geigen_est) < 0) Geigen_est <- -Geigen_est
  Geigen_mat[,i] <- Geigen_est
  
  Eeigen_est <- eigen(CE_fun_mat[,,i])$vectors[,1]
  if(trapz(t14, eigen_true14 * Eeigen_est) < 0) Eeigen_est <- -Eeigen_est
  Eeigen_mat[,i] <- Eeigen_est
}
## load bootstrap estimation results
bs_result_41 <- readRDS("sim_results/bs_result_41.rds")
SCB_result_41 <- readRDS("sim_results/SCB_result_41.rds")
fef_sap41_300 <- bs_result_41$fef_fin
CG_fun_sap41_300 <- bs_result_41$CG_fun_fin
CE_fun_sap41_300 <- bs_result_41$CE_fun_fin
#################################################################################
## SCB

## FE
par(mfrow = c(1,1), bty = "l")
plot(c(0,1), c(-5,300), type = "n", xlab = "Time", ylab = "")
lines(t14, fef_list[[41]], lty = "solid")
lines(t14, fef_list[[41]] + SCB_result_41$SCB_values$FE_sup, lty = "dashed")
lines(t14, fef_list[[41]] - SCB_result_41$SCB_values$FE_sup, lty = "dashed")
lines(t14, f_true14, lty = "solid", col = "red")
mtext("Fixed Effect", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))

## Eigenfuncitons
par(mfrow = c(1,2), bty = "l")
plot(c(0,1), c(-1,0.5), type = "n", xlab = "Time", ylab = "")
lines(t14, Geigen_mat[,41], lty = "solid")
lines(t14, Geigen_mat[,41] + SCB_result_41$SCB_values$Geigen_sup, lty = "dashed")
lines(t14, Geigen_mat[,41] - SCB_result_41$SCB_values$Geigen_sup, lty = "dashed")
lines(t14, eigen_true14, lty = "solid", col = "red")
mtext("Genetic Eigenfunction", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))

plot(c(0,1), c(-1,0.5), type = "n", xlab = "Time", ylab = "")
lines(t14, Eeigen_mat[,41], lty = "solid")
lines(t14, Eeigen_mat[,41] + SCB_result_41$SCB_values$Eeigen_sup, lty = "dashed")
lines(t14, Eeigen_mat[,41] - SCB_result_41$SCB_values$Eeigen_sup, lty = "dashed")
lines(t14, eigen_true14, lty = "solid", col = "red")
mtext("Environmental Eigenfunction", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))

## Covariance functions
plt_CG_SCB <- plot_ly(showscale = FALSE, scene='scene1') 
plt_CG_SCB <- plt_CG_SCB %>% add_surface(z = ~C_fun_true14, x = ~t14, y = ~t14)
plt_CG_SCB <- plt_CG_SCB %>% add_surface(z = ~CG_fun_mat[,,41] + SCB_result_41$SCB_values$CG_sup, x = ~t14, y = ~t14, 
                                         opacity = 0.5,colorscale =list(c(0, 1), c("grey", "lightgrey")))
plt_CG_SCB <- plt_CG_SCB %>% add_surface(z = ~CG_fun_mat[,,41] - SCB_result_41$SCB_values$CG_sup, x = ~t14, y = ~t14, 
                                         opacity = 0.5,colorscale =list(c(0, 1), c("grey", "lightgrey")))

plt_CE_SCB <- plot_ly(showscale = FALSE, scene='scene2') 
plt_CE_SCB <- plt_CE_SCB %>% add_surface(z = ~C_fun_true14, x = ~t14, y = ~t14)
plt_CE_SCB <- plt_CE_SCB %>% add_surface(z = ~CE_fun_mat[,,41] + SCB_result_41$SCB_values$CE_sup, x = ~t14, y = ~t14, 
                                         opacity = 0.5,colorscale =list(c(0, 1), c("grey", "lightgrey")))
plt_CE_SCB <- plt_CE_SCB %>% add_surface(z = ~CE_fun_mat[,,41] - SCB_result_41$SCB_values$CE_sup, x = ~t14, y = ~t14, 
                                         opacity = 0.5,colorscale =list(c(0, 1), c("grey", "lightgrey")))

plt_SCB <- subplot(plt_CG_SCB, plt_CE_SCB)
plt_SCB <- plt_SCB %>% layout(scene = list(domain=list(x=c(0,0.5),y=c(0,1)),
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
                                text = "Bootstrap Simultaneous Confidence Bands",
                                x = 0.5,
                                xanchor = "center"
                              ))
plt_SCB
################################################################################

