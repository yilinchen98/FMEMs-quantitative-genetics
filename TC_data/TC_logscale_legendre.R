# We repeat Irwin and Carter (2013) estimation for TC covariance functions.
# Here we fit the genetic FEME for raw data with orthogonal Legendre polynomials as basis functions. 

## Load Data
setwd("D:/KCL_2023-2027_PhD/FMEM_QuantitativeGenetics_Project/PhD_Project_Contents/R_code/TC_data")
TRFUN25PUP4 = read.delim("TRFUN25PUP4.DAT",header = FALSE)
names(TRFUN25PUP4)<-c("id","sire","dam","trait","x")
df <- data.frame(TRFUN25PUP4)

## Calculate genetic relationship matrix
id <- df$id
FirstUniqueIdPos <- which(duplicated(id) == FALSE)
pos = id[FirstUniqueIdPos] # extract ids for all subjects
sire_id = df$sire[FirstUniqueIdPos] # extract ids for sire
dam_id = df$dam[FirstUniqueIdPos] # extract ids for dam

pede <- editPed(sire = sire_id, dam = dam_id, label = pos)
ped<- with(pede, pedigree(label=label, sire=sire, dam=dam))
A <- getA(ped)[163:1035,163:1035]

## Transform to log scale
N = length(FirstUniqueIdPos) # N = 873 subjects
n = length(id) # n = 6860 observations
age_list <- split(df$x,id)
trait_list <- split(df$trait,id)
df$logtrait <- log10(df$trait)
log_trait_list <- split(df$logtrait,id)

## Scale time
### x = (x -min(x))/(max(x) - min(x)) scale time between [0,1]
age_list_new <- list()
for (i in 1:N){
  age_list_new[[i]] = (age_list[[i]]-min(age_list[[i]]))/(max(age_list[[i]])-min(age_list[[i]]))
}

## Rescale time interval to [-1,1]
df$x_new <-  df$x*(1/13)-1
x_new <- split(df$x_new, id)

par(mfrow = c(1,1))
par(mar = c(5, 6, 4, 2) + 0.1, font.main = 1)
plot(c(-1,1), c(0,3), type = "n",
     xlab = "Time", 
     ylab =  expression(Log(Mass,~10^-5~g)),
     xlim = c(-1, 1), ylim = c(0,3), 
     xaxs = "i", yaxs = "i",
     axes = FALSE)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(x_new[[i]], log_trait_list[[i]],type = "l", col = i)
}
axis(side = 1, at = seq(-1, 1, by = 0.2), pos = 0) 
axis(side = 2, at = seq(0, 3, by = 0.5), pos = 0) 
mtext("Log Mass Curves", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, v=0)
#############################################################################

# Redo the analysis of Irwin and Carter 2013

## orthonormal Legendre polynomial

### warning: some predictor variables are on different scales: consider rescaling (lme4)
library(sommer)
legendrePoly3_1 <- leg(x_new[[1]], n = 3)

## Example of Legendre polynomial
par(mfrow = c(1,2))
plot(c(-1,1), c(-2,2), type = "n",
     xlab = "Time", 
     ylab = "",
     xlim = c(-1, 1), ylim = c(-2, 2), 
     xaxs = "i", yaxs = "i",
     axes = FALSE)
axis(side = 1, at = seq(-1, 1, by = 0.1), pos = 0) 
axis(side = 2, at = seq(-2, 2, by = 0.5), pos = 0) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:4){
  lines(x_new[[1]], legendrePoly3_1[,i], col = i)
}
mtext("Orthonormal Legendre Poly of order 3 (Subject 1)", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, v=0)

legendrePoly3_2 <- leg(x_new[[2]], n = 3)
plot(c(-1,1), c(-2,2), type = "n",
     xlab = "Time", 
     ylab = "",
     xlim = c(-1, 1), ylim = c(-2, 2), 
     xaxs = "i", yaxs = "i",
     axes = FALSE)
axis(side = 1, at = seq(-1, 1, by = 0.1), pos = 0) 
axis(side = 2, at = seq(-2, 2, by = 0.5), pos = 0) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:4){
  lines(x_new[[2]], legendrePoly3_2[,i], col = i)
}
mtext("Orthonomal Legendre Poly of order 3 (Subject 2)", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, v=0)

## Orthogonal Legendre Basis

library(orthopolynom)
legendre_polys3 <- legendre.polynomials(n = 3, normalized = FALSE) # 3rd order Legendre polynomial
values1 <- polynomial.values(polynomials = legendre_polys3, x = x_new[[1]])
values2 <- polynomial.values(polynomials = legendre_polys3, x = x_new[[2]])

## Example of legendre polynomial
par(mfrow = c(1,2))
plot(c(-1,1), c(-2,2), type = "n",
     xlab = "Time", 
     ylab = "",
     xlim = c(-1, 1), ylim = c(-2, 2), 
     xaxs = "i", yaxs = "i",
     axes = FALSE)
axis(side = 1, at = seq(-1, 1, by = 0.1), pos = 0) 
axis(side = 2, at = seq(-2, 2, by = 0.5), pos = 0) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:4){
  lines(x_new[[1]], values1[[i]], col = i)
}
mtext("Orthogonal Legendre Poly of order 3 (Subject 1)", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, v=0)

plot(c(-1,1), c(-2,2), type = "n",
     xlab = "Time", 
     ylab = "",
     xlim = c(-1, 1), ylim = c(-2, 2), 
     xaxs = "i", yaxs = "i",
     axes = FALSE)
axis(side = 1, at = seq(-1, 1, by = 0.1), pos = 0) 
axis(side = 2, at = seq(-2, 2, by = 0.5), pos = 0) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:4){
  lines(x_new[[2]], values2[[i]], col = i)
}
mtext("Orthogonal Legendre Poly of order 3 (Subject 2)", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, v=0)

######################################################################################

## Construct basis (function) matrix for 873 subjects
legendre_polys3 <- legendre.polynomials(n = 3, normalized = FALSE)
legBasis3 <- list() # 3rd order polynomial (non-normalised)
for (i in 1:N){
  values <- polynomial.values(polynomials = legendre_polys3, x = x_new[[i]])
  legBasis3[[i]] <- do.call(cbind, values)
}
basis <- do.call(rbind, legBasis3)

## Fit both genetic and environmental random effects using 3rd order Legendre polynomials
df_leg <-  data.frame(id = id, 
                      trait = df$logtrait, 
                      phi1 = basis[,1],
                      phi2 = basis[,2],
                      phi3 = basis[,3],
                      phi4 = basis[,4])

fformLegendre <- trait ~ -1 + df_leg$phi1 + df_leg$phi2 + df_leg$phi3 + df_leg$phi4 +
  (-1 + df_leg$phi1 + df_leg$phi2 + df_leg$phi3 + df_leg$phi4 | df_leg$id) + 
  (-1 + df_leg$phi1 + df_leg$phi2 + df_leg$phi3 + df_leg$phi4 | df_leg$id)

system.time(
  ffLeg <- fit_genetic_fmm(fformLegendre, df_leg, A, 4)
) # user   system  elapsed

isSingular(ffLeg) # FALSE

summary(ffLeg)

## sparse FPCA 
logFPCA <- FPCA(Ly = log_trait_list, Lt = x_new, optns = list(nRegGrid = 100)) # FPCA through conditional expectation
timegrid <- logFPCA$workGrid
sample_mean <- logFPCA$mu
fpca_funs <- logFPCA$phi

par(mfrow = c(2, 1))
par(font.main = 1)
plot(c(-1,1), c(-0.4, 1.4), type = "n", 
     xlab = "Time", 
     ylab = "",
     xlim = c(-1, 1), ylim = c(-0.4,1.4), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(-1, 1, by = 0.2), pos = 0) 
axis(side = 2, at = seq(-0.4, 1.2, by = 0.2), pos = 0) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(timegrid, fpca_funs[,1], type = "l", col = i)
}
mtext("Principal Component 1 of Log Curves", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, v=0)

par(font.main = 1)
plot(c(-1,1), c(-0.4, 1.4), type = "n", 
     xlab = "Time", 
     ylab = "",
     xlim = c(-1, 1), ylim = c(-0.4,1.4), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(-1, 1, by = 0.2), pos = 0) 
axis(side = 2, at = seq(-0.4, 1.2, by = 0.2), pos = 0) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(timegrid, fpca_funs[,2], type = "l", col = i)
}
mtext("Principal Component 2 of Log Curves", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, v=0)
#################################################################################

## Sample covariance function
sample_cov_process <- GetCovSurface(Ly = log_trait_list, Lt = x_new, optns = list(nRegGrid = 100)) # smooth kernel: Gauss
sam_orin_process <- GetCovSurface(Ly = log_trait_list, Lt = age_list_new, optns = list(nRegGrid = 100))
timefine <- sam_orin_process$workGrid
sample_cov <- sample_cov_process$cov
sample_orin <- sam_orin_process$cov

fig1 <- plot_ly(showscale=FALSE, scene='scene1') 
fig1 <- fig1 %>% add_surface(z = ~sample_orin, x = timefine, y = timefine)

fig2 <- plot_ly(showscale = FALSE, scene='scene2') 
fig2 <- fig2 %>% add_surface(z = ~sample_cov, x= timegrid, y = timegrid)

fig_RR <- subplot(fig1, fig2) 
fig_RR <- fig_RR %>% layout(scene = list(title = "",
                                         domain=list(x=c(0,0.45),y=c(0.25,1)),
                                         xaxis=list(title = "Time"),
                                         yaxis =list(title = "Time") , 
                                         zaxis=list(range=c(-0.01,0.1),title = "Sample Covariance (0 < t < 1 )"),
                                         aspectmode='cube'),
                            scene2 = list(domain=list(x=c(0.50,0.95),y=c(0.25,1)),
                                          xaxis=list(title = "Time"),
                                          yaxis =list(title = "Time"),
                                          zaxis=list(range=c(-0.01,0.1),title = "Sample Covariance (-1 < t < 1)"),
                                          aspectmode='cube'))

fig_RR
#########################################################################################

## Model results

### Fixed effects
valuesfine <- polynomial.values(polynomials = legendre_polys3, x = timegrid)
basisfine <- do.call(cbind, valuesfine)

beta_leg <- fixef(ffLeg)
FE_leg <- basisfine %*% beta_leg

par(mfrow = c(2,1))
par(mar = c(5, 6, 4, 2) + 0.1, font.main = 1)
plot(c(-1,1), c(0,3), type = "n",
     xlab = "Time", 
     ylab =  expression(Log(Mass,~10^-5~g)),
     xlim = c(-1, 1), ylim = c(0,3), 
     xaxs = "i", yaxs = "i",
     axes = FALSE)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(x_new[[i]], log_trait_list[[i]],type = "l", col = i)
}
lines(timegrid, sample_mean, type = "l", col = "red", lwd = 3.0, lty = "solid")
axis(side = 1, at = seq(-1, 1, by = 0.2), pos = 0) 
axis(side = 2, at = seq(0, 3, by = 0.5), pos = 0) 
mtext("Log Mass Curves", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, v=0)

par(font.main = 1)
plot(c(-1,1), c(0, 3), type = "n", 
     xlab = "Time", 
     ylab = "",
     xlim = c(-1, 1), ylim = c(0,3), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(-1, 1, by = 0.2), pos = 0) 
axis(side = 2, at = seq(0, 3, by = 0.5), pos = 0) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
lines(timegrid, FE_leg, type = "l", lty = "dashed", cpl = "black")
lines(timegrid, sample_mean, type = "l", col = "red")
mtext("Fixed Effect", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, v=0)
legend("bottomright", legend= c("estimated FE", "sample mean"), lty = c("dashed", "solid"), bty = "n", col = c("black", "red"))
#################################################################################################

### Random effects
vcLeg <- VarCorr(ffLeg)
CGLeg <- vcLeg[["df_leg.id"]] # genetic covariance
CELeg <- vcLeg[["df_leg.id.1"]] # environmental covariance

CG_funLeg <- basisfine%*% CGLeg %*% t(basisfine) # estimated gen cov function
CE_funLeg <- basisfine %*% CELeg %*% t(basisfine) # estimated env cov function

fig3 <- plot_ly(showscale=FALSE, scene='scene1') 
fig3 <- fig3 %>% add_surface(z = ~CG_funLeg, x = timegrid, y = timegrid)

fig4 <- plot_ly(showscale = FALSE, scene='scene2') 
fig4 <- fig4 %>% add_surface(z = ~CE_funLeg, x= timegrid, y = timegrid)

fig_RR2 <- subplot(fig3, fig4) 
fig_RR2 <- fig_RR2 %>% layout(scene = list(title = "",
                                           domain=list(x=c(0,0.45),y=c(0.25,1)),
                                           xaxis=list(title = "Time"),
                                           yaxis =list(title = "Time") , 
                                           zaxis=list(range = c(-0.05, 0.3),title = "Genetic Covariance (Log(Mass, 10<sup>-5</sup> g/t))"),
                                           aspectmode='cube'),
                              scene2 = list(domain=list(x=c(0.50,0.95),y=c(0.25,1)),
                                            xaxis=list(title = "Time"),
                                            yaxis =list(title = "Time"),
                                            zaxis=list(range = c(-0.05, 0.3), title = "Environmental Covariance (Log(Mass, 10<sup>-5</sup> g/t))"),
                                            aspectmode='cube'))

fig_RR2

P_funLeg <- CG_funLeg + CE_funLeg
fig5 <- plot_ly(showscale = FALSE, scene='scene2') 
fig5 <- fig5 %>% add_surface(z = ~P_funLeg, x= timefine, y = timefine)

fig_RR3 <- subplot(fig5, fig3) 
fig_RR3 <- fig_RR3 %>% layout(title = "",
                              scene2 = list(domain=list(x=c(0,0.45),y=c(0.25,1)),
                                            xaxis=list(title = "Time"),
                                            yaxis =list(title = "Time") , 
                                            zaxis=list(range = c(-0.09, 0.7),title = "Phenotypic Covariance (Log(Mass, 10<sup>-5</sup> g/t))"),
                                            aspectmode='cube'),
                              scene = list(domain=list(x=c(0.50,0.95),y=c(0.25,1)),
                                           xaxis=list(title = "Time"),
                                           yaxis =list(title = "Time"),
                                           zaxis=list(range = c(-0.09, 0.7),title = "Genetic Covariance (Log(Mass, 10<sup>-5</sup> g/t))"),
                                           aspectmode='cube'))

fig_RR3

##########################################################################################

## Fit 3rd-order for genetic and 2nd-order for environmental

fformLegendre2 <- trait ~ -1 + df_leg$phi1 + df_leg$phi2 + df_leg$phi3 + df_leg$phi4 +
  (-1 + df_leg$phi1 + df_leg$phi2 + df_leg$phi3 + df_leg$phi4 | df_leg$id) + 
  (-1 + df_leg$phi1 + df_leg$phi2 + df_leg$phi3  | df_leg$id)

system.time(
  ffLeg2 <- fit_genetic_fmm(fformLegendre2, df_leg, A, c(4,3))
) # user   system  elapsed

isSingular(ffLeg2)
summary(ffLeg2)

## Model results
### Fixed effect
beta_leg2 <- fixef(ffLeg2)
FE_leg2 <- basisfine %*% beta_leg2

par(mfrow = c(1,1))
par(font.main = 1)
plot(c(-1,1), c(0, 3), type = "n", 
     xlab = "Time", 
     ylab = "",
     xlim = c(-1, 1), ylim = c(0,3), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(-1, 1, by = 0.2), pos = 0) 
axis(side = 2, at = seq(0, 3, by = 0.5), pos = 0) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
lines(timegrid, FE_leg2, type = "l", lty = "dashed", cpl = "black")
lines(timegrid, sample_mean, type = "l", col = "red")
mtext("Fixed Effect", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, v=0)
legend("bottomright", legend= c("estimated FE", "sample mean"), lty = c("dashed", "solid"), bty = "n", col = c("black", "red"))

### Random effect
vcLeg2 <- VarCorr(ffLeg2)
CGLeg2 <- vcLeg2[["df_leg.id"]] # genetic covariance
CELeg2 <- vcLeg2[["df_leg.id.1"]] # environmental covariance

CG_funLeg2 <- basisfine%*% CGLeg2 %*% t(basisfine) # estimated gen cov function
CE_funLeg2 <- basisfine[,1:3] %*% CELeg2 %*% t(basisfine[,1:3]) # estimated env cov function

fig3 <- plot_ly(showscale=FALSE, scene='scene1') 
fig3 <- fig3 %>% add_surface(z = ~CG_funLeg2, x = timegrid, y = timegrid)

fig4 <- plot_ly(showscale = FALSE, scene='scene2') 
fig4 <- fig4 %>% add_surface(z = ~CE_funLeg2, x= timegrid, y = timegrid)

fig_RR2 <- subplot(fig3, fig4) 
fig_RR2 <- fig_RR2 %>% layout(scene = list(title = "",
                                           domain=list(x=c(0,0.45),y=c(0.25,1)),
                                           xaxis=list(title = "Time"),
                                           yaxis =list(title = "Time") , 
                                           zaxis=list(range = c(-0.05, 0.5),title = "Genetic Covariance (Log(Mass, 10<sup>-5</sup> g/t))"),
                                           aspectmode='cube'),
                              scene2 = list(domain=list(x=c(0.50,0.95),y=c(0.25,1)),
                                            xaxis=list(title = "Time"),
                                            yaxis =list(title = "Time"),
                                            zaxis=list(range = c(-0.05, 0.5), title = "Environmental Covariance (Log(Mass, 10<sup>-5</sup> g/t))"),
                                            aspectmode='cube'))

fig_RR2

P_funLeg2 <- CG_funLeg2 + CE_funLeg2
fig5 <- plot_ly(showscale = FALSE, scene='scene2') 
fig5 <- fig5 %>% add_surface(z = ~P_funLeg2, x= timegrid, y = timegrid)

fig_RR3 <- subplot(fig5, fig3) 
fig_RR3 <- fig_RR3 %>% layout(title = "",
                              scene2 = list(domain=list(x=c(0,0.45),y=c(0.25,1)),
                                            xaxis=list(title = "Time"),
                                            yaxis =list(title = "Time") , 
                                            zaxis=list(range = c(-0.09, 0.7),title = "Phenotypic Covariance (Log(Mass, 10<sup>-5</sup> g/t))"),
                                            aspectmode='cube'),
                              scene = list(domain=list(x=c(0.50,0.95),y=c(0.25,1)),
                                           xaxis=list(title = "Time"),
                                           yaxis =list(title = "Time"),
                                           zaxis=list(range = c(-0.09, 0.7),title = "Genetic Covariance (Log(Mass, 10<sup>-5</sup> g/t))"),
                                           aspectmode='cube'))

fig_RR3

##################################################################################
## Fit 3rd-order for genetic and 2nd-order for environmental with time scales
## all start from -1 and end at 1

## Resale time s.t. all starts at -1 and ends at 1
x1 <- list()
for (i in 1:N){
  x1[[i]] = 2*(age_list[[i]]-min(age_list[[i]]))/(max(age_list[[i]])-min(age_list[[i]])) -1 
}

legendre_polys3 <- legendre.polynomials(n = 3, normalized = FALSE)
legBasis3_new <- list() # 3rd order polynomial (non-normalised)
for (i in 1:N){
  values <- polynomial.values(polynomials = legendre_polys3, x = x1[[i]])
  legBasis3_new[[i]] <- do.call(cbind, values)
}
basis_new <- do.call(rbind, legBasis3_new)

par(mfrow = c(1,2))
plot(c(-1,1), c(-2,2), type = "n",
     xlab = "Time", 
     ylab = "",
     xlim = c(-1, 1), ylim = c(-2, 2), 
     xaxs = "i", yaxs = "i",
     axes = FALSE)
axis(side = 1, at = seq(-1, 1, by = 0.1), pos = 0) 
axis(side = 2, at = seq(-2, 2, by = 0.5), pos = 0) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:4){
  lines(x1[[1]], legBasis3_new[[1]][,i], col = i)
}
mtext("Orthogonal Legendre Poly of order 3 (Subject 1)", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, v=0)

plot(c(-1,1), c(-2,2), type = "n",
     xlab = "Time", 
     ylab = "",
     xlim = c(-1, 1), ylim = c(-2, 2), 
     xaxs = "i", yaxs = "i",
     axes = FALSE)
axis(side = 1, at = seq(-1, 1, by = 0.1), pos = 0) 
axis(side = 2, at = seq(-2, 2, by = 0.5), pos = 0) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:4){
  lines(x1[[2]], legBasis3_new[[2]][,i], col = i)
}
mtext("Orthogonal Legendre Poly of order 3 (Subject 2)", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, v=0)
################################################################################################

df_leg_new <-  data.frame(id = id, 
                      trait = df$logtrait, 
                      phi1 = basis_new[,1],
                      phi2 = basis_new[,2],
                      phi3 = basis_new[,3],
                      phi4 = basis_new[,4])

fformLegendre3 <- trait ~ -1 + df_leg_new$phi1 + df_leg_new$phi2 + df_leg_new$phi3 + df_leg_new$phi4 +
  (-1 + df_leg_new$phi1 + df_leg_new$phi2 + df_leg_new$phi3 + df_leg_new$phi4 | df_leg_new$id) + 
  (-1 + df_leg_new$phi1 + df_leg_new$phi2 + df_leg_new$phi3  | df_leg_new$id)

system.time(
  ffLeg_new <- fit_genetic_fmm(fformLegendre3, df_leg_new, A, c(4,3))
) # user   system  elapsed

isSingular(ffLeg_new)
summary(ffLeg_new)

## Model results
### Fixed effect
logFPCA_new <- FPCA(Ly = log_trait_list, Lt = x1, optns = list(nRegGrid = 100)) # FPCA through conditional expectation
timegrid_new <- logFPCA_new$workGrid
sample_mean_new <- logFPCA_new$mu

valuesfine_new <- polynomial.values(polynomials = legendre_polys3, x = timegrid_new)
basisfine_new <- do.call(cbind, valuesfine_new)
beta_leg_new <- fixef(ffLeg_new)
FE_leg_new <- basisfine_new %*% beta_leg_new

par(mfrow = c(1,1))
par(font.main = 1)
plot(c(-1,1), c(0, 3), type = "n", 
     xlab = "Time", 
     ylab = "",
     xlim = c(-1, 1), ylim = c(0,3), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(-1, 1, by = 0.2), pos = 0) 
axis(side = 2, at = seq(0, 3, by = 0.5), pos = 0) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
lines(timegrid_new, FE_leg_new, type = "l", lty = "dashed", cpl = "black")
lines(timegrid_new, sample_mean_new, type = "l", col = "red")
mtext("Fixed Effect", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, v=0)
legend("bottomright", legend= c("estimated FE", "sample mean"), lty = c("dashed", "solid"), bty = "n", col = c("black", "red"))

### Random effect
vcLeg_new <- VarCorr(ffLeg_new)
CGLeg_new <- vcLeg_new[["df_leg_new.id"]] # genetic covariance
CELeg_new <- vcLeg_new[["df_leg_new.id.1"]] # environmental covariance

CG_funLeg_new <- basisfine_new%*% CGLeg_new %*% t(basisfine_new) # estimated gen cov function
CE_funLeg_new <- basisfine_new[,1:3] %*% CELeg_new %*% t(basisfine_new[,1:3]) # estimated env cov function

fig3 <- plot_ly(showscale=FALSE, scene='scene1') 
fig3 <- fig3 %>% add_surface(z = ~CG_funLeg_new, x = timegrid_new, y = timegrid_new)

fig4 <- plot_ly(showscale = FALSE, scene='scene2') 
fig4 <- fig4 %>% add_surface(z = ~CE_funLeg_new, x= timegrid_new, y = timegrid_new)

fig_RR2 <- subplot(fig3, fig4) 
fig_RR2 <- fig_RR2 %>% layout(scene = list(title = "",
                                           domain=list(x=c(0,0.45),y=c(0.25,1)),
                                           xaxis=list(title = "Time"),
                                           yaxis =list(title = "Time") , 
                                           zaxis=list(range = c(-0.01, 0.1),title = "Genetic Covariance (Log(Mass, 10<sup>-5</sup> g/t))"),
                                           aspectmode='cube'),
                              scene2 = list(domain=list(x=c(0.50,0.95),y=c(0.25,1)),
                                            xaxis=list(title = "Time"),
                                            yaxis =list(title = "Time"),
                                            zaxis=list(range = c(-0.001, 0.01), title = "Environmental Covariance (Log(Mass, 10<sup>-5</sup> g/t))"),
                                            aspectmode='cube'))

fig_RR2

P_funLeg_new <- CG_funLeg_new + CE_funLeg_new
fig5 <- plot_ly(showscale = FALSE, scene='scene2') 
fig5 <- fig5 %>% add_surface(z = ~P_funLeg_new, x= timegrid_new, y = timegrid_new)

fig_RR3 <- subplot(fig5, fig3) 
fig_RR3 <- fig_RR3 %>% layout(title = "",
                              scene2 = list(domain=list(x=c(0,0.45),y=c(0.25,1)),
                                            xaxis=list(title = "Time"),
                                            yaxis =list(title = "Time") , 
                                            zaxis=list(range = c(-0.01, 0.1),title = "Phenotypic Covariance (Log(Mass, 10<sup>-5</sup> g/t))"),
                                            aspectmode='cube'),
                              scene = list(domain=list(x=c(0.50,0.95),y=c(0.25,1)),
                                           xaxis=list(title = "Time"),
                                           yaxis =list(title = "Time"),
                                           zaxis=list(range = c(-0.01, 0.1),title = "Genetic Covariance (Log(Mass, 10<sup>-5</sup> g/t))"),
                                           aspectmode='cube'))

fig_RR3
