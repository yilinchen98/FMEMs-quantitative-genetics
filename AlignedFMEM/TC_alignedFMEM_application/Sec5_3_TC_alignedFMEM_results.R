# Fit the TC dataset (log10-transformed) with PCs as basis functions
source("PrepareR.R")
## Load Data
TRFUN25PUP4 <- read.delim("TRFUN25PUP4.DAT",header = FALSE)
names(TRFUN25PUP4) <- c("id","sire","dam","trait","x")
df <- data.frame(TRFUN25PUP4)

## calculate genetic relationship matrix
id <- df$id
FirstUniqueIdPos <- which(duplicated(id) == FALSE)
pos <- id[FirstUniqueIdPos] # extract ids for all subjects
sire_id <- df$sire[FirstUniqueIdPos] # extract ids for sire
dam_id <- df$dam[FirstUniqueIdPos] # extract ids for dam

pede <- editPed(sire = sire_id, dam = dam_id, label = pos)
ped<- with(pede, pedigree(label=label, sire=sire, dam=dam))
A <- getA(ped)[163:1035,163:1035]

## Rescale time interval to [0,1]
N <- length(FirstUniqueIdPos) # N = 873 subjects
n <- length(id) # n = 6860 observations
age_list <- split(df$x,id)
trait_list <- split(df$trait,id)

### x = (x -min(x))/(max(x) - min(x))
age_list_new <- list()
for (i in 1:N){
  age_list_new[[i]] <- (age_list[[i]]-min(age_list[[i]]))/(max(age_list[[i]])-min(age_list[[i]]))
}

df$x_rescaled <- unsplit(age_list_new,id)
df$logtrait <- log10(df$trait)
log_trait_list <- split(df$logtrait,id)

## plot growth curve
par(mfrow = c(1,2))
par(mar = c(5, 6, 4, 2) + 0.1)
plot(c(1, 25), c(0, 400), type = "n", 
     xlab = "Days", 
     ylab = expression(Mass~(10^-5~g)),
     xlim = c(0, 25), ylim = c(0, 400), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
for (i in 1:N){
  lines(age_list[[i]], trait_list[[i]], type = "l", col = i)
}
mtext("Growth Curves", side = 3, adj = 0, line = 1, font = 2)
axis(side = 1, at = seq(0, 25, by = 5), pos = 0, lwd = 2) 
axis(side = 2, at = seq(0, 400, by = 100), pos = 0, lwd = 2) 

par(mar = c(5, 6, 4, 2) + 0.1)
plot(c(0, 1), c(0, 3), type = "n", 
     xlab = "Standardised Time", 
     ylab = expression(Log(Mass~(10^-5~g))),
     xlim = c(0, 1), ylim = c(0, 3), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
for (i in 1:N){
  lines(age_list_new[[i]], log_trait_list[[i]], type = "l", col = i)
}
mtext("Growth Curves (log scale)", side = 3, adj = 0, line = 1, font = 2)
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0, lwd = 2) 
axis(side = 2, at = seq(0, 3, by = 1), pos = 0, lwd = 2) 
########################################################################################

## Smooth data, choose smoothing parameter <= 1e-4
timefine <- seq(0,1,length=100) # dense time grid
mass_smoothed <-list()
pred_mass_fine <- matrix(0,100,N) # store the smoothed mass predicted on the dense grid
pred_logmass_fine <- matrix(0,100,N) # store the smoothed logmass
lam <- rep(0,N) # smoothing parameter used for each log growth curve

for (i in 1:N) {
  ss_logmass <- smooth.spline(age_list_new[[i]], log10(trait_list[[i]]), cv=FALSE, all.knots=TRUE)
  
  # Check if lambda is greater than 1e-4
  if (ss_logmass$lambda > 1e-4) {
    # Redo smoothing with lambda set to 1e-4
    ss_logmass <- smooth.spline(age_list_new[[i]], log10(trait_list[[i]]), lambda=1e-4, all.knots=TRUE)
  }
  
  mass_smoothed[[i]] <- 10^(ss_logmass$y)
  pred_logmass_fine[,i] <- predict(ss_logmass, timefine)$y
  pred_mass_fine[,i] <- 10^(predict(ss_logmass, timefine)$y) ## predict smoothed data on a regular dense grid
  lam[i] <- ss_logmass$lambda
}

## Align Log Growth Curves
aligned_logmass_process <- time_warping(pred_logmass_fine, timefine)
aligned_logmass_curve <- aligned_logmass_process$fn
aligned_logmass_mean <- aligned_logmass_process$fmean
warping_logmass_funs <- aligned_logmass_process$warping_functions

## FPCA
fpcaobj_logmass <- prcomp(x=t(aligned_logmass_curve), retx = TRUE, center = TRUE, rank. = 4)
pcs_logmass <- fpcaobj_logmass$rotation # eigen vectors


par(mfrow = c(1,1))
pl <- fviz_eig(fpcaobj_logmass, 
               addlabels = TRUE, 
               main="Aligned Data", ncp = 6)
pl# scree plot

par(mfrow = c(1, 2))
par(font.main = 1)
plot(c(0,1), c(0, 1), type = "n", 
     xlab = "Standardised Time", 
     ylab = "",
     xlim = c(0, 1), ylim = c(-0.1,0.2), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
for (i in 1:N){
  lines(timefine, pcs_logmass[,1], type = "l", col = i)
}
mtext("Principal Component 1 of Aligned Curves", side = 3, adj = 0, line = 1, font = 2)
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0, lwd = 2) 
axis(side = 2, at = seq(-0.1, 0.2, by = 0.05), pos = 0, lwd = 2) 

par(font.main = 1)
plot(c(0,1), c(0, 1), type = "n", 
     xlab = "Standardised Time", 
     ylab = "",
     xlim = c(0, 1), ylim = c(-0.1,0.2), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
for (i in 1:N){
  lines(timefine, pcs_logmass[,2], type = "l", col = i)
}
mtext("Principal Component 2 of Aligned Curves", side = 3, adj = 0, line = 1, font = 2)
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0, lwd = 2) 
axis(side = 2, at = seq(-0.1, 0.2, by = 0.05), pos = 0, lwd =2 )
################################################################################
## Model Fitting (Log scale)

### Align raw data
aligned_logtrait <- list()
gamma_logmass <- list()
for (i in 1:N){
  gamma_logmass_inter <- convert_to_basisfunctions(timefine, warping_logmass_funs[,i], age_list_new[[i]])
  gamma_logmass[[i]] <- gamma_logmass_inter
  aligned_logtrait[[i]] <- warp_f_gamma(log10(trait_list[[i]]), age_list_new[[i]], gamma_logmass_inter)
}

par(mfrow = c(1, 2))
par(mar = c(5, 6, 4, 2) + 0.1, font.main = 1)
plot(c(0,1), c(0, 1), type = "n", 
     xlab = "Standardised Time", 
     ylab = "Warped Time",
     xlim = c(0, 1), ylim = c(0,1), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
for (i in 1:N){
  lines(age_list_new[[i]], gamma_logmass[[i]], type = "l", col = i)
}
mtext("Warping Functions", side = 3, adj = 0, line = 1, font = 2)
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0, lwd = 2) 
axis(side = 2, at = seq(0, 1, by = 0.1), pos = 0, lwd = 2) 

plot(c(0,1), c(0, 3), type = "n", 
     xlab = "Standardised Time", 
     ylab = expression(Log(Mass,~10^-5~g)),
     xlim = c(0, 1), ylim = c(0,3), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
for (i in 1:N){
  lines(age_list_new[[i]], aligned_logtrait[[i]], type = "l", col = i)
}
mtext("Aligned Raw Log Mass Curves", side = 3, adj = 0, line = 1, font = 2)
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0, lwd = 2) 
axis(side = 2, at = seq(0, 3, by = 0.5), pos = 0, lwd = 2) 
########################################################################################

### Prepare model basis
phi_logmass_list <- list() 
# create an empty list which stores eigenfunctions for 873 subjects
# evaluated at the original time points.

for (i in 1:N){
  phi_logmass <- convert_to_basisfunctions(t = timefine, eigenvecs = pcs_logmass,
                                           tout = age_list_new[[i]])
  phi_logmass_list[[i]] <- phi_logmass
}

phi_logmass <- do.call(rbind,phi_logmass_list)
colnames(phi_logmass) <- c("phi1", "phi2", "phi3", "phi4")

## Reform dataframe
df_logmass <- data.frame(id = id, trait = unsplit(aligned_logtrait,id), phi_logmass)

fmmFormL <- trait ~ -1 + df_logmass$phi1 + df_logmass$phi2 + df_logmass$phi3 + 
  (-1 + df_logmass$phi1 + df_logmass$phi2 | df_logmass$id) + 
  (-1 + df_logmass$phi1 + df_logmass$phi2 | df_logmass$id) 

system.time(
  ffL <- fit_genetic_fmm(fmmFormL, df_logmass, A, 2)
) # user   system  elapsed

summary(ffL)

isSingular(ffL) # original tol = 1e-4, not singular fit

## Model results
### Random effects
vcL <- VarCorr(ffL)
CGL <- vcL[["df_logmass.id"]] # genetic covariance
CEL <- vcL[["df_logmass.id.1"]] # environmental covariance

CG_funL <- pcs_logmass[,1:2] %*% CGL %*% t(pcs_logmass[,1:2]) # estimated gen cov function
CE_funL <- pcs_logmass[,1:2] %*% CEL %*% t(pcs_logmass[,1:2]) # estimated env cov function

P_funL <- CG_funL + CE_funL # estimated phenotypic covariance function
################################################################################

### SCB
CG_fun_sapL_300 <- readRDS("TC_CG_fun_bs_samples.rds")
mean_CGL <- apply(CG_fun_sapL_300, c(1, 2), mean)
CGL_diff <- sapply(seq_len(300), function(i){
  max(abs(CG_fun_sapL_300[,,i] - mean_CGL))
})
CGL_sup <- quantile(CGL_diff, probs = 0.95) # quantile
CGL_upper <- mean_CGL + CGL_sup
CGL_lower <- mean_CGL - CGL_sup

fig_P <- plot_ly(showscale = FALSE, scene='scene') 
fig_P <- fig_P %>% add_surface(z = ~P_funL, x= timefine, y = timefine)

fig_CG <- plot_ly(showscale=FALSE, scene='scene2') 
fig_CG <- fig_CG %>% add_surface(z = ~CG_funL, x = timefine, y = timefine)
fig_CG <- fig_CG %>% add_surface(z = ~CGL_upper, x = ~timefine, y = ~timefine, 
                             opacity = 0.5,colorscale = list(c(0, 1), c("grey", "lightgrey")))
fig_CG <- fig_CG %>% add_surface(z = ~CGL_lower, x = ~timefine, y = ~timefine, 
                             opacity = 0.5,colorscale = list(c(0, 1), c("grey", "lightgrey")))

fig_RR <- subplot(fig_P, fig_CG) 
fig_RR <- fig_RR %>% layout(title = "",
                              scene = list(domain=list(x=c(0,0.45),y=c(0.25,1)),
                                            xaxis=list(title = "Standardised Time"),
                                            yaxis =list(title = "Standardised Time") , 
                                            zaxis=list(range = c(-0.015,0.065),title = "Phenotypic Covariance (Log(Mass, 10<sup>-5</sup> g/t))"),
                                            aspectmode='cube'),
                              scene2 = list(domain=list(x=c(0.50,0.95),y=c(0.25,1)),
                                           xaxis=list(title = "Standardised Time"),
                                           yaxis =list(title = "Standardised Time"),
                                           zaxis=list(range = c(-0.015,0.065),title = "Genetic Covariance (Log(Mass, 10<sup>-5</sup> g/t))"),
                                           aspectmode='cube'))

fig_RR
################################################################################
### Genetic eigenfunction
eigenfunL_CG1 <- eigen(CG_funL)$vectors[,1]

GeigenL_300 <-  matrix(0,nrow = 100, ncol = 300)
for(i in 1:300){
  bs_GeigenL <- eigen(CG_fun_sapL_300[,,i])$vectors[,1]
  if(trapz(timefine, eigenfunL_CG1 * bs_GeigenL) < 0) bs_GeigenL <- -bs_GeigenL
  GeigenL_300[,i] <- bs_GeigenL
} # compute SCB
mean_GEL <- rowMeans(GeigenL_300)
GeigenL_diff <- sapply(seq_len(300), function(i){
  max(abs(GeigenL_300[,i] - mean_GEL))
})
GeigenL_sup <- quantile(GeigenL_diff, probs = 0.95) 

par(mfrow = c(1,1))
plot(c(0,1), c(-0.25, 0.1), type = "n", 
     xlab = "Standardised Time", 
     ylab = "",
     xlim = c(0, 1), ylim = c(-0.25,0.1), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
lines(timefine, eigenfunL_CG1, type = "l", lty = "solid", col = "red", lwd = 1.5)
lines(timefine, mean_GEL + GeigenL_sup, lty = "dashed")
lines(timefine, mean_GEL - GeigenL_sup, lty = "dashed")
mtext("First Genetic Eigenfunction", side = 3, adj = 0, line = 1, font = 2)
axis(side = 1, at = seq(0, 1, by = 0.1), pos = -0.25, lwd = 2) 
axis(side = 2, at = seq(-0.25, 0.1, by = 0.05), pos = 0, lwd = 2)
################################################################################

# Impact of time alignment (fit the mixed-effects model using misaligned data)

## FPCA of Misaligned curves
fpcaobj_lg <- prcomp(t(pred_logmass_fine), rank. = 8)
pcs_lg <- fpcaobj_lg$rotation
par(mfrow = c(1,1))
pcl <- fviz_eig(fpcaobj_lg, 
                addlabels = TRUE, 
                main="Misaligned Data", ncp = 10)
pcl
################################################################################

## Model fitting: fit original data

### prepare basis functions
phi_lg_list <- list() 

for (i in 1:N){
  phi_lg <- convert_to_basisfunctions(t = timefine, eigenvecs = pcs_lg,
                                      tout = age_list_new[[i]])
  phi_lg_list[[i]] <- phi_lg
}

phi_lg <- do.call(rbind,phi_lg_list)
colnames(phi_lg) <- c("phi1", "phi2", "phi3", "phi4", "phi5", "phi6", "phi7", "phi8")

### Reform dataframe
df_lg <- data.frame(id = id, trait = df$logtrait, phi_lg)

fmmFormlg <- trait ~ -1 + df_lg$phi1 + df_lg$phi2 + df_lg$phi3 + df_lg$phi4 + df_lg$phi5 + df_lg$phi6 + df_lg$phi7 + df_lg$phi8 +
  (-1 + df_lg$phi1 + df_lg$phi2 + df_lg$phi3 + df_lg$phi4 | df_lg$id) + 
  (-1 + df_lg$phi1 + df_lg$phi2 + df_lg$phi3 + df_lg$phi4 | df_lg$id) 

system.time(
  fflg <- fit_genetic_fmm(fmmFormlg, df_lg, A, 4)
) # user   system  elapsed

summary(fflg)

isSingular(fflg) # original tol = 1e-4, not singular

## Model results
### Random effects

vclg <- VarCorr(fflg)
CGlg <- vclg[["df_lg.id"]] # genetic covariance
CElg <- vclg[["df_lg.id.1"]] # environmental covariance

CG_funlg <- pcs_lg[,1:4] %*% CGlg %*% t(pcs_lg[,1:4]) # estimated gen cov function
CE_funlg <- pcs_lg[,1:4] %*% CElg %*% t(pcs_lg[,1:4]) # estimated env cov function

P_funlg <- CG_funlg + CE_funlg # estimated phenotypic covariance function

### Compare sample covariance for aligned and misaligned curves

sam_orin <- GetCovSurface(log_trait_list, age_list_new, optns = list(nRegGrid = 100)) # smooth kernel: Gauss
sam_registered <- GetCovSurface(aligned_logtrait, age_list_new, optns = list(nRegGrid = 100))

cov_sam_orin <- sam_orin$cov # sample cov of original data
cov_sam_registered <- sam_registered$cov # sample cov of the aligned data

fig1 <- plot_ly(showscale = FALSE, scene='scene') 
fig1 <- fig1 %>% add_surface(z = ~P_funlg, x= timefine, y = timefine)

fig2 <- plot_ly(showscale=FALSE, scene='scene2') 
fig2 <- fig2 %>% add_surface(z = ~cov_sam_orin, x = timefine, y = timefine)

fig3 <- plot_ly(showscale = FALSE, scene='scene3') 
fig3 <- fig3 %>% add_surface(z = ~cov_sam_registered, x= timefine, y = timefine)

fig_cov <- subplot(fig1, fig2, fig3)
fig_cov <- fig_cov %>% layout(scene = list(domain=list(x=c(0,0.33),y=c(0,1)),
                                            xaxis=list(title = "Standardised Time"),
                                            yaxis =list(title = "Standardised Time") , 
                                            zaxis=list(range=c(-0.001,0.1),title = "Phenotypic Covariance (Log(Mass, 10<sup>-5</sup> g/t))"),
                                            aspectmode='cube'),
                              scene2 = list(domain=list(x=c(0.33,0.67),y=c(0,1)),
                                            xaxis=list(title = "Standardised Time"),
                                            yaxis =list(title = "Standardised Time"),
                                            zaxis=list(range=c(-0.001,0.1),title = "Original (Log(Mass, 10<sup>-5</sup> g/t))"),
                                            aspectmode='cube'),
                              scene3 = list(domain=list(x=c(0.67,1),y=c(0,1)),
                                            xaxis=list(title = "Standardised Time"),
                                            yaxis =list(title = "Standardised Time"),
                                            zaxis=list(range=c(-0.001,0.1),title = "Aligned Raw (Log(Mass, 10<sup>-5</sup> g/t))"),
                                            aspectmode='cube')
  
)
fig_cov
