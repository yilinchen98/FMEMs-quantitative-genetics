# Here we find the joint model for analysis of amplitude and warping functions.
# We apply the combined FPCA proposed by Lee and Jung in 2016 for modeling genetic data, and
# in particular, we extend the combined FPCA into a functional mixed-effect framework for estimating
# genetic covariance function. 
source("PrepareR.R")
## Load data
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

## Transform to the log scale
N = length(FirstUniqueIdPos) # N = 873 subjects
n = length(id) # n = 6860 observations
age_list <- split(df$x,id)
trait_list <- split(df$trait,id)
df$logtrait <- log10(df$trait)
log_trait_list <- split(df$logtrait,id)

# In FDA, curves are often measured on the same time interval. Here we rescale
# time so that all curves are measured from 0 (hatching) to 1 (pupation). 

### x = (x -min(x))/(max(x) - min(x)) 
age_list_new <- list()
for (i in 1:N){
  age_list_new[[i]] = (age_list[[i]]-min(age_list[[i]]))/(max(age_list[[i]])-min(age_list[[i]]))
}
df$x_rescaled <- unsplit(age_list_new,id)

# The length of the larval period for each individual. We will later add the larval
# period for the genetic mixed-effects modelling.
x_length <- rep(0,length = N)
for (i in 1:N){
  x_length[i] <- max(age_list[[i]])-min(age_list[[i]])
}

## Data smoothing
grids <- 60
timefine <- seq(0,1,length = grids) # define a fine time grid
### choose smoothing parameter <= 1e-4
pred_logmass_fine <- matrix(0,grids,N) # store the smoothed logmass M x N on a grid of size M
lam <- rep(0,N) # smoothing parameter used for each log growth curve
for (i in 1:N) {
  ss_logmass <- smooth.spline(age_list_new[[i]], log_trait_list[[i]], cv=FALSE, all.knots=TRUE)
  
  # Check if lambda is greater than 1e-4
  if (ss_logmass$lambda > 1e-4){
    ss_logmass <- smooth.spline(age_list_new[[i]], log_trait_list[[i]], lambda=1e-4, all.knots=TRUE)
  }
  
  pred_logmass_fine[,i] <- predict(ss_logmass, timefine)$y # predict the smoothed curves on the regular dense grid
  lam[i] <- ss_logmass$lambda
}

## Registration based on Fisher-Rao Metric
aligned_FR <- time_warping(pred_logmass_fine, time = timefine)
y_FR <- aligned_FR$fn
y_FR_mean <- aligned_FR$fmean
h_FR <- aligned_FR$warping_functions

## Transform amplitude and warping functions
y_FR_centeredMean <- y_FR - y_FR_mean # centred the amplitude curves
clr_h_trans_FR <- gam_to_h(h_FR, smooth = FALSE) # centred log-ratio transform of warping functions

## Subsample
### On average there are 8 points per subject. Use linear interpolation to map
### to 8 uniform points in the unit interval.

## sampling scheme from the data
sampling_points <- sapply(log_trait_list, length)
max_points <- max(sampling_points) # 14 points per subject
min_points <- min(sampling_points) # 5 points per subject
ave_points <- round(n/N) # average points per subject = 8

t <- seq(0,1, length = 8)
amp <- convert_to_basisfunctions(t = timefine, eigenvecs = y_FR_centeredMean, tout = t)
gam <- convert_to_basisfunctions(t = timefine, eigenvecs = clr_h_trans_FR, tout = t)
################################################################################

# For the combined FPCA, we choose the scaling parameter C, which glues the amplitude functions
# and the transformed warping functions such that the ratio of the covariance of amplitude and warping 
# functions equals to 1.

## Compute the covariance ratio between amplitude and warping functions
## Choose the scaling parameter C: match the covariance ratio of amplitude and warping functions

cov_amp <- cov(t(amp))
cov_gam <- cov(t(gam))
norm_cov_amp <- norm(cov_amp, type ="M")
norm_cov_gam <- norm(cov_gam, type ="M")
C <- sqrt(norm_cov_amp/norm_cov_gam)

# The joint curves
## We combine the amplitude, warping functions with the length of the larval period (which is a scalar). 
## However, the scale of the larval period is not at the same order of the amplitude functions.
## We, therefore, scale the larval period length w.r.t. the variability of the amplitude functions.

amp_std <- apply(amp, 1, sd) # sd of the amplitude functions
gam_std <- apply(C*gam, 1, sd) # sd of the clr-transform warping functions
xLen <- x_length - mean(x_length) # zero-mean
xLen_std <- sd(xLen) # sd of the length vector
ave_amp_std <- mean(amp_std) # average sd of the amplitude
ave_gam_std <- mean(gam_std) # average sd of the clr-transform warping

C_xLen <- ave_amp_std/xLen_std # the scaling parameter for the larval period length
xLen_scaled <- C_xLen * xLen
xLen_scaled_std <- sd(xLen_scaled)

g <- rbind(amp, xLen_scaled, C*gam) # joint curves
g_mean <- rowMeans(g) # around zero
sample_cov <- cov(t(g))
#t_new <- seq(0,2,length = 17)
################################################################################
par(mfrow = c(1,3), bty = "n")
par(mar = c(5, 6, 4, 2) + 0.1)
plot(c(0,1), c(-1, 1.5), type = "n", 
     xlab = "Standardised Time", 
     ylab = expression(Log(Mass,~10^-5~g)),
     xlim = c(0, 1), ylim = c(-1,1.5), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
for (i in 1:N){
  lines(t, amp[,i], type = "l", col = i)
}
mtext("Amplitude", side = 3, adj = 0, line = 1, font = 2)
axis(side = 1, at = seq(0, 1, by = 0.1), pos = -1, lwd = 2) 
axis(side = 2, at = seq(-1, 1.5, by = 0.5), pos = 0, lwd = 2)

par(mar = c(5, 6, 4, 2) + 0.1)
plot(c(0,1), c(-1, 1.5), type = "n", 
     xlab = "Standardised Time", 
     ylab = "",
     xlim = c(0, 1), ylim = c(-1,1.5), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
for (i in 1:N){
  lines(t, C*gam[,i], type = "l", col = i)
}
mtext("Scaled CLR-Warping (8 Points)", side = 3, adj = 0, line = 1, font = 2)
axis(side = 1, at = seq(0, 1, by = 0.1), pos = -1, lwd = 2) 
axis(side = 2, at = seq(-1, 1.5, by = 0.5), pos = 0, lwd = 2)

par(mar = c(5, 6, 4, 2) + 0.1)
boxplot(xLen_scaled, ylim = c(-1,1.5),xaxs = "i", yaxs = "i")
axis(2, at = seq(-1, 1.5, by = 0.5),lwd = 2)
mtext("Scaled Larval Period", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
abline(h = 0, lwd = 2)
################################################################################

## Combined FPCA
fpca_obj <- prcomp(t(g), center = TRUE, retx = TRUE)
fpcs <- fpca_obj$rotation
################################################################################
par(mfrow = c(1,1))
screeplt <- fviz_eig(fpca_obj, 
                     addlabels = TRUE, 
                     main="Combined FPCA")
screeplt
################################################################################

par(mfrow = c(2,3))
for(k in 1:2){
  par(font.main = 1)
  plot(c(0,1), c(-1, 1), type = "n", 
       xlab = "Standardised Time", 
       ylab = "",
       xlim = c(0, 1), ylim = c(-1,1), 
       xaxs = "i", yaxs = "i",
       axes = FALSE) 
  grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
  lines(t, fpcs[1:8,k], type = "l")
  mtext(paste("PC", k, ": Amplitude", sep = ""), side = 3, adj = 0, line = 1, font = 2)
  axis(side = 1, at = seq(0, 1, by = 0.2), pos = 0, lwd = 2) 
  axis(side = 2, at = seq(-1, 1, by = 0.2), pos = 0, lwd = 2) 
  
  par(font.main = 1)
  plot(c(0,1), c(-1, 1), type = "n", 
       xlab = "Standardised Time", 
       ylab = "",
       xlim = c(0, 1), ylim = c(-1,1), 
       xaxs = "i", yaxs = "i",
       axes = FALSE) 
  grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
  lines(t, fpcs[10:17,k], type = "l")
  mtext(paste("PC", k, ": Phase", sep = ""), side = 3, adj = 0, line = 1, font = 2)
  axis(side = 1, at = seq(0, 1, by = 0.2), pos = 0, lwd = 2) 
  axis(side = 2, at = seq(-1, 1, by = 0.2), pos = 0, lwd = 2) 
  
  plot(c(0,1), c(-1, 1), type = "n", 
       xlab = "", 
       ylab = "",
       xlim = c(0, 1), ylim = c(-1,1), 
       xaxs = "i", yaxs = "i",
       axes = FALSE) 
  grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
  segments(x0 = 0.5, y0 = 0, x1 = 0.5, y1 = fpcs[9,k], col = "black", lwd = 1)  # Stick
  points(0.5, fpcs[9,k], pch = 16, col = "blue", cex = 1.5)
  mtext(paste("PC", k, ": Larval Period", sep = ""), side = 3, adj = 0, line = 1, font = 2)
  abline(h = 0,lwd = 2)
  axis(side = 2, at = seq(-1, 1, by = 0.2), pos = 0, lwd = 2) 
}
################################################################################

## Model Fitting
#df_combined <- data.frame(id = rep(pos,each = 17) , 
#trait = c(g), 
#phi1 = fpcs[,1],
#phi2 = fpcs[,2],
#phi3 = fpcs[,3],
#phi4 = fpcs[,4],
#phi5 = fpcs[,5],
#phi6 = fpcs[,6])

## Fit 5 FPCs
#fform5 <- trait ~ -1 + df_combined$phi1 + df_combined$phi2 + df_combined$phi3 + df_combined$phi4 + df_combined$phi5 + df_combined$phi6 +
#(-1 + df_combined$phi1 + df_combined$phi2 + df_combined$phi3 + df_combined$phi4 + df_combined$phi5 | df_combined$id) + 
#(-1 + df_combined$phi1 + df_combined$phi2 + df_combined$phi3 + df_combined$phi4 + df_combined$phi5 | df_combined$id) 
#system.time(
#ff5 <- fit_genetic_fmm(fform5, df_combined, A, 5)
#)

#summary(ff5)
#isSingular(ff5) # FALSE
################################################################################

## Model Results
CG_funComb_hat <- readRDS("JointFPCA_results/CG_funComb_hat.rds")
CE_funComb_hat <- readRDS("JointFPCA_results/CE_funComb_hat.rds")

P_funComb_hat <- CG_funComb_hat + CE_funComb_hat
################################################################################

## SCB for the joint model
CG_funcomb300 <- readRDS("JointFPCA_results/CG_funcomb300.rds")
CE_funcomb300 <- readRDS("JointFPCA_results/CE_funcomb300.rds")

Nbs <- 300 # use 300 non-singular bootstrap replicates.
t17 <- seq(0,1,length = 17)
GEign_funcomb_hat <- eigen(CG_funComb_hat)$vectors[,1]

GEign_funcomb300 <- matrix(0,17,Nbs)
for(i in 1:Nbs){
  Geigen_comb <- eigen(CG_funcomb300[,,i])$vectors[,1]
  if(trapz(t17, GEign_funcomb_hat * Geigen_comb) < 0) Geigen_comb <- -Geigen_comb
  GEign_funcomb300[,i] <- Geigen_comb
}

### covariance functions
mean_CG_funcomb <- apply(CG_funcomb300, c(1, 2), mean)
mean_CE_funcomb <- apply(CE_funcomb300, c(1, 2), mean)

CG_funcomb_diff <- rep(0,Nbs)
CE_funcomb_diff<- rep(0,Nbs)

for (i in 1:Nbs){
  CG_funcomb_diff[i] <- max(abs(CG_funcomb300[,,i] - mean_CG_funcomb))
  CE_funcomb_diff[i] <- max(abs(CE_funcomb300[,,i] - mean_CE_funcomb))
}

CG_funcomb_diff_sup <- quantile(CG_funcomb_diff, probs = 0.95)
CE_funcomb_diff_sup <- quantile(CE_funcomb_diff, probs = 0.95)

CG_funcomb_upper <- mean_CG_funcomb + CG_funcomb_diff_sup
CG_funcomb_lower <- mean_CG_funcomb - CG_funcomb_diff_sup

CE_funcomb_upper <- mean_CE_funcomb + CE_funcomb_diff_sup
CE_funcomb_lower <- mean_CE_funcomb - CE_funcomb_diff_sup

### Eigenfunctions of the genetic covariance function
mean_GEign_funcomb <- rowMeans(GEign_funcomb300)
GEign_funcomb_diff <- rep(0,Nbs)

for (i in 1:Nbs){
  GEign_funcomb_diff[i] <- max(abs(GEign_funcomb300[,i] - mean_GEign_funcomb))
}

GEign_funcomb_diff_sup <- quantile(GEign_funcomb_diff, probs = 0.95)

GEign_funcomb_upper <- mean_GEign_funcomb + GEign_funcomb_diff_sup
GEign_funcomb_lower <- mean_GEign_funcomb - GEign_funcomb_diff_sup

## Random effects

### covariance of amplitude functions

fig_amp1 <-  plot_ly(showscale=FALSE, scene='scene') 
fig_amp1 <- fig_amp1 %>% add_surface(z = ~sample_cov[1:8,1:8], x = t, y = t)

fig_amp2 <-  plot_ly(showscale=FALSE, scene='scene2') 
fig_amp2 <- fig_amp2 %>% add_surface(z = ~P_funComb_hat[1:8,1:8], x = t, y = t)

fig_amp3 <-  plot_ly(showscale=FALSE, scene='scene3') 
fig_amp3 <- fig_amp3 %>% add_surface(z = ~CG_funComb_hat[1:8,1:8], x = t, y = t)
fig_amp3 <- fig_amp3 %>% add_surface(z = ~CG_funcomb_upper[1:8,1:8], x = ~t, y = ~t, 
                                     opacity = 0.5, colorscale = list(c(0,1), c("grey", "lightgrey")))
fig_amp3 <- fig_amp3 %>% add_surface(z = ~CG_funcomb_lower[1:8,1:8], x = ~t, y = ~t, 
                                     opacity = 0.5, colorscale = list(c(0,1), c("grey", "lightgrey")))

fig_amp4 <-  plot_ly(showscale=FALSE, scene='scene4') 
fig_amp4 <- fig_amp4 %>% add_surface(z = ~CE_funComb_hat[1:8,1:8], x = t, y = t)
fig_amp4 <- fig_amp4 %>% add_surface(z = ~CE_funcomb_upper[1:8,1:8], x = ~t, y = ~t, 
                                     opacity = 0.5, colorscale = list(c(0,1), c("grey", "lightgrey")))
fig_amp4 <- fig_amp4 %>% add_surface(z = ~CE_funcomb_lower[1:8,1:8], x = ~t, y = ~t, 
                                     opacity = 0.5, colorscale = list(c(0,1), c("grey", "lightgrey")))


fig_ampcov <- subplot(fig_amp1, fig_amp2, fig_amp3, fig_amp4)
fig_ampcov <- fig_ampcov %>% layout(scene = list(title = "",
                                                 domain=list(x=c(0,0.45),y=c(0.55,1)),
                                                 xaxis=list(title = 'Standardised Time'),
                                                 yaxis =list(title = 'Standardised Time') , 
                                                 zaxis=list(range = c(-0.02,0.065),title = "Sample Cov-Amp (Log(Mass, 10<sup>-5</sup> g/t))"),
                                                 aspectmode='cube'),
                                    scene2 = list(domain=list(x=c(0.55,1),y=c(0.55,1)),
                                                  xaxis=list(title = 'Standardised Time'),
                                                  yaxis =list(title = 'Standardised Time'),
                                                  zaxis=list(range = c(-0.02,0.065),title = "P Cov-Amp (Log(Mass, 10<sup>-5</sup> g/t))"),
                                                  aspectmode='cube'),
                                    scene3 = list(domain=list(x=c(0,0.45),y=c(0,0.45)),
                                                  xaxis=list(title = 'Standardised Time'),
                                                  yaxis =list(title = 'Standardised Time'),
                                                  zaxis=list(range = c(-0.02,0.065),title = "Gen Cov-Amp (Log(Mass, 10<sup>-5</sup> g/t))"),
                                                  aspectmode='cube'),
                                    scene4 = list(domain=list(x=c(0.55,1),y=c(0,0.45)),
                                                  xaxis=list(title = 'Standardised Time'),
                                                  yaxis =list(title = 'Standardised Time'),
                                                  zaxis=list(range = c(-0.02,0.065),title = "Env Cov-Amp (Log(Mass, 10<sup>-5</sup> g/t))"),
                                                  aspectmode='cube'))
fig_ampcov
################################################################################

### covariance of warping functions
fig_warp1 <-  plot_ly(showscale=FALSE, scene='scene') 
fig_warp1 <- fig_warp1 %>% add_surface(z = ~sample_cov[10:17,10:17], x = t, y = t)

fig_warp2 <-  plot_ly(showscale=FALSE, scene='scene2') 
fig_warp2 <- fig_warp2 %>% add_surface(z = ~P_funComb_hat[10:17,10:17], x = t, y = t)

fig_warp3 <-  plot_ly(showscale=FALSE, scene='scene3') 
fig_warp3 <- fig_warp3 %>% add_surface(z = ~CG_funComb_hat[10:17,10:17], x = t, y = t)
fig_warp3 <- fig_warp3 %>% add_surface(z = ~CG_funcomb_upper[10:17,10:17], x = t, y = t,
                                       opacity = 0.5, colorscale = list(c(0,1), c("grey", "lightgrey")))
fig_warp3 <- fig_warp3 %>% add_surface(z = ~CG_funcomb_lower[10:17,10:17], x = t, y = t,
                                       opacity = 0.5, colorscale = list(c(0,1), c("grey", "lightgrey")))

fig_warp4 <-  plot_ly(showscale=FALSE, scene='scene4') 
fig_warp4 <- fig_warp4 %>% add_surface(z = ~CE_funComb_hat[10:17,10:17], x = t, y = t)
fig_warp4 <- fig_warp4 %>% add_surface(z = ~CE_funcomb_upper[10:17,10:17], x = t, y = t,
                                       opacity = 0.5, colorscale = list(c(0,1), c("grey", "lightgrey")))

fig_warp4 <- fig_warp4 %>% add_surface(z = ~CE_funcomb_lower[10:17,10:17], x = t, y = t,
                                       opacity = 0.5, colorscale = list(c(0,1), c("grey", "lightgrey")))

fig_warpcov <- subplot(fig_warp1, fig_warp2, fig_warp3, fig_warp4)
fig_warpcov <- fig_warpcov %>% layout(scene = list(title = "",
                                                   domain=list(x=c(0,0.45),y=c(0.55,1)),
                                                   xaxis=list(title = 'Standardised Time'),
                                                   yaxis =list(title = 'Standardised Time') , 
                                                   zaxis=list(range = c(-0.02,0.065),title = "Sample Cov-Phase"),
                                                   aspectmode='cube'),
                                      scene2 = list(domain=list(x=c(0.55,1),y=c(0.55,1)),
                                                    xaxis=list(title = 'Standardised Time'),
                                                    yaxis =list(title = 'Standardised Time'),
                                                    zaxis=list(range = c(-0.02,0.065),title = "P Cov-Phase"),
                                                    aspectmode='cube'),
                                      scene3 = list(domain=list(x=c(0,0.45),y=c(0,0.45)),
                                                    xaxis=list(title = 'Standardised Time'),
                                                    yaxis =list(title = 'Standardised Time'),
                                                    zaxis=list(range = c(-0.02,0.065),title = "Gen Cov-Phase"),
                                                    aspectmode='cube'),
                                      scene4 = list(domain=list(x=c(0.55,1),y=c(0,0.45)),
                                                    xaxis=list(title = 'Standardised Time'),
                                                    yaxis =list(title = 'Standardised Time'),
                                                    zaxis=list(range = c(-0.02,0.065),title = "Env Cov-Phase"),
                                                    aspectmode='cube'))
fig_warpcov
################################################################################

### cross covariance of amplitude and warping functions

fig_xcov1 <-  plot_ly(showscale=FALSE, scene='scene') 
fig_xcov1 <- fig_xcov1 %>% add_surface(z = ~sample_cov[1:8,10:17], x = t, y = t)

fig_xcov2 <-  plot_ly(showscale=FALSE, scene='scene2') 
fig_xcov2 <- fig_xcov2 %>% add_surface(z = ~P_funComb_hat[1:8,10:17], x = t, y = t)

fig_xcov3 <-  plot_ly(showscale=FALSE, scene='scene3') 
fig_xcov3 <- fig_xcov3 %>% add_surface(z = ~CG_funComb_hat[1:8,10:17], x = t, y = t)
fig_xcov3 <- fig_xcov3 %>% add_surface(z = ~CG_funcomb_upper[1:8,10:17], x = ~t, y = ~t,
                                       opacity = 0.5, colorscale = list(c(0,1), c("grey", "lightgrey")))
fig_xcov3 <- fig_xcov3 %>% add_surface(z = ~CG_funcomb_lower[1:8,10:17], x = ~t, y = ~t,
                                       opacity = 0.5, colorscale = list(c(0,1), c("grey", "lightgrey")))

fig_xcov4 <-  plot_ly(showscale=FALSE, scene='scene4') 
fig_xcov4 <- fig_xcov4 %>% add_surface(z = ~CE_funComb_hat[1:8,10:17], x = t, y = t)
fig_xcov4 <- fig_xcov4 %>% add_surface(z = ~CE_funcomb_upper[1:8,10:17], x = ~t, y = ~t,
                                       opacity = 0.5, colorscale = list(c(0,1), c("grey", "lightgrey")))
fig_xcov4 <- fig_xcov4 %>% add_surface(z = ~CE_funcomb_lower[1:8,10:17], x = ~t, y = ~t,
                                       opacity = 0.5, colorscale = list(c(0,1), c("grey", "lightgrey")))


fig_xcov <- subplot(fig_xcov1, fig_xcov2, fig_xcov3, fig_xcov4)
fig_xcov <- fig_xcov %>% layout(scene = list(title = "",
                                             domain=list(x=c(0,0.45),y=c(0.55,1)),
                                             xaxis=list(title = 'Standardised Time'),
                                             yaxis =list(title = 'Standardised Time') , 
                                             zaxis=list(range = c(-0.02,0.065),title = "Sample XCov"),
                                             aspectmode='cube'),
                                scene2 = list(domain=list(x=c(0.55,1),y=c(0.55,1)),
                                              xaxis=list(title = 'Standardised Time'),
                                              yaxis =list(title = 'Standardised Time'),
                                              zaxis=list(range = c(-0.02,0.065),title = "P XCov"),
                                              aspectmode='cube'),
                                scene3 = list(domain=list(x=c(0,0.45),y=c(0,0.45)),
                                              xaxis=list(title = 'Standardised Time'),
                                              yaxis =list(title = 'Standardised Time'),
                                              zaxis=list(range = c(-0.02,0.065),title = "Gen XCov"),
                                              aspectmode='cube'),
                                scene4 = list(domain=list(x=c(0.55,1),y=c(0,0.45)),
                                              xaxis=list(title = 'Standardised Time'),
                                              yaxis =list(title = 'Standardised Time'),
                                              zaxis=list(range = c(-0.02,0.065),title = "Env XCov"),
                                              aspectmode='cube'))
fig_xcov
################################################################################

## cross covariance of the larval period length with the amplitude and warping functions

par(mfrow = c(2,2))
par(font.main = 1)
plot(c(0,1), c(-0.03, 0.02), type = "n", 
     xlab = 'Standardised Time', 
     ylab = "",
     xlim = c(0, 1), ylim = c(-0.03,0.02), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0,lwd = 2) 
axis(side = 2, at = seq(-0.03, 0.02, by = 0.01), pos = 0,lwd = 2) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
lines(t, sample_cov[1:8,9], type = "l", col = "red", lwd = 1.5)
mtext("Sample XCov: Amplitude with Larval Period", side = 3, adj = 0, line = 1, font = 2)

par(font.main = 1)
plot(c(0,1), c(-0.03, 0.02), type = "n", 
     xlab = 'Standardised Time', 
     ylab = "",
     xlim = c(0, 1), ylim = c(-0.03,0.02), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0,lwd = 2) 
axis(side = 2, at = seq(-0.03, 0.02, by = 0.01), pos = 0,lwd = 2) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
lines(t, P_funComb_hat[1:8,9], type = "l", col = "red")
mtext("P XCov: Amplitude with Larval Period", side = 3, adj = 0, line = 1, font = 2)

par(font.main = 1)
par(font.main = 1)
plot(c(0,1), c(-0.03, 0.02), type = "n", 
     xlab = 'Standardised Time', 
     ylab = "",
     xlim = c(0, 1), ylim = c(-0.03,0.02), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0,lwd = 2) 
axis(side = 2, at = seq(-0.03, 0.02, by = 0.01), pos = 0,lwd = 2) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
lines(t, CG_funComb_hat[1:8,9], type = "l", col = "red")
lines(t, CG_funcomb_upper[1:8,9], type = "l", lty = "dashed")
lines(t, CG_funcomb_lower[1:8,9], type = "l", lty = "dashed")
#lines(t, mean_CG_funcomb[1:8,9], type = "l", col = "black")
mtext("Gen XCov: Amplitude with Larval Period", side = 3, adj = 0, line = 1, font = 2)


par(font.main = 1)
plot(c(0,1), c(-0.03, 0.02), type = "n", 
     xlab = 'Standardised Time', 
     ylab = "",
     xlim = c(0, 1), ylim = c(-0.03,0.02), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0,lwd = 2) 
axis(side = 2, at = seq(-0.03, 0.02, by = 0.01), pos = 0,lwd = 2) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
lines(t, CE_funComb_hat[1:8,9], type = "l", col = "red")
lines(t, CE_funcomb_upper[1:8,9], type = "l", lty = "dashed")
lines(t, CE_funcomb_lower[1:8,9], type = "l", lty = "dashed")
#lines(t, mean_CE_funcomb[1:8,9], type = "l", col = "black")
mtext("Env XCov: Amplitude with Larval Period", side = 3, adj = 0, line = 1, font = 2)
################################################################################

par(font.main = 1)
plot(c(0,1), c(-0.03,0.02), type = "n", 
     xlab = 'Standardised Time', 
     ylab = "",
     xlim = c(0, 1), ylim = c(-0.03,0.02), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0,lwd = 2) 
axis(side = 2, at = seq(-0.03, 0.02, by = 0.01), pos = 0,lwd = 2) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
lines(t, sample_cov[10:17,9], type = "l", col = "red")
mtext("Sample XCov: Phase with Larval Period", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, v=0)

par(font.main = 1)
plot(c(0,1), c(-0.03,0.02), type = "n", 
     xlab = 'Standardised Time', 
     ylab = "",
     xlim = c(0, 1), ylim = c(-0.03,0.02), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0,lwd = 2) 
axis(side = 2, at = seq(-0.03, 0.02, by = 0.01), pos = 0,lwd = 2) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
lines(t, P_funComb_hat[10:17,9], type = "l", col = "red")
mtext("P XCov: Phase with Larval Period", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, v=0)

par(font.main = 1)
plot(c(0,1), c(-0.03,0.02), type = "n", 
     xlab = 'Standardised Time', 
     ylab = "",
     xlim = c(0, 1), ylim = c(-0.03,0.02), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0,lwd = 2) 
axis(side = 2, at = seq(-0.03, 0.02, by = 0.01), pos = 0,lwd = 2) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
lines(t, CG_funComb_hat[10:17,9], type = "l", col = "red")
lines(t, CG_funcomb_upper[10:17,9], type = "l", lty = "dashed")
lines(t, CG_funcomb_lower[10:17,9], type = "l", lty = "dashed")
#lines(t, mean_CG_funcomb[10:17,9], type = "l", col = "black")
mtext("Gen XCov: Phase with Larval Period", side = 3, adj = 0, line = 1, font = 2)

par(font.main = 1)
plot(c(0,1), c(-0.03,0.02), type = "n", 
     xlab = 'Standardised Time', 
     ylab = "",
     xlim = c(0, 1), ylim = c(-0.03,0.02), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0,lwd = 2) 
axis(side = 2, at = seq(-0.03, 0.02, by = 0.01), pos = 0,lwd = 2) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
lines(t, CE_funComb_hat[10:17,9], type = "l", col = "red")
lines(t, CE_funcomb_upper[10:17,9], type = "l", lty = "dashed")
lines(t, CE_funcomb_lower[10:17,9], type = "l", lty = "dashed")
#lines(t, mean_CE_funcomb[10:17,9], type = "l", col = "black")
mtext("Env XCov: Phase with Larval Period Length", side = 3, adj = 0, line = 1, font = 2)
################################################################################

## variance of Larval Period Length
par(mfrow = c(2,2))
par(font.main = 1)
plot(c(0,1), c(-0.01, 0.03), type = "n", 
     xlab = "", 
     ylab = "",
     xlim = c(0, 1), ylim = c(-0.01,0.03), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
segments(x0 = 0.5, y0 = 0, x1 = 0.5, y1 = sample_cov[9,9], col = "black", lwd = 2.5)  # Stick
points(0.5, sample_cov[9,9], pch = 16, col = "red", cex = 1.5)
mtext("Sample Var: Larval Period", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0,lwd = 2)
axis(side = 2, at = seq(-0.01, 0.03, by = 0.01), pos = 0, lwd = 2) 

par(font.main = 1)
plot(c(0,1), c(-0.01, 0.03), type = "n", 
     xlab = "", 
     ylab = "",
     xlim = c(0, 1), ylim = c(-0.01,0.03), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
segments(x0 = 0.5, y0 = 0, x1 = 0.5, y1 = P_funComb_hat[9,9], col = "black", lwd = 2.5)  # Stick
points(0.5, P_funComb_hat[9,9], pch = 16, col = "red", cex = 1.5)
mtext("P Var: Larval Period", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0,lwd = 2)
axis(side = 2, at = seq(-0.01, 0.03, by = 0.01), pos = 0, lwd = 2) 

par(font.main = 1)
plot(c(0,1), c(-0.01, 0.03), type = "n", 
     xlab = "", 
     ylab = "",
     xlim = c(0, 1), ylim = c(-0.01,0.03), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
abline(h = 0,lwd = 1.5)
segments(x0 = 0.5, y0 = 0, x1 = 0.5, y1 = CG_funcomb_upper[9,9], col = "black", lwd = 2.5)
segments(x0 = 0.45, y0 = 0, x1 = 0.55, y1 = 0, col = "black", lwd = 2.5) 
segments(x0 = 0.45, y0 = CG_funcomb_upper[9,9], x1 = 0.55, y1 = CG_funcomb_upper[9,9], col = "black", lwd = 2.5) # Stick
points(0.5, CG_funComb_hat[9,9], pch = 16, col = "red", cex = 1.5)
mtext("Gen Var: Larval Period", side = 3, adj = 0, line = 1, font = 2)
axis(side = 2, at = seq(-0.01, 0.03, by = 0.01), pos = 0, lwd = 2)
#legend("bottomright", legend= "Initial estimate", pch = 16, bty = "n", col = "red")

par(font.main = 1)
plot(c(0,1), c(-0.01, 0.03), type = "n", 
     xlab = "", 
     ylab = "",
     xlim = c(0, 1), ylim = c(-0.01,0.03), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
segments(x0 = 0.5, y0 = 0, x1 = 0.5, y1 = CE_funcomb_upper[9,9], col = "black", lwd = 2.5)
segments(x0 = 0.45, y0 = 0, x1 = 0.55, y1 = 0, col = "black", lwd = 2.5) 
segments(x0 = 0.45, y0 = CE_funcomb_upper[9,9], x1 = 0.55, y1 = CE_funcomb_upper[9,9], col = "black", lwd = 2.5) # Stick
points(0.5, CE_funComb_hat[9,9], pch = 16, col = "red", cex = 1.5)
mtext("Env Var: Larval Period", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0,lwd = 1.5)
axis(side = 2, at = seq(-0.01, 0.03, by = 0.01), pos = 0, lwd = 2)
#legend("bottomright", legend= "Initial estimate", pch = 16, bty = "n", col = "red")
################################################################################

## Predicting Evolutionary Response to Selection

### eigen decomposition of the genetic covariance
par(mfrow = c(1,3))
par(font.main = 1)
plot(c(0,1), c(-1,1), type = "n", 
     xlab = 'Standardised Time', 
     ylab = "",
     xlim = c(0, 1), ylim = c(-1,1), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0,lwd = 2) 
axis(side = 2, at = seq(-1, 1, by = 0.2), pos = 0,lwd = 2) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
lines(t, GEign_funcomb_hat[1:8], type = "l", col = "red", lwd = 1.5)
lines(t, GEign_funcomb_upper[1:8], lty = "dashed", col = "black")
lines(t, GEign_funcomb_lower[1:8], lty = "dashed", col = "black")
mtext("G: 1st Eigenfunction (Amplitude)", side = 3, adj = 0, line = 1, font = 2)

par(font.main = 1)
plot(c(0,1), c(-1,1), type = "n", 
     xlab = 'Standardised Time', 
     ylab = "",
     xlim = c(0, 1), ylim = c(-1,1), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0,lwd = 2) 
axis(side = 2, at = seq(-1, 1, by = 0.2), pos = 0,lwd = 2) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
lines(t, GEign_funcomb_hat[10:17], type = "l", col = "red", lwd = 1.5)
lines(t, GEign_funcomb_upper[10:17], lty = "dashed", col = "black")
lines(t, GEign_funcomb_lower[10:17], lty = "dashed", col = "black")
mtext("G: 1st Eigenfunction (Warping)", side = 3, adj = 0, line = 1, font = 2)

par(font.main = 1)
plot(c(0,1), c(-1,1), type = "n", 
     xlab = "", 
     ylab = "",
     xlim = c(0, 1), ylim = c(-1,1), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 2, at = seq(-1, 1, by = 0.2), pos = 0,lwd = 2) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
segments(x0 = 0.45, y0 = GEign_funcomb_upper[9], x1 = 0.55, y1 = GEign_funcomb_upper[9], col = "black", lwd = 2.5) # Stick
segments(x0 = 0.45, y0 = GEign_funcomb_lower[9], x1 = 0.55, y1 = GEign_funcomb_lower[9], col = "black", lwd = 2.5) # Stick
segments(x0 = 0.5, y0 = GEign_funcomb_lower[9], x1 = 0.5, y1 = GEign_funcomb_upper[9], col = "black", lwd = 2.5)
points(0.5, GEign_funcomb_hat[9], pch = 16, col = "red", cex = 1.5)
mtext("G: 1st Eigenfunction (Length)", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, lwd = 2)
################################################################################

## inter-generation change in the mean 
GEigen_eval <- eigen(CG_funComb_hat)$values[1]

dy1 <- GEigen_eval * GEign_funcomb_hat[1:8]

dgam1 <- GEigen_eval * GEign_funcomb_hat[10:17]

gam_mean <- rowMeans(gam)
h <- h_to_gam(gam)
h_mean <- rowMeans(h) # mean of the warping functions, identity map
Y_mean <- convert_to_basisfunctions(timefine, y_FR_mean, t)

## predicted evolutionary response to selection
y_new1 <- Y_mean + dy1

clrwarping_new1 <- gam_mean + dgam1

h_new1 <- h_to_gam(clrwarping_new1)


par(mfrow = c(1,2))
par(mar = c(5, 6, 4, 2) + 0.1)
plot(c(0,1), c(0, 3), type = "n", 
     xlab = 'Standardised Time', 
     ylab = expression(Log(Mass,~10^-5~g)),
     xlim = c(0, 1), ylim = c(0,3), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
lines(t, Y_mean, type = "l", lty = "solid", col = "black")
lines(t, y_new1, type = "l", lty = "dashed", col = "red")
legend("bottomright", legend= c("Mean: current generation", "Mean: next generation"), lty = c("solid", "dashed"), bty = "n", col = c("black","red"))
mtext("Response to Seletction (1st Eigenfunction, Mass)", side = 3, adj = 0, line = 1, font = 2)
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0, lwd = 2) 
axis(side = 2, at = seq(-1, 3, by = 0.5), pos = 0, lwd = 2)

plot(c(0,1), c(0, 1), type = "n", 
     xlab = 'Standardised Time', 
     ylab = "Warped Time",
     xlim = c(0, 1), ylim = c(0,1), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
lines(t, h_mean, type = "l", lty = "solid", col = "black")
lines(t, h_new1, type = "l", lty = "dashed", col = "red")
legend("bottomright", legend= c("Mean: current generation", "Mean: next generation"), lty = c("solid", "dashed"), bty = "n", col = c("black","red"))
mtext("Response to Seletction (1st Eigenfunction, Warping)", side = 3, adj = 0, line = 1, font = 2)
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0, lwd = 2) 
axis(side = 2, at = seq(0, 1, by = 0.1), pos = 0, lwd = 2)
################################################################################
