# Fit the TC dataset with PCs as basis functions

## Load Data
TRFUN25PUP4 = read.delim("TRFUN25PUP4.DAT",header = FALSE)
names(TRFUN25PUP4)<-c("id","sire","dam","trait","x")
df <- data.frame(TRFUN25PUP4)

## calculate genetic relationship matrix
id <- df$id
FirstUniqueIdPos <- which(duplicated(id) == FALSE)
pos = id[FirstUniqueIdPos] # extract ids for all subjects
sire_id = df$sire[FirstUniqueIdPos] # extract ids for sire
dam_id = df$dam[FirstUniqueIdPos] # extract ids for dam

pede <- editPed(sire = sire_id, dam = dam_id, label = pos)
ped<- with(pede, pedigree(label=label, sire=sire, dam=dam))
A <- getA(ped)[163:1035,163:1035]

## Rescale time interval to [0,1]
N = length(FirstUniqueIdPos) # N = 873 subjects
n = length(id) # n = 6860 observations
age_list <- split(df$x,id)
trait_list <- split(df$trait,id)

### x = (x -min(x))/(max(x) - min(x))
age_list_new <- list()
for (i in 1:N){
  age_list_new[[i]] = (age_list[[i]]-min(age_list[[i]]))/(max(age_list[[i]])-min(age_list[[i]]))
}

df$x_rescaled <- unsplit(age_list_new,id)
df$logtrait <- log10(df$trait)
id <- df$id

## plot growth curve
par(mfrow = c(1,1))
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
########################################################################################

## Smooth data, choose smoothing parameter <= 1e-4

timefine <- seq(0,1,length=100)
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

## plot raw growth curves vs plot smoothed growth curve
par(mfrow = c(1,2))
par(mar = c(5, 6, 4, 2) + 0.1)
plot(c(0,1), c(0, 3), type = "n", 
     xlab = "Standardised Time", 
     ylab = expression(Mass~(10^-5~g)),
     xlim = c(0, 1), ylim = c(0,400), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
for (i in 1:N){
  lines(age_list_new[[i]], trait_list[[i]], type = "l", col = i)
}
mtext("Raw Growth Curves", side = 3, adj = 0, line = 1, font = 2)
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0, lwd = 2) 
axis(side = 2, at = seq(0, 400, by = 100), pos = 0, lwd = 2) 

par(mar = c(5, 6, 4, 2) + 0.1)
plot(c(0,1), c(0, 3), type = "n", 
     xlab = "Standardised Time", 
     ylab = expression(Mass~(10^-5~g)),
     xlim = c(0, 1), ylim = c(0,400), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
for (i in 1:N){
  lines(timefine, pred_mass_fine[,i], type = "l", col = i)
}
mtext("Smoothed Growth Curves", side = 3, adj = 0, line = 1, font = 2)
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0, lwd = 2) 
axis(side = 2, at = seq(0, 400, by = 100), pos = 0, lwd = 2) 

## Align Growth Curves

aligned_mass_process <- time_warping(pred_mass_fine, timefine)
aligned_mass_curve <- aligned_mass_process$fn
aligned_mass_mean <- aligned_mass_process$fmean
warping_mass_funs <- aligned_mass_process$warping_functions

par(mfrow = c(1, 2))
par(mar = c(5, 6, 4, 2) + 0.1, font.main = 1)
plot(c(0,1), c(0, 1), type = "n", 
     xlab = "Standardised Time", 
     ylab = "",
     xlim = c(0, 1), ylim = c(0,1), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
for (i in 1:N){
  lines(timefine, warping_mass_funs[,i], type = "l", col = i)
}
mtext("Warping Functions", side = 3, adj = 0, line = 1, font = 2)
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0, lwd = 2) 
axis(side = 2, at = seq(0, 1, by = 0.1), pos = 0, lwd = 2) 

plot(c(0,1), c(0, 3), type = "n", 
     xlab = "Standardised Time", 
     ylab = expression(Mass~(10^-5~g)),
     xlim = c(0, 1), ylim = c(0,400), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
for (i in 1:N){
  lines(timefine, aligned_mass_curve[,i], type = "l", col = i)
}
mtext("Aligned Mass Curves", side = 3, adj = 0, line = 1, font = 2)
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0, lwd = 2) 
axis(side = 2, at = seq(0, 400, by = 100), pos = 0, lwd = 2)

## FPCA

fpcaobj_mass <- prcomp(x=t(aligned_mass_curve), retx = TRUE, center = TRUE, rank. = 2)
pcs_mass <- fpcaobj_mass$rotation # eigen vectors

par(mfrow = c(1,1))
fviz_eig(fpcaobj_mass, 
         addlabels = TRUE, 
         main="", ncp = 6)

par(mfrow = c(1, 2))
par(font.main = 1)
plot(c(0,1), c(0, 1), type = "n", 
     xlab = "Standardised Time", 
     ylab = "",
     xlim = c(0, 1), ylim = c(-0.2,0.2), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
for (i in 1:N){
  lines(timefine, pcs_mass[,1], type = "l", col = i)
}
mtext("Principal Component 1", side = 3, adj = 0, line = 1, font = 2)
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0, lwd = 2) 
axis(side = 2, at = seq(-0.2, 0.2, by = 0.05), pos = 0, lwd = 2) 

par(font.main = 1)
plot(c(0,1), c(0, 1), type = "n", 
     xlab = "Standardised Time", 
     ylab = "",
     xlim = c(0, 1), ylim = c(-0.2,0.2), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
for (i in 1:N){
  lines(timefine, pcs_mass[,2], type = "l", col = i)
}
mtext("Principal Component 2", side = 3, adj = 0, line = 1, font = 2)
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0, lwd = 2) 
axis(side = 2, at = seq(-0.2, 0.2, by = 0.05), pos = 0, lwd = 2) 


## Model Fitting (Original scale)

### Align raw data
aligned_trait <- list()
gamma_mass <- list()
for (i in 1:N){
  gamma_mass_inter <- convert_to_basisfunctions(timefine, warping_mass_funs[,i], age_list_new[[i]])
  gamma_mass[[i]] <- gamma_mass_inter
  aligned_trait[[i]] <- warp_f_gamma(trait_list[[i]], age_list_new[[i]], gamma_mass_inter)
}

par(mfrow = c(1, 2))
par(mar = c(5, 6, 4, 2) + 0.1, font.main = 1)
plot(c(0,1), c(0, 1), type = "n", 
     xlab = "Standardised Time", 
     ylab = "",
     xlim = c(0, 1), ylim = c(0,1), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
for (i in 1:N){
  lines(age_list_new[[i]], gamma_mass[[i]], type = "l", col = i)
}
mtext("Warping Functions", side = 3, adj = 0, line = 1, font = 2)
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0, lwd = 2) 
axis(side = 2, at = seq(0, 1, by = 0.1), pos = 0, lwd = 2) 

plot(c(0,1), c(0, 400), type = "n", 
     xlab = "Standardised Time", 
     ylab = expression(Mass~(10^-5~g)),
     xlim = c(0, 1), ylim = c(0,400), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
for (i in 1:N){
  lines(age_list_new[[i]], aligned_trait[[i]], type = "l", col = i)
}
mtext("Aligned Raw Mass Curves", side = 3, adj = 0, line = 1, font = 2)
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0, lwd = 2) 
axis(side = 2, at = seq(0, 400, by = 100), pos = 0, lwd = 2) 

### Prepare model basis
phi_mass_list <- list() 
# create an empty list which stores eigenfunctions for 873 subjects
# evaluated at the original time points.

for (i in 1:N){
  phi_mass <- convert_to_basisfunctions(t = timefine, eigenvecs = pcs_mass,
                                           tout = age_list_new[[i]])
  phi_mass_list[[i]] <- phi_mass
}

phi_mass <- do.call(rbind,phi_mass_list)
colnames(phi_mass) <- c("phi1", "phi2")

## Reform dataframe
df_mass <- data.frame(id = id, trait = unsplit(aligned_trait,id), phi_mass)

fmmForm <- trait ~ -1 + df_mass$phi1 + df_mass$phi2 + 
  (-1 + df_mass$phi1 + df_mass$phi2 | df_mass$id) + 
  (-1 + df_mass$phi1 + df_mass$phi2 | df_mass$id) 

system.time(
  ff <- fit_genetic_fmm(fmmForm, df_mass, A, 2)
) # user   system  elapsed

summary(ff)

isSingular(ff) # TRUE singular fit, tol = 1e-4

## fixed effect
beta <- fixef(ff)
fef <- pcs_mass %*% beta

par(mfrow = c(1,1))
par(mar = c(5, 6, 4, 2) + 0.1)
plot(c(0,1), c(0, 400), type = "n", 
     xlab = "Standardised Time", 
     ylab =  expression(Mass~(10^-5~g)),
     xlim = c(0, 1), ylim = c(0,400), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
lines(timefine, fef, type = "l", col = "black")
lines(timefine, aligned_mass_mean, type = "l", col = "red")
mtext("Fixed Effect", side = 3, adj = 0, line = 1, font = 2)
legend("bottomright", legend= c("Estimated FE", "Mean"), lty = c("solid", "dashed"), bty = "n", col = c("black","red"))
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0, lwd = 2) 
axis(side = 2, at = seq(0, 400, by = 100), pos = 0, lwd = 2) 
###########################
vc <- VarCorr(ff)
CG <- vc[["df_mass.id"]] # genetic covariance
CE <- vc[["df_mass.id.1"]] # environmental covariance

CG_fun <- pcs_mass %*% CG %*% t(pcs_mass) # estimated gen cov function
CE_fun <- pcs_mass %*% CE %*% t(pcs_mass) # estimated env cov function

fig3 <- plot_ly(showscale=FALSE, scene='scene1') 
fig3 <- fig3 %>% add_surface(z = ~CG_fun, x = timefine, y = timefine)

fig4 <- plot_ly(showscale = FALSE, scene='scene2') 
fig4 <- fig4 %>% add_surface(z = ~CE_fun, x= timefine, y = timefine)

fig_RR2 <- subplot(fig3, fig4) 
fig_RR2 <- fig_RR2 %>% layout(scene = list(title = "",
                                           domain=list(x=c(0,0.45),y=c(0.25,1)),
                                           xaxis=list(title = "Standardised Time"),
                                           yaxis =list(title = "Standardised Time") , 
                                           zaxis=list(range=c(-10,600),title = "Genetic Covariance (Mass, 10<sup>-5</sup> g/t)"),
                                           aspectmode='cube'),
                              scene2 = list(domain=list(x=c(0.50,0.95),y=c(0.25,1)),
                                            xaxis=list(title = "Standardised Time"),
                                            yaxis =list(title = "Standardised Time"),
                                            zaxis=list(range=c(-10,600),title = "Environmental Covariance (Mass, 10<sup>-5</sup> g/t)"),
                                            aspectmode='cube'))

fig_RR2


P_fun <- CG_fun + CE_fun
fig5 <- plot_ly(showscale = FALSE, scene='scene2') 
fig5 <- fig5 %>% add_surface(z = ~P_fun, x= timefine, y = timefine)

fig_RR3 <- subplot(fig3, fig5) 
fig_RR3 <- fig_RR3 %>% layout(title = "",
                              scene = list(domain=list(x=c(0,0.45),y=c(0.25,1)),
                                           xaxis=list(title = "Standardised Time"),
                                           yaxis =list(title = "Standardised Time") , 
                                           zaxis=list(range=c(-10,1200),title = "Genetic Covariance (Mass, 10<sup>-5</sup> g/t)"),
                                           aspectmode='cube'),
                              scene2 = list(domain=list(x=c(0.50,0.95),y=c(0.25,1)),
                                            xaxis=list(title = "Standardised Time"),
                                            yaxis =list(title = "Standardised Time"),
                                            zaxis=list(range=c(-10,1200),title = "Phenotypic Covariance (Mass, 10<sup>-5</sup> g/t)"),
                                            aspectmode='cube'))

fig_RR3

####################################################################################
## Estimate sample covariance 

### Estimate sample covariance function based on aligned raw data
sample_cov_obj_orin<- GetCovSurface(aligned_trait, age_list_new, optns = list(nRegGrid = 100))
#Warning message:
#In CheckData(Ly, Lt) :
  #There is a time gap of at least 10% of the observed range across subjects

sample_cov_function_orin <- sample_cov_obj_orin$cov
timegrid_orin <- sample_cov_obj_orin$workGrid

### Estimate sample covariance function based on aligned smoothed data
Ly_orin <- list()
for (i in 1:N){
  Ly_orin[[i]] <- aligned_mass_curve[,i]
}
Lt_orin <- replicate(N, timefine, simplify = FALSE)

sample_cov_obj_smooth_orin<- GetCovSurface(Ly_orin, Lt_orin, optns = list(nRegGrid = 100))
sample_cov_function_smooth_orin <- sample_cov_obj_smooth_orin$cov

fig_sam1 <-  plot_ly(showscale=FALSE, scene='scene') 
fig_sam1 <- fig_sam1 %>% add_surface(z = ~sample_cov_function_orin, x = timegrid_orin, y = timegrid_orin)

fig_sam2 <-  plot_ly(showscale=FALSE, scene='scene2') 
fig_sam2 <- fig_sam2 %>% add_surface(z = ~sample_cov_function_smooth_orin, x = timefine, y = timefine)

fig_sam3 <- plot_ly(showscale=FALSE, scene='scene3') 
fig_sam3 <- fig_sam3 %>% add_surface(z = ~P_fun, x = timefine, y = timefine)

fig_sam <- subplot(fig_sam1, fig_sam2, fig_sam3)
fig_sam <- fig_sam %>% layout(scene = list(title = "",
                                           domain=list(x=c(0,0.32),y=c(0.25,1)),
                                           xaxis=list(title = "Standardised Time"),
                                           yaxis =list(title = "Standardised Time") , 
                                           zaxis=list(range=c(-10,1200),title = "Sample Covariance (Mass, 10<sup>-5</sup> g/t)"),
                                           aspectmode='cube'),
                              scene2 = list(domain=list(x=c(0.33,0.65),y=c(0.25,1)),
                                            xaxis=list(title = "Standardised Time"),
                                            yaxis =list(title = "Standardised Time"),
                                            zaxis=list(range = c(-10, 1200),title = "Smooth Sample Covariance (Mass, 10<sup>-5</sup> g/t)"),
                                            aspectmode='cube'),
                              scene3 = list(domain=list(x=c(0.66,1),y=c(0.25,1)),
                                            xaxis=list(title = "Standardised Time"),
                                            yaxis =list(title = "Standardised Time"),
                                            zaxis=list(range = c(-10,1200),title = "Phenotypic Covariance (Mass, 10<sup>-5</sup> g/t)"),
                                            aspectmode='cube'))
fig_sam

#########################################################################################################

## Eigenfunctions

eigenfun_CG1 <- eigen(CG_fun)$vectors[,1]
eigenfun_CG2 <- eigen(CG_fun)$vectors[,2]

eigenfun_CE1 <- eigen(CE_fun)$vectors[,1]
eigenfun_CE2 <- eigen(CE_fun)$vectors[,2]

par(mfrow = c(1,2))
plot(c(0,1), c(-0.25, 0.2), type = "n", 
     xlab = "Standardised Time", 
     ylab = "",
     xlim = c(0, 1), ylim = c(-0.25,0.2), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
lines(timefine, eigenfun_CG1, type = "l", lty = "solid")
lines(timefine, eigenfun_CG2, type = "l", lty = "dashed")
mtext("Genetic Eigenfunction (Original Scale)", side = 3, adj = 0, line = 1, font = 2)
legend("bottomright", legend= c("First EF", "Second EF"), lty = c("solid", "dashed"), bty = "n")
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0, lwd = 2) 
axis(side = 2, at = seq(-0.25, 0.2, by = 0.05), pos = 0, lwd = 2)

plot(c(0,1), c(-0.25, 0.2), type = "n", 
     xlab = "Standardised Time", 
     ylab = "",
     xlim = c(0, 1), ylim = c(-0.25,0.2), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
lines(timefine, eigenfun_CE1, type = "l", lty = "solid")
lines(timefine, eigenfun_CE2, type = "l", lty = "dashed")
mtext("Environmental Eigenfunction (Original Scale)", side = 3, adj = 0, line = 1, font = 2)
legend("bottomright", legend= c("First EF", "Second EF"), lty = c("solid", "dashed"), bty = "n")
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0, lwd = 2) 
axis(side = 2, at = seq(-0.25, 0.2, by = 0.05), pos = 0, lwd = 2) 
###################################################################################
