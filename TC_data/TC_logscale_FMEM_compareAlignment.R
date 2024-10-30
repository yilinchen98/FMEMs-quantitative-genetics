# We fit the TC growth data without alignment

## load data

setwd("D:/KCL_2023-2027_PhD/Genetics_FMEM_Project/FMEMs-quantitative-genetics/R_code/TC_data")
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
log_trait_list <- split(df$logtrait,id)

## plot oringinal growth curve
par(mar = c(5, 6, 4, 2) + 0.1)
plot(c(1, 25), c(0, 400), type = "n", 
     xlab = "Days", 
     ylab = expression(Mass~(10^-5~g)),
     xlim = c(0, 25), ylim = c(0, 400), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 25, by = 5), pos = 0) 
axis(side = 2, at = seq(0, 400, by = 100), pos = 0) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(age_list[[i]], trait_list[[i]], type = "l", col = i)
}
abline(h = 0, v=0)
#########################################################################################

## Data smoothing, set smoothing parameter <= 1e-4

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

## plot raw log growth curves vs plot smoothed log growth curves
par(mfrow = c(1,2))
par(mar = c(5, 6, 4, 2) + 0.1)
plot(c(0,1), c(0, 3), type = "n", 
     xlab = "Time", 
     ylab = expression(Log(Mass,~10^-5~g)),
     xlim = c(0, 1), ylim = c(0,3), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0) 
axis(side = 2, at = seq(0, 3, by = 0.5), pos = 0) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(age_list_new[[i]], log10(trait_list[[i]]), type = "l", col = i)
}
abline(h = 0, v=0)
mtext("Unsmoothed Log Growth Cruves", side = 3, adj = 0, line = 1, font = 2)

par(mar = c(5, 6, 4, 2) + 0.1)
plot(c(0,1), c(0, 3), type = "n", 
     xlab = "Time", 
     ylab = expression(Log(Mass,~10^-5~g)),
     xlim = c(0, 1), ylim = c(0,3), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0) 
axis(side = 2, at = seq(0, 3, by = 0.5), pos = 0) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(timefine, pred_logmass_fine[,i], type = "l", col = i)
}
abline(h = 0, v=0)
mtext("Smoothed Log Growth Cruves", side = 3, adj = 0, line = 1, font = 2)

#######################################################################################################

## Compare with registered curves
aligned_logmass_process <- time_warping(pred_logmass_fine, timefine)
aligned_logmass_curve <- aligned_logmass_process$fn
aligned_logmass_mean <- aligned_logmass_process$fmean
warping_logmass_funs <- aligned_logmass_process$warping_functions

## plot raw log growth curves vs plot smoothed registered log growth curves
par(mfrow = c(1,2))
par(mar = c(5, 6, 4, 2) + 0.1)
plot(c(0,1), c(0, 3), type = "n", 
     xlab = "Time", 
     ylab = expression(Log(Mass,~10^-5~g)),
     xlim = c(0, 1), ylim = c(0,3), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0) 
axis(side = 2, at = seq(0, 3, by = 0.5), pos = 0) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(age_list_new[[i]], log10(trait_list[[i]]), type = "l", col = i)
}
abline(h = 0, v=0)
mtext("Unsmoothed Log Growth Cruves", side = 3, adj = 0, line = 1, font = 2)

par(mar = c(5, 6, 4, 2) + 0.1)
plot(c(0,1), c(0, 3), type = "n", 
     xlab = "Time", 
     ylab = expression(Log(Mass,~10^-5~g)),
     xlim = c(0, 1), ylim = c(0,3), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0) 
axis(side = 2, at = seq(0, 3, by = 0.5), pos = 0) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(timefine, aligned_logmass_curve[,i], type = "l", col = i)
}
abline(h = 0, v=0)
mtext("Smoothed Registered Log Growth Cruves", side = 3, adj = 0, line = 1, font = 2)

####################################################################################

## Register unsmoothed curves

aligned_logtrait <- list()
gamma_logmass <- list()
for (i in 1:N){
  gamma_logmass_inter <- convert_to_basisfunctions(timefine, warping_logmass_funs[,i], age_list_new[[i]])
  gamma_logmass[[i]] <- gamma_logmass_inter
  aligned_logtrait[[i]] <- warp_f_gamma(log10(trait_list[[i]]), age_list_new[[i]], gamma_logmass_inter)
}

## plot raw log growth curves vs plot smoothed registered log growth curves
par(mfrow = c(1,2))
par(mar = c(5, 6, 4, 2) + 0.1)
plot(c(0,1), c(0, 3), type = "n", 
     xlab = "Time", 
     ylab = expression(Log(Mass,~10^-5~g)),
     xlim = c(0, 1), ylim = c(0,3), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0) 
axis(side = 2, at = seq(0, 3, by = 0.5), pos = 0) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(age_list_new[[i]], log10(trait_list[[i]]), type = "l", col = i)
}
abline(h = 0, v=0)
mtext("Unsmoothed Unregistered Log Growth Cruves", side = 3, adj = 0, line = 1, font = 2)

par(mar = c(5, 6, 4, 2) + 0.1)
plot(c(0,1), c(0, 3), type = "n", 
     xlab = "Time", 
     ylab = expression(Log(Mass,~10^-5~g)),
     xlim = c(0, 1), ylim = c(0,3), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0) 
axis(side = 2, at = seq(0, 3, by = 0.5), pos = 0) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(age_list_new[[i]], aligned_logtrait[[i]], type = "l", col = i)
}
abline(h = 0, v=0)
mtext("Unsmoothed Registered Log Growth Cruves", side = 3, adj = 0, line = 1, font = 2)

###########################################################################################

## Compare sample covariance for registered and unregistered curvevs

sam_orin <- GetCovSurface(log_trait_list, age_list_new, optns = list(nRegGrid = 100)) # smooth kernel: Gauss
sam_registered <- GetCovSurface(aligned_logtrait, age_list_new, optns = list(nRegGrid = 100))

cov_sam_orin <- sam_orin$cov
cov_sam_registered <- sam_registered$cov

fig1 <- plot_ly(showscale=FALSE, scene='scene1') 
fig1 <- fig1 %>% add_surface(z = ~cov_sam_orin, x = timefine, y = timefine)

fig2 <- plot_ly(showscale = FALSE, scene='scene2') 
fig2 <- fig2 %>% add_surface(z = ~cov_sam_registered, x= timefine, y = timefine)

fig_RR <- subplot(fig1, fig2) 
fig_RR <- fig_RR %>% layout(scene = list(title = "",
                                           domain=list(x=c(0,0.45),y=c(0.25,1)),
                                           xaxis=list(title = "Time"),
                                           yaxis =list(title = "Time") , 
                                           zaxis=list(range=c(-0.001,0.1),title = "Oringial (Log(Mass, 10<sup>-5</sup> g/t))"),
                                           aspectmode='cube'),
                              scene2 = list(domain=list(x=c(0.50,0.95),y=c(0.25,1)),
                                            xaxis=list(title = "Time"),
                                            yaxis =list(title = "Time"),
                                            zaxis=list(range=c(-0.001,0.1),title = "Aligned Raw (Log(Mass, 10<sup>-5</sup> g/t))"),
                                            aspectmode='cube'))

fig_RR


#########################################################################################


## FPCA

fpcaobj_lg <- prcomp(t(pred_logmass_fine), rank. = 8)
pcs_lg <- fpcaobj_lg$rotation
par(mfrow = c(1,1))
pcl <- fviz_eig(fpcaobj_lg, 
                addlabels = TRUE, 
                main="", ncp = 10)
pcl

par(mfrow = c(2, 1))
par(font.main = 1)
plot(c(0,1), c(0, 1), type = "n", 
     xlab = "Time", 
     ylab = "",
     xlim = c(0, 1), ylim = c(-0.2,0.2), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0) 
axis(side = 2, at = seq(-0.2, 0.25, by = 0.05), pos = 0) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(timefine, pcs_lg[,1], type = "l", col = i)
}
mtext("Principal Component 1 for Non-aligned Curves", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, v=0)

par(font.main = 1)
plot(c(0,1), c(0, 1), type = "n", 
     xlab = "Time", 
     ylab = "",
     xlim = c(0, 1), ylim = c(-0.2,0.25), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0) 
axis(side = 2, at = seq(-0.2, 0.25, by = 0.05), pos = 0) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(timefine, pcs_lg[,2], type = "l", col = i)
}
mtext("Principal Component 2 for Non-aligned Curves", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, v=0)

par(font.main = 1)
plot(c(0,1), c(0, 1), type = "n", 
     xlab = "Time", 
     ylab = "",
     xlim = c(0, 1), ylim = c(-0.2,0.25), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0) 
axis(side = 2, at = seq(-0.2, 0.25, by = 0.05), pos = 0) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(timefine, pcs_lg[,3], type = "l", col = i)
}
mtext("Principal Component 3 for Non-aligned Curves", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, v=0)

par(font.main = 1)
plot(c(0,1), c(0, 1), type = "n", 
     xlab = "Time", 
     ylab = "",
     xlim = c(0, 1), ylim = c(-0.2,0.25), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0) 
axis(side = 2, at = seq(-0.2, 0.25, by = 0.05), pos = 0) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(timefine, pcs_lg[,4], type = "l", col = i)
}
mtext("Principal Component 4 for Non-aligned Curves", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, v=0)
###############################################################################################

## Model fitting: fit original data

### prepare basis functions
phi_lg_list <- list() 
# create an empty list which stores eigenfunctions for 873 subjects
# evaluated at the original time points.

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

isSingular(fflg) # original tol = 1e-4 FALSE

theta <- fflg@theta
theta
## Model results

## fixed effect
smoothed_mean <- rowMeans(pred_logmass_fine)
betalg <- fixef(fflg)
feflg <- pcs_lg %*% betalg

df_FElg <- data.frame(
  Time = rep(timefine, 2),
  Value = c(smoothed_mean, feflg),
  FE = rep(c("mean", "estimated fixed effect"), each = length(timefine))
)

pfelg <- ggplot(df_FElg, aes(x = Time, y = Value, linetype = FE, color = FE)) +
  geom_line() +  
  labs(x = "Time", y = expression(Log(Mass,~10^-5~g)), title = "") +
  scale_color_manual(values = c("mean" = "red", "estimated fixed effect" = "black")) +  # Assign colors to lines
  scale_linetype_manual(values = c("mean" = "solid", "estimated fixed effect" = "dashed")) +  # Assign line types
  theme_minimal() +
  theme(legend.position = c(0.8,0.2))
pfelg
#ggsave("FElogmass.png", plot = pfe, width = 6, height = 3, dpi = 300)

## Random effects
vclg <- VarCorr(fflg)
CGlg <- vclg[["df_lg.id"]] # genetic covariance
CElg <- vclg[["df_lg.id.1"]] # environmental covariance

CG_funlg <- pcs_lg[,1:4] %*% CGlg %*% t(pcs_lg[,1:4]) # estimated gen cov function
CE_funlg <- pcs_lg[,1:4] %*% CElg %*% t(pcs_lg[,1:4]) # estimated env cov function

fig3 <- plot_ly(showscale=FALSE, scene='scene1') 
fig3 <- fig3 %>% add_surface(z = ~CG_funlg, x = timefine, y = timefine)

fig4 <- plot_ly(showscale = FALSE, scene='scene2') 
fig4 <- fig4 %>% add_surface(z = ~CE_funlg, x= timefine, y = timefine)

fig_RR2 <- subplot(fig3, fig4) 
fig_RR2 <- fig_RR2 %>% layout(scene = list(title = "",
                                           domain=list(x=c(0,0.45),y=c(0.25,1)),
                                           xaxis=list(title = "Time"),
                                           yaxis =list(title = "Time") , 
                                           zaxis=list(range=c(-0.001,0.1),title = "Genetic Covariance (Log(Mass, 10<sup>-5</sup> g/t))"),
                                           aspectmode='cube'),
                              scene2 = list(domain=list(x=c(0.50,0.95),y=c(0.25,1)),
                                            xaxis=list(title = "Time"),
                                            yaxis =list(title = "Time"),
                                            zaxis=list(range=c(-0.001,0.01),title = "Environmental Covariance (Log(Mass, 10<sup>-5</sup> g/t))"),
                                            aspectmode='cube'))

fig_RR2

P_funlg <- CG_funlg + CE_funlg
fig5 <- plot_ly(showscale = FALSE, scene='scene2') 
fig5 <- fig5 %>% add_surface(z = ~P_funlg, x= timefine, y = timefine)

fig_RR3 <- subplot(fig5, fig3) 
fig_RR3 <- fig_RR3 %>% layout(title = "",
                              scene2 = list(domain=list(x=c(0,0.45),y=c(0.25,1)),
                                            xaxis=list(title = "Time"),
                                            yaxis =list(title = "Time") , 
                                            zaxis=list(range=c(-0.001,0.1),title = "Phenotypic Covariance (Log(Mass, 10<sup>-5</sup> g/t))"),
                                            aspectmode='cube'),
                              scene = list(domain=list(x=c(0.50,0.95),y=c(0.25,1)),
                                           xaxis=list(title = "Time"),
                                           yaxis =list(title = "Time"),
                                           zaxis=list(range=c(-0.001,0.1),title = "Genetic Covariance (Log(Mass, 10<sup>-5</sup> g/t))"),
                                           aspectmode='cube'))

fig_RR3
################################################################################################################################

## FPCA of alinged smoothed curves

fpcaobj_logmass <- prcomp(x=t(aligned_logmass_curve), retx = TRUE, center = TRUE, rank. = 3)
pcs_logmass <- fpcaobj_logmass$rotation # eigen vectors

## Model Fitting (aligned raw data at Log scale)

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
colnames(phi_logmass) <- c("phi1", "phi2", "phi3")

## Reform dataframe
df_logmass <- data.frame(id = id, trait = unsplit(aligned_logtrait,id), phi_logmass)

fmmFormL <- trait ~ -1 + df_logmass$phi1 + df_logmass$phi2 + df_logmass$phi3 + 
  (-1 + df_logmass$phi1 + df_logmass$phi2 | df_logmass$id) + 
  (-1 + df_logmass$phi1 + df_logmass$phi2 | df_logmass$id) 

system.time(
  ffL <- fit_genetic_fmm(fmmFormL, df_logmass, A, 2)
) # user   system  elapsed

summary(ffL)

isSingular(ffL) # original tol = 1e-4

## Model results
### fixed effect
betaL <- fixef(ffL)
fefL <- pcs_logmass %*% betaL

df_FE <- data.frame(
  Time = rep(timefine, 2),
  Value = c(aligned_logmass_mean, fefL),
  FE = rep(c("mean", "estimated fixed effect"), each = length(timefine))
)

pfe <- ggplot(df_FE, aes(x = Time, y = Value, linetype = FE, color = FE)) +
  geom_line() +  
  labs(x = "Time", y = expression(Log(Mass,~10^-5~g)), title = "") +
  scale_color_manual(values = c("mean" = "red", "estimated fixed effect" = "black")) +  # Assign colors to lines
  scale_linetype_manual(values = c("mean" = "solid", "estimated fixed effect" = "dashed")) +  # Assign line types
  theme_minimal() +
  theme(legend.position = c(0.8,0.2))
pfe

## Random effects
vcL <- VarCorr(ffL)
CGL <- vcL[["df_logmass.id"]] # genetic covariance
CEL <- vcL[["df_logmass.id.1"]] # environmental covariance

CG_funL <- pcs_logmass[,1:2] %*% CGL %*% t(pcs_logmass[,1:2]) # estimated gen cov function
CE_funL <- pcs_logmass[,1:2] %*% CEL %*% t(pcs_logmass[,1:2]) # estimated env cov function
P_funL <- CG_funL + CE_funL
## 

fig_rawP <-  plot_ly(showscale=FALSE, scene='scene1') 
fig_rawP <- fig_rawP %>% add_surface(z = ~P_funlg, x= timefine, y = timefine)

fig_arawP <- plot_ly(showscale=FALSE, scene='scene2') 
fig_arawP <- fig_arawP %>% add_surface(z = ~P_funL, x= timefine, y = timefine)

fig_P <- subplot(fig_rawP, fig_arawP) 
fig_P <- fig_P %>% layout(scene = list(title = "",
                                           domain=list(x=c(0,0.45),y=c(0.25,1)),
                                           xaxis=list(title = "Time"),
                                           yaxis =list(title = "Time") , 
                                           zaxis=list(range=c(-0.001,0.1),title = "Phenotypic Original (Log(Mass, 10<sup>-5</sup> g/t))"),
                                           aspectmode='cube'),
                              scene2 = list(domain=list(x=c(0.50,0.95),y=c(0.25,1)),
                                            xaxis=list(title = "Time"),
                                            yaxis =list(title = "Time"),
                                            zaxis=list(range=c(-0.001,0.1),title = "Phenotypic Aligned (Log(Mass, 10<sup>-5</sup> g/t))"),
                                            aspectmode='cube'))
fig_P
