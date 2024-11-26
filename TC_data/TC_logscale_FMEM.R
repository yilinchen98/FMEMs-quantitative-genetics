# Fit the TC dataset (log10-transformed) with PCs as basis functions

## Load Data
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

## Rescale time interval to [0,1]
### x = (x -min(x))/(max(x) - min(x))
age_list_new <- list()
for (i in 1:N){
  age_list_new[[i]] = (age_list[[i]]-min(age_list[[i]]))/(max(age_list[[i]])-min(age_list[[i]]))
}

df$x_rescaled <- unsplit(age_list_new,id)


## Plot growth curves
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
mtext("Growth Curves", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, v=0)

## Plot log growth curves
par(mar = c(5, 6, 4, 2) + 0.1)
plot(c(1, 25), c(0, 3), type = "n", 
     xlab = "Days", 
     ylab = expression(Log(Mass,~10^-5~g)),
     xlim = c(0, 25), ylim = c(0, 3), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 25, by = 5), pos = 0) 
axis(side = 2, at = seq(0, 3, by = 0.5), pos = 0) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(age_list[[i]], log_trait_list[[i]], type = "l", col = i)
}
mtext("Log Growth Curves", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, v=0)
########################################################################################

## Problem with smoothing
timefine <- seq(0,1,length=100) # dense time grid
demo <- c(16,534,5,46)
mass_demo <- matrix(0,100,4)
for (i in 1:4){
  ss_demo <-  smooth.spline(age_list_new[[demo[i]]], trait_list[[demo[i]]], cv=FALSE, all.knots=TRUE)
  mass_demo[,i] <- predict(ss_demo, timefine)$y
}

## Smoothe data, choose smoothing parameter <= 1e-4

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

#############################################################################################################
# Original data points
original_df <- data.frame(
  Time = unlist(lapply(demo, function(i) age_list_new[[i]])),
  Mass = unlist(lapply(demo, function(i) trait_list[[i]])),
  Subject = factor(rep(demo, times = sapply(demo, function(i) length(age_list_new[[i]]))))
)

# Smoothed mass
smoothed_df <- data.frame(
  Time = rep(timefine, 4),
  Smoothed_Mass = c(mass_demo[,1], mass_demo[,2], mass_demo[,3], mass_demo[,4]),
  Subject = factor(rep(demo, each = length(timefine)))
)

# Predicted mass
predicted_df <- data.frame(
  Time = rep(timefine, 4),
  Predicted_Mass = c(pred_mass_fine[,demo[1]], pred_mass_fine[,demo[2]], pred_mass_fine[,demo[3]], pred_mass_fine[,demo[4]]),
  Subject = factor(rep(demo, each = length(timefine)))
)

# Create a function to generate the plot for each subject
create_subject_plot <- function(subject_index) {
  subject <- demo[subject_index]
  
  # Original data plot
  p_original <- ggplot(original_df[original_df$Subject == subject, ], aes(x = Time, y = Mass)) +
    geom_point(color = "black") +
    labs(x = "Time", y =expression(Mass~(10^-5~g)), title = paste("Subject", subject)) +
    theme_minimal()
  
  # Smoothed mass plot
  p_smoothed <- ggplot(smoothed_df[smoothed_df$Subject == subject, ], aes(x = Time, y = Smoothed_Mass)) +
    geom_line(color = "blue") +
    theme_minimal()
  
  # Predicted mass plot
  p_predicted <- ggplot(predicted_df[predicted_df$Subject == subject, ], aes(x = Time, y = Predicted_Mass)) +
    geom_line(color = "red") +
    theme_minimal()
  
  # Combine all plots 
  p <- p_original +
    geom_line(data = smoothed_df[smoothed_df$Subject == subject, ], aes(y = Smoothed_Mass), color = "blue") +
    geom_line(data = predicted_df[predicted_df$Subject == subject, ], aes(y = Predicted_Mass), color = "red")
  
  return(p)
}

plots <- lapply(1:4, create_subject_plot)

combined_plot <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + 
  plot_layout(ncol = 2, nrow = 2)
print(combined_plot)

##############################################################################################################

## Example of two growth curves with data ponits
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
lines(age_list_new[[1]],log_trait_list[[1]], type = "l")
lines(age_list_new[[1]],log_trait_list[[1]], type = "p")
lines(age_list_new[[25]],log_trait_list[[25]], type = "l")
lines(age_list_new[[25]],log_trait_list[[25]], type = "p")
abline(h = 0, v=0)
                       
## plot raw log growth curves vs plot smoothed log growth curves
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

#################################################################################################

## Align Log Growth Curves

aligned_logmass_process <- time_warping(pred_logmass_fine, timefine)
aligned_logmass_curve <- aligned_logmass_process$fn
aligned_logmass_mean <- aligned_logmass_process$fmean
warping_logmass_funs <- aligned_logmass_process$warping_functions

par(mfrow = c(1, 2))
par(mar = c(5, 6, 4, 2) + 0.1, font.main = 1)
plot(c(0,1), c(0, 1), type = "n", 
     xlab = "Time", 
     ylab = "Warping Functions",
     xlim = c(0, 1), ylim = c(0,1), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0) 
axis(side = 2, at = seq(0, 1, by = 0.1), pos = 0) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(timefine, warping_logmass_funs[,i], type = "l", col = i)
}
mtext("Warping Functions", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, v=0)

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
mtext("Aligned Log Mass Curves", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, v=0)

############################################################################################

## FPCA

fpcaobj_logmass <- prcomp(x=t(aligned_logmass_curve), retx = TRUE, center = TRUE, rank. = 4)
pcs_logmass <- fpcaobj_logmass$rotation # eigen vectors

par(mfrow = c(1,1))
pl <- fviz_eig(fpcaobj_logmass, 
         addlabels = TRUE, 
         main="", ncp = 6)
pl
#ggsave("pclogmass.png", plot = pl, width = 6, height = 4, dpi = 300)

par(mfrow = c(2, 1))
par(font.main = 1)
plot(c(0,1), c(0, 1), type = "n", 
     xlab = "Time", 
     ylab = "",
     xlim = c(0, 1), ylim = c(-0.1,0.2), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0) 
axis(side = 2, at = seq(-0.1, 0.2, by = 0.05), pos = 0) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(timefine, pcs_logmass[,1], type = "l", col = i)
}
mtext("Principal Component 1 of Aligned Curves", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, v=0)

par(font.main = 1)
plot(c(0,1), c(0, 1), type = "n", 
     xlab = "Time", 
     ylab = "",
     xlim = c(0, 1), ylim = c(-0.1,0.2), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0) 
axis(side = 2, at = seq(-0.1, 0.2, by = 0.05), pos = 0) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(timefine, pcs_logmass[,2], type = "l", col = i)
}
mtext("Principal Component 2 of Aligned Curves", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, v=0)

par(font.main = 1)
plot(c(0,1), c(0, 1), type = "n", 
     xlab = "Time", 
     ylab = "",
     xlim = c(0, 1), ylim = c(-0.3,0.15), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0) 
axis(side = 2, at = seq(-0.3, 0.15, by = 0.05), pos = 0) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(timefine, pcs_logmass[,3], type = "l", col = i)
}
mtext("Principal Component 3 of Aligned Curves", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, v=0)

par(font.main = 1)
plot(c(0,1), c(0, 1), type = "n", 
     xlab = "Time", 
     ylab = "",
     xlim = c(0, 1), ylim = c(-0.2,0.4), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0) 
axis(side = 2, at = seq(-0.2, 0.4, by = 0.05), pos = 0) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(timefine, pcs_logmass[,4], type = "l", col = i)
}
mtext("Principal Component 4 of Aligned Curves", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, v=0)
################################################################################################

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
     xlab = "Time", 
     ylab = "Warping Functions",
     xlim = c(0, 1), ylim = c(0,1), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0) 
axis(side = 2, at = seq(0, 1, by = 0.1), pos = 0) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(age_list_new[[i]], gamma_logmass[[i]], type = "l", col = i)
}
mtext("Warping Functions", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, v=0)

plot(c(0,1), c(0, 3), type = "n", 
     xlab = "Time", 
     ylab = expression(Mass~(10^-5~g)),
     xlim = c(0, 1), ylim = c(0,3), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0) 
axis(side = 2, at = seq(0, 3, by = 0.5), pos = 0) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(age_list_new[[i]], aligned_logtrait[[i]], type = "l", col = i)
}
mtext("Aligned Raw Log Mass Curves", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, v=0)

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
colnames(phi_logmass) <- c("phi1", "phi2", "phi3","phi4")

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

###############################################################################

## Model Results

## fixed effect
betaL <- fixef(ffL)
fefL <- pcs_logmass[,1:3] %*% betaL

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
#ggsave("FElogmass.png", plot = pfe, width = 6, height = 3, dpi = 300)

## Random effects
vcL <- VarCorr(ffL)
CGL <- vcL[["df_logmass.id"]] # genetic covariance
CEL <- vcL[["df_logmass.id.1"]] # environmental covariance

CG_funL <- pcs_logmass[,1:2] %*% CGL %*% t(pcs_logmass[,1:2]) # estimated gen cov function
CE_funL <- pcs_logmass[,1:2] %*% CEL %*% t(pcs_logmass[,1:2]) # estimated env cov function

fig3 <- plot_ly(showscale=FALSE, scene='scene1') 
fig3 <- fig3 %>% add_surface(z = ~CG_funL, x = timefine, y = timefine)

fig4 <- plot_ly(showscale = FALSE, scene='scene2') 
fig4 <- fig4 %>% add_surface(z = ~CE_funL, x= timefine, y = timefine)

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

P_funL <- CG_funL + CE_funL
fig5 <- plot_ly(showscale = FALSE, scene='scene2') 
fig5 <- fig5 %>% add_surface(z = ~P_funL, x= timefine, y = timefine)

fig_RR3 <- subplot(fig3, fig5) 
fig_RR3 <- fig_RR3 %>% layout(title = "",
                              scene = list(domain=list(x=c(0,0.45),y=c(0.25,1)),
                                           xaxis=list(title = "Time"),
                                           yaxis =list(title = "Time") , 
                                           zaxis=list(range=c(-0.001,0.1),title = "Genetic Covariance (Log(Mass, 10<sup>-5</sup> g/t))"),
                                           aspectmode='cube'),
                              scene2 = list(domain=list(x=c(0.50,0.95),y=c(0.25,1)),
                                            xaxis=list(title = "Time"),
                                            yaxis =list(title = "Time"),
                                            zaxis=list(range=c(-0.001,0.1),title = "Phenotypic Covariance (Log(Mass, 10<sup>-5</sup> g/t))"),
                                            aspectmode='cube'))

fig_RR3

####################################################################################


## Estimate sample covariance 

### Estimate sample covariance function based on aligned raw data
sample_cov_obj<- GetCovSurface(aligned_logtrait, age_list_new)
#Warning message:
#In CheckData(Ly, Lt) :
  #There is a time gap of at least 10% of the observed range across subjects

sample_cov_function <- sample_cov_obj$cov
timegrid <- sample_cov_obj$workGrid

### Estimate sample covariance function based on aligned smoothed data
Ly <- list()
for (i in 1:N){
  Ly[[i]] <- aligned_logmass_curve[,i]
}
Lt <- replicate(N, timefine, simplify = FALSE)
sample_cov_obj_smooth<- GetCovSurface(Ly, Lt)
sample_cov_function_smooth <- sample_cov_obj_smooth$cov

fig_sam1 <-  plot_ly(showscale=FALSE, scene='scene') 
fig_sam1 <- fig_sam1 %>% add_surface(z = ~sample_cov_function, x = timegrid, y = timegrid)

fig_sam2 <-  plot_ly(showscale=FALSE, scene='scene2') 
fig_sam2 <- fig_sam2 %>% add_surface(z = ~sample_cov_function_smooth, x = timefine, y = timefine)

fig_sam3 <- plot_ly(showscale=FALSE, scene='scene3') 
fig_sam3 <- fig_sam3 %>% add_surface(z = ~P_funL, x = timefine, y = timefine)

fig_sam <- subplot(fig_sam1, fig_sam2, fig_sam3)
fig_sam <- fig_sam %>% layout(scene = list(title = "",
                                           domain=list(x=c(0,0.32),y=c(0.25,1)),
                                           xaxis=list(title = "Time"),
                                           yaxis =list(title = "Time") , 
                                           zaxis=list(range=c(-0.001,0.1),title = "Sample Covariance (Log(Mass, 10<sup>-5</sup> g/t))"),
                                           aspectmode='cube'),
                              scene2 = list(domain=list(x=c(0.33,0.65),y=c(0.25,1)),
                                            xaxis=list(title = "Time"),
                                            yaxis =list(title = "Time"),
                                            zaxis=list(range = c(-0.001, 0.1),title = "Smooth Sample Covariance (Log(Mass, 10<sup>-5</sup> g/t))"),
                                            aspectmode='cube'),
                              scene3 = list(domain=list(x=c(0.66,1),y=c(0.25,1)),
                                            xaxis=list(title = "Time"),
                                            yaxis =list(title = "Time"),
                                            zaxis=list(range = c(-0.001, 0.1),title = "Phenotypic Covariance (Log(Mass, 10<sup>-5</sup> g/t))"),
                                            aspectmode='cube'))
fig_sam

#########################################################################################################

## Eigenfunctions

eigenfunL_CG1 <- eigen(CG_funL)$vectors[,1]
eigenfunL_CG2 <- eigen(CG_funL)$vectors[,2]

eigenfunL_CE1 <- eigen(CE_funL)$vectors[,1]
eigenfunL_CE2 <- eigen(CE_funL)$vectors[,2]

par(mfrow = c(1,2))
plot(c(0,1), c(-0.25, 0.2), type = "n", 
     xlab = "Time", 
     ylab = "",
     xlim = c(0, 1), ylim = c(-0.25,0.2), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = -0.25) 
axis(side = 2, at = seq(-0.25, 0.2, by = 0.05), pos = 0)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
lines(timefine, eigenfunL_CG1, type = "l", lty = "solid")
lines(timefine, eigenfunL_CG2, type = "l", lty = "dashed")
mtext("Genetic Eigenfunction (Log Scale)", side = 3, adj = 0, line = 1, font = 2)
legend("bottomright", legend= c("First EF", "Second EF"), lty = c("solid", "dashed"), bty = "n")
abline(v=0)
abline(h=0, lty = "dashed")

plot(c(0,1), c(-0.25, 0.2), type = "n", 
     xlab = "Time", 
     ylab = "",
     xlim = c(0, 1), ylim = c(-0.25,0.2), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = -0.25) 
axis(side = 2, at = seq(-0.25, 0.2, by = 0.05), pos = 0) 
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
lines(timefine, eigenfunL_CE1, type = "l", lty = "solid")
lines(timefine, eigenfunL_CE2, type = "l", lty = "dashed")
mtext("Environmental Eigenfunction (Log Scale)", side = 3, adj = 0, line = 1, font = 2)
legend("bottomright", legend= c("First EF", "Second EF"), lty = c("solid", "dashed"), bty = "n")
abline(v=0)
abline(h=0, lty = "dashed")
################################################################################################

## Selection to response

eigenval1 <- eigen(CG_funL)$values[1] # 1.352322
eigenval2 <- eigen(CG_funL)$values[2] # 0.04941301
deltaY1 <- eigenval1 * eigenfunL_CG1
deltaY2 <- eigenval2 * eigenfunL_CG2
Y_newgen1 <- aligned_logmass_mean + deltaY1
Y_newgen2 <- aligned_logmass_mean + deltaY2
                       
par(mfrow = c(1,1))
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
lines(timefine, aligned_logmass_mean, type = "l", lty = "solid", col = "red")
lines(timefine, Y_newgen1, type = "l", lty = "dashed", col = "black")
legend("bottomright", legend= c("Mean: current generation", "Mean: next generation"), lty = c("solid", "dashed"), bty = "n", col = c("red","black"))
mtext("Response to Seletction (1st Eigenfunction)", side = 3, adj = 0, line = 1, font = 2)
abline(v=0)
abline(h = 0)

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
lines(timefine, aligned_logmass_mean, type = "l", lty = "solid", col = "red")
lines(timefine, Y_newgen2, type = "l", lty = "dashed", col = "black")
legend("bottomright", legend= c("Mean: current generation", "Mean: next generation"), lty = c("solid", "dashed"), bty = "n", col = c("red","black"))
mtext("Response to Seletction (2nd Eigenfunction)", side = 3, adj = 0, line = 1, font = 2)
abline(v=0)
abline(h = 0)
