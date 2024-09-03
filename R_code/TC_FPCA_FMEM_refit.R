# Fit FMEM with smoothed data at original measurement points

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
df$logx <- log10(df$trait)

## plot growth curve
par(mar = c(5, 6, 4, 2) + 0.1)
plot(c(1, 25), c(0, 400), type = "n", 
     xlab = "Days", 
     ylab = expression(Mass~(10^-5~g)),
     xlim = c(0, 25), ylim = c(0, 400), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 25, by = 5), pos = 0) 
axis(side = 2, at = seq(0, 400, by = 100), pos = 0) 
box()
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(age_list[[i]], trait_list[[i]], type = "l", col = i)
}


##################################################################################################

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
box()
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(age_list_new[[i]], log10(trait_list[[i]]), type = "l", col = i)
}

par(mar = c(5, 6, 4, 2) + 0.1)
plot(c(0,1), c(0, 3), type = "n", 
     xlab = "Time", 
     ylab = expression(Log(Mass,~10^-5~g)),
     xlim = c(0, 1), ylim = c(0,3), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0) 
axis(side = 2, at = seq(0, 3, by = 0.5), pos = 0) 
box()
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(timefine, pred_logmass_fine[,i], type = "l", col = i)
}

## plot raw growth curves vs plot smoothed growth curve
par(mar = c(5, 6, 4, 2) + 0.1)
plot(c(0,1), c(0, 3), type = "n", 
     xlab = "Time", 
     ylab = expression(Mass~(10^-5~g)),
     xlim = c(0, 1), ylim = c(0,400), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0) 
axis(side = 2, at = seq(0, 400, by = 100), pos = 0) 
box()
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(age_list_new[[i]], trait_list[[i]], type = "l", col = i)
}

par(mar = c(5, 6, 4, 2) + 0.1)
plot(c(0,1), c(0, 3), type = "n", 
     xlab = "Time", 
     ylab = expression(Mass~(10^-5~g)),
     xlim = c(0, 1), ylim = c(0,400), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0) 
axis(side = 2, at = seq(0, 400, by = 100), pos = 0) 
box()
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(timefine, pred_mass_fine[,i], type = "l", col = i)
}


#######################################################################################################

## curve alignment

### Align mass curves
aligned_mass_process <- time_warping(f=pred_mass_fine, time=timefine)
aligned_mass_curve <- aligned_mass_process$fn
aligned_mean <- aligned_mass_process$fmean
warping_funs <- aligned_mass_process$warping_functions

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
box()
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(timefine, warping_funs[,i], type = "l", col = i)
}
mtext("Warping Functions", side = 3, adj = 0, line = 1, font = 2)

plot(c(0,1), c(0, 3), type = "n", 
     xlab = "Time", 
     ylab = expression(Mass~(10^-5~g)),
     xlim = c(0, 1), ylim = c(0,400), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0) 
axis(side = 2, at = seq(0, 400, by = 100), pos = 0) 
box()
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(timefine, aligned_mass_curve[,i], type = "l", col = i)
}
mtext("Aligned Mass Curves", side = 3, adj = 0, line = 1, font = 2)

###########################################################################################

# Align logmass curves
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
box()
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(timefine, warping_logmass_funs[,i], type = "l", col = i)
}
mtext("Warping Functions", side = 3, adj = 0, line = 1, font = 2)

plot(c(0,1), c(0, 3), type = "n", 
     xlab = "Time", 
     ylab = expression(Log(Mass,~10^-5~g)),
     xlim = c(0, 1), ylim = c(0,3), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0) 
axis(side = 2, at = seq(0, 3, by = 0.5), pos = 0) 
box()
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(timefine, aligned_logmass_curve[,i], type = "l", col = i)
}
mtext("Aligned Log Mass Curves", side = 3, adj = 0, line = 1, font = 2)

####################################################################################################

## FPCA

## growth curves
fpcaobj_mass <- prcomp(x=t(aligned_mass_curve), retx = TRUE, center = TRUE, rank. = 4)
pcs_mass <- fpcaobj_mass$rotation # eigen vectors

par(mfrow = c(1,1))
fviz_eig(fpcaobj_mass, 
         addlabels = TRUE, 
         main="", ncp = 6)

# Create a data frame for the first two principal components
df_pcs_mass <- data.frame(
  Time = rep(timefine, 2),  # Repeat the time vector for two components
  PC_Value = c(pcs_mass[, 1], pcs_mass[, 2]),  # Stack the first two principal components
  Component = factor(rep(1:2, each = length(timefine)))  # Indicate which principal component each value belongs to
)

p1 <- ggplot(df_pcs_mass[df_pcs_mass$Component == 1, ], aes(x = Time, y = PC_Value)) +
  geom_line() +
  labs(x = "Time", y = "", title = "Principal Component 1") +
  theme_minimal()

p2 <- ggplot(df_pcs_mass[df_pcs_mass$Component == 2, ], aes(x = Time, y = PC_Value)) +
  geom_line() +
  labs(x = "Time", y = "", title = "Principal Component 2") +
  theme_minimal()
combined_plot <- p1 / p2 
print(combined_plot)

## log-growth curvs
fpcaobj_logmass <- prcomp(x=t(aligned_logmass_curve), retx = TRUE, center = TRUE, rank. = 4)
pcs_logmass <- fpcaobj_logmass$rotation # eigen vectors

par(mfrow = c(1,1))
fviz_eig(fpcaobj_logmass, 
         addlabels = TRUE, 
         main="", ncp = 6)

df_pcs_logmass <- data.frame(
  Time = rep(timefine, 2),  # Repeat the time vector for two components
  PC_Value = c(pcs_logmass[, 1], pcs_logmass[, 2]),  # Stack the first two principal components
  Component = factor(rep(1:2, each = length(timefine)))  # Indicate which principal component each value belongs to
)

p1 <- ggplot(df_pcs_logmass[df_pcs_logmass$Component == 1, ], aes(x = Time, y = PC_Value)) +
  geom_line() +
  labs(x = "Time", y = "", title = "Principal Component 1") +
  theme_minimal()

p2 <- ggplot(df_pcs_logmass[df_pcs_logmass$Component == 2, ], aes(x = Time, y = PC_Value)) +
  geom_line() +
  labs(x = "Time", y = "", title = "Principal Component 2") +
  theme_minimal()
combined_plot <- p1 / p2 
print(combined_plot)

## Model fitting (original scale)

### Align raw data
aligned_trait <- list()
gamma_mass <- list()
for (i in 1:N){
  gamma_mass_inter <- convert_to_basisfunctions(timefine, warping_funs[,i], age_list_new[[i]])
  gamma_mass[[i]] <- gamma_mass_inter
  aligned_trait[[i]] <- warp_f_gamma(trait_list[[i]], age_list_new[[i]], gamma_mass_inter)
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
box()
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(age_list_new[[i]], gamma_mass[[i]], type = "l", col = i)
}
mtext("Warping Functions", side = 3, adj = 0, line = 1, font = 2)

plot(c(0,1), c(0, 3), type = "n", 
     xlab = "Time", 
     ylab = expression(Mass~(10^-5~g)),
     xlim = c(0, 1), ylim = c(0,400), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0) 
axis(side = 2, at = seq(0, 400, by = 100), pos = 0) 
box()
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(age_list_new[[i]], aligned_trait[[i]], type = "l", col = i)
}
mtext("Aligned Raw Mass Curves", side = 3, adj = 0, line = 1, font = 2)


#############################################################################################

phi_list <- list() 
# create an empty list which stores eigenfunctions for 873 subjects
# evaluated at the original time points.

for (i in 1:N){
  phi <- convert_to_basisfunctions(t = timefine, eigenvecs = pcs_mass,
                                   tout = age_list_new[[i]])
  phi_list[[i]] <- phi
}

phi <- do.call(rbind,phi_list)
colnames(phi) <- c("phi1", "phi2", "phi3", "phi4")

## Reform dataframe
df_mass <- data.frame(id =id, trait = unsplit(aligned_trait,id), phi)

fmmForm <- trait ~ -1 + df_mass$phi1 + df_mass$phi2 + df_mass$phi3 + df_mass$phi4 +
  (-1 + df_mass$phi1 + df_mass$phi2 | df_mass$id) + 
  (-1 + df_mass$phi1 + df_mass$phi2 | df_mass$id) 

system.time(
  ff <- fit_genetic_fmm(fmmForm, df_mass, A, 2)
) # user   system  elapsed

summary(ff)

isSingular(ff, tol = 1e-5) # original tol = 1e-4

## fixed effect
beta <- fixef(ff)
fef <- pcs_mass %*% beta

df_FEmass <- data.frame(
  Time = rep(timefine, 2),
  Value = c(aligned_mean, fef),
  FE = rep(c("mean", "estimated fixed effect"), each = length(timefine))
)

p<-ggplot(df_FEmass, aes(x = Time, y = Value, linetype = FE, color = FE)) +
  geom_line() +  
  labs(x = "Time", y = expression(Mass~(10^-5~g)), title = "") +
  scale_color_manual(values = c("mean" = "red", "estimated fixed effect" = "black")) +  # Assign colors to lines
  scale_linetype_manual(values = c("mean" = "solid", "estimated fixed effect" = "dashed")) +  # Assign line types
  theme_minimal() +
  theme(legend.position = c(0.8,0.2))

## random effects

## Extract covariance
vc <- VarCorr(ff)
CG <- vc[["df_mass.id"]] # genetic covariance
CE <- vc[["df_mass.id.1"]] # environmental covariance

CG_fun <- pcs_mass[,1:2] %*% CG %*% t(pcs_mass[,1:2]) # estimated gen cov function
CE_fun <- pcs_mass[,1:2] %*% CE %*% t(pcs_mass[,1:2]) # estimated env cov function

P_fun <- CG_fun + CE_fun

fig1 <- plot_ly(showscale=FALSE, scene='scene1') 
fig1 <- fig1 %>% add_surface(z = ~CG_fun, x = timefine, y = timefine)

fig2 <- plot_ly(showscale = FALSE, scene='scene2') 
fig2 <- fig2 %>% add_surface(z = ~CE_fun, x= timefine, y = timefine)

fig_RR1 <- subplot(fig1, fig2) 
fig_RR1 <- fig_RR1 %>% layout(title = "",
                              scene = list(domain=list(x=c(0,0.45),y=c(0.25,1)),
                                           xaxis=list(title = "Time"),
                                           yaxis =list(title = "Time") , 
                                           zaxis=list(range=c(-10,600),title = "Genetic Covariance (Mass, 10<sup>-5</sup> g/t)"),
                                           aspectmode='cube'),
                              scene2 = list(domain=list(x=c(0.50,0.95),y=c(0.25,1)),
                                            xaxis=list(title = "Time"),
                                            yaxis =list(title = "Time"),
                                            zaxis=list(range=c(-10,600),title = "Environmental Covariance (Mass, 10<sup>-5</sup> g/t)"),
                                            aspectmode='cube'))

fig_RR1

fig1 <- plot_ly(showscale=FALSE, scene='scene1') 
fig1 <- fig1 %>% add_surface(z = ~CG_fun, x = timefine, y = timefine)

fig2 <- plot_ly(showscale = FALSE, scene='scene2') 
fig2 <- fig2 %>% add_surface(z = ~P_fun, x= timefine, y = timefine)

fig_RR1 <- subplot(fig1, fig2) 
fig_RR1 <- fig_RR1 %>% layout(title = "",
                              scene = list(domain=list(x=c(0,0.45),y=c(0.25,1)),
                                           xaxis=list(title = "Time"),
                                           yaxis =list(title = "Time") , 
                                           zaxis=list(range=c(-10,1200),title = "Genetic Covariance (Mass, 10<sup>-5</sup> g/t)"),
                                           aspectmode='cube'),
                              scene2 = list(domain=list(x=c(0.50,0.95),y=c(0.25,1)),
                                            xaxis=list(title = "Time"),
                                            yaxis =list(title = "Time"),
                                            zaxis=list(range=c(-10,1200),title = "Phenotypic Covariance (Mass, 10<sup>-5</sup> g/t)"),
                                            aspectmode='cube'))

fig_RR1



##########################################################################################

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
box()
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(age_list_new[[i]], gamma_logmass[[i]], type = "l", col = i)
}
mtext("Warping Functions", side = 3, adj = 0, line = 1, font = 2)

plot(c(0,1), c(0, 3), type = "n", 
     xlab = "Time", 
     ylab = expression(Mass~(10^-5~g)),
     xlim = c(0, 1), ylim = c(0,3), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0) 
axis(side = 2, at = seq(0, 3, by = 0.5), pos = 0) 
box()
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
for (i in 1:N){
  lines(age_list_new[[i]], aligned_logtrait[[i]], type = "l", col = i)
}
mtext("Aligned Raw Log Mass Curves", side = 3, adj = 0, line = 1, font = 2)

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

fmmFormL <- trait ~ -1 + df_logmass$phi1 + df_logmass$phi2 + df_logmass$phi3 + df_logmass$phi4 +
  (-1 + df_logmass$phi1 + df_logmass$phi2 | df_logmass$id) + 
  (-1 + df_logmass$phi1 + df_logmass$phi2 | df_logmass$id) 

system.time(
  ffL <- fit_genetic_fmm(fmmFormL, df_logmass, A, 2)
) # user   system  elapsed

summary(ffL)

isSingular(ffL) # original tol = 1e-4

## fixed effect
betaL <- fixef(ffL)
fefL <- pcs_logmass %*% betaL

# Create a data frame
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
ggsave("FElogmass.png", plot = pfe, width = 6, height = 2, dpi = 300)
## Extract covariance
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
                                            zaxis=list(range=c(-0.002,0.002),title = "Environmental Covariance (Log(Mass, 10<sup>-5</sup> g/t))"),
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

######################################################################
eigenfun1_CGfunL <- Eigen(CG_funL)$vectors[,1]
eigenfun1_CEfunL <- Eigen(CE_funL)$vectors[,1]

par(mfrow = c(1,1))
par(mar = c(5, 6, 4, 2) + 0.1)
plot(c(0,1), c(-0.5,0.5), type = "n", 
     xlab = "Time", 
     ylab = "",
     xlim = c(0, 1), ylim = c(-0.5, 0.5), 
     xaxs = "i", yaxs = "i",
     axes = FALSE) 
axis(side = 1, at = seq(0, 1, by = 0.1), pos = 0) 
axis(side = 2, at = seq(-0.5, 0.5, by = 0.1), pos = 0) 
box()
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "solid")
lines(timefine, eigenfun1_CEfunL)


pc_CGL <- prcomp(CG_funL)
