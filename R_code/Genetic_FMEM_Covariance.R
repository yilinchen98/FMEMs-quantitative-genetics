# Fit FMEM with smoothed data at original measurement points

## load data
setwd("D:/KCL_2023-2027_PhD/Genetics_FMEM_Project/FMEMs-quantitative-genetics/R_code")
TRFUN25PUP4 = read.delim("TRFUN25PUP4.DAT",header = FALSE)
names(TRFUN25PUP4)<-c("id","sire","dam","trait","x")
df <- data.frame(TRFUN25PUP4)

## calculate genetic relationship matrix
FirstUniqueIdPos <- which(duplicated(df$id) == FALSE)
pos = df$id[FirstUniqueIdPos] # extract ids for all subjects
sire_id = df$sire[FirstUniqueIdPos] # extract ids for sire
dam_id = df$dam[FirstUniqueIdPos] # extract ids for dam

pede <- editPed(sire = sire_id, dam = dam_id, label = pos)
ped<- with(pede, pedigree(label=label, sire=sire, dam=dam))
A <- getA(ped)[163:1035,163:1035]

## Rescale time interval to [0,1]
N = length(FirstUniqueIdPos) # N = 873 subjects
n = length(df$id) # n = 6860 observations
age_list <- split(df$x,df$id)
trait_list <- split(df$trait,df$id)

### x = (x -min(x))/(max(x) - min(x))
age_list_new <- list()
for (i in 1:N){
  age_list_new[[i]] = (age_list[[i]]-min(age_list[[i]]))/(max(age_list[[i]])-min(age_list[[i]]))
}

df$x_rescaled <- unsplit(age_list_new,df$id)

## Smoothe data, choose smoothing parameter < 10^-4
agefine <- seq(0,1,length=100) # dense time grid
logmass <- matrix(0,100,N) # store smooth log growth curves
pred_mass <- matrix(0,100,N) # store the smoothed mass recovered by taking exp
lam <- rep(0,N) # smoothing parameter used for each log growth curve

for (i in 1:N){
  ss_logmass <- smooth.spline(age_list_new[[i]], log(trait_list[[i]]), cv=FALSE,
                              all.knots=TRUE)
  logmass[,i] <- predict(ss_logmass, agefine)$y
  pred_mass[,i] <- exp(predict(ss_logmass,agefine)$y)
  lam[i] <- ss_logmass$lambda
}

large_lam <- which(lam > 1e-4)
mass_adjusted <- matrix(0,100, length(large_lam))
for ( i in 1: length(large_lam)){
  ss_new <- smooth.spline(age_list_new[[large_lam[i]]], log(trait_list[[large_lam[i]]]), 
                          all.knots = TRUE, lambda = 1e-4) # restrict lambda <=1e-4
  mass_adjusted[,i] <- exp(predict(ss_new, agefine)$y)
}
colnames(mass_adjusted) <- as.character(large_lam)

## Let reform the smoothed growth curves into a matrix
trait_hat <- pred_mass 
for (i in 1: length(large_lam)){
  trait_hat[,large_lam[i]] <- mass_adjusted[,i]
}

## curve alignment
aligned_mass_process <- time_warping(f=trait_hat, time=agefine)
aligned_mass_curve <- aligned_mass_process$fn
aligned_mean <- aligned_mass_process$fmean
warping_funs <- aligned_mass_process$warping_functions

## FPCA
fpcaobj_mass <- prcomp(x=t(aligned_mass_curve), retx = TRUE, center = TRUE, rank. = 3)
eigen_mass <- fpcaobj_mass$rotation # eigen vectors

phi_list <- list() 
# create an empty list which stores eigenfunctions for 873 subjects
# evaluated at the original time points.

for (i in 1:N){
  phi <- convert_to_basisfunctions(t = agefine, eigenvecs = eigen_mass[,1:3],
                                   tout = age_list_new[[i]])
  phi_list[[i]] <- phi
}

phi <- do.call(rbind,phi_list)
colnames(phi) <- c("phi1", "phi2", "phi3")

## Interpolate curve
trait_hat_list <- list()

for (i in 1:N){
  inter_trait <- convert_to_basisfunctions(t = agefine, eigenvecs = trait_hat[,i],
                                           tout = age_list_new[[i]])
  trait_hat_list[[i]] <- inter_trait
}

## Reform dataframe
df_test <- data.frame(id = df$id, trait = unsplit(trait_hat_list, df$id), phi)

## Extract FMEM
fmmForm <- trait ~ -1 + df_test$phi1 + df_test$phi2 + df_test$phi3 + (-1 + df_test$phi1 + df_test$phi2 + df_test$phi3 | df_test$id) + 
  (-1 + df_test$phi1 + df_test$phi2 + df_test$phi3 | df_test$id) 

system.time(
  ff <- fit_genetic_fmm(fmmForm, df_test, A, eigen_mass)
) # user   system  elapsed

summary(ff)

## Extract fixed effect
fcoefs <- fixef(ff)
fixed_effect <- fcoefs[1] * eigen_mass[,1] + fcoefs[2] * eigen_mass[,2] + fcoefs[3] * eigen_mass[,3]
par(mfrow = c(1,1))
plot(agefine, fixed_effect, type = "l", xlab = "t", ylab="mass")

## Extract covariance 
vc <- VarCorr(ff)
CG <- vc[["df_test.id"]] ## estimated genetic covariance
CE <- vc[["df_test.id.1"]] ## estimated environmental covariance

CG_fun <- eigen_mass %*% CG %*% t(eigen_mass)
CE_fun <- eigen_mass %*% CE %*% t(eigen_mass)

fig1 <- plot_ly(showscale=FALSE, scene='scene1') 
fig1 <- fig1 %>% add_surface(z = ~CG_fun)

fig2 <- plot_ly(showscale = FALSE, scene='scene2') 
fig2 <- fig2 %>% add_surface(z = ~CE_fun)

fig_RR1 <- subplot(fig1, fig2) 
fig_RR1 <- fig_RR1 %>% layout(title = "Covariance Plot",
                              scene = list(domain=list(x=c(0,0.45),y=c(0.25,1)),
                                           xaxis=list(title = "time"),
                                           yaxis =list(title = "time") , 
                                           zaxis=list(range=c(-10,900),title = "Gen"),
                                           aspectmode='cube'),
                              scene2 = list(domain=list(x=c(0.50,0.95),y=c(0.25,1)),
                                            xaxis=list(title = "time"),
                                            yaxis =list(title = "time"),
                                            zaxis=list(range=c(-10,900),title = "Env"),
                                            aspectmode='cube'))

fig_RR1