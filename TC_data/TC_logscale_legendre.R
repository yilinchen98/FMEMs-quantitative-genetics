# Load Data
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

N = length(FirstUniqueIdPos) # N = 873 subjects
n = length(id) # n = 6860 observations
age_list <- split(df$x,id)
trait_list <- split(df$trait,id)

## Rescale time interval to [-1,1]
### x_recaled = -1 + 2 * (x - x_min)/(x_max - x_min)
age_list_new <- list()
for (i in 1:N){
  age_list_new[[i]] = -1 + 2*(age_list[[i]]-min(age_list[[i]]))/(max(age_list[[i]])-min(age_list[[i]]))
}

df$x_rescaled <- unsplit(age_list_new,id)
df$logtrait <- log10(df$trait)
log_trait_list <- split(df$logtrait,id)

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
  lines(age_list_new[[i]], log_trait_list[[i]],type = "l", col = i)
}
axis(side = 1, at = seq(-1, 1, by = 0.2), pos = 0) 
axis(side = 2, at = seq(0, 3, by = 0.5), pos = 0) 
mtext("Log Mass Curves", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, v=0)
#############################################################################

# Redo the analysis of Irwin and Carter 2013

## Orthogonal Legendre Basis
library(sommer)
legendrePoly3_1 <- leg(age_list_new[[1]], n = 3)

## Example of legendre ploynomial
par(mfrow = c(1,1))
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
  lines(age_list_new[[1]], legendrePoly3_1[,i], col = i)
}
mtext("Legendre Poly of order 3", side = 3, adj = 0, line = 1, font = 2)
abline(h = 0, v=0)

########################################################################

legBasis3 <- list() # 3rd order Legendre polynomial
for (i in 1:N){
  legBasis3[[i]] <- leg(age_list_new[[i]], n = 3)
}

basis <- do.call(rbind, legBasis3)

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

isSingular(ffLeg) # TRUE

summary(ffLeg)
theta <- ffLeg@theta
theta
################################################################

## Model results
timefine <- seq(-1,1,length = 100)
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

## fixed effect

basisfine <- leg(timefine, n = 3)
betaLeg <- fixef(ffLeg)
fefLeg <-  basisfine %*% betaLeg
smoothed_mean <- rowMeans(pred_logmass_fine)
df_FE <- data.frame(
  Time = rep(timefine, 2),
  Value = c(smoothed_mean, fefLeg),
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
vcLeg <- VarCorr(ffLeg)
CGLeg <- vcLeg[["df_leg.id"]] # genetic covariance
CELeg <- vcLeg[["df_leg.id.1"]] # environmental covariance

CG_funLeg <- basisfine%*% CGLeg %*% t(basisfine) # estimated gen cov function
CE_funLeg <- basisfine %*% CELeg %*% t(basisfine) # estimated env cov function

fig3 <- plot_ly(showscale=FALSE, scene='scene1') 
fig3 <- fig3 %>% add_surface(z = ~CG_funLeg, x = timefine, y = timefine)

fig4 <- plot_ly(showscale = FALSE, scene='scene2') 
fig4 <- fig4 %>% add_surface(z = ~CE_funLeg, x= timefine, y = timefine)

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

P_funLeg <- CG_funLeg + CE_funLeg
fig5 <- plot_ly(showscale = FALSE, scene='scene2') 
fig5 <- fig5 %>% add_surface(z = ~P_funLeg, x= timefine, y = timefine)

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

##################################################################################

## Use 2-order Legendre poly
fformLegendre2 <- trait ~ -1 + df_leg$phi1 + df_leg$phi2 + df_leg$phi3 + df_leg$phi4 + 
  (-1 + df_leg$phi1 + df_leg$phi2 + df_leg$phi3 | df_leg$id) + 
  (-1 + df_leg$phi1 + df_leg$phi2 + df_leg$phi3 | df_leg$id)

system.time(
  ffLeg2 <- fit_genetic_fmm(fformLegendre2, df_leg, A, 3)
) # user   system  elapsed

isSingular(ffLeg2) # TRUE

summary(ffLeg2)

## Random effects
vcLeg2 <- VarCorr(ffLeg2)
CGLeg2 <- vcLeg2[["df_leg.id"]] # genetic covariance
CELeg2 <- vcLeg2[["df_leg.id.1"]] # environmental covariance

CG_funLeg2 <- basisfine[,1:3]%*% CGLeg2 %*% t(basisfine[,1:3]) # estimated gen cov function
CE_funLeg2 <- basisfine[,1:3] %*% CELeg2 %*% t(basisfine[,1:3]) # estimated env cov function

fig3 <- plot_ly(showscale=FALSE, scene='scene1') 
fig3 <- fig3 %>% add_surface(z = ~CG_funLeg2, x = timefine, y = timefine)

fig4 <- plot_ly(showscale = FALSE, scene='scene2') 
fig4 <- fig4 %>% add_surface(z = ~CE_funLeg2, x= timefine, y = timefine)

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

P_funLeg2 <- CG_funLeg2 + CE_funLeg2
fig5 <- plot_ly(showscale = FALSE, scene='scene2') 
fig5 <- fig5 %>% add_surface(z = ~P_funLeg2, x= timefine, y = timefine)

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
