# Bootstrap Simulations for TC dataset

## load data
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
df$logx <- log10(df$trait)
id <- df$id
############################################################################################

# Log-growth Data
## Remove singular fit
sfit_ID_log <- which(sfit_sapL) # 34/800
fef_sapL_nonsing <- fef_sapL[,-sfit_ID_log]
CG_fun_sapL_nonsing <- CG_fun_sapL[,,-sfit_ID_log]
CE_fun_sapL_nonsing <- CE_fun_sapL[,,-sfit_ID_log]

## Functional Boxplot

### FE
timefine <- seq(0,1, length = 100)
par(mfrow = c(1,1), bty = "l")
FEL_nonsing <- fda::fbplot(fef_sapL_nonsing, x = timefine, xlab = "Time", ylab = "", xlim = c(0,1),
                      ylim = c(0,3), color = "grey", barcol = "lightblue")
lines(timefine,fefL_hat, lty = "solid", col = "red")
mtext("Fixed Effect (Log Scale)", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

### Eigenfunctions
eigenfun1_CGfunL <- eigen(CG_funL_hat)$vectors[,1]
eigenfun2_CGfunL <- eigen(CG_funL_hat)$vectors[,2]
eigenfun1_CEfunL <- eigen(CE_funL_hat)$vectors[,1]
eigenfun2_CEfunL <- eigen(CE_funL_hat)$vectors[,2]


GeigenL1_nonsing <- matrix(0,100,766)
EeigenL1_nonsing <- matrix(0,100,766)

GeigenL2_nonsing <- matrix(0,100,766)
EeigenL2_nonsing <- matrix(0,100,766)

for (i in 1:766){
  GeigenL1_nonsing[,i] <- eigen(CG_fun_sapL_nonsing[,,i])$vectors[,1]
  EeigenL1_nonsing[,i] <- eigen(CE_fun_sapL_nonsing[,,i])$vectors[,1]
  GeigenL2_nonsing[,i] <- eigen(CG_fun_sapL_nonsing[,,i])$vectors[,2]
  EeigenL2_nonsing[,i] <- eigen(CE_fun_sapL_nonsing[,,i])$vectors[,2]
}

par(mfrow = c(1,2), byt = "l")
GEL1_nonsing <- fda::fbplot(GeigenL1_nonsing, x = timefine, xlab = "Time", ylab = "", 
                      xlim = c(0,1), ylim = c(-0.5,0.5),
                      color = "grey", barcol = "lightblue",outliercol="orange")
lines(timefine,eigenfun1_CGfunL, lty = "solid", col = "red")
mtext("Genetic Eigenfunction 1 (Log)", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

EEL1_nonsing <- fda::fbplot(EeigenL1_nonsing, x = timefine, xlab = "Time", ylab = "", 
                         xlim = c(0,1),ylim = c(-0.5,0.5),
                         color = "grey", barcol = "lightblue",outliercol="orange")
lines(timefine,eigenfun1_CEfunL, lty = "solid", col = "red")
mtext("Environmental Eigenfunction 1 (Log)", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

GEL2_nonsing <- fda::fbplot(GeigenL2_nonsing, x = timefine, xlab = "Time", ylab = "", 
                   xlim = c(0,1), ylim = c(-0.5,0.5),
                   color = "grey", barcol = "lightblue",outliercol="orange")
lines(timefine,eigenfun2_CGfunL, lty = "solid", col = "red")
mtext("Genetic Eigenfunction 2 (Log)", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

EEL2_nonsing <- fda::fbplot(EeigenL2_nonsing, x = timefine, xlab = "Time", ylab = "", 
                   xlim = c(0,1),ylim = c(-0.5,0.5),
                   color = "grey", barcol = "lightblue",outliercol="orange")
lines(timefine,eigenfun2_CEfunL, lty = "solid", col = "red")
mtext("Environmental Eigenfunction 2 (Log)", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

################################################################################
## remove outliers
par(mfrow = c(1,1), bty = "l")
FE_300 <- fda::fbplot(fef_sapL_300, x = timefine, xlab = "Time", ylab = "", xlim = c(0,1),
                          ylim = c(0,3), color = "grey", barcol = "lightblue")
lines(timefine,fefL_hat, lty = "solid", col = "red")
mtext("Fixed Effect (Log Scale)", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

GeigenL1_300 <- matrix(0,100,300)
EeigenL1_300 <- matrix(0,100,300)
GeigenL2_300 <- matrix(0,100,300)
EeigenL2_300 <- matrix(0,100,300)

for (i in 1:300){
  GeigenL1_300[,i] <- eigen(CG_fun_sapL_300[,,i])$vectors[,1]
  EeigenL1_300[,i] <- eigen(CE_fun_sapL_300[,,i])$vectors[,1]
  GeigenL2_300[,i] <- eigen(CG_fun_sapL_300[,,i])$vectors[,2]
  EeigenL2_300[,i] <- eigen(CE_fun_sapL_300[,,i])$vectors[,2]
}

par(mfrow = c(1,2), bty = "l")
GEL1_300 <- fda::fbplot(GeigenL1_300, x = timefine, xlab = "Time", ylab = "", 
                   xlim = c(0,1), ylim = c(-0.25,0.1),
                   color = "grey", barcol = "lightblue",outliercol="orange")
lines(timefine,eigenfun1_CGfunL, lty = "solid", col = "red")
mtext("Genetic Eigenfunction 1 (Log Scale)", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

EEL1_300 <- fda::fbplot(EeigenL1_300, x = timefine, xlab = "Time", ylab = "", 
                   xlim = c(0,1),ylim = c(-0.25,0.1),
                   color = "grey", barcol = "lightblue",outliercol="orange")
lines(timefine,eigenfun1_CEfunL, lty = "solid", col = "red")
mtext("Environmental Eigenfunction 1 (Log Scale)", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

GEL2_300 <- fda::fbplot(GeigenL2_300, x = timefine, xlab = "Time", ylab = "", 
                   xlim = c(0,1), ylim = c(-0.15,0.2),
                   color = "grey", barcol = "lightblue",outliercol="orange")
lines(timefine,eigenfun2_CGfunL, lty = "solid", col = "red")
mtext("Genetic Eigenfunction 2 (Log Scale)", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

EEL2_300 <- fda::fbplot(EeigenL2_300, x = timefine, xlab = "Time", ylab = "", 
                   xlim = c(0,1),ylim = c(-0.15,0.2),
                   color = "grey", barcol = "lightblue",outliercol="orange")
lines(timefine,eigenfun2_CEfunL, lty = "solid", col = "red")
mtext("Environmental Eigenfunction 2 (Log Scale)", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

##########################################################

## SCB

## FE
par(mfrow = c(1,1), bty = "l")
FEmeanL <- rowMeans(fef_sapL_300) ## mean curve
FEL_diff<- rep(0,300)
for (i in 1:300){
  FEL_diff[i] <- max(abs(fef_sapL_300[,i] - FEmeanL))
}

FEL_bound <- quantile(FEL_diff, probs = 0.95)

plot(c(0,1), c(0,3), type = "n", xlab = "Time", ylab = "")
lines(timefine, FEmeanL, lty = "solid")
lines(timefine, FEmeanL + FEL_bound, lty = "dashed")
lines(timefine, FEmeanL - FEL_bound, lty = "dashed")
lines(timefine, fefL_hat, lty = "solid", col = "red")
mtext("Fixed Effect (Log Scale)", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

## Eigenfunctions
mean_GEL1 <- rowMeans(GeigenL1_300)
mean_EEL1 <- rowMeans(EeigenL1_300)
mean_GEL2 <- rowMeans(GeigenL2_300)
mean_EEL2 <- rowMeans(EeigenL2_300)

GeigenL_diff1 <- rep(0, 300)
EeigenL_diff1 <- rep(0, 300)
GeigenL_diff2 <- rep(0, 300)
EeigenL_diff2 <- rep(0, 300)
for (i in 1:300){
  GeigenL_diff1[i] <-  max(abs(GeigenL1_300[,i] - mean_GEL1))
  EeigenL_diff1[i] <-  max(abs(EeigenL1_300[,i] - mean_EEL1))
  GeigenL_diff2[i] <-  max(abs(GeigenL2_300[,i] - mean_GEL2))
  EeigenL_diff2[i] <-  max(abs(EeigenL2_300[,i] - mean_EEL2))
}

GeigenL1_sup1 <- quantile(GeigenL_diff1, probs = 0.95) 
EeigenL1_sup1 <- quantile(EeigenL_diff1, probs = 0.95) 
GeigenL2_sup2 <- quantile(GeigenL_diff2, probs = 0.95) 
EeigenL2_sup2 <- quantile(EeigenL_diff2, probs = 0.95)

par(mfrow = c(1,2), bty = "l")
plot(c(0,1), c(-0.25,0.1), type = "n", xlab = "Time", ylab = "")
lines(timefine, mean_GEL1, lty = "solid")
lines(timefine, mean_GEL1 + GeigenL1_sup1, lty = "dashed")
lines(timefine, mean_GEL1 - GeigenL1_sup1, lty = "dashed")
lines(timefine, eigenfun1_CGfunL, lty = "solid", col = "red")
mtext("Genetic Eigenfunction 1 (Log Scale)", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))

plot(c(0,1), c(-0.25,0.1), type = "n", xlab = "Time", ylab = "")
lines(timefine, mean_EEL1, lty = "solid")
lines(timefine, mean_EEL1 + EeigenL1_sup1, lty = "dashed")
lines(timefine, mean_EEL1 - EeigenL1_sup1, lty = "dashed")
lines(timefine, eigenfun1_CEfunL, lty = "solid", col = "red")
mtext("Environmental Eigenfunction 1 (log Scale)", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))


plot(c(0,1), c(-0.15,0.2), type = "n", xlab = "Time", ylab = "")
lines(timefine, mean_GEL2, lty = "solid")
lines(timefine, mean_GEL2 + GeigenL2_sup2, lty = "dashed")
lines(timefine, mean_GEL2 - GeigenL2_sup2, lty = "dashed")
lines(timefine, eigenfun2_CGfunL, lty = "solid", col = "red")
mtext("Genetic Eigenfunction 2 (log Scale)", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))

plot(c(0,1), c(-0.15,0.2), type = "n", xlab = "Time", ylab = "")
lines(timefine, mean_EEL2, lty = "solid")
lines(timefine, mean_EEL2 + EeigenL2_sup2, lty = "dashed")
lines(timefine, mean_EEL2 - EeigenL2_sup2, lty = "dashed")
lines(timefine, eigenfun2_CEfunL, lty = "solid", col = "red")
mtext("Environmental Eigenfunction 2 (Log Scale)", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, lty = "solid",col = rgb(0.8, 0.8, 0.8, alpha = 0.5))
#######################################################################

## Covariance functions
mean_CGL <- apply(CG_fun_sapL_300, c(1, 2), mean)
mean_CEL <- apply(CE_fun_sapL_300, c(1, 2), mean)

CGL_diff <- rep(0,300)
CEL_diff<- rep(0,300)

for (i in 1:300){
  CGL_diff[i] <- max(abs(CG_fun_sapL_300[,,i] - mean_CGL))
  CEL_diff[i] <- max(abs(CE_fun_sapL_300[,,i] - mean_CEL))
}

CGL_sup <- quantile(CGL_diff, probs = 0.95)
CEL_sup <- quantile(CEL_diff, probs = 0.95)

CGL_upper <- mean_CGL + CGL_sup
CGL_lower <- mean_CGL - CGL_sup

CEL_upper <- mean_CEL + CEL_sup
CEL_lower <- mean_CEL - CEL_sup

fig <- plot_ly(showscale = FALSE)
fig <- fig %>% add_surface(z = ~CG_funL_hat, x = ~timefine, y = ~timefine, colorscale = list(c(0, 'rgb(255, 0, 0)'), c(1, 'rgb(255, 150, 150)')))
fig <- fig %>% add_surface(z = ~CGL_upper, x = ~timefine, y = ~timefine, 
                           opacity = 0.98,colorscale = list(c(0, 'rgb(169, 169, 169)'), c(1, 'rgb(211, 211, 211)')))
fig <- fig %>% add_surface(z = ~CGL_lower, x = ~timefine, y = ~timefine, opacity = 0.98,colorscale = list(c(0, 'rgb(169, 169, 169)'), c(1, 'rgb(211, 211, 211)')))
fig <- fig %>% layout(scene = list(
  xaxis = list(title = 'Time'),
  yaxis = list(title = 'Time'),
  zaxis = list(range = c(-0.015,0.06),title = "Genetic Covariance (Log(Mass, 10<sup>-5</sup> g/t))")
))

fig

fig <- plot_ly(showscale = FALSE)
fig <- fig %>% add_surface(z = ~CE_funL_hat, x = ~timefine, y = ~timefine, colorscale = list(c(0, 'rgb(255, 0, 0)'), c(1, 'rgb(255, 150, 150)')))
fig <- fig %>% add_surface(z = ~CEL_upper, x = ~timefine, y = ~timefine, 
                           opacity = 0.98,colorscale = list(c(0, 'rgb(169, 169, 169)'), c(1, 'rgb(211, 211, 211)')))
fig <- fig %>% add_surface(z = ~CEL_lower, x = ~timefine, y = ~timefine, opacity = 0.98,colorscale = list(c(0, 'rgb(169, 169, 169)'), c(1, 'rgb(211, 211, 211)')))
fig <- fig %>% layout(scene = list(
  xaxis = list(title = 'Time'),
  yaxis = list(title = 'Time'),
  zaxis = list(range = c(-0.015,0.02),title = "Environmental Covariance (Log(Mass, 10<sup>-5</sup> g/t))")
))

fig

##########################################################################################################

par(mfrow = c(2,4))
for (k in sfit_ID_log[1:8]){
  Y_list <- split(Y_estL[,k], id) 
  plot(c(0,1), c(-0.5,3), type = "n", xlab = "Time", ylab = "Logmass", main = paste("Sample ID", k))
  for (i in 1:N){
    lines(age_list_new[[i]], Y_list[[i]], type = "l", col = i)
  }
}

for( k in sfit_ID_log[1:8]){
  fig1 <- plot_ly(showscale=FALSE, scene='scene1') 
  fig1 <- fig1 %>% add_surface(z = ~CG_fun_sapL[,,k], x = timefine, y = timefine)
  
  fig2 <- plot_ly(showscale = FALSE, scene='scene2') 
  fig2 <- fig2 %>% add_surface(z = ~CE_fun_sapL[,,k], x= timefine, y = timefine)
  
  fig_RR <- subplot(fig1, fig2) 
  fig_RR <- fig_RR %>% layout(title = paste("Sample ID", k),
                              scene = list(domain=list(x=c(0,0.45),y=c(0.25,1)),
                                           xaxis=list(title = "Time"),
                                           yaxis =list(title = "Time") , 
                                           zaxis=list(title = "Genetic Covariance (Log(Mass, 10<sup>-5</sup> g/t))"),
                                           aspectmode='cube'),
                              scene2 = list(domain=list(x=c(0.50,0.95),y=c(0.25,1)),
                                            xaxis=list(title = "Time"),
                                            yaxis =list(title = "Time"),
                                            zaxis=list(title = "Environmental Covariance (Log(Mass, 10<sup>-5</sup> g/t))"),
                                            aspectmode='cube'))
  
  print(fig_RR)
}

####################################################################################################

## Non singular
fig1 <- plot_ly(showscale=FALSE, scene='scene1') 
fig1 <- fig1 %>% add_surface(z = ~CG_fun_sapL[,,19], x = timefine, y = timefine)
  
fig2 <- plot_ly(showscale = FALSE, scene='scene2') 
fig2 <- fig2 %>% add_surface(z = ~CE_fun_sapL[,,10], x= timefine, y = timefine)
  
fig_RR <- subplot(fig1, fig2) 
fig_RR <- fig_RR %>% layout(title ="Sample 19 (Nonsingular)",
                              scene = list(domain=list(x=c(0,0.45),y=c(0.25,1)),
                                           xaxis=list(title = "Time"),
                                           yaxis =list(title = "Time") , 
                                           zaxis=list(title = "Genetic Covariance (Log(Mass, 10<sup>-5</sup> g/t))"),
                                           aspectmode='cube'),
                              scene2 = list(domain=list(x=c(0.50,0.95),y=c(0.25,1)),
                                            xaxis=list(title = "Time"),
                                            yaxis =list(title = "Time"),
                                            zaxis=list(title = "Environmental Covariance (Log(Mass, 10<sup>-5</sup> g/t))"),
                                            aspectmode='cube'))
  
print(fig_RR)

###########################################################################################

sfit_ID <- which(sfit_sap)
fef_sap_nonsing <- fef_sap[,-sfit_ID]
CG_fun_sap_nonsing <- CG_fun_sap[,,-sfit_ID]
CE_fun_sap_nonsing <- CE_fun_sap[,,-sfit_ID]

## Functional Boxplot
eigenfun1_CGfun <- eigen(CG_fun_hat)$vectors[,1]
eigenfun2_CGfun <- eigen(CG_fun_hat)$vectors[,2]

eigenfun1_CEfun <- eigen(CE_fun_hat)$vectors[,1]
eigenfun2_CEfun <- eigen(CE_fun_hat)$vectors[,2]

Geigen1 <- matrix(0,100,127)
Geigen2 <- matrix(0,100,127)
Eeigen1 <- matrix(0,100,127)
Eeigen2 <- matrix(0,100,127)
for (i in 1:127){
  Geigen1[,i] <- eigen(CG_fun_sap_nonsing[,,i])$vectors[,1]
  Geigen2[,i] <- eigen(CG_fun_sap_nonsing[,,i])$vectors[,2]
  Eeigen1[,i] <- eigen(CE_fun_sap_nonsing[,,i])$vectors[,1]
  Eeigen2[,i] <- eigen(CE_fun_sap_nonsing[,,i])$vectors[,2]
}

par(mfrow = c(1,2), byt = "l")
GE1 <- fda::fbplot(Geigen1, x = timefine, xlab = "Time", ylab = "", 
                            xlim = c(0,1), ylim = c(-0.2,0),
                            color = "grey", barcol = "lightblue",outliercol="orange")
lines(timefine,eigenfun1_CGfun, lty = "solid", col = "red")
mtext("Genetic Eigenfunction 1", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

EE1 <- fda::fbplot(Eeigen1, x = timefine, xlab = "Time", ylab = "", 
                            xlim = c(0,1),ylim = c(-0.2,0),
                            color = "grey", barcol = "lightblue",outliercol="orange")
lines(timefine,eigenfun1_CEfun, lty = "solid", col = "red")
mtext("Environmental Eigenfunction", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

GE2 <- fda::fbplot(Geigen2, x = timefine, xlab = "Time", ylab = "", 
                            xlim = c(0,1),
                            color = "grey", barcol = "lightblue",outliercol="orange")
lines(timefine,eigenfun2_CGfun, lty = "solid", col = "red")
mtext("Genetic Eigenfunction 2", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")

EE2 <- fda::fbplot(Eeigen2, x = timefine, xlab = "Time", ylab = "", 
                            xlim = c(0,1),
                            color = "grey", barcol = "lightblue",outliercol="orange")
lines(timefine,eigenfun2_CEfun, lty = "solid", col = "red")
mtext("Environmental Eigenfunction 2", side = 3, adj = 0, line = 1, font = 2)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
