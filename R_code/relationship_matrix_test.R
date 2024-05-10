library(pedigreemm)
library(lme4)
library(Matrix)

# Import data
setwd("D:/KCL_2023-2027_PhD/Year 1/Genetics/FMEMs-quantitative-genetics/R_code")
TRFUN25PUP4 = read.delim("TRFUN25PUP4.DAT",header = FALSE)
names(TRFUN25PUP4)<-c("id","sire","dam","trait","x")
df <- data.frame(TRFUN25PUP4) 

indices <- df$id
trait <- df$trait
time <- df$x
FirstUniqueIdPos <- which(duplicated(indices) == FALSE)

pos = indices[FirstUniqueIdPos] # extract ids for all subjects
sire_id = df$sire[FirstUniqueIdPos] # extract ids for sire
dam_id = df$dam[FirstUniqueIdPos] # extract ids for dam

pede <- editPed(sire = sire_id, dam = dam_id, label = pos)
ped<- with(pede, pedigree(label=label, sire=sire, dam=dam))
A <- getA(ped)[163:1035,163:1035]

######################################################################
posFirstUniqueId = which((!duplicated(df$id))==TRUE)
sr = c( rep(NA, 100),
        rep(NA, 9900),
        rep(NA, 1500))
sr[unique(df$id)] = (df$sire)[posFirstUniqueId]

dm = c( rep(NA, 100),
        rep(NA, 9900),
        rep(NA, 1500))
dm[unique(df$id)] = (df$dam)[posFirstUniqueId]

lbl = c(1:11500)
pedigree = pedigree(sr, dm, lbl)

position = unique(df$id)[order(unique(df$id))]
A1 = getA(pedigree)[position,position]

all(A == A1)


