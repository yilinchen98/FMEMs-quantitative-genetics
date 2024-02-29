library(refund)
library(sparseFLMM)
library(data.table)

data("acoustic_subset") # example of speach data in the package

# Modify our dataset to used sparseFLMM()

setwd("D:/KCL_2023-2027_PhD/Year 1/Genetics/FMEMs-quantitative-genetics/R_code")
TRFUN25PUP4 = read.delim("TRFUN25PUP4.DAT",header = FALSE)
names(TRFUN25PUP4)<-c("id","sire","dam","trait","x")
TC_df <- data.frame(TRFUN25PUP4)

FirstUniqueIdPos <- which(duplicated(TC_df$id) == FALSE)
N = length(FirstUniqueIdPos) # N = 873 subjects
n = length(TC_df$id) # n = 6860 observations
age_list <- split(TC_df$x,TC_df$id)
trait_list <- split(TC_df$trait,TC_df$id)

## Rescale time interval to [0,1]
## x = (x -min(x))/(max(x) - min(x))
age_list_new <- list()
for (i in 1:N){
  age_list_new[[i]] = (age_list[[i]]-min(age_list[[i]]))/(max(age_list[[i]])-min(age_list[[i]]))
}

TC_df$x <- unsplit(age_list_new,df$id)

parsedCurveInfo <- TC_df[,c(1,4,5)]
col_names <- c("subject_long", "y_vec", "t")
colnames(parsedCurveInfo) <- col_names
parsedCurveInfo$n_long <- parsedCurveInfo[,1] 
parsedCurveInfo <- as.data.table(parsedCurveInfo)

fit_sparseFLMM <- sparseFLMM(curve_info = parsedCurveInfo, use_RI = TRUE, 
                             use_simple = TRUE,
                             num_covariates = 3, covariate_form = rep("by",3), 
                             interaction = FALSE,
                             bf_covs = c(4),m_covs=list(c(2,3)))

# Visualise mean_hat
mean_hat <- fit_sparseFLMM$mean_hat

agefine <- seq(0,1,length=100)

par(mfrow=c(1,1))
plot(agefine, mean_hat$mean_pred, type="l", col="blue", xlab="time", ylab="mass")
lines(mean_mass, type="l", col="red")
legend(x="bottomright", legend=c("mean_fda", "mean_sparseFLMM"), 
       col=c("red","blue"), lwd = 1)
