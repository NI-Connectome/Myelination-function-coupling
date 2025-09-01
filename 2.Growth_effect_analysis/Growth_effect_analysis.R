# this code compute the growth effect of MFC
# use the "install.packages()" function to add the required packages before running this code
# Please set "YourPath before" running this code

library(mgcv)
library(R.matlab)
library(gratia)

MFC <- readMat("YourPath/data/MFC.mat")  #MFC：Vertex-level gMFC/sMFC/MFC
Age <- readMat("YourPath/data/Age_364.mat") 
GAM_data <- list()
GAM_data$PMA <- Age$PMA364     # PMA at scan
GAM_data$MFC_gs <- MFC$gMFC.v  # gMFC
GAM_data$MFC_ss <- MFC$sMFC.v  # sMFC
GAM_data$MFC_gss <- MFC$MFC.v  # MFC
fit_GAM <- list()              # create an empty list
current_GAM_data <- list()     # create an empty list
current_GAM_data$PMA <- GAM_data$PMA
Pvalue_all <- matrix(1:8589, nrow = 8589, ncol = 1)
derivs <- matrix(1:8589, nrow = 8589, ncol = 1)


### the growth effect of gMFC
for (i in 1:nrow(GAM_data$MFC_gs)) {
  current_GAM_data$current_MFC <- GAM_data$MFC_gs[i, ]
  fit_GAM <- gam(current_MFC ~ s(PMA, bs="tp", k=3), data=current_GAM_data, method="REML")  # GAM model: thin plate regression splined
  summary_model <- summary(fit_GAM)  # GAM model fitting result
  Pvalue_all[i,1] = summary_model$s.table[1,4] 
  deriv <- derivatives(fit_GAM, term='s(PMA)', order=1)  # the first derivative of th eage smooth function
  derivs[i,1] <- mean(deriv$.derivative)   #average derivative: growth rate
}
write.table(Pvalue_all,"YourPath/result/Pvalue_gs.csv",row.names=F, col.names=F)
write.table(derivs,"YourPath/result/derivs_gs.csv",row.names=F, col.names=F)


### the growth effect of sMFC
for (i in 1:nrow(GAM_data$MFC_ss)) {
  current_GAM_data$current_MFC <- GAM_data$MFC_ss[i, ]
  fit_GAM <- gam(current_MFC ~ s(PMA, bs="tp", k=3), data=current_GAM_data, method="REML")   # GAM model: thin plate regression splined
  summary_model <- summary(fit_GAM) # GAM model fitting result
  Pvalue_all[i,1] = summary_model$s.table[1,4]
  deriv <- derivatives(fit_GAM, term='s(PMA)', order=1) # the first derivative of th eage smooth function
  derivs[i,1] <- mean(deriv$.derivative)   # average derivative: growth rate
}

write.table(Pvalue_all,"YourPath/result/Pvalue_ss.csv",row.names=F, col.names=F)
write.table(derivs,"YourPath/result/derivs_ss.csv",row.names=F, col.names=F)


### the growth effect of MFC
for (i in 1:nrow(GAM_data$MFC_gss)) {
  current_GAM_data$current_MFC <- GAM_data$MFC_gss[i, ]
  fit_GAM <- gam(current_MFC ~ s(PMA, bs="tp", k=3), data=current_GAM_data, method="REML")   # GAM model: thin plate regression splined
  summary_model <- summary(fit_GAM) # GAM model fitting result
  Pvalue_all[i,1] = summary_model$s.table[1,4]
  deriv <- derivatives(fit_GAM, term='s(PMA)', order=1) # the first derivative of th eage smooth function
  derivs[i,1] <- mean(deriv$.derivative)   # average derivative: growth rate
}

write.table(Pvalue_all,"YourPath/result/Pvalue_gss.csv",row.names=F, col.names=F)
write.table(derivs,"YourPath/result/derivs_gss.csv",row.names=F, col.names=F)
