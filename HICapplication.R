# Application HIC data

rm(list=ls())
require("betareg")
require("xtable")
source("RobustEstimators.R")

data <- read.table("HICdata.txt")
head(data)

data$GDP <- data$GDP/1000
# Complete model: Urb and GDP in both mean and precision submodels 

y <- data$HIC
X <- as.matrix(cbind(rep(1,80), data[,c(1, 2)]))
Z <- as.matrix(cbind(rep(1,80), data[,c(1, 2)]))

fit_MLE <- betareg(y ~ X[,-1] | Z[,-1])
summary(fit_MLE)

fit_SMLE <- SMLE_BETA(y, X, Z) #0.06
fit_SMLE

fit_MDPDE <- MDPDE_BETA(y, X, Z) #0.42
fit_MDPDE

fit_LSMLE <- LSMLE_BETA(y, X, Z, weights = TRUE) #0.06
fit_LSMLE

fit_LMDPDE <- LMDPDE_BETA(y, X, Z) #0.06
fit_LMDPDE

# All the fits indicate that Urb may be excluded from the precision submodel
# Reduced model: Urb out of the precision submodel

Z <- as.matrix(cbind(rep(1,80), data[,c(2)]))

fit_MLE <- betareg(HIC ~ . | GDP, data = data)
summary(fit_MLE)

fit_SMLE <- SMLE_BETA(y, X, Z)
fit_SMLE

fit_MDPDE <- MDPDE_BETA(y, X, Z)
fit_MDPDE

# SMLE and MDPDE: Lack of stability

fit_LSMLE <- LSMLE_BETA(y, X, Z, weights = TRUE)
fit_LSMLE

fit_LMDPDE <- LMDPDE_BETA(y, X, Z)
fit_LMDPDE

# Similar results for the LMDPDE and the LSMLE
# Continue analysis with the LSMLE

# Excluding outlier observation # 1

fit_MLEwo1 <- betareg(HIC ~ . | GDP, data = data[-1,])
summary(fit_MLEwo1)

fit_LSMLEwo1 <- LSMLE_BETA(y[-1], X[-1,], Z[-1,], weights = TRUE)
fit_LSMLEwo1

fit_LMDPDEwo1 <- LMDPDE_BETA(y[-1], X[-1,], Z[-1,], weights = TRUE)
fit_LMDPDEwo1

# Residuals
source("Resfunction.r")

weights <- fit_LSMLE$weights
RP2_smle <- residuals_beta(y, X, Z, c(fit_LSMLE$beta,fit_LSMLE$gama), linkmu="logit", linkphi="log") 

plot(RP2_smle, weights, xlab="Residuals", ylab="Weights",main=""
     , pch = "+", xlim=c(-3,10), ylim=c(0,1.3), las=1, cex=1.5, cex.lab=1.5, cex.axis=1.4, cex.main=2.0)
identify(RP2_smle,weights,cex=1.3)

#zoomed versions
plot(RP2_smle, weights, xlab="Residuals", ylab="Weights",main=""
     , pch = "+", xlim=c(-3,3), ylim=c(0.8,1.2), las=1, cex=1.5, cex.lab=1.5, cex.axis=1.4, cex.main=2.0)
identify(RP2_smle,weights,cex=1.3)


# Envelope
source("envelope_function.r")
envelope_LSMLE(y=y, X=X, Z=Z, theta=c(fitMLE$coefficients$mean[,1],fitMLE$coefficients$precision[,1]),
               linkmu="logit", linkphi="log", LSMLE=F,
              main.title = " ", faixa.fixed = c(-6,10), labels.fixed =1:77)

envelope_LSMLE(y=y, X=X, Z=Z, theta=c(fit_LSMLE$beta,fit_LSMLE$gama), linkmu="logit", 
               linkphi="log", LSMLE=T,
               main.title = " ", faixa.fixed = c(-6,10), labels.fixed =1:77)

#zoomed versions

envelope_LSMLE(y=y, X=X, Z=Z, theta=c(fitMLE$coefficients$mean[,1],fitMLE$coefficients$precision[,1]),
               linkmu="logit", linkphi="log", LSMLE=F,
               main.title = "", faixa.fixed = c(-4,3), labels.fixed =1:77)

envelope_LSMLE(y=y, X=X, Z=Z, theta=c(fit_LSMLE$beta,fit_LSMLE$gama), linkmu="logit", 
               linkphi="log", LSMLE=T,
               main.title = "", faixa.fixed = c(-4,3), labels.fixed =1:77)
