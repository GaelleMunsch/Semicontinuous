#!/bin/bash

module load R/3.6.1

R --vanilla << END

################ Load packages
##############################

library(MASS) # for Negative binomial model
library(cplm) # for Compound Poisson-Gamma model
library(tweedie) # for Compound Poisson-Gamma model parameters estimation
library(VGAM) # for Tobit model

DATA <- read.delim("Phenotypic_file.txt")

DATA[,"NET1000"] <- DATA[,"NET"]*1000



######## Model estimations

model <- glm.nb(NET1000 ~ Age + Sex + tabac + status, data = DATA)
summary(model)

model_cpg <- cpglm(NET1000 ~ Age + Sex + tabac + status, data=DATA)
summary(model_cpg)

model_tobit <- vglm(NET1000 ~ Age + Sex + tabac + status, tobit(Lower = 0), data=DATA)
summary(model_tobit)



#-------------------------------------------------------------------#
######## Simulations from estimations -- Compound Poisson-Gamma model

CPG_ALL1000 <- data.frame(mat = matrix(ncol = 1000, nrow = 657))

for (simu in 1:1000) {
  set.seed(simu)
  CPG_ALL1000[,simu] <- rtweedie(n=dim(DATA)[1], mu=fitted(model_cpg), phi=as.numeric(model_cpg["phi"]), power=as.numeric(model_cpg["p"]))
  colnames(CPG_ALL1000)[simu] <- paste0("SIMU_",simu)
}

CPG_ALL1000 <- as.data.frame(cbind(DATA,CPG_ALL1000))

## Calculate RMSE
res <- c()
for (i in 1:1000){
res[i] <- sqrt(mean((CPG_ALL1000[,i+11] - CPG_ALL1000$NET1000)^2))
}

## print mean RMSE and 95% confidence interval
mean(res)
mean(res)-1.96*sd(res)
mean(res)+1.96*sd(res)
#-------------------------------------------------------------------#


#-------------------------------------------------------------------#
######## Simulations from estimations -- Negatif Binomial model

NB_ALL1000 <- data.frame(mat = matrix(ncol = 1000, nrow = 657))

for (simu in 1:1000) {
  set.seed(simu)
  NB_ALL1000[,simu] <- rnegbin(fitted(model), theta=model$theta)
  colnames(NB_ALL1000)[simu] <- paste0("SIMU_",simu)
}

NB_ALL1000 <- as.data.frame(cbind(DATA,NB_ALL1000))

## Calculate RMSE
res <- c()
for (i in 1:1000){
res[i] <- sqrt(mean((NB_ALL1000[,i+11] - NB_ALL1000$NET1000)^2))
}

## print RMSE and 95% confidence interval
mean(res)
mean(res)-1.96*sd(res)
mean(res)+1.96*sd(res)
#-------------------------------------------------------------------#




#-------------------------------------------------------------------#
######## Simulations from estimations -- Tobit model

TOBIT_ALL1000 <- data.frame(mat = matrix(ncol = 1000, nrow = 657))

for (simu in 1:1000) {
  set.seed(simu)
  TOBIT_ALL1000[,simu] <- rnorm(n=dim(DATA)[1], mean=model_tobit@predictors[,1], sd=exp(model_tobit@predictors[,2]))
  # TOBIT_ALL1000[which(TOBIT_ALL1000[,simu] < 0 )] <- 0  ##### possible to truncate simulations at zero
  colnames(TOBIT_ALL1000)[simu] <- paste0("SIMU_",simu)
}

TOBIT_ALL1000 <- as.data.frame(cbind(DATA,TOBIT_ALL1000))

## Calculate RMSE
res <- c()
for (i in 1:1000){
res[i] <- sqrt(mean((TOBIT_ALL1000[,i+11] - TOBIT_ALL1000$NET1000)^2))
}

## print RMSE and 95% confidence interval
mean(res)
mean(res)-1.96*sd(res)
mean(res)+1.96*sd(res)
#-------------------------------------------------------------------#



END