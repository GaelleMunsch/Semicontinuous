#!/bin/bash

module load R/3.6.1

R --vanilla << END



################ Load packages
##############################

library (plyr)
library(cplm)
library(VGAM)
library(MASS)
library(tweedie)



DATA <- read.delim("Phenotypic_file.txt")

DATA[,"NET1000"] <- DATA[,"NET"]*1000


#-----------------------------------------------------------------------#
###################### Simulation of the datasets #######################


###################### Creation of 10 000 bootstraps from NETs values #######################

DTBY10000 <- data.frame(mat = matrix(ncol = 10000, nrow = 657))

for (simu in 1:10000) {
  set.seed(simu)
  
  DTBY10000[,simu] <- sample(DATA[,"NET1000"], size=dim(DATA)[1], replace=T)
  colnames(DTBY10000)[simu] <- paste0("SIMU_",simu)
}


###################### Simulation of 4 genetic variants for each 10 000 bootstraps #######################

##############################
### Effect Allele frequency 1%
p=0.01
q=1-p

DTBSNP01 <- data.frame(mat = matrix(ncol = 10000, nrow = 657))

for (simu in 1:10000) {
  set.seed(simu)
  
  random <- runif(n=657,min=0,max=1)
  DTBSNP01[,simu] <- ifelse(random <= p*p, 2, ifelse(random <= 2*q*p,1, 0))
  colnames(DTBSNP01)[simu] <- paste0("SNP01_",simu)
}


##############################
### Effect Allele frequency 5%
p=0.05
q=1-p

DTBSNP05 <- data.frame(mat = matrix(ncol = 10000, nrow = 657))

for (simu in 1:10000) {
  set.seed(simu*2)
  
  random <- runif(n=657,min=0,max=1)
  DTBSNP05[,simu] <- ifelse(random <= p*p, 2, ifelse(random <= 2*q*p,1, 0))
  colnames(DTBSNP05)[simu] <- paste0("SNP05_",simu)
}


###############################
### Effect Allele frequency 10%
p=0.1
q=1-p

DTBSNP10 <- data.frame(mat = matrix(ncol = 10000, nrow = 657))

for (simu in 1:10000) {
  set.seed(simu*3)
  
  random <- runif(n=657,min=0,max=1)
  DTBSNP10[,simu] <- ifelse(random <= p*p, 2, ifelse(random <= 2*q*p,1, 0))
  colnames(DTBSNP10)[simu] <- paste0("SNP10_",simu)
}


###############################
### Effect Allele frequency 20%
p=0.2
q=1-p

DTBSNP20 <- data.frame(mat = matrix(ncol = 10000, nrow = 657))

for (simu in 1:10000) {
  set.seed(simu*4)
  
  random <- runif(n=657,min=0,max=1)
  DTBSNP20[,simu] <- ifelse(random <= p*p, 2, ifelse(random <= 2*q*p,1, 0))
  colnames(DTBSNP20)[simu] <- paste0("SNP20_",simu)
}


#-----------------------------------------------------------------------#











#-----------------------------------------------------------------------#
#################### Analyses with Compound Poisson-Gamma model #########


############################## FONCTION pour le CPG
CPGLM <- function(Y,DTBSNP01,DTBSNP05,DTBSNP10,DTBSNP20,data){
  print(which(colnames(data) == Y))
  SNP01 <- DTBSNP01[,which(colnames(data) == Y)]
  SNP05 <- DTBSNP05[,which(colnames(data) == Y)]
  SNP10 <- DTBSNP10[,which(colnames(data) == Y)]
  SNP20 <- DTBSNP20[,which(colnames(data) == Y)]
  
  DTB <- as.data.frame(cbind(cbind(cbind(cbind(data[,Y],SNP01),SNP05),SNP10),SNP20))
  
  assign("last.warning", NULL, envir = baseenv())

  #### Compound Poisson-Gamma model
  ##########################################
  Try_Error <- try(model <- cpglm(DTB[,"V1"] ~ DTB[,"SNP01"] + DTB[,"SNP05"] + DTB[,"SNP10"] + DTB[,"SNP20"]))
 
 
  if(length(warnings()) > 0){
    print(warnings())
          res <- list(Y,FALSE,NA,NA,NA,NA,NA,NA,NA,NA)
          names(res) <- c("SIMU","CONV","BETA_SNP01", "SE_SNP01","BETA_SNP05", "SE_SNP05","BETA_SNP10", "SE_SNP10","BETA_SNP20", "SE_SNP20")  
          return(res)
  }
  assign("last.warning", NULL, envir = baseenv())
  
  if(class(Try_Error)[1] == "try-error"){
          res <- list(Y,FALSE,NA,NA,NA,NA,NA,NA,NA,NA)
          names(res) <- c("SIMU","CONV","BETA_SNP01", "SE_SNP01","BETA_SNP05", "SE_SNP05","BETA_SNP10", "SE_SNP10","BETA_SNP20", "SE_SNP20")  
          return(res)
  }

  res <- list(Y, model[["converged"]], summary(model)[["coefficients"]][2,1], summary(model)[["coefficients"]][2,2], summary(model)[["coefficients"]][3,1], summary(model)[["coefficients"]][3,2], summary(model)[["coefficients"]][4,1], summary(model)[["coefficients"]][4,2], summary(model)[["coefficients"]][5,1], summary(model)[["coefficients"]][5,2])
    
  names(res) <- c("SIMU","CONV","BETA_SNP01", "SE_SNP01","BETA_SNP05", "SE_SNP05","BETA_SNP10", "SE_SNP10","BETA_SNP20", "SE_SNP20")  

  return(res)
}


CPG_YBOOT <- lapply(colnames(DTBY10000)[1:dim(DTBY10000)[2]],
                   FUN              = CPGLM,
                   data             = DTBY10000,
                   DTBSNP01         = DTBSNP01,
                   DTBSNP05         = DTBSNP05,
                   DTBSNP10         = DTBSNP10,
                   DTBSNP20         = DTBSNP20)

CPG_YBOOT <- ldply(CPG_YBOOT, data.frame)







######## Analyses of the results
RES <- CPG_YBOOT
table(RES$CONV)

############## for the genetic variant with an effect allele frequency at 1%
RES$STAT <- RES$BETA_SNP01 / RES$SE_SNP01

########### Number of significant statistical tests
#### alpha 0.1 %
dim(RES[which(RES$STAT < -3.291 | RES$STAT > 3.291),])[1]
#### alpha 0.5 %
dim(RES[which(RES$STAT < -3.09 | RES$STAT > 3.09),])[1]
#### alpha 1 %
dim(RES[which(RES$STAT < -2.576 | RES$STAT > 2.576),])[1]
#### alpha 5 %
dim(RES[which(RES$STAT < -1.96 | RES$STAT > 1.96),])[1]
#### alpha 10 %
dim(RES[which(RES$STAT < -1.645 | RES$STAT > 1.645),])[1]



############## for the genetic variant with an effect allele frequency at 5%
RES$STAT <- RES$BETA_SNP05 / RES$SE_SNP05

########### Number of significant statistical tests
#### alpha 0.1 %
dim(RES[which(RES$STAT < -3.291 | RES$STAT > 3.291),])[1]
#### alpha 0.5 %
dim(RES[which(RES$STAT < -3.09 | RES$STAT > 3.09),])[1]
#### alpha 1 %
dim(RES[which(RES$STAT < -2.576 | RES$STAT > 2.576),])[1]
#### alpha 5 %
dim(RES[which(RES$STAT < -1.96 | RES$STAT > 1.96),])[1]
#### alpha 10 %
dim(RES[which(RES$STAT < -1.645 | RES$STAT > 1.645),])[1]



############## for the genetic variant with an effect allele frequency at 10%
RES$STAT <- RES$BETA_SNP10 / RES$SE_SNP10

########### Number of significant statistical tests
#### alpha 0.1 %
dim(RES[which(RES$STAT < -3.291 | RES$STAT > 3.291),])[1]
#### alpha 0.5 %
dim(RES[which(RES$STAT < -3.09 | RES$STAT > 3.09),])[1]
#### alpha 1 %
dim(RES[which(RES$STAT < -2.576 | RES$STAT > 2.576),])[1]
#### alpha 5 %
dim(RES[which(RES$STAT < -1.96 | RES$STAT > 1.96),])[1]
#### alpha 10 %
dim(RES[which(RES$STAT < -1.645 | RES$STAT > 1.645),])[1]


############## for the genetic variant with an effect allele frequency at 20%
RES$STAT <- RES$BETA_SNP20 / RES$SE_SNP20

########### Number of significant statistical tests
#### alpha 0.1 %
dim(RES[which(RES$STAT < -3.291 | RES$STAT > 3.291),])[1]
#### alpha 0.5 %
dim(RES[which(RES$STAT < -3.09 | RES$STAT > 3.09),])[1]
#### alpha 1 %
dim(RES[which(RES$STAT < -2.576 | RES$STAT > 2.576),])[1]
#### alpha 5 %
dim(RES[which(RES$STAT < -1.96 | RES$STAT > 1.96),])[1]
#### alpha 10 %
dim(RES[which(RES$STAT < -1.645 | RES$STAT > 1.645),])[1]

#-----------------------------------------------------------------------#








#-----------------------------------------------------------------------#
#################### Analyses with Negatif Binomial model #########

NBGLM <- function(Y,DTBSNP01,DTBSNP05,DTBSNP10,DTBSNP20,data){
  print(which(colnames(data) == Y))
  
  SNP01 <- DTBSNP01[,which(colnames(data) == Y)]
  SNP05 <- DTBSNP05[,which(colnames(data) == Y)]
  SNP10 <- DTBSNP10[,which(colnames(data) == Y)]
  SNP20 <- DTBSNP20[,which(colnames(data) == Y)]
  
  DTB <- as.data.frame(cbind(cbind(cbind(cbind(data[,Y],SNP01),SNP05),SNP10),SNP20))
  
  assign("last.warning", NULL, envir = baseenv())
  
  #### Negatif Binomial model
  ##########################################
  Try_Error <- try(model <- glm.nb(DTB[,"V1"] ~ DTB[,"SNP01"] + DTB[,"SNP05"] + DTB[,"SNP10"] + DTB[,"SNP20"]))
  
  if(length(warnings()) > 0){
    print(warnings())
          res <- list(Y,FALSE,NA,NA,NA,NA,NA,NA,NA,NA)
          names(res) <- c("SIMU","CONV","BETA_SNP01", "SE_SNP01","BETA_SNP05", "SE_SNP05","BETA_SNP10", "SE_SNP10","BETA_SNP20", "SE_SNP20")  
          return(res)
  }
  assign("last.warning", NULL, envir = baseenv())
  
  if(class(Try_Error)[1] == "try-error"){
          res <- list(Y,FALSE,NA,NA,NA,NA,NA,NA,NA,NA)
          names(res) <- c("SIMU","CONV","BETA_SNP01", "SE_SNP01","BETA_SNP05", "SE_SNP05","BETA_SNP10", "SE_SNP10","BETA_SNP20", "SE_SNP20")  
          return(res)
  }

  res <- list(Y, model[["converged"]], summary(model)[["coefficients"]][2,1], summary(model)[["coefficients"]][2,2], summary(model)[["coefficients"]][3,1], summary(model)[["coefficients"]][3,2], summary(model)[["coefficients"]][4,1], summary(model)[["coefficients"]][4,2], summary(model)[["coefficients"]][5,1], summary(model)[["coefficients"]][5,2])
    
  names(res) <- c("SIMU","CONV","BETA_SNP01", "SE_SNP01","BETA_SNP05", "SE_SNP05","BETA_SNP10", "SE_SNP10","BETA_SNP20", "SE_SNP20")  

  return(res)
}


NB_YBOOT <- lapply(colnames(DTBY10000)[1:dim(DTBY10000)[2]],
                   FUN              = NBGLM,
                   data             = DTBY10000,
                   DTBSNP01         = DTBSNP01,
                   DTBSNP05         = DTBSNP05,
                   DTBSNP10         = DTBSNP10,
                   DTBSNP20         = DTBSNP20)

NB_YBOOT <- ldply(NB_YBOOT, data.frame)



######## Analyses of the results
RES <- NB_YBOOT
table(RES$CONV)

############## for the genetic variant with an effect allele frequency at 1%
RES$STAT <- RES$BETA_SNP01 / RES$SE_SNP01

########### Number of significant statistical tests
#### alpha 0.1 %
dim(RES[which(RES$STAT < -3.291 | RES$STAT > 3.291),])[1]
#### alpha 0.5 %
dim(RES[which(RES$STAT < -3.09 | RES$STAT > 3.09),])[1]
#### alpha 1 %
dim(RES[which(RES$STAT < -2.576 | RES$STAT > 2.576),])[1]
#### alpha 5 %
dim(RES[which(RES$STAT < -1.96 | RES$STAT > 1.96),])[1]
#### alpha 10 %
dim(RES[which(RES$STAT < -1.645 | RES$STAT > 1.645),])[1]



############## for the genetic variant with an effect allele frequency at 5%
RES$STAT <- RES$BETA_SNP05 / RES$SE_SNP05

########### Number of significant statistical tests
#### alpha 0.1 %
dim(RES[which(RES$STAT < -3.291 | RES$STAT > 3.291),])[1]
#### alpha 0.5 %
dim(RES[which(RES$STAT < -3.09 | RES$STAT > 3.09),])[1]
#### alpha 1 %
dim(RES[which(RES$STAT < -2.576 | RES$STAT > 2.576),])[1]
#### alpha 5 %
dim(RES[which(RES$STAT < -1.96 | RES$STAT > 1.96),])[1]
#### alpha 10 %
dim(RES[which(RES$STAT < -1.645 | RES$STAT > 1.645),])[1]



############## for the genetic variant with an effect allele frequency at 10%
RES$STAT <- RES$BETA_SNP10 / RES$SE_SNP10

########### Number of significant statistical tests
#### alpha 0.1 %
dim(RES[which(RES$STAT < -3.291 | RES$STAT > 3.291),])[1]
#### alpha 0.5 %
dim(RES[which(RES$STAT < -3.09 | RES$STAT > 3.09),])[1]
#### alpha 1 %
dim(RES[which(RES$STAT < -2.576 | RES$STAT > 2.576),])[1]
#### alpha 5 %
dim(RES[which(RES$STAT < -1.96 | RES$STAT > 1.96),])[1]
#### alpha 10 %
dim(RES[which(RES$STAT < -1.645 | RES$STAT > 1.645),])[1]


############## for the genetic variant with an effect allele frequency at 20%
RES$STAT <- RES$BETA_SNP20 / RES$SE_SNP20

########### Number of significant statistical tests
#### alpha 0.1 %
dim(RES[which(RES$STAT < -3.291 | RES$STAT > 3.291),])[1]
#### alpha 0.5 %
dim(RES[which(RES$STAT < -3.09 | RES$STAT > 3.09),])[1]
#### alpha 1 %
dim(RES[which(RES$STAT < -2.576 | RES$STAT > 2.576),])[1]
#### alpha 5 %
dim(RES[which(RES$STAT < -1.96 | RES$STAT > 1.96),])[1]
#### alpha 10 %
dim(RES[which(RES$STAT < -1.645 | RES$STAT > 1.645),])[1]

#-----------------------------------------------------------------------#



END