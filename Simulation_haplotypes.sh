################ Load packages
##############################
library (plyr)
library(data.table)

DATA <- read.delim("Phenotypic_file.txt")

DATA[,"NET1000"] <- DATA[,"NET"]*1000

#-----------------------------------------------------------------------#
###################### Simulation of the datasets #######################

## See Simulation.sh to create this file
DTBY10000 <- read.delim("DTBY10000.txt")


#--------------------------------------------------------------------------------------#
###################### Simulation of the haplotypes using 4 SNPs #######################

## Haplotypes frequencies
p1=0.89
p2=0.01
p3=0.05
p4=0.04
p5=0.01


#########################################################################################################
####################### Simulation of 2 databases (1 for each allele) according to haplotypes frequencies

DTBAlleles1 <- data.frame(mat = matrix(ncol = 10000, nrow = 657))
for (simu in 1:10000) {
  set.seed(simu)
  random <- runif(n=657,min=0,max=1)
  DTBAlleles1[,simu] <- ifelse(random <= p1, "H1", ifelse(random <= p1+p2,"H2", ifelse(random <= p1+p2+p3,"H3",ifelse(random <= p1+p2+p3+p4,"H4",ifelse(random <= p1+p2+p3+p4+p5,"H5",NA)))))
  colnames(DTBAlleles1)[simu] <- paste0("Haplo_",simu)
}

DTBAlleles2 <- data.frame(mat = matrix(ncol = 10000, nrow = 657))
for (simu in 1:10000) {
  set.seed(simu*1.5)
  random <- runif(n=657,min=0,max=1)
  DTBAlleles2[,simu] <- ifelse(random <= p1, "H1", ifelse(random <= p1+p2,"H2", ifelse(random <= p1+p2+p3,"H3",ifelse(random <= p1+p2+p3+p4,"H4",ifelse(random <= p1+p2+p3+p4+p5,"H5",NA)))))
  colnames(DTBAlleles2)[simu] <- paste0("Haplo_",simu)
}



####################################################################################################
####################### New database that represents the diplotypes (from DTBAlleles1 & DTBAlleles2)

DTBDiplo <- data.frame(mat = matrix(ncol = 10000, nrow = 657))
for (i in 1:dim(DTBDiplo)[1]){
  for (j in 1:dim(DTBDiplo)[2]){
    DTBDiplo[i,j] <- paste0(DTBAlleles1[i,j],DTBAlleles2[i,j])
  }
}



##################################################################
########## Creation of 1 database per SNP composing the haplotypes
DTBSNP1 <- data.frame(mat = matrix(ncol = 10000, nrow = 657))
DTBSNP2 <- data.frame(mat = matrix(ncol = 10000, nrow = 657))
DTBSNP3 <- data.frame(mat = matrix(ncol = 10000, nrow = 657))
DTBSNP4 <- data.frame(mat = matrix(ncol = 10000, nrow = 657))


######################################################################
########## From simulated diplotypes, we find the values of the 4 SNPs 

SNP1 <- function(DTBSNP, DTBDiplo, SIMU){

  list_SNP_0 <- c("H1H1","H1H2","H2H1","H1H3","H3H1","H2H2","H2H3","H3H2", "H3H3")
  list_SNP_1 <- c("H1H4","H4H1","H5H1","H1H5","H2H4","H4H2","H2H5","H5H2", "H3H4", "H4H3", "H3H5", "H5H3")
  list_SNP_2 <- c("H4H4","H4H5","H5H4","H5H5")
  
  DTBSNP[,SIMU] <- ifelse(DTBDiplo[,SIMU] %in%  list_SNP_0, 0, ifelse(DTBDiplo[,SIMU] %in%  list_SNP_1, 1, ifelse(DTBDiplo[,SIMU] %in%  list_SNP_2,2,NA)))

  return(DTBSNP[,SIMU])
}


DTBSNP1 <- lapply(colnames(DTBDiplo)[1:dim(DTBDiplo)[2]],
                   FUN            = SNP1,
                   DTBDiplo       = DTBDiplo,
                   DTBSNP         = DTBSNP1)
DTBSNP1 <- as.data.frame(DTBSNP1)
colnames(DTBSNP1) <- colnames(DTBDiplo)



SNP2 <- function(DTBSNP, DTBDiplo, SIMU){

  list_SNP_0 <- c("H1H1","H1H2","H2H1","H2H2")
  list_SNP_1 <- c("H2H3","H3H2","H1H3","H3H1","H1H4","H4H1","H5H1","H1H5","H2H4","H4H2","H2H5","H5H2")
  list_SNP_2 <- c("H4H4","H4H5","H5H4","H5H5", "H3H3", "H3H4", "H4H3", "H3H5", "H5H3")
  
  DTBSNP[,SIMU] <- ifelse(DTBDiplo[,SIMU] %in%  list_SNP_0, 0, ifelse(DTBDiplo[,SIMU] %in%  list_SNP_1, 1, ifelse(DTBDiplo[,SIMU] %in%  list_SNP_2,2,NA)))

  return(DTBSNP[,SIMU])
}

DTBSNP2 <- lapply(colnames(DTBDiplo)[1:dim(DTBDiplo)[2]],
                   FUN            = SNP2,
                   DTBDiplo       = DTBDiplo,
                   DTBSNP         = DTBSNP2)

DTBSNP2 <- as.data.frame(DTBSNP2)
colnames(DTBSNP2) <- colnames(DTBDiplo)




SNP3 <- function(DTBSNP, DTBDiplo, SIMU){

  list_SNP_0 <- c("H1H1","H1H3","H3H1","H4H4","H1H4","H4H1", "H3H3", "H3H4", "H4H3")
  list_SNP_1 <- c("H1H2","H2H1","H2H3","H3H2","H5H1","H1H5","H4H5","H5H4","H2H4","H4H2", "H3H5", "H5H3")
  list_SNP_2 <- c("H5H5","H2H2","H2H5","H5H2")
  
  DTBSNP[,SIMU] <- ifelse(DTBDiplo[,SIMU] %in%  list_SNP_0, 0, ifelse(DTBDiplo[,SIMU] %in%  list_SNP_1, 1, ifelse(DTBDiplo[,SIMU] %in%  list_SNP_2,2,NA)))

  return(DTBSNP[,SIMU])
}

DTBSNP3 <- lapply(colnames(DTBDiplo)[1:dim(DTBDiplo)[2]],
                   FUN            = SNP3,
                   DTBDiplo       = DTBDiplo,
                   DTBSNP         = DTBSNP3)

DTBSNP3 <- as.data.frame(DTBSNP3)
colnames(DTBSNP3) <- colnames(DTBDiplo)



SNP4 <- function(DTBSNP, DTBDiplo, SIMU){

  list_SNP_0 <- c("H1H1","H1H3","H3H1", "H3H3")
  list_SNP_1 <- c("H2H3","H3H2","H1H2","H2H1","H1H4","H4H1","H5H1","H1H5", "H3H4", "H4H3", "H3H5", "H5H3")
  list_SNP_2 <- c("H2H5","H5H2","H2H4","H4H2","H2H2","H4H4","H4H5","H5H4","H5H5")
  
  DTBSNP[,SIMU] <- ifelse(DTBDiplo[,SIMU] %in%  list_SNP_0, 0, ifelse(DTBDiplo[,SIMU] %in%  list_SNP_1, 1, ifelse(DTBDiplo[,SIMU] %in%  list_SNP_2,2,NA)))

  return(DTBSNP[,SIMU])
}

DTBSNP4 <- lapply(colnames(DTBDiplo)[1:dim(DTBDiplo)[2]],
                   FUN            = SNP4,
                   DTBDiplo       = DTBDiplo,
                   DTBSNP         = DTBSNP4)

DTBSNP4 <- as.data.frame(DTBSNP4)
colnames(DTBSNP4) <- colnames(DTBDiplo)

###################################################################################################
###################################################################################################
############### Then the same procedure described in Simulation.sh can be applied using those files



END

