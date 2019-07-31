#Hormone_jags


# Tomo Eguchi
# 20 February 2015

rm(list=ls())

library(rjags)
library(coda)

Tbegin <- Sys.time
save <- TRUE
## set MCMC parameters
n.adapt <- 10000
n.update <- 20000
n.iter <- 20000
n.chains <- 5
n.thin <- 1

sysInfo <- Sys.info()
if (sysInfo[1] == 'Linux') {
  source('~/Documents/R/TomosFunctions.R') 
} else if (sysInfo[1] == 'Windows'){
  if (sysInfo[4] == 'SWC-TEGUCHI-D'){ 
    source('~/R/TomosFunctions.R')
  }
}
#source('~/R/R_Work/TomosFunctions.R')
D <- dirSelector()
runDate <- Sys.Date()

## parameters to monitor - when this is changed, make sure to change
## summary statistics indices at the end of this script. 
parameters <- c("beta0", "beta_Sex",  
                "s_H1", "deviance",
                "Sex1", "Sex2") 

if (save == TRUE) {
  if (length(grep("Sex1", parameters)) > 0){
    saveFname <- paste(D$Rdir, 
                       "Allen_SexDetermination/RData/Hormone_", 
                       runDate, "M3_SexWithSex.RData", sep = "")
  } else {
    saveFname <- paste(D$Rdir, 
                       "Allen_SexDetermination/RData/Hormone_", 
                       runDate, "M3_SexWithoutSex.RData", sep = "") 
  }
}  



# define model and data file names
modelName <- "Model_Hormone_3_SexOnly.txt"

fileSDB <- paste(D$Rdir, 
                 "Allen_SexDetermination/Testo_SDB5.csv", sep = "")
fileBDAPAN <- paste(D$Rdir, 
                    "Allen_SexDetermination/Testo_BDAPAN2.csv", sep = "")

## get data files:
datSDB <- read.table(file = fileSDB, header = TRUE, sep = ",")
datBDAPAN <- read.table(file = fileBDAPAN, header = TRUE, sep = ",")

# Compute DOY for each.
datSDB$DOY <- as.numeric(strftime(paste(datSDB$YR, '-', 
                                        datSDB$MO, '-', 
                                        datSDB$DA, sep = ""), 
                                  format = "%j"))
datBDAPAN$DOY <- as.numeric(strftime(paste(datBDAPAN$YR, '-', 
                                           datBDAPAN$MO, '-', 
                                           datBDAPAN$DA, sep = ""), 
                                     format = "%j"))

# make 0-1 sex column
datSDB$Sex01 <- rep(0, times = dim(datSDB)[1])
datSDB$Sex01[datSDB$Sex == 'M'] <- 1
datSDB$Sex01[datSDB$Sex == 'U'] <- NA

datBDAPAN$Sex01 <- rep(0, times = dim(datBDAPAN)[1])
datBDAPAN$Sex01[datBDAPAN$Sex == 'M'] <- 1
datBDAPAN$Sex01[datBDAPAN$Sex == 'U'] <- NA

# log then de-mean per lab
datSDB$logTesto1_mean0 <- (log(datSDB$Testo1) - 
                          mean(log(c(datSDB$Testo1, datSDB$Testo2))))/(sd(log(c(datSDB$Testo1, datSDB$Testo2))))
datSDB$logTesto2_mean0 <- (log(datSDB$Testo2) - 
                          mean(log(c(datSDB$Testo1, datSDB$Testo2))))/(sd(log(c(datSDB$Testo1, datSDB$Testo2))))

datBDAPAN$logTesto1_mean0 <- (log(datBDAPAN$Testo1) - 
                             mean(log(c(datBDAPAN$Testo1, datBDAPAN$Testo2))))/(sd(log(c(datBDAPAN$Testo1, datBDAPAN$Testo2))))
datBDAPAN$logTesto2_mean0 <- (log(datBDAPAN$Testo2) - 
                             mean(log(c(datBDAPAN$Testo1, datBDAPAN$Testo2))))/(sd(log(c(datBDAPAN$Testo1, datBDAPAN$Testo2))))


# data for SDB
datSDB$SCL_mean0 <- (datSDB$SCL - mean(datSDB$SCL))/sd(datSDB$SCL)
datSDB$DOY_mean0 <- (datSDB$DOY - mean(datSDB$DOY))/sd(datSDB$DOY)
datSDB$Tail_mean0 <- (datSDB$Tail - mean(datSDB$Tail, na.rm=T))/sd(datSDB$Tail, na.rm=T)

uniqIDs <- unique(datSDB$ID)
H1 <- array(data = NA, dim = c(length(uniqIDs), 
                               max(datSDB$Rep), 2))
Tail1 <- DOY1 <- SCL1 <- matrix(data=NA, nrow = length(uniqIDs), ncol = max(datSDB$Rep))
n1 <- Sex1 <- vector(mode = 'numeric', length = length(uniqIDs))

for (j in 1:length(uniqIDs)){
  tmp <- datSDB[datSDB$ID == uniqIDs[j],]
  for (k in 1:max(tmp$Rep)){
    H1[j,k,1] <- tmp[tmp$Rep == k, "logTesto1_mean0"]
    H1[j,k,2] <- tmp[tmp$Rep == k, "logTesto2_mean0"]
    SCL1[j,k] <- tmp[tmp$Rep == k, "SCL_mean0"]
    DOY1[j,k] <- tmp[tmp$Rep == k, "DOY_mean0"]
    Tail1[j,k] <- tmp[tmp$Rep == k, "Tail_mean0"]
  }
  Sex1[j] <- tmp$Sex01[1]
  n1[j] <- max(tmp$Rep)
}

# data for BDAPAN - no recaptures
datBDAPAN$SCL_mean0 <- (datBDAPAN$SCL - mean(datBDAPAN$SCL))/sd(datBDAPAN$SCL)
datBDAPAN$DOY_mean0 <- (datBDAPAN$DOY - mean(datBDAPAN$DOY))/sd(datBDAPAN$DOY)
datBDAPAN$Tail_mean0 <- (datBDAPAN$Tail - mean(datBDAPAN$Tail, na.rm=T))/sd(datBDAPAN$Tail, na.rm=T)
uniqIDsBP <- unique(datBDAPAN$ID)
H2 <- matrix(data = NA, nrow = length(uniqIDsBP), 
             ncol = 2)
Tail2 <- DOY2 <- SCL2 <- n2 <- Sex2 <- vector(mode = "numeric", 
                                              length = length(uniqIDsBP))

for (j in 1:length(uniqIDsBP)){
  tmp <- datBDAPAN[datBDAPAN$ID == uniqIDsBP[j],]
  H2[j,1] <- tmp$logTesto1_mean0
  H2[j,2] <- tmp$logTesto2_mean0
  SCL2[j] <- tmp$SCL_mean0
  DOY2[j] <- tmp$DOY_mean0
  Tail2[j] <- tmp$Tail_mean0
  Sex2[j] <- tmp$Sex01
}

bugs.data <- list(H1 = H1,
                  H2 = H2,
                  Sex1 = Sex1,
                  Sex2 = Sex2,
                  N1 = length(uniqIDs), 
                  N2 = length(uniqIDsBP),
                  n1 = n1)

# initial values
# input n is the number of unknown sex
initsFunction <- function(ns){
  #mu_H1 <- matrix(rnorm(ns[1]*ns[2], 0, 10), 
  #                nrow = ns[1], ncol = ns[2])
  #mu_H2 <- rnorm(ns[3], 0, 10)
  
  tau_H1 <- rgamma(0.1, 0.01)
  
  beta0 <- rnorm(1, 0, 10)
  beta_Sex <- rnorm(1, 0, 10)
  
  A <- list(tau_H1 = tau_H1, beta0 = beta0, beta_Sex = beta_Sex)
  
  return(A)
}

# Create a list of length nChains with z in initsFunction
#inList <- list(CH, f)
inits <- lapply(replicate(n.chains, 
                          c(length(uniqIDs), max(datSDB$Rep),
                            length(uniqIDsBP)), 
                          simplify = FALSE), 
                initsFunction)

# use the following line without z in initsFunction
#inits <- apply(as.matrix(1:n.chains), 1, initsFunction)

jm <- jags.model(modelName, 
                 data = bugs.data, 
                 inits, 
                 n.chains = n.chains, 
                 n.adapt = n.adapt)

if (save == TRUE) save(list = ls(all = TRUE), file = saveFname)

parallel.seeds("base::BaseRNG", length(inits))
update(jm, n.iter = n.update)
#update(jm, n.iter = n.iter)
load.module("dic")

# ##############################################
# Using coda.samples - 
zm <- coda.samples(jm,
                   variable.names = parameters, 
                   n.iter = n.iter,
                   thin = n.thin)
if (save == TRUE) save(list = ls(all = TRUE), file = saveFname)

dicOut <- dic.samples(jm, 
                      n.iter=n.iter,
                      type = "pD")
if (save == TRUE) save(list = ls(all = TRUE), file = saveFname)

# use raftery.diag to figure out the required chain lengths:
if (length(grep("Sex1", parameters)) < 1){
  # use raftery.diag to figure out the required chain lengths:
  r.diag <- raftery.diag(zm)
  
  # test also with other convergence diagnostic tools:
  g.diag <- gelman.diag(zm)
  h.diag <- heidel.diag(zm)
}

meanDev <- mean(unlist(zm[, varnames(zm) == 'deviance']))
varDev  <- var(unlist(zm[, varnames(zm) == 'deviance']))
# uses Gelman's approximation because monitoring pD through
# jags.samples ended up in a lot of Inf's... 
DIC <- meanDev + 0.5 * varDev  # Gelman's approx. 

if (save == TRUE) save(list = ls(all = TRUE), file = saveFname)

