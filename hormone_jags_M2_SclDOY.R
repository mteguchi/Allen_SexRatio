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

# define model and data file names
modelName <- "Model_Hormone_2_SclDOY.txt"

fileSDB <- paste(D$Rdir, 
                 "Allen_SexDetermination/Testo_SDB3.csv", sep = "")
fileBDAPAN <- paste(D$Rdir, 
                    "Allen_SexDetermination/Testo_BDAPAN.csv", sep = "")

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
datSDB$logTesto1_mean0 <- log(datSDB$Testo1) - 
  mean(log(c(datSDB$Testo1, datSDB$Testo2)))
datSDB$logTesto2_mean0 <- log(datSDB$Testo2) - 
  mean(log(c(datSDB$Testo1, datSDB$Testo2)))

datBDAPAN$logTesto1_mean0 <- log(datBDAPAN$Testo1) - 
  mean(log(c(datBDAPAN$Testo1, datBDAPAN$Testo2)))
datBDAPAN$logTesto2_mean0 <- log(datBDAPAN$Testo2) - 
  mean(log(c(datBDAPAN$Testo1, datBDAPAN$Testo2)))

# or just take a log of each
datSDB$logTesto1 <- log(datSDB$Testo1)
datSDB$logTesto2 <- log(datSDB$Testo2)
datBDAPAN$logTesto1 <- log(datBDAPAN$Testo1)
datBDAPAN$logTesto2 <- log(datBDAPAN$Testo2)

if (save == TRUE) 
  saveFname <- paste(D$Rdir, 
                     "Allen_SexDetermination/RData/Hormone_", 
                     runDate, "M2.RData", sep = "")

## parameters to monitor - when this is changed, make sure to change
## summary statistics indices at the end of this script. 
parameters <- c("beta0", "beta_H", "beta_SCL", "beta_DOY",
                "s_H1", "s_H2", "m1", "m2", "s1", "s2",
                "mu_H1", "mu_H2", "deviance")
                #"Sex1", "Sex2", 

# data
bugs.data <- list(H1 = cbind(datSDB$logTesto1_mean0, 
                             datSDB$logTesto2_mean0),
                  H2 = cbind(datBDAPAN$logTesto1_mean0, 
                             datBDAPAN$logTesto2_mean0),
                  SCL1 = datSDB$SCL-mean(datSDB$SCL), 
                  SCL2 = datBDAPAN$SCL - mean(datBDAPAN$SCL),
                  DOY1 = datSDB$DOY - mean(datSDB$DOY), 
                  DOY2 = datBDAPAN$DOY - mean(datBDAPAN$DOY),
                  Sex1 = datSDB$Sex01, Sex2 = datBDAPAN$Sex01,
                  N1 = dim(datSDB)[1], N2 = dim(datBDAPAN)[1])

# initial values
# input n is the number of unknown sex
initsFunction <- function(n){
  
  #p1 <- rbeta(n[1], 2, 2)
  #p2 <- rbeta(n[2], 2, 2)
  # the parameters of gammas were computed using
  # moment matching for the two datasets
  # a = m^2/v, b = m/v
  mu_H1 <- rgamma(n[1], 3.78, 0.82)
  mu_H2 <- rgamma(n[2], 7.09, 1.24)
  tau_H1 <- rgamma(0.1, 0.01)
  tau_H2 <- rgamma(0.1, 0.01)
  
  beta0 <- rnorm(1, 0, 0.1)
  beta_H <- rnorm(1, 0, 0.1)
  beta_DOY <- rnorm(1, 0, 0.1)
  beta_SCL <- rnorm(1, 0, 0.1)
  
  #Sex1 <- rbinom(n[1], 1, 0.5)
  #Sex2 <- rbinom(n[2], 1, 0.5)
  
  A <- list(mu_H1 = mu_H1, mu_H2 = mu_H2,
            beta0 = beta0, beta_H = beta_H,
            beta_DOY = beta_DOY, beta_SCL = beta_SCL)
  
  return(A)
}

# Create a list of length nChains with z in initsFunction
#inList <- list(CH, f)
inits <- lapply(replicate(n.chains, 
                          c(dim(datSDB)[1], dim(datBDAPAN)[1]), 
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
r.diag <- raftery.diag(zm)

# test also with other convergence diagnostic tools:
g.diag <- gelman.diag(zm)
h.diag <- heidel.diag(zm)

meanDev <- mean(unlist(zm[, varnames(zm) == 'deviance']))
varDev  <- var(unlist(zm[, varnames(zm) == 'deviance']))
# uses Gelman's approximation because monitoring pD through
# jags.samples ended up in a lot of Inf's... 
DIC <- meanDev + 0.5 * varDev  # Gelman's approx. 

if (save == TRUE) save(list = ls(all = TRUE), file = saveFname)

