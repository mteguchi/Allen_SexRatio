#Hormone_jags


# Tomo Eguchi
# 20 February 2015

rm(list=ls())

library(rjags)
library(coda)

Tbegin <- Sys.time
save <- TRUE
## set MCMC parameters
n.adapt <- 5000
n.update <- 10000
n.iter <- 10000
n.chains <- 5
#n.thin <- 5

source('~/Documents/R/TomosFunctions.R')
#source('~/R/R_Work/TomosFunctions.R')
D <- dirSelector()
runDate <- Sys.Date()

# define model and data file names
modelName <- "Model_Hormone_1.txt"

fileSDB <- paste(D$Rdir, "Allen_SexDetermination/Testo_SDB2.csv", sep = "")
fileBDAPAN <- paste(D$Rdir, "Allen_SexDetermination/Testo_BDAPAN.csv", sep = "")

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
datSDB$Testo1_mean0 <- log(datSDB$Testo1) - 
  mean(log(c(datSDB$Testo1, datSDB$Testo2)))
datSDB$Testo2_mean0 <- log(datSDB$Testo2) - 
  mean(log(c(datSDB$Testo1, datSDB$Testo2)))

datBDAPAN$Testo1_mean0 <- log(datBDAPAN$Testo1) - 
  mean(log(c(datBDAPAN$Testo1, datBDAPAN$Testo2)))
datBDAPAN$Testo2_mean0 <- log(datBDAPAN$Testo2) - 
  mean(log(c(datBDAPAN$Testo1, datBDAPAN$Testo2)))

# Standardize SCL and DOY
datSDB$SCL1 <- datSDB$SCL - mean(datSDB$SCL)
datBDAPAN$SCL1 <- datBDAPAN$SCL - mean(datBDAPAN$SCL)

if (save == TRUE) saveFname <- paste(D$Rdir, 
                                     "Allen_SexDetermination/RData/Hormone_", 
                                     runDate, ".RData", sep = "");

## parameters to monitor - when this is changed, make sure to change
## summary statistics indices at the end of this script. 
parameters <- c("beta01", "beta02", 
                "beta_H1", "beta_H2",
                "beta_DOY1", "beta_DOY2", 
                "beta_SCL1", "beta_SCL2",
                "mu_H1", "mu_H2", "tau_H1", "tau_H2",
                "Sex1", "Sex2", "p1", "p2")

# data
bugs.data <- list(H1 = cbind(datSDB$Testo1_mean0, datSDB$Testo2_mean0),
                  H2 = cbind(datBDAPAN$Testo1_mean0, datBDAPAN$Testo2_mean0),
                  DOY1 = datSDB$DOY, DOY2 = datBDAPAN$DOY,
                  SCL1 = datSDB$SCL1, SCL2 = datBDAPAN$SCL1,
                  Sex1 = datSDB$Sex01, Sex2 = datBDAPAN$Sex01,
                  N1 = dim(datSDB)[1], N2 = dim(datBDAPAN)[1])

# initial values
# input n is the number of unknown sex
initsFunction <- function(n){
  
  #p1 <- rbeta(n[1], 2, 2)
  #p2 <- rbeta(n[2], 2, 2)
  mu_H1 <- rnorm(n[1], 0, 100)
  mu_H2 <- rnorm(n[2], 0, 100)
  tau_H1 <- rgamma(0.1, 0.01)
  tau_H2 <- rgamma(0.1, 0.01)
  
  beta01 <- rnorm(1, 0, 0.01)
  beta02 <- rnorm(1, 0, 0.01)
  beta_H1 <- rnorm(1, 0, 0.01)
  beta_H2 <- rnorm(1, 0, 0.01)
  beta_DOY1 <- rnorm(1, 0, 0.01)
  beta_DOY2 <- rnorm(1, 0, 0.01)
  beta_SCL1 <- rnorm(1, 0, 0.01)
  beta_SCL2 <- rnorm(1, 0, 0.01)
    
  #Sex1 <- rbinom(n[1], 1, 0.5)
  #Sex2 <- rbinom(n[2], 1, 0.5)
  
  A <- list(mu_H1 = mu_H1, mu_H2 = mu_H2,
            beta01 = beta01, beta02 = beta02,
            beta_H1 = beta_H1, beta_H2 = beta_H2,
            beta_DOY1 = beta_DOY1, beta_DOY2=beta_DOY2,
            beta_SCL1 = beta_SCL1, beta_SCL2 = beta_SCL2)
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
                   n.iter = n.iter)
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


