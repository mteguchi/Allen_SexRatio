#summary_Hormone_M2


# Tomo Eguchi
# 20 February 2015

rm(list=ls())
library(coda)
library(ggplot2)
Mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

gamfit <- function(x){
  a <- (mean(x)^2)/var(x)
  b <- (mean(x))/var(x)
  estim <- list(a = a, b = b)
  return(estim)
}

gamstats <- function(a, b){
  mean <- a/b
  var <- a/(b^2)
  stats <- list(mean = mean, variance = var)
  return(stats)
}
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
RdataFname <- paste(D$Rdir, 
                   "Allen_SexDetermination/RData/Hormone_2015-02-24M2.RData", 
                   sep = "")
load(RdataFname)
beta0 <- unlist(zm[,varnames(zm) == 'beta0'])
betaH <- unlist(zm[,varnames(zm) == 'beta_H'])
a1 <- unlist(zm[,varnames(zm) == 'a1'])
a2 <- unlist(zm[,varnames(zm) == 'a2'])
b1 <- unlist(zm[,varnames(zm) == 'b1'])
b2 <- unlist(zm[,varnames(zm) == 'b2'])

# marginal quantiles are not useful beause of correlation
# between the two - even after the centering. 
# for ellipsoid, use the following
library(MASS)
library(cluster)
samples <- as.matrix(cbind(beta0, betaH))
fit <- cov.mve(samples, quantile.used = nrow(samples) * 0.9)
points_in_ellipse <- samples[fit$best,]
ellipse_boundary <- predict(ellipsoidhull(points_in_ellipse))
plot(beta0, betaH)

lines(ellipse_boundary, col = "lightgreen", lwd = 3)
##############################

# for non-elliptical posterior use the following
#library(emdbook)

# minB0 <- min(ellipse_boundary[,1])
# maxB0 <- max(ellipse_boundary[,1])
# beta0_CI <- c(minB0, maxB0)
# betaH_CI <- c(ellipse_boundary[ellipse_boundary[,1] == minB0, 2], 
#               ellipse_boundary[ellipse_boundary[,1] == maxB0, 2])
# 
# plot(beta0, betaH)
# lines(ellipse_boundary, col = "blue", lwd = 2)
# points(beta0_CI[1], betaH_CI[1], col = "red", pch = 16)
# points(beta0_CI[2], betaH_CI[2], col = "red", pch = 16)

p50 <- vector(length = dim(datSDB)[1], mode = 'numeric')
p05 <- vector(length = dim(datSDB)[1], mode = 'numeric')
p95 <- vector(length = dim(datSDB)[1], mode = 'numeric')
m50 <- vector(length = dim(datSDB)[1], mode = 'numeric')
m05 <- vector(length = dim(datSDB)[1], mode = 'numeric')
m95 <- vector(length = dim(datSDB)[1], mode = 'numeric')
alpha <- vector(length = dim(datSDB)[1], mode = 'numeric')
beta <- vector(length = dim(datSDB)[1], mode = 'numeric')

for (i in 1:dim(datSDB)[1]){
  tmp <- unlist(zm[,varnames(zm) == paste0('mu_H1[', i, ']')])
  #t0 <- gamfit(tmp)
  m50[i] <- quantile(tmp, 0.5) #qgamma(0.5, t0$a, t0$b)
  #mmedian[i] <- median(unlist(zm[,varnames(zm) == paste0('mu_H1[', i, ']')]))
  m05[i] <- quantile(tmp, 0.05) #qgamma(0.05, t0$a, t0$b)
  m95[i] <- quantile(tmp, 0.95) #qgamma(0.95, t0$a, t0$b)
  #t <- exp(mean(beta0) + mean(betaH) * mmean[i])
  t05 <- exp(quantile(beta0 + betaH * m05[i], 0.05))
  t50 <- exp(quantile(beta0 + betaH * m50[i], 0.50))
  t95 <- exp(quantile(beta0 + betaH * m95[i], 0.95))
  
  p50[i] <- t50/(1+t50)
  p05[i] <- t05/(1+t05)
  p95[i] <- t95/(1+t95)
  
  #alpha[i] <- t0$a
  #beta[i] <- t0$b
}

all_df <- data.frame(ID = datSDB$ID, median_H0 = m50,
                     median_H = m50 + mean(log(c(datSDB$Testo1, datSDB$Testo2))),
                     low_p = p05, median_p = p50, high_p = p95,
                     Sex = datSDB$Sex, SCL = datSDB$SCL)

outfilename <- paste0(D$Rdir, 
                      "Allen_SexDetermination/out_2015-02-24.csv") 
write.table(all_df, file = outfilename,
            quote = FALSE, sep = ",", row.names = FALSE)

# plot the theoretical mean vs observed:
breaks <- log(c(2,10,100,1000,10000,100000))
labels <- as.character(format(exp(breaks), 
                              digits = 2))
p1 <- ggplot(data = all_df,
             aes(y = median_p, 
                 x = median_H)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin=low_p, ymax=high_p), width=.1) + 
  ylab("Pr(male)") + 
  xlab("Testosterone (pg/ml)") +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  scale_x_continuous(breaks = breaks,
                     labels = labels)
p1

breaks2 <- log(c(2,10,100,1000,10000))
labels <- as.character(format(exp(breaks2), 
                              digits = 2))
p2 <- ggplot(data = all_df[datSDB$Sex == "U",], 
             aes(y = median_p, 
                 x = median_H)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin=low_p, 
                    ymax=high_p), width=.1) + 
  ylab("Pr(male)") + 
  xlab("Testosterone (pg/ml)") + 
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  scale_x_continuous(breaks = breaks2,
                     labels = labels)
p2

plot(all_df$median_H, (datSDB$logTesto1_mean0+ datSDB$logTesto2_mean0)/2,
     xlab = "Estimated median testosterone",
     ylab = "Observed mean")

p3 <- ggplot(data = all_df, 
             aes(y = median_p, 
                 x = exp(median_H))) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin=low_p, 
                    ymax=high_p), width=.1) + 
  ylab("Pr(male)") + 
  xlab("testosterone")
p3

breaks4 <- log(c(2,10,100,1000,10000, 100000))
labels <- as.character(format(exp(breaks4), 
                              digits = 2))
DF_p4 <- datSDB[order(datSDB$Testo1),] #all_df[order(all_df$mean_H),] 
DF_p4$ID2 <- factor(c(1:dim(DF_p4)[1]), ordered = TRUE)
p4 <- ggplot(data = DF_p4, 
             aes(x = ID2, y = log(Testo1), fill = Sex)) + # change log(Testo1) to mean_H
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("white", "black", "gray")) + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12)) + 
  ylab("Testosterone (pg/ml)") + 
  xlab("") +
  scale_y_continuous(breaks = breaks4,
                     labels = labels)
p4

DF_p5 <- all_df[order(all_df$median_H),]
DF_p5$ID2 <- factor(c(1:dim(DF_p5)[1]), ordered = TRUE)
p5 <- ggplot(data = DF_p5[DF_p5$Sex == "U",], 
             aes(x = ID2, y = median_p)) + 
  geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin=low_p, ymax=high_p), width=.1) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12)) + 
  ylab("Testosterone (pg/ml)") + 
  xlab("")
p5

plot(all_df$median_p, datSDB$Sex01)
plot(all_df$median_H, datSDB$Sex01)

plot(exp(all_df$median_H), all_df$median_p, xlab = 'testosterone',
     ylab = 'Pr(male)', main = '');
points(exp(all_df$median_H[datSDB$Sex == "U"]), 
       all_df$median_p[datSDB$Sex == "U"], pch = 19)


# plot(m0[datSDB$Sex == "U"], p1[datSDB$Sex == "U"])
# points(m0[datSDB$Sex != "U"], p1[datSDB$Sex != "U"], pch=19)

# plot(log(datSDB$Testo1[datSDB$Sex == "U"]), p1[datSDB$Sex == "U"])
# points(log(datSDB$Testo1[datSDB$Sex != "U"]), 
#        p1[datSDB$Sex != "U"], pch = 19)

# get 95% PI for unknown sexes

# mean_b0 <- mean(beta0)
# mean_bH <- mean(betaH)
# idxU <- grep("U", datSDB$Sex)
# 
#   matU <- matrix(data = NA, nrow = length(idxU), ncol = 5)
# for (k in 1:length(idxU)){
#   tmp <- unlist(zm[,varnames(zm) == paste0('mu_H1[', idxU[k], ']')])
#   tmp2 <- mean_b0 + mean_bH * tmp
#   p <- exp(tmp2)/(1 + exp(tmp2))
#   matU[k, ] <- c(mean(p), quantile(p, probs = c(0.05, 0.5, 0.9)), 
#                  datSDB$SCL[idxU[k]])
# }
# 
# US_df <- data.frame(ID = datSDB$ID[idxU],
#                     mean = matU[,1], low = matU[,2],
#                     median = matU[,3], high = matU[,4],
#                     SCL = matU[,5], meanH = exp(m0[idxU]))
# 
