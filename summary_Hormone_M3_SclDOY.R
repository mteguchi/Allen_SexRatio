#summary_Hormone_M2


# Tomo Eguchi
# 20 February 2015

rm(list=ls())
library(coda)
library(ggplot2)
library(grid)

# from http://wresch.github.io/2012/11/02/subplots-with-ggplot2.html
subplot <- function(r, c) {
  # select viewport from layout
  #   used to do subplots
  # Args:
  #   r:  row
  #   c:  column
  # Returns:
  #   viewport for plotting
  
  viewport(layout.pos.col=c, layout.pos.row=r)
}

# from http://wresch.github.io/2012/11/02/subplots-with-ggplot2.html
vplayout <- function(r, c) {
  # Set up grid layout for creating subplots
  # Args:
  #   r:  row
  #   c:  column
  # Returns:
  #   viewport
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(r, c)))
}

sysInfo <- Sys.info()
if (sysInfo[1] == 'Linux') {
  source('~/Documents/R/TomosFunctions.R') 
} else if (sysInfo[1] == 'Windows'){
  if (sysInfo[4] == 'SWC-TEGUCHI-D'){ 
    source('~/R/TomosFunctions.R')
  } else if (sysInfo[4] == "SWC-TEGUCHI-L"){
    source("C:/Users/t_e/Documents/R/TomosFunctions.R")
  }
}
#source('~/R/R_Work/TomosFunctions.R')
D <- dirSelector()
RdataFname <- paste(D$Rdir, 
                   "Allen_SexDetermination/RData/Hormone_2015-04-22M3_SclDoyWithSex.RData", 
                   sep = "")
load(RdataFname)

outfilename <- paste0(D$Rdir, 
                      "Allen_SexDetermination/out_M3_SclDOY_",
                      runDate, ".csv") 

# posterior simulations for model check:
beta0 <- unlist(zm[,varnames(zm) == 'beta0'])
betaSex <- unlist(zm[,varnames(zm) == 'beta_Sex'])
betaDOY1 <- unlist(zm[,varnames(zm) == 'beta_DOY1'])
betaDOY2 <- unlist(zm[,varnames(zm) == 'beta_DOY2'])
betaSCL <- unlist(zm[,varnames(zm) == 'beta_SCL'])
#betaTail <- unlist(zm[,varnames(zm) == 'beta_Tail'])
s_H1 <- unlist(zm[,varnames(zm) == 's_H1'])
#mu_tail <- unlist(zm[,varnames(zm) == "mu_tail"])

pSex <- vector(mode = "numeric", length = length(uniqIDs))
for (i in 1:length(pSex)){
  sex1 <- unlist(zm[,varnames(zm) == paste0('Sex1[', i, ']')])
  pSex[i] <- sum(sex1)/length(sex1)
}

Sex2 <- pSex
Sex2[pSex > 0.5] <- 1
Sex2[pSex < 0.5] <- 0

H1a <- matrix(data = NA, nrow = dim(datSDB)[1],
                ncol = 5)
pMale <- vector(mode = "numeric", length = dim(datSDB)[1])
for (i in 1:dim(datSDB)[1]){
  tmpSex <- Sex2[uniqIDs == datSDB$ID[i]]
  pMale[i] <- pSex[uniqIDs == datSDB$ID[i]]
  ifelse(is.na(datSDB$Tail_mean0[i]),
         mu_H1 <- beta0 + betaSex * tmpSex + 
           (betaDOY1) * (datSDB$DOY_mean0[i]) + 
           (betaSCL) * (datSDB$SCL_mean0[i]),
         mu_H1 <- beta0 + betaSex * tmpSex + 
           (betaDOY1) * (datSDB$DOY_mean0[i]) + 
           (betaSCL) * (datSDB$SCL_mean0[i]))
  H1a[i, 1:3] <- quantile(rnorm(length(mu_H1), mu_H1, s_H1),
                         prob = c(0.025, 0.5, 0.975))
  H1a[i, 4:5] <- c(datSDB$logTesto1_mean0[i], 
                  datSDB$logTesto2_mean0[i])
}

# plot the results to see if data fall within predicted 95% PI
df1 <- data.frame(low = H1a[,1], high = H1a[,3], med = H1a[,2],
                  obs1 = H1a[,4], obs2 = H1a[,5], 
                  ID = as.factor(datSDB$ID))
df1$ID2 <- factor(c(1:dim(df1)[1]), ordered = TRUE)
p1 <- ggplot(data = df1, aes(y = med, x = ID)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin=low, ymax=high), width=.1) +
  xlab("ID") + 
  ylab("log Testosterone (pg/ml)") +
  theme(axis.text.x = element_text(size=12, angle=90, hjust=1),
        axis.text.y = element_text(size=12)) + 
  geom_point(aes(y = obs1, x = ID), size = 2, shape = 4) + 
  geom_point(aes(y = obs2, x = ID), size = 2, shape = 4)

p1

write.table(df1, file = paste0("posteriorSimulations_", runDate, ".csv"),
            quote = FALSE, sep = ",", row.names = FALSE)
# Looks ok. 

# marginal quantiles are not useful beause of correlation
# between the two - even after the centering. 
# for ellipsoid, use the following
library(MASS)
library(cluster)
samples <- as.matrix(cbind(beta0, betaSex))
fit <- cov.mve(samples, quantile.used = nrow(samples) * 0.9)
points_in_ellipse <- samples[fit$best,]
ellipse_boundary <- predict(ellipsoidhull(points_in_ellipse))
plot(beta0, betaSex)

lines(ellipse_boundary, col = "lightgreen", lwd = 3)
##############################

all_df <- data.frame(ID = as.factor(datSDB$ID), pMale = pMale,
                     #median_H = H1[,2] + mean(log(c(datSDB$Testo1, datSDB$Testo2))), 
                     Sex = datSDB$Sex, SCL = datSDB$SCL,
                     DOY = datSDB$DOY, Tail = datSDB$Tail,
                     #logTesto1 = datSDB$logTesto1, 
                     #logTesto2 = datSDB$logTesto2,
                     #median_T = exp(H1[,2] + mean(log(c(datSDB$Testo1, datSDB$Testo2)))),
                     obsMeanTesto = datSDB$TestoMean)


write.table(all_df, file = outfilename,
            quote = FALSE, sep = ",", row.names = FALSE)

# p0 <- ggplot(data = all_df, 
#              aes(x = (logTesto1+logTesto2)/2, y = median_H)) +
#   geom_point(size = 3) +
#   xlab("Estimated median testosterone") + 
#   ylab("Mean Log testosterone") + 
#   theme(axis.text.x = element_text(size=12),
#         axis.text.y = element_text(size=12))
# p0

# plot the theoretical mean vs observed:
breaks <- log(c(2,10,100,1000,10000,100000))
labels <- as.character(format(exp(breaks), 
                              digits = 2))
p2 <- ggplot() + 
  geom_point(data = all_df[all_df$Sex != "U",], 
             aes(y = pMale, 
                 x = log(obsMeanTesto)), 
             size = 4) +
  geom_point(data = all_df[all_df$Sex == "U",], 
             aes(y = pMale, 
                 x = log(obsMeanTesto)), 
             size = 4, color = "black", shape = 1) +
  ylab("Pr(male)") + 
  xlab("Observed mean testosterone (pg/ml)") +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  scale_x_continuous(breaks = breaks,
                     labels = labels)
p2

breaks2 <- log(c(2,10,100,1000,10000))
labels <- as.character(format(exp(breaks2), 
                              digits = 2))
p3 <- ggplot(data = all_df[datSDB$Sex == "U",], 
             aes(y = pMale, 
                 x = log(obsMeanTesto))) +
  geom_point(size = 4) +
  ylab("Pr(male)") + 
  xlab("Testosterone (pg/ml)") + 
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  scale_x_continuous(breaks = breaks2,
                     labels = labels)
p3

breaks4 <- log(c(2,10,100,1000,10000, 100000))
labels <- as.character(format(exp(breaks4), 
                              digits = 2))
DF_p4 <- datSDB[order(datSDB$Testo1),] #all_df[order(all_df$mean_H),] 
DF_p4$ID2 <- factor(c(1:dim(DF_p4)[1]), ordered = TRUE)
p4 <- ggplot(data = DF_p4, 
             aes(x = ID2, y = log(TestoMean), fill = Sex)) + # change log(Testo1) to mean_H
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c("white", "black", "gray")) + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=12)) + 
  ylab("Testosterone (pg/ml)") + 
  xlab("") +
  scale_y_continuous(breaks = breaks4,
                     labels = labels)
p4

DF_p5 <- data.frame(beta0 = beta0, betaSex = betaSex,
                    betaSCL = betaSCL,
                    betaDOY = betaDOY1, s_H1 = s_H1)
# for beta0, betaSex
p5_1 <- ggplot(data = DF_p5) + 
  geom_bar(aes(x = beta0), binwidth = 0.02) + 
  xlab((expression(beta [0]))) + 
  ylab("Count") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))


p5_2 <- ggplot(data = DF_p5) + 
  geom_bar(aes(x = betaSex), binwidth = 0.05) + 
  xlab((expression(beta [Sex]))) + 
  ylab("Count") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))


# for tail and SCL
p5_3 <- ggplot(data = DF_p5) + 
  geom_bar(aes(x = betaSCL), binwidth = 0.02) + 
  xlab((expression(beta [SCL]))) + 
  ylab("Count") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))


# for DOY
p5_4 <- ggplot(data = DF_p5) + 
  geom_bar(aes(x = betaDOY1), binwidth = 0.02)+
  xlab((expression(beta [DOY]))) + 
  ylab("Count") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12))


plots <- c(p5_1, p5_2, p5_3, p5_4)


vplayout(2, 2)
plot(p5_1, vp = subplot(1, 1))
plot(p5_2, vp = subplot(1, 2))
plot(p5_3, vp = subplot(2, 1))
plot(p5_4, vp = subplot(2, 2))

pMale <- vector(mode = "numeric", length=length(uniqIDs))
for (i in 1:length(uniqIDs)){
  tmp <- all_df[all_df$ID == uniqIDs[i],]
  pMale[i] <- max(tmp$pMale)
  
}

pMale_U <- pMale[all_df$Sex == "U"]
length(pMale_U[pMale_U > 0.5])


p6 <- ggplot(data = all_df, 
             aes(y = (obsMeanTesto), 
                 x = DOY)) +
  geom_point(size = 4) +
  ylab("T (pg/ml))") + 
  xlab("Day of year") + 
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))
p6
