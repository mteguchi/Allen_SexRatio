#Allen testosterone analysis


Tbegin <- Sys.time
sysInfo <- Sys.info()
runDate <- Sys.Date()
runTime <- Sys.time()
# Change the directory structure according to the computer
ifelse(sysInfo[1] == 'Linux',
       source('~/Documents/R/TomosFunctions.R'),
       source('~/R/TomosFunctions.R'))

D00 <- dirSelector()

fileSDB <- paste(D00$Rdir, "Allen_SexDetermination/Testo_SDB2.csv", sep = "")
fileBDAPAN <- paste(D00$Rdir, "Allen_SexDetermination/Testo_BDAPAN.csv", sep = "")
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
                            

# transform testosterone values to standardized to the maximum value of each lab
# then take the log (log(value/maxValue))
datSDB$Testo2X1 <- log(datSDB$Testo2) - log(max(datSDB$Testo2))
datBDAPAN$Testo2X1 <- log(datBDAPAN$Testo2) - log(max(datBDAPAN$Testo2))
datBDAPAN$Testo1X1 <- log(datBDAPAN$Testo1) - log(max(datBDAPAN$Testo1))

# Convert sex into a factor
datSDB$fSex <- factor(datSDB$Sex)
datBDAPAN$fSex <- factor(datBDAPAN$Sex)

# try CART - sex as the response but that doesn't make sense with respect
# to the functional responses; testosterone doesn't determine sex but sex
# determines the testosterone levels. Also, repeated measures are not 
# properly treated.
library(rpart)
fitCart_SDB <- rpart(Testo2X1 ~ Sex + DOY + SCL, data = datSDB)
print(fitCart_SDB)
par(mar = rep(0.1, 4), mai = rep(0.1, 4))
plot(fitCart_SDB)
text(fitCart_SDB)

fitCart_BDAPAN <- rpart(Testo2X1 ~ Sex + DOY + SCL, data = datBDAPAN)
print(fitCart_BDAPAN)
par(mar = rep(2, 4), mai = rep(0.1, 4))
plot(fitCart_BDAPAN)
text(fitCart_BDAPAN)

fitCart_BDAPAN1 <- rpart(Testo1X1 ~ Sex + DOY + SCL, data = datBDAPAN)
print(fitCart_BDAPAN1)
#par(mar = rep(2, 4), mai = rep(0.1, 4))
plot(fitCart_BDAPAN1)
text(fitCart_BDAPAN1)

# combine datasets:

fID_SDB <- as.factor(datSDB$ID)
allData_SD <- data.frame(ID = c(as.character(fID_SDB[datSDB$Sex != "U"]), as.character(datBDAPAN$ID)),
                      Sex =  (c(datSDB$Sex[datSDB$Sex != "U"], (datBDAPAN$Sex)) - 1),
                      Test2X1 = c(datSDB$Testo2X1[datSDB$Sex != "U"], datBDAPAN$Testo2X1),
                      DOY = c(datSDB$DOY[datSDB$Sex != "U"], datBDAPAN$DOY),
                      Tail = c(datSDB$Tail[datSDB$Sex != "U"], datBDAPAN$Tail))

# logistic regression analysis with repeated measures?
library(gtools)   # inv.logit function

# treat all independent.
fit.1 <- glm(Testo2X1 ~ Sex + DOY, 
             family = gaussian,
             data = allData)
fit.2 <- glm(Testo2X1 ~ Sex, family = gaussian,
             data = allData)
fit.3 <- glm(Testo2X1 ~ Sex + DOY + Tail,
             family = gaussian,
             data = allData)

# plot the logistic regression
# x <- seq(from = min(Testo2X1), to = max(Testo2X1), by = 0.01)
# y <- 1 - inv.logit(coef(fit.2)[1] + coef(fit.2)[2] * x)
#x2 <- exp(log(max(datSDB$Testo2)) + x)
# plot(x, y)

# predict the unknowns
# unknown_Testo <- datSDB$Testo2X1[datSDB$Sex == "U"]
# unknown_ID <- datSDB$ID[datSDB$Sex == "U"]
# unknown_Testo0 <- datSDB$Testo2[datSDB$Sex == "U"]
# p_f <- 1 - inv.logit(coef(fit.2)[1] + coef(fit.2)[2] * unknown_Testo)
# 
# out.dat <- as.data.frame(cbind(unknown_ID, unknown_Testo0, p_f))
# 
# write.csv(out.dat, 'Pr_Female1.csv', quote = FALSE, row.names = FALSE)

# but they are not all independent - individuals are repeatedly
# measured. Also, want to use multiple measurements of each
# sample. 




