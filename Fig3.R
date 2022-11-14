
# Figure 3A
library(tsiR)
library("doBy")
library("zoo")
load("TX.RData")
# plot the data
plot(data$time,data$cases,type="l")

tsir <- runtsir(data=data, IP = 1, xreg ="cumcases",regtype= "spline",userYhat = NULL,
                alpha=0.97, family='poisson',link='log',method='negbin') 

# checking model fit against data
plot(tsir$res$time, tsir$res$cases, type="l", lwd = 2, col="black")
lines(tsir$res$time, tsir$res$mean, type="l", lwd = 2, col="red")

####### Part 2
source('predtsirMultiControlVax.R', encoding = 'UTF-8')
# now we are going to predict forward, we need a dataset with the seasonal transmission rate and average births/pop 
data2 <- data.frame(seas_beta_order <- seq(1,52,1), sea_beta = tsir$beta, births = 368190/52, pop = 29145505)

controlWeekStart = 11 # week in year that NPIs go into place

times <- seq(1,100,1/52) # when simulating disease models sometimes we have to simulate for a few years to remove transient dynamics 
controlStart = 43*52 + controlWeekStart
controlEnd  = length(times) # controls in place forever

pred <- predtsirMultiControlVax(times = times, births = rep(data2$births, length = length(times)), beta = data2$sea_beta, alpha = 0.97, 
                         S0 =floor(0.8*data2$pop[1]), I0 = floor(0.2*data2$pop[1]), nsim = 10, stochastic = F, 
                         controlStart = controlStart, controlEnd =controlEnd, betachange = 0.8,
                         p2controlStart =  -1, p2controlEnd = -1, p2betachange = 1,
                         vacc.cov = 0.8,  vaccStart = controlStart +156, vaccEnd = length(times))
#This function does the forward simulation and allows you to put in a control period

subsetI <- pred$I$mean[times > 40]
subsetS <- pred$S$mean[times > 40]
timeuse <- seq(2016,2100,1/52)[1:length(subsetI)]
max(subsetI/data2$pop[1])
max(subsetS/data2$pop[1])

# Plot result
pdf("Texas_Vacc.pdf", width=6, height = 4)
par(mar=c(3,3,1,3))
par(mar=c(3,3,1,3))
cols <- c(dichromat_pal("DarkRedtoBlue.12")(12)[1:4],dichromat_pal("DarkRedtoBlue.12")(12)[9:12])
cols2 <- c(dichromat_pal("BluetoOrange.10")(10))

popuse <- mean(data$pop)
timeuse <- seq(2016,2100,1/52)[1:length(subsetI)]
plot(timeuse,subsetI/popuse,type="n",col="#F64740",lwd =2,ylim=c(0,0.005),xlab="",ylab=
       "",bty = "n",xaxs="i",xlim=c(2016,2035), main=stateslist[i], yaxt ="n", xaxt= "n")
polygon(c(timeuse[52*4+controlWeekStart],timeuse[52*4+controlWeekStart],max(timeuse),
                    max(timeuse)), c(-1,1,1,-1), border = NA,col="grey88")
polygon(c(timeuse[52*4+controlWeekStart + 156],timeuse[52*4+controlWeekStart + 156],max(timeuse),
          max(timeuse)), c(-1,1,1,-1), border = NA,col="bisque")
axis(side = 2, col.axis= cols[7], at = c(0.000,0.001,0.002,0.003), label = c("0.000","0.001","0.002","0.003"))
arrows(x0 = c(timeuse[52*4+controlWeekStart + 156]), x1 = c(timeuse[52*4+controlWeekStart + 156]), y0 = -0.2, y1 = 0.0052,
       lty = 3, length = 0.1, col="grey40")
arrows(x0 = c(timeuse[52*4+controlWeekStart]), x1 = c(timeuse[52*4+controlWeekStart]), y0 = -0.2, y1 = 0.0052,
       lty = 3, length = 0.1, col="grey40")
text(x = c(timeuse[52*4+controlWeekStart + 156]), y = 0.005,  substitute(paste(italic("Vacc. introduced"))), cex = 0.6, col="grey40")
text(x = c(timeuse[52*4+controlWeekStart]), y = 0.005,  substitute(paste(italic("NPI. introduced"))), cex = 0.6, col="grey40")
lines(timeuse,subsetI/popuse,col=cols[7],lwd =2,ylim=c(0,0.01),xlab="",ylab=
        "")
subsetI2 <- pred$I$mean[times > 39.5]/data2$pop[1]
rm <- rollmean(subsetI2, k = 52 )[1:length(timeuse)]
lines(timeuse, rm, col="grey30", lwd = 2, lty = 2)
mtext(side= 2, "I/N", line = 2,col = cols[7], at = 0.001)
mtext(side= 2, "mean(I/N)", line = 2,col="grey30", at =0.002)
axis(side = 1, at = seq(2015,2035,5), labels =  seq(2015,2035,5))

par(new = TRUE)
plot(timeuse, subsetS/popuse, col=cols[2],type="l",xlab="",ylab="", 
     axes = F,xaxs="i", ylim=c(0.1,0.3), lwd = 2, lty = 1,xlim=c(2016,2035))
axis(side = 4, col.axis=cols[2], at = c(0.15,0.20,0.25), label = c(0.15,"0.20",0.25))
title(xlab="Year", line = 2)
mtext(side = 4, line = 2, 'S/N',col=cols[2])
abline(v = seq(2016,2100,1),lty=2,col="gray")

par(new = TRUE)
subsetS2 <- pred$S$mean[times > 39.5]
beta2 = c(data2$sea_beta[26:52] , rep(data2$sea_beta, length = length(subsetS2)))[1:length(subsetS2)]
beta2[245:length(beta2)] <- beta2[245:length(beta2)]*0.8
reff = beta2*subsetS2
rmreff <- rollmean(reff, k = 52 )[1:length(timeuse)]
plot(timeuse, rmreff, ylab = " ", xlab = "", yaxt = "n", xaxt = "n", ylim=c(-2.5,1.5), type="l", lwd =2, col =cols2[9],xlim=c(2016,2035),bty = "n",
     xaxs="i", lty = 2)
axis(side = 2, col.axis=cols2[9], at = c(0.5,1,1.5), label = c(0.5,1,1.5))
mtext(side = 2, line = 2, expression('mean(R'[e]*')'),col=cols2[9], at= 1, lty  = 1)
dev.off()




##### Figure 3B
require("deSolve")
require("scales")
mu = 1/(50 * 52)
gamma = 0.5


R0 = seq(1,15,0.01)
Istar = (mu*(R0-1))/(R0*(gamma + mu))
Sstar = 1/R0
Reff  = R0*Sstar

### Plot 
pdf("Fig3B.pdf", width = 5, height =4)
cols <- c(dichromat_pal("DarkRedtoBlue.12")(12))
cols <- c(cols[1:5],cols[8:12])
par(mar=c(3,3,1,3))
plot(R0,Istar, type="l", lwd =2, bty = "n",ylim=c(0,max(Istar) + 0.00003))
pvec <- seq(0,1,0.1)
Istore <- NULL
for(i in 1:length(pvec)){
  IstarVac =  (mu*(R0*(1-pvec[i]) -1))/(R0*(gamma + mu))
  IstarVac[IstarVac <= 0] <- NA
  lines(R0, IstarVac, lwd = 2, col=cols[i])
  Istore <- c(Istore, IstarVac[1001])
}
text(11,Istore+0.00003, c("0%","10%","20%","30%","40%","50%","60%","70%","80%","90%"), pos = 4, col = c("black",cols[2:10]),cex = 0.8)
text(9, Istore[1] + 0.00003, substitute(paste(italic("Vacc. coverage ="))), cex = 0.8)
lines(R0,Istar, type="l", lwd =2, bty = "n")
par(new = TRUE)
plot(R0,Sstar, type="l", lwd =2, bty = "n", ylab="", yaxt="n", col="black", lty = 3)
axis(side=4, at = pretty(range(Sstar)), col="black")
mtext("S*", side=4, line=2, col="black")
title(xlab= expression('R'[0]*''), line = 2)
title(ylab="I*", line = 2)
dev.off()





