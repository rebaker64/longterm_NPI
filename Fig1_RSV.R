library(tsiR)
library("doBy")
library("scales")
library("zoo")
load("TX.RData")

####### Part 1
## Fit TSIR model
tsir <- runtsir(data=data, IP = 1, xreg ="cumcases",regtype= "spline",userYhat = NULL,
                alpha=0.97, family='poisson',link='log',method='negbin') 

# plot results
pdf("texas_tsir_fit.pdf", width = 5, height =3)
par(mar=c(3,3,1,1))
plot(tsir$res$time, tsir$res$cases, pch = 16, col="grey64", xlab="", ylab = "",ylim=c(0,1200))
polygon(c(tsir$res$time, rev(tsir$res$time))  ,c(tsir$res$mean + 1.96*tsir$res$sd, rev(tsir$res$mean - 1.96*tsir$res$sd)), border = NA, col=rgb(0,0,1,0.2))
title(xlab = "Year", line = 2)
title(ylab = "Cases", line = 2)
lines(tsir$res$time, tsir$res$mean, col="navy", lwd = 2)
dev.off()

#plot beta
pdf("texas_beta.pdf", width = 5, height =3)
par(mar=c(3,3,1,1))
plot(tsir$beta*mean(data$pop), type="n", xlab="", ylab = "", bty="n")
abline(h = mean(tsir$beta*mean(data$pop)),col="grey1", lty = 3)
polygon(c(seq(1,52,1), rev(seq(1,52,1)))  ,c(tsir$contact$betahigh*mean(data$pop), rev(tsir$contact$betalow)*mean(data$pop)), border = NA, col="grey64")
title(xlab="Week", line = 2)
title(ylab="Beta (t)", line = 2)
lines(tsir$beta*mean(data$pop), col="black", lwd = 2)
dev.off()

#estimate for R0
mean(tsir$beta*mean(data$pop))

####### Part 2
source('predtsirMultiControlVax.R', encoding = 'UTF-8')
# now we are going to predict forward, we need a dataset with the seasonal transmission rate and average births/pop 
# we can also use as seasonal birth rate here 
data2 <- data.frame(seas_beta_order <- seq(1,52,1), sea_beta = tsir$beta, births = 368190/52, pop = 29145505)

controlWeekStart = 11 # week in year that NPIs go into place

times <- seq(1,100,1/52) 
controlStart = 43*52 + controlWeekStart # week in simulation controls go in place #NB for endemic infections important to remove burn in period/transient dynamics
controlEnd =  controlStart + 52 # first controls in place for a year i.e. COVID-19 controls
p2controlStart = controlEnd + 52 # add controls back in after a year
p2controlEnd = length(times)
#p2controlEnd = p2controlStart + 8

pred <- predtsirMultiControlVax(times = times, births = rep(data2$births, length = length(times)), beta = data2$sea_beta, alpha = 0.97, 
                         S0 =floor(0.8*data2$pop[1]), I0 = floor(0.2*data2$pop[1]), nsim = 10, stochastic = F, 
                         controlStart = controlStart, controlEnd =controlEnd, betachange = 0.8,
                         p2controlStart =  p2controlStart, p2controlEnd = p2controlEnd , p2betachange = 0.8)
# I wrote this function, but it just does the forward simulation and allows you to put in a control period

subsetI <- pred$I$mean[times > 40]
subsetS <- pred$S$mean[times > 40]
timeuse <- seq(2016,2100,1/52)[1:length(subsetI)]

pdf("Texas_Main.pdf", width=6, height = 4)
par(mar=c(3,3,1,3))
cols <- c(dichromat_pal("DarkRedtoBlue.12")(12)[1:4],dichromat_pal("DarkRedtoBlue.12")(12)[9:12])
cols2 <- c(dichromat_pal("BluetoOrange.10")(10))

popuse <- mean(data$pop)
timeuse <- seq(2016,2100,1/52)[1:length(subsetI)]
plot(timeuse,subsetI/popuse,type="n",col="#F64740",lwd =2,ylim=c(0,0.005),xlab="",ylab=
       "",bty = "n",xaxs="i",xlim=c(2016,2035), main="Texas", yaxt ="n", xaxt= "n")
polygon(c(timeuse[52*4+controlWeekStart],timeuse[52*4+controlWeekStart],timeuse[52*4+controlWeekStart + 52],
          timeuse[52*4+controlWeekStart + 52]), c(-1,1,1,-1), border = NA,col="grey88")
polygon(c(timeuse[52*4+controlWeekStart + 104],timeuse[52*4+controlWeekStart + 104],max(timeuse),
          max(timeuse)), c(-1,1,1,-1), border = NA,col="grey88")
axis(side = 2, col.axis= cols[7], at = c(0.000,0.001,0.002,0.003), label = c("0.000","0.001","0.002","0.003"))
lines(timeuse,subsetI/popuse,col=cols[7],lwd =2,ylim=c(0,0.01),xlab="",ylab=
        "")
subsetI2 <- pred$I$mean[times > 39.5]/data2$pop[1]
rm <- rollmean(subsetI2, k = 52 )[1:length(timeuse)]
lines(timeuse, rm, col="grey30", lwd = 2, lty = 2)
mtext(side= 2, "I/N", line = 2,col = cols[7], at = 0.001)
mtext(side= 2, "mean(I/N)", line = 2,col="grey30", at =0.002)
axis(side = 1, at = seq(2015,2035,5), labels =  seq(2015,2035,5))
abline(v = seq(2016,2100,1),lty=2,col="gray")

text(2030, 0.005,expression('"RSV like": High R'[0]*', SIR'))

par(new = TRUE)
plot(timeuse, subsetS/popuse, col=cols[2],type="l",xlab="",ylab="", 
     axes = F,xaxs="i", ylim=c(0.1,0.3), lwd = 2, lty = 1,xlim=c(2016,2035))
axis(side = 4, col.axis=cols[2], at = c(0.15,0.20,0.25), label = c(0.15,"0.20",0.25))
title(xlab="Year", line = 2)
mtext(side = 4, line = 2, 'S/N',col=cols[2])

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

