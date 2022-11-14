##### seasonality - FLU
require("deSolve")
require("scales")
require('zoo')
sirsmod_cs = function(t, y, parameters) {
  #browser()
  S = y[1]
  I = y[2]
  R = y[3]
  mu = parameters[["mu"]]
  gamma = parameters[["gamma"]]
  omega = parameters[["omega"]]
  N = parameters[["N"]]
  controlstart = parameters[["controlstart"]]
  controlend = parameters[["controlend"]]
  controlstart2 = parameters[["controlstart2"]]
  controlend2 = parameters[["controlend2"]]
  betachange = parameters[["betachange"]]
  
  beta=R0.list[t]*gamma
  if(t >= controlstart & t < controlend){
    beta = betachange*beta
  }
  if(t >= controlstart2 & t < controlend2){
    beta = betachange*beta
  }
  dS = omega*R + mu * (N - S) - beta * S * I/N
  dI = beta * S * I/N - (mu + gamma) * I
  dR = gamma*I - omega*R - mu * R
  res = c(dS, dI, dR)
  list(res)
}

week <- seq(1,52,1)
qsim <-(0.006/2)*sin(2*pi*week/52 + 10.5) + 0.006 # based on fitting to NY q data
qsim[qsim < 0] <- 0
R0min = 1.2
R0max = 3
R0 = exp(-180*qsim + log(R0max - R0min)) + R0min
mean(R0)
times = seq(1, 52 * 100,by = 1)
R0.list <- rep(R0,length=length(times))
out = as.data.frame(ode(y = c(S = 0.19, I = 0.01, R = 0.8), times = times, 
                        func = sirsmod_cs, parms = list(mu = 1/(50 * 52), N = 1, R0.list = R0.list, gamma = 1,
                                                        omega = 1/40,
                                                        controlstart = (52*43 + 11), 
                                                        controlend =(52*43 + 11 + 52), betachange = 0.8, 
                                                        controlstart2 = (52*43 + 11 + 104) ,  controlend2 = (52*400))))



subsetI <- out$I[times > (39*52)]
subsetS <- out$S[times > (39*52)]

timeuse <- seq(2016,2100,1/52)[1:length(subsetI)]
controlWeekStart = 11

setwd("~/Dropbox/Controls/Plots/Texas_RSV_SIMS")
pdf("Flu_Main.pdf", width=6, height = 4)
par(mar=c(3,3,1,3))
cols <- c(dichromat_pal("DarkRedtoBlue.12")(12)[1:4],dichromat_pal("DarkRedtoBlue.12")(12)[9:12])
cols2 <- c(dichromat_pal("BluetoOrange.10")(10))

timeuse <- seq(2016,2100,1/52)[1:length(subsetI)]
plot(timeuse,subsetI,type="n",col="#F64740",lwd =2,xlab="",ylab=
       "",bty = "n",xaxs="i",xlim=c(2016,2035),  yaxt ="n", xaxt= "n", ylim=c(0,0.3))
polygon(c(timeuse[52*4+controlWeekStart],timeuse[52*4+controlWeekStart],timeuse[52*4+controlWeekStart + 52],
          timeuse[52*4+controlWeekStart + 52]), c(-1,1,1,-1), border = NA,col="grey88")
polygon(c(timeuse[52*4+controlWeekStart + 104],timeuse[52*4+controlWeekStart + 104],max(timeuse),
          max(timeuse)), c(-1,1,1,-1), border = NA,col="grey88")
axis(side = 2, col.axis= cols[7], at = c(0.000,0.04, 0.08), label = c("0.00","0.04","0.08"))
lines(timeuse,subsetI,col=cols[7],lwd =2,xlab="",ylab=
        "")
subsetI2 <- out$I[times > (38.5*52)]
rm <- rollmean(subsetI2, k = 52 )[1:length(timeuse)]
lines(timeuse, rm, col="grey30", lwd = 2, lty = 2)
mtext(side= 2, "I/N", line = 2,col = cols[7], at = 0.001)
mtext(side= 2, "mean(I/N)", line = 2,col="grey30", at =0.06)
axis(side = 1, at = seq(2015,2035,5), labels =  seq(2015,2035,5))
abline(v = seq(2016,2100,1),lty=2,col="gray")

text(2030, 0.3,expression('"Influenza like": Low R'[0]*', SIRS'))


par(new = TRUE)
plot(timeuse, subsetS, col=cols[2],type="l",xlab="",ylab="", 
     axes = F,xaxs="i", ylim=c(0,1), lwd = 2, lty = 1,xlim=c(2016,2035))
axis(side = 4, col.axis=cols[2], at = c(0.4, 0.6), label = c(0.4, 0.6))
title(xlab="Year", line = 2)
mtext(side = 4, line = 2, 'S/N',col=cols[2])

par(new = TRUE)
subsetS2 <- out$S[times > (38.5*52)]
beta2 = c(R0[26:52] , rep(R0, length = length(subsetS2)))[1:length(subsetS2)]
beta2[245:length(beta2)] <- beta2[245:length(beta2)]*0.8

reff = beta2*subsetS2

rmreff <- rollmean(reff, k = 52 )[1:length(timeuse)]
plot(timeuse, rmreff, ylab = " ", xlab = "", yaxt = "n", xaxt = "n", ylim=c(-2.5,1.5), type="l", lwd =2, col =cols2[9],xlim=c(2016,2035),bty = "n",
     xaxs="i", lty = 2)
axis(side = 2, col.axis=cols2[9], at = c(0.5,1,1.5), label = c(0.5,1,1.5))

mtext(side = 2, line = 2, expression('mean(R'[e]*')'),col=cols2[9], at= 1, lty  = 1)



dev.off()

