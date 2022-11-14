require("deSolve")
require("scales")

######################################################### Fig 2a
mu =  1/(50 * 52)
gamma = 1/2

R0 = seq(1,15,0.25)
Istar = (mu*(R0-1))/(R0*(gamma + mu))
Sstar = 1/R0
Reff  = R0*Sstar
R01 = 8
Istar1 = (mu*(R01-1))/(R01*(gamma + mu))
Sstar1 = 1/R01
R02 = (0.8*8)
Istar2 = (mu*(R02-1))/(R02*(gamma + mu))
Sstar2 = 1/R02

cols <- rev(dichromat_pal("LightBluetoDarkBlue.10")(10))
#setwd("~/Dropbox/Controls/Plots")
pdf("Fig1AR0.pdf", width = 4, height =2.5)
par(mar=c(3,3,1,1))
par(mar=c(3,3,3,3))
plot(R0,Istar, type="l", lwd =2)
segments(x0 = R01,y0=0, x1 = R01, y1 =Istar1,col="grey64", lty=3, lwd = 2)
segments(x0 = 0,y0=Istar1, x1 = R01, y1 =Istar1,col="grey64", lty=3, lwd = 2)
segments(x0 = R02,y0=0, x1 = R02, y1 =Istar2,col="grey64", lty=2, lwd = 2)
segments(x0 = 0,y0=Istar2, x1 = R02, y1 =Istar2,col="grey64", lty=2, lwd = 2)

par(new = TRUE)
plot(R0,Sstar, type="l", lwd =2, bty = "n", ylab="", yaxt="n", col=cols[1], ylim=c(0,1))
axis(side=4, at = pretty(range(Sstar)), col=cols[1])
#axis(side=4, col=cols[1])
mtext("S*", side=4, line=2, col=cols[1])
title(xlab=expression('R'[0]), line = 2)
title(ylab="I*", line = 2)
dev.off()

simod_c = function(t, y, parameters) {
  S = y[1]
  I = y[2]
  beta = parameters["beta"]
  mu = parameters["mu"]
  gamma = parameters["gamma"]
  N = parameters["N"]
  controlstart = parameters["controlstart"]
  controlend = parameters["controlend"]
  betachange = parameters["betachange"]
  if(t >= controlstart & t < controlend){
    beta = betachange*beta
  }
  dS = mu * (N - S) - beta * S * I/N
  dI = beta * S * I/N - (mu + gamma) * I
  res = c(dS, dI)
  list(res)
}

out = as.data.frame(ode(y = c(S = 0.08, I = 0.01), 
                        times = seq(0, 52 * 400,by = 1), 
                        func = simod_c, 
                        parms = c(mu = 1/(50 * 52), N = 1, beta = 4, 
                                  gamma = 1/2, controlstart = 52*200, 
                                  controlend =52*400, betachange = 0.8  )))

pdf("Fig1ASI.pdf", width = 4, height =4)
par(mar=c(3,3,1,1))
cols <- c(dichromat_pal("DarkRedtoBlue.12")(12))
cols <- c(cols[1:4],cols[9:12],"black")
plot(out$S, out$I, xlab = "", ylab = "", col= cols[1], type="n")
segments(x0 = head(out$S,-1),
         y0 = head(out$I,-1),
         x1 = tail(out$S,-1),
         y1 = tail(out$I,-1),
         lwd = 2,
         col = gradient_n_pal(cols)(seq(0, 1, length.out = length(out$I))))
mu = 1/(50 * 52)
gamma = 0.5
R0 = 8
Istar1 = (mu*(R0-1))/(R0*(gamma + mu))
Sstar1 = 1/R0
R0 = 0.8*8
Istar2 = (mu*(R0-1))/(R0*(gamma + mu))
Sstar2 = 1/R0
abline(h = Istar1, col = "grey64", lwd=2, lty=3)
abline(v = Sstar1, col = "grey64", lwd=2, lty=3)
abline(h = Istar2, col = "grey64", lwd=2, lty=2)
abline(v = Sstar2, col = "grey64", lwd=2, lty=2)
title(xlab="S/N", line = 2)
title(ylab="I/N", line = 2)
dev.off()


pdf("Fig1Ats.pdf", width = 4, height =2)
par(mar=c(3,3,1,1))
cols <- c(dichromat_pal("DarkRedtoBlue.12")(12))
cols <- c(cols[1:4],cols[9:12],"black")
plot(seq(1, length(out$S),1), out$I, xlab = "", ylab = "", col= cols[1], type="n", bty = "n")
polygon( c(52*200, 52*200, 52*400,52*400), c(-1,1,1,-1), col="grey81", border=NA)

segments(x0 = head(seq(1, length(out$S),1),-1),
         y0 = head(out$I,-1),
         x1 = tail(seq(1, length(out$S),1),-1),
         y1 = tail(out$I,-1),
         lwd = 2,
         col = gradient_n_pal(cols)(seq(0, 1, length.out = length(out$I))))
mu = 1/(50 * 52)
gamma = 0.5
R0 = 8
Istar1 = (mu*(R0-1))/(R0*(gamma + mu))
Sstar1 = 1/R0
R0 = 0.8*8
Istar2 = (mu*(R0-1))/(R0*(gamma + mu))
Sstar2 = 1/R0
abline(h = Istar1, col = "grey64", lwd=2, lty=3)
abline(h = Istar2, col = "grey64", lwd=2, lty=2)
title(xlab="Time (weeks)", line = 2)
title(ylab="I/N", line = 2)
dev.off()


######################################################### Fig 2b
R0start <- seq(2,15,1)
R0change = seq(0,0.99,0.001)
matout <- matrix(NA, ncol=length(R0change), nrow=length(R0start))

require("fields")
cols <- c(dichromat_pal("BluetoOrange.10")(10))
cols <- c(cols[1:4],cols[7:10])
cols = gradient_n_pal(cols)(seq(0, 1, length = 14))

pdf("Fig2B.pdf", width = 4, height =4)
par(mar=c(3,3,1,1))
plot(R0change, seq(0,1,length= length(R0change)), type="n", ylab ="",
     xlab = "",xaxt="n", yaxt ="n", bty = "n")
abline(1,-1, lwd = 3, lty = 3, col="grey64")
for(i in 1:length(R0start)){
  Ipre  =  (mu*(R0start[i]-1))/(R0start[i]*(gamma + mu))
  R0post = (1-R0change)*R0start[i]
  Ipost =(mu*(R0post-1))/(R0post*(gamma + mu))
  Ipost[Ipost <0] <- 0
  Ipct = Ipost/Ipre
  lines(R0change,Ipct, col=cols[i], lwd = 2)
  matout[i,] <- Ipct
}
axis(1, at = seq(0,1,0.2), labels=c("0%","20%","40%","60%","80%","100%"))
axis(2, at = seq(0,1,0.2), labels=c("100%","80%","60%","40%","20%","0%"))
arrows(0.3,0.4, 0.8, 0.8,
       code = 3, col = "black", cex = 0.8, length = 0.1)
title(xlab=expression('Reduction in R'[0]), line = 2)
title(ylab="Reduction in I*", line = 2)
text( 0.86, 0.80,expression('High R'[0]*'(= 15)'),cex = 0.8, pos = 3)
text(0.3,0.4,expression('Low R'[0]*'(= 2)'), cex =0.8, pos = 2)
dev.off()

i = 2
Ipre  =  (mu*(R0start[i]-1))/(R0start[i]*(gamma + mu))
R0post = (1-R0change)*R0start[i]
Ipost =(mu*(R0post-1))/(R0post*(gamma + mu))
Ipost[Ipost <0] <- 0
Ipct = Ipost/Ipre
out <- 1 - Ipct[101]
######################################################### Fig 2c

cols <- rev(dichromat_pal("LightBluetoDarkBlue.10")(10))

mu =  1/(50 * 52)
gamma = 1/2
durImm = 52
delta = 1/durImm
R0 = seq(1,15,0.25)
R01 = 2.5
R02 = 0.8*2.5
Sstar = 1/R0
Istar = ((mu + delta)*(1 - (1/R0)))/(gamma + mu + delta)
Sstar1 = 1/R01
Istar1 = ((mu + delta)*(1 - (1/R01)))/(gamma + mu + delta)
Sstar2 = 1/R02
Istar2 = ((mu + delta)*(1 - (1/R02)))/(gamma + mu + delta)

pdf("Fig2CR0.pdf", width = 4, height =2.5)
par(mar=c(3,3,1,1))
par(mar=c(3,3,3,3))
#par(mfrow=c(1,2))
plot(R0,Istar, type="l", lwd =2, bty = "n", ylab = "", xlab ="")
segments(x0 = R01,y0=0, x1 = R01, y1 =Istar1,col="grey64", lty=3, lwd = 2)
segments(x0 = 0,y0=Istar1, x1 = R01, y1 =Istar1,col="grey64", lty=3, lwd = 2)
segments(x0 = R02,y0=0, x1 = R02, y1 =Istar2,col="grey64", lty=2, lwd = 2)
segments(x0 = 0,y0=Istar2, x1 = R02, y1 =Istar2,col="grey64", lty=2, lwd = 2)

par(new = TRUE)
plot(R0,Sstar, lwd = 2, col=cols[1], type="l", ylim=c(0,1), ylab = "", xlab = "", yaxt= "n")
axis(side=4, at = pretty(range(Sstar)), col=cols[1])
mtext("S*", side=4, line=2, col=cols[1])
title(xlab=expression('R'[0]), line = 2)
title(ylab="I*", line = 2)
dev.off()



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
  betachange = parameters[["betachange"]]
  
  beta=R0.list[t]*gamma
  if(t >= controlstart & t < controlend){
    beta = betachange*beta
  }
  dS = omega*R + mu * (N - S) - beta * S * I/N
  dI = beta * S * I/N - (mu + gamma) * I
  dR = gamma*I - omega*R - mu * R
  res = c(dS, dI, dR)
  list(res)
}

times = seq(1, 52 * 400,by = 1)
R0 = 2.5
R0.list <- rep(R0,length=length(times))
out = as.data.frame(ode(y = c(S = 0.19, I = 0.01, R = 0.8), times = times, 
                        func = sirsmod_cs, 
                        parms = list(mu = 1/(50 * 52),
                                     N = 1, R0.list = R0.list, 
                                     gamma = 1/2,
                                     omega = 1/52,
                                     controlstart = 52*200,
                                     controlend =52*400, betachange = 0.8  )))




pdf("Fig2Cts.pdf", width = 4, height =2)
par(mar=c(3,3,1,1))
cols <- c(dichromat_pal("DarkRedtoBlue.12")(12))
cols <- c(cols[1:4],cols[9:12],"black")
plot(seq(1, length(out$S),1), out$I, xlab = "", ylab = "", col= cols[1], type="n", bty = "n")
polygon( c(52*200, 52*200, 52*400,52*400), c(-1,1,1,-1), col="grey81", border=NA)

segments(x0 = head(seq(1, length(out$S),1),-1),
         y0 = head(out$I,-1),
         x1 = tail(seq(1, length(out$S),1),-1),
         y1 = tail(out$I,-1),
         lwd = 2,
         col = gradient_n_pal(cols)(seq(0, 1, length.out = length(out$I))))
mu = 1/(50 * 52)
gamma = 0.5
R01 = 2.5
R02 = 0.8*2.5
Sstar = 1/R0
Istar = ((mu + delta)*(1 - (1/R0)))/(gamma + mu + delta)
Sstar1 = 1/R01
Istar1 = ((mu + delta)*(1 - (1/R01)))/(gamma + mu + delta)
Sstar2 = 1/R02
Istar2 = ((mu + delta)*(1 - (1/R02)))/(gamma + mu + delta)
abline(h = Istar1, col = "grey64", lwd=2, lty=3)
abline(h = Istar2, col = "grey64", lwd=2, lty=2)
title(xlab="Time (weeks)", line = 2)
title(ylab="I/N", line = 2)
dev.off()

pdf("Fig2CSI.pdf", width = 4, height =4)
par(mar=c(3,3,1,1))
cols <- c(dichromat_pal("DarkRedtoBlue.12")(12))
cols <- c(cols[1:4],cols[9:12],"black")
plot(out$S, out$I, xlab = "", ylab = "", col= cols[1], type="n")
segments(x0 = head(out$S,-1),
         y0 = head(out$I,-1),
         x1 = tail(out$S,-1),
         y1 = tail(out$I,-1),
         lwd = 2,
         col = gradient_n_pal(cols)(seq(0, 1, length.out = length(out$I))))

abline(h = Istar1, col = "grey64", lwd=2, lty=3)
abline(v = Sstar1, col = "grey64", lwd=2, lty=3)

abline(h = Istar2, col = "grey64", lwd=2, lty=2)
abline(v = Sstar2, col = "grey64", lwd=2, lty=2)
title(xlab="S/N", line = 2)
title(ylab="I/N", line = 2)

transientmin <- min(out$I[10400:20800])
transientX = out$S[which(out$I == transientmin)]
segments(x0 = transientX - 0.03 , y0 = transientmin, x1 = transientX + 0.03, y1 = transientmin, col="black")
text(x = transientX-0.035, y = transientmin-0.001,substitute(paste(italic("Transient min. I/N"))), col="black", cex = 0.8)

dev.off()




######################################################### Fig 2d
R0 = seq(1,15,0.25)
Istar = (mu*(R0-1))/(R0*(gamma + mu))
Sstar = 1/R0
pdf("Fig2D.pdf", width = 4, height =4)
par(mar=c(3,3,1,3))
#par(mfrow=c(1,2))
plot(R0,Istar, lwd =2, bty = "n", ylim=c(0,1), type="n")
cols <- c(dichromat_pal("DarkRedtoBlue.12")(12))
cols <- c(cols[1:4],cols[9:12])
col = gradient_n_pal(cols)(seq(0, 1, length = 20))
delta = seq(2,0.005, length=20)
for(i in 1:length(delta)){
  Sstar = 1/R0
  Istar = ((mu + delta[i])*(1 - (1/R0)))/(gamma + mu + delta[i])
  lines(R0,Istar, lwd = 2, col=col[i])
}
par(new = TRUE)
plot(R0,Sstar, type="l", lwd =2, bty = "n", ylab="", yaxt="n", col="black", ylim=c(0,1), lty = 3)
axis(side=4, at = pretty(range(Sstar)), col="black")
mtext("S*", side=4, line=2, col="black")
title(xlab=expression('R'[0]), line = 2)
title(ylab="I*", line = 2)
dev.off()



