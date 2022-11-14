##### seasonality - FLU
require("deSolve")
require("scales")

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

week <- seq(1,52,1)
qsim <-(0.006/2)*sin(2*pi*week/52 + 10.5) + 0.006 # based on fitting to NY q data
qsim[qsim < 0] <- 0
R0min = 1.2
R0max = 3
R0 = exp(-180*qsim + log(R0max - R0min)) + R0min
mean(R0)
times = seq(1, 52 * 400,by = 1)
R0.list <- rep(R0,length=length(times))
out = as.data.frame(ode(y = c(S = 0.3, I = 0.01, R = 0.8), times = times, 
                        func = sirsmod_cs, parms = list(mu = 1/(50 * 52), N = 1, R0.list = R0.list, gamma = 1,
                                                        omega = 1/40,
                                                        controlstart = 52*200, 
                                                        controlend =52*400, betachange = 0.8  )))
out <- out[(170*52):(230*52),]
cols <- c(dichromat_pal("DarkRedtoBlue.12")(12))


pdf("Fig4Alineplot.pdf", width = 2, height =1.5)
par(mar=c(3,3,1,1))
plot(out$I[1:52],type="l", lwd = 2, bty = "o", ylab = "", xlab = "", yaxt = "n", xaxt = "n",
     col = cols[3])
lines(out$I[3069:3121], lwd = 2, col = cols[12])
dev.off()




pdf("Fig4ASI.pdf", width = 4, height =4)
par(mar=c(3,3,1,1))
cols <- c(dichromat_pal("DarkRedtoBlue.12")(12))
cols <- c(cols[1],cols[1:4],cols[9:12])
plot(out$S, out$I, xlab = "", ylab = "", col= cols[1], type="n")
segments(x0 = head(out$S,-1),
         y0 = head(out$I,-1),
         x1 = tail(out$S,-1),
         y1 = tail(out$I,-1),
         lwd = 2,
         col = gradient_n_pal(cols)(seq(0, 1, length.out = length(out$I))))
mu = 1/(50 * 52)
gamma = 0.5
R0 = 2/(mu + gamma)
Istar1 = (mu*(R0-1))/(R0*(gamma + mu))
Sstar1 = 1/R0
R0 = (0.8*2)/(mu + gamma)
Istar2 = (mu*(R0-1))/(R0*(gamma + mu))
Sstar2 = 1/R0
title(xlab="S/N", line = 2)
title(ylab="I/N", line = 2)
dev.off()


betachangelist <- seq(0,1,0.01)
entout <- rep(NA, length = length(betachangelist))
meanIout <- rep(NA, length = length(betachangelist))
period  <- rep(1, length = length(betachangelist))
require("lomb")
for(i in 1:length(betachangelist)){
  out = as.data.frame(ode(y = c(S = 0.19, I = 0.01, R = 0.8), times = times, 
                          func = sirsmod_cs, parms = list(mu = 1/(50 * 52), N = 1, R0.list = R0.list, gamma = 1,
                                                          omega = 1/40,
                                                          controlstart = 52*200, 
                                                          controlend =52*400, betachange = betachangelist[i] )))
  subI <- out$I[15000:15208]
  meanIout[i] <-  mean(out$I[15000:15208])
  

  if(mean(subI) >0.0000001){
    ps <- subI/sum(subI)
    ps <- ps[ps != 0]
    entropy <- 1- (-1*sum(ps*log(ps)))/(-208*sum(1/208*log(1/208)))
    
    entout[i] <- entropy
    
    ls <- lsp(x = subI, times = seq(1,length(subI),1), type="period", from = 60 , to =200)
    if(ls$p.value < 0.05){
      period[i] <- 2
    }
  }
  
  
}


pdf("Fig4B.pdf", width = 4.1, height =4)
par(mar=c(3,3,1,3))
plot(1-betachangelist[period==1], entout[period==1], pch = 16, ylab = "", xlab = "", 
     bty = "n", xaxt = "n",ylim=c(min(entout, na.rm=T),max(entout, na.rm=T)),xlim=c(0,1))
points(1-betachangelist[period==2], entout[period==2], pch = 1, ylab = "", xlab = "", bty = "n", xaxt = "n",col="black")
axis(1, at = seq(0,1,0.1), c("0%", "10%", "20%" , "30%", "40%" ,"50%","60%","70%","80%","90%","100%"))
title(xlab=expression('Reduction in R'[0]), line = 2)
title(ylab="Intensity", line = 2)
par(new = TRUE)
plot(1-betachangelist, meanIout,ylab = "", xlab = "", bty = "n", xaxt = "n",  yaxt="n", 
     type="l", lwd = 2, col="blue", ylim=c(0,0.014))
axis(4)
mtext("mean I/N", side=4, line=2, col="blue")
dev.off()
