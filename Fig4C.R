require("lubridate")
require("doBy")
require("pals")#Function
library(geosphere)
require("plotrix")
require("raster")
require("rgdal")
require("dplyr")
require("doBy")
require("ncdf4")
require("lubridate")
library("harrypotter")
require("deSolve")
load(file="IntensityFlu.RData")


pdf("Fig4C.pdf", width = 4, height = 4)
par(mar=c(3,3,1,1))
plot(mg$lat_abs, mg$intensity, pch = 16, col="grey64", bty ="n", ylab="", xlab = "")
lm.out <- lm(intensity ~ lat_abs, data = mg)
newx = seq(min(mg$lat_abs),max(mg$lat_abs),by = 0.05)
conf_interval <- predict(lm.out, newdata=data.frame(lat_abs=newx), interval="confidence",
                         level = 0.95)
abline(lm.out, col="black", lwd = 2)
lines(newx, conf_interval[,2], col="grey64", lty=2, lwd = 2)
lines(newx, conf_interval[,3], col="grey64", lty=2, lwd = 2)
mg$mask <- 0
mg$mask[mg$newnames %in% c("China", "Japan", "South Korea", "Taiwan", "Hong Kong")] <- 1
sum(mg$mask)
points(mg$lat_abs[mg$mask==1], mg$intensity[mg$mask==1], pch = 16, col="blue")
points(mg$lat_abs[mg$newnames=="United Kingdom"], mg$intensity[mg$newnames=="United Kingdom"],col="firebrick", pch = 16)
points(mg$lat_abs[mg$newnames=="Germany"], mg$intensity[mg$newnames=="Germany"],col="firebrick", pch = 16)
points(mg$lat_abs[mg$newnames=="United States"], mg$intensity[mg$newnames=="United States"],col="firebrick", pch = 16)

text(mg$lat_abs[mg$mask==1], as.numeric(mg$intensity[mg$mask==1]), labels = mg$ISO[mg$mask==1], pos= 4, col="blue")
text(mg$lat_abs[mg$newnames=="United Kingdom"], as.numeric(mg$intensity[mg$newnames=="United Kingdom"]), 
     labels = mg$ISO[mg$newnames=="United Kingdom"], pos= 4, col="firebrick")
text(mg$lat_abs[mg$newnames=="Germany"], as.numeric(mg$intensity[mg$newnames=="Germany"]), 
     labels = mg$ISO[mg$newnames=="Germany"], pos= 4, col="firebrick")
text(mg$lat_abs[mg$newnames=="United States"], as.numeric(mg$intensity[mg$newnames=="United States"]), 
     labels = mg$ISO[mg$newnames=="United States"], pos= 4, col="firebrick")

title(xlab = "Latitude", line = 2)
title(ylab = "Intensity", line = 2)

dev.off()

lm.out <- lm(intensity~ mask + lat_abs, data = mg)
summary(lm.out)


