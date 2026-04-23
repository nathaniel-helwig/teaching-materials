##########   Introduction to Nonparametric Regression
##########   Nathaniel E. Helwig (helwig@umn.edu)
##########   Updated: 04-Jan-2017


#####   DEFINE PATHS AND PACKAGES   #####

# load fANCOVA package (for loess.as function)
if(!require(fANCOVA)){
  install.packages("fANCOVA")
  library("fANCOVA")
}

# load np package (for npreg function)
if(!require(np)){
  install.packages("np")
  library("np")
}


#####   NEED FOR NONPARAMETRIC REGRESSION   #####

# print data
anscombe

# summary statistics
summary(anscombe)

# fit four regression models
lm1 = lm(y1~x1, data=anscombe)
lm2 = lm(y2~x2, data=anscombe)
lm3 = lm(y3~x3, data=anscombe)
lm4 = lm(y4~x4, data=anscombe)

# plot four regression lines
dev.new(width=9,height=6,noRStudioGD=TRUE)
par(mfrow=c(2,2), mar=par()$mar+c(0,1,0,0))
xseq = seq(4,14,length=20)
plot(xseq,lm1$coef[1] + lm1$coef[2]*xseq,type="l",ylim=c(0,15),
     xlab=expression(italic(x)), ylab=expression(hat(italic(y))),
     main=expression("Study 1:   "*hat(y)*" = 3 + 0.5"*italic(x)),
     lwd=2, col="blue",cex.axis=1.25, cex.lab=1.5, cex.main=2)
text(5,13,expression(R^2*" = 0.67"))
plot(xseq,lm2$coef[1] + lm2$coef[2]*xseq,type="l",ylim=c(0,15),
     xlab=expression(italic(x)), ylab=expression(hat(italic(y))),
     main=expression("Study 2:   "*hat(y)*" = 3 + 0.5"*italic(x)),
     lwd=2, col="blue",cex.axis=1.25, cex.lab=1.5, cex.main=2)
text(5,13,expression(R^2*" = 0.67"))
plot(xseq,lm3$coef[1] + lm3$coef[2]*xseq,type="l",ylim=c(0,15),
     xlab=expression(italic(x)), ylab=expression(hat(italic(y))),
     main=expression("Study 3:   "*hat(y)*" = 3 + 0.5"*italic(x)),
     lwd=2, col="blue",cex.axis=1.25, cex.lab=1.5, cex.main=2)
text(5,13,expression(R^2*" = 0.67"))
plot(xseq,lm4$coef[1] + lm4$coef[2]*xseq,type="l",ylim=c(0,15),
     xlab=expression(italic(x)), ylab=expression(hat(italic(y))),
     main=expression("Study 4:   "*hat(y)*" = 3 + 0.5"*italic(x)),
     lwd=2, col="blue",cex.axis=1.25, cex.lab=1.5, cex.main=2)
text(5,13,expression(R^2*" = 0.67"))

# plot four regression lines
dev.new(width=9,height=6,noRStudioGD=TRUE)
par(mfrow=c(2,2), mar=par()$mar+c(0,1,0,0))
xseq = seq(4,14,length=20)

# plot four regression lines with data
plot(anscombe$x1,anscombe$y1,pch=19,ylim=c(0,15),
     xlab=expression(italic(x)), ylab=expression(hat(italic(y))),
     main=expression("Study 1:   "*hat(y)*" = 3 + 0.5"*italic(x)),
     lwd=2, col="red",cex.axis=1.25, cex.lab=1.5, cex.main=2)
abline(lm1,lwd=2,col="blue")
text(5,13,expression(R^2*" = 0.67"))
plot(anscombe$x2,anscombe$y2,pch=19,ylim=c(0,15),
     xlab=expression(italic(x)), ylab=expression(hat(italic(y))),
     main=expression("Study 2:   "*hat(y)*" = 3 + 0.5"*italic(x)),
     lwd=2, col="red",cex.axis=1.25, cex.lab=1.5, cex.main=2)
abline(lm2,lwd=2,col="blue")
text(5,13,expression(R^2*" = 0.67"))
plot(anscombe$x3,anscombe$y3,pch=19,ylim=c(0,15),
     xlab=expression(italic(x)), ylab=expression(hat(italic(y))),
     main=expression("Study 3:   "*hat(y)*" = 3 + 0.5"*italic(x)),
     lwd=2, col="red",cex.axis=1.25, cex.lab=1.5, cex.main=2)
abline(lm3,lwd=2,col="blue")
text(5,13,expression(R^2*" = 0.67"))
plot(anscombe$x4,anscombe$y4,pch=19,ylim=c(0,15),
     xlab=expression(italic(x)), ylab=expression(hat(italic(y))),
     main=expression("Study 4:   "*hat(y)*" = 3 + 0.5"*italic(x)),
     lwd=2, col="red",cex.axis=1.25, cex.lab=1.5, cex.main=2)
abline(lm4,lwd=2,col="blue")
text(5,13,expression(R^2*" = 0.67"))



#####   LOCAL AVERAGING   #####

# local averaging of sunspots data
data(sunspots)
yrs=start(sunspots)
yre=end(sunspots)
years=seq(yrs[1]+yrs[2]/12,yre[1]+yre[2]/12,by=1/12)
dev.new(width=8,height=6,noRStudioGD=TRUE)
par(mfrow=c(2,2))
plot(years,sunspots,type="l",main="Raw Data")
sunlocavg=supsmu(years,sunspots,span=0.01)
plot(sunlocavg,type="l",main="span=0.01")
sunlocavg=supsmu(years,sunspots,span=0.05)
plot(sunlocavg,type="l",main="span=0.05")
sunlocavg=supsmu(years,sunspots)
plot(sunlocavg,type="l",main="span=cv")

# local averaging of simulated data
set.seed(55455)
x=seq(0,1,length=50)
y=sin(2*pi*x)+rnorm(50,sd=0.5)
locavg=supsmu(x,y)
dev.new(width=6,height=6,noRStudioGD=TRUE)
plot(x,y)
lines(locavg$x,locavg$y)
lines(locavg$x,sin(2*pi*x),lty=2)
legend("bottomleft",c(expression(hat(eta)),expression(eta)),lty=1:2,bty="n")


#####   LOCAL REGRESSION   #####

# loess of sunspots data
data(sunspots)
yrs=start(sunspots)
yre=end(sunspots)
years=seq(yrs[1]+yrs[2]/12,yre[1]+yre[2]/12,by=1/12)
dev.new(width=8,height=6,noRStudioGD=TRUE)
par(mfrow=c(2,2))
plot(years,sunspots,type="l",main="Raw Data")
sunlocreg=loess(sunspots~years,span=0.01)
plot(sunlocreg$x,sunlocreg$fitted,type="l",main="span=0.01")
sunlocreg=loess(sunspots~years,span=0.05)
plot(sunlocreg$x,sunlocreg$fitted,type="l",main="span=0.05")
sunlocreg=loess.as(years,sunspots,criterion="gcv")
plot(sunlocreg$x,sunlocreg$fitted,type="l",main="span=gcv")

# loess of simulated data
set.seed(55455)
x=seq(0,1,length=50)
y=sin(2*pi*x)+rnorm(50,sd=0.5)
locreg=loess.as(x,y,criterion="gcv")
dev.new(width=6,height=6,noRStudioGD=TRUE)
plot(x,y)
lines(locreg$x,locreg$fitted)
lines(locreg$x,sin(2*pi*x),lty=2)
legend("bottomleft",c(expression(hat(eta)),expression(eta)),lty=1:2,bty="n")


#####   KERNEL REGRESSION   #####

# kernel regression of sunspots data
data(sunspots)
yrs=start(sunspots)
yre=end(sunspots)
sunspots=as.vector(sunspots)
years=seq(yrs[1]+yrs[2]/12,yre[1]+yre[2]/12,by=1/12)
dev.new(width=8,height=6,noRStudioGD=TRUE)
par(mfrow=c(2,2))
plot(years,sunspots,type="l",main="Raw Data")
sunkerreg=npreg(bws=0.01,txdat=years,tydat=sunspots)
plot(years,fitted(sunkerreg),type="l",main="bw=0.01")
sunkerreg=npreg(bws=0.05,txdat=years,tydat=sunspots)
plot(years,fitted(sunkerreg),type="l",main="bws=0.05")
sunkerreg=npreg(txdat=years,tydat=sunspots)
plot(years,fitted(sunkerreg),type="l",main="bws=cv (0.09)")

# kernel regression of simulated data
set.seed(55455)
x=seq(0,1,length=50)
y=sin(2*pi*x)+rnorm(50,sd=0.5)
kerreg=npreg(txdat=x,tydat=y)
dev.new(width=6,height=6,noRStudioGD=TRUE)
plot(x,y)
lines(x,fitted(kerreg))
lines(x,sin(2*pi*x),lty=2)
legend("bottomleft",c(expression(hat(eta)),expression(eta)),lty=1:2,bty="n")

