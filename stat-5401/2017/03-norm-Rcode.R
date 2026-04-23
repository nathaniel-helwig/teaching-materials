##########   Introduction to Normal Distribution
##########   Nathaniel E. Helwig (helwig@umn.edu)
##########   Updated: 16-Jan-2017


#####   DEFINE PATHS AND PACKAGES   #####

# load mvtnorm package (for multivariate normal calculations)
if(!require(mvtnorm)){
  install.packages("mvtnorm")
  library(mvtnorm) 
}

# load lattice package (for wireframe function)
if(!require(lattice)){
  install.packages("lattice")
  library(lattice)
}

# load bigsplines package (for imagebar function)
if(!require(bigsplines)){
  install.packages("bigsplines")
  library(bigsplines)
}


#####   UNIVARIATE NORMAL   #####

# plot standard normal
dev.new(height=5,width=8,noRStudioGD=TRUE)
par(mar=c(5,5.4,4,2)+0.1)
plot(seq(-4,4,l=200),dnorm(seq(-4,4,l=200)),type="l",lwd=2,
     xlab=expression(italic(x)),ylab=expression(italic(f(x))),
     cex.lab=2,cex.axis=2)


#####   BIVARIATE NORMAL   #####

# define bivariate normal pdf
dnormb <- function(x,y,mean=c(0,0),sd=c(0,0),rho=0){
  if(any(sd<=0L)){ stop('Input "sd" must be positive.') }
  if(abs(rho)>1L){ stop('Input "rho" must be <= 1 in magnitude.') }
  sfac = -1/(2*(1-rho^2))
  ex1p = ((x-mean[1])^2)/(sd[1]^2)
  ex2p = ((y-mean[2])^2)/(sd[2]^2)
  ex3p = 2*rho*(x-mean[1])*(y-mean[2])/prod(sd)
  top = exp(sfac*(ex1p+ex2p-ex3p))
  bot = 2*pi*prod(sd)*sqrt(1-rho^2)
  top/bot
}

# make some bivariate normal data
x = y = seq(-5, 5, length=50)
xy = expand.grid(x,y)
z = dnormb(xy[,1],xy[,2],sd=c(1,sqrt(2)),rho=.6/sqrt(2))

# plot bivariate normal wireframe
dev.new(height=8,width=8,noRStudioGD=TRUE)
wireframe(z~xy[,1]*xy[,2],xlab=list(label=expression(italic(x)),cex=2),
          ylab=list(label=expression(italic(y)),cex=2),
          zlab=list(label=expression(italic(f(x,y))),cex=2,rot=90),
          screen=list(z=-25,x=-60),scales=list(arrows=FALSE,cex=1.25,distance=c(1,1,1.25)),
          par.settings=list(axis.line=list(col="transparent")))

# plot bivariate normal using heat map
dev.new(height=8,width=8,noRStudioGD=TRUE)
imagebar(x,y,matrix(z,50,50),xlab=expression(italic(x)),ylab=expression(italic(y)),
         zlab=expression(italic(f(x,y))),zcex.lab=2,zcex.axis=1.25,cex.lab=2,cex.axis=1.25,zline=4)

# plot bivariate normal: different means
dev.new(height=5,width=15,noRStudioGD=TRUE)
par(mfrow=c(1,3))
x = y = seq(-5,5,l=50)
xy = expand.grid(x,y)
z = dnormb(xy[,1],xy[,2],mean=c(0,0),sd=c(1,sqrt(2)),rho=.6/sqrt(2))
imagebar(x,y,matrix(z,50,50),xlab=expression(italic(x)),ylab=expression(italic(y)),
         zlab=expression(italic(f(x,y))),zcex.lab=2,zcex.axis=1.25,zline=4,
         drawbar=T,cex.lab=2,cex.axis=1.25,cex.main=2.5,main=expression(italic(mu[x])*" = 0,  "*italic(mu[y])*" = 0"))
z = dnormb(xy[,1],xy[,2],mean=c(1,2),sd=c(1,sqrt(2)),rho=.6/sqrt(2))
imagebar(x,y,matrix(z,50,50),xlab=expression(italic(x)),ylab=expression(italic(y)),
         zlab=expression(italic(f(x,y))),zcex.lab=2,zcex.axis=1.25,zline=4,
         drawbar=T,cex.lab=2,cex.axis=1.25,cex.main=2.5,main=expression(italic(mu[x])*" = 1,  "*italic(mu[y])*" = 2"))
z = dnormb(xy[,1],xy[,2],mean=c(-1,-1),sd=c(1,sqrt(2)),rho=.6/sqrt(2))
imagebar(x,y,matrix(z,50,50),xlab=expression(italic(x)),ylab=expression(italic(y)),
         zlab=expression(italic(f(x,y))),zcex.lab=2,zcex.axis=1.25,zline=4,
         cex.lab=2,cex.axis=1.25,cex.main=2.5,main=expression(italic(mu[x])*" = -1,  "*italic(mu[y])*" = -1"))

# plot bivariate normal: different correlations
dev.new(height=5,width=15,noRStudioGD=TRUE)
par(mfrow=c(1,3))
x = y = seq(-5,5,l=50)
xy = expand.grid(x,y)
z = dnormb(xy[,1],xy[,2],mean=c(0,0),sd=c(1,sqrt(2)),rho=-.6/sqrt(2))
imagebar(x,y,matrix(z,50,50),xlab=expression(italic(x)),ylab=expression(italic(y)),
         zlab=expression(italic(f(x,y))),zcex.lab=2,zcex.axis=1.25,zline=4,
         drawbar=T,cex.lab=2,cex.axis=1.25,cex.main=2.5,main=expression(italic(rho)*" = -0.6/"*sqrt(2)))
z = dnormb(xy[,1],xy[,2],mean=c(0,0),sd=c(1,sqrt(2)),rho=0)
imagebar(x,y,matrix(z,50,50),xlab=expression(italic(x)),ylab=expression(italic(y)),
         zlab=expression(italic(f(x,y))),zcex.lab=2,zcex.axis=1.25,zline=4,
         drawbar=T,cex.lab=2,cex.axis=1.25,cex.main=2.5,main=expression(italic(rho)*" = 0"))
z = dnormb(xy[,1],xy[,2],mean=c(0,0),sd=c(1,sqrt(2)),rho=1.2/sqrt(2))
imagebar(x,y,matrix(z,50,50),xlab=expression(italic(x)),ylab=expression(italic(y)),
         zlab=expression(italic(f(x,y))),zcex.lab=2,zcex.axis=1.25,zline=4,
         cex.lab=2,cex.axis=1.25,cex.main=2.5,main=expression(italic(rho)*" = 1.2/"*sqrt(2)))

# plot bivariate normal: different variances
dev.new(height=5,width=15,noRStudioGD=TRUE)
par(mfrow=c(1,3))
x = y = seq(-5,5,l=50)
xy = expand.grid(x,y)
mysd = c(1,1)
z = dnormb(xy[,1],xy[,2],mean=c(0,0),sd=mysd,rho=.6/prod(mysd))
imagebar(x,y,matrix(z,50,50),xlab=expression(italic(x)),ylab=expression(italic(y)),
         zlab=expression(italic(f(x,y))),zcex.lab=2,zcex.axis=1.25,zline=4,
         drawbar=T,cex.lab=2,cex.axis=1.25,cex.main=2.5,main=expression(italic(sigma[y])*" = 1"))
mysd = c(1,sqrt(2))
z = dnormb(xy[,1],xy[,2],mean=c(0,0),sd=mysd,rho=.6/prod(mysd))
imagebar(x,y,matrix(z,50,50),xlab=expression(italic(x)),ylab=expression(italic(y)),
         zlab=expression(italic(f(x,y))),zcex.lab=2,zcex.axis=1.25,zline=4,
         drawbar=T,cex.lab=2,cex.axis=1.25,cex.main=2.5,main=expression(italic(sigma[y])*" = "*sqrt(2)))
mysd = c(1,2)
z = dnormb(xy[,1],xy[,2],mean=c(0,0),sd=mysd,rho=.6/prod(mysd))
imagebar(x,y,matrix(z,50,50),xlab=expression(italic(x)),ylab=expression(italic(y)),
         zlab=expression(italic(f(x,y))),zcex.lab=2,zcex.axis=1.25,zline=4,
         cex.lab=2,cex.axis=1.25,cex.main=2.5,main=expression(italic(sigma[y])*" = 2"))

# get bivariate normal CDF values
x = y = seq(-5,5,l=50)
xy = as.matrix(expand.grid(x,y))
sigma = matrix(c(1,.6,.6,2),2,2,byrow=T)
z = rep(NA,2500)
for(jj in 1:2500){
  z[jj] = pmvnorm(lower=-Inf,upper=xy[jj,],sigma=sigma)
}

# plot bivariate normal CDF
dev.new(height=8,width=8,noRStudioGD=TRUE)
imagebar(x,y,matrix(z,50,50),xlab=expression(italic(x)),ylab=expression(italic(y)),zlim=c(0,1),
         zlab=expression(italic(F(x,y))),zcex.lab=2,zcex.axis=1.25,cex.lab=2,cex.axis=1.25,zline=4)

# Example 1a
pnorm(1,lower=F)
pnorm(75,mean=60,sd=15,lower=F)

# Example 1b
pnorm(0.5,lower=F)
pnorm(75,mean=69,sd=12,lower=F)

# Example 1c
pnorm(20/sqrt(505),lower=F)
pnorm(150,mean=130,sd=sqrt(505),lower=F)

# Example 1d
pnorm(-10/sqrt(145),lower=F)
pnorm(0,mean=10,sd=sqrt(145),lower=F)

# Example 1e
pnorm(0.8,lower=F)
pnorm(150,mean=110,sd=50,lower=F)

# Example 2a
pnorm(1.5,lower=F)
pnorm(8,mean=5,sd=2,lower=F)

# Example 2b
pnorm(2.25/sqrt(119/32),lower=F)
pnorm(8,mean=5.75,sd=sqrt(119/32),lower=F)

# Example 2c
pnorm(1)
pnorm(63,mean=46,sd=17)
