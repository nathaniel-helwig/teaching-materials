##########   Density and Distribution Estimation
##########   Nathaniel E. Helwig (helwig@umn.edu)
##########   Updated: 04-Jan-2017


#####   DEFINE PATHS AND PACKAGES   #####

# load quantreg package (for akj function)
if(!require(quantreg)){
  install.packages("quantreg")
  library("quantreg")
}


#####   EMPIRICAL CDF   #####

# plot uniform ecdf
dev.new(width=6,height=6,noRStudioGD=TRUE)
set.seed(1)
x = runif(20)
xseq = seq(0,1,length.out=100)
Fhat = ecdf(x)
plot(Fhat)
lines(xseq,punif(xseq),col="blue",lwd=2)

# plot normal ecdf
dev.new(width=6,height=6,noRStudioGD=TRUE)
set.seed(1)
x=rnorm(20)
xseq=seq(-3,3,length.out=100)
Fhat=ecdf(x)
plot(Fhat)
lines(xseq,pnorm(xseq),col="blue",lwd=2)


#####   HISTOGRAM ESTIMATES   #####

# plot uniform histograms
dev.new(width=12,height=4,noRStudioGD=TRUE)
par(mfrow=c(1,3))
set.seed(1)
x = runif(20)
hist(x,main="Sturges")
hist(x,breaks="FD",main="FD")
hist(x,breaks="scott",main="scott")

# plot normal histograms
dev.new(width=12,height=4,noRStudioGD=TRUE)
par(mfrow=c(1,3))
set.seed(1)
x = rnorm(20)
hist(x,main="Sturges")
hist(x,breaks="FD",main="FD")
hist(x,breaks="scott",main="scott")


#####   KERNEL DENSITY ESTIMATION   #####

# plot uniform kde
dev.new(width=12,height=4,noRStudioGD=TRUE)
par(mfrow=c(1,3))
set.seed(1)
x = runif(20)
kde = density(x)
plot(kde)
kde = density(x,kernel="epanechnikov")
plot(kde)
kde = density(x,kernel="rectangular")
plot(kde)

# plot normal kde
dev.new(width=12,height=4,noRStudioGD=TRUE)
par(mfrow=c(1,3))
set.seed(1)
x = rnorm(20)
kde = density(x)
plot(kde)
kde = density(x,kernel="epanechnikov")
plot(kde)
kde = density(x,kernel="rectangular")
plot(kde)

# write a kde function with standard normal kernel
kdenorm <- function(x,bw,q=NULL){
  if(is.null(q)) {
    q = seq(min(x)-3*bw, max(x)+3*bw, length.out=512)
  }
  nx = length(x)
  n = length(q)
  xmat = matrix(q,n,nx) - matrix(x,n,nx,byrow=TRUE)
  denall = dnorm(xmat/bw) / bw
  denhat = apply(denall,1,mean)
  list(x=q, y=denhat, bw=bw)
}

# try our function
dev.new(width=6,height=6,noRStudioGD=TRUE)
set.seed(1)
x = rnorm(100)
plot(density(x,bw=0.4),ylim=c(0,0.5))
kde = kdenorm(x,bw=0.4)
lines(kde,col="red")
lines(seq(-4,4,l=500),dnorm(seq(-4,4,l=500)),lty=2)

# plot normal kde different bandwidths
dev.new(width=12,height=4,noRStudioGD=TRUE)
par(mfrow=c(1,3))
set.seed(1)
x = rnorm(20)
kde = density(x,bw=0.1)
plot(kde)
kde = density(x)
plot(kde)
kde = density(x,0.7)
plot(kde)

# functions for CV evaluation
intfun <- function(ix,x,bw) kdenorm(x,bw,ix)$y^2
kdecv <- function(bw,x){
  lo = min(x) - 3*bw
  up = max(x) + 3*bw
  ival = integrate(intfun,x=x,bw=bw,lower=lo,upper=up)$value
  nx = length(x)
  ival - (2/(nx-1)) * sum( kdenorm(x,bw,x)$y - dnorm(0)/(nx*bw) )
}

# CV example
dev.new(width=12,height=6,noRStudioGD=TRUE)
par(mfrow=c(1,2))
set.seed(1)
xseq = seq(-4,4,length.out=500)
x = rnorm(100)
cvhat = rep(0,101)
htest = seq(0.05,1,length.out=101)
for(j in 1:101) cvhat[j] = kdecv(htest[j],x)
plot(htest,cvhat)
bwhat = htest[which.min(cvhat)]
kde = kdenorm(x,bw=bwhat)
plot(kde,main=bquote(bw==.(kde$bw)),type="l")
lines(xseq,dnorm(xseq),lty=2)

# plot normal kde different bandwidth rules
dev.new(width=12,height=12,noRStudioGD=TRUE)
par(mfrow=c(2,2))
set.seed(1)
x = rnorm(20)
kde = density(x,bw="nrd")
plot(kde)
kde = density(x,bw="nrd0")
plot(kde)
kde = density(x,bw="ucv")
plot(kde)
kde = density(x,bw="bcv")
plot(kde)

# plot normal kde different bandwidth rules
dev.new(width=6,height=6,noRStudioGD=TRUE)
set.seed(1)
x = rnorm(20)
xs = sort(x)
xseq = seq(-5,5,length=100)
kde = akj(xs,xseq)
plot(xseq,kde$dens,type="l")
