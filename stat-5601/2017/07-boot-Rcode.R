##########   Bootstrap Resampling
##########   Nathaniel E. Helwig (helwig@umn.edu)
##########   Updated: 04-Jan-2017


#####   DEFINE PATHS AND PACKAGES   #####

# load MASS package (for rlm function)
if(!require(MASS)){
  install.packages("MASS")
  library("MASS")
}


#####   BACKGROUND INFORMATION   #####

# hypothetical ideal: example 1 (normal mean)
set.seed(1)
n = 100
B = c(200,500,1000,2000,5000,10000)
xseq = seq(-0.4,0.4,length=200)
dev.new(width=12,height=8,noRStudioGD=TRUE)
par(mfrow=c(2,3))
for(k in 1:6){
  X = replicate(B[k], rnorm(n))
  xbar = apply(X, 2, mean)
  hist(xbar, freq=F, xlim=c(-0.4,0.4), ylim=c(0,5),
       main=paste("Sampling Distribution: B =",B[k]))
  lines(xseq, dnorm(xseq, sd=1/sqrt(n)))
  legend("topright",expression(bar(x)*" pdf"),lty=1,bty="n")
}

# hypothetical ideal: example 2 (normal median)
set.seed(1)
n = 100
B = c(200,500,1000,2000,5000,10000)
xseq = seq(-0.4,0.4,length=200)
dev.new(width=12,height=8,noRStudioGD=TRUE)
par(mfrow=c(2,3))
for(k in 1:6){
  X = replicate(B[k], rnorm(n))
  xmed = apply(X, 2, median)
  hist(xmed, freq=F, xlim=c(-0.4,0.4), ylim=c(0,5),
       main=paste("Sampling Distribution: B =",B[k]))
  lines(xseq, dnorm(xseq, sd=1/sqrt(n)))
  legend("topright",expression(bar(x)*" pdf"),lty=1,bty="n")
}

# hypothetical ideal: example 3 (uniform mean)
set.seed(1)
n = c(3,5,10,20,50,100)
B = 10000
dev.new(width=12,height=8,noRStudioGD=TRUE)
par(mfrow=c(2,3))
for(k in 1:6){
  X = replicate(B, runif(n[k]))
  xbar = apply(X, 2, mean)
  xseq = seq(min(xbar)-0.05,max(xbar)+0.05,length=200)
  dhat = dnorm(xseq, mean=0.5, sd=1/sqrt(12*n[k]))
  hist(xbar, freq=F, ylim=c(0,max(dhat)+0.05),
       main=paste("Sampling Distribution: n =",n[k]))
  lines(xseq, dhat)
  legend("topright","asymp pdf",lty=1,bty="n")
}


#####   BOOTSTRAP BASICS   #####

# normal ecdf
dev.new(width=12,height=4,noRStudioGD=TRUE)
set.seed(1)
par(mfrow=c(1,3))
n = c(100,500,1000)
xseq = seq(-4,4,length=100)
for(j in 1:3){
  x = rnorm(n[j])
  plot(ecdf(x),main=paste("n =",n[j]))
  lines(xseq,pnorm(xseq),col="blue")
}

# t(10) ecdf
dev.new(width=12,height=4,noRStudioGD=TRUE)
set.seed(1)
par(mfrow=c(1,3))
n = c(100,500,1000)
xseq = seq(-4,4,length=100)
for(j in 1:3){
  x = rt(n[j],df=10)
  plot(ecdf(x),main=paste("n =",n[j]))
  lines(xseq,pt(xseq,df=10),col="blue")
}

# chi-square ecdf
dev.new(width=12,height=4,noRStudioGD=TRUE)
set.seed(1)
par(mfrow=c(1,3))
n = c(100,500,1000)
xseq = seq(0,40,length=100)
for(j in 1:3){
  x = rchisq(n[j],df=10)
  plot(ecdf(x),main=paste("n =",n[j]))
  lines(xseq,pchisq(xseq,df=10),col="blue")
}

# F(5,10) ecdf
dev.new(width=12,height=4,noRStudioGD=TRUE)
set.seed(1)
par(mfrow=c(1,3))
n = c(100,500,1000)
xseq = seq(0,20,length=100)
for(j in 1:3){
  x = rf(n[j],df1=5,df2=10)
  plot(ecdf(x),main=paste("n =",n[j]))
  lines(xseq,pf(xseq,df1=5,df2=10),col="blue")
}

# uniform ecdf
dev.new(width=12,height=4,noRStudioGD=TRUE)
set.seed(1)
par(mfrow=c(1,3))
n = c(100,500,1000)
xseq = seq(0,1,length=100)
for(j in 1:3){
  x = runif(n[j])
  plot(ecdf(x),main=paste("n =",n[j]))
  lines(xseq,punif(xseq),col="blue")
}



#####   BOOTSTRAP IN PRACTICE   #####

# bootstrap sampling function
bootsamp <- function(x,nsamp=10000){
  x = as.matrix(x)
  nx = nrow(x)
  bsamp = replicate(nsamp,x[sample.int(nx,replace=TRUE),])
}

# bootstrap statistic and standard error calculation
bootse <- function(bsamp,myfun,...){
  if(is.matrix(bsamp)){
    theta = apply(bsamp,2,myfun,...)
  } else {
    theta = apply(bsamp,3,myfun,...)
  }
  if(is.matrix(theta)){
    return(list(theta=theta,cov=cov(t(theta))))
  } else{
    return(list(theta=theta,se=sd(theta)))
  }
}

# ex 1: mean
dev.new(width=5,height=5,noRStudioGD=TRUE)
set.seed(1)
x = rnorm(500,mean=1)
bsamp = bootsamp(x)
bse = bootse(bsamp,mean)
mean(x)
sd(x)/sqrt(500)
bse$se
hist(bse$theta)

# ex 2: median
dev.new(width=5,height=5,noRStudioGD=TRUE)
set.seed(1)
x = rnorm(500,mean=1)
bsamp = bootsamp(x)
bse = bootse(bsamp,median)
median(x)
bse$se
hist(bse$theta)

# ex 3: variance
dev.new(width=5,height=5,noRStudioGD=TRUE)
set.seed(1)
x = rnorm(500,sd=2)
bsamp = bootsamp(x)
bse = bootse(bsamp,var)
var(x)
bse$se
hist(bse$theta)

# ex 4: mean difference
dev.new(width=5,height=5,noRStudioGD=TRUE)
set.seed(1)
x = rnorm(500,mean=3)
y = rnorm(500)
z = cbind(x,y)
bsamp = bootsamp(z)
myfun = function(z) mean(z[,1]) - mean(z[,2])
bse = bootse(bsamp,myfun)
myfun(z)
sqrt( (var(z[,1])+var(z[,2]))/nrow(z) )
bse$se
hist(bse$theta)

# ex 5: median difference
dev.new(width=5,height=5,noRStudioGD=TRUE)
set.seed(1)
x = rnorm(500,mean=3)
y = rnorm(500)
z = cbind(x,y)
bsamp = bootsamp(z)
myfun = function(z) median(z[,1]) - median(z[,2])
bse = bootse(bsamp,myfun)
myfun(z)
bse$se
hist(bse$theta)

# ex 6: correlation
dev.new(width=5,height=5,noRStudioGD=TRUE)
set.seed(1)
x = rnorm(500)
y = rnorm(500)
Amat = matrix(c(1,-0.25,-0.25,1),2,2)
Aeig = eigen(Amat,symmetric=TRUE)
evec = Aeig$vec
evalsqrt = diag(Aeig$val^0.5)
Asqrt = evec %*% evalsqrt %*% t(evec)
z = cbind(x,y)%*%Asqrt
bsamp = bootsamp(z)
myfun = function(z) cor(z[,1],z[,2])
bse = bootse(bsamp,myfun)
myfun(z)
(1-myfun(z)^2)/sqrt(nrow(z)-3)
bse$se
hist(bse$theta)

# ex 7: uniform max
dev.new(width=5,height=5,noRStudioGD=TRUE)
set.seed(1)
x = runif(500)
bsamp = bootsamp(x)
myfun = function(x) max(x)
bse = bootse(bsamp,myfun)
myfun(x)
bse$se
hist(bse$theta)

# bootstrap bias calculation
bootbias <- function(bse,theta,...){
  if(is.matrix(bse$theta)){
    return(apply(bse$theta,1,mean) - theta)
  } else{
    return(mean(bse$theta) - theta)
  }
}

# sample mean unbiased
set.seed(1)
x = rnorm(500,mean=1)
bsamp = bootsamp(x)
bse = bootse(bsamp,mean)
mybias = bootbias(bse,mean(x))
mybias
mean(x)
sd(x)/sqrt(500)
bse$se

# toy example bias
set.seed(1)
x = rnorm(500,mean=1)
bsamp = bootsamp(x)
bse = bootse(bsamp,function(x) mean(x)+10)
mybias = bootbias(bse,mean(x))
mybias
mean(x)
sd(x)/sqrt(500)
bse$se

# toy example mse
set.seed(1)
x = rnorm(500,mean=1)
bsamp = bootsamp(x)
bse = bootse(bsamp,function(x) mean(x)+10)
mybias = bootbias(bse,mean(x))
c(bse$se,mybias)
c(bse$se,mybias)^2
mse = (bse$se^2) + (mybias^2)
mse
100 + 1/length(x)

# plot accuracy vs precision
dev.new(width=8,height=8,noRStudioGD=TRUE)
par(mfrow=c(2,2))
sds = c(0.5,0.15)
mus = c(1,0)
ptitle = c("low precision","high precision")
atitle = c("low accuracy","high accuracy")
set.seed(1)
for(j in 1:2){
  for(k in 1:2){
    x = rnorm(500, mean=mus[j], sd=sds[k])
    y = rnorm(500, mean=1, sd=sds[k])
    plot(x,y,xlim=c(-3,3),ylim=c(-3,3),pch=19,cex=0.25,col="gray",
         main=paste0(atitle[j]," and ",ptitle[k]))
    points(0,1,col="red",pch=8,cex=2,lwd=2)
    legend("bottomright",c("Truth","Estimates"),col=c("red","gray"),
           pch=c(8,19),pt.cex=c(1,1),bty="n")
  }
}

# mean is smooth
meanfun <- function(x,z) mean(c(x,z))
set.seed(1)
z = rnorm(100)
x = seq(-4,4,length=200)
meanval = rep(0,200)
for(j in 1:200) meanval[j] = meanfun(x[j],z)
dev.new(width=6,height=6,noRStudioGD=TRUE)
plot(x,meanval,xlab=expression(x[1]),main="mean")

# median is unsmooth
medfun <- function(x,z) median(c(x,z))
set.seed(1)
z = rnorm(100)
x = seq(-4,4,length=200)
medval = rep(0,200)
for(j in 1:200) medval[j] = medfun(x[j],z)
dev.new(width=6,height=6,noRStudioGD=TRUE)
plot(x,medval,xlab=expression(x[1]),main="median")

# jackknife functions
jacksamp <- function(x){
  nx = length(x)
  jsamp = matrix(0,nx-1,nx)
  for(j in 1:nx) jsamp[,j] = x[-j]
  jsamp
}
jackse <- function(jsamp,myfun,...){
  nx = ncol(jsamp)
  theta = apply(jsamp,2,myfun,...)
  se = sqrt( ((nx-1)/nx) * sum( (theta-mean(theta))^2 ) )
  list(theta=theta,se=se)
}

# ex 1: mean (revisited)
dev.new(width=5,height=5,noRStudioGD=TRUE)
set.seed(1)
x = rnorm(500,mean=1)
jsamp = jacksamp(x)
jse = jackse(jsamp,mean)
mean(x)
sd(x)/sqrt(500)
jse$se
hist(jse$theta)

# ex 2: median (revisited)
dev.new(width=5,height=5,noRStudioGD=TRUE)
set.seed(1)
x = rnorm(500,mean=1)
jsamp = jacksamp(x)
jse = jackse(jsamp,median)
median(x)
jse$se
hist(jse$theta)


#####   BOOTSTRAPPING REGRESSION   #####

# bootstrap regression residuals
set.seed(1)
n = 500
x = rexp(n)
e = runif(n,min=-2,max=2)
y = 3 + 2*x + e
linmod = lm(y~x)
linmod$coef
yhat = linmod$fitted.values
bsamp = bootsamp(linmod$residuals)
bsamp = matrix(yhat,n,ncol(bsamp))+bsamp
myfun = function(y,x) lm(y~x)$coef
bse = bootse(bsamp,myfun,x=x)
bse$cov
sigsq = mean(linmod$residuals^2)
solve(crossprod(cbind(1,x))) * sigsq
dev.new(height=6,width=4,noRStudioGD=TRUE)
par(mfcol=c(2,1))
hist(bse$theta[1,],main=expression(hat(b)[0]))
hist(bse$theta[2,],main=expression(hat(b)[1]))

# bootstrap regression pairs
set.seed(1)
n = 500
x = rexp(n)
e = runif(n,min=-2,max=2)
y = 3 + 2*x + e
linmod = lm(y~x)
linmod$coef
z = cbind(y,x)
bsamp = bootsamp(z)
myfun = function(z) lm(z[,1]~z[,2])$coef
bse = bootse(bsamp,myfun)
bse$cov
sigsq = mean(linmod$residuals^2)
solve(crossprod(cbind(1,x)))*sigsq
dev.new(height=6,width=4,noRStudioGD=TRUE)
par(mfcol=c(2,1))
hist(bse$theta[1,],main=expression(hat(b)[0]))
hist(bse$theta[2,],main=expression(hat(b)[1]))

# bootstrap robust regression pairs
set.seed(1)
n = 500
x = rnorm(n)
e = rnorm(n,sd=x^2)
y = 3 + 2*x + e
linmod = lm(y~x)
linmod$coef
rlinmod = rlm(y~x)
rlinmod$coef
z = cbind(y,x)
bsamp = bootsamp(z)
myfun = function(z) rlm(z[,1]~z[,2])$coef
bse = bootse(bsamp,myfun)
bse$cov
sigsq = mean(rlinmod$residuals^2)
solve(crossprod(cbind(1,x)))*sigsq
dev.new(height=6,width=4,noRStudioGD=TRUE)
par(mfcol=c(2,1))
hist(bse$theta[1,],main=expression(hat(b)[0]))
hist(bse$theta[2,],main=expression(hat(b)[1]))

