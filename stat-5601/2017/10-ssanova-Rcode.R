##########   Smoothing Spline ANOVA
##########   Nathaniel E. Helwig (helwig@umn.edu)
##########   Updated: 04-Jan-2017


#####   DEFINE PATHS AND PACKAGES   #####

# load bigsplines package (for various functions)
if(!require(bigsplines)){
  install.packages("bigsplines")
  library("bigsplines")
}

# load gss package (for ssanova function)
if(!require(gss)){
  install.packages("gss")
  library("gss")
}


#####   INTRODUCTION   #####

# example 1: continuous and nominal
addfun = function(x1,x2){
  funval = sin(2*pi*x1)
  idx = which(x2=="a")
  funval[idx] = funval[idx] + 2
  funval
}
intfun = function(x1,x2){
  funval = sin(2*pi*x1)
  idx = which(x2=="a")
  funval[idx] = funval[idx] + 2 + sin(4*pi*x1[idx])
  funval
}
dev.new(width=12,height=6,noRStudioGD=TRUE)
par(mfrow=c(1,2))
x1 = seq(0,1,length=200)
plot(x1,addfun(x1,rep("a",200)),type="l",ylim=c(-2,4),main="Additive",
     ylab="y",cex.axis=1.25,cex.lab=1.5,cex.main=3)
lines(x1,addfun(x1,rep("b",200)),lty=2)
legend("bottomleft",legend=c(expression(x[2]*" = "*a),expression(x[2]*" = "*b)),
       lty=1:2,bty="n",cex=1.5)
plot(x1,intfun(x1,rep("a",200)),type="l",ylim=c(-2,4),main="Interaction",
     ylab="y",cex.axis=1.25,cex.lab=1.5,cex.main=3)
lines(x1,intfun(x1,rep("b",200)),lty=2)
legend("bottomleft",legend=c(expression(x[2]*" = "*a),expression(x[2]*" = "*b)),
       lty=1:2,bty="n",cex=1.5)

# example 2: continuous and continuous
addfun = function(x1,x2){
  sin(2*pi*x1) + cos(4*pi*x2*(1-x2))
}
intfun = function(x1,x2){
  sin(2*pi*x1) + cos(4*pi*x2*(1-x2)) + 2*sin(pi*(x1-x2))
}
xs = seq(0,1,length=50)
xg = expand.grid(xs,xs)
dev.new(width=12,height=6,noRStudioGD=TRUE)
par(mfrow=c(1,2))
zmat = matrix(addfun(xg[,1],xg[,2]),50,50)
image(xs,xs,zmat,xlab="x1",ylab="x2",main="Additive",
      cex.axis=1.25,cex.lab=1.5,cex.main=3)
zmat = matrix(intfun(xg[,1],xg[,2]),50,50)
image(xs,xs,zmat,xlab="x1",ylab="x2",main="Interaction",
      cex.axis=1.25,cex.lab=1.5,cex.main=3)


# natural spline: set-up plot
dev.new(width=12,height=4,noRStudioGD=TRUE)
par(mfrow=c(1,3))

# natural spline: no noise
set.seed(1)
x = seq(0,1,length=21)
y = sin(2*pi*x)
mysp = spline(x,y,method="natural")
plot(x,y,main="No Noise",ylim=c(-1.5,1.5))
lines(x,sin(2*pi*x),lty=2)
lines(mysp)

# natural spline: some noise
set.seed(1)
x = seq(0,1,length=21)
y = sin(2*pi*x) + rnorm(21,sd=0.15)
mysp = spline(x,y,method="natural")
plot(x,y,main="Some Noise",ylim=c(-1.5,1.5))
lines(x,sin(2*pi*x),lty=2)
lines(mysp)

# natural spline: more noise
set.seed(1)
x = seq(0,1,length=21)
y = sin(2*pi*x) + rnorm(21,sd=0.25)
mysp = spline(x,y,method="natural")
plot(x,y,main="More Noise",ylim=c(-1.5,1.5))
lines(x,sin(2*pi*x),lty=2)
lines(mysp)

# smoothing spline: set-up plot
dev.new(width=12,height=4,noRStudioGD=TRUE)
par(mfrow=c(1,3))

# smoothing spline: no noise
set.seed(1)
x = seq(0,1,length=21)
y = sin(2*pi*x)
mysp = smooth.spline(x,y)
plot(x,y,main="No Noise",ylim=c(-1.5,1.5))
lines(x,sin(2*pi*x),lty=2)
lines(x,mysp$y)

# smoothing spline: some noise
set.seed(1)
x = seq(0,1,length=21)
y = sin(2*pi*x) + rnorm(21,sd=0.15)
mysp = smooth.spline(x,y)
plot(x,y,main="Some Noise",ylim=c(-1.5,1.5))
lines(x,sin(2*pi*x),lty=2)
lines(x,mysp$y)

# smoothing spline: more noise
set.seed(1)
x = seq(0,1,length=21)
y = sin(2*pi*x) + rnorm(21,sd=0.25)
mysp = smooth.spline(x,y)
plot(x,y,main="More Noise",ylim=c(-1.5,1.5))
lines(x,sin(2*pi*x),lty=2)
lines(x,mysp$y)


#####   SSANOVA IN PRACTICE   #####

##*##   one-way SSANOVA   ##*##

## smooth.spline overview

# smooth.spline: default
dev.new(width=6,height=6,noRStudioGD=TRUE)
set.seed(1)
x = seq(0,1,length=100)
eta = 2 + x + sin(2*pi*x)
y = eta + rnorm(100)
plot(x,y)
smsp = smooth.spline(x,y)
lines(x,smsp$y)
lines(x,eta,lty=2)

# smooth.spline: changing spar
dev.new(width=12,height=4,noRStudioGD=TRUE)
par(mfrow=c(1,3))
plot(x,y,main="spar=0.25")
smsp = smooth.spline(x,y,spar=0.25)
lines(x,smsp$y)
lines(x,eta,lty=2)
plot(x,y,main="spar=0.75")
smsp = smooth.spline(x,y,spar=0.75)
lines(x,smsp$y)
lines(x,eta,lty=2)
plot(x,y,main="spar=1")
smsp = smooth.spline(x,y,spar=1)
lines(x,smsp$y)
lines(x,eta,lty=2)

# smooth.spline: changing nknots (fixed spar)
dev.new(width=12,height=4,noRStudioGD=TRUE)
par(mfrow=c(1,3))
plot(x,y,main="spar=0.5, nknots=10")
smsp = smooth.spline(x,y,spar=0.5,nknots=10)
lines(x,smsp$y)
lines(x,eta,lty=2)
plot(x,y,main="spar=0.5, nknots=20")
smsp = smooth.spline(x,y,spar=0.5,nknots=20)
lines(x,smsp$y)
lines(x,eta,lty=2)
plot(x,y,main="spar=0.5, nknots=30")
smsp = smooth.spline(x,y,spar=0.5,nknots=30)
lines(x,smsp$y)
lines(x,eta,lty=2)

# smooth.spline: changing cv
dev.new(width=12,height=6,noRStudioGD=TRUE)
par(mfrow=c(1,2))
plot(x,y,main="nknots=20, cv=TRUE")
smsp = smooth.spline(x,y,nknots=20,cv=TRUE)
lines(x,smsp$y)
lines(x,eta,lty=2)
plot(x,y,main="nknots=20, cv=FALSE")
smsp = smooth.spline(x,y,nknots=20)
lines(x,smsp$y)
lines(x,eta,lty=2)

# smooth.spline: changing nknots (gcv spar)
dev.new(width=12,height=4,noRStudioGD=TRUE)
par(mfrow=c(1,3))
plot(x,y,main="cv=FALSE, nknots=10")
smsp = smooth.spline(x,y,nknots=10)
lines(x,smsp$y)
lines(x,eta,lty=2)
plot(x,y,main="cv=FALSE, nknots=20")
smsp = smooth.spline(x,y,nknots=20)
lines(x,smsp$y)
lines(x,eta,lty=2)
plot(x,y,main="cv=FALSE, nknots=30")
smsp = smooth.spline(x,y,nknots=30)
lines(x,smsp$y)
lines(x,eta,lty=2)

# smooth.spline: prediction
dev.new(width=6,height=6,noRStudioGD=TRUE)
plot(x,y,main="Prediction")
smsp = smooth.spline(x,y)
newdata = seq(0,1,length=200)
yhat = predict(smsp,newdata)
lines(yhat)
lines(x,eta,lty=2)


## bigspline overview

# bigspline: default
dev.new(width=6,height=6,noRStudioGD=TRUE)
set.seed(1)
x = seq(0,1,length=100)
eta = 2 + x + sin(2*pi*x)
y = eta + rnorm(100)
plot(x,y)
bigsp=bigspline(x,y)
lines(x,bigsp$fitted)
lines(x,eta,lty=2)

# bigspline: changing lambdas
dev.new(width=12,height=4,noRStudioGD=TRUE)
par(mfrow=c(1,3))
plot(x,y,main="lambdas=10^-9")
bigsp = bigspline(x,y,lambdas=10^-9)
lines(x,bigsp$fitted)
lines(x,eta,lty=2)
plot(x,y,main="lambdas=10^-5")
bigsp = bigspline(x,y,lambdas=10^-5)
lines(x,bigsp$fitted)
lines(x,eta,lty=2)
plot(x,y,main="lambdas=1")
bigsp = bigspline(x,y,lambdas=1)
lines(x,bigsp$fitted)
lines(x,eta,lty=2)

# bigspline: changing nknots (fixed lambda)
dev.new(width=12,height=4,noRStudioGD=TRUE)
par(mfrow=c(1,3))
plot(x,y,main="lambdas=10^-5, nknots=10")
bigsp = bigspline(x,y,lambdas=10^-5,nknots=10)
lines(x,bigsp$fitted)
lines(x,eta,lty=2)
plot(x,y,main="lambdas=10^-5, nknots=20")
bigsp = bigspline(x,y,lambdas=10^-5,nknots=20)
lines(x,bigsp$fitted)
lines(x,eta,lty=2)
plot(x,y,main="lambdas=10^-5, nknots=30")
bigsp = bigspline(x,y,lambdas=10^-5,nknots=30)
lines(x,bigsp$fitted)
lines(x,eta,lty=2)

# bigspline: changing nknots (gcv lambda)
dev.new(width=12,height=4,noRStudioGD=TRUE)
par(mfrow=c(1,3))
plot(x,y,main="GCV, nknots=10")
bigsp = bigspline(x,y,nknots=10)
lines(x,bigsp$fitted)
lines(x,eta,lty=2)
plot(x,y,main="GCV, nknots=20")
bigsp = bigspline(x,y,nknots=20)
lines(x,bigsp$fitted)
lines(x,eta,lty=2)
plot(x,y,main="GCV, nknots=30")
bigsp = bigspline(x,y,nknots=30)
lines(x,bigsp$fitted)
lines(x,eta,lty=2)

# bigspline: prediction
dev.new(width=6,height=6,noRStudioGD=TRUE)
plot(x,y,main="Prediction")
bigsp = bigspline(x,y)
newdata = seq(0,1,length=200)
yhat = predict(bigsp,newdata)
lines(newdata,yhat)
lines(x,eta,lty=2)

# bigspline: linear and non-linear
dev.new(width=12,height=4,noRStudioGD=TRUE)
par(mfrow=c(1,3))
bigsp = bigspline(x,y)
newdata = seq(0,1,length=200)
plot(x,y,main="Full Prediction")
yhat = predict(bigsp,newdata)
lines(newdata,yhat)
lines(x,eta,lty=2)
plot(x,2+x,main="Linear Effect",type="l",lty=2)
yhat = predict(bigsp,newdata,effect="0")+predict(bigsp,newdata,effect="lin")
lines(newdata,yhat)
plot(x,sin(2*pi*x),main="Non-Linear Effect",type="l",lty=2)
yhat = predict(bigsp,newdata,effect="non")
lines(newdata,yhat)

# bigspline: Bayesian CIs
dev.new(width=6,height=6,noRStudioGD=TRUE)
set.seed(1)
x = seq(0,1,length=100)
eta = 2 + x + sin(2*pi*x)
y = eta + rnorm(100)
bigsp = bigspline(x,y,se.fit=TRUE)
cilo = bigsp$fitted - qnorm(0.975)*bigsp$se.fit
cihi = bigsp$fitted + qnorm(0.975)*bigsp$se.fit
plot(x,y)
lines(x,eta)
lines(bigsp$xunique,cilo,lty=2)
lines(bigsp$xunique,cihi,lty=2)
sum(eta>=cilo & eta<=cihi)/length(x)


## smooth.spline versus bigspline

# mini simulation
nsamp = 10^c(2:6)
simresults = NULL
xnew = seq(0,1,length=200)
set.seed(1)
for(j in 1:5){
  for(k in 1:10){
    x = seq(0,1,length=nsamp[j])
    eta = 2 + x + sin(2*pi*x)
    y = eta + rnorm(nsamp[j])
    
    tic = proc.time()
    ssmod = smooth.spline(x,y,nknots=20)
    sstoc = proc.time() - tic
    tmse = sum( (ssmod$y - eta)^2 ) / nsamp[j]
    simsp = data.frame(method="smsp",n=nsamp[j],time=sstoc[3],tmse=tmse,row.names=k)
    
    tic = proc.time()
    ssmod = bigspline(x,y,nknots=20)
    sstoc = proc.time() - tic
    tmse = sum( (predict(ssmod) - eta)^2 ) / nsamp[j]
    simbig = data.frame(method="big",n=nsamp[j],time=sstoc[3],tmse=tmse,row.names=k+1)
    
    simresults = rbind(simresults,simsp,simbig)
  }
}

round(tapply(simresults$tmse,list(simresults$method,simresults$n),median),5)
round(tapply(simresults$time,list(simresults$method,simresults$n),median),5)

# bigspline: linear and non-linear effects (revisited)
dev.new(width=12,height=6,noRStudioGD=TRUE)
par(mfrow=c(2,4))
nsamp = c(100,1000,10000,100000)
xnew = seq(0,1,length=200)
set.seed(1)
for(j in 1:4){
  x = seq(0,1,length=nsamp[j])
  eta = 2 + x + sin(2*pi*x)
  y = eta + rnorm(nsamp[j])
  ssbig = bigspline(x,y,nknots=20)
  plot(xnew,2+xnew,type="l",lty=2,main=bquote("Linear: "*n==.(nsamp[j])))
  lines(xnew,predict(ssbig,effect="0",newdata=xnew)+predict(ssbig,effect="lin",newdata=xnew))
  plot(xnew,sin(2*pi*xnew),type="l",lty=2,main=bquote("Non-linear: "*n==.(nsamp[j])))
  lines(xnew,predict(ssbig,effect="non",newdata=xnew))
}



##*##   two-way SSANOVA   ##*##

## additive function

# define additive function
addfun = function(x1,x2){
  funval = sin(2*pi*x1)
  idx = which(x2=="a")
  funval[idx] = funval[idx] + 2
  funval
}

# plot additive function
dev.new(width=6,height=6,noRStudioGD=TRUE)
xseq = seq(0,1,length=100)
plot(xseq,addfun(xseq,rep("a",100)),type="l",ylim=c(-2,4),
     xlab=expression(x[1]),ylab=expression(eta(x[1],x[2])))
lines(xseq,addfun(xseq,rep("b",100)),lty=2)
legend("bottomleft",c(expression(x[2]*" = a"),expression(x[2]*" = b")),lty=1:2,bty="n")
dev.copy2pdf(file="/Users/Nate/Dropbox/Classes/MN_AP/2015b_Fall/class_notes/notes10_ssa/figs/addfun.pdf")

# bigssa model fitting
n = 100
set.seed(55455)
x1v = seq(0,1,length=n)
x2v = factor(sample(letters[1:2],n,replace=TRUE))
eta = addfun(x1v,x2v)
y = eta + rnorm(n)
idx = binsamp(cbind(x1v,x2v),nmbin=c(20,2))
ssint = bigssa(y~x1v*x2v, type=list(x1v="cub",x2v="nom"), nknots=idx)
sum((ssint$fitted-eta)^2) / length(eta)
ssadd = bigssa(y~x1v+x2v, type=list(x1v="cub",x2v="nom"), nknots=idx)
sum((ssadd$fitted-eta)^2) / length(eta)
fitstats = rbind(ssint$info,ssadd$info)
rownames(fitstats)=c("int","add")
fitstats

# bigssa model predictions
dev.new(width=12,height=6,noRStudioGD=TRUE)
par(mfrow=c(1,2))
newdata = expand.grid(x1v=seq(0,1,length=100),x2v=c("a","b"))
yint = predict(ssint,newdata)
yadd = predict(ssadd,newdata)
plot(newdata[1:100,1],yint[1:100],main="Interaction",type="l",ylim=c(-2,4))
lines(newdata[101:200,1],yint[101:200],lty=2)
plot(newdata[1:100,1],yadd[1:100],main="Additive",type="l",ylim=c(-2,4))
lines(newdata[101:200,1],yadd[101:200],lty=2)

# ssanova model fitting
n = 100
set.seed(55455)
x1v = seq(0,1,length=n)
x2v = factor(sample(letters[1:2],n,replace=TRUE))
eta = addfun(x1v,x2v)
y = eta + rnorm(n)
idx = binsamp(cbind(x1v,x2v),nmbin=c(20,2))
ssint = ssanova(y~x1v*x2v,type=list(x1v="cubic",x2v="nominal"),id.basis=idx)
newdata = data.frame(x1v=x1v,x2v=x2v)
sum((predict(ssint,newdata)-eta)^2) / length(eta)
ssadd = ssanova(y~x1v+x2v,type=list(x1v="cubic",x2v="nominal"),id.basis=idx)
sum((predict(ssadd,newdata)-eta)^2) / length(eta)

# ssanova model predictions
dev.new(width=12,height=6,noRStudioGD=TRUE)
par(mfrow=c(1,2))
newdata = expand.grid(x1v=seq(0,1,length=100),x2v=c("a","b"))
yint = predict(ssint,newdata)
yadd = predict(ssadd,newdata)
plot(newdata[1:100,1],yint[1:100],main="Interaction",type="l",ylim=c(-2,4))
lines(newdata[101:200,1],yint[101:200],lty=2)
plot(newdata[1:100,1],yadd[1:100],main="Additive",type="l",ylim=c(-2,4))
lines(newdata[101:200,1],yadd[101:200],lty=2)


## interaction function

# define interaction function
intfun = function(x1,x2){
  funval = sin(2*pi*x1)
  idx = which(x2=="a")
  funval[idx] = funval[idx] + 2 + sin(4*pi*x1[idx])
  funval
}

# plot interaction function
dev.new(width=6,height=6,noRStudioGD=TRUE)
xseq = seq(0,1,length=100)
plot(xseq,intfun(xseq,rep("a",100)),type="l",ylim=c(-2,4),
     xlab=expression(x[1]),ylab=expression(eta(x[1],x[2])))
lines(xseq,intfun(xseq,rep("b",100)),lty=2)
legend("bottomleft",c(expression(x[2]*" = a"),expression(x[2]*" = b")),lty=1:2,bty="n")

# bigssa model fitting
n = 100
set.seed(55455)
x1v = seq(0,1,length=n)
x2v = factor(sample(letters[1:2],n,replace=TRUE))
eta = intfun(x1v,x2v)
y = eta + rnorm(n)
idx = binsamp(cbind(x1v,x2v),nmbin=c(20,2))
ssint = bigssa(y~x1v*x2v,type=list(x1v="cub",x2v="nom"),nknots=idx)
sum((ssint$fitted-eta)^2) / length(eta)
ssadd = bigssa(y~x1v+x2v,type=list(x1v="cub",x2v="nom"),nknots=idx)
sum((ssadd$fitted-eta)^2) / length(eta)
fitstats = rbind(ssint$info,ssadd$info)
rownames(fitstats) = c("int","add")
fitstats

# bigssa model predictions
dev.new(width=12,height=6,noRStudioGD=TRUE)
par(mfrow=c(1,2))
newdata = expand.grid(x1v=seq(0,1,length=100),x2v=c("a","b"))
yint = predict(ssint,newdata)
yadd = predict(ssadd,newdata)
plot(newdata[1:100,1],yint[1:100],main="Interaction",type="l",ylim=c(-2,4))
lines(newdata[101:200,1],yint[101:200],lty=2)
plot(newdata[1:100,1],yadd[1:100],main="Additive",type="l",ylim=c(-2,4))
lines(newdata[101:200,1],yadd[101:200],lty=2)

# ssanova model fitting
n = 100
set.seed(55455)
x1v = seq(0,1,length=n)
x2v = factor(sample(letters[1:2],n,replace=TRUE))
eta = intfun(x1v,x2v)
y = eta + rnorm(n)
idx = binsamp(cbind(x1v,x2v),nmbin=c(20,2))
ssint = ssanova(y~x1v*x2v,type=list(x1v="cubic",x2v="nominal"),id.basis=idx)
newdata = data.frame(x1v=x1v,x2v=x2v)
sum((predict(ssint,newdata)-eta)^2) / length(eta)
ssadd = ssanova(y~x1v+x2v,type=list(x1v="cubic",x2v="nominal"),id.basis=idx)
sum((predict(ssadd,newdata)-eta)^2) / length(eta)

# ssanova model predictions
dev.new(width=12,height=6,noRStudioGD=TRUE)
par(mfrow=c(1,2))
newdata = expand.grid(x1v=seq(0,1,length=100),x2v=c("a","b"))
yint = predict(ssint,newdata)
yadd = predict(ssadd,newdata)
plot(newdata[1:100,1],yint[1:100],main="Interaction",type="l",ylim=c(-2,4))
lines(newdata[101:200,1],yint[101:200],lty=2)
plot(newdata[1:100,1],yadd[1:100],main="Additive",type="l",ylim=c(-2,4))
lines(newdata[101:200,1],yadd[101:200],lty=2)

# bigssa model fitting (more data)
n = 1000
set.seed(55455)
x1v = seq(0,1,length=n)
x2v = factor(sample(letters[1:2],n,replace=TRUE))
eta = intfun(x1v,x2v)
y = eta + rnorm(n)
idx = binsamp(cbind(x1v,x2v),nmbin=c(20,2))
ssint = bigssa(y~x1v*x2v,type=list(x1v="cub",x2v="nom"),nknots=idx)
sum((ssint$fitted-eta)^2) / length(eta)
ssadd = bigssa(y~x1v+x2v,type=list(x1v="cub",x2v="nom"),nknots=idx)
sum((ssadd$fitted-eta)^2) / length(eta)
fitstats = rbind(ssint$info,ssadd$info)
rownames(fitstats) = c("int","add")
fitstats

# bigssa model predictions (more data)
dev.new(width=12,height=6,noRStudioGD=TRUE)
par(mfrow=c(1,2))
newdata = expand.grid(x1v=seq(0,1,length=100),x2v=c("a","b"))
yint = predict(ssint,newdata)
yadd = predict(ssadd,newdata)
plot(newdata[1:100,1],yint[1:100],main="Interaction",type="l",ylim=c(-2,4))
lines(newdata[101:200,1],yint[101:200],lty=2)
plot(newdata[1:100,1],yadd[1:100],main="Additive",type="l",ylim=c(-2,4))
lines(newdata[101:200,1],yadd[101:200],lty=2)

