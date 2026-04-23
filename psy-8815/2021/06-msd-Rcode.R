##########   Model Selection and Diagnostics
##########   Nathaniel E. Helwig (helwig@umn.edu)
##########   Updated: 04-Jan-2017


#####   DEFINE PATHS AND PACKAGES   #####

# load leaps package (for branch and bound regression)
if(!require(leaps)){
  install.packages("leaps")
  library(leaps)
}

# load faraway package (for variance inflation factor)
if(!require(faraway)){
  install.packages("faraway")
  library(faraway)
}


#####   MODEL SELECTION   #####

# load packages and define functions
getpress <- function(ix,y,x){
  if(any(ix)){
    linmod=lm(y~.,data=as.data.frame(x[,ix]))
  } else{
    linmod=lm(y~1)
  }
  sum((linmod$residuals/(1-hatvalues(linmod)))^2)
}
presslm <- function(x,y){
  x=as.data.frame(x)
  np=ncol(x)
  xlist=vector("list",np)
  for(j in 1:np){xlist[[j]]=c(TRUE,FALSE)}
  xall=expand.grid(xlist)
  allpress=apply(xall,1,getpress,y=y,x=x)
  list(which=as.matrix(xall),press=allpress)
}

# get states data
states=data.frame(state.x77,row.names=state.abb)
states[1:10,]

# full model predicting murder rate
fullmod=lm(Murder~.,data=states)
summary(fullmod)

# adjusted R^2 selection
X=states[,-5]
arsqmod=leaps(x=X,y=states$Murder,method="adjr2")
widx=which.max(arsqmod$adjr2)
xidx=(1:ncol(X))[arsqmod$which[widx,]]
Xin=data.frame(X[,xidx])
arsqmod=lm(states$Murder~.,data=Xin)
summary(arsqmod)

# stepwise AIC selection
smod=lm(states$Murder~.,data=states)
aicmod=step(smod,trace=0)
summary(aicmod)

# stepwise BIC selection
smod=lm(states$Murder~.,data=states)
bicmod=step(smod,k=log(50),trace=0)
summary(bicmod)

# Mallow's Cp selection
X=states[,-5]
cpmod=leaps(x=X,y=states$Murder,method="Cp")
widx=which.min(cpmod$Cp)
xidx=(1:ncol(X))[cpmod$which[widx,]]
Xin=data.frame(X[,xidx])
cpmod=lm(states$Murder~.,data=Xin)
summary(cpmod)

# PRESS selection
X=states[,-5]
prmod=presslm(x=X,y=states$Murder)
widx=which.min(prmod$press)
xidx=(1:ncol(X))[prmod$which[widx,]]
Xin=as.data.frame(X[,xidx])
prmod=lm(states$Murder~.,data=Xin)
summary(prmod)

# summarize results
xnames=colnames(states)[-5]
xtab=matrix(0,5,7)
rownames(xtab)=c("Ra^2","AIC","BIC","Cp","PRESS")
colnames(xtab)=xnames
xlist=list(arsqmod,aicmod,bicmod,cpmod,prmod)
for(j in 1:5){
  ix=match(names(attr(xlist[[j]]$terms,"dataClasses"))[-1],xnames)
  xtab[j,ix]=1
}
xtab


#####   MODEL DIAGNOSTICS   #####

# get states data
states=data.frame(state.x77,row.names=state.abb)
states[1:10,]

# fit model suggested by AIC (and Cp and Ra^2)
amod=lm(Murder~Population+Illiteracy+Life.Exp+Frost+Area,data=states)

# normality assumption
shapiro.test(amod$resid)
dev.new(width=9,height=4.5,noRStudioGD=TRUE)
par(mfrow=c(1,2))
qqnorm(amod$resid)
qqline(amod$resid)
hist(amod$resid,freq=F)
xseq=seq(-5,5,length=200)
lines(xseq,dnorm(xseq,sd=summary(amod)$sigma))

# linearity assumption
yhat=amod$fit
ehat=amod$resid
dev.new(width=6,height=6,noRStudioGD=TRUE)
plot(yhat,ehat,
     xlab=expression(hat(y)[i]),
     ylab=expression(hat(e)[i]),
     main="Residual Plot")
lines(range(yhat),c(0,0))
smod=smooth.spline(yhat,ehat)
lines(smod,col="blue")

# homogeneity of variance assumption
BPtest=function(mymod){
  mymod$model[,1]=(mymod$resid)^2
  newmod=lm(formula(mymod),data=mymod$model)
  modsum=summary(newmod)
  Rsq=modsum$r.squared
  BPstat=Rsq*(dim(mymod$model)[1])
  pval=1-pchisq(BPstat,modsum$df[1]-1)
  list(BP=BPstat,df=modsum$df[1]-1,pval=pval)
}
BPtest(amod)

# equal influence assumption
cookplot<-function(mymod,k=NULL,alpha=0.1,ptext=TRUE,...){
  nx=dim(mymod$model)[1]
  np=length(mymod$coef)
  cdist=cooks.distance(mymod)
  if(is.null(k)){k=qf(alpha,np,nx-np)} 
  ylim=range(cdist)
  if(ylim[1]>k){ylim[1]=k} else if(ylim[2]<k){ylim[2]=k}
  if(ptext){
    plot(1:nx,cdist,type="n",xlab=expression(italic(i)),ylim=ylim,
         ylab=expression(italic(D[i])),main="Cook's Distance Plot")
    text(1:nx,cdist,1:nx)
  } else{plot(1:nx,cdist,xlab=expression(italic(i)),ylim=ylim,
              ylab=expression(italic(D[i])),main="Cook's Distance Plot")}
  lines(c(1,nx),c(k,k),...)
}
dev.new(width=6,height=6,noRStudioGD=TRUE)
cookplot(amod)
rownames(states)[c(2,11,28)]

# multicollinearity
X=model.matrix(amod)[,-1]
X[1:4,]
vif(X)
