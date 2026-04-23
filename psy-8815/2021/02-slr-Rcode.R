##########   Simple Linear Regression
##########   Nathaniel E. Helwig (helwig@umn.edu)
##########   Updated: 04-Jan-2017


#####   DEFINE PATHS AND PACKAGES   #####

# define data path (to load data)
datapath = "~/Desktop/notes/data/"


#####   OVERVIEW of SLR MODEL   #####

# simple linear regression visualization
dev.new(height=5,width=8,noRStudioGD=TRUE)
set.seed(123)
x=1:10
e=rnorm(10,sd=2)
y=1.1+2*x+e
plot(x,y,xlim=c(0,11),ylim=c(0,22),cex=1.5,cex.lab=1.5,cex.axis=1.5)
abline(a=1.1,b=2,lwd=2)
for(i in 1:10){lines(rep(x[i],2),c(y[i],1.1+2*x[i]),lty=2,lwd=1.5)}
legend("topleft",c(expression(italic(E(y))*" = 1.1 + 2"*italic(x)),
                   expression(italic(e[i]))),lty=c(1,2),cex=2,lwd=rep(1.5,2),bty="n")



#####   DRINKING EXAMPLE   #####

# read-in drinking data
drinking=read.table(paste(datapath,"drinking.txt",sep=""),header=TRUE)
drinking

# fit simple linear regression model
drinkmod=lm(cirrhosis~alcohol,data=drinking)
summary(drinkmod)

# recreate summary table by hand
X=cbind(1,drinking$alcohol)
y=drinking$cirrhosis
XtX=crossprod(X)
Xty=crossprod(X,y)
XtXi=solve(XtX)
bhat=XtXi%*%Xty
yhat=X%*%bhat
ehat=y-yhat
sigsq=sum(ehat^2)/(nrow(X)-2)
bhatse=sqrt(sigsq*diag(XtXi))
tval=bhat/bhatse
pval=2*(1-pt(abs(tval),nrow(X)-2))
data.frame(bhat=bhat,se=bhatse,t=tval,p=pval)

# plot regression line
dev.new(height=5,width=8,noRStudioGD=TRUE)
plot(drinking$alcohol,drinking$cirrhosis,type="n",
     xlab="yearly alcohol (liters/person)",ylab="cirrhosis deaths (per 100,000)")
text(drinking$alcohol,drinking$cirrhosis,drinking$country)
abline(drinkmod$coef[1],drinkmod$coef[2])

# read-in new data
drinknew=read.table(paste(datapath,"drinknew.txt",sep=""),header=TRUE)
drinknew

# get predictions, CI, and PI
predict(drinkmod,newdata=drinknew)
predict(drinkmod,newdata=drinknew,interval="confidence",level=.9)
predict(drinkmod,newdata=drinknew,interval="prediction",level=.9)

# plot line with CI, confidence bound (CB), and PI
drng=range(drinking$alcohol)
drinkseq=data.frame(alcohol=seq(drng[1],drng[2],length.out=100))
civals=predict(drinkmod,newdata=drinkseq,interval="confidence")
pivals=predict(drinkmod,newdata=drinkseq,interval="prediction")
sevals=predict(drinkmod,newdata=drinkseq,se.fit=T)
dev.new(height=5,width=8,noRStudioGD=TRUE)
plot(drinking$alcohol,drinking$cirrhosis,ylim=c(-5,55),
     xlab="yearly alcohol (liters/person)",
     ylab="cirrhosis deaths (per 100,000)")
abline(drinkmod$coef[1],drinkmod$coef[2])
W=sqrt(2*qf(.95,2,13))
lines(drinkseq$alcohol,civals[,2],lty=2,col="blue",lwd=2)
lines(drinkseq$alcohol,civals[,3],lty=2,col="blue",lwd=2)
lines(drinkseq$alcohol,sevals$fit+W*sevals$se.fit,lty=3,col="green3",lwd=2)
lines(drinkseq$alcohol,sevals$fit-W*sevals$se.fit,lty=3,col="green3",lwd=2)
lines(drinkseq$alcohol,pivals[,2],lty=4,col="red",lwd=2)
lines(drinkseq$alcohol,pivals[,3],lty=4,col="red",lwd=2)
legend("topleft",c("CI","CB","PI"),lty=2:4,cex=2,
       lwd=rep(2,3),col=c("blue","green3","red"),bty="n")


#####   GPA EXAMPLE   #####

# read-in GPA data
gpa=read.table(paste(datapath,"sat.csv",sep=""),header=TRUE,sep=",")
gpa[1:10,]

# fit model predicting university GPA from high school GPA
gpamod=lm(univ_GPA~high_GPA,data=gpa)
summary(gpamod)

# recreate summary table by hand
X=cbind(1,gpa$high_GPA)
y=gpa$univ_GPA
XtX=crossprod(X)
Xty=crossprod(X,y)
XtXi=solve(XtX)
bhat=XtXi%*%Xty
yhat=X%*%bhat
ehat=y-yhat
sigsq=sum(ehat^2)/(nrow(X)-2)
bhatse=sqrt(sigsq*diag(XtXi))
tval=bhat/bhatse
pval=2*(1-pt(abs(tval),nrow(X)-2))
data.frame(bhat=bhat,se=bhatse,t=tval,p=pval)

# plot regression line
dev.new(height=5,width=8,noRStudioGD=TRUE)
par(mar=c(5,5.4,4,2)+0.1)
plot(gpa$high_GPA,gpa$univ_GPA,xlab="High School GPA",
     ylab="University GPA",cex.lab=2,cex.axis=2)
abline(a=gpamod$coef[1],gpamod$coef[2])

# predict for new hypothetical data
gpanew=data.frame(high_GPA=c(2.4,3,3.1,3.3,3.9),
                  univ_GPA=rep(NA,5))
gpanew
predict(gpamod,newdata=gpanew)

# plot line with CI, confidence bound (CB), and PI
drng=range(gpa$high_GPA)
gpaseq=data.frame(high_GPA=seq(drng[1],drng[2],length.out=100))
civals=predict(gpamod,newdata=gpaseq,interval="confidence")
pivals=predict(gpamod,newdata=gpaseq,interval="prediction")
sevals=predict(gpamod,newdata=gpaseq,se.fit=T)
dev.new(height=5,width=8,noRStudioGD=TRUE)
plot(gpa$high_GPA,gpa$univ_GPA,ylim=c(2,5),
     xlab="High School GPA",
     ylab="University GPA")
abline(gpamod$coef[1],gpamod$coef[2])
W=sqrt(2*qf(.95,2,103))
lines(gpaseq$high_GPA,civals[,2],lty=2,col="blue",lwd=2)
lines(gpaseq$high_GPA,civals[,3],lty=2,col="blue",lwd=2)
lines(gpaseq$high_GPA,sevals$fit+W*sevals$se.fit,
      lty=3,col="green3",lwd=2)
lines(gpaseq$high_GPA,sevals$fit-W*sevals$se.fit,
      lty=3,col="green3",lwd=2)
lines(gpaseq$high_GPA,pivals[,2],lty=4,col="red",lwd=2)
lines(gpaseq$high_GPA,pivals[,3],lty=4,col="red",lwd=2)
legend("topleft",c("CI","CB","PI"),lty=2:4,cex=2,
       lwd=rep(2,3),col=c("blue","green3","red"),bty="n")
