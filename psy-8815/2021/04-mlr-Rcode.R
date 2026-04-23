##########   Multiple Linear Regression
##########   Nathaniel E. Helwig (helwig@umn.edu)
##########   Updated: 04-Jan-2017


#####   DEFINE PATHS AND PACKAGES   #####

# define data path (to load data)
datapath = "~/Desktop/notes/data/"

# load lattice package (for wireframe function)
if(!require(lattice)){
  install.packages("lattice")
  library(lattice)
}

# load car package (for confidence ellipses)
if(!require(car)){
  install.packages("car")
  library("car")
}


#####   OVERVIEW of MLR MODEL   #####

# plot regression surface
dev.new(height=6,width=6,noRStudioGD=TRUE)
x1=seq(0,10,length.out=50)
x2=seq(0,5,length.out=50)
mydata=expand.grid(x1,x2)
y=2+2*mydata[,1]-2*mydata[,2]
wireframe(y~mydata[,2]*mydata[,1],xlab=list(label=expression(italic(x)[i2]),cex=2),
          ylab=list(label=expression(italic(x)[i1]),cex=2),
          zlab=list(label=expression(2+2*italic(x)[i1]-2*italic(x)[i2]),cex=2,rot=90,vjust=0),
          scales=list(arrows=FALSE,cex=1.5),
          main=list(label="Multiple regression surface",cex=2,vjust=2),
          zlim=c(-10,25),screen=list(z=45,x=-60),
          par.settings=list(axis.line=list(col="transparent")))

# read-in GPA data
gpa=read.table(paste(datapath,"sat.csv",sep=""),header=TRUE,sep=",")
gpa[1:10,]

# define variables
uGPA=gpa$univ_GPA
hGPA=gpa$high_GPA
vSAT=gpa$verb_SAT/100

# fit model
gpamod=lm(uGPA~hGPA+vSAT)
summary(gpamod)

# plot points
dev.new(height=6,width=6,noRStudioGD=TRUE)
cloud(uGPA~vSAT*hGPA,xlab=list(label="Verbal SAT",cex=2,rot=55),
      ylab=list(label="HS GPA",cex=2,rot=-15),
      zlab=list(label="Univ GPA",cex=2,rot=90,vjust=0),
      scales=list(arrows=FALSE,cex=1.5),
      main=list(label="3D Scatterplot",cex=2,vjust=2),
      zlim=c(2,4),screen=list(z=65,x=-60),
      par.settings=list(axis.line=list(col="transparent")))



#####   INFERENCES IN MLR   #####

# plot confidence regions
dev.new(height=4,width=8,noRStudioGD=TRUE)
par(mfrow=c(1,2))
confidenceEllipse(gpamod,c(2,3),levels=.9,xlim=c(.2,1),ylim=c(-.1,.4),
                  main=expression(alpha*" = "*.1),cex.main=2)
confidenceEllipse(gpamod,c(2,3),levels=.99,xlim=c(.2,1),ylim=c(-.1,.4),
                  main=expression(alpha*" = "*.01),cex.main=2)

# rescale SAT scores
summary(gpa[,1:3])
gpa[,2:3]=gpa[,2:3]/100
summary(gpa[,1:3])

# fit full model
gpaFmod=lm(univ_GPA~high_GPA+verb_SAT+math_SAT,data=gpa)
summary(gpaFmod)

# reduced model (dropping math SAT)
gpaRmod=update(gpaFmod,~.-math_SAT)
summary(gpaRmod)

# ANOVA table comparing models
anova(gpaRmod,gpaFmod)

# ANOVA table with type I SS (oreder matters!!)
anova(gpaRmod)
gpa2mod=lm(univ_GPA~verb_SAT+high_GPA,data=gpa)
anova(gpa2mod)

# testing b1=b2
xvar=gpa$high_GPA+gpa$verb_SAT
gpaEmod=lm(univ_GPA~xvar,data=gpa)
anova(gpaEmod,gpaRmod)

# testing b0=b1
high_GPA1p=1+gpa$high_GPA
gpaImod=lm(univ_GPA~0+high_GPA1p+verb_SAT,data=gpa)
gpaImod$coef
gpaRmod$coef
anova(gpaImod,gpaRmod)

# testing b1=3*b2
wvar=gpa$high_GPA+gpa$verb_SAT/3
gpaLmod=lm(univ_GPA~wvar,data=gpa)
anova(gpaLmod,gpaRmod)

# get information for reduced (R) model
sumRmod=summary(gpaRmod)
sumRmod$coef
sumRmod$sigma
sumRmod$sigma^2
sumRmod$r.squared
sumRmod$adj.r.squared

# mean center GPA and SAT
mean(gpa$high_GPA)
hGPA=gpa$high_GPA-mean(gpa$high_GPA)
mean(gpa$verb_SAT)
vSAT=gpa$verb_SAT-mean(gpa$verb_SAT)
gpaSmod=lm(gpa$univ_GPA~hGPA+vSAT)
summary(gpaSmod)$coef

# model F: manual calculations
XF=cbind(1,gpa$high_GPA,gpa$verb_SAT,gpa$math_SAT)
y=gpa$univ_GPA
XtXF=crossprod(XF)
XtyF=crossprod(XF,y)
XtXiF=solve(XtXF)
bhatF=XtXiF%*%XtyF
yhatF=XF%*%bhatF
ehatF=y-yhatF
sigsqF=sum(ehatF^2)/(nrow(XF)-ncol(XF))
bhatseF=sqrt(sigsqF*diag(XtXiF))
tvalF=bhatF/bhatseF
pvalF=2*(1-pt(abs(tvalF),nrow(XF)-ncol(XF)))
RsqF=1-sum(ehatF^2)/sum((y-mean(y))^2)
aRsqF=1-(sum(ehatF^2)/(nrow(XF)-ncol(XF)))/(sum((y-mean(y))^2)/(nrow(XF)-1))
data.frame(bhat=bhatF,se=bhatseF,t=tvalF,p=pvalF)
cbind(RsqF,aRsqF)

# model R: manual calculations
XR=cbind(1,gpa$high_GPA,gpa$verb_SAT)
y=gpa$univ_GPA
XtXR=crossprod(XR)
XtyR=crossprod(XR,y)
XtXiR=solve(XtXR)
bhatR=XtXiR%*%XtyR
yhatR=XR%*%bhatR
ehatR=y-yhatR
sigsqR=sum(ehatR^2)/(nrow(XR)-ncol(XR))
bhatseR=sqrt(sigsqR*diag(XtXiR))
tvalR=bhatR/bhatseR
pvalR=2*(1-pt(abs(tvalR),nrow(XR)-ncol(XR)))
RsqR=1-sum(ehatR^2)/sum((y-mean(y))^2)
aRsqR=1-(sum(ehatR^2)/(nrow(XR)-ncol(XR)))/(sum((y-mean(y))^2)/(nrow(XR)-1))
data.frame(bhat=bhatR,se=bhatseR,t=tvalR,p=pvalR)
cbind(RsqR,aRsqR)

# model E: manual calculations
XE=cbind(1,gpa$high_GPA+gpa$verb_SAT)
y=gpa$univ_GPA
XtXE=crossprod(XE)
XtyE=crossprod(XE,y)
XtXiE=solve(XtXE)
bhatE=XtXiE%*%XtyE
yhatE=XE%*%bhatE
ehatE=y-yhatE
sigsqE=sum(ehatE^2)/(nrow(XE)-ncol(XE))
bhatseE=sqrt(sigsqE*diag(XtXiE))
tvalE=bhatE/bhatseE
pvalE=2*(1-pt(abs(tvalE),nrow(XE)-ncol(XE)))
RsqE=1-sum(ehatE^2)/sum((y-mean(y))^2)
aRsqE=1-(sum(ehatE^2)/(nrow(XE)-ncol(XE)))/(sum((y-mean(y))^2)/(nrow(XE)-1))
data.frame(bhat=bhatE,se=bhatseE,t=tvalE,p=pvalE)
cbind(RsqE,aRsqE)

# model I: manual calculations
XI=cbind(1+gpa$high_GPA,gpa$verb_SAT)
y=gpa$univ_GPA
XtXI=crossprod(XI)
XtyI=crossprod(XI,y)
XtXiI=solve(XtXI)
bhatI=XtXiI%*%XtyI
yhatI=XI%*%bhatI
ehatI=y-yhatI
sigsqI=sum(ehatI^2)/(nrow(XI)-ncol(XI))
bhatseI=sqrt(sigsqI*diag(XtXiI))
tvalI=bhatI/bhatseI
pvalI=2*(1-pt(abs(tvalI),nrow(XI)-ncol(XI)))
RsqI=1-sum(ehatI^2)/sum((y-mean(y))^2)
aRsqI=1-(sum(ehatI^2)/(nrow(XI)-ncol(XI)))/(sum((y-mean(y))^2)/(nrow(XI)-1))
data.frame(bhat=bhatI,se=bhatseI,t=tvalI,p=pvalI)
cbind(RsqI,aRsqI)

# model L: manual calculations
XL=cbind(1,gpa$high_GPA+gpa$verb_SAT/3)
y=gpa$univ_GPA
XtXL=crossprod(XL)
XtyL=crossprod(XL,y)
XtXiL=solve(XtXL)
bhatL=XtXiL%*%XtyL
yhatL=XL%*%bhatL
ehatL=y-yhatL
sigsqL=sum(ehatL^2)/(nrow(XL)-ncol(XL))
bhatseL=sqrt(sigsqL*diag(XtXiL))
tvalL=bhatL/bhatseL
pvalL=2*(1-pt(abs(tvalL),nrow(XL)-ncol(XL)))
RsqL=1-sum(ehatL^2)/sum((y-mean(y))^2)
aRsqL=1-(sum(ehatL^2)/(nrow(XL)-ncol(XL)))/(sum((y-mean(y))^2)/(nrow(XL)-1))
data.frame(bhat=bhatL,se=bhatseL,t=tvalL,p=pvalL)
cbind(RsqL,aRsqL)

