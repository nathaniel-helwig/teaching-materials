##########   Analysis of Covariance
##########   Nathaniel E. Helwig (helwig@umn.edu)
##########   Updated: 04-Jan-2017


#####   DEFINE PATHS AND PACKAGES   #####

# define data path (to load data)
datapath = "~/Desktop/notes/data/"

# load car package (for Levene Test)
if(!require(car)){
  install.packages("car")
  library("car")
}

# load multcomp package (for glht function)
if(!require(multcomp)){
  install.packages("multcomp")
  library("multcomp")
}

#####   OVERVIEW OF DATA   #####

# load audit data
audit=read.table(paste(datapath,"audit.txt",sep=""),header=TRUE)
head(audit)

# look at pre- and post-test means by group
tapply(audit$pre,audit$method,mean)
tapply(audit$post,audit$method,mean)

# mean-center covariate
audit=cbind(audit,cpre=(audit$pre-mean(audit$pre)))

# plot data with overall least-squares line
dev.new(width=6,heigh=6,noRStudioGD=TRUE)
plot(audit$pre,audit$post,type="n",xlab="Pre-test Audit Proficiency",ylab="Post-test Audit Proficiency")
text(audit$pre,audit$post,audit$method)
abline(lm(post~pre,data=audit))

# plot data with method-specific least-squares lines
dev.new(width=6,heigh=6,noRStudioGD=TRUE)
mycol=rep("blue",30)
mycol[audit$method=="A"]="red"
mycol[audit$method=="C"]="green3"
plot(audit$pre,audit$post,type="n",xlab="Pre-test Audit Proficiency",ylab="Post-test Audit Proficiency")
text(audit$pre,audit$post,audit$method,col=mycol)
mycol=c("red","blue","green3")
mymod=c("A","B","C")
for(k in 1:3){abline(lm(post~pre,data=subset(audit,method==mymod[k])),col=mycol[k])}


#####   TESTING ASSUMPTIONS   #####

# test linear relationship between pre and post
rmod=lm(post~pre,data=audit)
rmod$coef
anova(rmod)

# test parallel slopes assumption
contrasts(audit$method)=contr.sum(3)
imod=lm(post~pre*method,data=audit)
anova(imod)

# test relationship between treatment and covariate
anova(lm(pre~method,data=audit))

# initial test of homogeneity of variance
leveneTest(post~method,data=audit)


#####   FITTING ANOVA   #####

# fit one-way ANOVA model
amod=lm(post~method,data=audit)
summary(amod)$sigma
anova(amod)

# Tukey's HSD
TukeyHSD(aov(post~method,data=audit))


#####   FITTING ANCOVA   #####

# fit uncentered ANCOVA model
umod=lm(post~pre+method,data=audit)
umod$coef
summary(umod)$sigma
anova(umod)

# fit centered ANCOVA model
cmod=lm(post~cpre+method,data=audit)
cmod$coef
summary(cmod)$sigma
anova(cmod)


#####   MULTIPLE COMPARISONS   #####

# calculations by hand
postmean=tapply(audit$post,audit$method,mean)
cpremean=tapply(audit$cpre,audit$method,mean)
postmean - cmod$coef[2]*cpremean
newdata=data.frame(method=c("A","B","C"),cpre=rep(0,3))
yhat=predict(cmod,newdata=newdata)
yhat
BmA=yhat[2]-yhat[1]
CmA=yhat[3]-yhat[1]
CmB=yhat[3]-yhat[2]
difs=c(BmA,CmA,CmB)
names(difs)=c("B-A","C-A","C-B")
difs

# comparing treatment effects
pwc=glht(cmod,linfct=mcp(method="Tukey"))
summary(pwc)
