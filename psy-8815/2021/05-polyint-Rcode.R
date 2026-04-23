##########   Polynomials and Interactions
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


#####   POLYNOMIAL REGRESSION   #####

# SLR when relationship is quadratic
dev.new(width=6,height=6,noRStudioGD=TRUE)
x=seq(-1,1,l=50)
y=2+2*(x^2)
qmod=lm(y~x)
plot(x,y,main="Quadratic")
abline(qmod)

# SLR when relationship is cubic
dev.new(width=6,height=6,noRStudioGD=TRUE)
x=seq(-1,1,l=50)
y=2+2*(x^3)
qmod=lm(y~x)
plot(x,y,main="Cubic")
abline(qmod)

# multicollinearity: problem
set.seed(123)
x=runif(100)*2
X=cbind(x,xsq=x^2,xcu=x^3)
cor(X)

# multicollinearity: partial solution
set.seed(123)
x=runif(100)*2
x=x-mean(x)
X=cbind(x,xsq=x^2,xcu=x^3)
cor(X)

# Gram-Schmidt process
orthog<-function(X,normalize=FALSE){
  np=dim(X)
  Z=matrix(0,np[1],np[2])
  Z[,1]=X[,1]
  for(k in 2:np[2]){
    Z[,k]=X[,k]
    for(j in 1:(k-1)){
      Z[,k]=Z[,k]-Z[,j]*sum(Z[,k]*Z[,j])/sum(Z[,j]^2)
    }
  }
  if(normalize){Z=Z%*%diag(colSums(Z^2)^-0.5)}
  Z
}

# Gram-Schmidt example
set.seed(123)
X=cbind(1,runif(10),runif(10))
crossprod(X)
Z=orthog(X)
crossprod(Z)
Z=orthog(X,norm=TRUE)
crossprod(Z)

# compare to poly function
set.seed(123)
x=runif(10)
X=cbind(1,x,xsq=x^2,xcu=x^3)
Z=orthog(X,norm=TRUE)
z=poly(x,degree=3)
Z[,2:4]=Z[,2:4]%*%diag(colSums(z^2)^0.5)
Z[1:3,]
cbind(Z[1:3,1],z[1:3,])

# cars example
cars=read.table(paste(datapath,"auto-mpg.data",sep=""),na.strings="?")
names(cars)=c("mpg","cylinder","disp","hp","weight","accel","year","origin","name")
idx=which(is.na(cars$hp))
mpg=cars$mpg[-idx]
hp=cars$hp[-idx]

# SLR model
dev.new(width=6,height=6,noRStudioGD=TRUE)
plot(hp,mpg)
linmod=lm(mpg~hp)
abline(linmod)

# quadratic regression model
dev.new(width=6,height=6,noRStudioGD=TRUE)
quadmod=lm(mpg~hp+I(hp^2))
quadmod
hpseq=seq(50,250,by=5)
Xmat=cbind(1,hpseq,hpseq^2)
hphat=Xmat%*%quadmod$coef
plot(hp,mpg)
lines(hpseq,hphat)

# cubic regression model
cubmod=lm(mpg~hp+I(hp^2)+I(hp^3))
summary(cubmod)

# orthogonal vs raw polynomials
quadomod=lm(mpg~poly(hp,degree=2))
summary(quadomod)$coef
summary(quadmod)$coef

# try our orthog function
X=orthog(cbind(1,hp,hp^2))
quadomod=lm(mpg~X[,2]+X[,3])
summary(quadomod)$coef


#####   HOUSE PRICE EXAMPLE   #####

# read-in house price data
house=read.table(paste(datapath,"houseprice.txt",sep=""),header=TRUE)
house[1:8,]

# plot raw data
dev.new(width=6,height=6,noRStudioGD=TRUE)
plot(house$value,house$price,pch=ifelse(house$corner==1,0,16),
     xlab="Appraisal Value",ylab="Selling Price")
legend("topleft",c("Corner","Non-Corner"),pch=c(0,16),bty="n")

# fit interaction model
hmod=lm(price~value*corner,data=house)
summary(hmod)
hmod$coef

# plot regression lines
dev.new(width=6,height=6,noRStudioGD=TRUE)
plot(house$value,house$price,pch=ifelse(house$corner==1,0,16),
     xlab="Appraisal Value",ylab="Selling Price")
abline(hmod$coef[1],hmod$coef[2])
abline(hmod$coef[1]+hmod$coef[3],hmod$coef[2]+hmod$coef[4],lty=2)
legend("topleft",c("Corner","Non-Corner"),lty=2:1,pch=c(0,16),bty="n")



#####   DEPRESSION EXAMPLE   #####

# read-in depression data
depression=read.table(paste(datapath,"depression.txt",sep=""),header=TRUE)
depression[1:4,]

# plot raw data
dev.new(width=6,height=6,noRStudioGD=TRUE)
plot(depression$age,depression$effect,xlab="Age",ylab="Effect",type="n")
text(depression$age,depression$effect,depression$method)

# fit interaction model
dmod=lm(effect~age*method,data=depression)
dmod$coef
summary(dmod)

# plot regression lines
dev.new(width=6,height=6,noRStudioGD=TRUE)
plot(depression$age,depression$effect,xlab="Age",ylab="Effect",type="n")
text(depression$age,depression$effect,depression$method)
abline(dmod$coef[1],dmod$coef[2])
abline(dmod$coef[1]+dmod$coef[3],dmod$coef[2]+dmod$coef[5],lty=2)
abline(dmod$coef[1]+dmod$coef[4],dmod$coef[2]+dmod$coef[6],lty=3)
legend("bottomright",c("A","B","C"),lty=1:3,bty="n")


#####   PLOT REGRESSION SURFACES   #####

# define simple, multiple (additive), and multiple (interaction)
my1fun<-function(x){2+2*x}
my2fun<-function(x,z){2+2*x-2*z}
my3fun<-function(x,z){2+2*x-2*z+x*z/4}

# plot simple regression
x=seq(0,10,l=100)
y1=my1fun(x)
dev.new(width=6,height=6,noRStudioGD=TRUE)
par(mar=c(5,5,4,2))
plot(x,y1,ylab=expression(2+2*italic(x)[i]),cex.axis=1.5,cex.lab=2,cex.main=2,
     xlab=expression(italic(x)[i]),main="Simple regression")

# plot multiple regression (additive)
mydata=expand.grid(x=seq(0,10,l=30),z=seq(0,5,l=30))
y2=my2fun(mydata[,1],mydata[,2])
dev.new(width=6,height=6,noRStudioGD=TRUE)
par(mar=c(5,5,4,2))
wireframe(y2~mydata[,2]*mydata[,1],xlab=list(label=expression(italic(x)[i2]),cex=2),
          ylab=list(label=expression(italic(x)[i1]),cex=2),
          zlab=list(label=expression(2+2*italic(x)[i1]-2*italic(x)[i2]),cex=.5),
          scales=list(arrows=FALSE,cex=1.5),
          main=list(label="Multiple regression (additive)",cex=2),
          zlim=c(-10,25),screen=list(z=45,x=-60))

# plot multiple regression (interaction)
mydata=expand.grid(x=seq(0,10,l=30),z=seq(0,5,l=30))
y3=my3fun(mydata[,1],mydata[,2])
dev.new(width=6,height=6,noRStudioGD=TRUE)
par(mar=c(5,5,4,2))
wireframe(y3~mydata[,2]*mydata[,1],xlab=list(label=expression(italic(x)[i2]),cex=2),
          ylab=list(label=expression(italic(x)[i1]),cex=2),
          zlab=list(label=expression(2+2*italic(x)[i1]-2*italic(x)[i2]+italic(x)[i1]*italic(x)[i2]/4),cex=.5),scales=list(arrows=FALSE,cex=1.5),
          main=list(label="Multiple regression (interaction)",cex=2),
          zlim=c(-10,25),screen=list(z=45,x=-60))



#####   EL NINO EXAMPLE   #####

# load data ("El Nino" data from http://archive.ics.uci.edu/ml/datasets.html)
myclass=c(rep("integer",5),rep("numeric",2),rep("character",5))
elnino=read.table(paste(datapath,"tao-all2.dat",sep=""),colClasses=myclass)
colnames(elnino)=c("obs","year","month","day","date","latitude","longitude",
                   "zon.winds","mer.winds","humidity","air.temp","ss.temp")

# subset data
elnino=elnino[elnino$ss.temp!=".",]   # missing sea surface temperature
ilon=which(elnino$long<0)             # negative longitudes
elnino$longitude[ilon]=elnino$longitude[ilon]+360  # convert to positive longitudes
elnino=elnino[elnino$long>=156 & elnino$year>=94,] # removes far West points and data before 1994
elnino$ss.temp=as.numeric(elnino$ss.temp)

# look at data
elnino[1:4,]

# fit additive model
eladd=lm(ss.temp~latitude+longitude,data=elnino)
summary(eladd)

# get fitted values for additive model
newdata=expand.grid(longitude=seq(min(elnino$longitude),max(elnino$longitude),length=50),
                    latitude=seq(min(elnino$latitude),max(elnino$latitude),length=50))
yadd=predict(eladd,newdata)

# plot regression surface for additive model
dev.new(width=6,height=6,noRStudioGD=TRUE)
image(seq(min(elnino$longitude),max(elnino$longitude),l=50),
      seq(min(elnino$latitude),max(elnino$latitude),l=50),
      matrix(yadd,50,50),col=rev(rainbow(100,end=3/4)),
      xlab="Longitude",ylab="Latitude",main="Additve Prediction")

# fit interaction model
elint=lm(ss.temp~latitude*longitude,data=elnino)
summary(elint)

# get predictions for intercation model
newdata=expand.grid(longitude=seq(min(elnino$longitude),max(elnino$longitude),length=50),
                    latitude=seq(min(elnino$latitude),max(elnino$latitude),length=50))
yint=predict(elint,newdata)

# plot regression surface for interaction model
dev.new(width=6,height=6,noRStudioGD=TRUE)
image(seq(min(elnino$longitude),max(elnino$longitude),l=50),
      seq(min(elnino$latitude),max(elnino$latitude),l=50),
      matrix(yint,50,50),col=rev(rainbow(100,end=3/4)),
      xlab="Longitude",ylab="Latitude",main="Interaction Prediction")

