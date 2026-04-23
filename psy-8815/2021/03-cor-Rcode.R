##########   Correlation and Geometry
##########   Nathaniel E. Helwig (helwig@umn.edu)
##########   Updated: 04-Jan-2017


#####   DEFINE PATHS AND PACKAGES   #####

# define data path (to load data)
datapath = "~/Desktop/notes/data/"


#####   PEARSON'S CORRELATION   #####

# function to simulate correlation sampling distribution
rsim<-function(rho,n){
  x=rnorm(n)
  y=rho*x+rnorm(n,sd=sqrt(1-rho^2))
  cor(x,y)
}

# simulate different n with rho=0
set.seed(1234)
r10z=replicate(10000,rsim(rho=0,n=10))
r30z=replicate(10000,rsim(rho=0,n=30))
r60z=replicate(10000,rsim(rho=0,n=60))
r90z=replicate(10000,rsim(rho=0,n=90))
r120z=replicate(10000,rsim(rho=0,n=120))
r150z=replicate(10000,rsim(rho=0,n=150))
dev.new(height=6,width=9,noRStudioGD=TRUE)
par(mfrow=c(2,3))
hist(r10z,main=expression(italic(n)*"="*10))
hist(r30z,main=expression(italic(n)*"="*30))
hist(r60z,main=expression(italic(n)*"="*60))
hist(r90z,main=expression(italic(n)*"="*90))
hist(r120z,main=expression(italic(n)*"="*120))
hist(r150z,main=expression(italic(n)*"="*150))

# simulate different n with rho!=0
set.seed(1234)
r1p=replicate(10000,rsim(rho=0.1,n=100))
r5p=replicate(10000,rsim(rho=0.5,n=100))
r9p=replicate(10000,rsim(rho=0.9,n=100))
r1n=replicate(10000,rsim(rho=-0.1,n=100))
r5n=replicate(10000,rsim(rho=-0.5,n=100))
r9n=replicate(10000,rsim(rho=-0.9,n=100))
dev.new(height=6,width=9,noRStudioGD=TRUE)
par(mfrow=c(2,3))
hist(r1p,main=expression(italic(n)*"="*100*",  "*italic(r)*"="*0.1))
hist(r5p,main=expression(italic(n)*"="*100*",  "*italic(r)*"="*0.5))
hist(r9p,main=expression(italic(n)*"="*100*",  "*italic(r)*"="*0.9))
hist(r1n,main=expression(italic(n)*"="*100*",  "*italic(r)*"="*-0.1))
hist(r5n,main=expression(italic(n)*"="*100*",  "*italic(r)*"="*-0.5))
hist(r9n,main=expression(italic(n)*"="*100*",  "*italic(r)*"="*-0.9))



#####   INFERENCES WITH CORRELATIONS   #####

# read-in GPA data
gpa=read.table(paste(datapath,"sat.csv",sep=""),header=TRUE,sep=",")
gpa[1:10,]

# calculate correlation using cor function
X=gpa$high_GPA
Y=gpa$univ_GPA
cor(X,Y)

# calculate correlation using cov and sd functions
cov(X,Y)/(sd(X)*sd(Y))

# calculate correlation long way
mux=mean(X)
muy=mean(Y)
cxy=sum((X-mux)*(Y-muy))
sx=sqrt(sum((X-mux)^2))
sy=sqrt(sum((Y-muy)^2))
cxy/(sx*sy)

# perform correlation test using cor.test function
cor.test(X,Y)

# perform correlation test by hand
gpacr=cor(X,Y)
tstar=gpacr*sqrt(length(X)-2)/sqrt(1-gpacr^2)
tstar
2*(1-pt(tstar,103))
z=log((1+gpacr)/(1-gpacr))/2
z
zlo=z-qnorm(.975)/sqrt(102)
zhi=z+qnorm(.975)/sqrt(102)
c(zlo,zhi)
rlo=(exp(2*zlo)-1)/(exp(2*zlo)+1)
rhi=(exp(2*zhi)-1)/(exp(2*zhi)+1)
c(rlo,rhi)

# Fisher-Z transformation
fisherz=function(r,n,rho0=0){
  z=log((1+r)/(1-r))/2
  z0=log((1+rho0)/(1-rho0))/2
  zstar=(z-z0)*sqrt(n-3)
  pval=2*(1-pnorm(abs(zstar)))
  list(z=z,pval=pval)
}
fisherz(cor(X,Y),105,rho0=0.7)



#####   GEOMETRICAL INTERPRETATIONS   #####

pcor=function(x,y,z,type=c("partial","part")){
  rxy=cor(x,y)
  rxz=cor(x,z)
  ryz=cor(y,z)
  pc=(rxy-ryz*rxz)/sqrt(1-rxz^2)
  if(type[1]=="partial"){pc=pc/sqrt(1-ryz^2)}
  pc
}
