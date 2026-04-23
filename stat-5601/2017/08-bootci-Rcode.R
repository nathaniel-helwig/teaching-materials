##########   Bootstrap Confidence Intervals
##########   Nathaniel E. Helwig (helwig@umn.edu)
##########   Updated: 04-Jan-2017


#####   DEFINE PATHS AND PACKAGES   #####

# load bootstrap package (for various functions)
if(!require(bootstrap)){
  install.packages("bootstrap")
  library("bootstrap")
}


###*###   RESAMPLING FUNCTIONS   ###*###

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


#####   CONFIDENCE INTERVALS   #####

# proper interpretation: example
set.seed(1)
n = 100
B = 10^4
X = replicate(B, rnorm(n))
xbar = apply(X, 2, mean)
xsd = apply(X, 2, sd)
cilo = xbar - qt(.95, df=n-1)*(xsd/sqrt(n))
ciup = xbar - qt(.05, df=n-1)*(xsd/sqrt(n))
ci90 = (0>=cilo & 0<=ciup)
mean(ci90)
summary(ci90)


#####   BASIC BOOTSTRAP CIs   #####

# ex 1: mean ci
quartz(width=5,height=5)
set.seed(1)
n = 50
x = rnorm(n,mean=1)
bsamp = bootsamp(x)
bse = bootse(bsamp,mean)
mean(x)
c(mean(x)-qt(0.975,df=n-1)*sd(x)/sqrt(n),
  mean(x)-qt(0.025,df=n-1)*sd(x)/sqrt(n))
ci = quantile(bse$theta,c(0.025,0.975))
ci
hist(bse$theta)
lines(rep(ci[1],2),c(0,1500),lty=2)
lines(rep(ci[2],2),c(0,1500),lty=2)
standT = c(mean(x)-qt(0.975,df=n-1)*sd(x)/sqrt(n),
           mean(x)-qt(0.025,df=n-1)*sd(x)/sqrt(n))
prcnt.mean = ci

# ex 2: median ci
quartz(width=5,height=5)
set.seed(1)
n = 50
x = rnorm(n,mean=1)
bsamp = bootsamp(x)
bse = bootse(bsamp,median)
median(x)
ci = quantile(bse$theta,c(0.025,0.975))
ci
hist(bse$theta)
lines(rep(ci[1],2),c(0,4000),lty=2)
lines(rep(ci[2],2),c(0,4000),lty=2)
prcnt.med = ci

# ex 3: variance ci
quartz(width=5,height=5)
set.seed(1)
n = 50
x = rnorm(n,sd=2)
bsamp = bootsamp(x)
bse = bootse(bsamp,var)
var(x)
c((n-1)*var(x)/qchisq(0.975,df=n-1),
  (n-1)*var(x)/qchisq(0.025,df=n-1))
ci = quantile(bse$theta,c(0.025,0.975))
ci
hist(bse$theta)
lines(rep(ci[1],2),c(0,2000),lty=2)
lines(rep(ci[2],2),c(0,2000),lty=2)
prcnt.var = ci
standZ.var = c((n-1)*var(x)/qchisq(0.975,df=n-1),
               (n-1)*var(x)/qchisq(0.025,df=n-1))

# transformation respecting
set.seed(1)
n = 50
x = rnorm(n,mean=1)
bsamp = bootsamp(x)
bse = bootse(bsamp,function(x) exp(mean(x)))
exp(mean(x))
mean(x)
ci = quantile(bse$theta,c(0.025,0.975))
ci
bse = bootse(bsamp,mean)
quantile(bse$theta,c(0.025,0.975))
log(ci)


#####   BETTER BOOTSTRAP CIs   #####

# expanded percentile example 1 (mean)
set.seed(1)
n = 50
x = rnorm(n,mean=1)
bsamp = bootsamp(x)
bse = bootse(bsamp,mean)
mean(x)
c(mean(x)-qnorm(0.975)*sd(x)*sqrt((n-1)/n)/sqrt(n),
  mean(x)-qnorm(0.025)*sd(x)*sqrt((n-1)/n)/sqrt(n))
quantile(bse$theta,c(0.025,0.975))
alphaD2 = pnorm(sqrt(n/(n-1))*qt(.025,df=n-1))
alphaD2
c(mean(x)-qt(0.975,df=n-1)*sd(x)/sqrt(n),
  mean(x)-qt(0.025,df=n-1)*sd(x)/sqrt(n))
quantile(bse$theta,c(alphaD2,1-alphaD2))
standZ = c(mean(x)-qnorm(0.975)*sd(x)*sqrt((n-1)/n)/sqrt(n),
           mean(x)-qnorm(0.025)*sd(x)*sqrt((n-1)/n)/sqrt(n))
eprcnt.mean = quantile(bse$theta,c(alphaD2,1-alphaD2))

# expanded percentile example 2 (median)
set.seed(1)
n = 50
x = rnorm(n,mean=1)
bsamp = bootsamp(x)
bse = bootse(bsamp,median)
median(x)
quantile(bse$theta,c(0.025,0.975))
alphaD2 = pnorm(sqrt(n/(n-1))*qt(.025,df=n-1))
alphaD2
quantile(bse$theta,c(alphaD2,1-alphaD2))
eprcnt.med = quantile(bse$theta,c(alphaD2,1-alphaD2))

# expanded percentile example 3 (variance)
set.seed(1)
n = 50
x = rnorm(n,sd=2)
bsamp = bootsamp(x)
bse = bootse(bsamp,var)
var(x)
quantile(bse$theta,c(0.025,0.975))
alphaD2 = pnorm(sqrt(n/(n-1))*qt(.025,df=n-1))
alphaD2
quantile(bse$theta,c(alphaD2,1-alphaD2))
eprcnt.var = quantile(bse$theta,c(alphaD2,1-alphaD2))

# bootstrap t-table example 1 (mean)
set.seed(1)
n = 50
x = rnorm(n,mean=1)
bsamp = bootsamp(x)
bse = bootse(bsamp,mean)
theta = mean(x)
bsampSE = apply(bsamp, 2, sd) * sqrt((n-1)/n) / sqrt(n)
Z = (bse$theta - theta) / bsampSE
cval = quantile(Z, probs=c(0.025,0.975))
# bootstrap t-table
c(theta - cval[2]*bse$se, theta - cval[1]*bse$se)
bootT.mean = c(theta - cval[2]*bse$se, theta - cval[1]*bse$se)
# percentile
quantile(bse$theta,c(0.025,0.975))

# bootstrap t-table example 2 (median)
set.seed(1)
n = 50
x = rnorm(n,mean=1)
bsamp = bootsamp(x)
bse = bootse(bsamp,median)
theta = median(x)
nsamp = ncol(bsamp)
bsampSE = rep(0, nsamp)
for(k in 1:nsamp){
  cat("samp:",k,"\n")
  bsampSE[k] = bootse(bootsamp(bsamp[,k],nsamp=2000),median)$se
}
Z = (bse$theta - theta) / bsampSE
cval = quantile(Z, probs=c(0.025,0.975))
# bootstrap t-table
c(theta - cval[2]*bse$se, theta - cval[1]*bse$se)
bootT.med = c(theta - cval[2]*bse$se, theta - cval[1]*bse$se)
# percentile
quantile(bse$theta,c(0.025,0.975))

# bootstrap t-table example 3 (variance)
set.seed(1)
n = 50
x = rnorm(n,sd=2)
bsamp = bootsamp(x)
bse = bootse(bsamp,var)
theta = var(x)
nsamp = ncol(bsamp)
bsampSE = rep(0, nsamp)
for(k in 1:nsamp){
  cat("samp:",k,"\n")
  bsampSE[k] = bootse(bootsamp(bsamp[,k],nsamp=2000),var)$se
}
Z = (bse$theta - theta) / bsampSE
cval = quantile(Z, probs=c(0.025,0.975))
# bootstrap t-table
c(theta - cval[2]*bse$se, theta - cval[1]*bse$se)
bootT.var = c(theta - cval[2]*bse$se, theta - cval[1]*bse$se)
# percentile
quantile(bse$theta,c(0.025,0.975))

# BCa confidence interval function
bcafun <- function(x,nboot,theta,...,alpha=0.05){
  theta.hat = theta(x)
  nx = length(x)
  bse = bootse(bootsamp(x,nboot),theta,...)
  jse = jackse(jacksamp(x),theta,...)
  z0 = qnorm(sum(bse$theta<theta.hat)/nboot)
  atop = sum((mean(jse$theta)-jse$theta)^3)
  abot = 6*(((jse$se^2)*nx/(nx-1))^(3/2))
  ahat = atop/abot
  alpha1 = pnorm(z0+(z0+qnorm(alpha))/(1-ahat*(z0+qnorm(alpha))))
  alpha2 = pnorm(z0+(z0+qnorm(1-alpha))/(1-ahat*(z0+qnorm(1-alpha))))
  confpoint = quantile(bse$theta,probs=c(alpha1,alpha2))
  list(confpoint=confpoint,z0=z0,acc=ahat,u=(jse$theta-theta.hat),
       theta=bse$theta,se=bse$se)
}

# BCa example 1: mean
set.seed(1)
n = 50
x = rnorm(n,mean=1)
c(mean(x)-qt(0.975,df=n-1)*sd(x)/sqrt(n),
  mean(x)-qt(0.025,df=n-1)*sd(x)/sqrt(n))
mybca = bcafun(x,10000,mean,alpha=0.025)
quantile(mybca$theta,probs=c(0.025,0.975))
mybca$conf
bca = bcanon(x,10000,mean,alpha=c(0.025,0.975))
bca$conf
bca.mean = mybca$conf

# BCa example 2: median
set.seed(1)
n = 50
x = rnorm(n,mean=1)
mybca = bcafun(x,10000,median,alpha=0.025)
quantile(mybca$theta,c(0.025,0.975))
mybca$conf
bca = bcanon(x,10000,median,alpha=c(0.025,0.975))
bca$conf
bca.med = mybca$conf

# BCa example 3: variance
set.seed(1)
n = 50
x = rnorm(n,sd=2)
mybca = bcafun(x,10000,var,alpha=0.025)
quantile(mybca$theta,c(0.025,0.975))
mybca$conf
bca = bcanon(x,10000,var,alpha=c(0.025,0.975))
bca$conf
bca.var = mybca$conf

# table results for example 1 (mean)
tab.mean = rbind(standZ,standT,prcnt.mean,
                 eprcnt.mean,bootT.mean,bca.mean)
rownames(tab.mean) = c("standZ","standT","prcnt",
                       "eprcnt","bootT","bca")
round(tab.mean,4)

# table results for example 2 (median)
tab.med = rbind(prcnt.med,eprcnt.med,bootT.med,bca.med)
rownames(tab.med) = c("prcnt","eprcnt","bootT","bca")
round(tab.med,4)

# table results for example 3 (variance)
tab.var = rbind(standZ.var,prcnt.var,eprcnt.var,bootT.var,bca.var)
rownames(tab.var) = c("standZ","prcnt","eprcnt","bootT","bca")
round(tab.var,4)
