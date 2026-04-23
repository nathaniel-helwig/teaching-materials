##########   Effect Sizes and Power Analyses
##########   Nathaniel E. Helwig (helwig@umn.edu)
##########   Updated: 04-Jan-2017


#####   DEFINE PATHS AND PACKAGES   #####

# define data path (to load data)
datapath = "~/Desktop/notes/data/"

# load pwr package (for power calculations)
if(!require(pwr)){
  install.packages("pwr")
  library("pwr")
}


#####   EFFECT SIZES   #####

# example 1: one-way ANOVA
omega.sq <- function(mod){
  atab = anova(mod)
  ssb = atab[["Sum Sq"]][1]
  ssw = atab[["Sum Sq"]][2]
  dfb = atab[["Df"]][1]
  msw = atab[["Mean Sq"]][2]
  (ssb - dfb*msw) / (ssb + ssw + msw)
}
set.seed(1)
g = factor(sample(c(1,2,3),100,replace=TRUE))
e = rnorm(100)
mu = rbind(c(0,0.05,0.1),c(0,0.5,1),c(0,5,10))
eta = omega = rep(NA,3)
for(k in 1:3){
  y = 2 + mu[k,g] + e
  mod = lm(y~g)
  eta[k] = summary(mod)$r.squared
  omega[k] = omega.sq(mod)
}
eta
omega

# example 2: two-way ANOVA
eta.sq <- function(mod,k=NULL){
  atab = anova(mod)
  if(is.null(k)){k = 1:(nrow(atab)-1)}
  sum(atab[k,2]) / sum(atab[,2])
}
A = factor(rep(c("male","female"),each=12))
B = factor(rep(c("a","b","c"),8))
set.seed(1)
e = rnorm(24)
muA = c(0,2)
muB = c(0,1,2)
y = 2 + muA[A] + muB[B] + e
mod = aov(y~A+B)
eta.sq(mod)
eta.sq(mod,k=1)
eta.sq(mod,k=2)

# example 3: simple regression
set.seed(1)
x = rnorm(100)
e = rnorm(100)
bs = c(0.05,0.5,5)
R = Ra = rep(NA,3)
for(k in 1:3){
  y = 2 + bs[k]*x + e
  smod = summary(lm(y~x))
  R[k] = smod$r.squared
  Ra[k] = smod$adj.r.squared
}
R
Ra


# example 4: multiple regression
gpa = read.csv(paste(datapath,"sat.csv",sep=""),header=TRUE)

g1mod = lm(univ_GPA~high_GPA,data=gpa)
Rsq1 = summary(g1mod)$r.squared
g2mod = lm(univ_GPA~high_GPA+math_SAT,data=gpa)
Rsq2 = summary(g2mod)$r.squared
(Rsq2-Rsq1)/(1-Rsq2)   # f^2 (math_SAT given high_GPA)

g1mod = lm(univ_GPA~math_SAT,data=gpa)
Rsq1 = summary(g1mod)$r.squared
g2mod = lm(univ_GPA~math_SAT+high_GPA,data=gpa)
Rsq2 = summary(g2mod)$r.squared
(Rsq2-Rsq1)/(1-Rsq2)   # f^2 (high_GPA given math_SAT)


# example 5: t-test
diffES <- 
  function(x,y,type=c("gs","g","d","D")){
    md = mean(x) - mean(y)
    nx = length(x)
    ny = length(y)
    if(type[1]=="gs"){
      m = nx + ny - 2
      cm = gamma(m/2)/(sqrt(m/2)*gamma((m-1)/2))
      spsq = ((nx-1)*sd(x)+(ny-1)*sd(y))/m
      theta = cm*md/sqrt(spsq)
    } else if(type[1]=="g"){
      spsq = ((nx-1)*sd(x)+(ny-1)*sd(y))/(nx+ny-2)
      theta = md/sqrt(spsq)
    } else if(type[1]=="d"){
      spsq = ((nx-1)*sd(x)+(ny-1)*sd(y))/(nx+ny)
      theta = md/sqrt(spsq)
    } else { theta = md/sd(x) }
}
set.seed(1)
e = rnorm(100)
mu = rbind(c(0,0.05),c(0,0.5),c(0,1))
gs = g = d = D = rep(NA,3)
for(k in 1:3){
  x = rnorm(100,mean=mu[k,1])
  y = rnorm(100,mean=mu[k,2])
  gs[k] = diffES(x,y)
  g[k] = diffES(x,y,type="g")
  d[k] = diffES(x,y,type="d")
  D[k] = diffES(x,y,type="D")
}
rtab = rbind(gs,g,d,D)
rownames(rtab) = c("gs","g","d","D")
colnames(rtab) = c("small","medium","large")
rtab


# example 6: one-way ANOVA (revisited)
rmsse <- function(x,g){
  mx = tapply(x,g,mean)
  ng = nlevels(g)
  nx = length(x)
  msd = sum((mx-mean(x))^2)/(ng-1)
  mse = sum((mx[g]-x)^2)/(nx-ng)
  sqrt(msd/mse)
}
set.seed(1)
g = factor(sample(c(1,2,3),100,replace=TRUE))
e = rnorm(100)
mu = rbind(c(0,0.05,0.1),c(0,0.5,1),c(0,5,10))
rvec = rep(NA,3)
for(k in 1:3){
  y = 2 + mu[k,g] + e
  rvec[k] = rmsse(y,g)
}
rvec


#####   POWER ANALYSES   #####
plotpower <- function(m0,m1,sd=1,alpha=0.05){
  x0seq = seq(m0-3*sd,m0+3*sd,length=500)
  x1seq = seq(m1-3*sd,m1+3*sd,length=500)
  cval = qnorm(1-alpha,m0,sd)
  power = round(pnorm(cval,m1,sd,lower=F),2)
  plot(x0seq,dnorm(x0seq,m0,sd),xlim=c(m0-3*sd-1,m1+3*sd+1),type="l",lwd=2,
       xlab=expression(italic(x)),ylab=expression(italic(f(x))),
       main=bquote(1-beta==.(power)),cex.main=2,cex.axis=1.25,cex.lab=1.5)
  lines(x1seq,dnorm(x1seq,m1,sd),lwd=2)
  px=c(rep(cval,2),seq(cval-0.01,m1-3*sd,length=50),m1-3*sd,cval)
  py=c(0,dnorm(cval,m1,sd),dnorm(seq(cval-0.01,m1-3*sd,length=50),m1,sd),
       rep(dnorm(m1-3*sd,m1,sd),2))
  polygon(px,py,col="lightblue",border=NA)
  px=c(rep(cval,2),seq(cval+0.1,m0+3*sd,length=50),m0+3*sd,cval)
  py=c(0,dnorm(cval,m0,sd),dnorm(seq(cval+0.1,m0+3*sd,length=50),m0,sd),
       rep(dnorm(m0+3*sd,m0,sd),2))
  polygon(px,py,col="pink",border=NA)
  legend("topleft",legend=c(as.expression(bquote(alpha==.(alpha))),
                            as.expression(bquote(beta==.(1-power)))),
         fill=c("pink","lightblue"),bty="n")
}

dev.new(width=6,height=4,noRStudioGD=TRUE)
par(mar=c(5.1,5.1,4.1,2.1))
plotpower(0,2)
text(1.95,0.025,expression(alpha))
text(0.8,0.1,expression(beta))
text(2.75,0.1,expression(1-beta))
text(-1.25,0.35,expression(H[0]))
text(3.25,0.35,expression(H[1]))

dev.new(width=9,height=6,noRStudioGD=TRUE)
par(mfrow=c(2,3),mar=c(5.1,5.1,4.1,2.1))
for(j in seq(0.5,3,by=0.5)){
  plotpower(0,j)
}

dev.new(width=9,height=6,noRStudioGD=TRUE)
par(mfrow=c(2,3),mar=c(5.1,5.1,4.1,2.1))
for(j in 1:6){
  plotpower(0,1,sd=1/sqrt(j))
}

dev.new(width=9,height=6,noRStudioGD=TRUE)
par(mfrow=c(2,3),mar=c(5.1,5.1,4.1,2.1))
for(j in c(0.01,0.02,0.05,0.1,0.15,0.2)){
  plotpower(0,2,alpha=j)
}


# one-sample t-test
power.t.test(n=NULL,delta=-1,sd=1,sig.level=0.05,power=0.80,
             type="one.sample",alternative="one.sided")
power.t.test(n=NULL,delta=1,sd=1,sig.level=0.05,power=0.80,
             type="one.sample",alternative="one.sided")
power.t.test(n=NULL,delta=1,sd=1,sig.level=0.05,power=0.80,
             type="one.sample",alternative="two.sided")
power.t.test(n=NULL,delta=0.2,sd=1,sig.level=0.05,power=0.80,
             type="one.sample",alternative="two.sided")
power.t.test(n=NULL,delta=0.02,sd=1,sig.level=0.05,power=0.80,
             type="one.sample",alternative="two.sided")


# two-sample t-test
power.t.test(n=NULL,delta=-1,sd=1,sig.level=0.05,power=0.80,
             alternative="one.sided")
power.t.test(n=NULL,delta=1,sd=1,sig.level=0.05,power=0.80,
             alternative="one.sided")
power.t.test(n=NULL,delta=1,sd=1,sig.level=0.05,power=0.80,
             alternative="two.sided")
power.t.test(n=NULL,delta=0.2,sd=1,sig.level=0.05,power=0.80,
             alternative="two.sided")
power.t.test(n=NULL,delta=0.02,sd=1,sig.level=0.05,power=0.80,
             alternative="two.sided")


# one-way ANOVA
power.anova.test(groups=3,n=NULL,between.var=1,within.var=2,
                 sig.level=0.05,power=0.80)
power.anova.test(groups=3,n=NULL,between.var=1,within.var=4,
                 sig.level=0.05,power=0.80)
power.anova.test(groups=3,n=NULL,between.var=1,within.var=100,
                 sig.level=0.05,power=0.80)
power.anova.test(groups=3,n=NULL,between.var=1,within.var=1000,
                 sig.level=0.05,power=0.80)


# multiple regression (3 predictors)
pwr.f2.test(u=2,v=NULL,f2=4,sig.level=0.05,power=0.80)
pwr.f2.test(u=2,v=NULL,f2=1,sig.level=0.05,power=0.80)
pwr.f2.test(u=2,v=NULL,f2=0.1,sig.level=0.05,power=0.80)
pwr.f2.test(u=2,v=NULL,f2=0.01,sig.level=0.05,power=0.80)
pwr.f2.test(u=4,v=NULL,f2=0.01,sig.level=0.05,power=0.80)

