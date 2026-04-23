##########   Nonparametric Independence Tests
##########   Nathaniel E. Helwig (helwig@umn.edu)
##########   Updated: 04-Jan-2017


#####   DEFINE PATHS AND PACKAGES   #####

# load NSM3 package (for various functions)
if(!require(NSM3)){
  install.packages("NSM3")
  library("NSM3")
}


#####   KENDALL'S TAU   #####

# distribution of test statistic under H0
x = matrix(rep(1:3,6),6,3,byrow=T)
y = matrix(c(1,2,3,
             1,3,2,
             2,1,3,
             2,3,1,
             3,1,2,
             3,2,1),6,3,byrow=T)
K=rep(0,6)
for(k in 1:6){
  for(i in 1:2){
    for(j in (i+1):3){
      xx = ifelse((y[k,j]-y[k,i])*(x[k,j]-x[k,i])>0,1,-1)
      K[k] = K[k] + xx
    }
  }
}

# example data
x = c(277,169,157,139,108,213,232,229,114,232,161,149,128)
y = c(256,118,137,144,146,221,184,188,97,231,114,187,230)
#data.frame(and1=rep("&",13),x=x,and2=rep("&",13),y=y,and3=rep("\\",13))

# make Q for all pairs
n = 13
Qmat = matrix(0,n-1,n-1)
colnames(Qmat) = 1:(n-1)
rownames(Qmat) = 2:n
for(i in 1:(n-1)){
  for(j in (i+1):n){
    qval = (y[j]-y[i])*(x[j]-x[i])
    if(qval>0){
      Qmat[j-1,i] = 1
    } else if(qval<0){
      Qmat[j-1,i] = -1
    }
  }
}


K = sum(Qmat)    # or #    K = sum(Qmat[Qmat>0]) + sum(Qmat[Qmat<0])
tauhat = K/(n*(n-1)/2)
taub = K/sqrt(77*78)   # Kendall's tau-b
K
tauhat
taub

cor.test(x,y,method="kendall",alternative="greater")
cor(x,y,method="kendall")
kendall.ci(x,y)

# Kendall CI
Qfun = function(i,j){
  qval = (y[j]-y[i])*(x[j]-x[i])
  q = 0
  if(qval>0) { q = 1 } else if(qval<0) { q = -1 }
  return(q)
}
Cvec = rep(0,n)
idx = 1:n
for(i in idx){
  for(j in idx[idx!=i]){
    Cvec[i] = Cvec[i] + Qfun(i,j)
  }
}
Cbar = mean(Cvec)
const = 2/(n*(n-1))
sigsq = const*( const*((n-2)/(n-1))*sum((Cvec-Cbar)^2) + 1 - tauhat^2 )
c(tauhat-qnorm(.975)*sqrt(sigsq), tauhat+qnorm(.975)*sqrt(sigsq))
c(tauhat-qnorm(.95)*sqrt(sigsq), 1)

# easy way
x = c(277,169,157,139,108,213,232,229,114,232,161,149,128)
y = c(256,118,137,144,146,221,184,188,97,231,114,187,230)
cor(x,y,method="kendall")
cor.test(x,y,method="kendall",alternative="greater")
kendall.ci(x,y,type="t")
kendall.ci(x,y,type="l")


#####   SPEARMAN'S RHO   #####

# distribution of test statistic under H0
x = matrix(rep(1:3,6),6,3,byrow=T)
y = matrix(c(1,2,3,
             1,3,2,
             2,1,3,
             2,3,1,
             3,1,2,
             3,2,1),6,3,byrow=T)
rs = rep(0,6)
for(k in 1:6){
  rs[k] = cor(x[k,],y[k,],method="spearman")
}

# example data
x = c(277,169,157,139,108,213,232,229,114,232,161,149,128)
y = c(256,118,137,144,146,221,184,188,97,231,114,187,230)
#data.frame(and1=rep("&",13),x=x,and2=rep("&",13),rx=rank(x),and3=rep("&",13),y=y,and4=rep("&",13),ry=rank(y),and5=rep("\\",13))

# example:  by hand (with R)
rx = rank(x)
ry = rank(y)
cor(rx,ry)
sum((rx-7)*(ry-7))/sqrt(sum((rx-7)^2)*sum((ry-7)^2))

# example:  using R (hard way)
x = c(277,169,157,139,108,213,232,229,114,232,161,149,128)
y = c(256,118,137,144,146,221,184,188,97,231,114,187,230)
rx = rank(x)
ry = rank(y)
mx = mean(rx)
my = mean(ry)
sum((rx-mx)*(ry-my))/sqrt(sum((rx-mx)^2)*sum((ry-my)^2))
cor(rx,ry)

# example:  using R (easy way)
x = c(277,169,157,139,108,213,232,229,114,232,161,149,128)
y = c(256,118,137,144,146,221,184,188,97,231,114,187,230)
cor(x,y,method="spearman")
cor.test(x,y,method="spearman")
