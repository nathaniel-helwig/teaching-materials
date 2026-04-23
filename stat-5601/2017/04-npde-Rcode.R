##########   Nonparametric Dispersion and Equality Tests
##########   Nathaniel E. Helwig (helwig@umn.edu)
##########   Updated: 04-Jan-2017


#####   DEFINE PATHS AND PACKAGES   #####

# load NSM3 package (for various functions)
if(!require(NSM3)){
  install.packages("NSM3")
  library("NSM3")
}


#####   DISPERSION TEST (ANSARI-BRADLEY)   #####

# example 1:  using R (hard way)
set.seed(1)
x = round(rnorm(11),2)
y = round(rnorm(10,sd=2),2)
m = length(x)
n = length(y)
N = m + n
z = sort(c(x,y),index=TRUE)
rz = seq(1,(N-1)/2)
rz = c(rz,(N+1)/2,rev(rz))
r = rz[sort(z$ix,index=TRUE)$ix]
sum(r[1:11])
sum(r[12:21])

# example 1:  using R (easy way)
set.seed(1)
x = round(rnorm(11),2)
y = round(rnorm(10,sd=2),2)
ansari.test(x,y)
ansari.test(x,y,alternative="less")
ansari.test(x,y,alternative="less",conf.int=0.95)  ## TYPO!!


#####   DISPERSION/LOCATION TEST (LEPAGE)   #####

# example 2:  get Wilcoxon rank sum test statistic
set.seed(1)
x = round(rnorm(11),2)
y = round(rnorm(10,sd=2),2)
m = length(x)
n = length(y)
N = m + n
rk = rank(c(x,y))
sum(rk[1:11])
sum(rk[12:21])

# example 2:  using R (hard way)
set.seed(1)
x = round(rnorm(11),2)
y = round(rnorm(10,sd=2),2)
m = length(x)
n = length(y)
N = m + n
z = sort(c(x,y),index=TRUE)
rz = seq(1,(N-1)/2)
rz = c(rz,(N+1)/2,rev(rz))
r = rz[sort(z$ix,index=TRUE)$ix]
C = sum(r[12:21])
rk = rank(c(x,y))
W = sum(rk[12:21])
Wstar = (W-n*(N+1)/2)/sqrt(m*n*(N+1)/12)
Cstar = (C-n*((N+1)^2)/(4*N))/sqrt(m*n*(N+1)*(3+N^2)/(48*(N^2)))
D = Wstar^2 + Cstar^2
D
1 - pchisq(D,2)

# example 2:  using R (easy way)
set.seed(1)
x = round(rnorm(11),2)
y = round(rnorm(10,sd=2),2)
pLepage(x,y)


#####   EQUALITY TEST (KOLMOGOROV-SMIRNOV)   #####

# empirical CDFs of X with m=3 and n=2
F1mat = matrix(c(0,0,1/3,2/3,1,
                 0,1/3,1/3,2/3,1,
                 0,1/3,2/3,2/3,1,
                 0,1/3,2/3,1,1,
                 1/3,1/3,1/3,2/3,1,
                 1/3,1/3,2/3,2/3,1,
                 1/3,1/3,2/3,1,1,
                 1/3,2/3,2/3,2/3,1,
                 1/3,2/3,2/3,1,1,
                 1/3,2/3,1,1,1),10,5,byrow=TRUE)

# empirical CDFs of Y with m=3 and n=2
F2mat = matrix(c(1/2,1,1,1,1,
                 1/2,1/2,1,1,1,
                 1/2,1/2,1/2,1,1,
                 1/2,1/2,1/2,1/2,1,
                 0,1/2,1,1,1,
                 0,1/2,1/2,1,1,
                 0,1/2,1/2,1/2,1,
                 0,0,1/2,1,1,
                 0,0,1/2,1/2,1,
                 0,0,0,1/2,1),10,5,byrow=TRUE)

# get omega-hat for each possibility
apply(abs(F1mat-F2mat),1,max)

# example 3:  using R (hard way)
set.seed(1)
x = round(rnorm(11),2)
y = round(rnorm(10,sd=2),2)
z = sort(c(x,y),index=TRUE)
zlab = c(rep("x",11),rep("y",10))
j = ifelse(zlab[z$ix]=="x",1L,2L)
F1vec = F2vec = 0
for(k in 2:22){
  if(j[k-1]==1L){
    F1vec = c(F1vec,F1vec[k-1]+1)
    F2vec = c(F2vec,F2vec[k-1]+0)
  } else{
    F1vec = c(F1vec,F1vec[k-1]+0)
    F2vec = c(F2vec,F2vec[k-1]+1)
  }
}
F1vec = F1vec[2:22]/11
F2vec = F2vec[2:22]/10
omega = abs(F1vec-F2vec)
max(omega)
11*10*max(omega)
data.frame(z=z$x, j=j, F1=F1vec, F2=F2vec, omega=omega)

# example 3:  using R (easy way)
set.seed(1)
x = round(rnorm(11),2)
y = round(rnorm(10,sd=2),2)
ks.test(x,y)
ks.test(x,y,alternative="less")
max(F2vec-F1vec)
ks.test(x,y,alternative="greater")
max(F1vec-F2vec)

# more data (and ties)
set.seed(1)
x = round(rnorm(100),2)
y = round(rnorm(100,sd=2),2)
ks.test(x,y)
