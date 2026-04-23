##########   Permutation Tests
##########   Nathaniel E. Helwig (helwig@umn.edu)
##########   Updated: 04-Jan-2017


#####   INTRODUCTION TO PERMUTATIONS   #####

# permutations function
# from http://stackoverflow.com/questions/11095992/generating-all-distinct-permutations-of-a-list-in-r
permutations <- function(n){
  if(n==1){
    return(matrix(1))
  } else {
    sp <- permutations(n-1)
    p <- nrow(sp)
    A <- matrix(nrow=n*p,ncol=n)
    for(i in 1:n){
      A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
    }
    return(A)
  }
}

# permutation examples
permutations(2)
permutations(3)
permutations(4)
permutations(5)
permutations(6)

# random permutations
set.seed(1)
n = 5
x = seq(0,20,length=n)
x
x[sample.int(n)]
x[sample.int(n)]


#####   ONE-SAMPLE PERMUTATION TESTS   #####

# one-sample permutation test function (Monte Carlo)
perm1samp <- function(x,myfun=mean,mu=0,nsamp=10000,
                      alternative=c("two.sided","less","greater")){  
  x = x - mu
  n = length(x)
  theta.hat = myfun(x)
  gmat = replicate(nsamp,sample(x=c(1,-1),size=n,replace=TRUE))
  theta.mc = apply(gmat*abs(x),2,myfun)
  if(alternative[1]=="less"){
    aslperm = sum(theta.mc <= theta.hat) / nsamp
  } else if(alternative[1]=="greater"){
    aslperm = sum(theta.mc >= theta.hat) / nsamp
  } else{
    aslperm = sum(abs(theta.mc) >= abs(theta.hat)) / nsamp
  }
  list(theta.hat=theta.hat,theta.mc=theta.mc,asl=aslperm)
}

# one-sample permutation test: statistic 1 (mean)
set.seed(1)
n = 50
x = rnorm(n,mean=1)
mean(x)
se = (sd(x)/sqrt(n))
cv = qt(.975,df=n-1)
c(mean(x)-cv*se, mean(x)+cv*se)
mseq = seq(0.5,1.5,by=0.1)
pvals = rep(0,length(mseq))
for(k in 1:length(mseq)){
  pvals[k] = perm1samp(x,mu=mseq[k])$asl
}
pval.stat1 = pvals
dev.new(width=4,height=4,noRStudioGD=TRUE)
plot(mseq,pvals,type="b",xlab=expression("Median under "*H[0]),
     ylab="permutation ASL",main="Statistic 1: Sample Mean",ylim=c(0,1))
lines(rep(mean(x),2),c(0,1),lty=1)
lines(rep(median(x),2),c(0,1),lty=2)
lines(c(0.5,1.5),rep(0.05,2),lty=3)
legend("topleft",c("mean","median","p=0.05"),lty=1:3,bty="n")

# one-sample permutation test: statistic 2 (positive signed rank)
set.seed(1)
n = 50
x = rnorm(n,mean=1)
median(x)
myfun <- function(x) {
  n = length(x)
  rx = rank(abs(x))
  sum(rx[x>0]) - n*(n+1)/4
}
mseq = seq(0.5,1.5,by=0.1)
pvals = rep(0,length(mseq))
for(k in 1:length(mseq)){
  pvals[k] = perm1samp(x,myfun,mu=mseq[k])$asl
}
pval.stat2 = pvals
dev.new(width=4,height=4,noRStudioGD=TRUE)
plot(mseq,pvals,type="b",xlab=expression("Median under "*H[0]),
     ylab="permutation ASL",main="Statistic 2: Signed Rank",ylim=c(0,1))
lines(rep(mean(x),2),c(0,1),lty=1)
lines(rep(median(x),2),c(0,1),lty=2)
lines(c(0.5,1.5),rep(0.05,2),lty=3)
legend("topleft",c("mean","median","p=0.05"),lty=1:3,bty="n")

# one-sample permutation test: statistic 3 (sign statistic)
set.seed(1)
n = 50
x = rnorm(n,mean=1)
myfun <- function(x) {
  n = length(x)
  sum(x>0) - n/2
}
mseq = seq(0.5,1.5,by=0.1)
pvals = rep(0,length(mseq))
for(k in 1:length(mseq)){
  pvals[k] = perm1samp(x,myfun,mu=mseq[k])$asl
}
pval.stat3 = pvals
dev.new(width=4,height=4,noRStudioGD=TRUE)
plot(mseq,pvals,type="b",xlab=expression("Median under "*H[0]),
     ylab="permutation ASL",main="Statistic 3: Sign",ylim=c(0,1))
lines(rep(mean(x),2),c(0,1),lty=1)
lines(rep(median(x),2),c(0,1),lty=2)
lines(c(0.5,1.5),rep(0.05,2),lty=3)
legend("topleft",c("mean","median","p=0.05"),lty=1:3,bty="n")


#####   TWO-SAMPLE PERMUTATION TESTS   #####

# two-sample permutation test function (Monte Carlo)
meandif = function(x,y){ mean(x) - mean(y) }
perm2samp <- function(x,y,myfun=meandif,nsamp=10000,
                      alternative=c("two.sided","less","greater")){
  theta.hat = myfun(x,y)
  m = length(x)
  n = length(y)
  N = m + n
  z = c(x,y)
  gmat = replicate(nsamp,sample.int(N,m))
  theta.mc = apply(gmat,2,function(g,z){myfun(z[g],z[-g])},z=z)
  if(alternative[1]=="less"){
    aslperm = sum(theta.mc <= theta.hat) / nsamp
  } else if(alternative[1]=="greater"){
    aslperm = sum(theta.mc >= theta.hat) / nsamp
  } else{
    aslperm = sum(abs(theta.mc) >= abs(theta.hat)) / nsamp
  }
  list(theta.hat=theta.hat,theta.mc=theta.mc,asl=aslperm)
}

# two-sample permutation test: statistic 1 (mean difference)
dev.new(width=6,height=6,noRStudioGD=TRUE)
set.seed(1)
x = rnorm(15)
y = rnorm(20,mean=1)
choose(35,15)
myfun = function(x,y) mean(x) - mean(y)
myfun(x,y)
mean(x) - mean(y)
ptest = perm2samp(x,y,myfun)
ptest$theta.hat
ptest$asl
hist(ptest$theta.mc)
lines(rep(ptest$theta.hat,2),c(0,2000),col="red",lty=2)
t.test(x,y)

# two-sample permutation test: statistic 2 (rank sum)
dev.new(width=6,height=6,noRStudioGD=TRUE)
set.seed(1)
x = rnorm(15)
y = rnorm(20,mean=1)
choose(35,15)
myfun = function(x,y){
  m = length(x)
  n = length(y)
  rx = rank(c(x,y))
  sum(rx[seq(along=x)]) - m*(m+n+1)/2
}
myfun(x,y)
ptest = perm2samp(x,y,myfun)
ptest$theta.hat
ptest$asl
hist(ptest$theta.mc)
lines(rep(ptest$theta.hat,2),c(0,2000),col="red",lty=2)

# two-sample permutation test: statistic 3 (variance)
dev.new(width=6,height=6,noRStudioGD=TRUE)
set.seed(1)
x = rnorm(15)
y = rnorm(20,sd=3)
choose(35,15)
myfun = function(x,y) log(var(x)/var(y))
myfun(x,y)
log(var(x)/var(y))
ptest = perm2samp(x,y,myfun)
ptest$theta.hat
ptest$asl
hist(ptest$theta.mc)
lines(rep(ptest$theta.hat,2),c(0,2000),col="red",lty=2)


#####   CORRELATION PERMUTATION TESTS   #####

# correlation permutation test function (Monte Carlo)
permcor <- function(x,y,method="pearson",nsamp=10000,
                    alternative=c("two.sided","less","greater")){
  n = length(x)
  if(n!=length(y)) stop("lengths of x and y must match")
  theta.hat = cor(x,y,method=method)
  gmat = replicate(nsamp,sample.int(n))
  theta.mc = apply(gmat,2,function(g){cor(x,y[g],method=method)})
  if(alternative[1]=="less"){
    aslperm = sum(theta.mc <= theta.hat) / nsamp
  } else if(alternative[1]=="greater"){
    aslperm = sum(theta.mc >= theta.hat) / nsamp
  } else{
    aslperm = sum(abs(theta.mc) >= abs(theta.hat)) / nsamp
  }
  list(theta.hat=theta.hat,theta.mc=theta.mc,asl=aslperm)
}

# correlation permutation test: statistic 1 (Pearson)
dev.new(width=6,height=6,noRStudioGD=TRUE)
set.seed(1)
n = 50
x = rnorm(n)
y = rnorm(n)
rho = -0.2
Amat = matrix(c(1,rho,rho,1),2,2)
Aeig = eigen(Amat,symmetric=TRUE)
evec = Aeig$vec
evalsqrt = diag(Aeig$val^0.5)
Asqrt = evec %*% evalsqrt %*% t(evec)
z = cbind(x,y)%*%Asqrt
x = z[,1]
y = z[,2]
ptest = permcor(x,y)
ptest$asl
hist(ptest$theta.mc)
lines(rep(ptest$theta.hat,2),c(0,2000),col="red",lty=2)

# correlation permutation test: statistic 2 (Spearman)
dev.new(width=6,height=6,noRStudioGD=TRUE)
set.seed(1)
n = 50
x = rnorm(n)
y = rnorm(n)
rho = -0.2
Amat = matrix(c(1,rho,rho,1),2,2)
Aeig = eigen(Amat,symmetric=TRUE)
evec = Aeig$vec
evalsqrt = diag(Aeig$val^0.5)
Asqrt = evec %*% evalsqrt %*% t(evec)
z = cbind(x,y)%*%Asqrt
x = z[,1]
y = z[,2]
ptest = permcor(x,y,method="spearman")
ptest$asl
hist(ptest$theta.mc)
lines(rep(ptest$theta.hat,2),c(0,2000),col="red",lty=2)

# correlation permutation test: statistic 3 (Kendall)
dev.new(width=6,height=6,noRStudioGD=TRUE)
set.seed(1)
n = 50
x = rnorm(n)
y = rnorm(n)
rho = -0.2
Amat = matrix(c(1,rho,rho,1),2,2)
Aeig = eigen(Amat,symmetric=TRUE)
evec = Aeig$vec
evalsqrt = diag(Aeig$val^0.5)
Asqrt = evec %*% evalsqrt %*% t(evec)
z = cbind(x,y)%*%Asqrt
x = z[,1]
y = z[,2]
ptest = permcor(x,y,method="kendall")
ptest$asl
hist(ptest$theta.mc)
lines(rep(ptest$theta.hat,2),c(0,2000),col="red",lty=2)

