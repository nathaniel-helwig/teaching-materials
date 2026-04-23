##########   Inferences about Multivariate Means
##########   Nathaniel E. Helwig (helwig@umn.edu)
##########   Updated: 16-Jan-2017


#####   DEFINE PATHS AND PACKAGES   #####

# load car package (for ellipse functions)
if(!require(car)){
  install.packages("car")
  library("car")
}

# load lsmeans package (for lsmeans and test functions)
if(!require(lsmeans)){
  install.packages("lsmeans")
  library("lsmeans")
}

# load biotools package (for boxM function)
if(!require(biotools)){
  install.packages("biotools")
  library("biotools")
}


#####   SINGLE MEAN VECTOR   #####

# define Hotelling T^2 function
T.test <- function(X, mu=0){
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  df2 <- n - p
  if(df2 < 1L) stop("Need nrow(X) > ncol(X).")
  if(length(mu) != p) mu <- rep(mu[1], p)
  xbar <- colMeans(X)
  S <- cov(X)
  T2 <- n * t(xbar - mu) %*% solve(S) %*% (xbar - mu)
  Fstat <- T2 / (p * (n-1) / df2)
  pval <- 1 - pf(Fstat, df1=p, df2=df2)
  data.frame(T2=as.numeric(T2), Fstat=as.numeric(Fstat), 
             df1=p, df2=df2, p.value=as.numeric(pval), row.names="")
}

# get data matrix
data(mtcars)
X <- mtcars[,c("mpg","disp","hp","wt")]
xbar <- colMeans(X)

# try Hotelling T^2 function
xbar
T.test(X)
T.test(X, mu=c(20,200,150,3))
T.test(X, mu=xbar)

# using lm function
y <- as.matrix(X)
anova(lm(y ~ 1))
y <- as.matrix(X) - matrix(c(20,200,150,3),nrow(X),ncol(X),byrow=T)
anova(lm(y ~ 1))
y <- as.matrix(X) - matrix(xbar,nrow(X),ncol(X),byrow=T)
anova(lm(y ~ 1))

# plot 2D confidence region (ellipse)
X <- mtcars[,c("mpg","disp","hp","wt")]
n <- nrow(X)
p <- ncol(X)
xbar <- colMeans(X)
S <- cov(X)
xnames <- names(xbar)
dev.new(width=8,height=4,noRStudioGD=TRUE)
par(mfrow=c(1,2))
levels <- c(0.99,0.95,0.90)
for(ii in 1:2){
  id <- c(ii,3)
  for(jj in 1:3){
    tconst <- sqrt((p/n)*((n-1)/(n-p)) * qf(levels[jj],p,n-p))
    if(jj>1){
      points(ellipse(center=xbar[id], shape=S[id,id], radius=tconst, draw=F), col=jj, pch=jj)
    } else {
      plot(ellipse(center=xbar[id], shape=S[id,id], radius=tconst, draw=F), 
           col=jj, pch=jj, xlab=xnames[id[1]], ylab=xnames[id[2]])
    }
    points(xbar[id[1]],xbar[id[2]],pch="*")
  }
  legend(ifelse(ii==1,"topright","topleft"),legend=levels,col=1:3,pch=1:3)
}

# T.ci function
T.ci <- function(mu, Sigma, n, avec=rep(1,length(mu)), level=0.95){
  p <- length(mu)
  if(nrow(Sigma)!=p) stop("Need length(mu) == nrow(Sigma).")
  if(ncol(Sigma)!=p) stop("Need length(mu) == ncol(Sigma).")
  if(length(avec)!=p) stop("Need length(mu) == length(avec).")
  if(level <=0 | level >= 1) stop("Need 0 < level < 1.")
  cval <- qf(level, p, n-p) * p * (n-1) / (n-p)
  zhat <- crossprod(avec, mu)
  zvar <- crossprod(avec, Sigma %*% avec) / n
  const <- sqrt(cval * zvar)
  c(lower = zhat - const, upper = zhat + const)
}

# T.ci example
X <- mtcars[,c("mpg","disp","hp","wt")]
n <- nrow(X)
p <- ncol(X)
xbar <- colMeans(X)
S <- cov(X)
xbar
T.ci(mu=xbar, Sigma=S, n=n, avec=c(1,0,0,0))
T.ci(mu=xbar, Sigma=S, n=n, avec=c(0,1,0,0))
T.ci(mu=xbar, Sigma=S, n=n, avec=c(0,0,1,0))
T.ci(mu=xbar, Sigma=S, n=n, avec=c(0,0,0,1))

# compare simultaneous vs one-at-a-time
TCI <- tCI <- NULL
for(k in 1:4){
  avec <- rep(0, 4)
  avec[k] <- 1
  TCI <- c(TCI, T.ci(xbar, S, n, avec))
  tCI <- c(tCI,
           xbar[k] - sqrt(S[k,k]/n) * qt(0.975, df=n-1),
           xbar[k] + sqrt(S[k,k]/n) * qt(0.975, df=n-1))
}
cnames <- paste(rep(names(X),each=2),sep=".",names(TCI))
rtab <- rbind(TCI, tCI)
colnames(rtab) <- cnames
round(rtab, 2)

# Bonferroni method
TCI <- tCI <- bon <- NULL
alpha <- 1 - 0.05/(2*4)
for(k in 1:4){
  avec <- rep(0, 4)
  avec[k] <- 1
  TCI <- c(TCI, T.ci(xbar, S, n, avec))
  tCI <- c(tCI,
           xbar[k] - sqrt(S[k,k]/n) * qt(0.975, df=n-1),
           xbar[k] + sqrt(S[k,k]/n) * qt(0.975, df=n-1))
  bon <- c(bon, 
           xbar[k] - sqrt(S[k,k]/n) * qt(alpha, df=n-1),
           xbar[k] + sqrt(S[k,k]/n) * qt(alpha, df=n-1))
}
cnames <- paste(rep(names(X),each=2),sep=".",names(TCI))
rtab <- rbind(TCI, tCI, bon)
colnames(rtab) <- cnames
round(rtab, 2)

# update T.test function with 'asymp' option
T.test <- function(X, mu=0, asymp=FALSE){
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  df2 <- n - p
  if(df2 < 1L) stop("Need nrow(X) > ncol(X).")
  if(length(mu) != p) mu <- rep(mu[1], p)
  xbar <- colMeans(X)
  S <- cov(X)
  T2 <- n * t(xbar - mu) %*% solve(S) %*% (xbar - mu)
  Fstat <- T2 / (p * (n-1) / df2)
  if(asymp){
    pval <- 1 - pchisq(T2, df=p)
  } else {
    pval <- 1 - pf(Fstat, df1=p, df2=df2)
  }
  data.frame(T2=as.numeric(T2), Fstat=as.numeric(Fstat), 
             df1=p, df2=df2, p.value=as.numeric(pval), 
             asymp=asymp, row.names="")
}

# compare finite sample and large sample p-values (n=10)
set.seed(1)
XX <- matrix(rnorm(10*4), 10, 4)
T.test(XX)
T.test(XX, asymp=TRUE)

# compare finite sample and large sample p-values (n=50)
set.seed(1)
XX <- matrix(rnorm(50*4), 50, 4)
T.test(XX)
T.test(XX, asymp=TRUE)

# compare finite sample and large sample p-values (n=100)
set.seed(1)
XX <- matrix(rnorm(100*4), 100, 4)
T.test(XX)
T.test(XX, asymp=TRUE)

# large sample (chi-square) method
TCI <- tCI <- bon <- chi <- NULL
alpha <- 1 - 0.05/(2*4)
for(k in 1:4){
  avec <- rep(0, 4)
  avec[k] <- 1
  TCI <- c(TCI, T.ci(xbar, S, n, avec))
  tCI <- c(tCI,
           xbar[k] - sqrt(S[k,k]/n) * qt(0.975, df=n-1),
           xbar[k] + sqrt(S[k,k]/n) * qt(0.975, df=n-1))
  bon <- c(bon, 
           xbar[k] - sqrt(S[k,k]/n) * qt(alpha, df=n-1),
           xbar[k] + sqrt(S[k,k]/n) * qt(alpha, df=n-1))
  chi <- c(chi, 
           xbar[k] - sqrt(S[k,k]/n) * sqrt(qchisq(0.95, df=p)),
           xbar[k] + sqrt(S[k,k]/n) * sqrt(qchisq(0.95, df=p)))
}
cnames <- paste(rep(names(X),each=2),sep=".",names(TCI))
rtab <- rbind(TCI, tCI, bon, chi)
colnames(rtab) <- cnames
round(rtab, 2)

# plot 2D prediction region (ellipse)
X <- mtcars[,c("mpg","disp","hp","wt")]
n <- nrow(X)
p <- ncol(X)
xbar <- colMeans(X)
S <- cov(X)
dev.new(width=8,height=4,noRStudioGD=TRUE)
par(mfrow=c(1,2))
levels <- c(0.99,0.95,0.90)
for(ii in 1:2){
  id <- c(ii,3)
  for(jj in 1:3){
    mconst <- sqrt((p/n)*((n-1)/(n-p)) * qf(levels[jj],p,n-p))
    pconst <- sqrt((p/n)*((n^2-1)/(n-p)) * qf(levels[jj],p,n-p))
    if(jj>1){
      points(ellipse(center=xbar[id], shape=S[id,id], radius=pconst, draw=F), col=jj, pch=jj)
      lines(ellipse(center=xbar[id], shape=S[id,id], radius=mconst, draw=F), col=jj, lty=jj)
    } else {
      plot(ellipse(center=xbar[id], shape=S[id,id], radius=pconst, draw=F), 
           col=jj, pch=jj, xlab=xnames[id[1]], ylab=xnames[id[2]])
      lines(ellipse(center=xbar[id], shape=S[id,id], radius=mconst, draw=F), col=jj, lty=jj)
    }
    points(xbar[id[1]],xbar[id[2]],pch="*")
  }
  legend(ifelse(ii==1,"topright","topleft"),legend=levels,col=1:3,pch=1:3)
}


#####   MULTIPLE MEAN VECTORS   #####

# define RM.test function
RM.test <- function(X, mu=0, C=NULL){
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  df2 <- n - p + 1
  if(df2 < 1L) stop("Need nrow(X) > ncol(X).")
  if(length(mu) != p) mu <- rep(mu[1], p)
  xbar <- colMeans(X)
  S <- cov(X)
  if(is.null(C)){
    C <- matrix(0, p-1, p)
    for(k in 1:(p-1)) C[k, 1:2 + 1*(k-1)] <- c(1, -1)
  } else {
    if(nrow(C) != (p-1)) stop("Need [ncol(X)-1] == nrow(C).")
    if(ncol(C) != p) stop("Need ncol(X) == ncol(C).")
    if(any(rowSums(C)>0L)) stop("Need rowSums(C) == rep(0, nrow(C)).")
  }
  T2 <- n * t(C %*% (xbar - mu)) %*% solve(C %*% S %*% t(C)) %*% (C %*% (xbar - mu))
  Fstat <- T2 / ((p-1) * (n-1) / df2)
  pval <- 1 - pf(Fstat, df1=p-1, df2=df2)
  data.frame(T2=as.numeric(T2), Fstat=as.numeric(Fstat), 
             df1=p-1, df2=df2, p.value=as.numeric(pval), row.names="")
}

# RM.test example w/ H0 true (10 data points)
set.seed(1)
XX <- matrix(rnorm(10*4), 10, 4)
RM.test(XX)

# RM.test example w/ H0 true (100 data points)
set.seed(1)
XX <- matrix(rnorm(100*4), 100, 4)
RM.test(XX)

# RM.test example w/ H0 true (500 data points)
set.seed(1)
XX <- matrix(rnorm(500*4), 500, 4)
RM.test(XX)

# RM.test example w/ H0 false (10 data points)
set.seed(1)
XX <- matrix(rnorm(10*4), 10, 4)
XX <- XX + matrix(c(0,0,0,0.25), 10, 4, byrow=TRUE)
RM.test(XX)

# RM.test example w/ H0 false (100 data points)
set.seed(1)
XX <- matrix(rnorm(100*4), 100, 4)
XX <- XX + matrix(c(0,0,0,0.25), 100, 4, byrow=TRUE)
RM.test(XX)

# RM.test example w/ H0 false (500 data points)
set.seed(1)
XX <- matrix(rnorm(500*4), 500, 4)
XX <- XX + matrix(c(0,0,0,0.25), 500, 4, byrow=TRUE)
RM.test(XX)

# define two-sample T.test function
T.test <- function(X, Y=NULL, mu=0, paired=FALSE, asymp=FALSE){
  if(is.null(Y)){
    # one-sample T^2 test
    X <- as.matrix(X)
    nx <- nrow(X)
    p <- ncol(X)
    df2 <- nx - p
    if(df2 < 1L) stop("Need nrow(X) > ncol(X).")
    if(length(mu) != p) mu <- rep(mu[1], p)
    xbar <- colMeans(X)
    S <- cov(X)
    T2 <- nx * t(xbar - mu) %*% solve(S) %*% (xbar - mu)
    Fstat <- T2 / (p * (nx-1) / df2)
    if(asymp){
      pval <- 1 - pchisq(T2, df=p)
    } else {
      pval <- 1 - pf(Fstat, df1=p, df2=df2)
    }
    return(data.frame(T2=as.numeric(T2), Fstat=as.numeric(Fstat), 
                      df1=p, df2=df2, p.value=as.numeric(pval), 
                      type="one-sample", asymp=asymp, row.names=""))
  } else {
    if(paired){
      # dependent two-sample T^2 test
      X <- as.matrix(X)
      Y <- as.matrix(Y)
      if(!identical(dim(X),dim(Y))) stop("Need dim(X) == dim(Y).")
      xx <- T.test(X-Y, mu=mu, asymp=asymp)
      xx$type <- "dep-sample"
      return(xx)
    } else {
      # independent two-sample T^2 test
      X <- as.matrix(X)
      Y <- as.matrix(Y)
      nx <- nrow(X)
      ny <- nrow(Y)
      p <- ncol(X)
      df2 <- nx + ny - p - 1
      if(p != ncol(Y)) stop("Need ncol(X) == ncol(Y).")
      if(min(nx,ny) <= p) stop("Need min(nrow(X),nrow(Y)) > ncol(X).")
      Sp <- ((nx-1)*cov(X) + (ny-1)*cov(Y)) / (nx + ny - 2)
      dbar <- colMeans(X) - colMeans(Y)
      T2 <- (1/((1/nx) + (1/ny))) * t(dbar - mu) %*% solve(Sp) %*% (dbar - mu)
      Fstat <- T2 / ((nx + ny - 2) * p / df2)
      if(asymp){
        pval <- 1 - pchisq(T2, df=p)
      } else {
        pval <- 1 - pf(Fstat, df1=p, df2=df2)
      }
      return(data.frame(T2=as.numeric(T2), Fstat=as.numeric(Fstat), 
                        df1=p, df2=df2, p.value=as.numeric(pval), 
                        type="ind-sample", asymp=asymp, row.names=""))
    } # end if(paired)
  } # end if(is.null(Y))
} # end T.test function

# independent samples T^2 test
X4 <- subset(mtcars, cyl==4)[,c("mpg","disp","hp","wt")]
X6 <- subset(mtcars, cyl==6)[,c("mpg","disp","hp","wt")]
X8 <- subset(mtcars, cyl==8)[,c("mpg","disp","hp","wt")]
T.test(X4, X6)
T.test(X4, X8)
T.test(X6, X8)

# compare finite sample and large sample p-values (n=10)
set.seed(1)
n <- 10
XX <- matrix(rnorm(n*4), n, 4)
YY <- matrix(rnorm(n*4), n, 4)
T.test(XX,YY)
T.test(XX,YY,asymp=T)

# compare finite sample and large sample p-values (n=50)
set.seed(1)
n <- 50
XX <- matrix(rnorm(n*4), n, 4)
YY <- matrix(rnorm(n*4), n, 4)
T.test(XX,YY)
T.test(XX,YY,asymp=T)

# compare finite sample and large sample p-values (n=100)
set.seed(1)
n <- 100
XX <- matrix(rnorm(n*4), n, 4)
YY <- matrix(rnorm(n*4), n, 4)
T.test(XX,YY)
T.test(XX,YY,asymp=T)

# redefine T.ci function
T.ci <- function(mu, Sigma, n, avec=rep(1,length(mu)), level=0.95){
  p <- length(mu)
  if(nrow(Sigma)!=p) stop("Need length(mu) == nrow(Sigma).")
  if(ncol(Sigma)!=p) stop("Need length(mu) == ncol(Sigma).")
  if(length(avec)!=p) stop("Need length(mu) == length(avec).")
  if(level <=0 | level >= 1) stop("Need 0 < level < 1.")
  zhat <- crossprod(avec, mu)
  if(length(n)==1L){
    cval <- qf(level, p, n-p) * p * (n-1) / (n-p)
    zvar <- crossprod(avec, Sigma %*% avec) / n
  } else {
    df2 <- n[1] + n[2] - p - 1
    cval <- qf(level, p, df2) * p * (n[1]+n[2]-2) / df2
    zvar <- crossprod(avec, Sigma %*% avec) * ( (1/n[1]) + (1/n[2]) )
  }
  const <- sqrt(cval * zvar)
  c(lower = zhat - const, upper = zhat + const)
}

# try T2ci with independent samples
X4 <- subset(mtcars, cyl==4)[,c("mpg","disp","hp","wt")]
X6 <- subset(mtcars, cyl==6)[,c("mpg","disp","hp","wt")]
n4 <- nrow(X4)
n6 <- nrow(X6)
dbar <- colMeans(X4) - colMeans(X6)
Sp <- ((n4-1)*cov(X4) + (n6-1)*cov(X6)) / (n4 + n6 - 2)
dbar
T.ci(dbar, Sp, c(n4,n6), c(1,0,0,0))
T.ci(dbar, Sp, c(n4,n6), c(0,1,0,0))
T.ci(dbar, Sp, c(n4,n6), c(0,0,1,0))
T.ci(dbar, Sp, c(n4,n6), c(0,0,0,1))

# redefine two-sample T.test function (with var.equal option)
T.test <- function(X, Y=NULL, mu=0, paired=FALSE, asymp=FALSE, var.equal=TRUE){
  if(is.null(Y)){
    # one-sample T^2 test
    X <- as.matrix(X)
    nx <- nrow(X)
    p <- ncol(X)
    df2 <- nx - p
    if(df2 < 1L) stop("Need nrow(X) > ncol(X).")
    if(length(mu) != p) mu <- rep(mu[1], p)
    xbar <- colMeans(X)
    S <- cov(X)
    T2 <- nx * t(xbar - mu) %*% solve(S) %*% (xbar - mu)
    Fstat <- T2 / (p * (nx-1) / df2)
    if(asymp){
      pval <- 1 - pchisq(T2, df=p)
    } else {
      pval <- 1 - pf(Fstat, df1=p, df2=df2)
    }
    return(data.frame(T2=as.numeric(T2), Fstat=as.numeric(Fstat), 
                      df1=p, df2=df2, p.value=as.numeric(pval), 
                      type="one-sample", asymp=asymp, row.names=""))
  } else {
    if(paired){
      # dependent two-sample T^2 test
      X <- as.matrix(X)
      Y <- as.matrix(Y)
      if(!identical(dim(X),dim(Y))) stop("Need dim(X) == dim(Y).")
      xx <- T.test(X-Y, mu=mu, asymp=asymp)
      xx$type <- "dep-sample"
      return(xx)
    } else {
      # independent two-sample T^2 test
      X <- as.matrix(X)
      Y <- as.matrix(Y)
      nx <- nrow(X)
      ny <- nrow(Y)
      p <- ncol(X)
      if(p != ncol(Y)) stop("Need ncol(X) == ncol(Y).")
      if(min(nx,ny) <= p) stop("Need min(nrow(X),nrow(Y)) > ncol(X).")
      dbar <- colMeans(X) - colMeans(Y)
      if(var.equal){
        df2 <- nx + ny - p - 1
        Sp <- ((nx-1)*cov(X) + (ny-1)*cov(Y)) / (nx + ny - 2)
        T2 <- (1/((1/nx) + (1/ny))) * t(dbar - mu) %*% solve(Sp) %*% (dbar - mu)
        Fstat <- T2 / ((nx + ny - 2) * p / df2)
      } else {
        Sx <- cov(X)
        Sy <- cov(Y)
        Sp <- (Sx/nx) + (Sy/ny)
        T2 <- t(dbar - mu) %*% solve(Sp) %*% (dbar - mu)
        SpInv <- solve(Sp)
        SxSpInv <- (1/nx) * Sx %*% SpInv
        SySpInv <- (1/ny) * Sy %*% SpInv
        nudx <- (sum(diag(SxSpInv %*% SxSpInv)) + (sum(diag(SxSpInv)))^2) / nx
        nudy <- (sum(diag(SySpInv %*% SySpInv)) + (sum(diag(SySpInv)))^2) / ny
        nu <- (p + p^2) / (nudx + nudy)
        df2 <- nu - p + 1
        Fstat <- T2 / (nu * p / df2)
      }
      if(asymp){
        pval <- 1 - pchisq(T2, df=p)
      } else {
        pval <- 1 - pf(Fstat, df1=p, df2=df2)
      }
      return(data.frame(T2=as.numeric(T2), Fstat=as.numeric(Fstat), 
                        df1=p, df2=df2, p.value=as.numeric(pval), 
                        type="ind-sample", asymp=asymp, var.equal=var.equal, row.names=""))
    } # end if(paired)
  } # end if(is.null(Y))
} # end T.test function

# independent samples T^2 test w/ var.equal
X4 <- subset(mtcars, cyl==4)[,c("mpg","disp","hp","wt")]
X6 <- subset(mtcars, cyl==6)[,c("mpg","disp","hp","wt")]
T.test(X4, X6)
T.test(X4, X6, var.equal=FALSE)

# independent samples T^2 test w/ var.equal
set.seed(1)
n <- 100
XX <- matrix(rnorm(n*4), n, 4)
YY <- matrix(rnorm(n*4), n, 4)
T.test(XX,YY)
T.test(XX,YY,var.equal=F)

# one-way manova
data(mtcars)
X <- as.matrix(mtcars[,c("mpg","disp","hp","wt")])
cylinder <- factor(mtcars$cyl)
mod <- lm(X ~ cylinder)
summary(manova(mod))
Manova(mod, test.statistic="Pillai")
Manova(mod, test.statistic="Wilks")
Manova(mod, test.statistic="Roy")
Manova(mod, test.statistic="Hotelling-Lawley")
  
# get least-squares means for each variable
p <- ncol(X)
lsm <- vector("list", p)
names(lsm) <- colnames(X)
for(j in 1:p){
  wts <- rep(0, p)
  wts[j] <- 1
  lsm[[j]] <- lsmeans(mod, "cylinder", weights=wts)
}
lsm[[1]]
lsm[[3]]

# get alpha level for Bonferroni correction
q <- p * 3 * (3-1) / 2
alpha <- 0.05 / (2*q)

# Bonferroni pairwise CIs for "mpg"
confint(contrast(lsm[[1]], "pairwise"), level=1-alpha, adj="none")

# Bonferroni pairwise CIs for "hp"
confint(contrast(lsm[[3]], "pairwise"), level=1-alpha, adj="none")

# Box's M test (for homogeneity of covariances)
boxM(X, cylinder)

# two-way manova with interaction
data(mtcars)
X <- as.matrix(mtcars[,c("mpg","disp","hp","wt")])
cylinder <- factor(mtcars$cyl)
transmission <- factor(mtcars$am)
mod <- lm(X ~ cylinder * transmission)
Manova(mod, test.statistic="Wilks")

# refit additive model
mod <- lm(X ~ cylinder + transmission)
Manova(mod, test.statistic="Wilks")

# get LS means for cylinder and transmission
lsmeans(mod, "cylinder", weights="show.levels")
lsmeans(mod, "transmission", weights="show.levels")
p <- ncol(X)
lsm.cyl <- lsm.trn <- vector("list", p)
names(lsm) <- colnames(X)
for(j in 1:p){
  wts <- rep(0, p*2)
  wts[1:2 + (j-1)*2] <- 1
  lsm.cyl[[j]] <- lsmeans(mod, "cylinder", weights=wts)
  wts <- rep(0, p*3)
  wts[1:3 + (j-1)*3] <- 1
  lsm.trn[[j]] <- lsmeans(mod, "transmission", weights=wts)
}

# print mpg LS mean for cylinder effect
lsm.cyl[[1]]

# print mpg LS mean for transmission effect
lsm.trn[[1]]

# get alpha level for Bonferroni correction
q <- p * (3 * (3-1) / 2 + 2 * (2-1) / 2)
alpha <- 0.05 / (2*q)

# Bonferroni pairwise CIs for "mpg" (cylinder effect)
confint(contrast(lsm.cyl[[1]], "pairwise"), level=1-alpha, adj="none")

# Bonferroni pairwise CIs for "mpg" (tranmission effect)
confint(contrast(lsm.trn[[1]], "pairwise"), level=1-alpha, adj="none")

# Bonferroni pairwise CIs for "hp"  (cylinder effect)
confint(contrast(lsm.cyl[[3]], "pairwise"), level=1-alpha, adj="none")

# Bonferroni pairwise CIs for "hp"  (transmission effect)
confint(contrast(lsm.trn[[3]], "pairwise"), level=1-alpha, adj="none")

