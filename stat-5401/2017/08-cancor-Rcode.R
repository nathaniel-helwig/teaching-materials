##########   Canonical Correlation Analysis
##########   Nathaniel E. Helwig (helwig@umn.edu)
##########   Updated: 16-Mar-2017


#####   DECATHLON EXAMPLE: GET DATA   #####

# read-in decathlon data
decathlon <- read.table("~/Downloads/OLYMPIC.dat",row.names=1)
colnames(decathlon) <- c("run100","long.jump","shot","high.jump","run400","hurdle",
                         "discus","pole.vault","javelin","run1500","score")

# resign running events (so higher score means better performance)
head(decathlon)
decathlon[,c(1,5,6,10)] <- (-1)*decathlon[,c(1,5,6,10)]
head(decathlon)

# check variable scales
apply(decathlon,2,mean)
apply(decathlon,2,sd)

# separate into running/jumping vs throwing/arm events
X <- as.matrix(decathlon[,c("shot","discus","javelin","pole.vault")])
Y <- as.matrix(decathlon[,c("run100","run400","run1500",
                            "hurdle","long.jump","high.jump")])
n <- nrow(X)
p <- ncol(X)
q <- ncol(Y)


#####   DECATHLON EXAMPLE: COVARIANCE CCA   #####

# canonical correlations of covariance (unstandardized data)
cca <- cancor(X, Y)

# cca (the normal way)
Sx <- cov(X)
Sy <- cov(Y)
Sxy <- cov(X,Y)
Sxeig <- eigen(Sx, symmetric=TRUE)
Sxisqrt <- Sxeig$vectors %*% diag(1/sqrt(Sxeig$values)) %*% t(Sxeig$vectors)
Syeig <- eigen(Sy, symmetric=TRUE)
Syisqrt <- Syeig$vectors %*% diag(1/sqrt(Syeig$values)) %*% t(Syeig$vectors)
Xmat <- Sxisqrt %*% Sxy %*% solve(Sy) %*% t(Sxy) %*% Sxisqrt
Ymat <- Syisqrt %*% t(Sxy) %*% solve(Sx) %*% Sxy %*% Syisqrt
Xeig <- eigen(Xmat, symmetric=TRUE)
Yeig <- eigen(Ymat, symmetric=TRUE)

# compare correlations (same)
cca$cor
rho <- sqrt(Xeig$values)
rho
sqrt(Yeig$values[1:p])

# compare linear combinations (different!)
Ahat <- Sxisqrt %*% Xeig$vectors
Bhat <- Syisqrt %*% Yeig$vectors
sum((cca$xcoef - Ahat)^2)
sum((cca$ycoef[,1:p] - Bhat[,1:p])^2)

# NOTE: you need to multiply R's xcoef and ycoef by sqrt(n-1)
#       to obtain the results we are expecting...

# compare linear combinations (same!)
Ahat <- Sxisqrt %*% Xeig$vectors
Bhat <- Syisqrt %*% Yeig$vectors
sum((cca$xcoef * sqrt(n-1) - Ahat)^2)
sum((cca$ycoef[,1:p] * sqrt(n-1) - Bhat[,1:p])^2)

# plot coefficients
dev.new(width=10, height=5, noRStudioGD=TRUE)
par(mfrow=c(1,2))
plot(Ahat[,1:2], xlab="A1 Coefficients", ylab="A2 Coefficients",
     type="n", main="X Coefficients", xlim=c(-2, 0.1), ylim=c(-1.1, 0.5))
text(Ahat[,1:2], labels=colnames(X))
abline(0,0,lty=3)
abline(v=0,lty=3)
plot(Bhat[,1:2], xlab="B1 Coefficients", ylab="B2 Coefficients",
     type="n", main="Y Coefficients", xlim=c(-2, 0.2), ylim=c(-2, 6))
text(Bhat[,1:2], labels=colnames(Y))
abline(0,0,lty=3)
abline(v=0,lty=3)

# define canonical variates
U <- X %*% Ahat
V <- Y %*% Bhat

# canonical variable covariances
round(cov(U),4)
round(cov(V),4)
round(cov(U,V),4)

# covariance of original and canonical variables (U and X)
Ainv <- solve(Ahat)
sum( ( cov(U, X) - crossprod(Ahat, Sx) )^2 )
sum( ( Sx - crossprod(Ainv) )^2 )
Sxhat <- matrix(0, p, p)
for(j in 1:p) Sxhat <- Sxhat + outer(Ainv[j,], Ainv[j,])
sum( (Sx - Sxhat)^2 )

# covariance of original and canonical variables (V and Y)
Binv <- solve(Bhat)
sum( ( cov(V, Y) - crossprod(Bhat, Sy) )^2 )
sum( ( Sy - crossprod(Binv) )^2 )
Syhat <- matrix(0, q, q)
for(j in 1:q) Syhat <- Syhat + outer(Binv[j,], Binv[j,])
sum( (Sy - Syhat)^2 )

# covariance of original and canonical variables (U and Y)
sum( (cov(U, Y) - crossprod(Ahat, Sxy))^2 )

# covariance of original and canonical variables (V and X)
sum( (cov(V, X) - crossprod(Bhat, t(Sxy)))^2 )

# covariance of canonical variables (U and V)
rhomat <- cbind(diag(rho), matrix(0, p, q-p))
sum( (cov(U, V) - rhomat)^2 )
sum( (Sxy - crossprod(Ainv, rhomat) %*% Binv)^2 )
Sxyhat <- matrix(0, p, q)
for(j in 1:p) Sxyhat <- Sxyhat + rho[j] * outer(Ainv[j,], Binv[j,])
sum( (Sxy - Sxyhat)^2 )

# error of approximation matrices (with r=2)
Ainv <- solve(Ahat)
Binv <- solve(Bhat)
r <- 2
Ex <- Sx - crossprod(Ainv[1:r,])
Ey <- Sy - crossprod(Binv[1:r,])
Exy <- Sxy - crossprod(diag(rho[1:r]) %*% Ainv[1:r,], Binv[1:r,])

# get norms of error matrices
sqrt(mean(Ex^2))
sqrt(mean(Ey^2))
sqrt(mean(Exy^2))


#####   DECATHLON EXAMPLE: CORRELATION CCA   #####

# standardize data
Xs <- scale(X)
Ys <- scale(Y)

# canonical correlations of correlations (standardized data)
ccas <- cancor(Xs, Ys)

# cca (the normal way)
Sx <- cov(Xs)
Sy <- cov(Ys)
Sxy <- cov(Xs,Ys)
Sxeig <- eigen(Sx, symmetric=TRUE)
Sxisqrt <- Sxeig$vectors %*% diag(1/sqrt(Sxeig$values)) %*% t(Sxeig$vectors)
Syeig <- eigen(Sy, symmetric=TRUE)
Syisqrt <- Syeig$vectors %*% diag(1/sqrt(Syeig$values)) %*% t(Syeig$vectors)
Xmat <- Sxisqrt %*% Sxy %*% solve(Sy) %*% t(Sxy) %*% Sxisqrt
Ymat <- Syisqrt %*% t(Sxy) %*% solve(Sx) %*% Sxy %*% Syisqrt
Xeig <- eigen(Xmat, symmetric=TRUE)
Yeig <- eigen(Ymat, symmetric=TRUE)

# compare correlations (same)
cca$cor
sqrt(Xeig$values)
sqrt(Yeig$values[1:p])

# compare linear combinations (different?)
Ahat <- Sxisqrt %*% Xeig$vectors
Bhat <- Syisqrt %*% Yeig$vectors
sum((ccas$xcoef * sqrt(n-1) - Ahat)^2)
sum((ccas$ycoef[,1:p] * sqrt(n-1) - Bhat[,1:p])^2)

# note that the signing is arbitary!!
ccas$ycoef[,1:p] * sqrt(n-1)
Bhat[,1:p]
Bhat[,1:p] <- Bhat[,1:p] %*% diag(c(-1,1,-1,1))
sum((ccas$ycoef[,1:p] * sqrt(n-1) - Bhat[,1:p])^2)

# plot coefficients
dev.new(width=10, height=5, noRStudioGD=TRUE)
par(mfrow=c(1,2))
plot(Ahat[,1:2], xlab="A1 Coefficients", ylab="A2 Coefficients",
     type="n", main="X Coefficients", xlim=c(-2, 0.1), ylim=c(-1.5, 2))
text(Ahat[,1:2], labels=colnames(X))
abline(0,0,lty=3)
abline(v=0,lty=3)
plot(Bhat[,1:2], xlab="B1 Coefficients", ylab="B2 Coefficients",
     type="n", main="Y Coefficients", xlim=c(-0.5, 0.5), ylim=c(-1.1, 0.8))
text(Bhat[,1:2], labels=colnames(Y))
abline(0,0,lty=3)
abline(v=0,lty=3)

# define canonical variates
U <- Xs %*% Ahat
V <- Ys %*% Bhat

# canonical variable covariances
round(cov(U),4)
round(cov(V),4)
round(cov(U,V),4)

# covariance of original and canonical variables (U and Xs)
Ainv <- solve(Ahat)
sum( ( cov(U, Xs) - crossprod(Ahat, Sx) )^2 )
sum( ( Sx - crossprod(Ainv) )^2 )
Sxhat <- matrix(0, p, p)
for(j in 1:p) Sxhat <- Sxhat + outer(Ainv[j,], Ainv[j,])
sum( (Sx - Sxhat)^2 )

# covariance of original and canonical variables (V and Ys)
Binv <- solve(Bhat)
sum( ( cov(V, Ys) - crossprod(Bhat, Sy) )^2 )
sum( ( Sy - crossprod(Binv) )^2 )
Syhat <- matrix(0, q, q)
for(j in 1:q) Syhat <- Syhat + outer(Binv[j,], Binv[j,])
sum( (Sy - Syhat)^2 )

# covariance of original and canonical variables (U and Ys)
sum( (cov(U, Ys) - crossprod(Ahat, Sxy))^2 )

# covariance of original and canonical variables (V and Xs)
sum( (cov(V, Xs) - crossprod(Bhat, t(Sxy)))^2 )

# covariance of canonical variables (U and V)
rhomat <- cbind(diag(rho), matrix(0, p, q-p))
sum( (cov(U, V) - rhomat)^2 )
sum( (Sxy - crossprod(Ainv, rhomat) %*% Binv)^2 )
Sxyhat <- matrix(0, p, q)
for(j in 1:p) Sxyhat <- Sxyhat + rho[j] * outer(Ainv[j,], Binv[j,])
sum( (Sxy - Sxyhat)^2 )

# error of approximation matrices (with r=2)
Ainv <- solve(Ahat)
Binv <- solve(Bhat)
r <- 2
Ex <- Sx - crossprod(Ainv[1:r,])
Ey <- Sy - crossprod(Binv[1:r,])
Exy <- Sxy - crossprod(diag(rho[1:r]) %*% Ainv[1:r,], Binv[1:r,])

# get norms of error matrices
sqrt(mean(Ex^2))
sqrt(mean(Ey^2))
sqrt(mean(Exy^2))
