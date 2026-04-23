##########   Data, Covariance, and Correlation Matrix
##########   Nathaniel E. Helwig (helwig@umn.edu)
##########   Updated: 16-Jan-2017


#####   DEFINE PATHS AND PACKAGES   #####

# load car package (for scatterplot functions)
if(!require(car)){
  install.packages("car")
  library("car")
}

# load scatterplot3d package (for scatterplot3d function)
if(!require(scatterplot3d)){
  install.packages("scatterplot3d")
  library("scatterplot3d")
}

# load bigsplines package (for imagebar function)
if(!require(bigsplines)){
  install.packages("bigsplines")
  library("bigsplines")
}

# load corrplot package (for corrplot function)
if(!require(corrplot)){
  install.packages("corrplot")
  library("corrplot")
}

# load RColorBrewer package (for brewer.pal function)
if(!require(RColorBrewer)){
  install.packages("RColorBrewer")
  library("RColorBrewer")
}


#####   THE DATA MATRIX   #####

# load mtcars data
data(mtcars)
class(mtcars)
dim(mtcars)
head(mtcars)
X <- as.matrix(mtcars)
class(X)

# get row means (3 ways)
rowMeans(X)[1:3]
c(mean(X[1,]),mean(X[2,]),mean(X[3,]))
apply(X,1,mean)[1:3]

# get column means (3 ways)
colMeans(X)[1:3]
c(mean(X[,1]),mean(X[,2]),mean(X[,3]))
apply(X,2,mean)[1:3]

# get column medians
apply(X,2,median)[1:3]
c(median(X[,1]),median(X[,2]),median(X[,3]))

# get column ranges
apply(X,2,range)[,1:3]
cbind(range(X[,1]),range(X[,2]),range(X[,3]))


#####   THE COVARIANCE MATRIX   #####

# calculate covariance matrix
n <- nrow(X)
C <- diag(n) - matrix(1/n, n, n)
Xc <- C %*% X
S <- t(Xc) %*% Xc / (n-1)
S[1:3,1:6]
# or #
Xc <- scale(X, center=TRUE, scale=FALSE)
S <- t(Xc) %*% Xc / (n-1)
S[1:3,1:6]

# calculate covariance matrix
S <- cov(X)
dim(S)

# check variance
S[1,1]
var(X[,1])
sum((X[,1]-mean(X[,1]))^2) / (n-1)

# check covariance
S[1:3,1:6]


#####   THE CORRELATION MATRIX   #####

# calculate correlation matrix
n <- nrow(X)
C <- diag(n) - matrix(1/n, n, n)
D <- diag(apply(X, 2, sd))
Xs <- C %*% X %*% solve(D)
R <- t(Xs) %*% Xs / (n-1)
R[1:3,1:6]
# or #
Xs <- scale(X, center=TRUE, scale=TRUE)
R <- t(Xs) %*% Xs / (n-1)
R[1:3,1:6]

# calculate correlation matrix
R <- cor(X)
dim(R)

# check correlation of mpg and cyl
R[1,2]
cor(X[,1],X[,2])
cov(X[,1],X[,2]) / (sd(X[,1]) * sd(X[,2]))

# check correlations
R[1:3,1:6]


#####   MISCELLANEOUS TOPICS   #####

# crossproduct calculations
X <- matrix(rnorm(2*3),2,3)
Y <- matrix(rnorm(2*3),2,3)
t(X) %*% Y
crossprod(X, Y)
X %*% t(Y)
tcrossprod(X, Y)

# vectorizing a matrix
X <- matrix(1:6,2,3)
X
c(X)
c(t(X))

# vec operator: property 2
set.seed(1)
a <- runif(3)
b <- runif(4)
c(tcrossprod(a,b))
kronecker(b,a)

# vec operator: property 3
set.seed(1)
A <- matrix(runif(10),5,2)
B <- matrix(runif(10),5,2)
crossprod(c(A), c(B))
sum(diag(crossprod(A, B)))

# vec operator: property 4
set.seed(1)
A <- matrix(runif(12),3,4)
B <- matrix(runif(20),4,5)
C <- matrix(runif(20),5,4)
c(A %*% B %*% C)
kronecker(t(C), A) %*% c(B)

# Kronecker products
X <- matrix(1:4,2,2)
Y <- matrix(5:10,2,3)
kronecker(X, Y)

# Kronecker product: property 2
set.seed(1)
A <- matrix(runif(12),3,4)
B <- matrix(runif(20),4,5)
sum(( t(kronecker(A,B)) - kronecker(t(A),t(B)) )^2)

# Kronecker product: property 3
set.seed(1)
a <- matrix(runif(3),3,1)
b <- matrix(runif(4),4,1)
sum(( kronecker(t(a),b) - tcrossprod(b,a) )^2)
sum(( kronecker(b,t(a)) - tcrossprod(b,a) )^2)

# Kronecker product: property 4
set.seed(1)
A <- matrix(runif(9),3,3)
B <- matrix(runif(16),4,4)
sum(diag(kronecker(A,B)))
sum(diag(A)) * sum(diag(B))

# Kronecker product: property 5
set.seed(1)
A <- matrix(runif(9),3,3)
B <- matrix(runif(16),4,4)
sum(( solve(kronecker(A,B)) - kronecker(solve(A), solve(B)) )^2)

# Kronecker product: property 6
set.seed(1)
A <- cbind(matrix(runif(9),3,3),1,1)
B <- cbind(matrix(runif(16),4,4),2,2)
pinv <- function(x){
  xsvd <- svd(x)
  xsvd$v %*% diag(1/xsvd$d) %*% t(xsvd$u)
}
sum(( pinv(kronecker(A,B)) - kronecker(pinv(A), pinv(B)) )^2)

# Kronecker product: property 7
set.seed(1)
A <- crossprod(matrix(runif(9),3,3))
B <- crossprod(matrix(runif(16),4,4))
sum(log(eigen(kronecker(A,B))$values))
ncol(B) * sum(log(eigen(A)$values)) + ncol(A) * sum(log(eigen(B)$values))

# Kronecker product: property 8
set.seed(1)
A <- cbind(matrix(runif(9),3,3),1,1)
B <- cbind(matrix(runif(16),4,4),2,2)
qr(kronecker(A,B))$rank
qr(A)$rank * qr(B)$rank

# Kronecker product: property 9
set.seed(1)
A <- matrix(runif(9),3,3)
B <- matrix(runif(16),4,4)
C <- matrix(runif(16),4,4)
sum(( kronecker(A, B + C) - (kronecker(A,B) + kronecker(A,C)) )^2)

# Kronecker product: property 10
set.seed(1)
A <- matrix(runif(9),3,3)
B <- matrix(runif(9),3,3)
C <- matrix(runif(16),4,4)
sum(( kronecker(A + B, C) - (kronecker(A,C) + kronecker(B,C)) )^2)

# Kronecker product: property 11
set.seed(1)
A <- matrix(runif(9),3,3)
B <- matrix(runif(16),4,4)
C <- matrix(runif(25),5,5)
sum(( kronecker(kronecker(A, B), C) - kronecker(A, kronecker(B, C)) )^2)

# Kronecker product: property 12
set.seed(1)
A <- matrix(runif(9),3,3)
B <- matrix(runif(16),4,4)
C <- matrix(runif(12),3,4)
D <- matrix(runif(20),4,5)
sum(( kronecker(A, B) %*% kronecker(C, D) - kronecker(A %*% C, B %*% D) )^2)

# two versions of scatterplot
dev.new(width=6,height=6,noRStudioGD=TRUE)
plot(mtcars$hp, mtcars$mpg, xlab="HP", ylab="MPG")
dev.new(width=6,height=6,noRStudioGD=TRUE)
scatterplot(mtcars$hp, mtcars$mpg, xlab="HP", ylab="MPG")

# two versions of scatterplot matrix
cylint <- as.integer(factor(mtcars$cyl))
dev.new(width=6,height=6,noRStudioGD=TRUE)
pairs(~mpg+disp+hp+wt, data=mtcars, col=cylint, pch=cylint)
dev.new(width=6,height=6,noRStudioGD=TRUE)
scatterplotMatrix(~mpg+disp+hp+wt|cyl, data=mtcars)

# three-dimensional scatterplot
dev.new(width=6,height=6,noRStudioGD=TRUE)
sp3d <- scatterplot3d(mtcars$hp, mtcars$wt, mtcars$mpg, 
                      type="h", color=cylint, pch=cylint, 
                      xlab="HP", ylab="WT", zlab="MPG")
fitmod <- lm(mpg ~ hp + wt, data=mtcars)
sp3d$plane3d(fitmod)

# color image (heat map)
fitmod <- lm(mpg ~ hp + wt, data=mtcars)
hpseq <- seq(50, 330, by=20)
wtseq <- seq(1.5, 5.4, length=15)
newdata <- expand.grid(hp=hpseq, wt=wtseq)
fit <- predict(fitmod, newdata)
fitmat <- matrix(fit, 15, 15)
dev.new(width=6,height=6,noRStudioGD=TRUE)
image(hpseq, wtseq, fitmat, xlab="HP", ylab="WT")
dev.new(width=6,height=6,noRStudioGD=TRUE)
imagebar(hpseq, wtseq, fitmat, xlab="HP", ylab="WT", zlab="MPG", col=heat.colors(12), ncolor=12)

# visualizing correlation matrices (corrplot)
cmat <- cor(mtcars)
dev.new(width=6,height=6,noRStudioGD=TRUE)
corrplot(cmat, method="circle")
dev.new(width=6,height=6,noRStudioGD=TRUE)
corrplot.mixed(cmat, lower="number", upper="ellipse")

# visualizing correlation matrices (heat map)
cmat <- cor(mtcars)
p <- nrow(cmat)
dev.new(width=6,height=6,noRStudioGD=TRUE)
imagebar(1:p, 1:p, cmat[,p:1], axes=F, zlim=c(-1,1), xlab="", ylab="",
         col=brewer.pal(7, "Spectral"))
axis(1, 1:p, labels=rownames(cmat))
axis(2, p:1, labels=colnames(cmat))
dev.new(width=6,height=6,noRStudioGD=TRUE)
imagebar(1:p, 1:p, cmat[,p:1], axes=F, zlim=c(-1,1), xlab="", ylab="",
         col=brewer.pal(7, "RdBu"))
axis(1, 1:p, labels=rownames(cmat))
axis(2, p:1, labels=colnames(cmat))
for(k in 1:p){
  for(j in 1:k){
    if(j < k) text(j, p+1-k, labels=round(cmat[j,k],2), cex=0.75)
  }
}
