##########   Introduction to Linear Algebra
##########   Nathaniel E. Helwig (helwig@umn.edu)
##########   Updated: 04-Jan-2017


# create matrix columnwise and rowwise
x = 1:9
x
matrix(x,nrow=3,ncol=3)
matrix(x,nrow=3,ncol=3,byrow=TRUE)

# replication warning
x = 1:9
x
matrix(x,nrow=3,ncol=4)

# scalar and matrix multiplication
x = 1:9
y = 9:1
X = matrix(x,nrow=3,ncol=3)
Y = matrix(y,nrow=3,ncol=3)
X
Y
X * Y
X %*% Y

# scalar multiplication error
x = 1:6
y = 6:1
X = matrix(x,2,3)
Y = matrix(y,3,2)
X
Y
X*Y
X%*%Y

# matrix multiplication error
x = 1:6
y = 6:1
X = matrix(x,2,3)
Y = matrix(y,2,3)
X
Y
X * Y
X %*% Y

# transpose
X = matrix(1:6,2,3)
X
t(X)

# matrix dimensions
X = matrix(1:6,2,3)
X
dim(X)
dim(t(X))

# crossproduct
X = matrix(1:6,3,2)
Y = matrix(1:9,3,3)
crossprod(X,Y)
t(X) %*% Y

# transpose crossproduct
X = matrix(1:6,2,3)
Y = matrix(1:9,3,3)
tcrossprod(X,Y)
X %*% t(Y)

# row and column sums
X = matrix(1:6,2,3)
X
rowSums(X)
colSums(X)

# row and column means
X = matrix(1:6,2,3)
X
rowMeans(X)
colMeans(X)

# matrix diagonal to vector, vector to diagonal matrix, and identity matrix
X = matrix(1:4,2,2)
X
diag(X)
diag(1:3)
diag(2)

# eigenvalue decomposition
X = matrix(1:9,3,3)
X = crossprod(X)
xeig = eigen(X,symmetric=TRUE)
xeig$val
xeig$vec
Xhat = xeig$vec %*% diag(xeig$val) %*% t(xeig$vec)
sum( (X - Xhat)^2 )

# Cholesky decomposition
set.seed(1)
X = matrix(runif(9),3,3)
X = crossprod(X)
xchol = chol(X)
t(xchol)
Xhat = crossprod(xchol)
sum( (X - Xhat)^2 )

# singular value decomposition
X = matrix(1:6,3,2)
xsvd = svd(X)
xsvd$d
xsvd$u
xsvd$v
Xhat = xsvd$u %*% diag(xsvd$d) %*% t(xsvd$v)
sum( (X - Xhat)^2 )

# QR decomposition
X = matrix(1:6,3,2)
xqr = qr(X)
Q = qr.Q(xqr)
Q
R = qr.R(xqr)
R
Xhat = Q %*% R[,sort(xqr$pivot,index=TRUE)$ix]
sum( (X - Xhat)^2 )
