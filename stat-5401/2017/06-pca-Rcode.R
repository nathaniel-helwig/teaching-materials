##########   Principal Components Analysis
##########   Nathaniel E. Helwig (helwig@umn.edu)
##########   Updated: 16-Mar-2017


#####   BACKGROUND   #####

# 2D rotate points
rotmat2d <- function(theta){
  matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)),2,2)
}
x <- seq(-2,2,length=11)
y <- 4*exp(-x^2) - 2
xy <- cbind(x,y)
rang <- c(30,45,60,90,180)
dev.new(width=12,height=8,noRStudioGD=TRUE)
par(mfrow=c(2,3))
plot(x,y,xlim=c(-3,3),ylim=c(-3,3),main="No Rotation")
text(x,y,labels=letters[1:11],cex=1.5)
for(j in 1:5){
  rmat <- rotmat2d(rang[j]*2*pi/360)
  xyrot <- xy%*%rmat
  plot(xyrot,xlim=c(-3,3),ylim=c(-3,3))
  text(xyrot,labels=letters[1:11],cex=1.5)
  title(paste(rang[j]," degrees"))
}


#####   PCA IN PRACTICE   #####

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

# PCA on covariance matrix (default)
n <- nrow(decathlon)
pcaCOV <- princomp(x=decathlon[,1:10])
names(pcaCOV)

# resign PCA solution
pcsign <- sign(colSums(pcaCOV$loadings^3))
pcaCOV$loadings <- pcaCOV$loadings %*% diag(pcsign)
pcaCOV$scores <- pcaCOV$scores %*% diag(pcsign)

# Note: R uses MLE of covariance matrix
sum((pcaCOV$sdev - sqrt(eigen((n-1)/n*cov(decathlon[,1:10]), symmetric=TRUE)$values))^2)

# PCA on correlation matrix
pcaCOR <- princomp(x=decathlon[,1:10], cor=TRUE)
names(pcaCOR)

# resign PCA solution
pcsign <- sign(colSums(pcaCOR$loadings^3))
pcaCOR$loadings <- pcaCOR$loadings %*% diag(pcsign)
pcaCOR$scores <- pcaCOR$scores %*% diag(pcsign)

# Note: PC standard deviations are square roots of correlation matrix eigenvalues
sum((pcaCOR$sdev - sqrt(eigen(cor(decathlon[,1:10]), symmetric=TRUE)$values))^2)

# plot linear combination weights (PC loadings)
dev.new(width=10, height=5, noRStudioGD=TRUE)
par(mfrow=c(1,2))
plot(pcaCOV$loadings[,1:2], xlab="PC1 Loaings", ylab="PC2 Loadings",
     type="n", main="PCA of Covariance Matrix", xlim=c(-0.15, 1.15), ylim=c(-0.15, 1.15))
text(pcaCOV$loadings[,1:2], labels=colnames(decathlon))
plot(pcaCOR$loadings[,1:2], xlab="PC1 Loaings", ylab="PC2 Loadings",
     type="n", main="PCA of Correlation Matrix", xlim=c(0.05, 0.45), ylim=c(-0.5,0.6))
text(pcaCOR$loadings[,1:2], labels=colnames(decathlon))

# check correlation of PC scores w/ overall decathlon score
round(cor(decathlon$score, pcaCOV$scores),3)
round(cor(decathlon$score, pcaCOR$scores),3)
dev.new(width=10, height=5, noRStudioGD=TRUE)
par(mfrow=c(1,2))
plot(pcaCOV$scores[,2], decathlon$score, xlab="PC2 Score", 
     ylab="Decathlon Score", main="PCA of Covariance Matrix")
plot(pcaCOR$scores[,1], decathlon$score, xlab="PC1 Score", 
     ylab="Decathlon Score", main="PCA of Correlation Matrix")

# scree plot for covariance and correlation solutions
dev.new(width=10, height=5, noRStudioGD=TRUE)
par(mfrow=c(1,2))
plot(1:10, pcaCOV$sdev^2, type="b", xlab="# PCs", ylab="Variance of PC", 
     main="PCA of Covariance Matrix")
plot(1:10, pcaCOR$sdev^2, type="b", xlab="# PCs", ylab="Variance of PC", 
     main="PCA of Correlation Matrix")
