##########   Factor Analysis
##########   Nathaniel E. Helwig (helwig@umn.edu)
##########   Updated: 16-Mar-2017


#####   FACTOR ROTATION   #####

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


#####   DECATHLON EXAMPLE   #####

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

# scree plot for factor analysis results
famods <- vector("list", 6)
for(k in 1:6){
  famods[[k]] <- factanal(x=decathlon[,1:10], factors=k)
}
vafs <- sapply(famods, function(x) sum(x$loadings^2)) / nrow(famods[[1]]$loadings)
vaf.scree <- vafs - c(0, vafs[1:5])
dev.new(width=10, height=5, noRStudioGD=TRUE)
plot(1:6, vaf.scree, type="b", xlab="# Factors", ylab="Proportion of Variance", 
     main="FA Scree Plot")

# FA on correlation matrix (w/ varimax rotation)
famod <- factanal(x=decathlon[,1:10], factors=2)
names(famod)

# plot factor loadings and uniquenesses
dev.new(width=10, height=5, noRStudioGD=TRUE)
par(mfrow=c(1,2))
plot(famod$loadings, xlab="F1 Loadings", ylab="F2 Loadings",
     type="n", main="Factor Loadings", xlim=c(-0.35, 1), ylim=c(-0.35, 1))
text(famod$loadings, labels=colnames(decathlon))
plot(famod$uniquenesses, xlab=expression("Variable ("*X[j]*")"), 
     ylab=expression("Uniqueness ("*psi[j]*")"),
     type="n", main="Factor Uniquenesses", xlim=c(0.5,10.5))
text(famod$uniquenesses, labels=colnames(decathlon))

# refit model and get FA scores (NOT GOOD IDEA!!)
famodW <- factanal(x=decathlon[,1:10], factors=2, scores="Bartlett")
famodR <- factanal(x=decathlon[,1:10], factors=2, scores="regression")
round(cor(famodW$scores, famodR$scores), 4)

# check correlation of FA scores w/ overall decathlon score
round(cor(decathlon$score, famodR$scores), 4)
round(cor(decathlon$score, famodW$scores), 4)
dev.new(width=10, height=5, noRStudioGD=TRUE)
par(mfrow=c(1,2))
plot(famodW$scores[,1], decathlon$score, xlab="F1 Score", 
     ylab="Decathlon Score", main="Weighted Least Squares Method")
plot(famodR$scores[,1], decathlon$score, xlab="F1 Score", 
     ylab="Decathlon Score", main="Regression Method")

# try promax rotation
famod.promax <- promax(famod$loadings)
dev.new(width=10, height=5, noRStudioGD=TRUE)
par(mfrow=c(1,2))
plot(famod$loadings, xlab="F1 Loadings", ylab="F2 Loadings",
     type="n", main="Varimax Factor Loadings", xlim=c(-0.35, 1), ylim=c(-0.35, 1))
text(famod$loadings, labels=colnames(decathlon))
abline(0,0,lty=3)
abline(v=0,lty=3)
plot(famod.promax$loadings, xlab="F1 Loadings", ylab="F2 Loadings",
     type="n", main="Promax Factor Loadings", xlim=c(-0.35, 1), ylim=c(-0.35, 1))
text(famod.promax$loadings, labels=colnames(decathlon))
abline(0,0,lty=3)
abline(v=0,lty=3)

# compare loadings
oldFAloadings <- famod$loadings
newFAloadings <- famod$loadings %*% famod.promax$rotmat
sum((newFAloadings - famod.promax$loadings)^2)

# compare reproduced data before and after rotation
oldFAscores <- famodR$scores
newFAscores <- oldFAscores %*% t(solve(famod.promax$rotmat))
Xold <- tcrossprod(oldFAscores, oldFAloadings)
Xnew <- tcrossprod(newFAscores, newFAloadings)
sum((Xold - Xnew)^2)

# population and sample factor score covariance matrix (after rotation)
tcrossprod(solve(famod.promax$rotmat))    # population
cor(newFAscores)                          # sample
