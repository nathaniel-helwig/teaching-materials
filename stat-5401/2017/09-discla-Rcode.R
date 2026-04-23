##########   Discrimination and Classification
##########   Nathaniel E. Helwig (helwig@umn.edu)
##########   Updated: 16-Mar-2017


### load 'MASS' package
if(!require(MASS)){
  install.packages("MASS")
  library(MASS)
}

### load 'klaR' package
if(!require(klaR)){
  install.packages("klaR")
  library(klaR)
}


#####   LINEAR DISCRIMINANT ANALYSIS   #####

# look at Fisher's iris data
head(iris)
colMeans(iris[iris$Species=="setosa",1:4])
colMeans(iris[iris$Species=="versicolor",1:4])
colMeans(iris[iris$Species=="virginica",1:4])
p <- 4L
g <- 3L

# make pooled covariances matrix
Sp <- matrix(0, p, p)
nx <- rep(0, g)
lev <- levels(iris$Species)
for(k in 1:g){
  x <- iris[iris$Species==lev[k],1:p]
  nx[k] <- nrow(x)
  Sp <- Sp + cov(x) * (nx[k] - 1)
}
Sp <- Sp / (sum(nx) - g)
round(Sp,3)

# fit lda model
ldamod <- lda(Species ~ ., data=iris, prior=rep(1/3, 3))
names(ldamod)

# check the LDA coefficients/scalings
ldamod$scaling
crossprod(ldamod$scaling, Sp) %*% ldamod$scaling

# create the (centered) discriminant scores
mu.k <- ldamod$means
mu <- colMeans(mu.k)
dscores <- scale(iris[,1:p], center=mu, scale=F) %*% ldamod$scaling
sum((dscores - predict(ldamod)$x)^2)

# plot the scores and coefficients
spid <- as.integer(iris$Species)
dev.new(width=10, height=5, noRStudioGD=TRUE)
par(mfrow=c(1,2))
plot(dscores, xlab="LD1", ylab="LD2", pch=spid, col=spid,
     main="Discriminant Scores", xlim=c(-10, 10), ylim=c(-3, 3))
abline(h=0, lty=3)
abline(v=0, lty=3)
legend("top",lev,pch=1:3,col=1:3,bty="n")
plot(ldamod$scaling, xlab="LD1", ylab="LD2", type="n",
     main="Discriminant Coefficients", xlim=c(-4, 3), ylim=c(-1, 3))
text(ldamod$scaling, labels=rownames(ldamod$scaling))
abline(h=0, lty=3)
abline(v=0, lty=3)

# visualize the LDA paritions
dev.new(width=5, height=5, noRStudioGD=TRUE)
species <- factor(rep(c("s","c","v"), each=50))
partimat(x=dscores[,2:1], grouping=species, method="lda")

# visualize the LDA paritions (for all pairs)
dev.new(width=10, height=5, noRStudioGD=TRUE)
species <- factor(rep(c("s","c","v"), each=50))
partimat(x=iris[,1:4], grouping=species, method="lda")

# make confusion matrix (and APER)
confusion <- table(iris$Species, predict(ldamod)$class)
confusion
n <- sum(confusion)
aper <- (n - sum(diag(confusion))) / n
aper

# use CV to get expected AER
ldamodCV <- lda(Species ~ ., data=iris, prior=rep(1/3, 3), CV=TRUE)
confusionCV <- table(iris$Species, ldamodCV$class)
confusionCV
eaer <- (n - sum(diag(confusionCV))) / n
eaer

# split into separate matrices for each flower
Xs <- subset(iris, Species=="setosa")
Xc <- subset(iris, Species=="versicolor")
Xv <- subset(iris, Species=="virginica")

# split into training and testing
set.seed(1)
sid <- sample.int(n=50, size=35)
cid <- sample.int(n=50, size=35)
vid <- sample.int(n=50, size=35)
Xtrain <- rbind(Xs[sid,], Xc[cid,], Xv[vid,])
Xtest <- rbind(Xs[-sid,], Xc[-cid,], Xv[-vid,])

# fit lda to training and evaluate on testing
ldatrain <- lda(Species ~ ., data=Xtrain, prior=rep(1/3, 3))
confusionTest <- table(Xtest$Species, predict(ldatrain, newdata=Xtest)$class)
confusionTest
n <- sum(confusionTest)
aer <- (n - sum(diag(confusionTest))) / n
aer

# split into training and testing (100 splits)
nrep <- 100
aer <- rep(0, nrep)
set.seed(1)
for(k in 1:nrep){
  cat("rep:",k,"\n")
  sid <- sample.int(n=50, size=35)
  cid <- sample.int(n=50, size=35)
  vid <- sample.int(n=50, size=35)
  Xtrain <- rbind(Xs[sid,], Xc[cid,], Xv[vid,])
  Xtest <- rbind(Xs[-sid,], Xc[-cid,], Xv[-vid,])
  ldatrain <- lda(Species ~ ., data=Xtrain, prior=rep(1/3, 3))
  confusionTest <- table(Xtest$Species, predict(ldatrain, newdata=Xtest)$class)
  confusionTest
  n <- sum(confusionTest)
  aer[k] <- (n - sum(diag(confusionTest))) / n
}
mean(aer)


#####   QUADRATIC DISCRIMINANT ANALYSIS   #####

# fit qda model
qdamod <- qda(Species ~ ., data=iris, prior=rep(1/3, 3))
names(qdamod)

# check the QDA coefficients/scalings
dim(qdamod$scaling)
dnames <- dimnames(qdamod$scaling)
dnames
round(crossprod(qdamod$scaling[,,1], cov(Xs[,1:p])) %*% qdamod$scaling[,,1], 4)
round(crossprod(qdamod$scaling[,,2], cov(Xc[,1:p])) %*% qdamod$scaling[,,2], 4)
round(crossprod(qdamod$scaling[,,3], cov(Xv[,1:p])) %*% qdamod$scaling[,,3], 4)

# visualize the QDA paritions
dev.new(width=10, height=5, noRStudioGD=TRUE)
partimat(Species ~ ., data=iris, method="qda")

# make confusion matrix (and APER)
confusion <- table(iris$Species, predict(qdamod)$class)
confusion
n <- sum(confusion)
aper <- (n - sum(diag(confusion))) / n
aper

# use CV to get expected AER
qdamodCV <- qda(Species ~ ., data=iris, prior=rep(1/3, 3), CV=TRUE)
confusionCV <- table(iris$Species, qdamodCV$class)
confusionCV
eaer <- (n - sum(diag(confusionCV))) / n
eaer

# split into training and testing
set.seed(1)
sid <- sample.int(n=50, size=35)
cid <- sample.int(n=50, size=35)
vid <- sample.int(n=50, size=35)
Xtrain <- rbind(Xs[sid,], Xc[cid,], Xv[vid,])
Xtest <- rbind(Xs[-sid,], Xc[-cid,], Xv[-vid,])

# fit qda to training and evaluate on testing
qdatrain <- qda(Species ~ ., data=Xtrain, prior=rep(1/3, 3))
confusionTest <- table(Xtest$Species, predict(qdatrain, newdata=Xtest)$class)
confusionTest
n <- sum(confusionTest)
aer <- (n - sum(diag(confusionTest))) / n
aer

# split into training and testing (100 splits)
nrep <- 100
aer <- rep(0, nrep)
set.seed(1)
for(k in 1:nrep){
  cat("rep:",k,"\n")
  sid <- sample.int(n=50, size=35)
  cid <- sample.int(n=50, size=35)
  vid <- sample.int(n=50, size=35)
  Xtrain <- rbind(Xs[sid,], Xc[cid,], Xv[vid,])
  Xtest <- rbind(Xs[-sid,], Xc[-cid,], Xv[-vid,])
  qdatrain <- qda(Species ~ ., data=Xtrain, prior=rep(1/3, 3))
  confusionTest <- table(Xtest$Species, predict(qdatrain, newdata=Xtest)$class)
  confusionTest
  n <- sum(confusionTest)
  aer[k] <- (n - sum(diag(confusionTest))) / n
}
mean(aer)

# visualize LDA and QDA results via PCA
ldaid <- as.integer(predict(ldamod)$class)
qdaid <- as.integer(predict(qdamod)$class)
pcamod <- princomp(iris[,1:4])
dev.new(width=10, height=5, noRStudioGD=TRUE)
par(mfrow=c(1,2))
plot(pcamod$scores[,1:2], xlab="PC1", ylab="PC2", pch=ldaid, col=ldaid,
     main="LDA Results", xlim=c(-4, 4), ylim=c(-2, 2))
legend("topright",lev,pch=1:3,col=1:3,bty="n")
abline(h=0,lty=3)
abline(v=0,lty=3)
plot(pcamod$scores[,1:2], xlab="PC1", ylab="PC2", pch=qdaid, col=qdaid,
     main="QDA Results", xlim=c(-4, 4), ylim=c(-2, 2))
legend("topright",lev,pch=1:3,col=1:3,bty="n")
abline(h=0,lty=3)
abline(v=0,lty=3)

