##########   Multivariate Linear Regression
##########   Nathaniel E. Helwig (helwig@umn.edu)
##########   Updated: 16-Jan-2017


#####   DEFINE PATHS AND PACKAGES   #####

# load car package (for various functions)
if(!require(car)){
  install.packages("car")
  library("car")
}


#####   MULTIPLE LINEAR REGRESSION   #####

# fit multiple linear regression model
data(mtcars)
head(mtcars)
mtcars$cyl <- factor(mtcars$cyl)
mod <- lm(mpg ~ cyl + am + carb, data=mtcars)
coef(mod)

# sum-of-squares information
anova(mod)
Anova(mod, type=2)
Anova(mod, type=3)

# summarize results
smod <- summary(mod)
names(smod)
summary(mod)$r.squared
summary(mod)$adj.r.squared

# get mean-squared error in 3 ways
n <- length(mtcars$mpg)
p <- length(coef(mod)) - 1
smod$sigma^2
sum((mod$residuals)^2) / (n - p - 1)
sum((mtcars$mpg - mod$fitted.values)^2) / (n - p - 1)

# get MLE of error variance
smod$sigma^2 * (n-p-1) / n

# inference information
summary(mod)
confint(mod)

# confidence interval
newdata <- data.frame(cyl=factor(6, levels=c(4,6,8)), am=1, carb=4)
predict(mod, newdata, interval="confidence")

# prediction interval
newdata <- data.frame(cyl=factor(6, levels=c(4,6,8)), am=1, carb=4)
predict(mod, newdata, interval="prediction")

# confidence region
dev.new(height=4,width=8,noRStudioGD=TRUE)
par(mfrow=c(1,2))
confidenceEllipse(mod,c(2,3),levels=.9,xlim=c(-11,3),ylim=c(-14,0),
                  main=expression(alpha*" = "*.1),cex.main=2)
confidenceEllipse(mod,c(2,3),levels=.99,xlim=c(-11,3),ylim=c(-14,0),
                  main=expression(alpha*" = "*.01),cex.main=2)


#####   MULTIVARIATE LINEAR REGRESSION   #####

# fit multivariate multiple linear regression model
data(mtcars)
head(mtcars)
mtcars$cyl <- factor(mtcars$cyl)
Y <- as.matrix(mtcars[,c("mpg","disp","hp","wt")])
mvmod <- lm(Y ~ cyl + am + carb, data=mtcars)
coef(mvmod) 

# get SSCP matrices
ybar <- colMeans(Y)
n <- nrow(Y)
m <- ncol(Y)
Ybar <-  matrix(ybar, n, m, byrow=TRUE)
SSCP.T <- crossprod(Y - Ybar)
SSCP.R <- crossprod(mvmod$fitted.values - Ybar)
SSCP.E <- crossprod(Y - mvmod$fitted.values)
SSCP.T
SSCP.R + SSCP.E

# estimated error covariance matrix
n <- nrow(Y)
p <- nrow(coef(mvmod)) - 1
SSCP.E <- crossprod(Y - mvmod$fitted.values)
SigmaHat <- SSCP.E / (n - p - 1)
SigmaTilde <- SSCP.E / n
SigmaHat
SigmaTilde

# inference for coefficients
mvsum <- summary(mvmod)
mvsum
mvsum[[1]]
mvsum[[3]]

# testing a reduced model (easy way)
mvmod0 <- lm(Y ~ am + carb, data=mtcars)
anova(mvmod, mvmod0, test="Wilks")
anova(mvmod, mvmod0, test="Pillai")
anova(mvmod, mvmod0, test="Hotelling-Lawley")
anova(mvmod, mvmod0, test="Roy")

# testing a reduced model (by hand)
Etilde <- n * SigmaTilde
SigmaTilde1 <- crossprod(Y - mvmod0$fitted.values) / n
Htilde <- n * (SigmaTilde1 - SigmaTilde)
HEi <- Htilde %*% solve(Etilde)
HEi.values <- eigen(HEi)$values
c(Wilks=prod(1/(1+HEi.values)), Pillai=sum(HEi.values / (1 + HEi.values)))

# confidence interval
newdata <- data.frame(cyl=factor(6, levels=c(4,6,8)), am=1, carb=4)
predict(mvmod, newdata, interval="confidence")

# prediction interval
newdata <- data.frame(cyl=factor(6, levels=c(4,6,8)), am=1, carb=4)
predict(mvmod, newdata, interval="prediction")

# confidence and prediction interval function
pred.mlm <- function(object, newdata, level=0.95,
                     interval = c("confidence", "prediction")){
  form <- as.formula(paste("~",as.character(formula(object))[3]))
  xnew <- model.matrix(form, newdata)
  fit <- predict(object, newdata)
  Y <- model.frame(object)[,1]
  X <- model.matrix(object)
  n <- nrow(Y)
  m <- ncol(Y)
  p <- ncol(X) - 1
  sigmas <- colSums((Y - object$fitted.values)^2) / (n - p - 1)
  fit.var <- diag(xnew %*% tcrossprod(solve(crossprod(X)), xnew))
  if(interval[1]=="prediction") fit.var <- fit.var + 1
  const <- qf(level, df1=m, df2=n-p-m) * m * (n - p - 1) / (n - p - m) 
  vmat <- (n/(n-p-1)) * outer(fit.var, sigmas)
  lwr <- fit - sqrt(const) * sqrt(vmat)
  upr <- fit + sqrt(const) * sqrt(vmat)
  if(nrow(xnew)==1L){
    ci <- rbind(fit, lwr, upr)
    rownames(ci) <- c("fit", "lwr", "upr")
  } else {
    ci <- array(0, dim=c(nrow(xnew), m, 3))
    dimnames(ci) <- list(1:nrow(xnew), colnames(Y), c("fit", "lwr", "upr") )
    ci[,,1] <- fit
    ci[,,2] <- lwr
    ci[,,3] <- upr
  }
  ci
}

# confidence interval
newdata <- data.frame(cyl=factor(6, levels=c(4,6,8)), am=1, carb=4)
pred.mlm(mvmod, newdata)

# prediction interval
newdata <- data.frame(cyl=factor(6, levels=c(4,6,8)), am=1, carb=4)
pred.mlm(mvmod, newdata, interval="prediction")

# confidence interval (multiple new observations)
newdata <- data.frame(cyl=factor(c(4,6,8), levels=c(4,6,8)), am=c(0,1,1), carb=c(2,4,6))
pred.mlm(mvmod, newdata)

# prediction interval (multiple new observations)
newdata <- data.frame(cyl=factor(c(4,6,8), levels=c(4,6,8)), am=c(0,1,1), carb=c(2,4,6))
pred.mlm(mvmod, newdata, interval="prediction")
