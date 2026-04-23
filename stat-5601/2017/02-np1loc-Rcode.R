##########   Nonparametric Location Tests: One-Sample
##########   Nathaniel E. Helwig (helwig@umn.edu)
##########   Updated: 04-Jan-2017


#####   DEFINE PATHS AND PACKAGES   #####

# load NSM3 package (for various functions)
if(!require(NSM3)){
  install.packages("NSM3")
  library("NSM3")
}

# load BSDA package (for SIGN.test function)
if(!require(BSDA)){
  install.packages("BSDA")
  library("BSDA")
}


#####   BACKGROUND INFORMATION   #####

# order and rank example (no ties)
x = c(3,12,11,18,14,10)
x
sort(x)
rank(x)

# order and rank example (with ties)
x = c(3,11,11,14,14,11)
x
sort(x)
rank(x)


#####   WILCOXON SIGNED RANK TEST   #####

# histogram of sum of independent variables (CLT)
set.seed(1)
dev.new(width=8,height=6,noRStudioGD=TRUE)
par(mfrow=c(2,2))
hist(replicate(5000,sum(runif(2))),main="n=2")
hist(replicate(5000,sum(runif(5))),main="n=5")
hist(replicate(5000,sum(runif(10))),main="n=10")
hist(replicate(5000,sum(runif(100))),main="n=100")

# normality test on sum of independent variables
set.seed(1)
shapiro.test(replicate(5000,sum(runif(2))))
shapiro.test(replicate(5000,sum(runif(100))))

# example 3.1:  hypothesis test
pre = c(1.83,0.50,1.62,2.48,1.68,1.88,1.55,3.06,1.30)
post = c(0.878,0.647,0.598,2.050,1.060,1.290,1.060,3.140,1.290)
z = post - pre
R = rank(abs(z))
wilcox.test(z,alternative="less")
wilcox.test(post,pre,alternative="less",paired=TRUE)

# example 3.1:  estimate theta
owa(pre,post)

# example 3.1:  CI for theta
wilcox.test(post,pre,alternative="less",paired=TRUE,conf.int=TRUE)


#####   FISHER SIGN TEST   #####

# example 3.1:  hypothesis test (revisted)
z = post - pre
SIGN.test(z,alternative="less")
SIGN.test(post,pre,alternative="less")

# example 3.1:  CI for theta (revisted)
zs = sort(z)
round(pbinom(0:9,9,1/2),4)
zs[7:8]
zs[7]+(zs[8]-zs[7])*(0.95-0.9102)/(0.9805-0.9102)


#####   UNIVARIATE SYMMETRY   #####

# symmetric data
set.seed(1)
x = rnorm(50)
dev.new(width=6,height=6,noRStudioGD=TRUE)
hist(x)
test = RFPW(x)
c(test$obs.stat, test$p.val)
2*(1 - pnorm(abs(test$obs.stat)))  # two-sided
pnorm(test$obs.stat)               # left-skew
1 - pnorm(test$obs.stat)           # right-skew

# asymmetric data
set.seed(1)
x = rchisq(50,df=3)
dev.new(width=6,height=6,noRStudioGD=TRUE)
hist(x)
test = RFPW(x)
c(test$obs.stat, test$p.val)
2*(1 - pnorm(abs(test$obs.stat)))  # two-sided
pnorm(test$obs.stat)               # left-skew
1 - pnorm(test$obs.stat)           # right-skew


#####   BIVARIATE SYMMETRY   #####

# example 3.1:  exchangeability test
pre = c(1.83,0.50,1.62,2.48,1.68,1.88,1.55,3.06,1.30)
post = c(0.878,0.647,0.598,2.050,1.060,1.290,1.060,3.140,1.290)
HollBivSym(pre,post)
set.seed(1)
pHollBivSym(pre,post)
