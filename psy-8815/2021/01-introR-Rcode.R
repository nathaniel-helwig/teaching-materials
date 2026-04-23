##########   Introduction to R and Programming
##########   Nathaniel E. Helwig (helwig@umn.edu)
##########   Updated: 04-Jan-2017


#####   BASIC CALCULATIONS   #####

# addition and subtraction
3+2
3-2

# multiplication and division
3/2
3*2

# exponents
3^2
2^3

# constants
pi
exp(1)

# infinity
Inf
1+Inf

# empty
NULL
1+NULL

# missing
NA
1+NA

# machine epsilon
.Machine$double.eps
0>.Machine$double.eps

# storing and manipulating objects
x=3
y=2
x+y
x*y
x^y


#####   USING R FUNCTIONS   #####

# combine
c(1,3,-2)
c("a","a","b","b","a")

# sum and mean
sum(c(1,3,-2))
mean(c(1,3,-2))

# variance and standard deviation
var(c(1,3,-2))
sd(c(1,3,-2))

# minimum and maximum
min(c(1,3,-2))
max(c(1,3,-2))

# correlation and covariance
x=c(1,3,4,6,8)
y=c(2,3,5,7,9)
cor(x,y)
cov(x,y)

# combine as columns and rows
cbind(x,y)
rbind(x,y)


#####   OBJECT CLASSES IN R   #####

# numeric class
x=c(1,3,-2)
x
class(x)

# integer class
x=c(1L,3L,-2L)
x
class(x)

# character class
x=c("a","a","b")
x
class(x)

# factor class
x=factor(c("a","a","b"))
x
class(x)

# matrix class
x=c(1,3,-2)
y=c(2,0,7)
z=cbind(x,y)
z
class(z)
class(z[,1])
class(z[,2])

# data frame class
x=c(1,3,-2)
y=c("a","a","b")
z=data.frame(x,y)
z
class(z)
class(z$x)
class(z$y)

# print and summary (class customized)
x=c(1,3,-2)
y=factor(c("a","a","b"))
print(x)
print(y)
summary(x)
summary(y)

# range (class customized)
x=c(1,3,-2)
y=factor(c("a","a","b"))
range(x)
range(y)


#####   NORMAL DISTRIBUTION   #####

# standard normal density function
dnorm(-4)
dnorm(-2)
dnorm(0)
dnorm(2)
dnorm(4)

# plot standard normal density function
dev.new(height=8,width=8,noRStudioGD=TRUE)
x=seq(-4,4,by=.1)
plot(x,dnorm(x),type="l")

# normal density with different mean and variance
dnorm(-3,mean=1,sd=sqrt(2))
dnorm(-1,mean=1,sd=sqrt(2))
dnorm(1,mean=1,sd=sqrt(2))
dnorm(3,mean=1,sd=sqrt(2))
dnorm(5,mean=1,sd=sqrt(2))

# standard normal distribution function
pnorm(-4)
pnorm(-2)
pnorm(0)
pnorm(2)
pnorm(4)

# plot standard normal distribution function
dev.new(height=8,width=8,noRStudioGD=TRUE)
x=seq(-4,4,by=.1)
plot(x,pnorm(x),type="l")

# normal distribution with different mean and variance
pnorm(-3,mean=1,sd=sqrt(2))
pnorm(-1,mean=1,sd=sqrt(2))
pnorm(1,mean=1,sd=sqrt(2))
pnorm(3,mean=1,sd=sqrt(2))
pnorm(5,mean=1,sd=sqrt(2))

# standard normal quantile function
qnorm(.005)
qnorm(.025)
qnorm(.5)
qnorm(.975)
qnorm(.995)

# plot standard normal quantiles
dev.new(height=8,width=8,noRStudioGD=TRUE)
x=seq(-4,4,by=.1)
plot(x,dnorm(x),type="l")
qx=qnorm(.025)
lines(x=rep(qx,2),
      y=c(0,dnorm(qx)))
lines(x=rep(-qx,2),
      y=c(0,dnorm(-qx)))

# normal quantiles with different mean and variance
qnorm(.005,mean=1,sd=sqrt(2))
qnorm(.025,mean=1,sd=sqrt(2))
qnorm(.5,mean=1,sd=sqrt(2))
qnorm(.975,mean=1,sd=sqrt(2))
qnorm(.995,mean=1,sd=sqrt(2))

# simulating data
set.seed(12345)
xvals=rnorm(1000,mean=0,sd=2)
xseq=seq(-5,7,l=100)
hist(xvals,freq=FALSE)
lines(xseq,dnorm(xseq,sd=2))

# testing normality
qqnorm(xvals)
qqline(xvals)
shapiro.test(xvals)


#####   STUDENT'S t DISTRIBUTION   #####

# t density function
dt(0,df=1)
dt(0,df=10)
dt(0,df=100)

# t distribution function
pt(0,df=1)
pt(0,df=10)
pt(0,df=100)

# t quantile function
qt(.975,df=1)
qt(.975,df=10)
qt(.975,df=100)
qt(.995,df=1)
qt(.995,df=10)
qt(.995,df=100)

# one-sample t test example
x=c(15.5, 16.2, 16.1, 15.8, 15.6, 16.0, 15.8, 15.9, 16.2)
mean(x)
sd(x)
var(x)
t.test(x)
t.test(x,mu=16,alternative="less",conf.level=.95)

# two-sample t test example
x=c(70, 82, 78, 74, 94, 82)
y=c(64, 72, 60, 76, 72, 80, 84, 68)
t.test(x,y,alternative="greater",var.equal=TRUE)


#####   COMMON DISTRIBUTIONS   #####

# chi-squared density function
dchisq(1,df=1)
dchisq(1,df=10)
dchisq(1,df=100)

# chi-squared distribution function
pchisq(1,df=1)
pchisq(1,df=10)
pchisq(1,df=100)

# chi-squared quantile function
qchisq(.975,df=1)
qchisq(.975,df=10)
qchisq(.975,df=100)
qchisq(.995,df=1)
qchisq(.995,df=10)
qchisq(.995,df=100)

# F density function
df(1,df1=1,df2=1)
df(1,df1=1,df2=10)
df(1,df1=10,df2=10)

# F distribution function
pf(1,df1=1,df2=1)
pf(1,df1=1,df2=10)
pf(1,df1=10,df2=10)

# F quantile function
qf(.975,df1=1,df2=1)
qf(.975,df1=1,df2=10)
qf(.975,df1=10,df2=10)
qf(.995,df1=1,df2=1)
qf(.995,df1=1,df2=10)
qf(.995,df1=10,df2=10)


#####   BASIC PROGRAMMING   #####

# less than, greater than, equal/not equal to
x=y=10
x<y
x<=y
x>y
x>=y
x==y
x!=y

# | = or,  & = and
x=10
y=11
(x<11 | y<11)
(x<11 & y<11)

# if/else statements (example 1)
x=10
if(x>5){
  x=x/2
  y=2*x
} else {
  x=x*2
  y=x
}
x
y

# if/else statements (example 2)
x=4
if(x>5){
  x=x/2
  y=2*x
} else {
  x=x*2
  y=x
}
x
y

# if/else statements (example 3)
myfun<-function(x){
  if(x>5){
    x=x/2
    y=2*x
  } else {
    x=x*2
    y=x
  }
  list(x=x,y=y)
}
myfun(10)
myfun(5)

# for loop
x=11:15
x
for(j in 1:5){
  x[j]=x[j]+1
}
x

# vectorized version of for loop
x=11:15
x
x=x+1
x

# while statement (example 1)
x=80
iter=0
while(x<100){
  x=x+sqrt(x)/10
  iter=iter+1
}
x
iter

# while statement (example 2)
x=80
iter=0
while(x<100 & iter<20){
  x=x+sqrt(x)/10
  iter=iter+1
}
x
iter

# improper while statement
x=80
iter=0
while(x<100){
  x=x-sqrt(x)/10
  iter=iter+1
}

# infinite while statement
x=80
iter=0
while(x<100){
  x=x-x/10
  iter=iter+1
}
