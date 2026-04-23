##########   Nonparametric Location Tests: k-Sample
##########   Nathaniel E. Helwig (helwig@umn.edu)
##########   Updated: 04-Jan-2017


#####   DEFINE PATHS AND PACKAGES   #####

# load NSM3 package (for various functions)
if(!require(NSM3)){
  install.packages("NSM3")
  library("NSM3")
}


#####   WILCOXON RANK SUM TEST   #####

# example 4.2:  load data
data(alcohol.intake)
alcohol.intake

# example 4.2:  using R (hard way)
r = rank(c(alcohol.intake$x,alcohol.intake$y))
sum(r[1:12])
sum(r[13:23])

# example 4.2:  using R's wilcox.test function
control = alcohol.intake$x
sst = alcohol.intake$y
wilcox.test(control,sst,alternative="greater")
12*11 - 81 + (11*12/2)

# example 4.1:  estimate delta
d = as.vector(outer(control,sst,"-"))
sort(d)
median(d)

# example 4.1:  CI for delta
wilcox.test(control,sst,alternative="greater",conf.int=TRUE)


#####   KRUSKAL-WALLIS ANOVA   #####

# example:  data
sync = c(23,27,23,22,28,24,18,33,21,15,
         19,25,29,25,19,30,29,24,23,36,
         22,16,30,17,19,26,20,17,21,23)
cond = factor(rep(c("fast","normal","slow"),10))

# example:  using R (hard way)
N = 30
Rj = tapply(rank(sync),cond,sum)
H = (12/(N*(N+1)))*sum(Rj^2)/10 - 3*(N+1)
H
tj = tapply(sync,sync,length)
tj = tj[tj>1]
tj
Hstar = H/(1-sum(tj^3-tj)/(N^3-N))
Hstar
1 - pchisq(Hstar,2)

# example:  using R (easy way)
kruskal.test(sync,cond)

# example:  number of ties
sort(sync)
tj = c(2,3,2,2,4,2,2,2,2)


#####   FRIEDMAN TEST   #####

# example 7.1:  data
rounding.times = matrix(c(5.40, 5.50, 5.55,
                          5.85, 5.70, 5.75,
                          5.20, 5.60, 5.50,
                          5.55, 5.50, 5.40,
                          5.90, 5.85, 5.70,
                          5.45, 5.55, 5.60,
                          5.40, 5.40, 5.35,
                          5.45, 5.50, 5.35,
                          5.25, 5.15, 5.00,
                          5.85, 5.80, 5.70,
                          5.25, 5.20, 5.10,
                          5.65, 5.55, 5.45,
                          5.60, 5.35, 5.45,
                          5.05, 5.00, 4.95,
                          5.50, 5.50, 5.40,
                          5.45, 5.55, 5.50,
                          5.55, 5.55, 5.35,
                          5.45, 5.50, 5.55,
                          5.50, 5.45, 5.25,
                          5.65, 5.60, 5.40,
                          5.70, 5.65, 5.55,
                          6.30, 6.30, 6.25),ncol=3,byrow=TRUE)

# example 7.1:  using R (hard way)
rtrank = t(apply(rounding.times,1,rank))
n = 22
k = 3
vrt = as.vector(rtrank)
tj = tapply(vrt,list(rep(1:n,k),vrt),length)
cval = 0
for(i in 1:n){
  tidx = which(is.na(tj[i,])==FALSE)
  tij = tj[i,tidx]
  if(length(tij)<k){cval=cval+sum(tij^3)-k}
}
top = 12*sum((colSums(rtrank)-n*(k+1)/2)^2)
bot = n*k*(k+1)-(1/(k-1))*cval
Sc = top/bot
Sc
1 - pchisq(Sc,2)

# example 7.1:  using R (easy way)
friedman.test(rounding.times)

# example 7.1:  combine times and ranks (for slides)
cbind(rounding.times[,1],rtrank[,1],
      rounding.times[,2],rtrank[,2],
      rounding.times[,3],rtrank[,3])
colSums(rtrank)

# example 7.1:  pretending there are no ties
n = 22
k = 3
(12/(n*k*(k+1)))*sum(colSums(rtrank)^2)-3*n*(k+1)

# example 7.1:  corrected S-star for ties
(12*((53-44)^2 + (47-44)^2 + (32-44)^2))/(22*3*4-0.5*(6*4))

# example 7.1:  find groups with ties
apply(rtrank,1,function(x) length(unique(x)))

