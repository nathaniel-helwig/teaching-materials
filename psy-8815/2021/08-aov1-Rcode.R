##########   One-Way Analysis of Variance
##########   Nathaniel E. Helwig (helwig@umn.edu)
##########   Updated: 04-Jan-2017


#####   CATEGORICAL PREDICTORS   #####

# dummy coding
x=factor(rep(letters[1:3],each=5))
x
contrasts(x)<-contr.treatment(nlevels(x))
contrasts(x)

# effect coding
x=factor(rep(letters[1:3],each=5))
x
contrasts(x)<-contr.sum(nlevels(x))
contrasts(x)


#####   ONE-WAY ANOVA MODEL   #####

# simulate data from one-way ANOVA
x=factor(sample(letters[1:3],30,replace=TRUE))
mus=c(-1,0,2)
y=mus[as.integer(x)]+rnorm(30)

# fit model with dummy coding (default)
mymod=lm(y~x)
summary(mymod)

# fit model with effect coding
contrasts(x)<-contr.sum(nlevels(x))
mymod=lm(y~x)
summary(mymod)

# define memory data
sync=c(23,27,23,22,28,24,18,33,21,15,
       19,25,29,25,19,30,29,24,23,36,
       22,16,30,17,19,26,20,17,21,23)
cond=factor(rep(c("fast","normal","slow"),10))

# get sums and sums-of-squares
tapply(sync,cond,sum)
sumsq=function(x){sum(x^2)}
tapply(sync,cond,sumsq)

# fit one-way ANOVA model
smod=lm(sync~cond)
summary(smod)$coef

# refit using normal group as baseline (reference)
contrasts(cond)
contrasts(cond)<-contr.treatment(3,base=2)
contrasts(cond)
smod=lm(sync~cond)
summary(smod)$coef

# refit using effect coding
contrasts(cond)<-contr.sum(nlevels(cond))
contrasts(cond)
smod=lm(sync~cond)
summary(smod)$coef

# SS, MS, F, and p-value by hand
sst=(4738+7742+4810) - (704^2)/30
sse=(4738+7742+4810) - ((212^2+274^2+218^2)/10)
ssr=sst-sse
msr=ssr/2
mse=sse/27
Fstat=msr/mse
1-pf(Fstat,2,27)

# ANOVA table in R
anova(smod)


#####   MULTIPLE COMPARISONS   #####

# define pairwise comparisons
L1 = 21.2 - 27.4   # fast - normal
L2 = 27.4 - 21.8   # normal - slow
L3 = 21.8 - 21.2   # slow - fast

# define variance of each pairwise comparison
VL=mse*2/10

# compare different possible critical values
qt(.975,27)
qt(1-(.05/3)/2,27)
qtukey(.95,3,27)/sqrt(2)
sqrt(2*qf(.95,2,27))

# pairwise comparison CIs using no correction
c(L1-qt(.975,27)*sqrt(VL),L1+qt(.975,27)*sqrt(VL))
c(L2-qt(.975,27)*sqrt(VL),L2+qt(.975,27)*sqrt(VL))
c(L3-qt(.975,27)*sqrt(VL),L3+qt(.975,27)*sqrt(VL))

# pairwise comparison CIs using Bonferroni correction
c(L1-qt(1-(.05/3)/2,27)*sqrt(VL),L1+qt(1-(.05/3)/2,27)*sqrt(VL))
c(L2-qt(1-(.05/3)/2,27)*sqrt(VL),L2+qt(1-(.05/3)/2,27)*sqrt(VL))
c(L3-qt(1-(.05/3)/2,27)*sqrt(VL),L3+qt(1-(.05/3)/2,27)*sqrt(VL))

# pairwise comparison CIs using Tukey correction
c(L1-(qtukey(.95,3,27)/sqrt(2))*sqrt(VL),L1+(qtukey(.95,3,27)/sqrt(2))*sqrt(VL))
c(L2-(qtukey(.95,3,27)/sqrt(2))*sqrt(VL),L2+(qtukey(.95,3,27)/sqrt(2))*sqrt(VL))
c(L3-(qtukey(.95,3,27)/sqrt(2))*sqrt(VL),L3+(qtukey(.95,3,27)/sqrt(2))*sqrt(VL))

# pairwise comparison CIs using Scheffe correction
c(L1-sqrt(2*qf(.95,2,27))*sqrt(VL),L1+sqrt(2*qf(.95,2,27))*sqrt(VL))
c(L2-sqrt(2*qf(.95,2,27))*sqrt(VL),L2+sqrt(2*qf(.95,2,27))*sqrt(VL))
c(L3-sqrt(2*qf(.95,2,27))*sqrt(VL),L3+sqrt(2*qf(.95,2,27))*sqrt(VL))

# adding another contrast:  normal - (slow+fast)/2
qt(1-(.05/4)/2,27)     # Bonferroni critical value
sqrt(2*qf(.95,2,27))   # Scheffe critical value
mse*(1+0.25+0.25)/10   # variance of contrast

# Bonferroni CI for new contrast
c(5.9-qt(1-(.05/4)/2,27)*sqrt(mse*(1+0.25+0.25)/10),
  5.9+qt(1-(.05/4)/2,27)*sqrt(mse*(1+0.25+0.25)/10))

# Scheffe CI for new contrast
c(5.9-sqrt(2*qf(.95,2,27))*sqrt(mse*(1+0.25+0.25)/10),
  5.9+sqrt(2*qf(.95,2,27))*sqrt(mse*(1+0.25+0.25)/10))

