##########   Linear Mixed-Effects Regression Models
##########   Nathaniel E. Helwig (helwig@umn.edu)
##########   Updated: 04-Jan-2017


#####   DEFINE PATHS AND PACKAGES   #####

# define data path (to load data)
datapath = "~/Desktop/notes/data/"

# load car package (for Anova function)
if(!require(car)){
  install.packages("car")
  library("car")
}

# load lme4 package (for lmer function)
if(!require(lme4)){
  install.packages("lme4")
  library("lme4")
}


#####   ONE-WAY RM ANOVA   #####

# load grocery data
url = "http://ww2.coastal.edu/kingw/statistics/R-tutorials/text/groceries.txt"
groceries = read.table(url, header=TRUE)
groceries

# convert to long data frame
grocery = data.frame(price = as.numeric(unlist(groceries[,2:5])), 
                     item = rep(groceries$subject,4), 
                     store = rep(LETTERS[1:4],each=10))
grocery[1:12,]

# check and reset contrasts
contrasts(grocery$store)
contrasts(grocery$store) <- contr.sum(4)
contrasts(grocery$store)

# define one-way repeated measures ANOVA function
aov1rm <- function(X){
  X = as.matrix(X)
  n = nrow(X)
  a = ncol(X)
  mu = mean(X)
  rhos = rowMeans(X) - mu
  alphas = colMeans(X) - mu
  ssa = n*sum(alphas^2)
  msa = ssa / (a - 1)
  mss = a*sum(rhos^2) / (n - 1)
  ehat = X - ( mu + matrix(rhos,n,a) + matrix(alphas,n,a,byrow=TRUE) )
  sse = sum(ehat^2)
  mse = sse / ( (a-1)*(n-1) )
  Fstat = msa / mse
  pval = 1 - pf(Fstat,a-1,(a-1)*(n-1))
  Cmat = cov(X)
  Jmat = diag(a) - matrix(1/a,a,a)
  Dmat = Jmat%*%Cmat%*%Jmat
  gg = ( sum(diag(Dmat))^2 ) / ( (a-1)*sum(Dmat^2) )
  hf = (n*(a-1)*gg - 2) / ( (a-1)*(n - 1 - (a-1)*gg) )
  pgg = 1 - pf(Fstat,gg*(a-1),gg*(a-1)*(n-1))
  phf = 1 - pf(Fstat,hf*(a-1),hf*(a-1)*(n-1))
  list(mu = mu, alphas = alphas, rhos = rhos,
       Fstat = c(F=Fstat,df1=(a-1),df2=(a-1)*(n-1)),
       pvals = c(pGG=pgg,pHF=phf,p=pval),
       epsilon = c(GG=gg,HF=hf),
       vcomps = c(sigsq.e=mse, sigsq.rho=((mss-mse)/a)) )
}

# fit model using aov with fixed-effects syntax
amod = aov(price ~ store + item, data=grocery)
summary(amod)

# fit model using aov with mixed-effects syntax
amod = aov(price ~ store + Error(item/store), data=grocery)
summary(amod)

# fit model using lmer (using ML approach)
nmod = lmer(price ~ 1 + (1 | item), data=grocery, REML=FALSE)
amod = lmer(price ~ store + (1 | item), data=grocery, REML=FALSE)
anova(amod,nmod)

# fit model using lm and Anova (car)
lmod = lm(as.matrix(groceries[,2:5]) ~ 1)
store = LETTERS[1:4]
almod = Anova(lmod, type="III", idata=data.frame(store=store), idesign=~store)
summary(almod,multivariate=FALSE)$univariate
summary(almod,multivariate=FALSE)$pval.adj

# fit model using our function
amod = aov1rm(groceries[,2:5])
amod$Fstat
amod$pvals
amod$eps


#####   LINEAR MIXED EFFECTS MODEL   #####

# load data and declare types
timss = read.table(paste(datapath,"timss1997.txt",sep=""),header=TRUE,
                   colClasses=c(rep("factor",4),rep("numeric",3)))
head(timss)

# get mean math and hoursTV info by school
grpMmath = with(timss,tapply(math,idschool,mean))
grpMhoursTV = with(timss,tapply(hoursTV,idschool,mean))
dev.new(width=8,height=4,noRStudioGD=TRUE)
par(mfrow=c(1,2))
hist(grpMmath,xlim=c(125,175))
hist(grpMhoursTV,xlim=c(2,4))

# merge school mean scores with timss data.frame
timss = merge(timss,data.frame(idschool=names(grpMmath),
                               grpMmath=as.numeric(grpMmath),
                               grpMhoursTV=as.numeric(grpMhoursTV)))
head(timss)

# define group-centered math and hoursTV
timss = cbind(timss,grpCmath=(timss$math-timss$grpMmath),
              grpChoursTV=(timss$hoursTV-timss$grpMhoursTV))
head(timss)

# note: 
# grpCmath and grpChoursTV are level-1 variables
# grpMmath and grpMhoursTV are level-2 variables

# plot level-1 vs level-2 variables (no linear relation)
dev.new(width=8,height=4,noRStudioGD=TRUE)
par(mfrow=c(1,2))
plot(timss$grpMmath,timss$grpCmath,xlab="grpMmath",ylab="grpCmath")
plot(timss$grpMhoursTV,timss$grpChoursTV,xlab="grpMhoursTV",ylab="grpChoursTV")
# note possible ceiling effect for hoursTV

# random one-way ANOVA (ANOVA II Model)
ramod = lmer(science ~ 1 + (1|idschool), data=timss, REML=FALSE)
ramod
fixef(ramod)
ranef(ramod)$idschool[1:4,]

# add math as fixed effect
rimod = lmer(science ~ 1 + math + (1|idschool), data=timss, REML=FALSE)
print(summary(rimod),correlation=FALSE)

# likelihood-ratio test for math
anova(rimod,ramod)

# add group-centered and group-mean math as fixed effects
ri2mod = lmer(science ~ 1 + grpCmath + grpMmath + (1|idschool), 
              data=timss, REML=FALSE)
print(summary(ri2mod),correlation=FALSE)

# likelihood-ratio test for grpMmath and grpCmath
anova(ri2mod,ramod)

# add grade as fixed effects
ri3mod = lmer(science ~ 1 + grpCmath + grpMmath + grade + (1|idschool), 
              data=timss, REML=FALSE)
print(summary(ri3mod),correlation=FALSE)

# likelihood-ratio test for grade given grpMmath and grpCmath
anova(ri3mod,ri2mod)

# add gender as fixed effects
ri4mod = lmer(science ~ 1 + grpCmath + grpMmath + grade + gender 
              + (1|idschool), data=timss, REML=FALSE)
print(summary(ri4mod),correlation=FALSE)

# likelihood-ratio test for gender given grpMmath, grpCmath, and grade
anova(ri4mod,ri3mod)

# add grpChoursTV and grpMhoursTV as fixed effects
ri5mod = lmer(science ~ 1 + grpCmath + grpMmath + grade + gender
              + grpChoursTV + grpMhoursTV + (1|idschool), 
              data=timss, REML=FALSE)
print(summary(ri5mod),correlation=FALSE)

# likelihood-ratio test for grpChoursTV and grpMhoursTV given others
anova(ri5mod,ri4mod)

# random intercept and random slopes (unconstrained covariance)
risucmod = lmer(science ~ 1 + grpCmath + grpMmath + grade + gender 
                + grpChoursTV + grpMhoursTV + (grpCmath+grpChoursTV|idschool), 
                data=timss, REML=FALSE)
print(summary(risucmod),correlation=FALSE)

# random intercept and random slopes (variance components)
risvcmod = lmer(science ~ 1 + grpCmath + grpMmath + grade + gender 
                + grpChoursTV + grpMhoursTV + (grpCmath+grpChoursTV||idschool), 
                data=timss, REML=FALSE)
print(summary(risvcmod),correlation=FALSE)

# likelihood-ratio test for all covariance terms
anova(risucmod,risvcmod)

# random intercept and random slopes (intercept constrained)
risicmod = lmer(science ~ 1 + grpCmath + grpMmath + grade + gender 
                + grpChoursTV + grpMhoursTV + (1|idschool) 
                + (0+grpCmath+grpChoursTV|idschool), data=timss, REML=FALSE)
print(summary(risicmod),correlation=FALSE)

# likelihood-ratio test for covariance terms
anova(risucmod,risicmod)
anova(risicmod,risvcmod)

# with alpha=0.05 we would choose unconstrained covariance
# with alpha=0.01 we would choose variance components covariance
