##########   Factorial Analysis of Variance
##########   Nathaniel E. Helwig (helwig@umn.edu)
##########   Updated: 04-Jan-2017


#####   DEFINE PATHS AND PACKAGES   #####

# define data path (to load data)
datapath = "~/Desktop/notes/data/"

# load car package (for Type II sums of squares)
if(!require(car)){
  install.packages("car")
  library("car")
}


#####   BALANCED TWO-WAY ANOVA   #####

# load hypertension data
bp=scan(file=paste(datapath,"hypertension.dat",sep=""))
diet=factor(rep(rep(c("no","yes"),each=6),6))
drug=factor(rep(rep(c("X","Y","Z"),each=12),2))
biof=factor(rep(c("present","absent"),each=36))
hyper=data.frame(bp=bp,diet=diet,drug=drug,biof=biof)

# get sums and sums-of-squares
tapply(bp,list(drug,diet),sum)
sumsq=function(x){sum(x^2)}
tapply(bp,list(drug,diet),sumsq)

# parameter calculations by hand
mu = mean(hyper$bp)
alphas = tapply(hyper$bp,drug,mean) - mu
betas = tapply(hyper$bp,diet,mean) - mu
ab11 = (2136/12) - (4188/24) - (6948/36) + mu
ab12 = (2052/12) - (4188/24) - (6336/36) + mu
ab21 = (2424/12) - (4578/24) - (6948/36) + mu
ab22 = (2154/12) - (4578/24) - (6336/36) + mu
ab31 = (2388/12) - (4518/24) - (6948/36) + mu
ab32 = (2130/12) - (4518/24) - (6336/36) + mu
alphabetas=c(ab11,ab12,ab21,ab22,ab31,ab32)

# fit model with dummy coding for diet
contrasts(hyper$drug) <- contr.sum(3)
contrasts(hyper$drug)
contrasts(hyper$diet) <- contr.treatment(2,base=1)
contrasts(hyper$diet)
mymod=lm(bp~drug*diet,data=hyper)
summary(mymod)

# fit model with effect coding
contrasts(hyper$drug) <- contr.sum(3)
contrasts(hyper$drug)
contrasts(hyper$diet) <- contr.sum(2)
contrasts(hyper$diet)
mymod=lm(bp~drug*diet,data=hyper)
summary(mymod)

# sums-of-squares by hand
sst=2473492 - (13284^2)/72
sse=2473492 - (2136^2 + 2052^2 + 2424^2 + 2154^2 + 2388^2 + 2130^2)/12
ssr=22594 - 12814
ssa=24*((-10)^2 + 6.25^2 + 3.75^2)
ssb=36*((-8.5)^2 + 8.5^2)
ssab=12*((-5)^2 + 5^2 + 2.75^2 + (-2.75)^2 + 2.25^2 + (-2.25)^2)

# mean squares by hand
msr=ssr/5
mse=sse/66
msa=ssa/2
msb=ssb
msab=ssab/2

# F-tests by hand
Fs=msr/mse
Fa=msa/mse
Fb=msb/mse
Fab=msab/mse
1-pf(Fs,5,66)
1-pf(Fa,2,66)
1-pf(Fb,1,66)
1-pf(Fab,2,66)

# anova table in R
anova(mymod)

# all pairwise comparisons for drug:diet
mymod=aov(bp~drug*diet,data=hyper)
TukeyHSD(mymod,"drug:diet")

# sample calculations for first row
(2424-2136)/12 - qtukey(.95,6,66)*sqrt((12814/66)*(2/12))/sqrt(2)
(2424-2136)/12 + qtukey(.95,6,66)*sqrt((12814/66)*(2/12))/sqrt(2)

# pariwise comparisons for drug by hand
(4578-4188)/24    # Y - X
(4518-4188)/24    # Z - X
(4518-4578)/24    # Z - Y
(12814+903)/68    # MSE of additive model
qtukey(.95,3,68)  # critical value from studentized range distribution

# pairwise comparison CI: drug Y - X
c((4578-4188)/24 - qtukey(.95,3,68)*sqrt((2/24)*(12814+903)/68)/sqrt(2),
  (4578-4188)/24 + qtukey(.95,3,68)*sqrt((2/24)*(12814+903)/68)/sqrt(2))

# pairwise comparison CI: drug X vs Z
c((4518-4188)/24 - qtukey(.95,3,68)*sqrt((2/24)*(12814+903)/68)/sqrt(2),
  (4518-4188)/24 + qtukey(.95,3,68)*sqrt((2/24)*(12814+903)/68)/sqrt(2))

# pairwise comparison CI: drug Y vs Z
c((4518-4578)/24 - qtukey(.95,3,68)*sqrt((2/24)*(12814+903)/68)/sqrt(2),
  (4518-4578)/24 + qtukey(.95,3,68)*sqrt((2/24)*(12814+903)/68)/sqrt(2))

# pariwise comparisons for diet by hand
(6336-6948)/36    # yes - no
qtukey(.95,2,68)  # critical value from studentized range distribution
(12814+903)/68    # MSE of additive model
c((6336-6948)/36 - qtukey(.95,2,68)*sqrt((2/36)*(12814+903)/68)/sqrt(2),
  (6336-6948)/36 + qtukey(.95,2,68)*sqrt((2/36)*(12814+903)/68)/sqrt(2))

# pairwise comparison of drug in R
mymod=aov(bp~drug+diet,data=hyper)
TukeyHSD(mymod,"drug")

# pairwise comparison of diet in R
mymod=aov(bp~drug+diet,data=hyper)
TukeyHSD(mymod,"diet")


#####   BALANCED THREE-WAY ANOVA   #####

# load hypertension data
bp=scan(file=paste(datapath,"hypertension.dat",sep=""))
diet=factor(rep(rep(c("no","yes"),each=6),6))
drug=factor(rep(rep(c("X","Y","Z"),each=12),2))
biof=factor(rep(c("present","absent"),each=36))
hyper=data.frame(bp=bp,diet=diet,drug=drug,biof=biof)

# assign effect coding scheme
contrasts(hyper$drug)<-contr.sum(3)
contrasts(hyper$diet)<-contr.sum(2)
contrasts(hyper$biof)<-contr.sum(2)

# fit model 3-way ANOVA with all interactions
mymod=lm(bp~drug*diet*biof,data=hyper)
anova(mymod)

# fit model 3-way ANOVA with all 2-way intactions
mymod=lm(bp~drug*diet+drug*biof+diet*biof,data=hyper)
anova(mymod)

# fit model 3-way ANOVA with only additive effects
mymod=lm(bp~drug+diet+biof,data=hyper)
anova(mymod)

# interaction plot
yhat=tapply(hyper$bp,list(hyper$drug,hyper$diet,hyper$biof),mean)
dev.new(width=12,height=6,noRStudioGD=TRUE)
par(mfrow=c(1,2))
mytitles=c("Biofeedback Absent","Biofeedback Present")
for(k in 1:2){
  plot(1:3,yhat[,1,k],ylim=c(165,215),xlab="Drug",
       ylab="Mean BP",main=mytitles[k],axes=FALSE,type="l")
  lines(1:3,yhat[,2,k],lty=2)
  legend("topleft",c("Diet No","Diet Yes"),lty=1:2,bty="n")
  axis(1,at=1:3,labels=c("X","Y","Z"))
  axis(2)
}

# pariwise comparisons on marginal means
mymod=aov(bp~drug+diet+biof,data=hyper)
TukeyHSD(mymod,"drug")
TukeyHSD(mymod,"diet")
TukeyHSD(mymod,"biof")


#####   UNBALANCED THREE-WAY ANOVA   #####

# load hypertension data
bp=scan(file=paste(datapath,"hypertension.dat",sep=""))
diet=factor(rep(rep(c("no","yes"),each=6),6))
drug=factor(rep(rep(c("X","Y","Z"),each=12),2))
biof=factor(rep(c("present","absent"),each=36))
hyper=data.frame(bp=bp,diet=diet,drug=drug,biof=biof)

# assign effect coding scheme
contrasts(hyper$drug)<-contr.sum(3)
contrasts(hyper$diet)<-contr.sum(2)
contrasts(hyper$biof)<-contr.sum(2)

# fit unbalanced 3-way ANOVA
mymod=lm(bp~drug*diet*biof,data=hyper[1:71,])

# Type I SS
anova(mymod)

# Type II SS
Anova(mymod,type=2)

# Type III SS
Anova(mymod,type=3)


#####   UNBALANCED TWO-WAY ANOVA   #####

# load hypertension data
bp=scan(file=paste(datapath,"hypertension.dat",sep=""))
diet=factor(rep(rep(c("no","yes"),each=6),6))
drug=factor(rep(rep(c("X","Y","Z"),each=12),2))
biof=factor(rep(c("present","absent"),each=36))
hyper=data.frame(bp=bp,diet=diet,drug=drug,biof=biof)

# fit unbalanced 2-way ANOVA
mymod=lm(bp~drug*diet,data=hyper[1:71,])

# check type I and type II
anova(mymod)
Anova(mymod)

# type I for drug
F1mod=lm(bp~drug,data=hyper[1:71,])
R1mod=lm(bp~1,data=hyper[1:71,])
anova(R1mod,F1mod)
anova(F1mod)
anova(mymod)     # note: same Sum Sq for drug

# type II for drug
F2mod=lm(bp~diet+drug,data=hyper[1:71,])
R2mod=lm(bp~diet,data=hyper[1:71,])
anova(R2mod,F2mod)
anova(F2mod)
Anova(mymod)     # note: same Sum Sq for drug

