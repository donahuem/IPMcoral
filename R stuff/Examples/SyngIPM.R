## Syng IPM: harvest of Golden grass scapes in Brazil
## Orou G. Gaoue (ogaoue@bio.miami.edu), March 11-12, 2010

rm(list=ls(all=TRUE)) ## Clear all

############################# PART (1) #######################
# Read in the data MODIFY THE PATH FOR YOUR FILE SYSTEM
dat<-read.table("/Users/ogaoue/Desktop/IPM Syng/dat.syng.txt", header=T); 
attach(dat)

### I - Data manip and exploration #####

str(dat) ## To check the structure of the dataset
## Redefining the class for each variable
for (i in 1:4){
	dat[,i]<-as.factor(dat[,i])
	}
	
for (i in 7:8){
	dat[,i]<-as.factor(dat[,i])
	}	
	
for (i in 11:16){
	dat[,i]<-as.factor(dat[,i])
	}

## Creating a new variables in "dat' for survival in 07
dat$lived07<-as.factor(1*(!is.na(size.06) & !is.na(size.07)))

## We are subseting the data to remove the big ind (>=80 mm diam), which seem to be outiliers for the survival model.
## dat<-subset(dat, size.06<80)
	
attach(dat) ## attach the data
str(dat) 
head(dat, 5)
summary(dat)

## Will be using log transformed size for all analyses
s0<-log(size.06); s1<-log(size.07)

######################### (PART 2) ############################
###  Fiting the statistical models 

### A. Survival (logistic regression for size and harvest)
surv<-glm(lived07~s0*treat, family=binomial, data=dat)

## surv2<-glm(lived07~s0*fam.count06*treat, family=binomial, data=dat)

## Survival analysis with a mixed model, with site as a random effect
library(lme4)
(surv1<-glmer(lived07~s0*treat+(1|Site), family=binomial, data=dat))
summary(surv)

## The 2 models are equally supported (AIC = 2863 and 2863.9). We'll select the simplest glm model. Let's save our regression coefficients for harvest (H) and unharvested (UH) pops: 

							## *********
s.itH<-surv$coef[1]			## intercept for harvest
s.slH<-surv$coef[2]			## slope for harvest

s.itUH<-surv$coef[3]+s.itH	## intercept for unharvest
s.slUH<-surv$coef[4]+s.slH	## slope for unharvest


### B. GROWTH
## How does the data look?
op<-par(mfrow=c(1,2))
plot(s0, s1, col=as.numeric(treat), xlab="Diameter at t (mm)", ylab="Diameter at t+1 (mm)")
plot(size.06, size.07, col=as.numeric(treat), ylim=c(0,100))
par(op)

## Running linear model of sizes with or without exponential variance of size, and with or without variance partitionned by treatment (harvest)
library(nlme)
g1<-gls(s1~s0*treat, na.action=na.omit, data=dat, method="ML")

g2p<-gls(s1~s0*treat, na.action=na.omit, weight=varPower(form=~fitted(.)|treat), data=dat, method="ML")

g2<-gls(s1~s0*treat, na.action=na.omit, weight=varExp(form=~fitted(.)), data=dat, method="ML")

anova(g2, g2p)
## Modeling variance as a power, rather than exponential function of size did not improve our model. so we will stick with the varExp. And now simplify the variance model.

g3<-gls(s1~s0*treat, na.action=na.omit, weight=varExp(form=~fitted(.)|treat), data=dat, method="ML")

## Comparing the 3 models to find the best and simplest!!!
anova(g1, g2, g3) ## Model g2 is one of the best and the simplest

g4<-gls(s1~s0:treat, na.action=na.omit, weight=varExp(form=~fitted(.)), data=dat, method="ML")

g5<-gls(s1~s0+treat, na.action=na.omit, weight=varExp(form=~fitted(.)), data=dat, method="ML")

g6<-gls(s1~s0, na.action=na.omit, weight=varExp(form=~fitted(.)), data=dat, method="ML")

anova(g2, g4, g5, g6)
## The simplest model g6, harvest independent growth is as equally supported as the other models and we will consider that. We update that to use a REML method instead of a ML for our true coefficients.

g7<-update(g6, method="REML")
summary(g7)

## Let save the coefficients        ## *********
g.it<-summary(g7)$coef[1]   		## common intercept
g.sl<-summary(g7)$coef[2]			## common slope
g.sig2<-summary(g7)$sigma^2 		## overall variance
varExp<-g7$modelStruct$varStruct	## size-dependent variance


## C. FERTILITY - SEED ORIGINATED SEEDLINGS

## C.1. Number of fruits, modeled as a Poisson (?)
plot(size.06, esc.06, log="x")

## The plot shows that reproducing ind are the one with size >= 20 mm
## We chop off the individual that are not reproducing (size < 20), removed individuals that has NA in 2006 (new born), and replaced the NAs for adults by 0 because they did not produce scape, not that they can't produce any.
size0620<-size.06[size.06>=20 & !is.na(size.06)]
esc20<-esc.06[size.06>=20 &!is.na(size.06)]
esc20[is.na(esc20)]<-0
treat20<-treat[size.06>=20 & !is.na(size.06)]

## And run a Poisson model of size dependent scapes production
scape<-glm(esc20~log(size0620)*treat20, poisson, data=dat)
summary(scape)

## Looking at the ratio of the deviance/df, 3432/2036 = 1.685658 >> 1, we estimated that there is an overdispersion and opted to use a negative binomial distribution instead.

library(MASS)
scape.nb<-glm.nb(esc20~log(size0620)*treat20, data=dat)
scape.nb2<-glm.nb(esc20~log(size0620), data=dat)
anova(scape.nb, scape.nb2)

## We fit the full model including treatment (harvest), and realized that there is not effect of harvest on scape number. We removed treatment from the model and tested if this may make the model worst. It did not! We keep the simplest model 'scape.nb2'.

summary(scape.nb2)

## Let's save the regression coefficient for the neg. binomial
f.int<-summary(scape.nb2)$coef[1]   ## *********
f.sl<- summary(scape.nb2)$coef[2]

## C.2 Estimating the scape/seedling ratio "p.est". 
## We used the total number of scape produced in 2006 (nscape06), and the number of seedling produced in 2007 (nsdl07) as the length of the vector of individuals that were NA in 2006 (not born) and NOT NA in 2007 (just born) and we also considered seedlings that originated from scape (ori.07=0) rather than clonal.

nscape06<-sum(esc.06, na.rm=T)
nsdl07<-length(size.07[is.na(size.06) & !is.na(size.07) & ori.07=="0"])
p.est<-nsdl07/nscape06

## The value for nsdl07 (199) is different from the raw value (168). We need to figure out what is wrong but in the meantime, we will use the 168 value to get the p.est calculated
nsdl.sx<-168
p.est.cal<-nsdl.sx/nscape06   ## *********

## C.3. Size distribution of seedlings produced in 2007
## Checking the distribution of size (histogram), raw data. We will log tranform the sizes of the seedling to be consistent with other functions
hist(size.07[ori.07=="0"])

## We calculated the mean and standard deviation of the sizes of seedling that originated from seeds to build the normal dist of size
mu.sdl<-mean(log(size.07[ori.07=="0"]), na.rm=T) ## *********
sig.sdl<-sqrt(var(log(size.07[ori.07=="0"]), na.rm=T))


## D. FERTILITY - CLONAL 
## D.1. Probability of producing clonal 
clonality<-cl.rep.06.07[!is.na(cl.rep.06.07)]
size.clon<-size.06[!is.na(cl.rep.06.07)]
treat.clon<-treat[!is.na(cl.rep.06.07)]

clon<-glm(clonality~log(size.clon)*treat.clon, binomial, data=dat)
clon1<-glm(clonality~log(size.clon), binomial, data=dat)
anova(clon, clon1, test="Chi")

## There is no size dependent probability of producing clonals. We select model "clon1" that is the simplest with equal surpport as "clon"

## Let's save the regression coefficients
c.int<-summary(clon1)$coef[1]  ## *********
c.sl<-summary(clon1)$coef[2]

## D.2. Number of clonals seedlings producted by each adults. It is going to be non size dependent.

## Number of seedling clonal originated
nsdl07cl<-length(size.07[is.na(size.06) & !is.na(size.07) & ori.07=="1"])
## This formula (gives 439) is not giving what we should get (= 407)
nsdl07cl<-407

## The total number of adults that produced clonals should be 
totclon<-length(clonality[clonality=="1"])
## but there is a problem with the data and the true number is 339
totclon<-339

## The mean number of clonal produced per adults
mu.cl07<-nsdl07cl/totclon   ## *********

## D.3. Size distribution of the clonals
## We will calculate the mean and standard deviation of the the size of the clonals to build the normal distribution

mu.sdlcl<-mean(log(size.07[ori.07=="1"]), na.rm=T) ## *********
sig.sdlcl<-sqrt(var(log(size.07[ori.07=="1"]), na.rm=T))


########################## PART (3) ###########################
### Matrix of all the regression coefficients
ntreat<-length(unique(treat))
ncoef<-16
p.vec<-matrix(0, ncoef, ntreat)

p.vec[1,]<-c(s.itH, s.itUH)
p.vec[2,]<-c(s.slH, s.slUH)
p.vec[3,]<-rep(g.it,ntreat)
p.vec[4,]<-rep(g.sl,ntreat)
p.vec[5,]<-rep(g.sig2,ntreat)
p.vec[6,]<-rep(varExp,ntreat)
p.vec[7,]<-rep(f.int,ntreat)
p.vec[8,]<-rep(f.sl,ntreat)
p.vec[9,]<-rep(p.est.cal,ntreat)
p.vec[10,]<-rep(mu.sdl,ntreat)
p.vec[11,]<-rep(sig.sdl,ntreat)
p.vec[12,]<-rep(c.int,ntreat)
p.vec[13,]<-rep(c.sl,ntreat)
p.vec[14,]<-rep(mu.cl07,ntreat)
p.vec[15,]<-rep(mu.sdlcl,ntreat)
p.vec[16,]<-rep(sig.sdlcl,ntreat)

dimnames(p.vec)<-list(c("int.survival", "slope.survival", "int.growth", "slope.growth", "var.growth", "varExp", "int.nscape", "slope.nscape", "p.establish", "mean.sdl.size", "sigma.sdl.size", "int.clonality", "slope.clonality", "mean.nclonal", "mean.ram.size", "sigma.ram.size"), c("Harvest", "Unharvest"))

############################# PART (4)#########################
## FIGURES

op<-par(mfrow=c(1,3))
## PLOT: Survival
os<-order(s0) # ordering size 
os.size<-(s0)[os] # reorganizing size based on that order
os.surv<-with(dat, (lived07)[os]) # reoganizing surv based on order os

## We chop off the last 605 obs of os.size and os.surv to remove the new born that were not included in the survival column.
os.size<-os.size[!is.na(os.size)]
os.surv<-os.surv[1:length(os.size)]

psz<-tapply(os.size,as.numeric(cut(os.size,100)),mean,na.rm=T)
ps<-tapply((as.numeric(os.surv)-1),as.numeric(cut(os.size,100)),mean,na.rm=T)
sizef<-seq(min(s0,na.rm=T),max(s0,na.rm=T),length=100)

plot(psz,ps,pch=19, main="Survival", xlab=expression("Diameter at t (log scale)"), ylab=expression("Probability of surviving"), cex.lab=1.5)

points(sizef, exp(s.itH + s.slH*sizef)/(1 + exp(s.itH +s.slH*sizef)),type="l", col="red", lty=3)

points(sizef,exp(s.itUH + s.slUH*sizef)/(1 + exp(s.itUH +s.slUH*sizef)),type="l", col="blue", lty=3)


## PLOT for Growth
## The issue here is how to fit the regression line when the plot is on log-log scale. I presented here two approaches (use log10 or ggplot2)

## 1st method:log10
## quartz() ## use x11() for PC!!!
plot(size.06, size.07, log="xy", main="Growth", col=as.numeric(treat), xlab="Diameter at t (mm, log scale)", ylab="Diameter at t+1 (mm, log scale)", cex.lab=1.5)
abline(lm(log10(size.07)~log10(size.06)), lty=2, lwd=2, col="blue")

## 2d method: ggplot2
## install.packages("ggplot2", dependencies=T)
## library(ggplot2) 
## qplot(size.06, size.07, log="xy") + geom_smooth(method = lm)


## PLOT Fecundity
## Reorganizing, re-ordering (mainly the fitted values) the data to get a graph that makes sense!!
os20<-order(size0620)
size0620n<-size0620[os20]
fit20<-scape.nb2$fitted.values
fit20<-fit20[os20]

## There is the plot
plot(size.06,esc.06, log="x", main="Fecundity", xlab="Diameter at t (mm, log scale)", ylab="Number of scapes at t", cex.lab=1.5)
points(size0620n,fit20, type="l", lwd=3, lty=3, col="red");

par(op)


############################# PART (5) ##########################
## The IPM FUNCTIONS: s(x), g(y,x), f(y,x), c(y,x), p(y,x), K(y,x)

### A. SURVIVAL function s(x)
sx<-function(x, pvec, treat) {
	xbeta<-pvec[1,treat] + pvec[2, treat]*x;
	s<-exp(xbeta)/(1+exp(xbeta))
	return(s);
}

### B. GROWTH function g(y,x)
gyx<-function(y, x, pvec, treat) {
	mux<-pvec[3, treat] + pvec[4, treat]*x
	sigmax2<-pvec[5, treat]*exp(2*pvec[6, treat]*mux)
	sigmax<-sqrt(sigmax2)
	g<-dnorm(y, mux, sigmax)
	return(g)
}

### C.The SURVIVAL-GROWTH function P(y, x)
pyx<-function(y,x, pvec, treat) { 
	p<-sx(x, pvec, treat)*gyx(y, x, pvec, treat)
	return(p) 
	}

### D. FERTILITY (SEEDS) function f(x,y)

fyx<-function(y, x, pvec, treat) {
	nscape<-exp(pvec[7, treat] + pvec[8, treat]*x); ## log(mu)=a + bx
	scd<-dnorm(y, pvec[10, treat], pvec[11, treat])
	f<-nscape*pvec[9, treat]*scd
	return(f)
}

### D. CLONAL FERTILITY
cyx<-function(y, x, pvec, treat) {
	xbeta<-pvec[12, treat] + pvec[13, treat]*x;
	p_clon<-exp(xbeta)/(1+exp(xbeta))
	nclone<-pvec[14, treat]
	scd_c<-dnorm(y, pvec[15, treat], pvec[16, treat])
	c<-p_clon*nclone*scd_c
	return(c)
}

### The (master) KERNEL: K(y,x)= p(y,x) + f(y,x) + c(y,x)
Kyx<-function(y, x, pvec, treat) {
	k<-pyx(y, x, pvec, treat)+fyx(y, x, pvec, treat)+cyx(y, x, pvec, treat)
	return(k) 
}

######################### PART (6) #############################
### NUMERICAL INTEGRATION 

## The big matrix function

bigmat<-function(bigM, pvec, treat){
	## Set matrix size and convergence tolerance 
	min.sz<-0.9*min((range(s0,na.rm=T)))
	max.sz<-1.1*max((range(s0,na.rm=T)))

	# Compute meshpoints iteration matrix KD 
	h=(max.sz-min.sz)/(bigM+1); 
	y=(seq(min.sz,max.sz,length=bigM)+seq(min.sz+h,max.sz+h,length=bigM))/2;  

	## Apply Kyx funct for y and y, pvec (=p.vec), treat (choose!!)
	K=outer(y,y,Kyx, pvec, treat);
	KD=h*K;
	return(KD);  		  ## This is your have the matrix
}

############################# PART (7) ##########################
## Eigenvalues and eigenvectors analysis, sensitivity,...

library(popbio)
library(Biodem)

## Checks: Ploting lambda for various big matrix size to see when it level off:

bigm<-c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
lam<-matrix(0, length(bigm), 2)
ntreat<-2

for (i in 1:length(bigm)){
	for (j in 1:ntreat){
	library(popbio)
	A<-array(0, c(bigm[i], bigm[i], length(bigm), ntreat))
	A[,,i,j]<-bigmat(bigm[i], pvec=p.vec, treat=j)
	lam[i,j]<-eigen.analysis(A[,,i,j])$lambda1
	}
}	

## PLOT FOR THE BIG MATRIX
op<-par(mfrow=c(1,2))

plot(bigm, lam[,1], main="Harvest", type="l", xlab="", ylab="", col="red", lwd=3)

plot(bigm, lam[,2], main="Unharvest", type="l", xlab="", ylab="", col="blue", lty=2, lwd=3)

par(oma=c(5, 4.5, 0, 0))
title(xlab="Size of the big matrix", ylab=expression(paste("Population growth rate, ", lambda)), outer=T, cex.lab=1.2)

par(op)

## The big matrix is 300 x 300, more than that, lambda does not change much

## Lambda, SSD, RV, Damping ratio, Elasticity matrix
big<-array(0, c(300, 300, ntreat))
eig<-vector("numeric", ntreat)

for (i in 1:ntreat){
	big[,,i]<-bigmat(300, pvec=p.vec, treat=i)
	eig[i]<-eigen.analysis(big[,,i])$lambda1
}

## To get lambda for each treatment
eig

############################# END OF #############################
