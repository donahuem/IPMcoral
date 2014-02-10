#Model reproduction from fission using mixed-effects logistic and poisson regression
library(lme4) 

#make a centered size variable
MC$ctrlnarea_t <- MC$lnarea_t-mean(MC$lnarea_t)

#Subset to those colonies that underwent fission
MCf <- MC[!is.na(MC$numfiss),]
MCf$clnareafiss <- MCf$lnareafiss - mean(MCf$lnareafiss)
MCf$clnarea_prefiss <- MCf$lnarea_prefiss - mean(MCf$lnarea_prefiss)

#What are we modeling here?
#1. Probability of fission (Ft+1) as a function of size (At) (color by fate_t1), by location

#2. Probability distribution of NUMBER of fission products GIVEN area of colony fissed (area_prefiss) by location

#3. Size distribution of fission products

#####################################################################################
#(1)  Test for effect of size on fission probability
#What are we modeling?
plot(MC$lnarea_t,MC$fiss,col=MC$fate_t1)
legend("right",legend=levels(MC$fate_t1),col=1:6,pch=1)
xyplot(MC$fiss~MC$lnarea_t|MC$utrans, group=MC$fate_t1,auto.key=list(columns=2),xlab="Log(Area)[t]",ylab="Fission",main=list(label=sp))
histogram(~MC$fiss|MC$utrans,type="count")

attach(MC)

#first, just look at FE model
fiss.AU <- glm(fiss~ctrlnarea_t*utrans,family=binomial)
fiss.A_U <- glm(fiss~ctrlnarea_t+utrans,family=binomial)
AIC(fiss.A_U,fiss.AU)
#fiss.A_U has lower AIC and none of the lnarea*utrans slopes are sig, indicating that there is no site * size interaction

#Next, add RE for quad and year  #since utrans is fixed effect, let year RE be simply the mean effect of year on the intercept (year_t)
fiss.A_U.Q0Y0 <- glmer(fiss~ctrlnarea_t+utrans +(1|uquad/uyear),family=binomial)
#estimated variance for uquad is so small as to be meaningless
fiss.A_U.uY0 <- glmer(fiss~ctrlnarea_t+utrans +(1|uyear),family=binomial)
fiss.A_U.Y0 <- glmer(fiss~ctrlnarea_t+utrans +(1|year_t),family=binomial)
fiss.A_U.TY0 <- glmer(fiss~ctrlnarea_t+utrans +(1|utrans/uyear),family=binomial)
#compare 
AIC(fiss.A_U,fiss.AU,fiss.A_U.Q0Y0,fiss.A_U.Y0)
#simpler is better fiss.A_U.Y0: I don't understand why, but no variance is attributed to uquad
#For PC, fiss.A_U is best

#RE only model
fiss.A.U0Q0Y0 <- glmer(fiss~ctrlnarea_t +(1|utrans/uquad/uyear),family=binomial)

#Refit on natural (not-centered) scale
fiss.A_U.Y0 <- glmer(fiss~lnarea_t+utrans +(1|year_t),family=binomial)

#Plot the fission probability function & data
x <-seq(min(lnarea_t),max(lnarea_t),.1)
ns <- length(unique(utrans))
b <-attr(fiss.A_U.Y0,"fixef")[1]+x*attr(fiss.A_U.Y0,"fixef")[2]
pfs <- exp(b)/(1+exp(b))
cl<-rainbow[ns]
#plot(exp(x+mean(lnarea_t)),pxz,col=1,type="l",ylim=c(0,1))
plot(x,pfs,col=cl[1],type="l",ylim=c(0,1),xlab="Log(Area)[t]", ylab="Prob(Fission)")
#title("P(fission) ~ Logit(Log(Area)+Transect)")
points(lnarea_t[utrans==unique(utrans)[1]],fiss[utrans==unique(utrans)[1]],col=cl[1])
for (i in 1:(ns-1)) {
  b <- attr(fiss.A_U.Y0,"fixef")[1]+attr(fiss.A_U.Y0,"fixef")[i+2]+x*attr(fiss.A_U.Y0,"fixef")[2]
  pfs <- exp(b)/(1+exp(b))
  #lines(exp(x+mean(lnarea_t)),pxz,col=i+1)}
  lines(x,pfs,col=cl[i+1])
  points(lnarea_t[utrans==unique(utrans)[i+1]],fiss[utrans==unique(utrans)[i+1]],col=cl[i+1])
  }
legend("left",legend=unique(utrans),col=cl[1:ns],lty=1,pt.cex=0.3,bty="n")

detach(MC)

####################################################################################
#2. Number of Fission Products Given Area_prefiss

attach(MCf)

#plot the problem\
hist(MCf$numfiss)
histogram(~MCf$numfiss|MCf$utrans,type="count",xlab="Number of Fission Products", ylab="Frequency")
xyplot(MCf$numfiss~MCf$lnarea_prefiss|MCf$utrans,xlab="Log(Area) before Fission", ylab="Number of Fission Products",main=list(label=sp))

#First, I fit with Poisson [not correct -- need zero-truncated]
nfiss.glm0 <- glm(numfiss~1,family=poisson)
nfiss.glmU <- glm(numfiss~utrans,family=poisson)
nfiss.glmA <- glm(numfiss~clnarea_prefiss,family=poisson)
nfiss.glmUA <- glm(numfiss~utrans*clnarea_prefiss,family=poisson)
nfiss.glmU_A <- glm(numfiss~utrans+clnarea_prefiss,family=poisson)
drop1(nfiss.glmUA,test="Chi")
drop1(nfiss.glmU_A,test="Chi")
AIC(nfiss.glm0,nfiss.glmU,nfiss.glmA,nfiss.glmUA,nfiss.glmU_A)
#best model is a function of log_area_prefiss but not site
#  check for overdispersion: phi = 67.8/121 = 0.56 --> underdispersed
#  Zuur et al (p.225) suggest not to correct for underdispersion
#BEST FE Poisson MODEL: nfiss.glmA AIC=344.22, Res.Dev.=67, df=118

#what about year?
#First, look at Fixed Effect
nfiss.glmAY <- glm(numfiss~clnarea_prefiss+year_t,family=poisson)
nfiss.glmA_U_Y <- glm(numfiss~clnarea_prefiss+utrans+year_t,family=poisson)
drop1(nfiss.glmA_U_Y)
#suggests year is not important

nfiss.glmeUA.Y0 <- glmer(numfiss~clnarea_prefiss*utrans+(1|year_t),family=poisson)
nfiss.glmeU_A.Y0 <- glmer(numfiss~clnarea_prefiss+utrans+(1|year_t),family=poisson)
nfiss.glmeA.Y0 <- glmer(numfiss~clnarea_prefiss+(1|year_t),family=poisson)
nfiss.glmeA.Y0 <- glmer(numfiss~clnarea_prefiss+(1|year_t),family=poisson)
anova(nfiss.glmeA.Y0,nfiss.glmeUA.Y0,nfiss.glmeU_A.Y0)
#the variance estimates assoc with year are *tiny* and the deviance is 65.7 -not much better
#CONCLUSION: don't bother with REs for year
#also checked FE of year in ztPoisson - no improvement

#ZERO-TRUNCATED POISSON: Next, I realized that I really need a zero-truncated poisson bc there is no possibility of 0 fission products 
# given the structure of the data
library("VGAM")
nfiss.zt0 <- vglm(numfiss~1,family=pospoisson)
AIC(nfiss.zt0)
nfiss.ztU <- vglm(numfiss~utrans,family=pospoisson)
AIC(nfiss.ztU)
nfiss.ztA <- vglm(numfiss~clnarea_prefiss,family=pospoisson)
AIC(nfiss.ztA)
nfiss.ztUA <- vglm(numfiss~utrans*clnarea_prefiss,family=pospoisson)
AIC(nfiss.ztUA)
nfiss.ztU_A <- vglm(numfiss~utrans+clnarea_prefiss,family=pospoisson)
AIC(nfiss.ztU_A)
#once we use a proper model,nfiss.ztU_A is the best fit (delAIC=6)
#CONCLUSION: main effects of utrans and clnarea_t with zero-truncated

#Also test for effect of density
nfiss.ztU_A_D1 <- vglm(numfiss~utrans+clnarea_prefiss+spcoverD,family=pospoisson)
nfiss.ztU_A_D2 <- vglm(numfiss~utrans+clnarea_prefiss+spcover_t,family=pospoisson)
nfiss.ztU_A_D3 <- vglm(numfiss~utrans+clnarea_prefiss+allcover_t,family=pospoisson)
AIC(nfiss.ztU_A_D1)
#none improve the fit (for PC, D3 improved the fit a bit - come back to this)

#refit on natural (not centered) scale
nfiss.ztU_A <- vglm(numfiss~utrans+lnarea_prefiss,family=pospoisson)
summary(nfiss.ztU_A)

#now, do some model evaluation
par(mfrow=c(1,1))
plot(resid(nfiss.ztU_A,type="pearson")~clnarea_prefiss)
plot(resid(nfiss.ztU_A,type="pearson")~utrans)
plot(resid(nfiss.ztU_A,type="pearson")~year_t)
#plot(resid(nfiss.ztU_A,type="pearson")~predict(nfiss.ztU_A))
#there are a few outliers, but nothing systematic

ns <- length(unique(utrans))
cl<-rainbow(ns)
temp_lap <- rep(seq((min(lnarea_t)-.1),(max(lnarea_t)+.1),.05),9)
temp_ut<-{}
for (i in 1:ns) {
  temp_ut <- c(temp_ut,rep(unique(utrans)[i],length(seq((min(lnarea_t)-.1),(max(lnarea_t)+.1),.05))))
}
temp_u <- unique(utrans)[temp_ut]
MyData<-data.frame(lnarea_prefiss=temp_lap,utrans=temp_u)
G <- predictvglm(nfiss.ztU_A,newdata=MyData,type="link",se=TRUE)#with "link" can get s.e. bands; need to exp(G) + 1 for predictions
plot(lnarea_prefiss[utrans==unique(utrans)[1]],numfiss[utrans==unique(utrans)[1]],col=cl[1],xlab="Log(Area)[t]",ylim=c(0.5,10.5),xlim=c(-1,7),ylab="Number of Fission Products")
lines(MyData$lnarea_prefiss[MyData$utrans==unique(MyData$utrans)[1]],exp(G$fit[MyData$utrans==unique(MyData$utrans)[1]])+1,lty=1,col=cl[1])
#title("Number of Fission Products ~ Zero-Trunc-Poisson(Log(Area) + Transect)")
for (i in 2:ns) {
  points(lnarea_prefiss[utrans==unique(utrans)[i]],numfiss[utrans==unique(utrans)[i]],col=cl[i])
  lines(MyData$lnarea_prefiss[MyData$utrans==unique(MyData$utrans)[i]],exp(G$fit[MyData$utrans==unique(MyData$utrans)[i]])+1,lty=1,col=cl[i])
}
legend("topleft",legend=unique(utrans),col=cl[1:ns],lty=1,pt.cex=0.3,bty="n")

detach(MCf)

############################################################################
#3. Size distribution of fission products
#This could be a function of numfiss or area_prefiss or fissarea, but we will ignore this for now

MCfp <- MC[MC$fisser_t!="",]
attach(MCfp)
hist(lnarea_t)
histogram(~lnarea_t|utrans,xlab="Log(Area) of Fission Products", main=list(label="Log(Area) of Fission Products ~ N(mu_trans, var_trans)"))
#histogram(~lnarea_t|year_t)

#look at plot, it's clear that mean and var will vary by utrans.
#no effect of year
sfiss.glsU <- gls(lnarea_t~utrans)
sfiss.glsU_Y <- gls(lnarea_t~utrans+year_t)
anova(sfiss.glsU_Y)
summary(sfiss.glsU)
plot(resid(sfiss.glsU),utrans)  #variance differs by utrans


vU<- varIdent(form=~1|utrans)
sfiss.glsU.vU <- gls(lnarea_t~utrans,weights=vU)
anova(sfiss.glsU,sfiss.glsU.vU)

detach(MCfp)
