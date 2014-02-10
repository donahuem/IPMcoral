
library(lme4)
library(MASS)
library(pscl)
library(glmmADMB)

attach(MCr)
rcrt.glm0 <- glm(recruits~1,family=poisson)
rcrt.glmU <- glm(recruits~utrans,family=poisson)
drop1(rcrt.glmU,test="Chisq")
#Keep utrans p-value=6.257e-15

#But overdispersion remains, even when we consider the covariates, let's try a quasipoisson
rcrt.glmq0 <- glm(recruits~1,family=quasipoisson)
rcrt.glmqU <- glm(recruits~utrans,family=quasipoisson)
drop1(rcrt.glmqU,test="Chisq")
#keep utrans p-value=9.344e-05

#Could also try a negative binomial to deal with overdispersion
rcrt.nb0 <- glm.nb(recruits~1)
rcrt.nbU <- glm.nb(recruits~utrans)
drop1(rcrt.nbU,test="Chisq")
#Keep utrans p-value=1.585e-05

#ZERO-INFALATED MODELS
  
#another problem, however, is that the recruit data is zero-inflated.  That is, it has more zeros that you would 
#expect.  Think of it this way -- some years there is recruitment failure, but in the years where there are recruits,
#there can still be some quads that don't get any.  So, there are two processes:  probability of recruitment failure
#and probability of no reruits given there is recruitment in the area
#Try ZIP and ZINB models

#for zero-inflated models, before the | is the count model; after the | is the binomial presence/absence model
# dist is the distribution for the count data; link is the link function fro the binomial 0/1 part


rcrt.zip0 <- zeroinfl(recruits~1|1,dist="poisson", link="logit")  #this fits one parm for binom 0/1 and one parm for poisson count
rcrt.zip0.U <-zeroinfl(recruits~1|utrans,dist="poisson", link="logit")
AIC(rcrt.zip0,rcrt.zip0.U)
#keep utrans delta=12.9847
#ZERO-INFLATED NEGATIVE BINOMIAL FIT BETTER THAN POISSON

rcrt.zinb0 <- zeroinfl(recruits~1|1,dist="negbin", link="logit")  #this fits one parm for binom 0/1 and one parm for poisson count
rcrt.zinb0.U <-zeroinfl(recruits~1|utrans,dist="negbin", link="logit")
# keep utrans delta=11.985


#let's try adding year as a RE
rcrt.lme0.Y <- glmer(recruits~1+(1|year),family=poisson)
rcrt.lmeU.Y <- glmer(recruits~utrans+(1|year),family=poisson)


rcrt.glmm.nb0<-glmmadmb(recruits~1+(1|year),family="nbinom")
rcrt.glmm.nbU<-glmmadmb(recruits~utrans+(1|year),family="nbinom")

AIC(rcrt.glm0,rcrt.glmU,rcrt.glmq0,rcrt.glmqU,rcrt.nb0,rcrt.nbU,rcrt.zip0,rcrt.zip0.U,rcrt.zinb0,rcrt.zinb0.U,rcrt.lme0.Y,rcrt.lmeU.Y,rcrt.glmm.nb0,rcrt.glmm.nbU)
#rcrt.lmeU.Y <- glmer(recruits~utrans+(1|year),family=poisson) is the best fit by AIC


par(mfrow=c(3,3))
plot(resid(rcrt.glm0)~utrans,ylim=c(-2,8))
plot(resid(rcrt.glmU)~utrans,ylim=c(-2,8))
plot(resid(rcrt.glmq0)~utrans,ylim=c(-2,8))
plot(resid(rcrt.glmqU)~utrans,ylim=c(-2,8))
plot(resid(rcrt.nb0)~utrans,ylim=c(-2,8))
plot(resid(rcrt.nbU)~utrans,ylim=c(-2,8))
plot(resid(rcrt.zip0)~utrans,ylim=c(-2,8))
plot(resid(rcrt.zip0.U)~utrans,ylim=c(-2,8))  #MC: resids do look better 
plot(resid(rcrt.zinb0)~utrans,ylim=c(-2,8))
plot(resid(rcrt.zinb0.U)~utrans,ylim=c(-2,8))
plot(resid(rcrt.lme0.Y)~utrans)
plot(resid(rcrt.lmeU.Y)~utrans)
plot(resid(rcrt.glmm.nb0)~utrans)
plot(resid(rcrt.glmm.nbU)~utrans)

par(mfrow=c(3,4))
xyplot(resid(rcrt.glm0)~utrans,ylim=c(-2,8))
xyplot(resid(rcrt.glmU)~utrans,ylim=c(-2,8))
xyplot(resid(rcrt.glmq0)~utrans,ylim=c(-2,8))
xyplot(resid(rcrt.glmqU)~utrans,ylim=c(-2,8))
xyplot(resid(rcrt.nb0)~utrans,ylim=c(-2,8))
xyplot(resid(rcrt.nbU)~utrans,ylim=c(-2,8))
xyplot(resid(rcrt.zip0)~utrans,ylim=c(-2,8))
xyplot(resid(rcrt.zip0.U)~utrans,ylim=c(-2,8))  #MC: resids do look better 
xyplot(resid(rcrt.zinb0)~utrans,ylim=c(-2,8))
xyplot(resid(rcrt.zinb0.U)~utrans,ylim=c(-2,8))
xyplot(resid(rcrt.lme0.Y)~utrans)
xyplot(resid(rcrt.lmeU.Y)~utrans)
xyplot(resid(rcrt.glmm.nb0)~utrans)
xyplot(resid(rcrt.glmm.nbU)~utrans)

detach(MCr)
#################################################################################
# Size distribution of recruits

attach(MCrs)
rec.size<-lm(sizeNext~1)
rec.sizeU<-lm(sizeNext~utrans)
drop1(rec.sizeU,test="Chisq")
#keep utrans p-value=2.698e-15
detach(MCrs)

par(mfrow=c(3,4))
x<-seq(-2,3,by=0.01)
with(MCrs[MCrs$utrans==levels(MCrs$utrans)[1],],hist(sizeNext,freq=FALSE,ylim=c(0,1),xlim=c(-2,3)))
lines(x,dnorm(x,coefficients(rec.sizeU)[1],summary(rec.sizeU)$sigma))
for (i in 2:10){
with(MCrs[MCrs$utrans==levels(MCrs$utrans)[i],],hist(sizeNext,freq=FALSE,ylim=c(0,1),xlim=c(-2,3)))
lines(x,dnorm(x,coefficients(rec.sizeU)[1]+coefficients(rec.sizeU)[i],summary(rec.sizeU)$sigma))
}
##################################################################################

library(lme4)
library(MASS)
library(pscl)
library(glmmADMB)

attach(MCd)
dens.glm0 <- glm(colnum~1,family=poisson)
dens.glmU <- glm(colnum~utrans,family=poisson)
drop1(dens.glmU,test="Chisq")
#keep utrans p-value=2.2e-16

#But overdispersion remains, even when we consider the covariates, let's try a quasipoisson
dens.glmq0 <- glm(colnum~1,family=quasipoisson)
dens.glmqU <- glm(colnum~utrans,family=quasipoisson)
drop1(dens.glmqU,test="Chisq")
# keep utrans p-value=2.013e-09

#Could also try a negative binomial to deal with overdispersion

dens.nb0 <- glm.nb(colnum~1)
dens.nbU <- glm.nb(colnum~utrans)
drop1(dens.nbU,test="Chisq")
# keep utrans p-value=4.475e-11

dens.lme0.Y <- glmer(colnum~1+(1|year),family=poisson)
dens.lmeU.Y <- glmer(colnum~utrans+(1|year),family=poisson)


dens.glmm.nb0<-glmmadmb(colnum~1+(1|year),family="nbinom")
dens.glmm.nbU<-glmmadmb(colnum~utrans+(1|year),family="nbinom")

AIC(dens.glm0,dens.glmU,dens.glmq0,dens.glmqU,dens.nb0,dens.nbU,dens.lme0.Y,dens.lmeU.Y,dens.glmm.nb0,dens.glmm.nbU)
#best model is dens.lmeU.Y <- glmer(colnum~utrans+(1|year),family=poisson) AIC=689.5554

par(mfrow=c(3,3))
plot(resid(dens.glm0)~utrans,ylim=c(-3,10))
plot(resid(dens.glmU)~utrans,ylim=c(-3,10))
plot(resid(dens.glmq0)~utrans,ylim=c(-3,10))
plot(resid(dens.glmqU)~utrans,ylim=c(-3,10))
plot(resid(dens.nb0)~utrans,ylim=c(-3,10))
plot(resid(dens.nbU)~utrans,ylim=c(-3,10))
plot(resid(dens.lme0.Y)~utrans,ylim=c(-3,10))
plot(resid(dens.lmeU.Y)~utrans,ylim=c(-3,10))
plot(resid(dens.glmm.nb0)~utrans,ylim=c(-3,10))
plot(resid(dens.glmm.nbU)~utrans,ylim=c(-3,10))
detach(MCd)
#################################################################################
# Size distribution of colonies

attach(MCds)
col.size<-lm(sizeNext~1)
col.sizeU<-lm(sizeNext~utrans)
drop1(col.sizeU,test="Chisq")
#keep utrans p-value=2.2e-16
AIC(col.size,col.sizeU)
detach(MCds)

par(mfrow=c(3,3))
x<-seq(-2,10,by=0.01)
with(MCds[MCds$utrans==levels(MCds$utrans)[1],],hist(sizeNext,freq=FALSE,ylim=c(0,1),xlim=c(-2,10)))
lines(x,dnorm(x,coefficients(col.sizeU)[1],summary(col.sizeU)$sigma))
for (i in 2:10){
  with(MCds[MCds$utrans==levels(MCds$utrans)[i],],hist(sizeNext,freq=FALSE,ylim=c(0,1),xlim=c(-2,10)))
  lines(x,dnorm(x,coefficients(col.sizeU)[1]+coefficients(col.sizeU)[i],summary(col.sizeU)$sigma))
}



