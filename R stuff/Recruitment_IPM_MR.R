library(lme4)
attach(MCr)
hist(sizeNext)
histogram(~sizeNext|utrans,type="count")

MCr$lnspcover <- log(MCr$spcover_t+.1)
MCr$lnspcover2<- MCr$lnspcover^2
MCr$pspcover <- MCr$spcover_t/(0.33*100^2)  #percent cover based on 0.33m2 
MCr$pspcover2 <- MCr$pspcover^2           

#MCr$logitspcover <- log(MCr$pspcover/(1-MCr$pspcover))   #too many zeros for this to work
detach(MCr)
attach(MCr)
rcrt.glm0 <- glm(Freq~1,family=poisson)
#AIC=748, note overdispersion:  Res Deviance/df = 527/219 = 2.4. 
rcrt.glmU <- glm(Freq~utrans,family=poisson)
#AIC=677; lambda sig varies by utrans; overdisp:  440/219 = 2.09
rcrt.glmU_D <- glm(Freq~utrans + lnspcover,family=poisson)
#AIC=616; this is a bit better (del AIC is 7); overdispersion: 376/210 = 1.79
rcrt.glmUD <- glm(Freq~utrans*lnspcover,family=poisson)
#AIC=618; overdispersion:  363/202=1.99
rcrt.glmU_D2 <- glm(Freq~utrans + lnspcover + lnspcover2,family=poisson)
#AIC=589.15; overdispersion: 348/209=1.67
drop1(rcrt.glmUD,test="Chi")
drop1(rcrt.glmU_D,test="Chi")
drop1(rcrt.glmU_D2,test="Chi")
#using percent cover doesn't fit as well as lnspcover
rcrt.glmUP <- glm(Freq~utrans*pspcover,family=poisson)
#AIC=655, phi=401/202=1.99
rcrt.glmU_P <- glm(Freq~utrans+pspcover,family=poisson)
#AIC=670
rcrt.glmU_P2 <- glm(Freq~utrans+pspcover+pspcover2,family=poisson)
#for MC:AIC=654; phi= 413/209=1.92; 
#for PC: AIC=678; phi=438/212
#for PL: AIC=513, phi=318/234=1.35


#let's try adding year as a RE
rcrt.lme0.Y <- glmer(recruits~1+(1|year),family=poisson)
#AIC=512.8,ResDev=508
rcrt.lmeU.Y <- glmer(recruits~utrans+(1|year),family=poisson)
#AIC=436, ResDev=416  (PC: AIC = 638; PL: 415.3)
#can't include both spcover and year%in%utrans in the model bc not unique
#  looking at the ResDev (348-376 including lnspcover and no RE; 416 including only year),
#  it seems like including lnspcover is good
plot(resid(rcrt.lmeU.Y)~utrans)
plot(resid(rcrt.glmU_D2)~utrans)
xyplot(resid(rcrt.glmU_D2)~lnspcover|utrans)
#That said, the residuals look nicer for rcrt.U_D2 than rcrt.U.Y

#But overdispersion remains, even when we consider the covariates, let's try a quasipoisson
rcrt.glmq0 <- glm(recruits~1,family=quasipoisson)
rcrt.glmqU <- glm(recruits~utrans,family=quasipoisson)
  #phi-2.69
rcrt.glmqU_D <- glm(recruits~utrans + lnspcover,family=quasipoisson)
#phi=2.28
rcrt.glmqU_D2 <- glm(recruits~utrans + lnspcover + lnspcover2,family=quasipoisson)
#phi=1.66
rcrt.glmqUD <- glm(recruits~utrans*lnspcover,family=quasipoisson)
#phi=1.75
rcrt.glmqUD2 <-glm(recruits~utrans+lnspcover+lnspcover2+lnspcover2:utrans+lnspcover:utrans,family=quasipoisson)
#phi=1.575; ResDev=304.4 
drop1(rcrt.glmqUD,test="Chi")
drop1(rcrt.glmqU_D,test="Chi")
drop1(rcrt.glmqU_D2,test="Chi")
drop1(rcrt.glmqUD2,test="Chi")
#model comparisons are similar, but this scales the parms and p-values accordingly.
#PC: best is U_D2
#PL: best is UD2

#Could also try a negative binomial to deal with overdispersion
library(MASS)
rcrt.nb0 <- glm.nb(recruits~1)
rcrt.nbU <- glm.nb(recruits~utrans)
rcrt.nbU_D <- glm.nb(recruits~utrans+lnspcover)
rcrt.nbU_D2 <- glm.nb(recruits~utrans+lnspcover+lnspcover2)
rcrt.nbUD <- glm.nb(recruits~utrans*lnspcover)
rcrt.nbUD2i <- glm.nb(recruits~utrans+lnspcover+lnspcover2+lnspcover:utrans)
rcrt.nbUD2ii <- glm.nb(recruits~utrans+lnspcover+lnspcover2+lnspcover2:utrans)

anova(rcrt.nb0,rcrt.nbU,rcrt.nbU_D,rcrt.nbU_D2,rcrt.nbUD,rcrt.nbUD2i,rcrt.nbUD2ii)
anova(rcrt.nbU_D2,rcrt.nbUD2i)
#MC: best model is rcrt.nbU_D2; ResDev=165.8; AIC 518.84
#PC: best model is rcrt.nbU_D2, but marginal compared to ncUD2i; for ease, same as MC
#Neg Binomal model has lower residual deviance than Poisson models
summary(rcrt.nbU_D2)  #for all spp

#ERO-INFALATED MODELS
  library(pscl)
#another problem, however, is that the recruit data is zero-inflated.  That is, it has more zeros that you would 
#expect.  Think of it this way -- some years there is recruitment failure, but in the years where there are recruits,
#there can still be some quads that don't get any.  So, there are two processes:  probability of recruitment failure
#and probability of no reruits given there is recruitment in the area
#Try ZIP and ZINB models

#for zero-inflated models, before the | is the count model; after the | is the binomial presence/absence model
# dist is the distribution for the count data; link is the link function fro the binomial 0/1 part
rcrt.zip0 <- zeroinfl(recruits~1|1,dist="poisson", link="logit")  #this fits one parm for binom 0/1 and one parm for poisson count
rcrt.zip0.U <-zeroinfl(recruits~1|utrans,dist="poisson", link="logit")
rcrt.zipD.U <-zeroinfl(recruits~lnspcover|utrans,dist="poisson", link="logit")
rcrt.zipD2.U <-zeroinfl(recruits~lnspcover+lnspcover2|utrans,dist="poisson", link="logit")
rcrt.zipU_D.0 <-zeroinfl(recruits~utrans+lnspcover|1,dist="poisson", link="logit")
rcrt.zipU_D2.0 <-zeroinfl(recruits~utrans+lnspcover+lnspcover2|1,dist="poisson", link="logit")
   #AIC=527.5
rcrt.zipUD.0 <-zeroinfl(recruits~utrans+lnspcover+utrans:lnspcover|1,dist="poisson", link="logit")
rcrt.zipUD2.0 <-zeroinfl(recruits~utrans+lnspcover+lnspcover2+utrans:lnspcover|1,dist="poisson", link="logit")
rcrt.zipU_D.U <-zeroinfl(recruits~utrans+lnspcover|utrans,dist="poisson", link="logit")
rcrt.zipU_D2.U <-zeroinfl(recruits~utrans+lnspcover+lnspcover2|utrans,dist="poisson", link="logit")
   #AIC=524.8
rcrt.zipUD.U <-zeroinfl(recruits~utrans+lnspcover+utrans:lnspcover|utrans,dist="poisson", link="logit")
rcrt.zipUD2.U <-zeroinfl(recruits~utrans+lnspcover+lnspcover2+utrans:lnspcover|utrans,dist="poisson", link="logit")
AIC(rcrt.zip0,rcrt.zip0.U,rcrt.zipD.U,rcrt.zipD2.U,rcrt.zipU_D.0,rcrt.zipU_D2.0,rcrt.zipU_D.U,rcrt.zipU_D2.U,rcrt.zipUD.0,rcrt.zipUD2.0,rcrt.zipUD.U,rcrt.zipUD2.U)
summary(rcrt.zipU_D2.0)
summary(rcrt.zipU_D2.U)
plot(resid(rcrt.zipU_D2.0)~utrans)
plot(resid(rcrt.zipU_D2.U)~utrans)  #MC: resids do look better 
#for PC, rcrt.zipUD2.U is best AIC=521 (delAIC=3.5)
#for PL, zipU_D2 best
          
#ZERO-INFLATED NEGATIVE BINOMIAL FIT BETTER THAN POISSON
rcrt.zinb0 <- zeroinfl(recruits~1|1,dist="negbin", link="logit")  #this fits one parm for binom 0/1 and one parm for poisson count
rcrt.zinb0.U <-zeroinfl(recruits~1|utrans,dist="negbin", link="logit")
rcrt.zinbD.U <-zeroinfl(recruits~lnspcover|utrans,dist="negbin", link="logit")
rcrt.zinbD2.U <-zeroinfl(recruits~lnspcover+lnspcover2|utrans,dist="negbin", link="logit")
rcrt.zinbU_D.0 <-zeroinfl(recruits~utrans+lnspcover|1,dist="negbin", link="logit")
rcrt.zinbU_D2.0 <-zeroinfl(recruits~utrans+lnspcover+lnspcover2|1,dist="negbin", link="logit")
  #AIC=517.84
rcrt.zinbUD.0 <-zeroinfl(recruits~utrans+lnspcover+utrans:lnspcover|1,dist="negbin", link="logit")
rcrt.zinbUD2.0 <-zeroinfl(recruits~utrans+lnspcover+lnspcover2+utrans:lnspcover|1,dist="negbin", link="logit")
rcrt.zinbU_D.U <-zeroinfl(recruits~utrans+lnspcover|utrans,dist="negbin", link="logit")
rcrt.zinbU_D2.U <-zeroinfl(recruits~utrans+lnspcover+lnspcover2|utrans,dist="negbin", link="logit")
  #AIC=518.3
rcrt.zinbUD.U <-zeroinfl(recruits~utrans+lnspcover+utrans:lnspcover|utrans,dist="negbin", link="logit")
rcrt.zinbUD2.U <-zeroinfl(recruits~utrans+lnspcover+lnspcover2+utrans:lnspcover|utrans,dist="negbin", link="logit")
AIC(rcrt.zinb0,rcrt.zinb0.U,rcrt.zinbD.U,rcrt.zinbD2.U,rcrt.zinbU_D.0,rcrt.zinbU_D2.0,rcrt.zinbU_D.U,rcrt.zinbU_D2.U,rcrt.zinbUD.0,rcrt.zinbUD2.0,rcrt.zinbUD.U,rcrt.zinbUD2.U)
lrtest(rcrt.zinbU_D2.0,rcrt.zinbU_D2.U)
summary(rcrt.zinbU_D2.0)
summary(rcrt.zinbU_D2.U)
plot(resid(rcrt.zinbU_D2.0)~utrans)   #resids do look better
plot(resid(rcrt.zinbU_D2.U)~utrans) 
#For MC and PC and PL: NB fits better and simpler model is best now: nbU_D2.0


#Now, compare ZI and non-ZI negative binomial fits.
lrtest(rcrt.nbU_D2,rcrt.zinbU_D2.0)
summary(rcrt.zinbU_D2.0)
summary(rcrt.nbU_D2)
plot(resid(rcrt.zinbU_D2.0,type="pearson")~utrans,ylim=c(-1,8))
plot(resid(rcrt.nbU_D,type="pearson")~utrans,ylim=c(-1,8))
plot(resid(rcrt.zinbU_D2.0,type="pearson")~fitted(rcrt.zinbU_D2.0),ylim=c(-1,8),xlim=c(-.1,5))
plot(resid(rcrt.nbU_D2,type="pearson")~fitted(rcrt.nbU_D2),ylim=c(-1,8),xlim=c(-.1,5))
plot(resid(rcrt.zinbU_D2.0,type="pearson")~lnspcover,ylim=c(-1,8))
plot(resid(rcrt.nbU_D2,type="pearson")~lnspcover,ylim=c(-1,8))#  the extra ZI parameter does not significantly improve model fit for 
#Overall, the residuals look a bit better without Zero Inflation

#For PL, best model by AIC is rcrt.lmeU.Y:  AIC=415.3
xyplot(resid(rcrt.lmeU.Y,type="pearson")~utrans)
xyplot(resid(rcrt.lmeU.Y,type="pearson")~lnspcover|utrans)
plot(resid(rcrt.zinbU_D2.0,type="pearson")~lnspcover)
plot(resid(rcrt.nbU_D2,type="pearson")~lnspcover)#  the extra ZI parameter does not significantly improve model fit for 


#Let's plot the competing models

# #Plot quasipoisson spcover(^2) and transect
# ns <- length(unique(utrans))
# x <- seq(0,8,.1)
# fit <- exp(coef(rcrt.glmqU_D2)[1]+x*(coef(rcrt.glmqU_D2)[ns+1])+x^2*(coef(rcrt.glmqU_D2)[ns+2]))
# plot(lnspcover[utrans==unique(utrans)[1]],recruits[utrans==unique(utrans)[1]],col=1,xlab="Area Covered by Conspecifics",ylim=c(-0.5,11.5),xlim=c(0,8),ylab="Number of Recruits")
# lines(x,fit,lty=1,col=1)
# 
# for (i in 2:ns) {
#   fit<-exp(coef(rcrt.glmqU_D2)[1] + coef(rcrt.glmqU_D2)[i] +x*(coef(rcrt.glmqU_D2)[ns+1])+x^2*(coef(rcrt.glmqU_D2)[ns+2]))
#   points(lnspcover[utrans==unique(utrans)[i]],recruits[utrans==unique(utrans)[i]],col=i)
#   lines(x,fit,lty=1,col=i)
# }
# legend("topleft",legend=unique(utrans),col=1:ns,lty=1,pt.cex=0.3)
# title("Overdispersed Poisson fit")

#Plot rcrt.nbU_D2: negative binomial spcover and transect
ns <- length(unique(utrans))
x <- seq(range(lnspcover)[1],range(lnspcover)[2],.1)
fit <- exp(coef(rcrt.nbU_D2)[1]+x*(coef(rcrt.nbU_D2)[ns+1])+x^2*(coef(rcrt.nbU_D2)[ns+2]))
plot(lnspcover[utrans==unique(utrans)[1]],recruits[utrans==unique(utrans)[1]],col=1,xlab="Area Covered by Conspecifics",ylim=c(-0.5,(max(recruits)+0.5)),xlim=range(lnspcover),ylab="Number of Recruits")
lines(x,fit,lty=1,col=1)

for (i in 2:ns) {
  fit<-exp(coef(rcrt.nbU_D2)[1] + coef(rcrt.nbU_D2)[i] +x*coef(rcrt.nbU_D2)[ns+1]+x^2*coef(rcrt.nbU_D2)[ns+2])
  points(lnspcover[utrans==unique(utrans)[i]],recruits[utrans==unique(utrans)[i]],col=i)
  lines(x,fit,lty=1,col=i)
}
legend("topleft",legend=unique(utrans),col=1:ns,lty=1,pt.cex=0.3)
title("Negative Binomial fit")

detach(MCr)
#################################################################################
# Size distribution of recruits

attach(MCrs)
library(lmtest)
hist(lnarea_t[fate_t=="R"])
histogram(~lnarea_t[fate_t=="R"]|utrans[fate_t=="R"],type="count",xlab="Log(Area) of Recruit",main=list(label=sp))

srcrt.glsU <- gls(lnarea_t~utrans)
srcrt.glsUY <- gls(lnarea_t~utrans+year_t)
lrtest(srcrt.glsUY,srcrt.glsU)
vU<- varIdent(form=~1|utrans)
srcrt.glsU.vU <- gls(lnarea_t~utrans,weights=vU)
AIC(srcrt.glsU.vU,srcrt.glsU)
#MC: AIC=603 for varIdent; AIC=609.6 for one var

plot(resid(srcrt.glsU.vU,type="pearson")~utrans)
plot(resid(srcrt.glsU,type="pearson")~utrans)
plot(resid(srcrt.glsU.vU,type="pearson")~fitted(srcrt.glsU.vU))
plot(resid(srcrt.glsU,type="pearson")~fitted(srcrt.glsU))

summary(srcrt.glsU.vU)
detach(MCrs)
