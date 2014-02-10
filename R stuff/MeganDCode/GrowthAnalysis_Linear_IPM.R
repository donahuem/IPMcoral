#Model growth with linear mixed-effects model in nlme
library(nlme) 
library(lattice)

#cont<-nls.control(maxiter=10000,minFactor=1/(2^12))
lmc<-lmeControl(niter=10000,msMaxIter=10000)

#Montipora capitata - growth including S, fissers and fusees
attach(MCa)
#attach(PLa)
#******************************************
#Should I include maxWH and depth?  
#These values correspond to the transect and site scales - can I ignore them in the analysis and deal with
#them post-hoc by, for instance, looking at whether transect blups vary with wavex or depth?
#Ignore for now and analyze blups later


#*******************************************
#First step is a linear model of growth
grow <- lm(lnarea_t1 ~ clnarea_t)
summary(grow) 
x<-seq(min(lnarea_t),max(lnarea_t),.1)
plot(x,coef(grow)[1]+coef(grow)[2]*(x-mean(lnarea_t)),type="l")
points(lnarea_t,lnarea_t1,col=utrans)
lines(c(0,7),c(0,7),lty=2)
plot(lnarea_t,resid(grow),col=utrans)

#Next step is to look at effect of utrans (commented this out b/c requires a data= )
#grow.lisU <- lmList(lnarea_t1~clnarea_t|utrans,data=PLa)
#grow.lisU <- lmList(lnarea_t1~clnarea_t|utrans,data=PCa)
#grow.lisU <- lmList(lnarea_t1~clnarea_t|utrans)
#summary(grow.lisU)
coplot(lnarea_t1~lnarea_t|utrans,xlab="Log(Area)[t]",ylab="Log(Area)[t+1]",show.given=FALSE,axlabels=abbreviate(levels(utrans)))
xyplot(lnarea_t1~lnarea_t|utrans,xlab="Log(Area)[t]",ylab="Log(Area)[t+1]",auto.key=list(columns=2),main=list(label=sp))
#coplot(resid(grow.lisU)~clnarea_t|utrans)
#plot(intervals(grow.lisU))


#Now, look at effect of fyear_t
#grow.lisY <-lmList(lnarea_t1~clnarea_t|fyear_t,PLa)
#grow.lisY <-lmList(lnarea_t1~clnarea_t|fyear_t,PCa)
#grow.lisY <-lmList(lnarea_t1~clnarea_t|fyear_t,MCa)
#summary(grow.lisY)
#coplot(lnarea_t1~lnarea_t|fyear_t)
#coplot(resid(grow.lisY)~lnarea_t|fyear_t)
#plot(intervals(grow.lisY))

#Now, fit RE model for utrans
grow.U <- lme(lnarea_t1~clnarea_t,random=~clnarea_t|utrans)
summary(grow.U)
plot(lnarea_t,resid(grow.U))

#Now, fit RE model for fyear_t
grow.Y <- lme(lnarea_t1~clnarea_t,random=~clnarea_t|fyear_t,control=lmc)
#grow.Y <- lme(lnarea_t1~clnarea_t,random=~lnarea_t|fyear_t,MCa,control=cont)
#this doesnt' converge for MCa,PC or PL.  allow only intercept to vary by year
grow.Y0 <- lme(lnarea_t1~clnarea_t,random=~1|fyear_t,control=lmc)
summary(grow.Y0)
plot(lnarea_t,resid(grow.Y0))
grow.Y1 <- lme(lnarea_t1~clnarea_t,random=~clnarea_t-1|fyear_t,control=lmc)
summary(grow.Y1)
plot(lnarea_t,resid(grow.Y0))

#Include both Year and Utrans
#!!!NOTE that the order of the list matters -- the outermost level comes first
#this order implies that for every utrans %in% fyear_t
grow.YU <- lme(lnarea_t1~clnarea_t,random=list(fyear_t=~clnarea_t,utrans=~clnarea_t))
#this doesn't coverge for PL
grow.Y0U <- lme(lnarea_t1~clnarea_t,random=list(fyear_t=~1,utrans=~clnarea_t))
grow.Y1U <- lme(lnarea_t1~clnarea_t,random=list(fyear_t=~clnarea_t-1,utrans=~clnarea_t))
grow.Y0U0 <- lme(lnarea_t1~clnarea_t,random = list(fyear_t=~1,utrans=~1))
grow.Y0U1 <- lme(lnarea_t1~clnarea_t,random = list(fyear_t=~-1,utrans=~clnarea_t-1))
anova(grow.U,grow.Y0,grow.Y1,grow.YU,grow.Y0U,grow.Y1U,grow.Y0U0,grow.Y0U1)#MC & PC
#anova(grow.U,grow.Y0,grow.Y1,grow.Y0U,grow.Y1U,grow.Y0U0,grow.Y0U1)#PL
anova(grow.U,grow.Y0U0)

grow.UY <- lme(lnarea_t1~clnarea_t,random=list(utrans=~clnarea_t,fyear_t=~clnarea_t))
#for PC,PL, doesn't converge
grow.UY0 <- lme(lnarea_t1~clnarea_t,random=list(utrans=~clnarea_t,fyear_t=~1))
grow.UY1 <- lme(lnarea_t1~clnarea_t,random=list(utrans=~lnarea_t,fyear_t=~clnarea_t-1))
grow.U0Y0 <- lme(lnarea_t1~clnarea_t,random = list(utrans=~1,fyear_t=~1))
grow.U1Y0 <- lme(lnarea_t1~clnarea_t,random = list(utrans=~clnarea_t-1,fyear_t=~-1))
anova(grow.U,grow.Y0,grow.Y1,grow.UY,grow.UY0,grow.UY1,grow.U0Y0,grow.U1Y0)#for MC:
#anova(grow.U,grow.Y0,grow.Y1,grow.UY0,grow.UY1,grow.U0Y0,grow.U1Y0)#for PL PC
anova(grow.UY0,grow.U0Y0)
anova(grow.U,grow.UY0,grow.U0Y0)#for PL
#Best model for MCa is grow.U0Y0.  RE for utrans and year%in%utrans for intercept
#Best model for PCa is grow.U1Y0 or grow.U0Y0 - keep the same for simplicity?
#Best model for PLa is grow.U0Y0.
#difference between U0Y0 and Y0U0 is tricky -- Y0U0 envisions the sampling as going out every year to measure sites (and quads in sites)
#  while U0Y0 envisions returning to sites and doing repeated measures on those sites -- this is how I think of it.

#Quad? Yup for MC.  Yes but Less convincing for PC and PL.
grow.U0Q0Y0 <- lme(lnarea_t1~clnarea_t,random = list(utrans=~1,uquad=~1,fyear_t=~1))
anova(grow.U0Y0,grow.U0Q0Y0)

#Get BLUPs for utrans level
coef(grow.U0Y0,level=1) #BLUPs for utrans
coef(grow.U0Y0,level=2) #BLUPs for year within utrans

#Alternatively, consider a FE model for utrans with RE for fyear_t
grow.U0.Y0 <- lmer(lnarea_t1~clnarea_t+utrans +(1|fyear_t))
grow.U.Y0 <- lmer(lnarea_t1~clnarea_t*utrans +(1|fyear_t))

#To compare FE and RE, need to fit both with ML b/c different fixed effects
grow.U.Y0.ml <- update(grow.U.Y0,method="ML")
grow.U0.Y0.ml <- update(grow.U0.Y0,method="ML")
grow.U0Y0.ml <- update(grow.U0Y0,method="ML")
anova(grow.U.Y0.ml,grow.U0.Y0.ml,grow.U0Y0.ml)
anova(grow.U.Y0.ml,grow.U0Y0.ml)
#FOR MC: to the extent that this comparison is legitimate, grow.U.Y.ml with utrans as FE for intercept and slope
#  is somewhat better than nested REs
#FOR PC: the RE model grow.U0Y0 is better than the FE
#FOR ML: the FE model grow.U.Y0 is better than RE grow.U0Y0

#Add Quad as RE to FE model?
grow.U0.Q0Y0 <- lme(lnarea_t1~clnarea_t+utrans,random = ~1|uquad/fyear_t)
grow.U.Q0Y0 <- lme(lnarea_t1~clnarea_t*utrans,random = ~1|uquad/fyear_t)
#U.Q0Y0 is better (del AIC=3, LRT=0.0125)

#to compare with models with different REs, fit with ML
grow.U0Q0Y0.ml <- update(grow.U0Q0Y0,method="ML")
grow.U0.Q0Y0.ml <- update(grow.U0.Q0Y0,method="ML")
grow.U.Q0Y0.ml <- update(grow.U.Q0Y0,method="ML")
anova(grow.U0.Q0Y0.ml,grow.U.Q0Y0.ml,grow.U0Q0Y0.ml)
anova(grow.U.Q0Y0.ml,grow.U0Y0.ml)
#FOR MC & PL: Best model is grow.U.Q0Y0
summary(grow.U.Q0Y0)
#For PC: Best model is a tossup: grow.U0.Q0Y0, grow.U0Q0Y0, grow.U.Q0Y0
summary(grow.U0Q0Y0)


#Include variance decreasing with lnarea_t?  YES for MC, YES for PC, YES for PL
plot(predict(grow.U.Q0Y0),resid(grow.U.Q0Y0))
grow.U.Q0Y0.vP <- update(grow.U.Q0Y0,weights=varPower(form=~clnarea_t))
anova(grow.U.Q0Y0.vP,grow.U.Q0Y0)
#FINAL LINEAR MODEL for MC:  
#Lt1 ~ b0 + b1.Lt + b1.utrans + b3.Lt:utrans
# modeling var(b0) ~ quad + year%in%quad and var(err)~varP(Lt)

#Fit a model without centering and wo quad-var for the IPM
grow.U.Y0.vP <- lme(lnarea_t1~lnarea_t*utrans,random = ~1|fyear_t, weights=varPower(form=~lnarea_t))



#Add coral cover? (spcover_t MAYBE; allcover_t NO) - leave this for later
grow.ccU.Q0Y0.ml <- lme(lnarea_t1~clnarea_t*utrans+spcover_t,random = ~1|uquad/fyear_t,method="ML")
grow.acU.Q0Y0.ml <- lme(lnarea_t1~clnarea_t*utrans+allcover_t,random = ~1|uquad/fyear_t,method="ML")
anova(grow.U.Q0Y0.ml,grow.ccU.Q0Y0.ml)
anova(grow.U.Q0Y0.ml,grow.acU.Q0Y0.ml)

#Hmm, perhaps coral cover can account for some of the variation attributed at the quad scale
#use the fullest fixed model and test for REs
grow.ccAU.Y0 <- lmer(lnarea_t1~clnarea_t*utrans*log(spcover_t) +(1|fyear_t))
grow.cc_AU.Y0 <- lmer(lnarea_t1~clnarea_t*utrans+spcover_t +(1|fyear_t))
grow.ccAU.Q0Y0 <- lmer(lnarea_t1~clnarea_t*utrans*spcover_t + (1|uquad/fyear_t))
grow.cdAU.Q0Y0 <- lmer(lnarea_t1~clnarea_t*utrans*spcoverD +(1|uquad/fyear_t))
#Yup, grow.ccAU.Y0 (AIC=2449) fits better than grow.ccAU.Q0Y0 (AIC=2456)
#But, this may be driven by correlation between larger colonies having larger total cover in the quad.
xyplot(lnarea_t1~log(spcover_t)|utrans)
xyplot(lnarea_t1~log(spcover_t-area_t)|utrans)

#Use grow.ccAU.Y0 as best model
x1 <-seq(min(lnarea_t),max(lnarea_t),.1)
x2 <-seq(min(spcover_t),max(spcover_t),1)


#Plot growth curves
#For all 3 species use grow.U.Q0Y0.vP
#adjust axes so x is lnarea_t
x <-seq(min(lnarea_t),max(lnarea_t),.1)
ns <- length(unique(utrans))
cl<-rainbow(ns)
y <-fixef(grow.U.Q0Y0.vP)[1]+x*fixef(grow.U.Q0Y0.vP)[2]-mean(lnarea_t)*fixef(grow.U.Q0Y0.vP)[2]
plot(x,y,col=cl[1],type="l",xlab="Log(Area)[t]",ylab="Log(Area)[t+1]")
points(lnarea_t[utrans==unique(utrans)[1]],lnarea_t1[utrans==unique(utrans)[1]],col=cl[1])
lines(c(-2,10),c(-2,10),lty=2)
for (i in 2:ns) {
  y <- fixef(grow.U.Q0Y0.vP)[1]+fixef(grow.U.Q0Y0.vP)[i+1]+x*(fixef(grow.U.Q0Y0.vP)[2]+fixef(grow.U.Q0Y0.vP)[i+ns])-(fixef(grow.U.Q0Y0.vP)[2]+fixef(grow.U.Q0Y0.vP)[i+ns])*mean(lnarea_t)
  lines(x,y,col=cl[i],type="l")
  points(lnarea_t[utrans==unique(utrans)[i]],lnarea_t1[utrans==unique(utrans)[i]],col=cl[i])
  lines(c(-2,10),c(-2,10),lty=2)
  }
legend("bottomright",legend=unique(utrans),col=cl[1:ns],lty=1,pt.cex=0.3,bty="n")
title(paste(sp,": VarFunc=power(lnarea)",""))


detach(MCa)
