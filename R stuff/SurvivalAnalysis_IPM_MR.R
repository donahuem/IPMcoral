#Model survival using mixed-effects logistic regression
sp<-"Montipora capitata"
dat<-MC
library(lme4) 

#What are we modeling here?
plot(dat$size,dat$Psurv,col=c(1:8)[dat$utrans],xlab="Log(Area)[t]", ylab="Survival")
legend("right",legend=levels(dat$utrans),col=c(1:8),pch=1)
xyplot(Psurv~size|utrans,xlab="Log(Area)[t]",data=dat, ylab="Survival",auto.key=list(columns=2), main=list(label=sp))

#first, just look at FE model
surv.full<-glm(Psurv~size*utrans+size*fyear+utrans*fyear,data=dat,family=binomial)
drop1(surv.full,test="Chi")
surv.U0Y<-glm(Psurv~size*fyear+utrans*fyear,data=dat,family=binomial)
drop1(surv.U0Y,test="Chi")
surv.U0Y0<-glm(Psurv~size+utrans*fyear,data=dat,family=binomial)
drop1(surv.U0Y0,test="Chi")
##utrans*fyear interaction significantly improves fit but wouldn't include as a fixed effect so look at FE before RE## 
surv.U<-glm(Psurv~size*utrans,data=dat,family=binomial)
drop1(surv.U,test="Chi")
surv.U0<-glm(Psurv~size+utrans,data=dat,family=binomial)
drop1(surv.U0,test="Chi")
##Keep utrans main effects (P-value=7.028e-07)##

############################################################################################################################################
#Next, try RE for quad and year
surv.U.Q0Y0 <- glmer(Psurv~size*utrans +(1|uyear),data=dat,family=binomial)
surv.U.Q0 <- glmer(Psurv~size*utrans +(1|uquad),data=dat,family=binomial)
surv.U.Y0 <- glmer(Psurv~size*utrans+(1|fyear),data=dat,family=binomial)
#compare 
anova(surv.U.Q0Y0,surv.U.Q0,surv.U.Y0)
##surv.U.Q0Y0 better by AIC but do we want to use our concatenated year/quad?##
##If we don't want to use uyear, fyear is next best of RE but does not decrease AIC from FE model##
surv.U.Y<-glmer(Psurv~size*utrans+(size|fyear),data=dat,family=binomial)
anova(surv.U.Y,surv.U.Y0)
##surv.U.Y0 is a better fit than surv.U.Y)##
drop1(surv.U.Q0Y0,test="Chi")
##AIC decreases from 1485.6 to 1477.0 when 2-way interaction dropped##
surv.U0.Y0<-glmer(Psurv~size+utrans+(1|fyear),data=dat,family=binomial)
drop1(surv.U0.Y0,test="Chi")
surv.U0.Q0Y0<-glmer(Psurv~size+utrans+(1|uyear),data=dat,family=binomial)
drop1(surv.U0.Q0Y0,test="Chi")


#Plot the survivorship function
x <-seq(min(dat$size,na.rm=TRUE),max(dat$size,na.rm=TRUE),.1)
nx <- length(unique(dat$utrans))
b <-coefficients(surv.U0)[1]+x*coefficients(surv.U0)[2]
pxz <- exp(b)/(1+exp(b))
cl<-rainbow(nx)
plot(x,pxz,col=cl[1],type="l",ylim=c(0,1),xlab="Log(Area)[t]",ylab="Prob(Survival)",cex.label=3)
points(dat$size[dat$utrans==unique(dat$utrans)[1]],dat$Psurv[dat$utrans==unique(dat$utrans)[1]],col=cl[1])
for (i in 1:(nx-1)) {
  b <- coefficients(surv.U0)[1]+coefficients(surv.U0)[i+2]+x*coefficients(surv.U0)[2]
  pxz <- exp(b)/(1+exp(b))
  lines(x,pxz,col=cl[i+1])
  points(dat$size[dat$utrans==unique(dat$utrans)[i+1]],dat$Psurv[dat$utrans==unique(dat$utrans)[i+1]],col=cl[i+1])
  }
legend("right",legend=unique(dat$utrans),col=cl[1:nx],lty=1,pt.cex=0.3,bty="n")
title(sp)
