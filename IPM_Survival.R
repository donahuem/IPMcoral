#Model survival using mixed-effects logistic regression
library(lme4) 
sp<-"Montipora capitata"
dat<-MC
nx<-10 # number of sites for this species
x <-seq(min(dat$size,na.rm=TRUE),max(dat$size,na.rm=TRUE),.1)


#What are we modeling here?
plot(dat$size,dat$Psurv,col=c(1:8)[dat$utrans],xlab="Log(Area)[t]", ylab="Survival")
legend("right",legend=levels(dat$utrans),col=c(1:8),pch=1)
xyplot(Psurv~size|utrans,xlab="Log(Area)[t]",data=dat, ylab="Survival",auto.key=list(columns=2), main=list(label=sp))

#first, just look at FE model
surv.full<-glm(Psurv~size*utrans+size*fyear+utrans*fyear,data=dat,family=binomial)
drop1(surv.full,test="Chi")
#drop size *fyear p-value=0.3999
surv.U0Y<-glm(Psurv~size*utrans+utrans*fyear,data=dat,family=binomial)
drop1(surv.U0Y,test="Chi")
#drop size*utrans p-value=0.2507
surv.U0Y0<-glm(Psurv~size+utrans*fyear,data=dat,family=binomial)
drop1(surv.U0Y0,test="Chi")
#keep utrans*fyear p-value=6.904e-06 and size main effect p-value=2.2e-16
##utrans*fyear interaction significantly improves fit but wouldn't include as a fixed effect so look at FE before RE## 

surv.U<-glm(Psurv~size*utrans,data=dat,family=binomial)
drop1(surv.U,test="Chi")
#drop size*utrans p-value=0.5459
surv.U0<-glm(Psurv~size+utrans,data=dat,family=binomial)
drop1(surv.U0,test="Chi")
##Keep size p-value=2.2e-16 and utrans p-value=1.278e-06 main effects##

############################################################################################################################################
#Next, try RE for quad and year
surv.U.Q0Y0 <- glmer(Psurv~size*utrans +(1|uyear),data=dat,family=binomial)
surv.U.Q0 <- glmer(Psurv~size*utrans +(1|uquad),data=dat,family=binomial)
surv.U.Y0 <- glmer(Psurv~size*utrans+(1|fyear),data=dat,family=binomial)
#compare 
anova(surv.U.Q0Y0,surv.U.Q0,surv.U.Y0)
##surv.U.Q0Y0 better by AIC but we do not want to use our concatenated year/quad##
##If we don't want to use uyear, fyear is next best and decreases AIC from FE (both with size+utrans only) by delta=11.724##
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
