#Model survival using mixed-effects logistic regression
sp<-"Montipora capitata"
attach(MC)
library(lme4) 

#What are we modeling here?
plot(lnarea_t,surv,col=fate_t1,xlab="Log(Area)[t]", ylab="Survival")
  legend("bottomright",legend=levels(fate_t1)[2:6],col=2:6,pch=1)
xyplot(surv~lnarea_t|utrans,groups=fate_t1,xlab="Log(Area)[t]", ylab="Survival",auto.key=list(columns=2), main=list(label=sp))

#test for effect of size on survival 

#first, just look at FE model
surv.U <- glm(surv~size*utrans,family=binomial)
surv.U0 <- glm(surv~size+utrans,family=binomial)
AIC(surv.U0,surv.U)
#surv.U0 has lower AIC and none of the lnarea*utrans slopes are sig, indicating that there is no site * size interaction

#Next, add RE for quad and year
surv.U0.Q0Y0 <- glmer(surv~clnarea_t+utrans +(1|uquad/uyear),family=binomial)
#estimated variance for uquad is so small as to be meaningless
surv.U0.Y0 <- glmer(surv~clnarea_t+utrans +(1|uyear),family=binomial)
#compare 
anova(surv.U0.Q0Y0,surv.U0.Y0)
#simpler is better: I don't understand why, but no variance is attributed to uquad
#refit with non-centered size
surv.A_U.Y0<- glmer(surv~lnarea_t+utrans +(1|uyear),family=binomial)
AIC(surv.A_U.Y0)

#To look at variance components, use RE model for all
surv.U0Q0Y0 <- glmer(surv~clnarea_t +(1|utrans/uquad/fyear_t),family=binomial)
#wrong syntax: surv.U0Q0Y0.2 <- glmer(surv~clnarea_t +(1|utrans) + (1|uquad) + (1|fyear_t),family=binomial)
#equiv to correct: surv.U0Q0Y0.3 <- glmer(surv~clnarea_t +(1|utrans/uquad/uyear),family=binomial)
#equiv to correct: surv.U0Q0Y0.4 <- glmer(surv~clnarea_t +(1|utrans) + (1|uquad) + (1|uyear),family=binomial)
surv.UQ0Y0 <- glmer(surv~clnarea_t +(clnarea_t|utrans)+(1|uquad)+(1|uyear),family=binomial)
#compare
anova(surv.U0Q0Y0,surv.UQ0Y0)
#simpler is better: U0Q0Y0 for all RE model

#Plot the survivorship function
#first, refit on natural scale (not centered)

x <-seq(min(lnarea_t),max(lnarea_t),.1)
ns <- length(unique(utrans))
b <-attr(surv.A_U.Y0,"fixef")[1]+x*attr(surv.A_U.Y0,"fixef")[2]
pxz <- exp(b)/(1+exp(b))
cl<-rainbow(ns)
#plot(exp(x+mean(lnarea_t)),pxz,col=1,type="l",ylim=c(0,1))
plot(x,pxz,col=cl[1],type="l",ylim=c(0,1),xlab="Log(Area)[t]",ylab="Prob(Survival)",cex.label=3)
points(lnarea_t[utrans==unique(utrans)[1]],surv[utrans==unique(utrans)[1]],col=cl[1])
for (i in 1:(ns-1)) {
  b <- attr(surv.A_U.Y0,"fixef")[1]+attr(surv.A_U.Y0,"fixef")[i+2]+x*attr(surv.A_U.Y0,"fixef")[2]
  pxz <- exp(b)/(1+exp(b))
  #lines(exp(x+mean(lnarea_t)),pxz,col=i+1)}
  lines(x,pxz,col=cl[i+1])
  points(lnarea_t[utrans==unique(utrans)[i+1]],fiss[utrans==unique(utrans)[i+1]],col=cl[i+1])
  }
legend("right",legend=unique(utrans),col=cl[1:ns],lty=1,pt.cex=0.3,bty="n")
title(sp)

detach(MC)
