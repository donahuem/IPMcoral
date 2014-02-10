library(nlme)
library(lme4)#for glmer

#######################Survival##########################
surv.reg.re1<-glmer(Psurv~size+(size|utrans),data=MC,family=binomial)
summary(surv.reg.re1)
surv.int.re1=c(fixed.effects(surv.reg.re1)[1]+random.effects(surv.reg.re1)[[1]][1:10,1])
surv.slope.re1=c(fixed.effects(surv.reg.re1)[2]+random.effects(surv.reg.re1)[[1]][1:10,2])

surv.reg.re01<-glmer(Psurv~size+(size-1|utrans),data=MC,family=binomial)
summary(surv.reg.re01)
surv.int.re01=fixed.effects(surv.reg.re01)[1]
surv.slope.re01=c(fixed.effects(surv.reg.re01)[2]+random.effects(surv.reg.re01)[[1]][1:10,1])

surv.reg.re0<-glmer(Psurv~size+(1|utrans),data=MC,family=binomial)
summary(surv.reg.re0)
surv.int.re0=c(fixed.effects(surv.reg.re0)[1]+random.effects(surv.reg.re0)[[1]][1:10,1])
surv.slope.re0=fixed.effects(surv.reg.re0)[2]

par(mfrow=c(2,2))
for(i in 1:10){
  x<-subset(MC,MC$utrans==levels(MC$utrans)[i])
  plot(x$size,x$Psurv,xlim=c(-2,10),xlab="Size (t)",ylab="Probability of Survival")
  u1<-exp(surv.int[i]+surv.slope*y)
  #u2<-exp(surv.int.re1[i]+surv.slope.re1[i]*y)
  u3<-exp(surv.int.re0[i]+surv.slope.re0*y)
  #u4<-exp(surv.int.re01+surv.slope.re01[i]*y)
  lines(y,u1/(u1+1))
  #lines(y,u2/(u2+1),col=c("purple"))
  lines(y,u3/(u3+1),col=c("red"))
  #lines(y,u4/(u4+1),col=c("green"))
}
# intercept alone is a better fit than slope alone.  Adding slope to intercept model does not improve fit p-value=0.7955
#This agrees with fixed effect model

##############Growth#########################3
growth.reg.re.full<-nlme(sizeNext~a+b*size+c*exp(-d*size),data=GONGS,fixed=a+b+c+d~1,random=a+b+c+d~1|utrans,start=list(fixed=c(1,1,1,1)))

#a only has re by site
growth.reg.re.a0<-nlme(sizeNext~a+b*size+c*exp(-d*size),data=GONGS,fixed=a+b+c+d~1,random=a~1|utrans,start=list(fixed=c(1,1,1,1)),control=gnlsControl(nlsTol=0.1))
growth.int.a0<-c(fixed.effects(growth.reg.re.a0)[1]+random.effects(growth.reg.re.a0)[1:10,1])
growth.slope.a0<-fixed.effects(growth.reg.re.a0)[2]
growth.e1.a0<-fixed.effects(growth.reg.re.a0)[3]
growth.e2.a0<-fixed.effects(growth.reg.re.a0)[4]
growth.sd.a0<-summary(growth.reg.re.a0)$sigma

growth.reg.re.a1<-nlme(sizeNext~a+b*size+c*exp(-d*size),data=GONGS,fixed=a+b+c+d~1,random=a~size|utrans,start=list(fixed=c(1,1,1,1)),control=gnlsControl(nlsTol=0.1))
growth.int.a1<-c(fixed.effects(growth.reg.re.a1)[1]+random.effects(growth.reg.re.a1)[1:10,1]+random.effects(growth.reg.re.a1)[1:10,2]*y)
growth.slope.a1<-fixed.effects(growth.reg.re.a1)[2]
growth.e1.a1<-fixed.effects(growth.reg.re.a1)[3]
growth.e2.a1<-fixed.effects(growth.reg.re.a1)[4]
growth.sd.a1<-summary(growth.reg.re.a1)$sigma


growth.reg.re.a01<-nlme(sizeNext~a+b*size+c*exp(-d*size),data=GONGS,fixed=a+b+c+d~1,random=a~size-1|utrans,start=list(fixed=c(1,1,1,1)),control=gnlsControl(nlsTol=0.1))
growth.int.a01<-c(fixed.effects(growth.reg.re.a01)[1]+random.effects(growth.reg.re.a01)[1:10,1]*y)
growth.slope.a01<-fixed.effects(growth.reg.re.a01)[2]
growth.e1.a01<-fixed.effects(growth.reg.re.a01)[3]
growth.e2.a01<-fixed.effects(growth.reg.re.a01)[4]
growth.sd.a01<-summary(growth.reg.re.a01)$sigma

par(mfrow=c(2,2))
for(i in 1:10){
  x<-subset(GONGS,GONGS$utrans==levels(GONGS$utrans)[i])
  plot(x$size,x$sizeNext,xlim=c(-2,10))
  lines(y,growth.int[i]+growth.slope*y+growth.e1*exp(-1*growth.e2*y),col=c("black"),lwd=2)
  lines(y,growth.int.a0[i]+growth.slope.a0*y+growth.e1.a0*exp(-1*growth.e2.a0*y),col=c("red"))
  lines(y,growth.int.a1[i]+growth.slope.a1*y+growth.e1.a1*exp(-1*growth.e2.a1*y),col=c("purple"))
  lines(y,growth.int.a01[i]+growth.slope.a01*y+growth.e1.a01*exp(-1*growth.e2.a01*y),col=c("green"))
}
AIC(growth.reg.re.a0,growth.reg.re.a1,growth.reg.re.a01)
##growth.reg.re.a1 lowest by delta AIC=19.815
##growth.reg.re.a0 matches FE growth.reg most closely

#b only has random effect by site
growth.reg.re.b<-nlme(sizeNext~a+b*size+c*exp(-d*size),data=GONGS,fixed=a+b+c+d~1,random=b~1|utrans,start=list(fixed=c(1,1,1,1)))

growth.int.b<-fixed.effects(growth.reg.re.b)[1]
growth.slope.b<-c(fixed.effects(growth.reg.re.b)[2]+random.effects(growth.reg.re.b)[1:10,1])
growth.e1.b<-fixed.effects(growth.reg.re.b)[3]
growth.e2.b<-fixed.effects(growth.reg.re.b)[4]
growth.sd.b<-summary(growth.reg.re.b)$sigma
par(mfrow=c(2,2))
for(i in 1:10){
  x<-subset(GONGS,GONGS$utrans==levels(GONGS$utrans)[i])
  plot(x$size,x$sizeNext)
  lines(y,growth.int[i]+growth.slope*y+growth.e1*exp(-1*growth.e2*y),col=c("black"),lwd=2)
  lines(y,growth.int.b+growth.slope.b[i]*y+growth.e1.b*exp(-1*growth.e2.b*y),col=c("red"))
}
#c only has RE by site
growth.reg.re.c<-nlme(sizeNext~a+b*size+c*exp(-d*size),data=GONGS,fixed=a+b+c+d~1,random=c~1|utrans,start=list(fixed=c(1,1,1,1)))

growth.int.c<-fixed.effects(growth.reg.re.c)[1]
growth.slope.c<-fixed.effects(growth.reg.re.c)[2]
growth.e1.c<-c(fixed.effects(growth.reg.re.c)[3]+random.effects(growth.reg.re.c)[1:10,1])
growth.e2.c<-fixed.effects(growth.reg.re.c)[4]
growth.sd.c<-summary(growth.reg.re.c)$sigma
par(mfrow=c(2,2))
for(i in 1:10){
  x<-subset(GONGS,GONGS$utrans==levels(GONGS$utrans)[i])
  plot(x$size,x$sizeNext)
  lines(y,growth.int[i]+growth.slope*y+growth.e1*exp(-1*growth.e2*y),col=c("black"),lwd=2)
  lines(y,growth.int.c+growth.slope.c*y+growth.e1.c[i]*exp(-1*growth.e2.c*y),col=c("red"))
}
#d only RE by site
growth.reg.re.d<-nlme(sizeNext~a+b*size+c*exp(-d*size),data=GONGS,fixed=a+b+c+d~1,random=d~size-1|utrans,start=list(fixed=c(1,1,1,1)))
growth.int.d<-fixed.effects(growth.reg.re.d)[1]
growth.slope.d<-fixed.effects(growth.reg.re.d)[2]
growth.e1.d<-fixed.effects(growth.reg.re.d)[3]
growth.e2.d<-c(fixed.effects(growth.reg.re.d)[4]+random.effects(growth.reg.re.d)[1:10,1])
growth.sd.d<-summary(growth.reg.re.d)$sigma

par(mfrow=c(2,2))
for(i in 1:10){
  x<-subset(GONGS,GONGS$utrans==levels(GONGS$utrans)[i])
  plot(x$size,x$sizeNext)
  lines(y,growth.int[i]+growth.slope*y+growth.e1*exp(-1*growth.e2*y),col=c("black"),lwd=2)
  lines(y,growth.int.d+growth.slope.d*y+growth.e1.d*exp(-1*growth.e2.d[i]*y),col=c("red"))
}

AIC(growth.reg.re.a0,growth.reg.re.b,growth.reg.re.c,growth.reg.re.d)
AIC(growth.reg.re.a,growth.reg.re.b,growth.reg.re.c,growth.reg.re.d)
AIC(growth.reg.re.a,growth.reg.re.b,growth.reg.re.c,growth.reg.re.d)


par(mfrow=c(2,2))
for(i in 1:10){
  x<-subset(GONGS,GONGS$utrans==levels(GONGS$utrans)[i])
  plot(x$size,x$sizeNext,xlim=c(-2,10),ylim=c(-2,10),xlab="Size (t)",ylab="Sizt (t+1)")
  lines(y,growth.int[i]+growth.slope*y+growth.e1*exp(-1*growth.e2*y),col=c("black"),lwd=2)
  lines(y,growth.int.a0[i]+growth.slope.a0*y+growth.e1.a0*exp(-1*growth.e2.a0*y),col=c("red"))
  lines(y,growth.int.b+growth.slope.b[i]*y+growth.e1.b*exp(-1*growth.e2.b*y),col=c("purple"))
  lines(y,growth.int.c+growth.slope.c*y+growth.e1.c[i]*exp(-1*growth.e2.c*y),col=c("green"))
}


##########################ES Growth#################
es.reg.re0<-lme(sizeNext~size,data=GOES,random=~1|utrans)
es.int.re0<-c(fixed.effects(es.reg.re0)[1]+random.effects(es.reg.re0)[1:10,1])
es.slope.re0<-fixed.effects(es.reg.re0)[2]
es.sd.re0<-(summary(es.reg.re0)$sigma)

es.reg.re1<-lme(sizeNext~size,data=GOES,random=~size|utrans,control=lmeControl(maxIter=100))
es.int.re1<-c(fixed.effects(es.reg.re01)[1]+random.effects(es.reg.re0)[1:10,1])
es.slope.re1<-c(fixed.effects(es.reg.re1)[2]+random.effects(es.reg.re1)[1:10,2])
es.sd.re1<-(summary(es.reg.re1)$sigma)
##Error message iteration
es.reg.re01<-lme(sizeNext~size,data=GOES,random=~size-1|utrans)
es.int.re01<-fixed.effects(es.reg.re01)[1]
es.slope.re01<-c(fixed.effects(es.reg.re01)[2]+random.effects(es.reg.re01)[1:10,1])
es.sd.re01<-(summary(es.reg.re01)$sigma)

par(mfrow=c(2,2))
for(i in 1:10){
  x<-subset(GOES,GOES$utrans==levels(GOES$utrans)[i])
  plot(x$size,x$sizeNext,xlim=c(-2,10),ylim=c(-2,10))
  lines(y,es.int+es.slope*y,lwd=4)
  lines(y,es.int.re0[i]+es.slope.re0*y,col=c("red"))
  lines(y,es.int.re01+es.slope.re01[i]*y,col=c("green"))
}
AIC(es.reg.re1,es.reg.re01,es.reg.re0)
anova(es.reg,es.reg.re0)
# addition of RE does not improve the fit p-value=0.9999
#############################Fission Growth##############################
fiss.full<-gls(sizeNext~size,data=fissonly)
fiss.int<-coefficients(fiss.full)[1]
fiss.slope<-coefficients(fiss.full)[2]
fiss.sd<-summary(fiss.full)$sigma

fiss.full.re1<-lme(sizeNext~size,random=~size|utrans,data=fissonly)
fiss.int.re1<-c(fixef(fiss.full.re1)[1]+ranef(fiss.full.re1)[1:10,1])
fiss.slope.re1<-c(fixef(fiss.full.re1)[2]+ranef(fiss.full.re1)[1:10,2])
fiss.sd.re1<-summary(fiss.full.re1)$sigma

fiss.full.re01<-lme(sizeNext~size,random=~size-1|utrans,data=fissonly)
fiss.int.re01<-fixef(fiss.full.re01)[1]
fiss.slope.re01<-c(fixef(fiss.full.re01)[2]+ranef(fiss.full.re01)[1:10,1])
fiss.sd.re01<-summary(fiss.full.re01)$sigma

fiss.full.re0<-lme(sizeNext~size,random=~1|utrans,data=fissonly)
fiss.int.re0<-c(fixef(fiss.full.re0)[1]+ranef(fiss.full.re0)[1:10,1])
fiss.slope.re0<-fixef(fiss.full.re0)[2]
fiss.sd.re0<-summary(fiss.full.re0)$sigma

par(mfrow=c(2,2))
for(i in 1:10){
  x<-subset(fissonly,fissonly$utrans==levels(fissonly$utrans)[i])
  plot(x$size,x$sizeNext,xlim=c(-2,10),ylim=c(-2,10))
  lines(y,fiss.int+fiss.slope*y)
  lines(y,fiss.int.re1[i]+fiss.slope.re1[i]*y,col=c("red"))
  #lines(y,fiss.int.re0[i]+fiss.slope.re0*y,col=c("purple"))
  #lines(y,fiss.int.re01+fiss.slope.re01[i]*y,col=c("green"))
}
AIC(fiss.full.re1,fiss.full.re01,fiss.full.re0)

anova(fiss.full.re0,fiss.full)
#Adding RE does not improve model p-value= 0.9998
#######################fusion Growth####################

fuse.full<-gls(sizeNext~size,data=fuseonly)
summary(fuse.full)
fuse.int<-coefficients(fuse.full)[1]
fuse.slope<-coefficients(fuse.full)[2]
fuse.sd<-summary(fuse.full)$sigma

fuse.full.re1<-lme(sizeNext~size,random=~size|utrans,data=fuseonly)
summary(fuse.full.re1)
fuse.int.re1<-c(fixef(fuse.full.re1)[1]+ranef(fuse.full.re1)[1:10,1])
fuse.slope.re1<-c(fixef(fuse.full.re1)[2]+ranef(fuse.full.re1)[1:10,2])
fuse.sd.re1<-summary(fuse.full.re1)$sigma

fuse.full.re01<-lme(sizeNext~size,random=~size-1|utrans,data=fuseonly)
summary(fuse.full.re01)
fuse.int.re01<-fixef(fuse.full.re01)[1]
fuse.slope.re01<-c(fixef(fuse.full.re01)[2]+ranef(fuse.full.re01)[1:10,1])
fuse.sd.re01<-summary(fuse.full.re01)$sigma
                 
fuse.full.re0<-lme(sizeNext~size,random=~1|utrans,data=fuseonly)
summary(fuse.full.re0)
fuse.int.re0<-c(fixef(fuse.full.re0)[1]+ranef(fuse.full.re0)[1:10,1])
fuse.slope.re0<-fixef(fuse.full.re0)[2]
fuse.sd.re0<-summary(fuse.full.re0)$sigma

par(mfrow=c(2,2))
for(i in 1:10){
  x<-subset(fuseonly,fuseonly$utrans==levels(fuseonly$utrans)[i])
  plot(x$size,x$sizeNext,xlim=c(-2,10),ylim=c(-2,10))
  lines(y,fuse.int+fuse.slope*y)
  #lines(y,fuse.int.re1[i]+fuse.slope.re1[i]*y,col=c("red"))
  lines(y,fuse.int.re0[i]+fuse.slope.re0*y,col=c("purple"))
  lines(y,fuse.int.re01+fuse.slope.re01[i]*y,col=c("green"))
}
AIC(fuse.full.re01,fuse.full.re0)

anova(fuse.full.re0,fuse.full)
##RE does not improve fit p-value=0.9998

##################################Probability of Extreme shrinkage
Pes.re1<-glmer(PES~size+(size|utrans),family=binomial,data=growthonly)
pes.int.re1<-c(fixef(Pes.re1)[1]+ranef(Pes.re1)[[1]][1:10,1])
pes.slope.re1<-c(fixef(Pes.re1)[2]+ranef(Pes.re1)[[1]][1:10,2])

Pes.re01<-glmer(PES~size+(size-1|utrans),family=binomial,data=growthonly)
pes.int.re01<-fixef(Pes.re01)[1]
pes.slope.re01<-c(fixef(Pes.re01)[2]+ranef(Pes.re01)[[1]][1:10,1])

Pes.re0<-glmer(PES~size+(1|utrans),family=binomial,data=growthonly)
pes.int.re0<-c(fixef(Pes.re0)[1]+ranef(Pes.re0)[[1]][1:10,1])
pes.slope.re0<-fixef(Pes.re0)[2]

par(mfrow=c(2,2))
for(i in 1:10){
 x<-subset(growthonly,growthonly$utrans==levels(growthonly$utrans)[i])
 plot(x$size,x$Pes,xlim=c(-2,10),ylim=c(0,1))
 u1<-exp(pes.int)
 u2<-exp(pes.int.re1[i]+pes.slope.re1[i]*y)
 u3<-exp(pes.int.re01+pes.slope.re01[i]*y)
 u4<-exp(pes.int.re0[i]+pes.slope.re0*y)
 lines(y,rep(u1/(u1+1),300),col=c("black"),lwd=2)
 lines(y,u2/(u2+1),col=c("red"))
 lines(y,u3/(u3+1),col=c("purple"))
 lines(y,u4/(u4+1),col=c("green"))
}
AIC(Pes.re1,Pes.re01,Pes.re0)
anova(Pes.re1,Pes.re0)
anova(Pes.re1,Pes.re01)
anova(Pes.re0,Pes.re01)                

Pes.re01<-glmer(PES~1+(1|utrans),family=binomial,data=growthonly)
pes.int.re01<-fixef(Pes.re01)[1]
pes.slope.re01<-c(fixef(Pes.re01)[2]+ranef(Pes.re01)[[1]][1:10,1])

####RE does not appear to improve model.  n so low that breaking it up by utrans is silly.
#################################Probability of Fission########
fiss.re1=glmer(Pfiss~size+(size|utrans),data=MC,family=binomial)
summary(fiss.re1)
Pfiss.int.re1=c(fixed.effects(fiss.re1)[1]+random.effects(fiss.re1)[[1]][1:10,1])
Pfiss.slope.re1=c(fixed.effects(fiss.re1)[2]+random.effects(fiss.re1)[[1]][1:10,2])

fiss.re01=glmer(Pfiss~size+(size-1|utrans),data=MC,family=binomial)
summary(fiss.re01)
Pfiss.int.re01=fixed.effects(fiss.re01)[1]
Pfiss.slope.re01=c(fixed.effects(fiss.re01)[2]+random.effects(fiss.re01)[[1]][1:10,1])

fiss.re0=glmer(Pfiss~size+(1|utrans),data=MC,family=binomial)
summary(fiss.re0)
Pfiss.int.re0=c(fixed.effects(fiss.re0)[1]+random.effects(fiss.re0)[[1]][1:10,1])
Pfiss.slope.re0=fixed.effects(fiss.re0)[2]

par(mfrow=c(2,2))
for(i in 1:10){
  x<-subset(MC,MC$utrans==levels(MC$utrans)[i])
  plot(x$size,x$Pfiss,xlim=c(-2,10),ylim=c(0,1),xlab="Size (t)",ylab="Probability of fissLG")
  u1<-exp(Pfiss.int[i]+Pfiss.slope[i]*y)
  u2<-exp(Pfiss.int.re1[i]+Pfiss.slope.re1[i]*y)
  u3<-exp(Pfiss.int.re01+Pfiss.slope.re01[i]*y)
  u4<-exp(Pfiss.int.re0[i]+Pfiss.slope.re0*y)
  lines(y,u1/(u1+1),col=c("black"),lwd=2)
  lines(y,u2/(u2+1),col=c("red"))
  lines(y,u3/(u3+1),col=c("purple"))
  lines(y,u4/(u4+1),col=c("green"))
}
AIC(fiss.re1,fiss.re01,fiss.re0)
anova(fiss.re1,fiss.re0)
anova(fiss.re1,fiss.re01)
anova(fiss.re0,fiss.re01)

##intercept only is best fit
################################Probablility of Fusion
fuse.re1<-glmer(Pfuse~size+(size|utrans),data=MC,family=binomial)
summary(fuse.re1)
Pfuse.int.re1=c(fixed.effects(fuse.re1)[1]+random.effects(fuse.re1)[[1]][1:10,1])
Pfuse.slope.re1=c(fixed.effects(fuse.re1)[2]+random.effects(fuse.re1)[[1]][1:10,2])

fuse.re01<-glmer(Pfuse~size+(size-1|utrans),data=MC,family=binomial)
summary(fuse.re01)
Pfuse.int.re01=fixed.effects(fuse.re01)[1]
Pfuse.slope.re01=c(fixed.effects(fuse.re01)[2]+random.effects(fuse.re01)[[1]][1:10,1])

fuse.re0<-glmer(Pfuse~size+(1|utrans),data=MC,family=binomial)
summary(fuse.re0)
Pfuse.int.re0=c(fixed.effects(fuse.re0)[1]+random.effects(fuse.re0)[[1]][1:10,1])
Pfuse.slope.re0=fixed.effects(fuse.re0)[2]

par(mfrow=c(2,2))
for(i in 1:10){
  x<-subset(MC,MC$utrans==levels(MC$utrans)[i])
  plot(x$size,x$Pfuse,xlim=c(-2,10),ylim=c(0,1),xlab="size (t)",ylab="Probability of fusion")
  u1<-exp(Pfuse.int[i]+Pfuse.slope*y)
  u2<-exp(Pfuse.int.re1[i]+Pfuse.slope.re1[i]*y)
  u3<-exp(Pfuse.int.re01+Pfuse.slope.re01[i]*y)
  u4<-exp(Pfuse.int.re0[i]+Pfuse.slope.re0*y)
  lines(y,u1/(u1+1),col=c("black"),lwd=2)
  lines(y,u2/(u2+1),col=c("red"))
  lines(y,u3/(u3+1),col=c("purple"))
  lines(y,u4/(u4+1),col=c("green"))
}
AIC(fuse.re1,fuse.re01,fuse.re0)
anova(fuse.re1,fuse.re0)
anova(fuse.re1,fuse.re01)
anova(fuse.re0,fuse.re01)

##incercept varies with site is best model of these this goes along with FE model


####################### Probability of large fission#####################
fuseLG.re1=glmer(PfuseLG~size+(size|utrans),data=MC,family=binomial)
summary(fuseLG.re1)
PfuseLG.int.re1=c(fixed.effects(fuseLG.re1)[1]+random.effects(fuseLG.re1)[[1]][1:10,1])
PfuseLG.slope.re1=c(fixed.effects(fuseLG.re1)[2]+random.effects(fuseLG.re1)[[1]][1:10,2])

fuseLG.re01=glmer(PfuseLG~size+(size-1|utrans),data=MC,family=binomial)
summary(fuseLG.re01)
PfuseLG.int.re01=fixed.effects(fuseLG.re01)[1]
PfuseLG.slope.re01=c(fixed.effects(fuseLG.re01)[2]+random.effects(fuseLG.re01)[[1]][1:10,1])

fuseLG.re0=glmer(PfuseLG~size+(1|utrans),data=MC,family=binomial)
summary(fuseLG.re0)
PfuseLG.int.re0=c(fixed.effects(fuseLG.re0)[1]+random.effects(fuseLG.re0)[[1]][1:10,1])
PfuseLG.slope.re0=fixed.effects(fuseLG.re0)[2]


par(mfrow=c(2,2))
for(i in 1:10){
  x<-subset(MC,MC$utrans==levels(MC$utrans)[i])
  plot(x$size,x$PfuseLG,xlim=c(-2,10),ylim=c(0,1),xlab="Size (t)",ylab="Probability of FuseLG")
  u1<-exp(PfuseLG.int[i]+PfuseLG.slope*y)
  u2<-exp(PfuseLG.int.re1[i]+PfuseLG.slope.re1[i]*y)
  u3<-exp(PfuseLG.int.re01+PfuseLG.slope.re01[i]*y)
  u4<-exp(PfuseLG.int.re0[i]+PfuseLG.slope.re0*y)
  lines(y,u1/(u1+1),col=c("black"),lwd=2)
  lines(y,u2/(u2+1),col=c("red"))
  lines(y,u3/(u3+1),col=c("purple"))
  lines(y,u4/(u4+1),col=c("green"))
}
AIC(fuseLG.re1,fuseLG.re01,fuseLG.re0)
anova(fuseLG.re1,fuseLG.re0)
anova(fuseLG.re1,fuseLG.re01)

##Intercept only RE with site is best fit matches FE

################Probability of lo clonenum###############
hilo.re1 <- glmer(Philo~size+(size|utrans),data=MCf,family=binomial)
summary(hilo.re1)
Philo.int.re1<-c(fixed.effects(hilo.re1)[1]+random.effects(hilo.re1)[[1]][1:10,1])
Philo.slope.re1<-c(fixed.effects(hilo.re1)[2]+random.effects(hilo.re1)[[1]][1:10,2])

hilo.re01 <- glmer(Philo~size+(size-1|utrans),data=MCf,family=binomial)
summary(hilo.re01)
Philo.int.re01<-fixed.effects(hilo.re01)[1]
Philo.slope.re01<-c(fixed.effects(hilo.re01)[2]+random.effects(hilo.re01)[[1]][1:10,1])

hilo.re0 <- glmer(Philo~size+(1|utrans),data=MCf,family=binomial)
summary(hilo.re0)
Philo.int.re0<-c(fixed.effects(hilo.re0)[1]+random.effects(hilo.re0)[[1]][1:10,1])
Philo.slope.re0<-fixed.effects(hilo.re0)[2]

par(mfrow=c(2,2))
for(i in 1:10){
  x<-subset(MCf,MCf$utrans==levels(MCf$utrans)[i])
  plot(x$size,x$Philo,xlim=c(-2,10),ylim=c(0,1),xlab="Size (t)",ylab="Probability of lo fission #")
  u1<-exp(Philo.int+Philo.slope*y)
  u2<-exp(Philo.int.re1[i]+Philo.slope.re1[i]*y)
  u3<-exp(Philo.int.re01+Philo.slope.re01[i]*y)
  u4<-exp(Philo.int.re0[i]+Philo.slope.re0*y)
  lines(y,u1/(u1+1),col=c("black"),lwd=4)
  lines(y,u2/(u2+1),col=c("red"))
  lines(y,u3/(u3+1),col=c("purple"))
  lines(y,u4/(u4+1),col=c("green"))
}
AIC(hilo.re1,hilo.re01,hilo.re0)
anova(hilo.re1,hilo.re0)
anova(hilo.re1,hilo.re01)
anova(hilo.re0,hilo.re01)
####Of these intercept only is best but the glm without any FE is better.  
##Also, several sites do not have appear to have incidence of high clone number


############################number of clones################3
clonenumhi=glmer(numfiss~size+(size|utrans),data=hinum,family=pospoisson)
summary(clonenumhi)
clonenumhi.int=coefficients(clonenumhi)[1]
clonenumhi.slope=coefficients(clonenumhi)[2]


#####################clonesize######################33
clonesize.re1=lme(sizeNext~size,data=MCfs,random=~size|utrans)
#Iteration error message
clonesize=gls(sizeNext~size,data=MCfs)
summary(clonesize)
clonesize.int=coefficients(clonesize)[1]
clonesize.slope=coefficients(clonesize)[2]
clonesize.sd=summary(clonesize)$sigma

clonesize.re01=lme(sizeNext~size,data=MCfs,random=~size-1|utrans)
summary(clonesize.re01)
clonesize.int.re01=fixed.effects(clonesize.re01)[1]
clonesize.slope.re01=c(fixed.effects(clonesize.re01)[2]+random.effects(clonesize.re01)[1:10,1])
clonesize.sd2.re01=summary(clonesize.re01)$sigma

clonesize.re0=lme(sizeNext~size,data=MCfs,random=~1|utrans)
summary(clonesize.re0)
clonesize.int.re0=c(fixed.effects(clonesize.re0)[1]+random.effects(clonesize.re0)[1:10,1])
clonesize.slope.re0=fixed.effects(clonesize.re0)[2]
clonesize.sd2.re0=summary(clonesize.re0)$sigma

par(mfrow=c(2,2))
for(i in 1:10){
  x<-subset(MCfs,MCfs$utrans==levels(MCfs$utrans)[i])
  plot(x$size,x$sizeNext,xlim=c(-2,10),ylim=c(-2,10), xlab="Size (t)",ylab="Size (t+1)")
  lines(y,clonesize.int+clonesize.slope*y)
  lines(y,clonesize.int.re0[i]+clonesize.slope.re0*y,col=c("red"))
  lines(y,clonesize.int.re01+clonesize.slope.re01[i]*y,col=c("green"))
}
AIC(clonesize,clonesize.re01,clonesize.re0)
anova(clonesize,clonesize.re0)
##RE does not improve model p-va;ue=0.5442 (intercept) p-value=0.1308 (re slope only)
##Still need to deal with variance function


rcrt.lme.UY <- glmer(recruits~1+(1|utrans)+(1|year),family=poisson,data=MCr)
summary(rcrt.lmeU.Y)
recruit.int= c(fixef(rcrt.lmeU.Y)[1],fixef(rcrt.lmeU.Y)[1]+fixef(rcrt.lmeU.Y)[2:nx])

rcrt.lme.Y <- glmer(recruits~1+(1|year),family=poisson,data=MCr)


AIC(rcrt.lme.Y,rcrt.lme.UY)
anova(rcrt.lme.Y,rcrt.lme.UY)

##model with random effect of year and utrans is better


