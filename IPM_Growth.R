#Model growth with linear mixed-effects model in nlme
library(nlme) 
library(lme4)
library(lattice)
MCg<-within(MCg,fate<-relevel(fate,ref="growth"))
dat<-MC #data for this species
growdat<-MCg #growth data for this species
fissdat<-MCf #fission only data for this species
SMfissdat<-MCfs #small fission only data for this species
x<-seq(from=-2,to=10,by=0.01)#Range of values for lines added to plots
nx<-10 # number of sites for this species
nf<-3 #number of fates for this species
lmc<-lmeControl(niter=10000,msMaxIter=10000)

#Alternative transformations
par(mfrow=c(2,2))
hist(growdat$asize,main="Untransformed area")
hist(growdat$size,main="Natural log of area")
hist((growdat$asize)^(1/2),main="Square root area")
hist((growdat$asize)^(1/3),main="cube root area")
plot(growdat$asize,growdat$asizeNext,col=c(1:nx)[growdat$utrans],main="Growth functions by utrans untransformed")
plot(growdat$size,growdat$sizeNext,col=c(1:nx)[growdat$utrans],main="Growth functions by utrans log")
plot((growdat$asize)^(1/2),(growdat$asizeNext)^(1/2),col=c(1:nx)[growdat$utrans],main="Growth functions by utrans square root")
plot((growdat$asize)^(1/3),(growdat$asizeNext)^(1/3),col=c(1:nx)[growdat$utrans],main="Growth functions by utrans cubed root")
legend("bottomright", legend = sprintf(levels(growdat$utrans)), col = c(1:nx), pch = 1, xjust = 1, bg = "white")
#look at normality of residuals with square root transformation
grow.U.fullsqrt<-lm(sqrtsizeNext~sqrtsize*utrans,data=growdat)
plot(grow.U.fullsqrt)
## log normal transformation is still a better option for normalizing data##

###########################################################################################################################
#Linear model of growth using only size as explantory
grow <- lm(sizeNext ~ size,data=growdat)
summary(grow)
#look at residuals
plot(grow)
##Heterogeneity in residuals: decreases with size##

###########################################################################################################################
#Look at the effect of site = utrans
grow.lisU <- lmList(sizeNext~size|utrans,data=growdat)
summary(grow.lisU)
plot(intervals(grow.lisU))
xyplot(sizeNext~size|utrans,data=growdat,xlab="Log(Area)[t]",ylab="Log(Area)[t+1]",main="Growth Function by site")

grow.U<-lm(sizeNext~size*utrans,data=growdat)
plot(grow.U)

par(mfrow=c(1,1))
cl<-rainbow(nx)
plot(growdat$size,growdat$sizeNext,col=cl[growdat$utrans],main="Growth functions by utrans")
legend("bottomright", legend = sprintf(levels(growdat$utrans)), col = cl, pch = 1, xjust = 1, bg = "white")
lines(c(0,7),c(0,7),lty=2)
b<-coefficients(grow.U)[1]+x*coefficients(grow.U)[2]
lines(x,b,col=cl[1])
for (i in 1:(nx-1)) {
  b <- coefficients(grow.U)[1]+coefficients(grow.U)[i+2]+x*(coefficients(grow.U)[2]+coefficients(grow.U)[i+nx+1])
  lines(x,b,col=cl[i+1])
}

###########################################################################################################################             
#Should we include fate as a fixed effect (does fission or fusion change growth intercept or slope)?
grow.F<-lm(sizeNext~size*fate,data=growdat)
plot(growdat$size,growdat$sizeNext,col=cl[growdat$fate],main="Growth functions by utrans")
legend("bottomright", legend = sprintf(levels(growdat$fate)), col = cl, pch = 1, xjust = 1, bg = "white")
lines(c(0,7),c(0,7),lty=2)
b<-coefficients(grow.F)[1]+x*coefficients(grow.F)[2]
lines(x,b,col=cl[1])
for (i in 1:(nf-1)) {
  b <- coefficients(grow.F)[1]+coefficients(grow.F)[i+2]+x*(coefficients(grow.F)[2]+coefficients(grow.F)[i+4])
  lines(x,b,col=cl[i+1])
}
#are patterns associated with fate similar between sites?
xyplot(sizeNext~size|utrans,data=growdat,groups=growdat$fate,auto.key=TRUE,xlab="Log(Area)[t]",ylab="Log(Area)[t+1]",main="Growth Function all")
# Do we have enough fission and fusion transitions to estimate seperate parameters for each growth fate (fission, fusion, grow)?
xtabs(~utrans+fate,data=growdat)

###########################################################################################################################             
#Looks like we could use site and fate as fixed effects.  There is some heterogenity in the residuals.
#What if we add a varFunc? See Zurr chapter 4 protocol on page 90
#Select varFunc using fullest model that will converge = no 3way interaction
##w/ 3-way interaction get error in estimate : computed "gls" fit is singular, rank 48##
##removing HonoluaS and Papaula 10m (lowest ns) does not improve convergence##
grow.UF.nov<-gls(sizeNext~size*utrans+size*fate+utrans*fate,data=growdat)
grow.UF.vps<-gls(sizeNext~size*utrans+size*fate+utrans*fate,data=growdat,weights=varPower(form=~size))
grow.UF.ves<-gls(sizeNext~size*utrans+size*fate+utrans*fate,data=growdat,weights=varExp(form=~size))
anova(grow.UF.vps,grow.UF.ves)
##varExp better fit than varPower by delta=15.849 ##
anova(grow.UF.nov,grow.UF.ves)
##varExp is significanly better fit (L=39.04758,p=<0.001,delta=37.048) than no varFunc
grow.UF.vesu<-gls(sizeNext~size*utrans+size*fate+utrans*fate,data=growdat,weights=varExp(form=~size|utrans))
anova(grow.UF.ves,grow.UF.vesu)
##including utrans in varFunc improves fit (L=35.7443,p-value=<0.0001,delta=17.744)
grow.UF.vesf<-gls(sizeNext~size*utrans+size*fate+utrans*fate,data=growdat,weights=varExp(form=~size|fate))
anova(grow.UF.ves,grow.UF.vesf)
##including fate in varFunc improves fit over size alone (L=21.15661,p-value=<0.0001,delta=17.156) 
anova(grow.UF.vesu,grow.UF.vesf)
##replacing utrans with fate in varFunc does not improve the model = AIC for utrans lower (L=14.58769,p-value=0.0417,delta=0.588)
grow.UF.vesuf<-gls(sizeNext~size*utrans+size*fate+utrans*fate,data=growdat,weights=varExp(form=~size|utrans*fate))
anova(grow.UF.vesu,grow.UF.vesuf)
##including fate with utrans in varFunc improves the models (L=56.34169,p-value=<0.0001,delta=16.342)


#Select Fixed Effects using varExp(form=~size|utrans)
grow.UF.full<-gls(sizeNext~size*utrans+size*fate+utrans*fate,data=growdat,method="ML",weights=varExp(form=~size|utrans))
grow.UF.drop1<-gls(sizeNext~size*utrans+size*fate,data=growdat,method="ML",weights=varExp(form=~size|utrans))
anova(grow.UF.full,grow.UF.drop1)
#utrans*fate interaction does not improve fit (L= 16.21699, p-value= 0.5774, delta=19.783)
grow.UF.drop2<-gls(sizeNext~size*utrans+fate,data=growdat,method="ML",weights=varExp(form=~size|utrans))
anova(grow.UF.drop1,grow.UF.drop2)
#keep size:fate (L=8.849814, p-value=0.012, delta=5.415)
grow.UF.drop3<-gls(sizeNext~utrans+size*fate,data=growdat,method="ML",weights=varExp(form=~size|utrans))
anova(grow.UF.drop1,grow.UF.drop3)
##keep size:utrans (L=21.18752,p-value=0.0058, delta=5.187)

#refit best fixed effects model with REML
grow.UF.vesu<-gls(sizeNext~size*utrans+size*fate,data=growdat,weights=varExp(form=~size|utrans))

#look at normalized residuals to see if varFunc makes a difference
Esu<-resid(grow.UF.vesu,type="normalized")
coplot(Esu~size|utrans,data=growdat,ylab="normalized residuals")
Esf<-resid(grow.UF.vesf,taype="normalized")
coplot(Esf~size|utrans,data=growdat,ylab="normalized residuals")
###########################################################################################################################             
#Now, look at effect of year
grow.lisY <-lmList(sizeNext~size|fyear,growdat)
summary(grow.lisY)
xyplot(sizeNext~size|fyear,data=growdat,xlab="Log(Area)[t]",ylab="Log(Area)[t+1]",main="Growth Function all")
plot(intervals(grow.lisY))

#year as random
grow.UF.vesu<-gls(sizeNext~size*utrans+size*fate+utrans*fate,weights=varExp(form=~size|utrans),data=growdat,control=lmc)
grow.UF.vesu.Y <- lme(sizeNext~size*utrans+size*fate+utrans*fate,random=~size|fyear,weights=varExp(form=~size|utrans),data=growdat,control=lmc)
grow.UF.vesu.Y0 <- lme(sizeNext~size*utrans+size*fate+utrans*fate, random=~1|fyear,weights=varExp(form=~size|utrans),data=growdat,control=lmc)
grow.UF.vesu.Y1 <- lme(sizeNext~size*utrans+size*fate+utrans*fate, random=~size-1|fyear,weights=varExp(form=~size|utrans),data=growdat,control=lmc)
anova(grow.UF.vesu,grow.UF.vesu.Y)
anova(grow.UF.vesu,grow.UF.vesu.Y0)
anova(grow.UF.vesu,grow.UF.vesu.Y1)
##Year as random effect is not significantly different from fixed effects model##


###########################################################################################################################
##probability of fissing, fusing or growing given that I survive.  See survival analysis for probability of survivial##
fuse.full<-glm(Pfuse~size*utrans+size*fyear+utrans*fyear,data=dat,family=binomial)
drop1(fuse.full,test="Chi")
##drop size*fyear##
fuse.UY<-glm(Pfuse~size*utrans+utrans*fyear,data=dat,family=binomial)
drop1(fuse.UY,test="Chi")
##keep size*utrans and utrans*fyear##
#We wouldn't really keep year as FE so try just size and utrans
fuse.U<-glm(Pfuse~size*utrans,data=dat,family=binomial)
drop1(fuse.U,test="Chi")
##keep size*utrans interaction (p-value = 0.007868)## 
fuse.U0<-glm(Pfuse~size+utrans,data=dat,family=binomial)
drop1(fuse.U0,test="Chi")


#check fyear random effect
fuse.U.Q0Y0 <- glmer(Pfuse~size*utrans +(1|uyear),data=dat,family=binomial)
fuse.U.Q0 <- glmer(Pfuse~size*utrans +(1|uquad),data=dat,family=binomial)
fuse.U.Y0 <- glmer(Pfuse~size*utrans+(1|fyear),data=dat,family=binomial)
#compare 
anova(fuse.U.Q0Y0,fuse.U.Q0,fuse.U.Y0)
##Q0Y0 and Q0 improve AIC over fuse.U.  Q0Y0 is best but the concatenation does not make very much sense##
##Do not know if I can really compare FE only and RE models but AIC of random effect of year is higher by delta = 63 so we won't use it##

growth.full<-glm(Pgrowth~size*utrans+size*fyear+utrans*fyear,data=dat,family=binomial)
drop1(growth.full,test="Chi")
growth.UY<-glm(Pgrowth~size*utrans+size*fyear,data=dat,family=binomial)
drop1(growth.UY,test="Chi")
growth.UY0<-glm(Pgrowth~size*utrans+fyear,data=dat,family=binomial)
drop1(growth.UY0,test="Chi")
growth.U<-glm(Pgrowth~size*utrans,data=dat,family=binomial)
drop1(growth.U,test="Chi")
##Keep size*utrans interaction (p-value=0.003809)


fiss.full<-glm(Pfiss~size*utrans+size*fyear+utrans*fyear,data=dat,family=binomial)
drop1(fiss.full,test="Chi")
fiss.UY0<-glm(Pfiss~size*utrans+utrans*fyear,data=dat,family=binomial)
drop1(fiss.UY0,test="Chi")
fiss.U0Y0<-glm(Pfiss~size+utrans*fyear,data=dat,family=binomial)
drop1(fiss.U0Y0,test="Chi")
##keep utrans*fyear but we wouldn't use this as FE so try without year##

fiss.U<-glm(Pfiss~size*utrans,data=dat,family=binomial)
drop1(fiss.U,test="Chi")
fiss.U0<-glm(Pfiss~size+utrans,data=dat,family=binomial)
drop1(fiss.U0,test="Chi")
##keep utrans main effect##

#check fyear random effect
fiss.U.Q0Y0 <- glmer(Pfiss~size*utrans +(1|uyear),data=dat,family=binomial)
fiss.U.Q0 <- glmer(Pfiss~size*utrans +(1|uquad),data=dat,family=binomial)
fiss.U.Y0 <- glmer(Pfiss~size*utrans+(1|fyear),data=dat,family=binomial)
#compare 
anova(fiss.U.Q0Y0,fiss.U.Q0,fiss.U.Y0)
##Q0Y0 and Q0 improve AIC over fuse.U.  Q0Y0 is best but the concatenation does not make very much sense##
##Do not know if I can really compare FE only and RE models but AIC of random effect of year is higher by delta = 35 so we won't use it##

###########################################################################################################################
#probability of fuseLG given survival and fusion
fuseLG.full<-glm(PfuseLG~size*utrans+size*fyear+utrans*fyear,data=dat,family=binomial)
drop1(fuseLG.full,test="Chi")
fuseLG.UY<-glm(PfuseLG~size*utrans+size*fyear,data=dat,family=binomial)
drop1(fuseLG.UY,test="Chi")
fuseLG.U0Y<-glm(PfuseLG~utrans+size*fyear,data=dat,family=binomial)
drop1(fuseLG.U0Y,test="Chi")
fuseLG.U0Y0<-glm(PfuseLG~utrans+size+fyear,data=dat,family=binomial)
drop1(fuseLG.U0Y0,test="Chi")
fuseLG.U0<-glm(PfuseLG~size+utrans,data=dat,family=binomial)
drop1(fuseLG.U0,test="Chi")
##keep utrans main effect