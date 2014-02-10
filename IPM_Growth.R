#Model growth with linear mixed-effects model in nlme
library(nlme) 
library(lme4)
library(lattice)
MCg<-within(MCg,fate<-relevel(fate,ref="growth"))
dat<-MC #data for this species
growdat<-MCg #growth data for this species
#fissdat<-MCf #fission only data for this species
#SMfissdat<-MCfs #small fission only data for this species
x<-seq(from=-2,to=10,by=0.01)#Range of values for lines added to plots
nx<-10 # number of sites for this species
nf<-3 #number of fates for this species
lmc<-lmeControl(niter=10000,msMaxIter=10000)# use to test year as RE

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
##Considered removing HonS and Papaula 10m because of low n.  Removing HonoluaS and Papaula 10m (lowest ns) does not improve convergence##
#Variance component by size: no variance, Power structure, exponential, constant power 
grow.UF.nov<-gls(sizeNext~size*utrans+size*fate+utrans*fate,data=growdat)
grow.UF.vps<-gls(sizeNext~size*utrans+size*fate+utrans*fate,data=growdat,weights=varPower(form=~size))
grow.UF.ves<-gls(sizeNext~size*utrans+size*fate+utrans*fate,data=growdat,weights=varExp(form=~size))
grow.UF.vcps<-gls(sizeNext~size*utrans+size*fate+utrans*fate,data=growdat,weights=varConstPower(form=~size))
#Variance component by size and transect: no variance, Power structure, exponential, constant power, combination (similar to varExp by utrans.  Combination of varIdent and varExp)
grow.UF.vesu<-gls(sizeNext~size*utrans+size*fate+utrans*fate,data=growdat,weights=varExp(form=~size|utrans))
grow.UF.vpsu<-gls(sizeNext~size*utrans+size*fate+utrans*fate,data=growdat,weights=varPower(form=~size|utrans))
grow.UF.vcpsu<-gls(sizeNext~size*utrans+size*fate+utrans*fate,data=growdat,weights=varConstPower(form=~size|utrans))
grow.UF.vciesu<-gls(sizeNext~size*utrans+size*fate+utrans*fate,data=growdat,weights=varComb(varIdent(form=~1|utrans),varExp(form=~size)))
grow.UF.vcipsu<-gls(sizeNext~size*utrans+size*fate+utrans*fate,data=growdat,weights=varComb(varIdent(form=~1|utrans),varPower(form=~size)))
#Compare model AIC
AIC(grow.UF.nov,grow.UF.vps,grow.UF.ves,grow.UF.vcps,grow.UF.vesu,grow.UF.vpsu,grow.UF.vcpsu,grow.UF.vciesu,grow.UF.vcipsu)
##by AIC varExp is best model for size only.  varPower was best model for size|utrans.  All models with utrans included had lower AIC 
##varComb and varExp including utrans should be really similar but they are not.  varExp is a better fit so we will go with that one.

#Compare residuals
E<-resid(grow.UF.nov,type="normalized")
coplot(E~size|utrans,data=growdat,ylab="normalized residuals")
E<-resid(grow.UF.vps,type="normalized")
coplot(E~size|utrans,data=growdat,ylab="normalized residuals")
E<-resid(grow.UF.ves,type="normalized")
coplot(E~size|utrans,data=growdat,ylab="normalized residuals")
E<-resid(grow.UF.vcps,type="normalized")
coplot(E~size|utrans,data=growdat,ylab="normalized residuals")
E<-resid(grow.UF.vesu,type="normalized")
coplot(E~size|utrans,data=growdat,ylab="normalized residuals")
E<-resid(grow.UF.vpsu,type="normalized")
coplot(E~size|utrans,data=growdat,ylab="normalized residuals")
E<-resid(grow.UF.vcpsu,type="normalized")
coplot(E~size|utrans,data=growdat,ylab="normalized residuals")
E<-resid(grow.UF.vciesu,type="normalized")
coplot(E~size|utrans,data=growdat,ylab="normalized residuals")
E<-resid(grow.UF.vcipsu,type="normalized")
coplot(E~size|utrans,data=growdat,ylab="normalized residuals")

#Plot variance by size
growth.sd2=((summary(grow.UF.vps)$sigma))^2
growth.varFunc=coef(grow.UF.vps$modelStruct$varStruct,"exp")
sigmaxp2<-growth.sd2*abs(x)^(2*growth.varFunc)
plot(sigmaxp2~x)

growth.sd2=((summary(grow.UF.ves)$sigma))^2
growth.varFunc=coef(grow.UF.ves$modelStruct$varStruct,"exp")
sigmaxp2<-growth.sd2*exp(2*growth.varFunc*x)
plot(sigmaxp2~x)

growth.sd2=((summary(grow.UF.vcps)$sigma))^2
growth.varFunc=coef(grow.UF.vcps$modelStruct$varStruct)
sigmaxp2<-growth.sd2*(growth.varFunc[1]+abs(x)^(growth.varFunc[2])^2)
plot(sigmaxp2~x)

par(mfrow=c(3,4))
growth.sd2=((summary(grow.UF.vesu)$sigma))^2
growth.varFunc=coef(grow.UF.vesu$modelStruct$varStruct,"exp")
for(i in 1:10){
  sigmaxp2<-growth.sd2*exp(2*growth.varFunc[i]*x)
  plot(sigmaxp2~x)}

par(mfrow=c(3,4))
growth.sd2=((summary(grow.UF.vpsu)$sigma))^2
growth.varFunc=coef(grow.UF.vpsu$modelStruct$varStruct,"exp")
for(i in 1:10){
  sigmaxp2<-growth.sd2*abs(x)^(2*growth.varFunc[i])
  plot(sigmaxp2~x)
}

par(mfrow=c(3,4))
growth.sd2=((summary(grow.UF.vcpsu)$sigma))^2
growth.varFunc=coef(grow.UF.vcpsu$modelStruct$varStruct)
for (i in 2:11){
  sigmaxp2<-growth.sd2*(growth.varFunc[1]+abs(x)^(growth.varFunc[i])^2)
  plot(sigmaxp2~x)
  
}
## varExp is the only model that makes sense for our data.  The model is a better fit when utrans is included in variance structure.

#Select Fixed Effects using varExp(form=~size|utrans)
grow.UF.full<-gls(sizeNext~size*utrans+size*fate+utrans*fate,data=growdat,method="ML",weights=varExp(form=~size|utrans))
grow.UF.drop1<-gls(sizeNext~size*utrans+size*fate,data=growdat,method="ML",weights=varExp(form=~size|utrans))
anova(grow.UF.full,grow.UF.drop1)
#drop utrans*fate interaction. does not improve fit (L= 16.70579, p-value= 0.5434, delta=19.294)
grow.UF.full2<-gls(sizeNext~size*utrans+size*fate,data=growdat,method="ML",weights=varExp(form=~size|utrans))
grow.UF.drop2<-gls(sizeNext~size*utrans+fate,data=growdat,method="ML",weights=varExp(form=~size|utrans))
anova(grow.UF.full2,grow.UF.drop2)
#keep size:fate (L=10.63786, p-value=0.0049, delta=6.638)
grow.UF.drop3<-gls(sizeNext~utrans+size*fate,data=growdat,method="ML",weights=varExp(form=~size|utrans))
anova(grow.UF.full2,grow.UF.drop3)
##keep size:utrans (L=24.0147,p-value=0.0043, delta=6.015)

#refit best fixed effects model with REML
grow.UF.vesu<-gls(sizeNext~size*utrans+size*fate,data=growdat,weights=varExp(form=~size|utrans))
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
#No way to directly compare RE model with a model with no RE so we fit year as a fixed effect to see if it makes a difference
fuse.full<-glm(Pfuse~size*utrans+size*fyear+utrans*fyear,data=dat,family=binomial)
drop1(fuse.full,test="Chi")
##drop utrans*fyear##
fuse.UY<-glm(Pfuse~size*utrans+size*fyear,data=dat,family=binomial)
drop1(fuse.UY,test="Chi")
##keep size*utrans and size*fyear##

#We wouldn't really keep year as FE.  Year might make a difference as a RE.
#Here we look at FE only = size and utrans
fuse.U<-glm(Pfuse~size*utrans,data=dat,family=binomial)
drop1(fuse.U,test="Chi")
##drop size*utrans(p-value=0.1646) 
fuse.U0<-glm(Pfuse~size+utrans,data=dat,family=binomial)
drop1(fuse.U0,test="Chi")
#Keep utrans main effect (p-value=2e-16)


#Compare RE effects models: quad and year, quad only, year only
fuse.U.Q0Y0 <- glmer(Pfuse~size*utrans +(1|uyear),data=dat,family=binomial)
fuse.U.Q0 <- glmer(Pfuse~size*utrans +(1|uquad),data=dat,family=binomial)
fuse.U.Y0 <- glmer(Pfuse~size*utrans+(1|fyear),data=dat,family=binomial)
#compare 
anova(fuse.U.Q0Y0,fuse.U.Q0,fuse.U.Y0)
##Q0Y0 and Q0 improve AIC over fuse.U.  Q0Y0 is best but the concatenation does not make very much sense##
##Do not know if I can really compare FE only and RE models but AIC of random effect of year is higher by delta = 7.1 so we won't use it##

growth.full<-glm(Pgrowth~size*utrans+size*fyear+utrans*fyear,data=dat,family=binomial)
drop1(growth.full,test="Chi")
#drop size*fyear p-value=0.0547562
growth.UY<-glm(Pgrowth~utrans*fyear+size*utrans,data=dat,family=binomial)
drop1(growth.UY,test="Chi")
#keep utrans*size and utrans*fyear

#We wouldn't really keep year as FE.  Year might make a difference as a RE.
#Here we look at FE only = size and utrans
growth.full<-glm(Pgrowth~size*utrans,data=dat,family=binomial)
drop1(growth.full,test="Chi")
#keep size*utrans p-value=0.003599

#Compare RE effects models: quad and year, quad only, year only
growth.U.Q0Y0 <- glmer(Pgrowth~size*utrans +(1|uyear),data=dat,family=binomial)
growth.U.Q0 <- glmer(Pgrowth~size*utrans +(1|uquad),data=dat,family=binomial)
growth.U.Y0 <- glmer(Pgrowth~size*utrans+(1|fyear),data=dat,family=binomial)
#compare 
anova(growth.U.Q0Y0,growth.U.Q0,growth.U.Y0)
##Q0Y0 and Q0 improve AIC over fuse.U.  Q0Y0 is best but the concatenation does not make very much sense##
##Do not know if I can really compare FE only and RE models but AIC of random effect of year is higher by delta = 42 so we won't use it##

fiss.full<-glm(Pfiss~size*utrans+size*fyear+utrans*fyear,data=dat,family=binomial)
drop1(fiss.full,test="Chi")
#drop size*fyear p-value=0.0547562
fiss.UY0<-glm(Pfiss~size*utrans+utrans*fyear,data=dat,family=binomial)
drop1(fiss.UY0,test="Chi")
#keep utrans*size and utrans*fyear

#We wouldn't really keep year as FE.  Year might make a difference as a RE.
#Here we look at FE only = size and utrans
fiss.U<-glm(Pfiss~size*utrans,data=dat,family=binomial)
drop1(fiss.U,test="Chi")
#keep size*utrans p-value=0.003599


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
#drop utrans*fyear p-value=0.8650
fuseLG.UY<-glm(PfuseLG~size*utrans+size*fyear,data=dat,family=binomial)
drop1(fuseLG.UY,test="Chi")
#drop size*utrans p-value=0.6107
fuseLG.U0Y<-glm(PfuseLG~utrans+size*fyear,data=dat,family=binomial)
drop1(fuseLG.U0Y,test="Chi")
#drop size*fyear p-value=0.4139
fuseLG.U0Y0<-glm(PfuseLG~utrans+size+fyear,data=dat,family=binomial)
drop1(fuseLG.U0Y0,test="Chi")
#drop fyear p-value=0.2785

##year does not matter so do not need to test for RE
#Here we look at FE only = size and utrans
fuseLG.full<-glm(PfuseLG~size*utrans,data=dat,family=binomial)
drop1(fuseLG.full,test="Chi")
#drop size*utrans p-value=0.7149
fuseLG.U0<-glm(PfuseLG~size+utrans,data=dat,family=binomial)
drop1(fuseLG.U0,test="Chi")
##keep size p-value=2.2e-16 and utrans p-value=0.0001138 main effects