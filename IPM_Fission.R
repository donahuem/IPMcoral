#Model reproduction from fission using mixed-effects logistic and poisson regression
library(lme4) 
library(nlme) 
library(lattice)
library(VGAM)
dat<-MC
fissdat<-MCf #fission only data for this species
SMfissdat<-MCfs #small fission only data for this species
nx<-10 # number of sites for this species
x <-seq(min(dat$size,na.rm=TRUE),max(dat$size,na.rm=TRUE),.1)
#What are we modeling here?
#1. Probability of fission (Ft+1) as a function of size (At) (color by fate_t1), by location

#2. Probability distribution of NUMBER of fission products GIVEN area of colony fissed (area_prefiss) by location

#3. Size distribution of fission products

#####################################################################################
#(1)  Test for effect of size on fission probability
#What are we modeling?
plot(dat$size,dat$Pfiss,col=c(1:6)[dat$fateIPM])
legend("right",legend=levels(dat$fateIPM),col=1:6,pch=1)
xyplot(dat$Pfiss~dat$size|dat$utrans, group=dat$fateIPM,auto.key=list(columns=2),xlab="Log(Area)[t]",ylab="Fission",main=list(label="M. capitata"))
histogram(~dat$Pfiss|dat$utrans,type="count")



#first, just look at FE model
fiss.U <- glm(Pfiss~csize*utrans,data=dat,family=binomial)
drop1(fiss.U,test="Chi")
fiss.U0 <- glm(Pfiss~csize+utrans,data=dat,family=binomial)
drop1(fiss.U0,test="Chi")
##Keep utrans main effect

#Next, add RE for quad and year  #since utrans is fixed effect, let year RE be simply the mean effect of year on the intercept (fyear)
fiss.U.Q0Y0 <- glmer(Pfiss~csize*utrans +(1|uyear),data=dat,family=binomial)
fiss.U.Y0 <- glmer(Pfiss~csize*utrans +(1|fyear),data=dat,family=binomial)
fiss.U.Q0 <- glmer(Pfiss~csize*utrans +(1|uquad),data=dat,family=binomial)
#compare 
AIC(fiss.U,fiss.U.Q0Y0,fiss.U.Y0,fiss.U.Q0)
#fiss.U.Q0Y0 is best by AIC

#Select FE
drop1(fiss.U.Q0Y0)
fiss.U0.Q0Y0<-glmer(Pfiss~csize+utrans+(1|uyear),data=dat,family=binomial)
drop1(fiss.U0.Q0Y0)
##keep utrans main effect


#Plot the fission probability function & data

b <-attr(fiss.U0.Q0Y0,"fixef")[1]+x*attr(fiss.U0.Q0Y0,"fixef")[2]
pfx <- exp(b)/(1+exp(b))
cl<-rainbow(nx)
plot(x,pfx,col=cl[1],type="l",ylim=c(0,1),xlab="Log(Area)[t]", ylab="Prob(Fission)")
#title("P(fission) ~ Logit(Log(Area)+Transect)")
points(dat$size[dat$utrans==unique(dat$utrans)[1]],dat$Pfiss[dat$utrans==unique(dat$utrans)[1]],col=cl[1])
for (i in 1:(nx-1)) {
  b <- attr(fiss.U0.Q0Y0,"fixef")[1]+attr(fiss.U0.Q0Y0,"fixef")[i+2]+x*attr(fiss.U0.Q0Y0,"fixef")[2]
  pfx <- exp(b)/(1+exp(b))
  lines(x,pfx,col=cl[i+1])
  points(dat$size[dat$utrans==unique(dat$utrans)[i+1]],dat$Pfiss[dat$utrans==unique(dat$utrans)[i+1]],col=cl[i+1])
  }
legend("left",legend=unique(dat$utrans),col=cl[1:nx],lty=1,pt.cex=0.3,bty="n")


####################################################################################
#2. Number of Fission Products Given Area_prefiss
#plot the problem\
hist(MC$numfiss)
histogram(~MC$numfiss|MC$utrans,type="count",xlab="Number of Fission Products", ylab="Frequency")
xyplot(MC$numfiss~MC$size|MC$utrans,xlab="Log(Area) before Fission", ylab="Number of Fission Products",main=list(label="M. capitata"))

#First, I fit with Poisson [not correct -- need zero-truncated]
nfiss.UY <- glm(numfiss~csize*utrans+csize*fyear,data=fissdat,family=poisson)
drop1(nfiss.UY)
nfiss.U0Y <- glm(numfiss~utrans+csize*fyear,data=fissdat,family=poisson)
drop1(nfiss.U0Y)
nfiss.Y<- glm(numfiss~csize*fyear,data=fissdat,family=poisson)
drop1(nfiss.Y)
nfiss.Y0<-glm(numfiss~csize+fyear,data=fissdat,family=poisson)
drop1(nfiss.Y0)
nfiss<-glm(numfiss~csize,data=fissdat,family=poisson)
drop1(nfiss)

##Keep size

nfiss.S.Y0 <- glmer(numfiss~csize*utrans+(1|fyear),data=fissdat,family=poisson)
nfiss.glmeU_A.Y0 <- glmer(numfiss~csize+utrans+(1|fyear),data=fissdat,family=poisson)
nfiss.glmeA.Y0 <- glmer(numfiss~csize+(1|fyear),data=fissdat,family=poisson)
anova(nfiss.S.Y0,nfiss.glmeA.Y0,nfiss.glmeU_A.Y0)
#the variance estimates assoc with year are *tiny* and the deviance is 65.7 -not much better
#CONCLUSION: don't bother with REs for year
#also checked FE of year in ztPoisson - no improvement

#ZERO-TRUNCATED POISSON: Next, I realized that I really need a zero-truncated poisson bc there is no possibility of 0 fission products 
# given the structure of the data

nfiss.UY <- vglm(numfiss~utrans*csize+fyear,data=fissdat,family=pospoisson)
AIC(nfiss.UY)
#AIC=346.0597
nfiss.U <- vglm(numfiss~utrans*csize,data=fissdat,family=pospoisson)
AIC(nfiss.U)
#AIC=350.819
##fyear improves fit by AIC by delta 4.7593
nfiss.U0<- vglm(numfiss~utrans+csize,data=fissdat,family=pospoisson)
AIC(nfiss.U0)
#AIC=357.8996
##keep size*utrans interaction

#Refit best model uncentered
nfiss.U <- vglm(numfiss~utrans*size,data=fissdat,family=pospoisson)

#now, do some model evaluation
par(mfrow=c(1,1))
plot(resid(nfiss.U,type="pearson")~fissdat$csize)
plot(resid(nfiss.U,type="pearson")~fissdat$utrans)
plot(resid(nfiss.U,type="pearson")~fissdat$fyear)
#there are a few outliers, but nothing systematic

cl<-rainbow(nx)
temp_lap <- rep(seq((min(fissdat$size)-.1),(max(fissdat$size)+.1),.05),nx)
temp_ut<-{}
for (i in 1:nx) {
  temp_ut <- c(temp_ut,rep(unique(fissdat$utrans)[i],length(seq((min(fissdat$size)-.1),(max(fissdat$size)+.1),.05))))
}
temp_u <- unique(fissdat$utrans)[temp_ut]
MyData<-data.frame(size=temp_lap,utrans=temp_u)
G <- predictvglm(nfiss.U,newdata=MyData,data=fissdat,type="link",se=TRUE)#with "link" can get s.e. bands; need to exp(G) + 1 for predictions
plot(fissdat$size[fissdat$utrans==unique(fissdat$utrans)[1]],fissdat$numfiss[fissdat$utrans==unique(fissdat$utrans)[1]],col=cl[1],xlab="Log(Area)[t]",ylim=c(0.5,10.5),xlim=c(-1,7),ylab="Number of Fission Products")
lines(MyData$size[MyData$utrans==unique(MyData$utrans)[1]],exp(G$fit[MyData$utrans==unique(MyData$utrans)[1]])+1,lty=1,col=cl[1])
#title("Number of Fission Products ~ Zero-Trunc-Poisson(Log(Area) + Transect)")
for (i in 2:nx) {
  points(fissdat$size[fissdat$utrans==unique(fissdat$utrans)[i]],fissdat$numfiss[fissdat$utrans==unique(fissdat$utrans)[i]],col=cl[i])
  lines(MyData$size[MyData$utrans==unique(MyData$utrans)[i]],exp(G$fit[MyData$utrans==unique(MyData$utrans)[i]])+1,lty=1,col=cl[i])
}
legend("topleft",legend=unique(fissdat$utrans),col=cl[1:nx],lty=1,pt.cex=0.3,bty="n")


############################################################################
#3. Size distribution of fission products
#This could be a function of numfiss or size= size of parent colony or sumfiss


hist(SMfissdat$sizeNext)
histogram(~SMfissdat$sizeNext|SMfissdat$utrans,xlab="Log(Area) of Fission Products", main=list(label="Log(Area) of Fission Products ~ N(mu_trans, var_trans)"))
histogram(~SMfissdat$sizeNext|SMfissdat$fyear,xlab="Log(Area) of Fission Products", main=list(label="Log(Area) of Fission Products ~ N(mu_trans, var_trans)"))
plot(SMfissdat$size,SMfissdat$sizeNext)
plot(SMfissdat$fec2,SMfissdat$sizeNext)
plot(SMfissdat$fec2,SMfissdat$size)

#look at plot, it's clear that mean and var will vary by utrans.
#no effect of year
sfiss.UY<- gls(sizeNext~size*utrans+size*fyear,data=SMfissdat,method="ML")
sfiss.UY0<- gls(sizeNext~size*utrans+fyear,data=SMfissdat,method="ML")
anova(sfiss.UY,sfiss.UY0)
sfiss.U<- gls(sizeNext~size*utrans,data=SMfissdat,method="ML")
anova(sfiss.UY0,sfiss.U)
sfiss.U0<- gls(sizeNext~size+utrans,data=SMfissdat,method="ML")
anova(sfiss.U,sfiss.U0)
##drop size*utrans interaction
sfiss<- gls(sizeNext~size,data=SMfissdat,method="ML")
anova(sfiss.U0,sfiss)
##drop utrans

sfiss.U<-gls(sizeNext~size*utrans,data=SMfissdat) 
sfiss.U.viu <- gls(sizeNext~size*utrans,data=SMfissdat,weights=varIdent(form=~1|utrans))
sfiss.U.vesu<-gls(sizeNext~size*utrans,data=SMfissdat,weights=varExp(form=~size|utrans))
sfiss.U.vpsu<-gls(sizeNext~size*utrans,data=SMfissdat,weights=varPower(form=~size|utrans))
sfiss.U.ves<-gls(sizeNext~size*utrans,data=SMfissdat,weights=varExp(form=~size))
anova(sfiss.U,sfiss.U.viu)
anova(sfiss.U,sfiss.U.vesu)
anova(sfiss.U,sfiss.U.vpsu)
anova(sfiss.U.vesu,sfiss.U.ves)
##varExp better by AIC.  size alone better than size|utrans
sfiss.U.ves<-gls(sizeNext~size*utrans,method="ML",data=SMfissdat,weights=varExp(form=~size))
sfiss.U0.ves<-gls(sizeNext~size+utrans,method="ML",data=SMfissdat,weights=varExp(form=~size))
anova(sfiss.U.ves,sfiss.U0.ves)
sfiss.ves<-gls(sizeNext~size,method="ML",data=SMfissdat,weights=varExp(form=~size))
anova(sfiss.U0.ves,sfiss.ves)
sfiss.1.ves<-gls(sizeNext~1,method="ML",data=SMfissdat,weights=varExp(form=~size))
anova(sfiss.ves,sfiss.1.ves)
##do not use utrans keep size##

##refit best model with reml##
sfiss.ves<-gls(sizeNext~size,data=SMfissdat,weights=varExp(form=~size))
sfiss<-gls(sizeNext~size,data=SMfissdat)
Es<-resid(sfiss.ves,type="normalized")
coplot(Es~SMfissdat$size|SMfissdat$utrans,ylab="normalized residuals")

E<-resid(sfiss,type="normalized")
coplot(E~SMfissdat$size|SMfissdat$utrans,ylab="normalized residuals")