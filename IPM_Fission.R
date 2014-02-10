#Model reproduction from fission using mixed-effects logistic and poisson regression
library(lme4) 
library(nlme) 
library(lattice)
library(VGAM)#for zero truncated poisson distribution
library(MASS)# for negative binomial distribution
library(countreg) # for zerotrunc for geometric distribution
library(nnet) # for multinomial distribution
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
#(1)  Test for effect of size on fission probability (See IPM_growth)
#What are we modeling?
plot(dat$size,dat$Pfiss,col=c(1:6)[dat$fateIPM])
legend("right",legend=levels(dat$fateIPM),col=1:6,pch=1)
xyplot(dat$Pfiss~dat$size|dat$utrans, group=dat$fateIPM,auto.key=list(columns=2),xlab="Log(Area)[t]",ylab="Fission",main=list(label="M. capitata"))
histogram(~dat$Pfiss|dat$utrans,type="count")



#first, just look at FE model See growth
fiss.U <- glm(Pfiss~size*utrans,data=dat,family=binomial)
drop1(fiss.U,test="Chi")
#keep size*utrans p-value=0.003599

#Next, add RE for quad and year  #since utrans is fixed effect, let year RE be simply the mean effect of year on the intercept (fyear)
fiss.U.Q0Y0 <- glmer(Pfiss~size*utrans +(1|uyear),data=dat,family=binomial)
fiss.U.Y0 <- glmer(Pfiss~size*utrans +(1|fyear),data=dat,family=binomial)
fiss.U.Q0 <- glmer(Pfiss~size*utrans +(1|uquad),data=dat,family=binomial)
#compare 
AIC(fiss.U,fiss.U.Q0Y0,fiss.U.Y0,fiss.U.Q0)
#fiss.U.Q0Y0 is best by AIC
#year as RE does not improve fit

#Plot the fission probability function & data


####################################################################################
#2. Number of Fission Products Given Area_prefiss
#plot the problem
hist(MC$numfiss)
histogram(~MC$numfiss|MC$utrans,type="count",xlab="Number of Fission Products", ylab="Frequency")
xyplot(MC$numfiss~MC$size|MC$utrans,xlab="Log(Area) before Fission", ylab="Number of Fission Products",main=list(label="M. capitata"))

#First, I fit with Poisson [not correct -- need zero-truncated]
nfiss.UY <- glm(numfiss~size*utrans+size*fyear,data=fissdat,family=poisson)
drop1(nfiss.UY,test="Chisq")
#drop size*utrans p-value=0.3651
nfiss.U0Y <- glm(numfiss~utrans+size*fyear,data=fissdat,family=poisson)
drop1(nfiss.U0Y,test="Chisq")
#drop utrans p-value=0.2287
nfiss.Y<- glm(numfiss~size*fyear,data=fissdat,family=poisson)
drop1(nfiss.Y,test="Chisq")
#drop size*fyear p-value=0.2105
nfiss.Y0<-glm(numfiss~size+fyear,data=fissdat,family=poisson)
drop1(nfiss.Y0,test="Chisq")
#drop fyear p-value=0.1138
nfiss<-glm(numfiss~size,data=fissdat,family=poisson)
drop1(nfiss,test="Chisq")
#Keep size p-value=1.083e-07

#year does not seem to matter
#try with FE only
nfiss.U <- glm(numfiss~size*utrans,data=fissdat,family=poisson)
drop1(nfiss.U,test="Chisq")
#drop size*utrans p-value=0.1842
nfiss.U0 <- glm(numfiss~size+utrans,data=fissdat,family=poisson)
drop1(nfiss.U0,test="Chisq")
#drop utrans p-value=0.3153
nfiss <- glm(numfiss~size,data=fissdat,family=poisson)
drop1(nfiss,test="Chisq")
#keep size p-value=1.083e-07

nfiss.S.Y0 <- glmer(numfiss~size*utrans+(1|fyear),data=fissdat,family=poisson)
nfiss.glmeU_A.Y0 <- glmer(numfiss~size+utrans+(1|fyear),data=fissdat,family=poisson)
nfiss.glmeA.Y0 <- glmer(numfiss~size+(1|fyear),data=fissdat,family=poisson)
anova(nfiss.S.Y0,nfiss.glmeA.Y0,nfiss.glmeU_A.Y0)
#the variance estimates assoc with year are *tiny* and the deviance is 65.7 -not much better
#CONCLUSION: don't bother with REs for year
#also checked FE of year in ztPoisson - no improvement

#ZERO-TRUNCATED POISSON: Next, I realized that I really need a zero-truncated poisson bc there is no possibility of 0 fission products 
# given the structure of the data

nfiss.UY <- vglm(numfiss~utrans*size+fyear,data=fissdat,family=pospoisson)
AIC(nfiss.UY)
#AIC=346.0597
nfiss.U <- vglm(numfiss~utrans*size,data=fissdat,family=pospoisson)
AIC(nfiss.U)
#AIC=350.819
##fyear improves fit by AIC by delta 4.7593
nfiss.U0<- vglm(numfiss~utrans+size,data=fissdat,family=pospoisson)
AIC(nfiss.U0)
#AIC=357.8996
##keep size*utrans interaction
##Zero truncated Poisson gave really crazy estimates of numfiss at larger size values so we tried a few other things
# Zero truncated negative binomial
nfiss.nbz1<-zerotrunc(numfiss~size*utrans,data=fissdat,dist="negbin")
nfiss.nbz2<-zerotrunc(numfiss~size+utrans,data=fissdat,dist="negbin")
nfiss.nbz3<-zerotrunc(numfiss~size,data=fissdat,dist="negbin")
AIC(nfiss.nbz1,nfiss.nbz2,nfiss.nbz3)
# keep size main effect

#Negative Bionomial
nfiss.nb1 <- glm.nb(numfiss~size*utrans,data=fissdat)
nfiss.nb2 <- glm.nb(numfiss~size+utrans,data=fissdat)
nfiss.nb3 <- glm.nb(numfiss~size,data=fissdat)
AIC(nfiss.nb1,nfiss.nb2,nfiss.nb3)
# keep size main effect

#Zero truncated Geometric Distribution
#https://stat.ethz.ch/pipermail/r-help/2011-August/285557.html
nfiss.geo1<-zerotrunc(numfiss~size*utrans,data=fissdat,dist="geometric")
nfiss.geo2<-zerotrunc(numfiss~size+utrans,data=fissdat,dist="geometric")
nfiss.geo3<-zerotrunc(numfiss~size,data=fissdat,dist="geometric")
AIC(nfiss.geo1,nfiss.geo2,nfiss.geo3)
#keep size main effect


##Multinomial
#http://www.ats.ucla.edu/stat/r/dae/mlogit.htm
xtabs(~numfiss,fissdat)
index<-c(1,2,3,4,5,6,7,10)
values<-c("one","two","three",">three",">three",">three",">three",">three")
fissdat$fnumfiss<-values[match(fissdat$numfiss,index)]
xtabs(~fnumfiss,fissdat)

nfiss.multi1<-multinom(fnumfiss~size*utrans,fissdat)
nfiss.multi2<-multinom(fnumfiss~size+utrans,fissdat)
nfiss.multi3<-multinom(fnumfiss~size,fissdat)
AIC(nfiss.multi1,nfiss.multi2,nfiss.multi3)
#keep size main effect

#Number of fission products#
# Several ways: zero truncated poisson, negative binomial, geometric, negative binomial not truncated, multinomial probability of 1-4 products and 1-10 products
dsize<-data.frame(size=x)
clonenum<-vglm(numfiss~size,data=fissdat,family=pospoisson)
clonenb<- glm.nb(numfiss~size,data=fissdat)
clonenbz<-zerotrunc(numfiss~size,data=fissdat,dist="negbin")
clonegeo<-zerotrunc(numfiss~size,data=fissdat,dist="geometric")

plot(fissdat$size,fissdat$numfiss,col=cl[fissdat$utrans],type="p",pch=18,xlab="Size (t)",ylab="# of fission products/'parent colony'",ylim=c(0,10),xlim=c(-2,10))
lines(x,predict(clonenum,newdata=dsize, "response"), col="red")
lines(x,predict(clonenb,newdata=dsize, "response"),col="blue")#for negative binomial distribution
lines(x,predict(clonenbz,newdata=dsize, "response"))
lines(x,predict(clonegeo,newdata=dsize, "response"),col="green")# for geometric distribution
site=1#same for all sites
lines(x,p.c1.x(x,params)+p.c2.x(x,params)*2+p.c3.x(x,params)*3+(1-p.c1.x(x,params)-p.c2.x(x,params)-p.c3.x(x,params))*4,col="magenta")
lines(x,p2.c1.x(x,params)+p2.c2.x(x,params)*2+p2.c3.x(x,params)*3+p2.c4.x(x,params)*4+p2.c6.x(x,params)*6+p2.c7.x(x,params)*7+p2.c10.x(x,params)*10+(1-p2.c1.x(x,params)-p2.c2.x(x,params)-p2.c3.x(x,params)-p2.c4.x(x,params)-p2.c6.x(x,params)-p2.c7.x(x,params)-p2.c10.x(x,params))*5,col="purple")#for multinomial distribution
legend("topleft", legend = c("ztpoisson","negbin","ztnegbin","ztgeometric","multi10","multi4"), col = c("red","blue","black","green","magenta","purple"), lwd = 2, xjust = 1, bg = "white")

#now, do some model evaluation
par(mfrow=c(1,1))
plot(resid(nfiss.U,type="pearson")~fissdat$size)
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
drop1(sfiss.UY,test="Chisq")
#drop size*utrans p-value=0.2169
sfiss.U0Y<- gls(sizeNext~utrans+size*fyear,data=SMfissdat,method="ML")
drop1(sfiss.U0Y,test="Chisq")
#drop size*year p-value=0.10847
sfiss.U0Y0<- gls(sizeNext~size+utrans+fyear,data=SMfissdat,method="ML")
drop1(sfiss.U0Y0,test="Chisq")
#keep size p-value=1.947e-06, utrans p-value=0.03567, fyear p-value=0.04412
##year may matter but we wouldn't keep it as a main effect.  

##Look at FE only
sfiss.U<- gls(sizeNext~size*utrans,data=SMfissdat,method="ML")
drop1(sfiss.U,test="Chisq")
#drop size*utrans p-value=0.2392
sfiss.U0<- gls(sizeNext~size+utrans,data=SMfissdat,method="ML")
drop1(sfiss.U0,test="Chisq")
#drop utrans p-value=0.1417
sfiss<- gls(sizeNext~size,data=SMfissdat,method="ML")
drop1(sfiss,test="Chisq")
#Keep size main effect 2.297e-16


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

##year as RE
sfiss.U.ves<-gls(sizeNext~size*utrans,weights=varExp(form=~size),data=SMfissdat)
sfiss.U.ves.Y <- lme(sizeNext~size*utrans,random=~size|fyear,weights=varExp(form=~size),data=SMfissdat)
sfiss.U.ves.Y0 <- lme(sizeNext~size*utrans, random=~1|fyear,weights=varExp(form=~size),data=SMfissdat)
sfiss.U.ves.Y1 <- lme(sizeNext~size*utrans, random=~size-1|fyear,weights=varExp(form=~size),data=SMfissdat)
anova(sfiss.U.ves,sfiss.U.ves.Y)
anova(sfiss.U.ves,sfiss.U.ves.Y0)
anova(sfiss.U.ves,sfiss.U.ves.Y1)
#RE of year does not improve model p-value=0.4902

##refit best model with reml##
sfiss.ves<-gls(sizeNext~size,data=SMfissdat,weights=varExp(form=~size))
sfiss<-gls(sizeNext~size,data=SMfissdat)
Es<-resid(sfiss.ves,type="normalized")
coplot(Es~SMfissdat$size|SMfissdat$utrans,ylab="normalized residuals")

E<-resid(sfiss,type="normalized")
coplot(E~SMfissdat$size|SMfissdat$utrans,ylab="normalized residuals")