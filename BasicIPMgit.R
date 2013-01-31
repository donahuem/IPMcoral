######Basic IPMpack Olowalu, MC, all years, cheat Fecundity, colonal Matrix######
rm(list=ls(all=TRUE))
setwd("C://Users/Megan/IPMcoral")
library(IPMpack)
data<- read.delim("sortIPMpack.txt")

######Clean up data: remove M, add NA, make recruit fecs,add stage and stageNext######
data<-subset(data,fate_t!="M")
data$fate_prior[data$fate_prior==""]<-NA
data$offspringNext[data$offspringNext==""]<-NA
data$fuseLGname[data$fuseLGname==""]<-NA
data$size<-log(data$asize)
data$sizeNext<-log(data$asizeNext)
data$surv<-data$survIPM
index<-c("growth","dead","fissSM","fissLG","fuseSM","fuseLG","recruit")
values<-c(1,NA,NA,1,NA,1,NA)
data$rec<-values[match(data$fateIPM,index)]
values2<-c("continuous","continuous","continuous","continuous","continuous","continuous",NA)
data$stage<-as.factor(values2[match(data$fateIPM,index)])
values3<-c("continuous","dead","continuous","continuous","dead","continuous","continuous")
data$stageNext<-as.factor(values3[match(data$fateIPM,index)])
data$fec4<-ifelse(data$fec1==1,data$fec2,NA)
#######Using Olowalu MC because best n###########
OMC<-subset(data,sitename=="Olow"&species=="MC")
######growth,fuseLg, fissLg only for growth makeGrowthObj only removes values where size or sizeNext are NA so keeps fuseSM and fissSM. We don't want those for growth######
OMCg<-subset(OMC,fateIPM=="growth"|fateIPM=="fuseLG"|fateIPM=="fissLG")
######Construct growth+survival=Pmatrix using IPMpack Vignette 1.5#######
par(mfrow=c(1,2),bty="l",pty="m")##for growth and survival comp plots
growthcomp1 <- growthModelComp(OMC,expVars=c(sizeNext~1, sizeNext~size, sizeNext~size + size2),regressionType="constantVar",testType="AIC",makePlot=TRUE, ylim=c(-1,10))
gr1 <- makeGrowthObj(dataf = OMC, Formula = sizeNext~size)
survcomp1 <- survModelComp(OMC,expVars=c(surv~1, surv~size, surv~size + size2),testType="AIC",makePlot=TRUE,ylim=c(0.4,1))
sv1 <- makeSurvObj(OMC, Formula = surv~size)
par(mfrow = c(1, 3), bty = "l", pty = "m")##for growth, survival function plots and image of Pmatrix kernel
p1 <- picGrow(OMCg, gr1)
p2 <- picSurv(OMC, sv1, ncuts = 30)
minsize<-min(OMC$size,na.rm=T)
maxsize<-max(OMC$size,na.rm=T)+2######Added 2 to maxsize to improve diagnostics. Seems reasonable makes max size a little less than 1.5m####
Pmatrix <- createIPMPmatrix(nBigMatrix = 100,minSize = minsize, maxSize = maxsize,growObj = gr1, survObj = sv1)
#diagnosticsPmatrix(Pmatrix,gr1,sv1)
filled.contour(Pmatrix@meshpoints,Pmatrix@meshpoints,t(Pmatrix),color=heat.colors,nlevels=20,cex.lab=1.5,main="Pmatrix:survival and growth", xlab= "Size at t",ylab="Size at t+1")
LE1 <- meanLifeExpect(Pmatrix)
#plot(LE1,Pmatrix@meshpoints)
pTime1 <- passageTime(mean(OMC$size, na.rm = TRUE), Pmatrix)

######Fecundity Matrix######

#######How will we deal with recruitment?######
#######MeganD and I talked about adding some number of recruits to the size distribution vector with each iteration######
#######That would work but I wanted to try options to keep probability in matrix for eigenvalue/vector calculations######
#######Bruno aspergillosis paper used constant R*probability density.  We can force makeFecObj to do this.  Was really overthinking this one but Cory's code helped me understand IPMpack object######
######Taken from Jess' Monocarp exercise page5######
dummy.fit<-glm(rec~1,data=OMC,family=poisson)
######chose this intercept value based on average number of recruits per adult from Olowalu MC data######
######could use any value######
dummy.fit$coefficients[1]=0.0905897
fo<-new("fecObj")
fo@fitFec[[1]]<-dummy.fit
fo@Transform<-"none"
#fo@fecNames <- c("rec")
fo@fecConstants<-data.frame(NA)
fo@fecByDiscrete<-data.frame(NA)
dummy.fit1 <- lm(sizeNext~1, data=OMC)
dummy.fit1$coefficients[1] <- mean(OMC$sizeNext[OMC$fateIPM=="recruit"])
fo@offspringRel <- dummy.fit1
fo@sdOffspringSize <- sd(OMC$sizeNext[OMC$fateIPM=="recruit"])
fo@offspringSplitter <- data.frame(continuous=1)
fo@vitalRatesPerOffspringType <- data.frame(continuous=c(1,1))
rownames(fo@vitalRatesPerOffspringType) <- c("rec","NA")

Fmatrix<-createIPMFmatrix(fecObj=fo,minSize=minsize,maxSize=maxsize,nBigMatrix = 100)
par(mfrow = c(1, 1), bty = "l", pty = "m")
filled.contour(Fmatrix@meshpoints,Fmatrix@meshpoints,t(Fmatrix),color=heat.colors,nlevels=20,cex.lab=1.5,main="Fmatrix:Sexual Recruitment", xlab= "Size at t",ylab="Size at t+1")


####Clonal Matrix from Eelke's exercise#######

# whether corals fission or not
#No clonalModelComp() yet so used Eelke's code
co1<-makeClonalObj(OMC, Formula = fec1~1, Family = c("binomial"), Transform=c("none"))
co2<-makeClonalObj(OMC, Formula = fec1~size, Family = c("binomial"), Transform=c("none"))
co3<-makeClonalObj(OMC, Formula = fec1~size+size2, Family = c("binomial"), Transform=c("none"))

x<-seq(from=-2,to=10,length=1001)
x0<-data.frame(size=x,size2=x*x)
par(mfrow = c(1, 3), bty = "l", pty = "m")

# plot these models for comparison
os <- order(OMC$size)  
osClon <- (OMC$fec1)[os] 
osSize<-(OMC$size)[os] 
binnedSize <- tapply(osSize, as.numeric(cut(osSize, breaks=20)), mean, na.rm = TRUE); # bin Size data
binnedClon <- tapply(osClon, as.numeric(cut(osSize, breaks=20)), mean, na.rm = TRUE) #bin Survival probabilities
plot(binnedSize, binnedClon, pch = 19, cex=2, xlab = "Size at t", ylab = "Producing clones at t+1", main = "Whether or not fission")
y0<-predict(co1@fitFec[[1]],newdata=x0);y0<-exp(y0)/(exp(y0)+1);lines(x,y0,col="red")
y0<-predict(co2@fitFec[[1]],newdata=x0);y0<-exp(y0)/(exp(y0)+1);lines(x,y0,col="green")
y0<-predict(co3@fitFec[[1]],newdata=x0);y0<-exp(y0)/(exp(y0)+1);lines(x,y0,col="blue")
legend("topleft", legend = sprintf("%s: %s = %.1f",c("1","size","size+size2"), c("AIC"), c(AIC(co1@fitFec[[1]]),AIC(co2@fitFec[[1]]),AIC(co3@fitFec[[1]]))), col = c(2:4), lty = 1, xjust = 1, bg = "white")

# number of fission products
plot(OMC$size, OMC$fec2, pch = 19, xlab = "Size at t", ylab = "# fission products at t+1", main = "# fissSM when fissing", cex = 2)
co1 <- makeClonalObj(OMC, Formula = c(fec1~size, fec2~1), Family = c("binomial","poisson"), Transform=c("none","none"))
co2 <- makeClonalObj(OMC, Formula = c(fec1~size, fec2~size), Family = c("binomial","poisson"), Transform=c("none","none"))
y0<-predict(co1@fitFec[[2]],newdata=x0);y0<-exp(y0)+1;lines(x,y0,col="red")
y0<-predict(co2@fitFec[[2]],newdata=x0);y0<-exp(y0)+1;lines(x,y0,col="green")
legend("topleft", legend = sprintf("%s: %s = %.1f",c("1","size"), c("AIC"), c(AIC(co1@fitFec[[2]]),AIC(co2@fitFec[[2]]))), col = c(2:4), lty = 1, xjust = 1, bg = "white")

# fissSM size in relation to fissLG size
plot(OMC$size[OMC$offspringNext=="clonal"], OMC$sizeNext[OMC$offspringNext=="clonal"], pch = 19, xlab = "FissLG Size at t", ylab = "# FissSM size at t+1", main = "FissLG - FissSM size relationship", cex = 2)
co1<-makeClonalObj(OMC, offspringSizeExplanatoryVariables = "1", Formula = c(fec1~size, fec4~size), Family = c("binomial","poisson"), Transform=c("none","-1"))
co2<-makeClonalObj(OMC, offspringSizeExplanatoryVariables = "size", Formula = c(fec1~size, fec4~size), Family = c("binomial","poisson"), Transform=c("none","-1"))
co3<-makeClonalObj(OMC, offspringSizeExplanatoryVariables = "size+size2", Formula = c(fec1~size, fec4~size), Family = c("binomial","poisson"), Transform=c("none","-1"))
#plot these models for comparison
y0<-predict(co1@offspringRel,newdata=x0);lines(x,y0,col="red")
y0<-predict(co2@offspringRel,newdata=x0);lines(x,y0,col="green")
y0<-predict(co3@offspringRel,newdata=x0);lines(x,y0,col="blue")
legend("topleft", legend = sprintf("%s: %s = %.1f",c("1","size","size+size2"), c("AIC"), c(AIC(co1@offspringRel),AIC(co2@offspringRel),AIC(co3@offspringRel))), col = c(2:4), lty = 1, xjust = 1, bg = "white")
co<-makeClonalObj(OMC, offspringSizeExplanatoryVariables = "size", Formula = c(fec1~size, fec4~size), Family = c("binomial","poisson"), Transform=c("none","none"))
Cmatrix<-createIPMCmatrix(nBigMatrix=100,clonalObj=co,minSize=minsize,maxSize=maxsize)
par(mfrow = c(1, 1), bty = "l", pty = "m")
filled.contour(Cmatrix@meshpoints,Cmatrix@meshpoints,t(Cmatrix),color=heat.colors,nlevels=20,cex.lab=1.5,main="Cmatrix:Fission Recruitment", xlab= "Size at t",ylab="Size at t+1")
######Full IPM#######
IPM<-Pmatrix+Fmatrix+Cmatrix
Re(eigen(IPM)$value[1])
sensitivity1 <- sens(IPM)
elasticity1 <- elas(IPM)
filled.contour(Cmatrix@meshpoints,Cmatrix@meshpoints,t(Cmatrix),color=heat.colors,nlevels=20,cex.lab=1.5,main="Cmatrix:Fission Recruitment", xlab= "Size at t",ylab="Size at t+1")
persp(Pmatrix@meshpoints,Pmatrix@meshpoints,t(IPM),phi=30,theta=300,col="lightgrey",border=NA,shade=0.75,main="Pmatrix:survival and growth", xlab= "Size at t",ylab="Size at t+1")
sens<-sensParams(growObj=gr1, survObj=sv1, fecObj=fo, nBigMatrix=100, minSize=minsize, maxSize=maxsize)
barplot(sens$elam, main = expression("Parameter elasticity of population growth rate "* lambda), las = 2, cex.names = 0.5)

Pmatrix <- createIPMPmatrix(nBigMatrix = 100,minSize = minsize, maxSize = maxsize,growObj = gr1, survObj = sv1,correction="constant")
Fmatrix<-createIPMFmatrix(fecObj=fo,minSize=minsize,maxSize=maxsize,nBigMatrix = 100,correction="constant")
Cmatrix<-createIPMCmatrix(nBigMatrix=100,clonalObj=co,minSize=minsize,maxSize=maxsize,correction="constant")
IPMconstant<-Pmatrix+Fmatrix+Cmatrix
Re(eigen(IPMconstant)$value[1])