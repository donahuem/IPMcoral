library(lattice)
library(popbio)
#Data collected by EKB 50 colonies measured along each 10m transect in CRAMP transect
#measured length, width and height.  surface area calculated as spherical cap or ellipse
#percent partial mortality estimated for each colony
#size used here is the ln transformed surface area minus the estimated percentage partial mortality times total ln transformed surface area
#See EKB's dissertation for details
EKBSFD <- read.csv("C:/Users/Megan/IPMcoral/R stuff/EKBSFD.csv")
EKBSFD<-subset(EKBSFD,EKBSFD$Reef!="Puamana Deep")
EKBSFD<-droplevels(EKBSFD)
index<-levels(EKBSFD$Reef)
values<-c("Honolua.N","Honolua.S","Olowalu.03m","Olowalu.07m","Puamana.03m")
EKBSFD$utrans<-values[match(EKBSFD$Reef,index)]
EKBSFD$utrans<-as.factor(EKBSFD$utrans)
#subset for Montipora capitat only
EKBMC<-subset(EKBSFD,Species=="Mcap")
#Double check count for subset 
xtabs(~Species,data=EKBMC)
#what morphologies used for M. capitata?
#All ellipse
xtabs(~Morphology,data=EKBMC)
#Adjust for partial mortality size=lnarea-lnarea*%partial mortality/100
EKBMC$size<-EKBMC$Ln.Area-EKBMC$Ln.Area*(EKBMC$X..Mortality/100)
#plot stable.stage distribution estimated by ipm against proportional distribution of EKB's SFD data
par(mfrow=c(2,3))
ipm<-c(1,2,7,8,10)
for (i in 1:5){
  x<-EKBMC$size[EKBMC$utrans==levels(EKBMC$utrans)[i]]
  yfit<-dnorm(y,mean=mean(x),sd=sd(x))
  yfit<-yfit/sum(yfit)
  plot(y,yfit,col="red",type="l",xlim=c(-2,10),ylim=c(0,0.013),xlab= "ln(area)",ylab="proportional distribution",main=levels(EKBMC$utrans)[i])
  site=ipm[i]
  K<-Kernel(y,n,params)
  A<-eigen.analysis(K)
  lines(y,A$stable.stage)
}
legend("topright",legend=c("Eric","IPM"),col=c("red","black"),lty=1)

#plot Eric's SFD vs. SFD of data used to build IPM
par(mfrow=c(2,3))
for (i in 1:5){
  x1<-EKBMC$size[EKBMC$utrans==levels(EKBMC$utrans)[i]]
  x2<-MC$size[MC$utrans==levels(MC$utrans)[ipm[i]]]
  hist(x1,col="red",freq=F,main=levels(EKBMC$utrans)[i],xlim=c(-2,10),ylim=c(0,0.55),xlab="ln(area)")
  hist(x2,col="green",freq=F,add=T)
}
legend("topright", legend = c("Eric","IPM"), col=c("red","green"), pch=15, xjust = 1, bg = "white")

par(mfrow=c(3,4))
for (i in 1:10){
  x<-MC$size[MC$utrans==levels(MC$utrans)[i]]
  yfit<-dnorm(y,mean=mean(x,na.rm=T),sd=sd(x,na.rm=T))
  yfit<-yfit/sum(yfit)
  plot(y,yfit,col="blue",type="l",xlim=c(-2,10),xlab= "ln(area)",ylab="proportional distribution",main=levels(MC$utrans)[i])
  site=i
  K<-Kernel(y,n,params)
  A<-eigen.analysis(K)
  lines(y,A$stable.stage)
}
legend("topright",legend=c("IPMdata","IPMstablestage"),col=c("blue","black"),lty=1)

par(mfrow=c(2,3))
ipm<-c(1,2,7,8,10)
for (i in 1:5){
  x<-EKBMC$size[EKBMC$utrans==levels(EKBMC$utrans)[i]]
  yfit<-dnorm(y,mean=mean(x),sd=sd(x))
  yfit<-yfit/sum(yfit)
  plot(y,yfit,col="red",type="l",xlim=c(-2,10),ylim=c(0,0.013),xlab= "ln(area)",ylab="proportional distribution",main=levels(EKBMC$utrans)[i])
  x<-MC$size[MC$utrans==levels(MC$utrans)[ipm[i]]]
  yfit<-dnorm(y,mean=mean(x,na.rm=T),sd=sd(x,na.rm=T))
  yfit<-yfit/sum(yfit)
  lines(y,yfit,col="blue")
  site=ipm[i]
  K<-Kernel(y,n,params)
  A<-eigen.analysis(K)
  lines(y,A$stable.stage)
}
legend("topright",legend=c("Eric","IPMdata","IPMstablestage"),col=c("red","blue","black"),lty=1)

#MR collected data along 6 CRAMP transects in May 2009
#MR data collected at Kanahena Bay(not used in ipm), Maalaea (not used in ipm),Molokini, Olowalu, Kahekili, Honolua Bay
#1msq quadrats photographed and autonomous tissue traced in imageJ
MRSFD <- read.csv("C:/Users/Megan/IPMcoral/MauiSFD.csv")
MRSFD<-subset(MRSFD,MRSFD$Site!="Kanahena")
MRSFD<-subset(MRSFD,MRSFD$Site!="Maalaea")
MRSFD<-droplevels(MRSFD)
index<-c("HonoluaN","Kahekili","Molokini","Olowalu")
values<-c("Honolua.N","Kahekili.03m","Molokini.13m","Olowalu.03m")
MRSFD$utrans<-values[match(MRSFD$Site,index)]
MRSFD$utrans<-as.factor(MRSFD$utrans)
#ln transform surface area measurments
MRSFD$size<-log(MRSFD$Area)
MRMC<-subset(MRSFD, MRSFD$Species=="M. capitata")
MRMC<-subset(MRMC,MRMC$Area>0)# Not ideal but data only to 1 decimal place.  Would need to go back to image J.  Not sure that it is worth it given quality of photos.  Lose 5 observations all from HonN
par(mfrow=c(2,2))
ipm<-c(1,3,6,7)
for (i in 1:4){
  x<-MRMC$size[MRMC$utrans==levels(MRMC$utrans)[i]]
  yfit<-dnorm(y,mean=mean(x),sd=sd(x))
  yfit<-yfit/sum(yfit)
  plot(y,yfit,col="green",type="l",xlim=c(-2,10),ylim=c(0,0.013),xlab= "ln(area)",ylab="proportional distribution",main=levels(MRMC$utrans)[i])
  site=ipm[i]
  K<-Kernel(y,n,params)
  A<-eigen.analysis(K)
  lines(y,A$stable.stage)
}
legend("topright",legend=c("Megan","IPM"),col=c("green","black"),lty=1)

par(mfrow=c(1,2))
mr<-c(1,4)
ekb<-c(1,3)
ipm<-c(1,7)
for (i in 1:2){
  x<-MRMC$size[MRMC$utrans==levels(MRMC$utrans)[mr[i]]]
  yfit<-dnorm(y,mean=mean(x),sd=sd(x))
  yfit<-yfit/sum(yfit)
  plot(y,yfit,col="green",type="l",xlim=c(-2,10),ylim=c(0,0.013),xlab= "ln(area)",ylab="proportional distribution",main=levels(MRMC$utrans)[mr[i]])
  x<-EKBMC$size[EKBMC$utrans==levels(EKBMC$utrans)[ekb[i]]]
  yfit<-dnorm(y,mean=mean(x),sd=sd(x))
  yfit<-yfit/sum(yfit)
  lines(y,yfit,col="red")
  site=ipm[i]
  K<-Kernel(y,n,params)
  A<-eigen.analysis(K)
  lines(y,A$stable.stage)
}
legend("topright",legend=c("Megan","Eric","IPM"),col=c("green","red","black"),lty=1)


#plot SFD of data used for IPM vs. SFD data collected by MR
par(mfrow=c(2,2))
ipm<-c(1,3,6,7)
for (i in 1:4){
  x1<-MRMC$size[MRMC$utrans==levels(MRMC$utrans)[i]]
  x2<-MC$size[MC$utrans==levels(MC$utrans)[ipm[i]]]
  hist(x1,col="red",freq=F,xlim=c(-2,10),ylim=c(0,0.55),xlab="ln(area)",main=levels(MRMC$utrans)[i])
  hist(x2,col="green",freq=F,add=T)
}
legend("topright", legend = c("Megan","IPM"), col=c("red","green"), pch=15, xjust = 1, bg = "white")


#plot SFD from EKB vs. MR
par(mfrow=c(1,2))
mr<-c(1,4)
ekb<-c(1,3)
for (i in 1:2){
  x1<-MRMC$size[MRMC$utrans==levels(MRMC$utrans)[mr[i]]]
  x2<-EKBMC$size[EKBMC$utrans==levels(EKBMC$utrans)[ekb[i]]]
  hist(x1,col="green",freq=F,xlim=c(-2,10),ylim=c(0,0.55),xlab="ln(area)",main=levels(MRMC$utrans)[mr[i]])
  hist(x2,col="red",freq=F,add=T)
}
legend("topright", legend = c("Eric","Megan"), col=c("red","green"), pch=15, xjust = 1, bg = "white")



KaheSFD<- read.csv("C:/Users/Megan/IPMcoral/KahekiliSFD.csv")
KaheSFD$area<-pi*(0.5*KaheSFD$max.diameter)^2
KaheSFD$lnarea<-log(KaheSFD$area)
KaheMC<-subset(KaheSFD,KaheSFD$spp=="mc")
histogram(~KaheMC$lnarea|KaheMC$Transect)
par(mfrow=c(1,1))
hist(KaheMC$lnarea,freq=F,col="red",main="Kahekili",xlab="ln(area)")
hist(MRSFD$size[MRSFD$Site=="Kahekili"],freq=F,add=T,col="green")
hist(MC$size[MRSFD$utrans=="Kahekili.03m"],freq=F,add=T,col="blue")
legend("topright", legend = c("KaheMC","MRMC","IPM"), col=c("red","green","blue"), pch=15, xjust = 1, bg = "white")


par(mfrow=c(1,1))
x<-MRMC$size[MRMC$utrans==levels(MRMC$utrans)[2]]
yfit<-dnorm(y,mean=mean(x),sd=sd(x))
yfit<-yfit/sum(yfit)
plot(y,yfit,col="green",type="l",xlim=c(-2,10),ylim=c(0,0.013),xlab= "ln(area)",ylab="proportional distribution",main=levels(MRMC$utrans)[2])
yfit<-dnorm(y,mean=mean(KaheMC$lnarea),sd=sd(KaheMC$lnarea))
yfit<-yfit/sum(yfit)
lines(y,yfit,col="red")
site=3
K<-Kernel(y,n,params)
A<-eigen.analysis(K)
lines(y,A$stable.stage)
legend("topright",legend=c("MRMC","KaheMC","IPM"),col=c("green","red","black"),lty=1)


##################################################
#Eric's matrix models from page 151-153 of his disseratation
Olow03<-matrix(c(0.422,0.203,0.000,0.016,0.340,0.489,0.064,0.021,0.200,0.240,0.360,0.160,0.309,0.255,0.164,0.255),4,4)
Olow08<-matrix(c(0.385,0.096,0.019,0.019,0.267,0.200,0.233,0.067,0.300,0.100,0.300,0.250,0.390,0.183,0.061,0.354),4,4)
Puam03<-matrix(c(0.701,0.013,0.000,0.615,0.308,0.000,0.000,0.000,0.000),3,3)
HonoN<-matrix(c(0.538,0.114,0.008,0.438,0.425,0.096,0.387,0.129,0.484),3,3)
HonoS<-matrix(c(0.471,0.039,0.000,0.500,0.450,0.050,0.583,0.167,0.250),3,3)
sizeclasses<-c(0,2.5,10,25)
logsizeclasses<-c(log(0.5),log(2.5),log(10),log(25))
#Olowalu 3m
#use this code to have same scale as IPM
logsizeclasses<-c(-2,log(2.5),log(10),log(25),10)#use this code to have same scale as IPM
midpt<-0.5*(logsizeclasses[1:4]+logsizeclasses[2:(4+1)])
midpt<-midpt[1:4]
image(logsizeclasses,logsizeclasses,t(Olow03),col=topo.colors(300),main="Olowalu 3m",xlab="Size (t)",ylab="Size (t+1)",zlim=c(mi,0.7),xlim=c(-2,10),ylim=c(-2,10))
for (i in 1:4)
  for(j in 1:4)
    text(midpt[i],midpt[j],Olow03[j,i])
mtext("Transition Probabilities", side = 4, outer = TRUE, cex = 1, line = 1,col = "grey20")
vertical.image.legend(col=topo.colors(300),zlim=c(0,0.7))



image(logsizeclasses,logsizeclasses,t(Olow03),col=topo.colors(300),main="Olowalu.03m",xlab="size at t",ylab="size at t+1",zlim=c(mi,0.7))
for (i in 1:4)
  for(j in 1:4)
  text(logsizeclasses[i],logsizeclasses[j],Olow03[j,i])
vertical.image.legend(col=topo.colors(300),zlim=c(0,0.7)
#Olowalu 8m
image(logsizeclasses,logsizeclasses,t(Olow08),col=topo.colors(10),main="Olowalu.08m",xlab="size at t",ylab="size at t+1",zlim=c(0,0.701))
for (i in 1:4)
  for(j in 1:4)
    text(logsizeclasses[i],logsizeclasses[j],Olow08[j,i])
vertical.image.legend(col=topo.colors(10),zlim=c(0,0.701))
#Honolua N
# only 3 size classes
sizeclasses<-c(0,2.5,10)
logsizeclasses<-c(log(0.5),log(2.5),log(10))
image(logsizeclasses,logsizeclasses,t(HonoN),col=topo.colors(10),main="Honolua.N",xlab="size at t",ylab="size at t+1",zlim=c(0,0.701))
for (i in 1:3)
  for(j in 1:3)
    text(logsizeclasses[i],logsizeclasses[j],HonoN[j,i])
vertical.image.legend(col=topo.colors(10),zlim=c(0,0.701))

#Honolua S
image(logsizeclasses,logsizeclasses,t(HonoS),col=topo.colors(10),main="Honolua.S",xlab="size at t",ylab="size at t+1",zlim=c(0,0.701))
for (i in 1:3)
  for(j in 1:3)
    text(logsizeclasses[i],logsizeclasses[j],HonoS[j,i])
vertical.image.legend(col=topo.colors(10),zlim=c(0,0.701))

#Puamana 03m
image(logsizeclasses,logsizeclasses,t(Puam03),col=topo.colors(10),main="Puamana.03m",xlab="size at t",ylab="size at t+1",zlim=c(0,0.701))
for (i in 1:3)
  for(j in 1:3)
    text(logsizeclasses[i],logsizeclasses[j],Puam03[j,i])
vertical.image.legend(col=topo.colors(10),zlim=c(0,0.701))


EKBeig<-rep(0,5)
EKBeig[1]<-eigen.analysis(HonoN)$lambda1
EKBeig[2]<-eigen.analysis(HonoS)$lambda1
EKBeig[3]<-eigen.analysis(Olow03)$lambda1
EKBeig[4]<-eigen.analysis(Olow08)$lambda1
EKBeig[5]<-eigen.analysis(Puam03)$lambda1
EKBeig
## To get lambda for each site
big<-array(0, c(n, n, nx))
eig<-vector("numeric", nx)
for (site in 1:nx){
big[,,site]<-Kernel(y,n, params)
A<-eigen.analysis(big[,,site])
eig[site]<-A$lambda1
}
eig
IPMeig<-c(eig[1:2],eig[7:8],eig[10])
compeig<-matrix(c(EKBeig,IPMeig),5,2)
rownames(compeig)<-c("Honolua N","Honolua S","Olowalu 03m","Olowalu 08m","Puamana 03m")
colnames(compeig)<-c("EKB","IPM")
compeig<-t(compeig)
barplot(compeig,beside=T,col=c("red","green"),xlab="Site",ylab="lambda")
legend("topright",legend=c("EKB","IPM"),col=c("red","green"),pch=15)


par(mfrow=c(5,2))
B<-eigen.analysis(HonoN)$stable.stage
B<-as.matrix(B)
rownames(B)<-c("0-2.5","2.5-10","10-25")
barplot(t(B),xlab="area in square cm",ylab="proportional distribution",main="Honolua.N")
site=1
K<-Kernel(y,n,params)
A<-eigen.analysis(K)
plot(y,A$stable.stage,type="l",xlim=c(-2,10),xlab= "ln(area)",ylab="proportional distribution",main="Honolua.N")

B<-eigen.analysis(HonoS)$stable.stage
B<-as.matrix(B)
rownames(B)<-c("0-2.5","2.5-10","10-25")
barplot(t(B),xlab="area in square cm",ylab="proportional distribution",main="Honolua.S")
site=2
K<-Kernel(y,n,params)
A<-eigen.analysis(K)
plot(y,A$stable.stage,type="l",xlim=c(-2,10),xlab= "ln(area)",ylab="proportional distribution",main="Honolua.S")

B<-eigen.analysis(Olow03)$stable.stage
B<-as.matrix(B)
rownames(B)<-c("0-2.5","2.5-10","10-25",">25")
barplot(t(B),xlab="area in square cm",ylab="proportional distribution",main="Olowalu.03m")
site=7
K<-Kernel(y,n,params)
A<-eigen.analysis(K)
plot(y,A$stable.stage,type="l",xlim=c(-2,10),xlab= "ln(area)",ylab="proportional distribution",main="Olowalu.03m")

B<-eigen.analysis(Olow08)$stable.stage
B<-as.matrix(B)
rownames(B)<-c("0-2.5","2.5-10","10-25",">25")
barplot(t(B),xlab="area in square cm",ylab="proportional distribution",main="Olowalu.08m")
site=8
K<-Kernel(y,n,params)
A<-eigen.analysis(K)
plot(y,A$stable.stage,type="l",xlim=c(-2,10),xlab= "ln(area)",ylab="proportional distribution",main="Olowalu.08m")

B<-eigen.analysis(Puam03)$stable.stage
B<-as.matrix(B)
rownames(B)<-c("0-2.5","2.5-10","10-25")
barplot(t(B),xlab="area in square cm",ylab="proportional distribution",main="Puamana.03m")
site=10
K<-Kernel(y,n,params)
A<-eigen.analysis(K)
plot(y,A$stable.stage,type="l",xlim=c(-2,10),xlab= "ln(area)",ylab="proportional distribution",main="Puamana.03m")


par(mfrow=c(2,2))
sizeclasses<-c(0,2.5,10,25)
logsizeclasses<-c(log(0.5),log(2.5),log(10),log(25))
#Olowalu 3m
A<-eigen.analysis(Olow03)
m<-max(A$elasticities)
image(logsizeclasses,logsizeclasses,t(A$elasticities),col=topo.colors(10),main="Olowalu.03m",xlab="size at t",ylab="size at t+1",zlim=c(0,m))
for (i in 1:4)
  for(j in 1:4)
    text(logsizeclasses[i],logsizeclasses[j],A$elasticities[j,i])
vertical.image.legend(col=topo.colors(10),zlim=c(0,m))
site=7
K<-Kernel(y,n,params)
A<-eigen.analysis(K)
m<-max(A$elasticities)
image(y,y,t(A$elasticities),col=topo.colors(300),main=levels(MC$utrans)[site],xlab="size at t",ylab="size at t+1",zlim=c(0,m))
vertical.image.legend(col=topo.colors(300),zlim=c(0,m))

#Olowalu 8m
A<-eigen.analysis(Olow08)
m<-max(A$elasticities)
image(logsizeclasses,logsizeclasses,t(A$elasticities),col=topo.colors(10),main="Olowalu.08m",xlab="size at t",ylab="size at t+1",zlim=c(0,m))
for (i in 1:4)
  for(j in 1:4)
    text(logsizeclasses[i],logsizeclasses[j],A$elasticities[j,i])
vertical.image.legend(col=topo.colors(10),zlim=c(0,m))

site=8
K<-Kernel(y,n,params)
A<-eigen.analysis(K)
m<-max(A$elasticities)
image(y,y,t(A$elasticities),col=topo.colors(300),main=levels(MC$utrans)[site],xlab="size at t",ylab="size at t+1",zlim=c(0,m))
vertical.image.legend(col=topo.colors(300),zlim=c(0,m))

#Honolua N
# only 3 size classes
sizeclasses<-c(0,2.5,10)
logsizeclasses<-c(log(0.5),log(2.5),log(10))
A<-eigen.analysis(HonoN)
m<-max(A$elasticities)
image(logsizeclasses,logsizeclasses,t(A$elasticities),col=topo.colors(10),main="Honolua.N",xlab="size at t",ylab="size at t+1",zlim=c(0,m))
for (i in 1:3)
  for(j in 1:3)
    text(logsizeclasses[i],logsizeclasses[j],A$elasticities[j,i])
vertical.image.legend(col=topo.colors(10),zlim=c(0,m))


site=1
K<-Kernel(y,n,params)
A<-eigen.analysis(K)
m<-max(A$elasticities)
image(y,y,t(A$elasticities),col=topo.colors(300),main=levels(MC$utrans)[site],xlab="size at t",ylab="size at t+1",zlim=c(0,m))
vertical.image.legend(col=topo.colors(300),zlim=c(0,m))

#Honolua S
A<-eigen.analysis(HonoS)
m<-max(A$elasticities)
image(logsizeclasses,logsizeclasses,t(A$elasticities),col=topo.colors(10),main="Honolua.S",xlab="size at t",ylab="size at t+1",zlim=c(0,m))
for (i in 1:3)
  for(j in 1:3)
    text(logsizeclasses[i],logsizeclasses[j],A$elasticities[j,i])
vertical.image.legend(col=topo.colors(10),zlim=c(0,m))

site=2
K<-Kernel(y,n,params)
A<-eigen.analysis(K)
m<-max(A$elasticities)
image(y,y,t(A$elasticities),col=topo.colors(300),main=levels(MC$utrans)[site],xlab="size at t",ylab="size at t+1",zlim=c(0,m))
vertical.image.legend(col=topo.colors(300),zlim=c(0,m))

#Puamana 03m
A<-eigen.analysis(Puam03)
m<-max(A$elasticities)
image(logsizeclasses,logsizeclasses,t(A$elasticities),col=topo.colors(10),main="Puamana.03m",xlab="size at t",ylab="size at t+1",zlim=c(0,m))
for (i in 1:3)
  for(j in 1:3)
    text(logsizeclasses[i],logsizeclasses[j],A$elasticities[j,i])
vertical.image.legend(col=topo.colors(10),zlim=c(0,m))

site=10
K<-Kernel(y,n,params)
A<-eigen.analysis(K)
ma<-max(A$elasticities)
mi<-min(A$elasticities)
image(y,y,t(A$elasticities),col=topo.colors(300),main=levels(MC$utrans)[site],xlab="size at t",ylab="size at t+1",zlim=c(mi,ma))
vertical.image.legend(col=topo.colors(300),zlim=c(mi,ma))
