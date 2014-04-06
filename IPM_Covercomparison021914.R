cl<-c("cyan","mediumturquoise","dodgerblue","green3","purple","olivedrab2","orchid1","blue","darkblue","darkorchid4")

## To get lambda for each site
big<-array(0, c(n, n, nx))
eig<-vector("numeric", nx)
for (site in 1:nx){
  big[,,site]<-Kernel(y,n, params)
  A<-eigen.analysis(big[,,site])
  eig[site]<-A$lambda1
}
eig
##Lambda is multiplicative growth rate of number of colonies.  This may not be directly related to cover.
##What if we multiply lambda by the stable size distribution and the size in area of all of the size classes?
ss<-rep(0,10)
for (site in 1:10){
  K<-Kernel(y,n,params)
  A<-eigen.analysis(K)
  ss[site]<-sum(A$stable.stage*y)
}
ssl<-rep(0,10)
for(site in 1:10){
  K<-Kernel(y,n,params)
  A<-eigen.analysis(K)
  ssl[site]<-sum(A$stable.stage*y)*A$lambda1
}
ssp<-(ssl-ss)/ss
ssd<-ssl-ss

ss1<-rep(0,10)
for (site in 1:10){
  K<-Kernel(y,n,params)
  A<-eigen.analysis(K)
  ss[site]<-sum(A$stable.stage*exp(y))
}
ssl1<-rep(0,10)
for(site in 1:10){
  K<-Kernel(y,n,params)
  A<-eigen.analysis(K)
  ssl[site]<-sum(A$stable.stage*exp(y))*A$lambda1
}
ssp1<-(ssl1-ss1)/ss1
ssd1<-ssl1-ss1
AllC<- read.csv("C:/Users/Megan/IPMcoral/MauiCover.csv")
AllC<-na.omit(AllC)
AllC$trans<-paste(AllC$Site,AllC$Depth..m.)
AllC$trans<-as.factor(AllC$trans)
index<-levels(AllC$trans)
values<-c(levels(MC$utrans)[1:4],levels(MC$utrans)[6],levels(MC$utrans)[5],levels(MC$utrans)[7:10])
AllC$utrans<-values[match(AllC$trans,index)]
AllC$utrans<-as.factor(AllC$utrans)
AllC$fyear<-as.factor(AllC$year)
######Maui coral cover M. capitata only
MCC<- subset(AllC,AllC$Taxon=="Montipora capitata")
MCC<-na.omit(MCC)
MCC$pcover<-MCC$cover/100
MCC$asinsqrt<-asin(sqrt(MCC$pcover))
par(mfrow=c(3,4))
for(i in 1:10){
  x<-subset(MCC,MCC$utrans==levels(MCC$utrans)[i])
  plot(x$year,x$cover,pch=18,main=levels(MCC$utrans)[i],xlab="year",ylab="M. capitata cover")
}

###M. capitata cover slope over 1999-2012 vs. lambda
slopes<-matrix(0,10,4)
for (i in 1:10){
  x<-subset(MCC,MCC$utrans==levels(MCC$utrans)[i])
  overall<-lm(asinsqrt~year,data=x)
  slopes[i,1]<-coefficients(overall)[2]
  slopes[i,2]<-anova(overall)[1,5]
}
slopes<-as.data.frame(slopes)
colnames(slopes)<-c("slopes","p-value","lambda","utrans")
slopes$eig<-eig
slopes$utrans<-as.factor(levels(MC$utrans))

par(mfrow=c(1,1))
plot(slopes$slopes,slopes$eig,col=cl[slopes$utrans],pch=18,xlab="M. capitata %coverage slope 1999-2012",ylab="IPM lambda")
legend("bottomright",legend=sprintf(levels(slopes$utrans)),col=cl,pch=18)
cor.test(slopes$slopes,slopes$eig)
#r=0.4346596, p-value=0.2094

slopes$ssp<-ssp
plot(slopes$slopes,slopes$ssp,col=cl[slopes$utrans],pch=18,xlab="M. capitata % coverage slope transects 1999-2012",ylab="[(ssd*y*lambda)-(ssd*y)]/ssd")
cor.test(slopes$slopes,slopes$ssp)
legend("bottomright",legend=sprintf(levels(slopes$utrans)),col=cl,pch=18)
slopes$ssd<-ssl-ss
plot(slopes$slopes,slopes$ssd,col=cl[slopes$utrans],pch=18,xlab="M. capitata % coverage slope transects 1999-2012",ylab="(ssd*y*lambda)-(ssd*y)")
cor.test(slopes$slopes,slopes$ssd)
legend("topleft",legend=sprintf(levels(slopes$utrans)),col=cl,pch=18)

###Does transect cover data match quadrat cover data?
quadcov<-data.frame(MC$uyear,MC$spcover,MC$tyear)
quadcov<-unique(quadcov)
meancov<-tapply(quadcov$MC.spcover,c(quadcov$MC.tyear),mean)
meanquadcov<-data.frame(MC$utrans,MC$year)
meanquadcov<-unique(meanquadcov)
meanquadcov$mean<-meancov
meanquadcov$cover<-meanquadcov$mean/3290*100
colnames(meanquadcov)<-c("utrans","year","mean","cover")

par(mfrow=c(3,4))
for(i in 1:10){
  x<-subset(MCC,MCC$utrans==levels(MCC$utrans)[i])
  x1<-subset(meanquadcov,meanquadcov$utrans==levels(meanquadcov$utrans)[i])
  maxs<-c(max(x$cover),max(x1$cover))
  ma<-max(maxs)
  mins<-c(min(x$cover),min(x1$cover))
  mi<-min(mins)
  plot(x$year,x$cover,pch=18,main=levels(MCC$utrans)[i],xlab="year",ylab="M. capitata cover",ylim=c(mi,ma))
  points(x1$year,x1$cover,pch=18,col="red")
}


###Does quadrat M. capitata cover data match quadrat lambda data?
meanquadcov$pcover<-meanquadcov$cover/100
meanquadcov$asinsqrt<-asin(sqrt(meanquadcov$pcover))
par(mfrow=c(3,4))

###M. capitata cover from quadrate data slope over 1999-2012 vs. lambda
slopes<-rep(0,10)
for (i in 1:10){
  x<-subset(meanquadcov,meanquadcov$utrans==levels(meanquadcov$utrans)[i])
  overall<-lm(asinsqrt~year,data=x)
  slopes[i]<-coefficients(overall)[2]
}

slopes<-as.data.frame(slopes)
slopes$eig<-eig
slopes$utrans<-levels(MCC$utrans)
slopes$utrans<-as.factor(slopes$utrans)
par(mfrow=c(1,1))
plot(slopes$slopes,slopes$eig,col=cl[slopes$utrans],pch=18,xlab="M. capitata %coverage slope from quadrats 1999-2005",ylab="IPM lambda")
legend("bottomright",legend=sprintf(levels(slopes$utrans)),col=cl,pch=18)
slopes$ssp<-ssp
plot(slopes$slopes,slopes$ssp,col=cl[slopes$utrans],pch=18,xlab="M. capitata %coverage slope from quadrats",ylab="(ssd*y*lambda)-(ssd*y)/ssd")
legend("bottomright",legend=sprintf(levels(slopes$utrans)),col=cl,pch=18)

slopes$ssd<-ssl-ss
plot(slopes$slopes,slopes$ssd,col=cl[slopes$utrans],pch=18,xlab="M. capitatav% coverage slope from quadrats",ylab="(ssd*y*lambda)-(ssd*y)")
legend("bottomright",legend=sprintf(levels(slopes$utrans)),col=cl,pch=18)
cor.test(slopes$slopes,slopes$ssd)
#r=0.93001173,p-value=9.638e-05

##Truncated transect M. capitat coverage to match regression above 
MCCtrunc<-subset(MCC,MCC$fyear=="1999"|MCC$fyear=="2000"|MCC$fyear=="2001"|MCC$fyear=="2002"|MCC$fyear=="2003"|MCC$fyear=="2004"|MCC$fyear=="2005")
MCCtrunc<-droplevels(MCCtrunc)
par(mfrow=c(3,4))
for(i in 1:10){
  x<-subset(MCCtrunc,MCCtrunc$utrans==levels(MCCtrunc$utrans)[i])
  plot(x$year,x$cover,pch=18,main=levels(MCCtrunc$utrans)[i],xlab="year",ylab="M. capitata cover")
}

###M. capitata cover slope over 1999-2012 vs. lambda
slopes<-rep(0,10)
for (i in 1:10){
  x<-subset(MCCtrunc,MCCtrunc$utrans==levels(MCCtrunc$utrans)[i])
  overall<-lm(asinsqrt~year,data=x)
  slopes[i]<-coefficients(overall)[2]
}

slopes<-as.data.frame(slopes)
slopes$eig<-eig
slopes$utrans<-levels(MCCtrunc$utrans)
slopes$utrans<-as.factor(slopes$utrans)
par(mfrow=c(1,1))
plot(slopes$slopes,slopes$eig,col=cl[slopes$utrans],pch=18,xlab="M. capitata %coverage slope transects truncated 1999-2005",ylab="IPM lambda")
legend("bottomright",legend=sprintf(levels(slopes$utrans)),col=cl,pch=18)


slopes$ssp<-ssp
plot(slopes$slopes,slopes$ssp,col=cl[slopes$utrans],pch=18,xlab="M. capitata % coverage slope transects truncated 1999-2005",ylab="[(ssd*y*lambda)-(ssd*y)]/lambda")
legend("bottomright",legend=sprintf(levels(slopes$utrans)),col=cl,pch=18)
slopes$ssd<-ssl-ss
plot(slopes$slopes,slopes$ssd,col=cl[slopes$utrans],pch=18,xlab="M. capitata % coverage slope transects truncated 1999-2005",ylab="(ssd*y*lambda)-(ssd*y)")

cor.test(slopes$slopes,slopes$ssd)
#r=0.8033003, p=value=0.005127


###Does transect Total coral cover data match quadrat cover data?
quadcov<-data.frame(MC$uyear,MC$allcover,MC$tyear)
quadcov<-unique(quadcov)
meancov<-tapply(quadcov$MC.allcover,c(quadcov$MC.tyear),mean)
meanquadcov<-data.frame(MC$utrans,MC$year)
meanquadcov<-unique(meanquadcov)
meanquadcov$mean<-meancov
meanquadcov$cover<-meanquadcov$mean/3290*100
colnames(meanquadcov)<-c("utrans","year","mean","cover")

par(mfrow=c(3,4))
for(i in 1:10){
  x<-subset(TCC,TCC$utrans==levels(TCC$utrans)[i])
  x1<-subset(meanquadcov,meanquadcov$utrans==levels(meanquadcov$utrans)[i])
  maxs<-c(max(x$cover),max(x1$cover))
  ma<-max(maxs)
  mins<-c(min(x$cover),min(x1$cover))
  mi<-min(mins)
  plot(x$year,x$cover,pch=18,main=levels(MCC$utrans)[i],xlab="year",ylab="total coral cover",ylim=c(mi,ma))
  points(x1$year,x1$cover,pch=18,col="red")
}