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
####Maui coral cover by speices
AllC<- read.csv("C:/Users/Megan/IPMcoral/MauiCover.csv")
AllC<-na.omit(AllC)
AllC$trans<-paste(AllC$Site,AllC$Depth..m.)
AllC$trans<-as.factor(AllC$trans)
index<-levels(AllC$trans)
values<-c(levels(MC$utrans)[1:4],levels(MC$utrans)[6],levels(MC$utrans)[5],levels(MC$utrans)[7:10])
AllC$utrans<-values[match(AllC$trans,index)]
AllC$utrans<-as.factor(AllC$utrans)
AllC$fyear<-as.factor(AllC$year)
par(mfrow=c(3,4))
for(i in 1:10){
  x<-subset(AllC,AllC$utrans==levels(AllC$utrans)[i])
  m<-max(AllC$cover[AllC$utrans==levels(AllC$utrans)[i]],na.rm=T)
  plot(x$year[x$Taxon==levels(x$Taxon)[1]],x$cover[x$Taxon==levels(x$Taxon)[1]],type="l",ylim=c(0,m),col="red",main=levels(x$utrans)[i],xlab="year",ylab="% cover by spp")
  lines(x$year[x$Taxon==levels(x$Taxon)[2]],x$cover[x$Taxon==levels(x$Taxon)[2]],col="blue")
  lines(x$year[x$Taxon==levels(x$Taxon)[3]],x$cover[x$Taxon==levels(x$Taxon)[3]],col="green")
  lines(x$year[x$Taxon==levels(x$Taxon)[4]],x$cover[x$Taxon==levels(x$Taxon)[4]],col="purple")
  lines(x$year[x$Taxon==levels(x$Taxon)[5]],x$cover[x$Taxon==levels(x$Taxon)[5]],col="magenta")
  lines(x$year[x$Taxon==levels(x$Taxon)[6]],x$cover[x$Taxon==levels(x$Taxon)[6]],col="black")
}
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
legend("center", legend = sprintf(levels(AllC$Taxon)), col = c("red","blue","green","purple","magenta","black"), lwd = 1, xjust = 1, bg = "white")

ACtotal<-matrix(0,10,14)
for(i in 1:10){
  x<-subset(AllC,AllC$utrans==levels(AllC$utrans)[i])
  ACtotal[i,]<-tapply(x$cover,x$fyear,sum)
}
ACtotal<-as.data.frame(ACtotal)
ACtotal<-stack(ACtotal)
years<-rep(c(1999:2012),each=10)
ACtotal$year<-years
utrans1<-levels(AllC$utrans)
ACtotal$utrans<-rep(utrans1,14)
ACtotal$utrans<-as.factor(ACtotal$utrans)

par(mfrow=c(3,4))
for(i in 1:10){
  x<-subset(ACtotal,ACtotal$utrans==levels(ACtotal$utrans)[i])
  plot(x$year,x$values,type="l",xlab="year",ylab="% cover of coral (top 6 species)",main=levels(x$utrans)[i])
}

###Maui Total coral cover 
TCC<- read.csv("C:/Users/Megan/IPMcoral/MauiTCC.csv")
TCC<-na.omit(TCC)
TCC$pcover<-TCC$cover/100
TCC$asinsqrt<-asin(sqrt(TCC$pcover))
par(mfrow=c(3,4))
for(i in 1:10){
  x<-subset(TCC,TCC$utrans==levels(TCC$utrans)[i])
  plot(x$year,x$cover,type="l",xlab="year",ylab="% cover of coral (top 6 species)",main=levels(x$utrans)[i])
}


##Maui total coral cover slope over 1999-2012 vs lambda
slopes<-rep(0,10)
for (i in 1:10){
  x<-subset(TCC,TCC$utrans==levels(TCC$utrans)[i])
  overall<-lm(asinsqrt~year,data=x)
  slopes[i]<-coefficients(overall)[2]
}
slopes<-as.data.frame(slopes)
slopes$eig<-eig
slopes$utrans<-levels(TCC$utrans)
slopes$utrans<-as.factor(slopes$utrans)
par(mfrow=c(1,1))
plot(slopes$slopes,slopes$eig,col=cl[slopes$utrans],pch=18,xlab="Cramp %coverage slope",ylab="IPM lambda")
legend("bottomleft",legend=sprintf(levels(slopes$utrans)),col=cl,pch=18)
cor.test(slopes$slopes,slopes$eig)
#r=0.1826784, p-value=0.6135

#remove "outliers"
#identify(slopes$slopes,slopes$eig)
slopes1<-rbind(slopes[1,],slopes[3:8,],slopes[10,])
plot(slopes1$slopes,slopes1$eig,col=cl[slopes1$utrans],pch=18,xlab="Cramp %coverage slope",ylab="IPM lambda")
legend("topright",legend=sprintf(levels(slopes1$utrans)),col=cl,pch=18)
cor.test(slopes1$slopes,slopes1$eig)
#r=-0.2782645, p-value=0.5046


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
slopes<-rep(0,10)
for (i in 1:10){
  x<-subset(MCC,MCC$utrans==levels(MCC$utrans)[i])
  overall<-lm(asinsqrt~year,data=x)
  slopes[i]<-coefficients(overall)[2]
}
slopes<-as.data.frame(slopes)
slopes$eig<-eig
slopes$utrans<-levels(MCC$utrans)
slopes$utrans<-as.factor(slopes$utrans)
par(mfrow=c(1,1))
plot(slopes$slopes,slopes$eig,col=cl[slopes$utrans],pch=18,xlab="M. capitata %coverage slope 1999-2012",ylab="IPM lambda")
legend("bottomright",legend=sprintf(levels(slopes$utrans)),col=cl,pch=18)
cor.test(slopes$slopes,slopes$eig)
#r=0.5024241, p-value=0.1389


#remove "outliers"
slopes1<-rbind(slopes[1,],slopes[3:8,],slopes[10,])
plot(slopes1$slopes,slopes1$eig,col=cl[slopes1$utrans],pch=18,xlab="Cramp % MC coverage slope",ylab="IPM lambda")
legend("topright",legend=sprintf(levels(slopes1$utrans)),col=cl,pch=18)
cor.test(slopes1$slopes,slopes1$eig)
#r=0.4593588, p-value=0.2522

MCC$fyear<-as.factor(MCC$year)
year.slope<-matrix(0,nrow=10,ncol=13)
for (i in (1999:2011)){
  x<-subset(MCC,MCC$fyear==i|MCC$fyear==(i+1))
  for (j in 1:10){
    xs<-subset(x,x$utrans==levels(x$utrans)[j])
    lmxs<-lm(asinsqrt~year,data=xs)
    year.slope[j,(i-1998)]<-coefficients(lmxs)[2]
  }
}
year.slope<-as.data.frame(year.slope)
year.slope<-stack(year.slope)
year.slope$year<-rep(c(1999:2011),each=10)
year.slope$utrans<-as.factor(rep(levels(MCC$utrans),13))
par(mfrow=c(3,4))
for (i in 1:10){
plot(year.slope$year[year.slope$utrans==levels(year.slope$utrans)[i]],year.slope$values[year.slope$utrans==levels(year.slope$utrans)[i]],main=levels(MC$utrans)[i],xlab="year",ylab="year to year slope (MC%cover)")
}
slope.slope<-rep(0,10)
yr<-seq(1999,2011,1)
par(mfrow=c(3,4))
for (i in 1:10){
  x<-subset(year.slope,year.slope$utrans==levels(year.slope$utrans)[i])
  test<-lm(values~year,data=x)
  slope.slope[i]<-coefficients(test)[2]
  plot(x$year,x$values,pch=18)
  lines(yr,coefficients(test)[1]+coefficients(test)[2]*yr)
}

slope.slope<-as.data.frame(slope.slope)
slope.slope$utrans<-as.factor(levels(MCC$utrans))
slope.slope$eig<-slopes$eig
par(mfrow=c(1,1))
plot(slope.slope$slope.slope,slope.slope$eig,col=cl[slope.slope$utrans],pch=18,xlab="slope of year to year slopes",ylab="IPM lambda")
legend("bottomleft",legend=sprintf(levels(slope.slope$utrans)),col=cl,pch=18)
cor.test(slope.slope$slope.slope,slope.slope$eig)
#r=-0.5397051, p-value=0.1074



slopes$ssp<-ssp
slopes.slopes$ssp<-ssp
plot(slopes$slopes,slopes$ssp,col=cl[slopes$utrans],pch=18,xlab="M. capitata % coverage slope transects 1999-2012",ylab="[(ssd*y*lambda)-(ssd*y)]/ssd")
cor.test(slopes$slopes,slopes$ssp)
legend("bottomright",legend=sprintf(levels(slopes$utrans)),col=cl,pch=18)
slopes$ssd<-ssl-ss
plot(slopes$slopes,slopes$ssd,col=cl[slopes$utrans],pch=18,xlab="M. capitata % coverage slope transects 1999-2012",ylab="(ssd*y*lambda)-(ssd*y)")
cor.test(slopes$slopes,slopes$ssd)
legend("topleft",legend=sprintf(levels(slopes$utrans)),col=cl,pch=18)
##What if we don't take the values in intermeidate years (between start in 1999 or 2000 and end in 2012) into account
## like t-test of start vs. end
MC1999<-subset(MCC,MCC$fyear=="1999")
MC2000<-subset(MCC,MCC$fyear=="2000")
MCstart<-rbind(MC1999[1:8,],MC2000[3:4,])
MC2012<-subset(MCC,MCC$fyear=="2012")
MCSE<-MC2012$utrans
MCSE<-as.data.frame(MCSE)
colnames(MCSE)<-c("utrans")
MCSE$pcover2012<-MC2012$pcover
MCSE$pcoverstart<-c(MCstart$pcover[1:4],MCstart$pcover[9:10],MCstart$pcover[5:8])
MCSE$asinsqrtstart<-c(MCstart$asinsqrt[1:4],MCstart$asinsqrt[9:10],MCstart$asinsqrt[5:8])
MCSE$asinsqrt2012<-c(MC2012$asinsqrt[1:4],MC2012$asinsqrt[9:10],MC2012$asinsqrt[5:8])
MCSE$eig<-eig
MCSE$asindiff<-(MCSE$asinsqrt2012-MCSE$asinsqrtstart)
par(mfrow=c(1,1))
plot(MCSE$asindiff,MCSE$eig,pch=18,col=cl[MCSE$utrans])
legend("bottomright", legend = sprintf(levels(MCSE$utrans)), col = cl, pch=18, xjust = 1, bg = "white")
MCSE$asinprop<-(MCSE$asinsqrt2012-MCSE$asinsqrtstart)/MCSE$asinsqrtstart
plot(MCSE$asinprop,MCSE$eig,pch=18,col=cl[MCSE$utrans])


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
slope.slope$ssp<-ssp
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