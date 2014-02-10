####Observed vs. "predicted" data for Non-linear growth model with linear intercept varying by site.
####growth functions by fateIPM
gg.sim=function(n,x,params) {
  sigmaxp2<-params[7,site]*exp(2*params[8,site]*x)
  sigmaxp<-sqrt(sigmaxp2)
  g<-rnorm(n,params[3,site]+params[4,site]*x+params[5,site]*exp((-1)*params[6,site]*x),sigmaxp)
  
}
ig.sim=function(n,x,params) {
  g<-rnorm(n,params[9,site]+params[10,site]*x,params[11,site])
}
ug.sim=function(n,x,params) {
  sigmaxp2<-params[14,site]*exp(2*params[15,site]*x)
  sigmaxp<-sqrt(sigmaxp2) 
  g<-rnorm(n,params[12,site]+params[13,site]*x,sigmaxp)
}

clo.sim=function(x,params) {
  nfiss<-exp(params[26,site]+params[27,site]*x)
}

chi.sim=function(x,params) {
  nfiss<-exp(params[28,site]+params[29,site]*x)
}

csize.sim=function(xp,x,params) {
  sigmaxp2<-params[24,site]*exp(2*params[25,site]*x)
  sigmaxp<-sqrt(sigmaxp2)
  c<-rnorm(n,params[22,site]+params[23,site]*x,sigmaxp)
}
cl<-c("cyan","mediumturquoise","dodgerblue","green3","purple","olivedrab2","orchid1","blue","darkblue","darkorchid4")



#######################Observed x or Size at t.  Simulated y or Size at t+1#######################
##Plot observed vs. simulated stasis growth data by transect
gsim<-list(vector("numeric",149),vector("numeric",29),vector("numeric",62),vector("numeric",72),vector("numeric",65),vector("numeric",87),vector("numeric",116),vector("numeric",73),vector("numeric",40),vector("numeric",108))
for (i in 1:10){
  x<-subset(growthonly,growthonly$utrans==levels(growthonly$utrans)[i])
  ggx<-x$size
  site=i
  n=length(ggx)
  gsim[[i]]<-gg.sim(n,ggx,params)
}
par(mfrow=c(3,4))
for (i in 1:10){
x<-subset(growthonly,growthonly$utrans==levels(growthonly$utrans)[i])
plot(x$size,x$sizeNext,main=levels(x$utrans)[i],xlim=c(-2,10),ylim=c(-2,10),xlab="Size t",ylab="Size t+1")
points(x$size,gsim[[i]],col="red",pch=18)
}
legend("bottomright", legend = c("observed","simulated"), col = c("black","red"), pch =c(1,18), xjust = 1, bg = "white")

##Plot observed vs. simulated Fission and Fusion growth data separately (These functions do not vary by transect)
igall<-ig.sim(136,fissonly$size,params)
par(mfrow=c(1,1))
plot(fissonly$size,igall,main="fission growth",xlab="Size t",ylab="Size t+1")
points(fissonly$size,fissonly$sizeNext,pch=18,col="red")
legend("bottomright", legend = c("observed","simulated"), col = c("black","red"), pch =c(1,18), xjust = 1, bg = "white")

ugall<-ug.sim(length(fuseonly$size),fuseonly$size,params)
plot(fuseonly$size,ugall,main="fusion growth",xlab="Size t",ylab="Size t+1")
points(fuseonly$size,fuseonly$sizeNext,pch=18,col="red")
legend("bottomright", legend = c("observed","simulated"), col = c("black","red"), pch =c(1,18), xjust = 1, bg = "white")

##Plot observed vs. simulated growth (Stasis, Fission and fusion together)
plot(growthonly$size,growthonly$sizeNext,main="All growth: simulated data",xlab="Size t",ylab="Size t+1")
points(fissonly$size,fissonly$sizeNext)
points(fuseonly$size,fuseonly$sizeNext)
for (i in 1:10){
  x<-subset(growthonly,growthonly$utrans==levels(growthonly$utrans)[i])
  points(x$size,gsim[[i]],col="black",pch=18)
}
points(fissonly$size,fissonly$sizeNext,pch=18,col="red")
points(fuseonly$size,fuseonly$sizeNext,pch=18,col="green")
legend("bottomright", legend = c("observed","simulated"), col = c("black","red"), pch =c(1,18), xjust = 1, bg = "white")


##### fate probabilities 
#plot observed vs. simulated probability of fissing (varies by transect)     
Pfisssim<-list(vector("numeric",312),vector("numeric",71),vector("numeric",161),vector("numeric",240),vector("numeric",144),vector("numeric",194),vector("numeric",400),vector("numeric",260),vector("numeric",73),vector("numeric",275))

for(i in 1:10){
  site=i
  x<-subset(MC,MC$utrans==levels(MC$utrans)[i])
  ggx<-x$size
  site=i
  Pfisssim[[i]]<-rbinom(n=length(ggx),size=1,prob=p.fiss.x(ggx,params))
}

par(mfrow=c(3,4))
for (i in 1:10){
  x<-subset(MC,MC$utrans==levels(MC$utrans)[i])
  plot(x$size,x$Pfiss,main=levels(x$utrans)[i],xlim=c(-2,10),ylim=c(0,1),xlab="Size t",ylab="Size t+1")
  points(x$size,Pfisssim[[i]],col="red",pch=18)
}
legend("bottomright", legend = c("observed","simulated"), col = c("black","red"), pch =c(1,18), xjust = 1, bg = "white")

#plot observed vs. simulated probability of fuseing (varies by transect)   
Pfusesim<-list(vector("numeric",312),vector("numeric",71),vector("numeric",161),vector("numeric",240),vector("numeric",144),vector("numeric",194),vector("numeric",400),vector("numeric",260),vector("numeric",73),vector("numeric",275))

for(i in 1:10){
  site=i
  x<-subset(MC,MC$utrans==levels(MC$utrans)[i])
  ggx<-x$size
  site=i
  Pfusesim[[i]]<-rbinom(n=length(ggx),size=1,prob=p.fuse.x(ggx,params))
}

par(mfrow=c(3,4))
for (i in 1:10){
  x<-subset(MC,MC$utrans==levels(MC$utrans)[i])
  plot(x$size,x$Pfuse,main=levels(x$utrans)[i],xlim=c(-2,10),ylim=c(0,1),xlab="Size t",ylab="Size t+1")
  points(x$size,Pfusesim[[i]],col="red",pch=18)
}
legend("bottomright", legend = c("observed","simulated"), col = c("black","red"), pch =c(1,18), xjust = 1, bg = "white")

#plot observed vs. simulated probability of being a large fuser (varies by transect) 
PfuseLGsim<-list(vector("numeric",312),vector("numeric",71),vector("numeric",161),vector("numeric",240),vector("numeric",144),vector("numeric",194),vector("numeric",400),vector("numeric",260),vector("numeric",73),vector("numeric",275))

for(i in 1:10){
  site=i
  x<-subset(MC,MC$utrans==levels(MC$utrans)[i])
  ggx<-x$size
  site=i
  PfuseLGsim[[i]]<-rbinom(n=length(ggx),size=1,prob=p.fuseLG.x(ggx,params))
}

par(mfrow=c(3,4))
for (i in 1:10){
  x<-subset(MC,MC$utrans==levels(MC$utrans)[i])
  plot(x$size,x$PfuseLG,main=levels(x$utrans)[i],xlim=c(-2,10),ylim=c(0,1),xlab="Size t",ylab="Size t+1")
  points(x$size,PfuseLGsim[[i]],col="red",pch=18)
}

legend("bottomright", legend = c("observed","simulated"), col = c("black","red"), pch =c(1,18), xjust = 1, bg = "white")

#plot observed vs. simulated probability of stasis (varies by transect) 
Pgrowthsim<-list(vector("numeric",312),vector("numeric",71),vector("numeric",161),vector("numeric",240),vector("numeric",144),vector("numeric",194),vector("numeric",400),vector("numeric",260),vector("numeric",73),vector("numeric",275))

for(i in 1:10){
  site=i
  x<-subset(MC,MC$utrans==levels(MC$utrans)[i])
  ggx<-x$size
  site=i
  Pgrowthsim[[i]]<-rbinom(n=length(ggx),size=1,prob=1-p.fiss.x(y,params=params) - p.fuse.x(y,params=params) + p.fiss.x(y,params=params)*p.fuse.x(y,params=params))
}

par(mfrow=c(3,4))
for (i in 1:10){
  x<-subset(MC,MC$utrans==levels(MC$utrans)[i])
  plot(x$size,x$Pgrowth,main=levels(x$utrans)[i],xlim=c(-2,10),ylim=c(0,1),xlab="Size t",ylab="Size t+1")
  points(x$size,Pgrowthsim[[i]],col="red",pch=18)
}
legend("bottomright", legend = c("observed","simulated"), col = c("black","red"), pch =c(1,18), xjust = 1, bg = "white")


####Survival####

#plot observed vs. simulated probability of survival (varies by transect) 
Psurvsim<-list(vector("numeric",312),vector("numeric",71),vector("numeric",161),vector("numeric",240),vector("numeric",144),vector("numeric",194),vector("numeric",400),vector("numeric",260),vector("numeric",73),vector("numeric",275))

for(i in 1:10){
  site=i
  x<-subset(MC,MC$utrans==levels(MC$utrans)[i])
  ggx<-x$size
  site=i
  Psurvsim[[i]]<-rbinom(n=length(ggx),size=1,prob=s.x(y,params=params))
}
####Clone number

#probability of lo fiss product number  
par(mfrow=c(1,1))
Philosim<-rbinom(n=length(MCf$size),size=1,prob=p.hilo.x(MCf$size,params))
plot(MCf$size,MCf$Philo,main="Probability of Low # of fission Products",xlim=c(-2,10),ylim=c(0,1),xlab="Size t",ylab="Size t+1")
points(MCf$size,Philosim,col="red",pch=18)
legend("bottomright", legend = c("observed","simulated"), col = c("black","red"), pch =c(1,18), xjust = 1, bg = "white")


#Fission offspring using separate vglm pospoisson function for hi vs lo number of fission products


lonumsim<-rpospois(n=length(lonum$size),clo.sim(lonum$size,params))
plot(lonum$size,lonum$numfiss,main="# of clones lo",xlim=c(-2,10),xlab="Size t",ylab="# of fission products")
points(lonum$size,lonumsim,col="red",pch=18)
legend("bottomright", legend = c("observed","simulated"), col = c("black","red"), pch =c(1,18), xjust = 1, bg = "white")



hinumsim<-rpospois(n=length(hinum$size),chi.sim(hinum$size,params))
plot(hinum$size,hinum$numfiss,main="# of clones hi",xlim=c(-2,10),xlab="Size t",ylab="# of fission products")
points(hinum$size,hinumsim,col="red",pch=18)
legend("bottomright", legend = c("observed","simulated"), col = c("black","red"), pch =c(1,18), xjust = 1, bg = "white")

par(mfrow=c(1,1))
hist(fisssizedat$sizeNext,freq=F,main="Size Distribution of Clonal Recruits",xlab="Size t+1",xlim=c(-2,10),ylim=c(0,0.5))
hist(csize.sim(100000,fisssizedat$size,params),freq=F,add=T,col=rainbow(1,alpha=0.25))
legend("bottomright", legend = c("observed","simulated"), col = c("black","red"), pch =15, xjust = 1, bg = "white")

#recruitment = sexual offspring 
par(mfrow=c(1,1))
hist(recsizedat$sizeNext,freq=F,main="Size Distribution of Sexual Recruits",xlab="Size t+1",xlim=c(-2,10),ylim=c(0,0.5))
hist(rnorm(100000,params[32,site],params[33,site]),freq=F,add=T,col=rainbow(1,alpha=0.25))
legend("bottomright", legend = c("observed","simulated"), col = c("black","red"), pch =15, xjust = 1, bg = "white")


#######Simulated data across full range (i.e. y or -2 to 10) of IPM kernel

####growth functions by fateIPM
plot(0,0,xlim=c(-2,10),ylim=c(-2,10),type="n",main="Stasis Growth",xlab="Size t",ylab="Size t+1")

for (i in 1:10){
  site=i
  points(y,gg.sim(300,y,params),col=cl[i],pch=18,cex=0.5)
}
legend("bottomright", legend = levels(growthonly$utrans), col = cl, pch =18, xjust = 1, bg = "white")

plot(0,0,xlim=c(-2,10),ylim=c(-2,10),type="n",main="Fission&Fusion Growth",xlab="Size t",ylab="Size t+1")
points(y,ig.sim(300,y,params),col="red",pch=18,cex=0.7)
points(y,ug.sim(300,y,params),col="green",pch=18,cex=0.7)
legend("bottomright", legend = c("fission","fusion"), col = c("red","green"), pch =18, xjust = 1, bg = "white")

plot(MCg$size,MCg$sizeNext,col=c("red","green","black")[MCg$fate],main="Stasis,Fission&Fusion Growth: Observed Data",pch=18,cex=0.7,xlab="Size t",ylab="Size t+1")
legend("bottomright", legend = c("stasis","fission","fusion"), col = c("black","red","green"), pch =18, xjust = 1, bg = "white")

##simulated data across range of full kernel = y = -2:10 Plot stasis, fission fusion growth together
plot(0,0,xlim=c(-2,10),ylim=c(-2,10),type="n",main="Stasis,Fission&Fusion Growth: Simulated data",xlab="Size t",ylab="Size t+1")
for (i in 1:10){
  site=i
  points(y,gg.sim(300,y,params),pch=18,cex=0.5)
}
points(y,ig.sim(300,y,params),col="red",pch=18,cex=0.7)
points(y,ug.sim(300,y,params),col="green",pch=18,cex=0.7)
legend("bottomright", legend = c("stasis","fission","fusion"), col = c("black","red","green"), pch =18, xjust = 1, bg = "white")

##### fate probabilities 
par(mfrow=c(1,1))

plot(0,0,xlim=c(-2,10),ylim=c(0,1),type="n",ylab="Probability of Fission",xlab="Size t")

for(i in 1:10){
  site=i
  points(y,rbinom(n=length(y),size=1,prob=p.fiss.x(y,params)),pch=18,cex=0.7)
}


plot(0,0,xlim=c(-2,10),ylim=c(0,1),type="n",ylab="Probability of Fusion",xlab="Size t")

for(i in 1:10){
  site=i
  points(y,rbinom(n=length(y),size=1,prob=p.fuse.x(y,params)),pch=18,cex=0.5)
}

plot(0,0,xlim=c(-2,10),ylim=c(0,1),type="n",ylab="Probability of Large Fusion",xlab="Size t")
for(i in 1:10){
  site=i
  points(y,rbinom(n=length(y),size=1,prob=p.fuseLG.x(y,params)),pch=18,cex=0.5)
}

plot(0,0,xlim=c(-2,10),ylim=c(0,1),type="n",ylab="Probability of growth",xlab="Size t")
for(i in 1:10){
  site=i
  points(y,rbinom(n=length(y),size=1,prob=1-p.fiss.x(y,params=params) - p.fuse.x(y,params=params) + p.fiss.x(y,params=params)*p.fuse.x(y,params=params)),pch=18,cex=0.5)
}
 

####Clone number

#probability of lo fiss product number  
plot(0,0,xlim=c(-2,10),ylim=c(0,1),type="n",ylab="Probability of low # of fusion products",xlab="Size t")
points(y,rbinom(n=length(y),size=1,prob=p.hilo.x(y,params)))


#Fission offspring using separate vglm pospoisson function for hi vs lo number of fission products
plot(0,0,xlim=c(-2,10),ylim=c(0,10),type="n",ylab="# of fission products",xlab="Size t")
points(y,rpospois(n=length(y),clo.sim(y,params)),col="green",pch=18)
points(y,rpospois(n=length(y),chi.sim(y,params)),col="red",pch=18)
legend("topleft", legend = c("High","Low"), col = c("red","green"), pch =18, xjust = 1, bg = "white")

par(mfrow=c(1,1))
hist(fisssizedat$sizeNext,freq=F,main="Size Distribution of Clonal Recruits",xlab="Size t+1",xlim=c(-2,10),ylim=c(0,0.5))
hist(csize.sim(100000,y,params),freq=F,add=T,col=rainbow(1,alpha=0.25))
legend("bottomright", legend = c("observed","simulated"), col = c("black","red"), pch =15, xjust = 1, bg = "white")

################Combine Probability and Growth/Clonal ##############

####growth functions by fateIPM

PGsim<-rep.int(0,3000)
for (i in 1:10){
  site=i
  start<-(i-1)*300+1
  end<-start+299
  PGsim[start:end]<-gg.sim(300,y,params)*rbinom(n=length(y),size=1,prob=1-p.fiss.x(y,params=params) - p.fuse.x(y,params=params) + p.fiss.x(y,params=params)*p.fuse.x(y,params=params))
}
PGsim<-as.data.frame(PGsim)
PGsim$yval<-rep.int(y,10)
PGsim<-subset(PGsim,PGsim!=0)


plot(0,0,xlim=c(-2,10),ylim=c(-2,10),type="n",main="Stasis Growth",xlab="Size t",ylab="Size t+1")
points(PGsim$yval,PGsim$PGsim,col="black",pch=18,cex=0.5)


PIsim<-rep.int(0,300)
PIsim<-ig.sim(300,y,params)*rbinom(n=length(y),size=1,prob=p.fiss.x(y,params))
PIsim<-as.data.frame(PIsim)
PIsim$yval<-y
PIsim<-subset(PIsim,PIsim!=0)
plot(0,0,xlim=c(-2,10),ylim=c(-2,10),type="n",main="Fission&Fusion Growth",xlab="Size t",ylab="Size t+1")
points(PIsim$yval,PIsim$PIsim,col="red",pch=18,cex=0.7)


PUsim<-rep.int(0,3000)
for (i in 1:10){
  site=i
  start<-(i-1)*300+1
  end<-start+299
  PUsim[start:end]<-ug.sim(300,y,params)*rbinom(n=length(y),size=1,prob=p.fuse.x(y,params)*p.fuseLG.x(y,params))
}
PUsim<-as.data.frame(PUsim)
PUsim$yval<-rep.int(y,10)
PUsim<-subset(PUsim,PUsim!=0)

points(PUsim$yval,PUsim$PUsim,col="green",pch=18,cex=0.7)
legend("bottomright", legend = c("fission","fusion"), col = c("red","green"), pch =18, xjust = 1, bg = "white")

plot(0,0,xlim=c(-2,10),ylim=c(-2,10),type="n",main="Stasis,Fission&Fusion Growth:simulated data with fate prob ",xlab="Size t",ylab="Size t+1")
points(PIsim$yval,PIsim$PIsim,col="red",pch=18,cex=0.7)
points(PGsim$yval,PGsim$PGsim,col="black",pch=18,cex=0.5)
points(PUsim$yval,PUsim$PUsim,col="green",pch=18,cex=0.7)
legend("bottomright", legend =c("growth","fusion","fission") , col = c("black","green","red"), pch =18, xjust = 1, bg = "white")


####Clone number

#probability of lo fiss product number  
points(y,rbinom(n=length(y),size=1,prob=p.hilo.x(y,params)))


#Fission offspring using separate vglm pospoisson function for hi vs lo number of fission products
Closim<-rep.int(0,300)
Closim<-rpospois(n=length(y),clo.sim(y,params))*rbinom(n=length(y),size=1,prob=p.hilo.x(y,params))
Closim<-as.data.frame(Closim)
Closim$yval<-y
Closim<-subset(Closim,Closim!=0)

Chisim<-rep.int(0,300)
Chisim<-rpospois(n=length(y),chi.sim(y,params))*(1-rbinom(n=length(y),size=1,prob=p.hilo.x(y,params)))
Chisim<-as.data.frame(Chisim)
Chisim$yval<-y
Chisim<-subset(Chisim,Chisim!=0)
plot(0,0,xlim=c(-2,10),ylim=c(0,10),type="n",ylab="# of fission products",xlab="Size t")
points(Closim$yval,Closim$Closim,col="green",pch=18)
points(Chisim$yval,Chisim$Chisim,col="red",pch=18)
legend("topleft", legend = c("High","Low"), col = c("red","green"), pch =18, xjust = 1, bg = "white")



par(mfrow=c(1,1))
hist(growthonly$sizeNext,freq=F,main="Size t+1 Stasis growth",xlab="Size t+1",xlim=c(-2,10),ylim=c(0,0.5))
hist(PGsim$PGsim,freq=F,add=T,col=rainbow(1,alpha=0.25))
legend("right", legend = c("observed","simulated"), col = c("black","red"), pch =15, xjust = 1, bg = "white")

hist(fissonly$sizeNext,freq=F,main="Size t+1 Fission growth",xlab="Size t+1",xlim=c(-2,10),ylim=c(0,0.5))
hist(PIsim$PIsim,freq=F,add=T,col=rainbow(1,alpha=0.25))
legend("right", legend = c("observed","simulated"), col = c("black","red"), pch =15, xjust = 1, bg = "white")


hist(fuseonly$sizeNext,freq=F,main="Size t+1 Fusion growth",xlab="Size t+1",xlim=c(-2,10),ylim=c(0,0.5))
hist(PUsim$PUsim,freq=F,add=T,col=rainbow(1,alpha=0.25))
legend("right", legend = c("observed","simulated"), col = c("black","red"), pch =15, xjust = 1, bg = "white")

PSsim<-list(vector("numeric",14),vector("numeric",3),vector("numeric",7),vector("numeric",5),vector("numeric",1),vector("numeric",3),vector("numeric",7),vector("numeric",5),vector("numeric",4),vector("numeric",5))

for (i in 1:10){
  x<-subset(GOES,GOES$utrans==levels(GOES$utrans)[i])
  ggx<-x$size
  site=i
  start<-(i-1)*300+1
  end<-start+299
  PSsim[[i]]<-gg.sim(length(ggx),ggx,params)
}
PSsim<-c(PSsim[[1]],PSsim[[2]],PSsim[[3]],PSsim[[4]],PSsim[[5]],PSsim[[6]],PSsim[[7]],PSsim[[8]],PSsim[[9]],PSsim[[10]])

hist(GOES$sizeNext,freq=F,main="Size t+1 Extreme Shrinkage nonlinear",xlab="Size t+1",xlim=c(-2,10),ylim=c(0,0.5))
hist(PSsim,freq=F,add=T,col=rainbow(1,alpha=0.5))

####If you want to repeat
####growth functions by fateIPM
par(mfrow=c(3,4))
z<-0
plot(0,0,xlim=c(-2,10),ylim=c(-2,10),type="n")

for (i in 1:10){
  site=i
  repeat{
    points(y,gg.sim(300,y,params),col=cl[i])
    z<-z+1
    print(z)
    if(z>10) break()
  }
}

repeat{
  points(y,ig.sim(300,y,params),col="red")
  points(y,ug.sim(300,y,params),col="green")
  z<-z+1
  print(z)
  if(z>10) break()
}
plot(MCg$size,MCg$sizeNext,col=c("red","green","black")[MCg$fate])

##### fate probabilities 

z<-0
plot(0,0,xlim=c(-2,10),ylim=c(0,1),type="n")

for(i in 1:10){
  site=i
  repeat{
    points(y,rbinom(n=length(y),size=1,prob=p.fiss.x(y,params)))
    z<-z+1
    print(z)
    if(z>10) break()
  }
}

z<-0
plot(0,0,xlim=c(-2,10),ylim=c(0,1),type="n")

for(i in 1:10){
  site=i
  repeat{
    points(y,rbinom(n=length(y),size=1,prob=p.fuse.x(y,params)))
    z<-z+1
    print(z)
    if(z>10) break()
  }
}

for(i in 1:10){
  site=i
  repeat{
    points(y,rbinom(n=length(y),size=1,prob=p.fuseLG.x(y,params)))
    z<-z+1
    print(z)
    if(z>10) break()
  }
}

for(i in 1:10){
  site=i
  repeat{
    points(y,rbinom(n=length(y),size=1,prob=1-p.fiss.x(y,params=params) - p.fuse.x(y,params=params) + p.fiss.x(y,params=params)*p.fuse.x(y,params=params)))
    z<-z+1
    print(z)
    if(z>10) break()
  }
}
