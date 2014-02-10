###replacement for histogram plot for fission growth with fate probability
plot(0,0,type="n",xlab="size(t)",ylab="cummulative probability",main="Probability of fission",xlim=c(-2,10),ylim=c(0,0.2))
test<-matrix(0,nrow=300,ncol=2)
test[,1]<-y

prac<-subset(MC,MC$fateIPM=="growth"|MC$fateIPM=="fissLG"|MC$fateIPM=="fuseLG"|MC$fateIPM=="fuseSM"|MC$fateIPM=="M")
for (i in 1:300){
  lx<-subset(prac,prac$size<y[i])
  fn<-subset(lx,lx$fateIPM=="fissLG")
  test[i,2]<-length(fn$size)/1223
  
}

Pfisssim<-list(vector("numeric",197),vector("numeric",37),vector("numeric",88),vector("numeric",125),vector("numeric",100),vector("numeric",140),vector("numeric",219),vector("numeric",156),vector("numeric",46),vector("numeric",115))
Pfisssim2<-list(vector("numeric",197),vector("numeric",37),vector("numeric",88),vector("numeric",125),vector("numeric",100),vector("numeric",140),vector("numeric",219),vector("numeric",156),vector("numeric",46),vector("numeric",115))
for(i in 1:10){
  site=i
  x<-subset(prac,prac$utrans==levels(MC$utrans)[i])
  Pfisssim[[i]]<-rnorm(length(x$size),mean(x$size),sd(x$size))
  Pfisssim2[[i]]<-rbinom(n=length(x$size),size=1,prob=p.fiss.x(Pfisssim[[i]],params))
}

test2<-matrix(0,nrow=300,ncol=2)
test2[,1]<-y
bleh<-matrix(0,nrow=1223,ncol=2)
bleh[,1]<-c(Pfisssim[[1]],Pfisssim[[2]],Pfisssim[[3]],Pfisssim[[4]],Pfisssim[[5]],Pfisssim[[6]],Pfisssim[[7]],Pfisssim[[8]],Pfisssim[[9]],Pfisssim[[10]])
bleh[,2]<-c(Pfisssim2[[1]],Pfisssim2[[2]],Pfisssim2[[3]],Pfisssim2[[4]],Pfisssim2[[5]],Pfisssim2[[6]],Pfisssim2[[7]],Pfisssim2[[8]],Pfisssim2[[9]],Pfisssim2[[10]])

for (i in 1:300){
  lx<-subset(bleh,bleh[,1]<y[i])
  fn<-subset(lx,lx[,2]==1)
  test2[i,2]<-length(fn[,2])/1223
  
}

lines(test2[,1],test2[,2],col="red")
lines(test[,1],test[,2],lwd=2)
####################################################################################
plot(0,0,type="n",xlab="size(t)",ylab="cummulative probability",main="Probability of fusion",xlim=c(-2,10),ylim=c(0,0.2))
test<-matrix(0,nrow=300,ncol=2)
test[,1]<-y

prac<-subset(MC,MC$fateIPM=="growth"|MC$fateIPM=="fissLG"|MC$fateIPM=="fuseLG"|MC$fateIPM=="fuseSM"|MC$fateIPM=="M")
for (i in 1:300){
  lx<-subset(prac,prac$size<y[i])
  fn<-subset(lx,lx$fateIPM=="fuseLG")
  test[i,2]<-length(fn$size)/1223
  
}

Pfusesim<-list(vector("numeric",197),vector("numeric",37),vector("numeric",88),vector("numeric",125),vector("numeric",100),vector("numeric",140),vector("numeric",219),vector("numeric",156),vector("numeric",46),vector("numeric",115))
Pfusesim2<-list(vector("numeric",197),vector("numeric",37),vector("numeric",88),vector("numeric",125),vector("numeric",100),vector("numeric",140),vector("numeric",219),vector("numeric",156),vector("numeric",46),vector("numeric",115))
for(i in 1:10){
  site=i
  x<-subset(prac,prac$utrans==levels(MC$utrans)[i])
  Pfusesim[[i]]<-rnorm(length(x$size),mean(x$size),sd(x$size))
  Pfusesim2[[i]]<-rbinom(n=length(x$size),size=1,prob=p.fuse.x(Pfusesim[[i]],params)*p.fuseLG.x(Pfusesim[[i]],params))
}

test2<-matrix(0,nrow=300,ncol=2)
test2[,1]<-y
bleh<-matrix(0,nrow=1223,ncol=2)
bleh[,1]<-c(Pfusesim[[1]],Pfusesim[[2]],Pfusesim[[3]],Pfusesim[[4]],Pfusesim[[5]],Pfusesim[[6]],Pfusesim[[7]],Pfusesim[[8]],Pfusesim[[9]],Pfusesim[[10]])
bleh[,2]<-c(Pfusesim2[[1]],Pfusesim2[[2]],Pfusesim2[[3]],Pfusesim2[[4]],Pfusesim2[[5]],Pfusesim2[[6]],Pfusesim2[[7]],Pfusesim2[[8]],Pfusesim2[[9]],Pfusesim2[[10]])

for (i in 1:300){
  lx<-subset(bleh,bleh[,1]<y[i])
  fn<-subset(lx,lx[,2]==1)
  test2[i,2]<-length(fn[,2])/1223
  
}

lines(test2[,1],test2[,2],col="red")

lines(test[,1],test[,2],lwd=2)
###############################################################################3

plot(0,0,type="n",xlab="size(t)",ylab="cummulative probability",main="Probability of ES Growth",xlim=c(-2,10),ylim=c(0,0.5))
test<-matrix(0,nrow=300,ncol=2)
test[,1]<-y

prac<-subset(MC,MC$fateIPM=="growth"|MC$fateIPM=="fissLG"|MC$fateIPM=="fuseLG"|MC$fateIPM=="fuseSM"|MC$fateIPM=="M")
for (i in 1:300){
  lx<-subset(prac,prac$size<y[i])
  f<-subset(growthonly,growthonly$ES=="ES")
  fn<-subset(f,f$size<y[i])
  test[i,2]<-length(fn$size)/1223
  
}

PSsim<-matrix(0,1223,2)
PSsim[,1]<-rnorm(length(prac$size),mean(prac$size),sd(prac$size))
PSsim[,2]<-rbinom(n=length(prac$size),size=1,prob=(1-p.fiss.x(PSsim[[i]],params=params) - p.fuse.x(PSsim[[i]],params=params) + p.fiss.x(PSsim[[i]],params=params)*p.fuse.x(PSsim[[i]],params=params)*p.es.x(PSsim[[i]],params=params)))


test2<-matrix(0,nrow=300,ncol=2)
test2[,1]<-y
bleh<-matrix(0,nrow=1223,ncol=2)
bleh[,1]<-c(PSsim[[1]],PSsim[[2]],PSsim[[3]],PSsim[[4]],PSsim[[5]],PSsim[[6]],PSsim[[7]],PSsim[[8]],PSsim[[9]],PSsim[[10]])
bleh[,2]<-c(PSsim2[[1]],PSsim2[[2]],PSsim2[[3]],PSsim2[[4]],PSsim2[[5]],PSsim2[[6]],PSsim2[[7]],PSsim2[[8]],PSsim2[[9]],PSsim2[[10]])

for (i in 1:300){
  lx<-subset(PSsim,PSsim[,1]<y[i])
  fn<-subset(lx,lx[,2]==1)
  test2[i,2]<-length(fn[,2])/1223
  
}

lines(test2[,1],test2[,2],col="red")
lines(test[,1],test[,2],lwd=2)




E<-resid(growth.reg,type="normalized")
par(mfrow=c(1,1))
plot(GONGS$size,E,xlim=c(-2,10),xlab="Size (t)",ylab="normalized residuals")
site=1
x<-seq(-2,20,0.1)
sigmaxp2<-params[7,site]*exp(2*params[8,site]*x)
sigmaxp<-sqrt(sigmaxp2)
lines(x,sigmaxp)
lines(x,-1*sigmaxp)
