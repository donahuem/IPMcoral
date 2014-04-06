#####################################################################Plot Vital Rate Regressions#######################################################

#http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
cl<-c("cyan","mediumturquoise","dodgerblue","green3","purple","olivedrab2","orchid1","blue","darkblue","darkorchid4")
x<-seq(from=-2.4,to=10.5,by=0.01)#Range of values for lines added to plots

##Survival##
plot(dat$size,dat$Psurv,col=cl[dat$utrans],type="p",xlab="Size (t)",ylab="Survival Probability")
for (site in 1:(nx)) {
  lines(y,s.x(y,params=params),col=cl[site],lwd=2)
}
legend("right", legend = sprintf(levels(growdat$utrans)), col = cl, lwd = 2, xjust = 1, bg = "white")

##Growth##
#Just to show why we include fate as an explanatory variable
plot(growdat$size,growdat$sizeNext,col=c(cl[8],cl[3],cl[5])[growdat$fate],pch=18,cex=0.95,ylab="Size (t+1)",xlab="Size (t)")
fategrow<-lm(sizeNext~size+fate,data=growdat)
lines(x,coefficients(fategrow)[1]+coefficients(fategrow)[2]*x,col=cl[8],lwd=2)
lines(x,coefficients(fategrow)[1]+coefficients(fategrow)[3]+coefficients(fategrow)[2]*x,col=cl[3],lwd=2)
lines(x,coefficients(fategrow)[1]+coefficients(fategrow)[4]+coefficients(fategrow)[2]*x,col=cl[5],lwd=2)
legend("bottomright", legend = sprintf(levels(growdat$fate)), col=c(cl[8],cl[3],cl[5]), lty=1,lwd=2, xjust = 1, bg = "white")

# Growth functions for stasis, fission, and fusion by site
par(mfrow=c(3,1),mar=c(5,2,2,10)+0.1,xpd=F)
plot(growthonly$size,growthonly$sizeNext,col=cl[growthonly$utrans],type="p",pch=18,cex=0.65,main="Stasis",xlab="Size (t)",ylab="Size (t+1)",ylim=c(-2,10),xlim=c(-2,10))
for (site in 1:(nx)) {
  lines(x,params[3,site]+params[4,site]*x+params[5,site]*exp(-params[6,site]*x),col=cl[site],lwd=1)
}
legend("bottomright", legend = sprintf(levels(growdat$utrans)), col = cl, lwd=2, xjust = 1, bg = "white")

plot(fissonly$size,fissonly$sizeNext,col=cl[fissonly$utrans],type="p",pch=18,cex=0.65,main="Fission",xlab="Size (t)",ylab="Size (t+1)",ylim=c(-2,10),xlim=c(-2,10))
  lines(x,params[9,site]+params[10,site]*x,lwd=1)
legend("bottomright", legend = sprintf(levels(growdat$utrans)), col = cl, pch=18, xjust = 1, bg = "white")

plot(fuseonly$size,fuseonly$sizeNext,col=cl[fuseonly$utrans],type="p",pch=18,cex=0.65,main="Fusion",xlab="Size (t)",ylab="Size (t+1)",ylim=c(-2,10),xlim=c(-2,10))
lines(x,params[12,site]+params[13,site]*x,lwd=1)
legend("bottomright", legend = sprintf(levels(growdat$utrans)), col = cl, pch=18, xjust = 1, bg = "white")



#Density functions for growth by fate
par(mfrow=c(3,4))
for (site in 1:nx){
  plot(MC$Psurv~MC$size, type="n")
  lines(y,gg.yx(y,y,params),col="red")
  lines(y,ig.yx(y,y,params),col="blue")
  lines(y,ug.yx(y,y,params),col="green")
}

##Fate Probabilities##
par(mfrow=c(1,1))
#Probability of growth#
plot(dat$size,dat$Pgrowth,col=cl[dat$utrans],type="p",ylim=c(0,1),xlab="Size (t)",ylab="Prob(growth)",xlim=c(-2,10))
for (site in 1:(nx)) {
  lines(x,1-p.fiss.x(x,params=params) - p.fuse.x(x,params=params) + p.fiss.x(x,params=params)*p.fuse.x(x,params=params),col=cl[site],lwd=2)
}
legend(15,1), legend = sprintf(levels(growdat$utrans)), col = cl, lwd = 2, xjust = 1, bg = "white")
##Probability of fission##
plot(dat$size,dat$Pfiss,col=cl[dat$utrans],type="p",ylim=c(0,1),xlim=c(-2,10),xlab="Size (t)",ylab="Prob(fission)")
for (site in 1:(nx)) {
  lines(x,p.fiss.x(x,params=params)-p.fiss.x(x,params=params)*p.fuse.x(x,params=params),col=cl[site],lwd=2)
}
legend(13.2,1, legend = sprintf(levels(growdat$utrans)), col = cl, lwd=2, xjust = 1, bg = "white")

#Probability of Large fusion#
plot(dat$size,dat$Pfuse,col=cl[dat$utrans],type="p",ylim=c(0,1),xlim=c(-2,10),xlab="Size (t)",ylab="Prob(fusion)")
for (site in 1:(nx)) {
  lines(x,p.fuse.x(x,params=params)*p.fuseLG.x(x,params=params),col=cl[site],lwd=2)
}
legend(13.2,1, legend = sprintf(levels(growdat$utrans)), col = cl, lwd=2, xjust = 1, bg = "white")

#Fate probabilities using code from IPM.R
par(mfrow=c(3,4))
for (site in 1:nx){
  S=s.x(y,params=params) # survival 
  I=p.fiss.x(y,params=params)-p.fiss.x(y,params=params)*p.fuse.x(y,params=params) #fisser-only = I
  U=(p.fuse.x(y,params=params)*p.fuseLG.x(y,params=params))  #fusion of U AND Ms:  include Ms in fusion growth
  Gr=1-p.fiss.x(y,params=params) - p.fuse.x(y,params=params) + p.fiss.x(y,params=params)*p.fuse.x(y,params=params)
  plot(MC$Psurv~MC$size)
  lines(y,S,col="black")
  lines(y,I,col="blue")
  lines(y,U,col="green")
  lines(y,Gr,col="red")
}

#Number of fission products#
# Several ways: zero truncated poisson, negative binomial, geometric, negative binomial not truncated, multinomial probability of 1-4 products and 1-10 products
dsize<-data.frame(size=x)
clonenum<-vglm(numfiss~size,data=fissdat,family=pospoisson)
clonenb<- glm.nb(numfiss~size,data=MCf)
clonenbz<-vglm(numfiss~size,data=fissdat,family=posnegbinomial)
clonegeo<-zerotrunc(numfiss~size,data=fissdat,dist="geometric")

plot(fissdat$size,fissdat$numfiss,col=cl[fissdat$utrans],type="p",pch=18,xlab="Size (t)",ylab="# of fission products/'parent colony'",ylim=c(0,10),xlim=c(-2,10))
lines(x,predict(clonenum,newdata=dsize, "response"), col="red")
lines(x,predict(clonenb,newdata=dsize, "response"),col="blue")#for negative binomial distribution
lines(x,predict(clonenbz,newdata=dsize, "response"))
lines(x,predict(clonegeo,newdata=dsize, "response"),col="green")# for geometric distribution
site=1#same for all sites
lines(x,p.c1.x(x,params)+p.c2.x(x,params)*2+p.c3.x(x,params)*3+(1-p.c1.x(x,params)-p.c2.x(x,params)-p.c3.x(x,params))*4,col="magenta")
lines(x,p2.c1.x(x,params)+p2.c2.x(x,params)*2+p2.c3.x(x,params)*3+p2.c4.x(x,params)*4+p2.c6.x(x,params)*6+p2.c7.x(x,params)*7+p2.c10.x(x,params)*10+(1-p2.c1.x(x,params)-p2.c2.x(x,params)-p2.c3.x(x,params)-p2.c4.x(x,params)-p2.c6.x(x,params)-p2.c7.x(x,params)-p2.c10.x(x,params))*5,col="purple")#for multinomial distribution
legend("topleft", legend = c("ztpoisson","negbin","ztnegbin","ztgeometric","multi4","multi10"), col = c("red","blue","black","green","magenta","purple"), lwd = 2, xjust = 1, bg = "white")

#Plot from page 127 of notes
library(plyr)#for mean by group = utrans

#Plot one point per site with average number of fission products
clone=glm(Pfiss~utrans,data=MC,family=binomial)
Pfiss.mean<-c(coefficients(clone)[1],coefficients(clone)[2:10])
for (i in 1:10){
  Pfiss.mean[i]<-exp(Pfiss.mean[i])/(1+exp(Pfiss.mean[i]))}
Pfiss.mean
numfiss.mean<-ddply(MCf,~utrans,summarise,mean=mean(numfiss))
numfiss.mean
bleh<-cbind(numfiss.mean,Pfiss.mean)
plot(bleh$mean,bleh$Pfiss.mean,col=cl[bleh$utrans],pch=18,xlab="mean # of fission products",ylab="mean P(fiss)")
legend("topright", legend = sprintf(levels(bleh$utrans)), col = cl[bleh$utrans], pch=18, xjust = 1, bg = "white")

#plot of probability of fission by number of fission product.  Several points per site (number of points matches number of fission events per site)
bleh<-data.frame(MCf$utrans,MCf$numfiss)
index<-levels(MCf$utrans)
Pfiss.mean<-as.data.frame(Pfiss.mean)
values<-Pfiss.mean[,1]
bleh$Pfiss.mean<-values[match(bleh$MCf.utrans,index)]
plot(bleh$MCf.numfiss,bleh$Pfiss.mean,col=cl[bleh$MCf.utrans],pch=18,xlab="# of fission products",ylab="Mean P(fiss) by site")
legend("topright", legend = sprintf(levels(bleh$MCf.utrans)), col = cl, pch=18, xjust = 1, bg = "white")

##size distribution of clones##
hist(fisssizedat$sizeNext,ylim=c(0,1),freq=FALSE)
image(y,y,t(Cmulti), col=topo.colors(300),main=sn,axes=TRUE)
points(fisssizedat$size,fisssizedat$sizeNext,pch=16,col="white")

#Recruit size distribution
par(mfrow=c(3,4))
x<-seq(-2,3,by=0.01)
with(MCrs[MCrs$utrans==levels(MCrs$utrans)[1],],hist(sizeNext,freq=FALSE,ylim=c(0,1),xlim=c(-2,3)))
lines(x,dnorm(x,coefficients(rec.sizeU)[1],summary(rec.sizeU)$sigma))
for (i in 2:10){
  with(MCrs[MCrs$utrans==levels(MCrs$utrans)[i],],hist(sizeNext,freq=FALSE,ylim=c(0,1),xlim=c(-2,3)))
  lines(x,dnorm(x,coefficients(rec.sizeU)[1]+coefficients(rec.sizeU)[i],summary(rec.sizeU)$sigma))
}
#Colony density (for recruitment) size distribution
par(mfrow=c(3,3))
x<-seq(-2,10,by=0.01)
with(MCds[MCds$utrans==levels(MCds$utrans)[1],],hist(sizeNext,freq=FALSE,ylim=c(0,1),xlim=c(-2,10)))
lines(x,dnorm(x,coefficients(col.sizeU)[1],summary(col.sizeU)$sigma))
for (i in 2:10){
  with(MCds[MCds$utrans==levels(MCds$utrans)[i],],hist(sizeNext,freq=FALSE,ylim=c(0,1),xlim=c(-2,10)))
  lines(x,dnorm(x,coefficients(col.sizeU)[1]+coefficients(col.sizeU)[i],summary(col.sizeU)$sigma))
}

##Plot kernel
maxP<-matrix(nrow=1,ncol=nx)
par(mfrow=c(3,4))
for(site in 1:nx){
  K<-Kernel(y,n,params)
  sn<-levels(dat$utrans)[site]
  image(y,y,t(K), col=topo.colors(300),main=sn,axes=TRUE)
  vertical.image.legend(col=topo.colors(n),zlim=c(0,max(abs(K))))
  maxP[1,site]<-max(K)
}
mtext("Size (t)", side = 1, outer = TRUE, cex = 1, line = 2.2,col = "grey20")
mtext("Size (t+1)", side = 2, outer = TRUE, cex = 1, line = 2.2,col = "grey20")
par(mfrow=c(1,1))
par(mar=c(4,1,1,1),oma=c(1,1,1,1))
vertical.image.legend(col=topo.colors(n),zlim=c(0,0.05)

par(mfrow=c(3,4))
for (site in 1:10){
  sn<-levels(dat$utrans)[site]
  K<-Kernel(y,n,params)
  image(y,y,t(K), col=topo.colors(300),main=sn,axes=TRUE,xlab="Size (t)",ylab="Size (t+1)")
  vertical.image.legend(col=topo.colors(n),zlim=c(min(K),max(K)))
}
                      
mins<-rep(0,10)
maxs<-rep(0,10)
for (site in 1:nx){
  K<-Kernel(y,n,params)
  maxs[site]<-max(K)
  mins[site]<-min(K)
}
ma<-max(maxs)
mi<-min(mins)
par(mfrow=c(3,4))
par(mar=c(.5,.5,2,.5))
par(oma=c(6,6,1,10))

for (site in 1:10){
  sn<-c("Honolua N","Honolua S","Kahekili 3m","Kahekili 7m","Molokini 7m","Molokini 13m","Olowalu 3m","Olowalu 7m","Papaula 10m","Puamana 3m")
  K<-Kernel(y,n,params)
  image(y,y,t(K), col=topo.colors(300),zlim=c(mi,ma),main=sn[site],xlab="",ylab="",axes=FALSE)
  if (site %in% c(7,8,9,10)) axis(1)
  if (site %in% c(1,5,9)) axis(2)
}
mtext("Size (t)", side = 1, outer = TRUE, cex = 1, line = 3,col = "grey20")
mtext("Size (t+1)", side = 2, outer = TRUE, cex = 1, line = 3,col = "grey20")
mtext("Transition Probabilities", side = 4, outer = TRUE, cex = 1, line = 5,col = "grey20")

par(mfrow=c(1,1))
par(mar=c(1,1,1,1),oma=c(1,1,1,5))
library(aqfig)
vertical.image.legend(col=topo.colors(300),zlim=c(mi,ma))
                      
#Elemental elasticities                      

mins<-rep(0,10)
maxs<-rep(0,10)
for (site in 1:nx){
  K<-Kernel(y,n,params)
  A<-eigen.analysis(K)
  maxs[site]<-max(A$elasticities)
  mins[site]<-min(A$elasticities)
}
ma<-max(maxs)
mi<-min(mins)     
par(mfrow=c(3,4))
par(mar=c(.5,.5,2,.5))
par(oma=c(6,6,1,10))
for (site in 1:nx){
  sn<-c("Honolua N","Honolua S","Kahekili 3m","Kahekili 7m","Molokini 7m","Molokini 13m","Olowalu 3m","Olowalu 7m","Papaula 10m","Puamana 3m")
  K<-Kernel(y,n,params)
  A<-eigen.analysis(K)
  image(y,y,t(A$elasticities),col=topo.colors(300),main=sn[site],xlab="",ylab="",zlim=c(mi,ma),axes=FALSE)
  if (site %in% c(7,8,9,10)) axis(1)
  if (site %in% c(1,5,9)) axis(2)
}
mtext("Size (t)", side = 1, outer = TRUE, cex = 1, line = 2.2,col = "grey20")
mtext("Size (t+1)", side = 2, outer = TRUE, cex = 1, line = 2.2,col = "grey20")
mtext("Elasticities", side = 4, outer = TRUE, cex = 1, line = 5,col = "grey20")
par(mfrow=c(1,1))
par(mar=c(1,1,1,1),oma=c(1,1,1,5))
vertical.image.legend(col=topo.colors(300),zlim=c(mi,ma))                   
                      
                      
                      
                      
par(mfrow=c(3,4))
par(mar=c(.5,.5,2,.5))
par(oma=c(6,6,1,10))
for (i in 1:10){
  x<-MC$size[MC$utrans==levels(MC$utrans)[i]]
  
  site=i
  K<-Kernel(y,n,params)
  A<-eigen.analysis(K)
  plot(y,A$stable.stage,type="l",yaxt="n",xaxt="n",ylim=c(0,0.012),xlim=c(-2,10),main=levels(MC$utrans)[i])
  lines(y,A$stable.stage)
  if (site %in% c(7,8,9,10)) axis(1)
  if (site %in% c(1,5,9)) axis(2)
  }
mtext("Size (t)", side = 1, outer = TRUE, cex = 1, line = 2.2,col = "grey20")
mtext("Proportional Distribution ", side = 2, outer = TRUE, cex = 1, line = 2.2,col = "grey20")

                      
                      
                      
                      
par(mfrow=c(2,2))
for (i in c(2:3,5:7,9:11,16:18,20,21,28:31,33,34)){
plot(eig$eig,elasP[,i],ylab=rownames(params)[i],col=cl,pch=18)
}
