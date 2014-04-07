cl1<-topo.colors(10)

#Plot Kernel by site
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
mtext("Size t (log(Area))", side = 1, outer = TRUE, cex = 1, line = 3,col = "grey20")
mtext("Size t+1 (log (Area))", side = 2, outer = TRUE, cex = 1, line = 3,col = "grey20")
mtext("Transition Probabilities", side = 4, outer = TRUE, cex = 1, line = 5,col = "grey20")

par(mfrow=c(1,1))
par(mar=c(1,1,1,1),oma=c(1,1,1,5))
library(aqfig)
vertical.image.legend(col=topo.colors(300),zlim=c(mi,ma))
#Plot single Kernel
par(oma=c(6,6,1,10))
site=6
image(y,y,t(K), col=topo.colors(300),zlim=c(mi,ma),xlab="",ylab="")
mtext("Size (t)", side = 1, outer = TRUE, cex = 1, line = 3,col = "grey20")
mtext("Size (t+1)", side = 2, outer = TRUE, cex = 1, line = 3,col = "grey20")
mtext("Transition Probabilities", side = 4, outer = TRUE, cex = 1, line = 5,col = "grey20")
par(oma=c(6,6,1,8))
vertical.image.legend(col=topo.colors(300),zlim=c(mi,ma))

##########################################################################
#plot average curves for survival, growth, Pfiss/Pfuse/growth, # of fission offspring
par(mfrow=c(2,2))
par(oma=c(1,3,1,3))
par(mar=c(4,4,3,2))

#Growth
plot(GONGS$size,GONGS$sizeNext,pch=18,cex=0.5, xlim=c(-2,8),ylim=c(-2,8),xlab="Size t (log(Area)",ylab="Size t+1 (log(Area))",col=cl1[3])
points(fissonly$size,fissonly$sizeNext, pch=18, cex=0.5, col=cl1[6])
points(fuseonly$size,fuseonly$sizeNext, pch=18, cex=0.5, col=cl1[1])
points(GOES$size,GOES$sizeNext, pch=18, cex=0.5, col=cl1[4])
x0<-seq(from=-2.4,to=10.5,by=0.01)#Range of values for lines added to plots
lines(x0,fixed.effects(growth.reg.re.a0)[1]+fixed.effects(growth.reg.re.a0)[2]*x0+fixed.effects(growth.reg.re.a0)[3]*exp(-1*fixed.effects(growth.reg.re.a0)[4]*x0),col=cl1[3],lwd=2)
lines(x0,fiss.int+fiss.slope*x0,col=cl1[6],lwd=2)
lines(x0,fuse.int+fuse.slope*x0,col=cl1[1],lwd=2)
lines(x0,es.int+es.slope*x0,col=cl1[4],lwd=2)
legend("bottomright", legend = c("normal","fission","fusion","ES"), col = cl1[c(3,6,1,4)], lwd = 2, xjust = 1, bg = "white",bty="n")
mtext(letters[1],side=3,line=1,adj=0.01)

#Fate probability
plot(MC$size, MC$Pfiss,type="n",xlim=c(-2,8),xlab="Size t (log(Area))",ylab="Fate Probabilities")
pf<-exp(fixed.effects(fiss.re0)[1]+fixed.effects(fiss.re0)[2]*x0)
pfu<-exp(fixed.effects(fuse.re0)[1]+fixed.effects(fuse.re0)[2]*x0)
pful<-exp(fixed.effects(fuseLG.re0)[1]+fixed.effects(fuseLG.re0)[2]*x0)
pe<-exp(coefficients(Pes))
pe1<-(1-pf/(pf+1) - pfu/(pfu+1) + pf/(pf+1)*pfu/(pfu+1))*(pe/(pe+1))
pg<-(1-pf/(pf+1) - pfu/(pfu+1) + pf/(pf+1)*pfu/(pfu+1))*(1-pe/(pe+1))
lines(x0,pf/(pf+1),col=cl1[6],lwd=2)
lines(x0,pfu/(pfu+1)*pful/(pful+1),col=cl1[1],lwd=2)
lines(x0,pe1,col=cl1[4],lwd=2)
lines(x0,pg,col=cl1[3],lwd=2)
mtext(letters[2],side=3,line=1,adj=0.01)
#Survival
plot(MC$size,MC$Psurv,pch=18,cex=0.5,xlim=c(-2,8),xlab="Size t (log(Area))",ylab="Survival Probability")
u<-exp(fixed.effects(surv.reg.re0)[1]+surv.slope.re0*x0)
lines(x0,u/(u+1),lwd=2)
mtext(letters[3],side=3,line=1,adj=0.01)

# of fission products
plot(MCf$size,MCf$numfiss,pch=18,cex=0.5,xlim=c(-2,8),ylim=c(0,10),xlab="Size t (log(Area))",ylab="# of clonal offspring per fission event")
u1<-exp(clonenum1.int+clonenum1.slope*x0)
u2<-exp(clonenum2.int+clonenum2.slope*x0)
u3<-exp(clonenum3.int+clonenum3.slope*x0)
uall<-1+u1+u2+u3
p4<-1-u1/uall-u2/uall-u3/uall
nfiss<-u1/uall+u2/uall*2+u3/uall*3+p4*70/11#70/11 is weighted average of # of offspring over 4.
lines(x0,nfiss,lwd=2)
mtext(letters[4],side=3,line=1,adj=0.01)

####################################################################################
#Eigen Values
big<-array(0, c(n, n, nx))
eig<-vector("numeric", nx)

for (site in 1:nx){
  big[,,site]<-Kernel(y,n, params)
  A<-eigen.analysis(big[,,site])
  eig[site]<-A$lambda1
}

eig
par(mar = c(7, 4, 4, 2) + 0.1)
xvals<-barplot(eig,ylab="Population Growth Rate (lambda)")
text(xvals, par("usr")[3] - 0.05, srt = 45, adj = 1,labels = unique(MC$utrans), xpd = TRUE)
###############################################################################
##Parameter Elasticities

eigsensP<-array(0,c(nx,ncoef))
eigsensN<-array(0,c(nx,ncoef))

for (site in 1:nx){
  for (j in 1:ncoef) {
    svecP<-params
    svecP[j,]<-1.1*svecP[j,]
    sensP<-Kernel(y,n, params=svecP)
    eigsensP[site,j]<-eigen.analysis(sensP)$lambda1
    svecN<-params
    svecN[j,]<-0.9*svecN[j,]
    sensN<-Kernel(y,n, params=svecN)
    eigsensN[site,j]<-eigen.analysis(sensN)$lambda1
  }
}
#increase by 10%
elasP <- (eig-eigsensP)/eig/(0.1)
colnames(elasP)<-rownames(params)
rownames(elasP)<-levels(MC$utrans)
#decrease by 10%
#increase by 10%
elasN <- (eig-eigsensN)/eig/(0.1)
colnames(elasN)<-rownames(params)
rownames(elasN)<-levels(MC$utrans)


#Plots for All parameter Elasticities by site
par(mfrow = c(2, 2))
par(mar = c(7, 2, 2, 2) + 0.1)
for (i in c(1:nx)){
  barplot(elasP[i,],main=paste(levels(MC$utrans)[i],":Elasticity"),names.arg=rownames(params),las=2,ylim=c(min(elasP),max(elasP)))
  
}

par(mfrow = c(2, 2))
par(mar = c(7, 2, 2, 2) + 0.1)
for (i in c(1:nx)){
  barplot(elasN[i,],main=paste(levels(MC$utrans)[i],":Elasticity"),names.arg=rownames(params),las=2,ylim=c(min(elasN),max(elasN)))
  
}

#Plots for All parameter elasticities by parameter
par(mfrow = c(2, 2))
par(mar = c(7, 2, 2, 2) + 0.1)
for (i in c(1:ncoef)){
  barplot(elasP[,i],main=colnames(elasP)[i],names.arg=rownames(elasP),las=2,ylim=c(min(elasP),max(elasP)))
  
}

#############################################################################33
#Find top 5 parameter elasticities
rankelas<-matrix(0,32,10)
colnames(rankelas)<-levels(MC$utrans)
rownames(rankelas)<-rownames(params)
for (i in 1:10){
  rankelas[,i]<-rank(abs(elasP[i,]))
}

rankelasN<-matrix(0,32,10)
colnames(rankelasN)<-levels(MC$utrans)
rownames(rankelasN)<-rownames(params)
for (i in 1:10){
  rankelasN[,i]<-rank(abs(elasN[i,]))
}

#Plot for highest ranked parameters for all sites
par(mfrow=c(3,4))
par(mar=c(2.5,.5,2,.5))
par(oma=c(6,5,1,6))
for(i in 1:10){
  sn<-c("Honolua N","Honolua S","Kahekili 3m","Kahekili 7m","Molokini 7m","Molokini 13m","Olowalu 3m","Olowalu 7m","Papaula 10m","Puamana 3m")
  barplot(elasP[i,c(12,
corelas1<-matrix(0,32,2)
rownames(corelas1)<-rownames(params)
colnames(corelas1)<-c("r","pvalue")
par(mfrow=c(3,4))
for (i in 1:32){
  x<-cor.test(slopes$slopes,elasP[,i])
  corelas1[i,1]<-x$estimate
  corelas1[i,2]<-x$p.value
  plot(slopes$slopes,elasP[,i],main=rownames(corelas1)[i])
}

bleh<-subset(corelas1,corelas1[,2]<0.05)
par(mfrow=c(3,4))
for (i in 1:32){
  plot(slopes$posneg,elasP[,i],main=rownames(corelas1)[i],pch=18, xlim=c(-2,2))
}
#no significant relationship between4,15,25,14,18,7,13)],main=sn[i],las=2,xlim=c(min(elasP),max(elasP)),horiz=TRUE,col=cl,axisnames=FALSE)
  #if(i %in% c(1,5,9))barplot(elasP[i,c(12,4,15,25,14,18,7,13)],main=sn[i],las=2,xlim=c(min(elasP),max(elasP)),horiz=TRUE,axes=FALSE)
  #if(i %in% c(2,3,4,6,7,8,10))barplot(elasP[i,c(12,4,15,25,14,18,7,13)],main=sn[i],las=2,xlim=c(min(elasP),max(elasP)),horiz=TRUE,axisnames=FALSE,axes=FALSE)
  #if(i %in% c(7,8,9,10))axis(1)
}
subparam<-c(expression(italic(gi)[1]),expression(italic(gn)[sigma]),expression(italic("pi")[0]),expression(italic(gu)[sigma]),expression(italic(cn1)[0]),expression(italic(gu)[0]),expression(italic(gn)[1]),expression(italic(gi)[0]))
legend(3,10, legend = subparam, col = cl[8:1], pch=15, xjust = 1, bg = "white",bty="n",xpd=NA,cex=1.25)
mtext(expression(paste("Elasticity of ",lambda, " to parameter estimates")), side = 1, outer = TRUE, cex = 1, line = 3,col = "grey20")

###############################################################################################
#look for correlation between lambda and parameter elasticities
corelas<-matrix(0,32,2)
rownames(corelas)<-rownames(params)
colnames(corelas)<-c("r","pvalue")
par(mfrow=c(3,4))
for (i in 1:32){
  x<-cor.test(eig,elasP[,i])
  corelas[i,1]<-x$estimate
  corelas[i,2]<-x$p.value
  plot(eig,elasP[,i],col=cl1,main=rownames(corelas)[i],pch=18)
}
corelas<-as.data.frame(corelas)
sigonly.lam<-subset(corelas,corelas[,2]<0.05)
#plot statistically significant relationships only
par(mfrow=c(3,4))
for (i in c(12,4,15,14,18,7,13,3,11,16,19,22,23,31)){
  plot(eig,elasP[,i],col=cl1,main=rownames(corelas)[i],pch=18)
}

###############################################################################################
#correlation between parameter elasticities and cover
AllC<- read.csv("C:/Users/Megan/IPMcoral/MauiCover.csv")
AllC<-na.omit(AllC)
AllC$trans<-paste(AllC$Site,AllC$Depth..m.)
AllC$trans<-as.factor(AllC$trans)
index<-levels(AllC$trans)
values<-c(levels(MC$utrans)[1:4],levels(MC$utrans)[6],levels(MC$utrans)[5],levels(MC$utrans)[7:10])
AllC$utrans<-values[match(AllC$trans,index)]
AllC$utrans<-as.factor(AllC$utrans)
AllC$fyear<-as.factor(AllC$year)
MCC<- subset(AllC,AllC$Taxon=="Montipora capitata")
MCC<-na.omit(MCC)
MCC$pcover<-MCC$cover/100
MCC$asinsqrt<-asin(sqrt(MCC$pcover))

#Plot trajectories in M. capitata cover from 1999-2012
par(mfrow=c(3,4))
for(i in 1:10){
  x<-subset(MCC,MCC$utrans==levels(MCC$utrans)[i])
  plot(x$year,x$cover,pch=18,main=levels(MCC$utrans)[i],xlab="year",ylab="M. capitata cover")
  bleh<-lm(cover~year,data=x)
  dyear<-data.frame(year=x$year)
  lines(x$year,predict(bleh,newdata=dyear))
}

#Determine whether slope is significant at alpha=0.5
slopesig<-matrix(0,10,2)
rownames(slopesig)<-levels(MCC$utrans)
colnames(slopesig)<-c("slope","pvalue")
for(i in 1:10){
  x<-subset(MCC,MCC$utrans==levels(MCC$utrans)[i])
  bleh<-lm(cover~year,data=x)
  slopesig[i,1]<-coefficients(bleh)[2]
  slopesig[i,2]<-anova(bleh)[1,5]
}
slopesig<-as.data.frame(slopesig)
sigonly.slopes<-subset(slopesig,slopesig$pvalue<0.05)
#HonoluaS, Kahekili07, Molokini 13m, Olowalu 03 and Papula 10m have significant relationship with time

###M. capitata cover slope over 1999-2012 vs. lambda
slopes<-matrix(0,10,4)
for (i in 1:10){
  x<-subset(MCC,MCC$utrans==levels(MCC$utrans)[i])
  overall<-lm(asinsqrt~year,data=x)
  slopes[i,1]<-coefficients(overall)[2]
  slopes[i,2]<-anova(overall)[1,5]
}
slopes<-as.data.frame(slopes)
colnames(slopes)<-c("slopes","pvalue","lambda","utrans")
slopes$eig<-eig
slopes$utrans<-as.factor(levels(MC$utrans))

#no signiricant relationship between cover and eigen values

#Try again using only areas with significant changes in cover
corelas1<-matrix(0,32,2)
rownames(corelas1)<-rownames(params)
colnames(corelas1)<-c("r","pvalue")
par(mfrow=c(3,4))
for (i in 1:32){
  x<-cor.test(slopes1$slopes,elasP[c(2,4,6,7,9),i])
  corelas1[i,1]<-x$estimate
  corelas1[i,2]<-x$p.value
  plot(slopes1$slopes,elasP[c(2,4,6,7,9),i],main=rownames(corelas1)[i])
}
bleh<-subset(corelas1,corelas1[,2]<0.05)

#significant relationship with Pfiss.slope, PfuseLg.int, r
#####################################################################
#Try with truncated CRAMP transect data
#Plot trajectories in M. capitata cover from 1999-2005
MCCtrunc<-subset(MCC,MCC$fyear=="1999"|MCC$fyear=="2000"|MCC$fyear=="2001"|MCC$fyear=="2002"|MCC$fyear=="2003"|MCC$fyear=="2004"|MCC$fyear=="2005")
MCCtrunc<-droplevels(MCCtrunc)
par(mfrow=c(3,4))
for(i in 1:10){
  x<-subset(MCCtrunc,MCCtrunc$utrans==levels(MCCtrunc$utrans)[i])
  plot(x$year,x$cover,pch=18,main=levels(MCCtrunc$utrans)[i],xlab="year",ylab="M. capitata cover")
}

par(mfrow=c(3,4))
for(i in 1:10){
  x<-subset(MCCtrunc,MCCtrunc$utrans==levels(MCCtrunc$utrans)[i])
  plot(x$year,x$cover,pch=18,main=levels(MCCtrunc$utrans)[i],xlab="year",ylab="M. capitata cover")
  bleh<-lm(cover~year,data=x)
  dyear<-data.frame(year=x$year)
  lines(x$year,predict(bleh,newdata=dyear))
}

#Determine whether slope is significant at alpha=0.5
slopesig<-matrix(0,10,2)
rownames(slopesig)<-levels(MCCtrunc$utrans)
colnames(slopesig)<-c("slope","pvalue")
for(i in 1:10){
  x<-subset(MCCtrunc,MCCtrunc$utrans==levels(MCCtrunc$utrans)[i])
  bleh<-lm(cover~year,data=x)
  slopesig[i,1]<-coefficients(bleh)[2]
  slopesig[i,2]<-anova(bleh)[1,5]
}
slopesig<-as.data.frame(slopesig)
sigonly.slopes<-subset(slopesig,slopesig$pvalue<0.05)

#Only significant change for Molokini 13m

###M. capitata cover slope over 1999-2005 vs. lambda
slopes<-matrix(0,10,4)
for (i in 1:10){
  x<-subset(MCCtrunc,MCCtrunc$utrans==levels(MCCtrunc$utrans)[i])
  overall<-lm(asinsqrt~year,data=x)
  slopes[i,1]<-coefficients(overall)[2]
  slopes[i,2]<-anova(overall)[1,5]
}
slopes<-as.data.frame(slopes)
colnames(slopes)<-c("slopes","pvalue","lambda","utrans")
slopes$eig<-eig
slopes$utrans<-as.factor(levels(MC$utrans))

corelas1<-matrix(0,32,2)
rownames(corelas1)<-rownames(params)
colnames(corelas1)<-c("r","pvalue")
par(mfrow=c(3,4))
for (i in 1:32){
  x<-cor.test(slopes$slopes,elasP[,i])
  corelas1[i,1]<-x$estimate
  corelas1[i,2]<-x$p.value
  plot(slopes$slopes,elasP[,i],main=rownames(corelas1)[i])
}

bleh<-subset(corelas1,corelas1[,2]<0.05)
#Significant relationship between cover and es.sd, fiss.int, fiss.slope


#no significant relationship between cover and elasticities
#no signiricant relationship between cover and eigen values

#Try again using only areas with significant changes in cover
corelas1<-matrix(0,32,2)
rownames(corelas1)<-rownames(params)
colnames(corelas1)<-c("r","pvalue")
par(mfrow=c(3,4))
for (i in 1:32){
  x<-cor.test(slopes1$slopes,elasP[c(2,4,6,7,9),i])
  corelas1[i,1]<-x$estimate
  corelas1[i,2]<-x$p.value
  plot(slopes1$slopes,elasP[c(2,4,6,7,9),i],main=rownames(corelas1)[i])
}
bleh<-subset(corelas1,corelas1[,2]<0.05)

#significant relationship with Pfiss.slope, PfuseLg.int, r

##########################################################################
#Elasticities to recruitment by orders of magnitude
svecP<-params
svecP[32,]<-1.11*svecP[32,]
sensP<-Kernel(y,n, params=svecP)
eigsensPrec<-eigen.analysis(sensP)$lambda1
svecN<-params
svecN[32,]<-0.89*svecN[32,]
sensN<-Kernel(y,n, params=svecN)
eigsensNrec<-eigen.analysis(sensN)$lambda1

elasPrec<-(eig-eigsensPrec)/eig/(0.1)
elasRP<-elasP
elasRP[,32]<-elasPrec

elasNrec<-(eig-eigsensNrec)/eig/(0.1)
elasRN<-elasN
elasRN[,32]<-elasNrec
#Plot All parameter Elasticities
par(mfrow = c(2, 2))
par(mar = c(7, 2, 2, 2) + 0.1)
for (i in c(1:nx)){
  barplot(elasRP[i,],main=paste(levels(MC$utrans)[i],":Elasticity"),names.arg=rownames(params),las=2,ylim=c(min(elasRP),max(elasRP)))
  
}
par(mfrow = c(2, 2))
par(mar = c(7, 2, 2, 2) + 0.1)
for (i in c(1:nx)){
  barplot(elasRN[i,],main=paste(levels(MC$utrans)[i],":Elasticity"),names.arg=rownames(params),las=2,ylim=c(min(elasRN),max(elasRN)))
  
}
#Find top 5 parameter elasticities
rankelasRP<-matrix(0,32,10)
colnames(rankelasRP)<-levels(MC$utrans)
rownames(rankelasRP)<-rownames(params)
for (i in 1:10){
  rankelasRP[,i]<-rank(abs(elasRP[i,]))
}

rankelasRN<-matrix(0,32,10)
colnames(rankelasRN)<-levels(MC$utrans)
rownames(rankelasRN)<-rownames(params)
for (i in 1:10){
  rankelasRN[,i]<-rank(abs(elasRN[i,]))
}
