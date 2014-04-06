MCC1 <- read.csv("C:/Users/Megan/IPMcoral/MCC1.csv")

slopesig<-matrix(0,6,2)
rownames(slopesig)<-levels(MCC1$Site1)
colnames(slopesig)<-c("slope","p-value")
for(i in 1:6){
  x<-subset(MCC1,MCC1$Site1==levels(MCC1$Site1)[i])
  bleh<-lm(cover~year,data=x)
  slopesig[i,1]<-coefficients(bleh)[2]
  slopesig[i,2]<-anova(bleh)[1,5]
}
View(slopesig)
par(mfrow=c(3,4))
diff<-c(-0.079313826,  1.627486354,	-1.717392042,	-0.941133964,	-1.84334922,	-7.97070708,	2.385712354,	3.94997946,	6.623062202,	0.199279688)
diffelas<-matrix(0,32,2)
rownames(diffelas)<-rownames(params)
colnames(diffelas)<-c("r","p-value")
for (i in 1:32){
  x<-cor.test(slopes$slopes,elasP[,i])
  diffelas[i,1]<-x$estimate
  diffelas[i,2]<-x$p.value
 plot(diff,elasP[,i],main=rownames(diffelas)[i],pch=18)
}



MCCtrunc<-subset(MCC,MCC$fyear=="1999"|MCC$fyear=="2000"|MCC$fyear=="2001"|MCC$fyear=="2002"|MCC$fyear=="2003"|MCC$fyear=="2004"|MCC$fyear=="2005")
MCCtrunc<-droplevels(MCCtrunc)
par(mfrow=c(3,4))
for(i in 1:10){
  x<-subset(MCCtrunc,MCCtrunc$utrans==levels(MCCtrunc$utrans)[i])
  plot(x$year,x$cover,pch=18,main=levels(MCCtrunc$utrans)[i],xlab="year",ylab="M. capitata cover")
}

###M. capitata cover slope over 1999-2012 vs. lambda

slopes<-matrix(0,10,4)
for (i in 1:10){
  x<-subset(MCCtrunc,MCCtrunc$utrans==levels(MCCtrunc$utrans)[i])
  overall<-lm(asinsqrt~year,data=x)
  slopes[i,1]<-coefficients(overall)[2]
  slopes[i,2]<-anova(overall)[1,5]
}
slopes<-as.data.frame(slopes)
colnames(slopes)<-c("slopes","p-value","lambda","utrans")
slopes$eig<-eig
slopes$utrans<-as.factor(levels(MC$utrans))
corelas1<-matrix(0,32,2)
rownames(corelas1)<-rownames(params)
colnames(corelas1)<-c("r","p-value")
par(mfrow=c(3,4))
for (i in 1:32){
  x<-cor.test(slopes$slopes,elasP[,i])
  corelas1[i,1]<-x$estimate
  corelas1[i,2]<-x$p.value
  plot(slopes$slopes,elasP[,i],main=rownames(corelas1)[i])
}

bleh<-subset(corelas1,corelas1[,2]<0.05)



corparams<-matrix(0,32,2)
rownames(corparams)<-rownames(params)
colnames(corparams)<-c("r","pvalue")
for (i in 1:32){
plot(eig,params[i,],pch=18, col=cl,main=rownames(params)[i])
x<-cor.test(eig,params[i,])
corparams[i,1]<-x$estimate
corparams[i,2]<-x$p.value
}






par(mfrow=c(1,2))
plot(elasP[,12],slopes$slopes,xlab="",ylab="slope of change in coral cover 1999-2005",pch=18,cex=2)
x<-lm(slopes$slopes~elasP[,12])
x0<-seq(min(elasP[,12])-1,max(elasP[,12])+1,0.1)
lines(x0,coefficients(x)[1]+coefficients(x)[2]*x0,lwd=2)
mtext(letters[1],side=3,line=0.5,adj=0.01)

plot(elasP[,12],slopes$eig,xlab="",ylab=expression(paste("multiplicative growth rate (  ",lambda,")")),pch=18,cex=2)
x<-lm(slopes$eig~elasP[,12])
lines(x0,coefficients(x)[1]+coefficients(x)[2]*x0,lwd=2)
mtext(letters[2],side=3,line=0.5,adj=0.01)

mtext(expression(paste("Elasticity of ",lambda, " to 10% increase in ", italic(gi)[0])), side = 1, outer = TRUE, cex = 1, line = -2,col = "grey20")
