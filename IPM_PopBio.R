###############################################################Lambda,SSD,RV,Damping ratio, Elasticity matrix###############################################################
## To get lambda for each site
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

##to get elasticitie for each site

eigsensP<-array(0,c(nx,ncoef))
eigsensN<-array(0,c(nx,ncoef))
#sensN<-array(0, c(n, n, nx,ncoef))
#sensP<-array(0, c(n, n, nx,ncoef))

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

elasP <- (eig-eigsensP)/eig/(0.1)
par(mfrow=c(2,2))
par(mar = c(7, 2, 2, 2) + 0.1)
for (i in c(1:ncoef)){
  barplot(elasP[,i],horiz=TRUE,xlab="Elasticity",col=cl[1:nx],main=paste(rownames(params)[i],":Elasticity"))#,xlim=c(-2.5,0.3))
}
#may not want xlim.  Good to understand relative elasticities of different params but not so good for utrans comparison
plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
legend("center",legend=sprintf(levels(MC$utrans)),col=cl,pch=15)


par(mfrow = c(2, 2))
par(mar = c(7, 2, 2, 2) + 0.1)
for (i in c(1:nx)){
  barplot(elasP[i,],main=paste(levels(MC$utrans)[i],":Elasticity"),names.arg=rownames(params),las=2,ylim=c(min(elasP),max(elasP)))
}

par(mfrow = c(3, 1))
par(mar=c(.5,.5,2,.5))
par(oma=c(10,4,1,1))
barplot(elasP[2,],main=levels(MC$utrans)[2],las=2,ylim=c(min(elasP),max(elasP)),names.arg=rownames(params),axisnames=FALSE)
barplot(elasP[8,],main=levels(MC$utrans)[8],las=2,ylim=c(min(elasP),max(elasP)),names.arg=rownames(params),axisnames=FALSE)
barplot(elasP[9,],main=levels(MC$utrans)[9],las=2,ylim=c(min(elasP),max(elasP)),names.arg=rownames(params))

mtext("elasticity of lambda to Vital rate parameters", side = 2, outer = TRUE, cex = 1, line = 2.2,col = "grey20")


