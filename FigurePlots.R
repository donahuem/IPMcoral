library(aqfig)
cl1<-topo.colors(11)
cl<-c("cyan","mediumturquoise","dodgerblue","green3","purple","olivedrab2","orchid1","blue","darkblue","darkorchid4")
sn<-c("Mean","Honolua N","Honolua S","Kahekili 3m","Kahekili 7m","Molokini 7m","Molokini 13m","Olowalu 3m","Olowalu 7m","Papaula 10m","Puamana 3m")

#Plot single Kernel

site=1 #can change site for different site's plot
K<-Kernel(y,n,params)
par(mar=c(5,4,4,2)+0.1)
par(oma=c(1,1,1,7))
image(y,y,t(K), col=topo.colors(300),xlab="Size t (log(Area))",ylab="Size t+1 (log(Area))")
mtext("Transition Probabilities", side = 4, outer = TRUE, cex = 1, line = 4,col = "grey20")
par(oma=c(1,1,1,5))
vertical.image.legend(col=topo.colors(300),zlim=c(min(K),max(K)))

#Elemental elasticities                      
A<-eigen.analysis(K)
par(oma=c(1,1,1,7))
image(y,y,t(A$elasticities),col=topo.colors(300),xlab="Size t (log(Area))",ylab="Size t+1 (log(Area))")
mtext("Elasticities", side = 4, outer = TRUE, cex = 1, line = 4,col = "grey20")
par(oma=c(1,1,1,5))
vertical.image.legend(col=topo.colors(300),zlim=c(min(A$elasticities),max(A$elasticities)))     


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
#Get Eigen Values
big<-array(0, c(n, n, nx))
eig<-vector("numeric", nx)

for (site in 1:nx){
  big[,,site]<-Kernel(y,n, params)
  A<-eigen.analysis(big[,,site])
  eig[site]<-A$lambda1
}

###############################################################################
##Parameter Elasticities
#'P' is positive perturbation 'N' is negative perturbation
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
elasP <- (eig-eigsensP)/eig/(-0.1)#(a/da)*(deig/eig)= (a/(a-1.1*a))*(eig-eigsensP)/eig
colnames(elasP)<-rownames(params)
rownames(elasP)<-sn

#decrease by 10%
elasN <- (eig-eigsensN)/eig/(0.1)#10*(eig-eigsensN)/eig
colnames(elasN)<-rownames(params)
rownames(elasN)<-sn

#############################################################################33
#Rank parameter elasticities
rankelas<-matrix(0,32,11)
colnames(rankelas)<-sn
rownames(rankelas)<-rownames(params)
for (i in 1:nx){
  rankelas[,i]<-rank(abs(elasP[i,]))
}

rankelasN<-matrix(0,32,11)
colnames(rankelasN)<-sn
rownames(rankelasN)<-rownames(params)
for (i in 1:nx){
  rankelasN[,i]<-rank(abs(elasN[i,]))
}
##############Plot 10 parameter elasticities with highest absolute values for mean IPM############################
site=1 #can change site for different site's plot
K<-Kernel(y,n,params)
par(mfrow=c(1,1))
par(oma=c(3,2,1,6))
barplot(elasP[site,c(12,4,15,24,3,25,14,7,13,6)],las=2,xlim=c(min(elasP),max(elasP)),horiz=TRUE,col=cl,axisnames=FALSE)
subparam<-c(expression(italic(g)['3,'][italic(n)]),expression(italic(i)[1]),expression(italic(g)[sigma]),expression(italic(i)[sigma]),expression(italic(c)['1,1']),expression(italic(g)['0,'][italic(n)]),expression(italic(c)['0,1']),expression(italic(u)[0]),expression(italic(g)['1,'][italic(n)]),expression(italic(i)[0]))
legend(2.5,11, legend = subparam, col = cl[10:1], pch=15, xjust = 1, bg = "white",bty="n",xpd=NA,cex=1)
mtext(expression(paste("Elasticity of ",lambda, " to parameter estimates")), side = 1, cex = 1,line=4, col = "grey20")


#####################################################################################################################
#Recruitment sensitivity 
#increase and decrease order of magnitude above and below observed
##Parameter Elasticities
#'P' is positive perturbation 'N' is negative perturbation
recsensP<-array(0,c(nx,1))
recsensN<-array(0,c(nx,1))

for (site in 1:nx){
    svecP<-params
    svecP[32,]<-10*svecP[32,]
    sensP<-Kernel(y,n, params=svecP)
    recsensP[site]<-eigen.analysis(sensP)$lambda1
    svecN<-params
    svecN[32,]<-0.1*svecN[32,]
    sensN<-Kernel(y,n, params=svecN)
    recsensN[site]<-eigen.analysis(sensN)$lambda1
}

#increase by 10%
recelasP <- (eig-recsensP)/eig/(-9)#(a/da)*(deig/eig)= (a/(a-1.1*a))*(eig-eigsensP)/eig
rownames(recelasP)<-sn

#decrease by 10%
recelasN <- (eig-recsensN)/eig/(0.9)#10*(eig-eigsensN)/eig
rownames(recelasN)<-sn
############################################################
##Elasticity of lambda to perturbation of functions
#could not figure out how to loop this so I did it manually.  
#Increased and decreased function by 10% in Kernel function
#probably should not include the probability of NGS because it is just a combination of other parameters.
#Then ran this code to get elasticity

funcsensP<-array(0,c(nx,nx))
funcsensN<-array(0,c(nx,nx))

for (site in 1:nx){
    sensP<-Kernel(y,n, params=params)
    funcsensP[site,func]<-eigen.analysis(sensP)$lambda1
}
for (site in 1:nx){
  sensN<-Kernel(y,n, params=params)
  funcsensN[site,func]<-eigen.analysis(sensN)$lambda1
}
#increase by 10%
funcelasP <- (eig-funcsensP)/eig/(-0.1)#(a/da)*(deig/eig)= (a/(a-1.1*a))*(eig-eigsensP)/eig
colnames(funcelasP)<-c("S","I","U","GR","ES","GS","GG","GI","GU","C","F")
rownames(funcelasP)<-sn

#decrease by 10%
funcelasN <- (eig-funcsensN)/eig/(0.1)#10*(eig-eigsensN)/eig
colnames(funcelasN)<-c("S","I","U","GR","ES","GS","GG","GI","GU","C","F")
rownames(funcelasN)<-sn
funcname<-c("S","I","U","GR","ES","GS","GG","GI","GU","C","F")

#Use this if you want to plot function elasticities
funcelasP <- read.csv("C:/Users/Megan/IPMcoral/funcelasP.csv")
funcelasN <- read.csv("C:/Users/Megan/IPMcoral/funcelasN.csv")
Pnames<-funcelasP[,1]
Nnames<-funcelasN[,1]
funcelasP<-funcelasP[2:12]
rownames(funcelasP)<-Pnames
funcelasN<-funcelasN[2:12]
rownames(funcelasN)<-Nnames
funcelasP<-as.matrix(funcelasP)
barplot(funcelasP[1,],las=2,xlim=c(min(funcelasP[1,]),max(funcelasP[1,])),horiz=TRUE,col=cl1,axisnames=FALSE)
subparam<-c("Sexual","Clonal","FusionGrowth","FissionGrowth","NGS","ES","P(ES)","P(NGS)","P(Fusion)","P(Fission)","P(Survival)")
legend("topright", legend = subparam, col = cl1[11:1], pch=15, xjust = 1, bg = "white",bty="n",xpd=NA,cex=1)
mtext(expression(paste("Elasticity of ",lambda, " to parameter estimates")), side = 1, cex = 1,line=4, col = "grey20")

#Example of how function changes when multiplied by 1.1 or 0.9
plot(y,s.x(y,params)*1.1,type="l",col="red",xlab="Size at t",ylab="probability of survival")
lines(y,s.x(y,params))
lines(y,s.x(y,params)*0.9,col="blue")
legend("bottomright",legend=c("s.x","s.x+10%","s.x-10%"),col=c("black","red","blue"),lwd=1)

######################################################################################

