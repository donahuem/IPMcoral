###########################################################MeganR IPM Code #############################################################################
library(nlme)
library(lme4)
library(VGAM)#For vglmm function with zero truncated poisson
library(aqfig)# for vertical.image.legend
library(popbio)
MCg<-within(MCg,fate<-relevel(fate,ref="growth")) #makes fate "growth" the zero-level in regression
dat<-MC #data for this species
growdat<-MCg #growth data for this species
fissdat<-MCf #fission only data for this species
SMfissdat<-MCfs #small fission only data for this species
x<-seq(from=-2,to=10,by=0.01)#Range of values for lines added to plots
nx<-10 # number of sites for this species
nf<-3 #number of growth fates for this species
ncoef<-22 # number of coefficients
n<-300 # number of cuts for meshpoints/discretization into matrix
lmc<-lmeControl(niter=10000,msMaxIter=10000)
################################################################Vital Rate Functions###################################################################################

#survival regression
surv.reg=glm(Psurv~size+utrans,data=dat,family=binomial)
summary(surv.reg)
surv.int=c(coefficients(surv.reg)[1], coefficients(surv.reg)[1]+coefficients(surv.reg)[3:(3+nx-1-1)])
surv.slope=coefficients(surv.reg)[2]

# growth regression
growth.reg<-gls(sizeNext~size*utrans+size*fate,weights=varExp(form=~size|utrans),data=growdat,control=lmc)
summary(growth.reg)
growth.int=c(coefficients(growth.reg)[1], coefficients(growth.reg)[1]+coefficients(growth.reg)[3:(3+nx-1-1)])
growth.slope=c(coefficients(growth.reg)[2], coefficients(growth.reg)[2]+coefficients(growth.reg)[(nx+nf+1):(2*nx+nf-1)])
growth.sd2=((summary(growth.reg)$sigma))^2
growth.varFunc=coef(growth.reg$modelStruct$varStruct,"exp")

fiss.int=growth.int+coefficients(growth.reg)[2+nx]
fiss.slope=growth.slope+coefficients(growth.reg)[(2*nx+nf)]

fuse.int=growth.int+coefficients(growth.reg)[2+nx+1]
fuse.slope=growth.slope+coefficients(growth.reg)[(2*nx+nf+1)]

#probability of fission conditional on survival (for fission growth function and clonal matrix)
clone=glm(Pfiss~size+utrans,data=dat,family=binomial)
summary(clone)
Pfiss.int=c(coefficients(clone)[1],coefficients(clone)[1]+coefficients(clone)[3:(nx+1)])
Pfiss.slope=coefficients(clone)[2]

# probability of fusion conditional on survival(for fusion growth function)
fuse=glm(Pfuse~size*utrans,data=dat,family=binomial)
summary(fuse)
Pfuse.int=c(coefficients(fuse)[1],coefficients(fuse)[1]+coefficients(fuse)[3:(nx+1)])
Pfuse.slope=c(coefficients(fuse)[2],coefficients(fuse)[2]+coefficients(fuse)[(nx+2):(2*nx)])

fuseLG=glm(PfuseLG~size+utrans,data=dat,family=binomial)
summary(fuseLG)
PfuseLG.int=c(coefficients(fuseLG)[1],coefficients(fuseLG)[1]+coefficients(fuseLG)[3:(nx+1)])
PfuseLG.slope=coefficients(fuseLG)[2]

#number of clones per fission event: zero truncated poisson 
clonenum=vglm(numfiss~utrans*size,data=fissdat,family=pospoisson)
summary(clonenum)
clonenum.int=c(coefficients(clonenum)[1],coefficients(clonenum)[1]+coefficients(clonenum)[2:nx])
clonenum.slope=c(coefficients(clonenum)[nx+1],coefficients(clonenum)[nx+1]+coefficients(clonenum)[(nx+2):(2*nx)])

# size distribution of clones
clonesize=gls(sizeNext~size,data=SMfissdat,weights=varExp(form=~size))
summary(clonesize)
clonesize.int=coefficients(clonesize)[1]
clonesize.slope=coefficients(clonesize)[2]
clonesize.sd2=(summary(clonesize)$sigma)^2
clonesize.varFunc=coef(clonesize$modelStruct$varStruct,"exp")

#recruitment rate
#replace recruit.int with number of recruits per "adult" averaged quad over years then total of quads
#recruit.int=0.0905897
#recruit.int=0.01

#size distribution of recruits
#recruit.size.mean=mean(dat$sizeNext[dat$fateIPM=="recruit"])
#recruit.size.sd=sd(dat$sizeNext[dat$fateIPM=="recruit"])


#make slots for parameters
#params<-data.frame(surv.int=NA,surv.slope=NA,growth.int=NA,growth.slope=NA,growth.sd=NA,growth.varFunc=NA,recruit.int=NA,recruit.size.mean=NA,recruit.size.sd=NA,clone.int=NA,clone.slope=NA,clonenum.int=NA,clonenum.slope=NA,clonesize.int=NA,clonesize.slope=NA,clonesize.sd=NA,fuse.int=NA,fuse.slope=NA,grow.int=NA,grow.slope=NA)

params<-matrix(0, ncoef, nx)  #creates a matrix of zeros with ncoef rows and nx (num of sites) columns

params[1,]<-c(surv.int)                #survival intercept by site
params[2,]<-rep(surv.slope)            #survival slope wrt size, same across sites
params[3,]<-c(growth.int)              #growth intercept by site
params[4,]<-c(growth.slope)            #growth slope by site
params[5,]<-rep(growth.sd2)            #variance in growth, same across sites
params[6,]<-c(growth.varFunc)          #additional variance is exponential function of size|utrans
params[7,]<-c(fiss.int)                #fiss growth intercept by site
params[8,]<-c(fiss.slope)              #fiss growth slope by site
params[9,]<-c(fuse.int)                #fuse growth intercept by site
params[10,]<-c(fuse.slope)             #fuse growth slope by site
params[11,]<-c(Pfiss.int)              #Probability of fission intercept by site
params[12,]<-rep(Pfiss.slope)          #Probability of fission slope, same across sites
params[13,]<-c(Pfuse.int)              #Probability of fusion intercept by site
params[14,]<-c(Pfuse.slope)            #Probability of fusion slope by site
params[15,]<-c(PfuseLG.int)            #Probability of being large fusion intercept given that you survived and fused by site
params[16,]<-rep(PfuseLG.slope)        #probability of being a large fusion slope given that you survived and fused same across sites
params[17,]<-rep(clonesize.int)        #clone size intercept, same across sites
params[18,]<-rep(clonesize.slope)      #clone size slope, same across sites
params[19,]<-rep(clonesize.sd2)        #clone size sd, same across sites
params[20,]<-rep(clonesize.varFunc)    #clone size variance increases with size by site
params[21,]<-c(clonenum.int)           #number of clones intercept, same across sites
params[22,]<-c(clonenum.slope)         #number of clones slope same for fates and across sites
#params[23,]<-rep(recruit.int)          #number of recruits (kind of made up), same across sites     
#params[24,]<-rep(recruit.size.mean)    #mean size of recruits, same across sites
#params[25,]<-rep(recruit.size.sd)      #sd size of recruits, same across sites
rownames(params)<-c("surv.int","surv.slope","growth.int","growth.slope","growth.sd2","growth.varFunc","fiss.int","fiss.slope","fuse.int","fuse.slope","Pfiss.int","Pfiss.slope","Pfuse.int","Pfuse.slope","PfuseLG.int","PfuseLG.slope","clonesize.int","clonesize.slope","clonesize.sd2","clonesize.varFunc","clonenum.int","clonenum.slope")

## utility functions
#probability of surviving
  s.x=function(x,params) {
    u=exp(params[1,site]+params[2,site]*x)
    return(u/(1+u))
  }
#growth functions by fateIPM
  gg.yx=function(xp,x,params) {
    sigmaxp2<-params[5,site]*exp(2*params[6,site]*x)
    sigmaxp<-sqrt(sigmaxp2)
    g<-dnorm(xp,params[3,site]+params[4,site]*x,sigmaxp)
  }
  ig.yx=function(xp,x,params) {
    sigmaxp2<-params[5,site]*exp(2*params[6,site]*x)
    sigmaxp<-sqrt(sigmaxp2)
    g<-dnorm(xp,params[7,site]+params[8,site]*x,sigmaxp)
  }
  ug.yx=function(xp,x,params) {
    sigmaxp2<-params[5,site]*exp(2*params[6,site]*x)
    sigmaxp<-sqrt(sigmaxp2)
    g<-dnorm(xp,params[9,site]+params[10,site]*x,sigmaxp)
  }
  

# probability of fissing     
  p.fiss.x=function(x,params) {
    u=exp(params[11,site]+params[12,site]*x)
    return(u/(1+u))
  }
# probability of fusion
  p.fuse.x=function(x,params) {
    u=exp(params[13,site]+params[14,site]*x)
    return(u/(1+u))
  }
#probability of large fusion  
  p.fuseLG.x=function(x,params) {
    u=exp(params[15,site]+params[16,site]*x)
    return(u/(1+u))
  }

#Fission offspring 
  c.yx=function(xp,x,params) {
    sigmaxp2<-params[19,site]*exp(2*params[20,site]*x)
    sigmaxp<-sqrt(sigmaxp2)
    nfiss<-exp(params[21,site]+params[22,site]*x)
    #return: prob of fission * number of fission products * prob that a fisser of size x prodcues a product of size xp
    p.fiss.x(x,params)*nfiss*dnorm(xp,params[17,site]+params[18,site]*x,sigmaxp)
  }

# number of clones per adult for "constant correction" of clonal matrix (no size distribution of offspring)
  c.y=function(x,params) {
    nfiss<-exp(params[21,site]+params[22,site]*x)
    u=p.fiss.x(x,params)*nfiss
    return (u)
  }
#recruitment = sexual offspring 
  #f.yx=function(xp,x,params) {
   # u=dnorm(xp,params[23,site],params[24,site])*
    #exp(params[25,site])
    #return (u)
 # }
# recruits per adult for "constant correction" of fecundity matrix (no size distribution of offspring)
  #f.y=function(xp,params) {
   # u=exp(params$recruit.int)
    #return (u)
  #}

#####################################################################Plot Vital Rate Regressions#######################################################
##Survival##
cl<-c("cyan","mediumturquoise","dodgerblue","green3","purple","olivedrab2","orchid1","blue","darkblue")

#cl<-rainbow(nx)
plot(dat$size,dat$Psurv,col=cl[dat$utrans],type="p",yaxt="n",xlab="Size (t)",ylab="Survival Probability")
axis(2,at=seq(0,1, by=0.2))
for (site in 1:(nx)) {
  lines(x,s.x(x,params=params),col=cl[site],lwd=2)
}
legend("right", legend = sprintf(levels(growdat$utrans)), col = cl, lwd = 2, xjust = 1, bg = "white")

##Growth##
plot(growdat$size,growdat$sizeNext,col=c(cl[3],cl[8],cl[5],cl[1])[growdat$fateIPM],pch=18,cex=0.95,ylab="Size (t+1)",xlab="Size (t)")
fategrow<-lm(sizeNext~size+fate,data=growdat)
lines(x,coefficients(fategrow)[1]+coefficients(fategrow)[2]*x,col=cl[3],lwd=2)
lines(x,coefficients(fategrow)[1]+coefficients(fategrow)[3]+coefficients(fategrow)[2]*x,col=cl[8],lwd=2)
lines(x,coefficients(fategrow)[1]+coefficients(fategrow)[4]+coefficients(fategrow)[2]*x,col=cl[5],lwd=2)
legend("bottomright", legend = sprintf(levels(growdat$fate)), col=c(cl[3],cl[8],cl[5],cl[1]), lty=1,lwd=2, xjust = 1, bg = "white")

par(mfrow=c(2,2))

plot(growdat$size,growdat$sizeNext,col=cl[growdat$utrans],type="p",pch=18,cex=0.65,main="Stasis",xlab="Size (t)",ylab="Size (t+1)")
for (site in 1:(nx)) {
  lines(x,params[3,site]+params[4,site]*x,col=cl[site],lwd=2)
}

plot(growdat$size,growdat$sizeNext,col=cl[growdat$utrans],type="p",pch=18,cex=0.65,main="Fission",xlab="Size (t)",ylab="Size (t+1)")
for (site in 1:(nx)) {
  lines(x,params[7,site]+params[8,site]*x,col=cl[site],lwd=2)
}

plot(growdat$size,growdat$sizeNext,col=cl[growdat$utrans],type="p",pch=18,cex=0.65,main="Fusion",xlab="Size (t)",ylab="Size (t+1)")
for (site in 1:(nx)) {
  lines(x,params[9,site]+params[10,site]*x,col=cl[site],lwd=2)
}

legend("bottomright", legend = sprintf(levels(growdat$utrans)), col = cl, lwd=2, xjust = 1, bg = "white")

##Fate Probabilities##
par(mfrow=c(1,1))
#Probability of growth#
plot(dat$size,dat$Pgrowth,col=cl[dat$utrans],type="p",ylim=c(0,1),xlab="Size (t)",ylab="Prob(growth)")
for (site in 1:(nx)) {
  lines(x,1-p.fiss.x(x,params=params)-p.fiss.x(x,params=params)*p.fuse.x(x,params=params)*p.fuseLG.x(x,params=params)-p.fuse.x(x,params=params)*p.fuseLG.x(x,params=params)-p.fiss.x(x,params=params)*p.fuse.x(x,params=params)*p.fuseLG.x(x,params=params),col=cl[site],lwd=2)
}
legend("topright", legend = sprintf(levels(growdat$utrans)), col = cl, lwd = 2, xjust = 1, bg = "white")
##Probability of fission##
plot(dat$size,dat$Pfiss,col=cl[dat$utrans],type="p",ylim=c(0,1),xlab="Size (t)",ylab="Prob(fission)")
for (site in 1:(nx)) {
  lines(x,p.fiss.x(x,params=params)-p.fiss.x(x,params=params)*p.fuse.x(x,params=params)*p.fuseLG.x(x,params=params),col=cl[site],lwd=2)
}
legend("left", legend = sprintf(levels(growdat$utrans)), col = cl, lwd=2, xjust = 1, bg = "white")

#Probability of Large fusion#
plot(dat$size,dat$Pfuse,col=cl[dat$utrans],type="p",ylim=c(0,1),xlab="Size (t)",ylab="Prob(fusion)")
for (site in 1:(nx)) {
  lines(x,p.fuse.x(x,params=params)*p.fuseLG.x(x,params=params)-p.fiss.x(x,params=params)*p.fuse.x(x,params=params)*p.fuseLG.x(x,params=params),col=cl[site],lwd=2)
}
legend("left", legend = sprintf(levels(growdat$utrans)), col = cl, lwd=2, xjust = 1, bg = "white")



#Number of fission products#
par(mfrow=c(1,1))
plot(fissdat$size,fissdat$numfiss,col=cl[fissdat$utrans],type="p",pch=18,xlab="Size (t)",ylab="# of fission products/'parent colony'")
for (site in 1:(nx)) {
  lines(x,exp(params[21,site]+params[22,site]*x)+1,col=cl[site],lwd=2)
}
legend("topleft", legend = sprintf(levels(fissdat$utrans)), col = cl, lwd = 2, xjust = 1, bg = "white")

##Does same as above but predicts using vglmm instead of my function to get density.  See Megan D's note below about exp(G)+1##
temp_lap <- rep(x,nx)
temp_ut<-{}
for (i in 1:nx) {
  temp_ut <- c(temp_ut,rep(unique(fissdat$utrans)[i],length(x)))
}
temp_u <- unique(fissdat$utrans)[temp_ut]
MyData<-data.frame(size=temp_lap,utrans=temp_u)  
G <- predictvglm(clonenum,newdata=MyData,type="link",se=TRUE)#with "link" can get s.e. bands; need to exp(G) + 1 for predictions
plot(fissdat$size,fissdat$numfiss,col=cl[fissdat$utrans],type="p",xlab="Log(Area)[t]",ylim=c(0.5,10.5),xlim=c(-2,10),ylab="Number of Fission Products")
lines(MyData$size[MyData$utrans==unique(MyData$utrans)[1]],exp(G$fit[MyData$utrans==unique(MyData$utrans)[1]])+1,lty=1,col=cl[1])
#title("Number of Fission Products ~ Zero-Trunc-Poisson(Log(Area) + Transect)")
for (i in 2:nx) {
  lines(MyData$size[MyData$utrans==unique(MyData$utrans)[i]],exp(G$fit[MyData$utrans==unique(MyData$utrans)[i]])+1,lty=1,col=cl[i])
}
legend("topleft",legend=unique(fissdat$utrans),col=cl[1:nx],lty=1,pt.cex=0.3,bty="n")

##size distribution of clones##
hist(SMfissdat$sizeNext,ylim=c(0,1),freq=FALSE)
for (site in 1:nx){
  sigmaxp2<-params[19,site]*exp(2*params[20,site]*x)
  sigmaxp<-sqrt(sigmaxp2)
  lines(x,dnorm(x,(params[17,site]+params[18,site]*x),sigmaxp),col=cl[site])
}


##########################################################discretization/build component kernels##################################
# integration limits - these limits span the range of sizes observed in the data set, and then some.
min.size<-min(dat$size,na.rm=T)
max.size<-max(dat$size,na.rm=T)+2
# the function outer() evaluates the kernel at all pairwise combinations of the two vectors y and y. 
# P is the growth/suvival kernel.  Pi, Pu, Pm,and Pg are growth kernels for 4 growth fates (fission, fusion, and growth)
Kernel<-function(n,params,site){

  #boundary points (the edges of the cells defining the kernel)
  b=min.size+c(0:n)*(max.size-min.size)/n 
  
  # mesh points= midpoints of the cells (could also use cumulative)
  y=0.5*(b[1:n]+b[2:(n+1)])
  
  # width of the cells
  h=y[2]-y[1] 
  S=s.x(y,params=params) # survival 
  I=p.fiss.x(y,params=params)-p.fiss.x(y,params=params)*p.fuse.x(y,params=params) #fisser-only = I
  U=(p.fuse.x(y,params=params)*p.fuseLG.x(y,params=params)  #fusion of U AND Ms:  include Ms in fusion growth
  Gr=1-p.fiss.x(y,params=params) - p.fuse.x(y,params=params) + p.fiss.x(y,params=params)*p.fuse.x(y,params=params)
  
  GG=h*outer(y,y,gg.yx,params=params)
  GI=h*outer(y,y,ig.yx,params=params)
  GU=h*outer(y,y,ug.yx,params=params)
 
  
  Pg=GG
  for (i in 1:n) Pg[,i]=GG[,i]*S[i]*(Gr[i])
  
  Pi=GI
  for(i in 1:n) Pi[,i]=GI[,i]*S[i]*I[i]
  
  Pu=GU
  for(i in 1:n) Pu[,i]=GU[,i]*S[i]*U[i]
  
    
#C is the clonal offspring kernel
  C=h*outer(y,y,c.yx,params=params)
  
#F is the fecundity/sexually reproduced offspring kernel
  #F=h*outer(y,y,f.yx,params=params)
  
  ###############################################################IPM!###############################################################
#Full kernel with 3 growth functions.  One each for fission, fusion (including M) and growth
  P = Pg+Pi+Pu
  #####################################################################constant correction#######################################################
  #which will multiply every column of the IPM by a constant sufficient to adjust values to those predicted for survival at that size
  #For Pmatrix
  nvals <- colSums(P, na.rm = TRUE)
  loc0 <- which(nvals == 0, arr.ind = TRUE)
  if (length(loc0) > 0) {
    print("warnings - columns that sum to 0 or that have NAs - assuming survival is along the diagonal; plot your Pmatrix to check it")
    P[, loc0] <- 0
    P[cbind(loc0, loc0)] <- S
  }
  nvals <- colSums(P, na.rm = TRUE)
  P<- t((t(P)/nvals) * S*(I+U+Gr))  # should be:  t((t(P)/nvals)*(S*p.smfuse)
     #want to normalize to the survival and the prob of being a small fuser
  #for F matrix
  #f<-vector(length=n)
  #for(i in 1:n) f[i]<- f.y(y,params=params)
  #correction.here <- f/colSums(F)
  #F<- t(t(F) * correction.here)
  #for Cmatrix
  c<-c.y(y,params=params)
  correction.here <- c/colSums(C)
  C<- t(t(C) * correction.here)
    
  K=P+C
  
  return(K)
}  
###############################################################Lambda,SSD,RV,Damping ratio, Elasticity matrix###############################################################
par(mfrow = c(3, 3))
par(cex = 0.6)
par(mar = c(1,1,1,1), oma = c(4, 4, 2, 4))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))
n=300
#boundary points (the edges of the cells defining the kernel)
b=min.size+c(0:n)*(max.size-min.size)/n 

# mesh points= midpoints of the cells (could also use cumulative)
y=0.5*(b[1:n]+b[2:(n+1)])
for(site in 1:nx){
  K<-Kernel(n,params,site)
  sn<-levels(dat$utrans)[site]
  image(y,y,t(K), col=topo.colors(300),main=sn,axes=FALSE)
}
mtext("Size (t)", side = 1, outer = TRUE, cex = 1, line = 2.2,col = "grey20")
mtext("Size (t+1)", side = 2, outer = TRUE, cex = 1, line = 2.2,col = "grey20")
par(mfrow=c(1,1))
par(mar=c(4,1,1,1),oma=c(1,1,1,1))
vertical.image.legend(col=topo.colors(n),zlim=c(0,0.05))

## Lambda, SSD, RV, Damping ratio, Elasticity matrix
big<-array(0, c(n, n, nx))
eig<-vector("numeric", nx)

for (site in 1:nx){
  big[,,site]<-Kernel(n, params, site)
  A<-eigen.analysis(big[,,site])
  eig[site]<-A$lambda1
}

eigsensP<-array(0,c(nx,length(params[,1])))
eigsensN<-array(0,c(nx,length(params[,1])))
#sensN<-array(0, c(n, n, nx,length(params[,1])))
#sensP<-array(0, c(n, n, nx,length(params[,1])))

for (site in 1:nx){
  for (j in 1:length(params[,1])) {
    svecP<-params
    svecP[j,]<-1.1*svecP[j,]
    sensP<-Kernel(n, params=svecP, site)
    eigsensP[site,j]<-eigen.analysis(sensP)$lambda1
    svecN<-params
    svecN[j,]<-0.9*svecN[j,]
    sensN<-Kernel(n, params=svecN, site)
    eigsensN[site,j]<-eigen.analysis(sensN)$lambda1
  }
}

elasP <- (eig-eigsensP)/eig/(0.1)
par(mfrow=c(2,2))
par(mar = c(7, 2, 2, 2) + 0.1)
for (i in c(1,4,9,15)){
  barplot(elasP[,i],horiz=TRUE,xlab="Elasticity",col=cl[1:nx])
}
par(mfrow = c(2, 2))
par(mar = c(7, 2, 2, 2) + 0.1)
for (i in c(1:nx)){
barplot(elasP[i,],main="Elasticity",names.arg=rownames(params),las=2)
}

#legend("bottom",legend=sort(unique(utrans),decreasing=TRUE),horiz=TRUE,pch=1,col=cl[ns:1],cex=1,bty="n")

## To get lambda for each site
eig

par(mar = c(7, 4, 4, 2) + 0.1)
xvals<-barplot(eig)
text(xvals, par("usr")[3] - 0.25, srt = 45, adj = 1,labels = unique(dat$utrans), xpd = TRUE)


