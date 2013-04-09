###########################################################MeganR IPM Code #############################################################################
library(nlme)
library(lme4)
library(VGAM)
#Range of values for lines added to plots
x<-seq(from=-2,to=10,by=0.01)
x0<-data.frame(size=x,size2=x*x)
nx<-9 # number of sites for this species
ncoef<-32 # number of coefficients
lmc<-lmeControl(niter=10000,msMaxIter=10000)
MCg<-within(MCg,fateIPM<-relevel(fateIPM,ref="growth"))
################################################################Vital Rate Functions###################################################################################

#survival regression
surv.reg=glmer(Psurv~size+utrans+(1|uyear),data=MC,family=binomial)
summary(surv.reg)
surv.int=c(attr(surv.reg,"fixef")[1], attr(surv.reg,"fixef")[1]+attr(surv.reg,"fixef")[3:(nx+1)])
surv.slope=attr(surv.reg,"fixef")[2]

# growth regression
growth.reg<-lme(sizeNext~size*utrans+size*fateIPM,random=~1|uquad,weights=varExp(form=~size|utrans),data=MCg,control=lmc)
summary(growth.reg)
growth.int=c(fixef(growth.reg)[1], fixef(growth.reg)[1]+fixef(growth.reg)[3:(nx+1)])
growth.slope=c(fixef(growth.reg)[2], fixef(growth.reg)[2]+fixef(growth.reg)[14:21])
growth.sd2=((summary(growth.reg)$sigma))^2
growth.varFunc=coef(growth.reg$modelStruct$varStruct,"exp")

fiss.int=c(fixef(growth.reg)[1]+fixef(growth.reg)[2+nx], fixef(growth.reg)[1]+fixef(growth.reg)[2+nx]+fixef(growth.reg)[3:(nx+1)])
fiss.slope=c(fixef(growth.reg)[2]+fixef(growth.reg)[22], fixef(growth.reg)[2]+fixef(growth.reg)[22]+fixef(growth.reg)[14:21])

fuse.int=c(fixef(growth.reg)[1]+fixef(growth.reg)[2+nx+1], fixef(growth.reg)[1]+fixef(growth.reg)[2+nx+1]+fixef(growth.reg)[3:(nx+1)])
fuse.slope=c(fixef(growth.reg)[2]+fixef(growth.reg)[23], fixef(growth.reg)[2]+fixef(growth.reg)[23]+fixef(growth.reg)[14:21])

M.int=c(fixef(growth.reg)[1]+fixef(growth.reg)[2+nx+2], fixef(growth.reg)[1]+fixef(growth.reg)[2+nx+2]+fixef(growth.reg)[3:(nx+1)])
M.slope=c(fixef(growth.reg)[2]+fixef(growth.reg)[24], fixef(growth.reg)[2]+fixef(growth.reg)[24]+fixef(growth.reg)[14:21])

#probability of growth (for fission, fusion, growth with diff growth functions)
grow=glm(Pgrowth~size*utrans,data=MC,family=binomial)
summary(grow)
Pgrow.int=c(coefficients(grow)[1],coefficients(grow)[1]+coefficients(grow)[3:(nx+1)])
Pgrow.slope=c(coefficients(grow)[2],coefficients(grow)[2]+coefficients(grow)[11:18])

#probability of cloning (for fission growth function and clonal matrix)
clone=glmer(Pfiss~size+utrans +(1|uyear),data=MC,family=binomial)
summary(clone)
Pfiss.int=c(fixef(clone)[1],fixef(clone)[1]+fixef(clone)[3:(nx+1)])
Pfiss.slope=fixef(clone)[2]

# probability of fusion (for fusion growth function)
fuse=glmer(Pfuse~size*utrans +(1|uyear),data=MC,family=binomial)
summary(fuse)
Pfuse.int=c(fixef(fuse)[1],fixef(fuse)[1]+fixef(fuse)[3:(nx+1)])
Pfuse.slope=c(fixef(fuse)[2],fixef(fuse)[2]+fixef(fuse)[11:18])

fuseLG=glm(PfuseLG~size+utrans,data=MC,family=binomial)
summary(fuseLG)
PfuseLG.int=c(coefficients(fuseLG)[1],coefficients(fuseLG)[1]+coefficients(fuseLG)[3:(nx+1)])
PfuseLG.slope=coefficients(fuseLG)[2]

#probability of M (fission and fusion in same time step.  For M growth function)
M<-glm(PM~size+utrans,data=MC,family=binomial)
summary(M)
PM.int=c(coefficients(M)[1],coefficients(M)[1]+coefficients(M)[3:(nx+1)])
PM.slope=coefficients(M)[2]

#number of clones per fission event 
clonenum=vglm(numfiss~utrans*size+fateIPM,data=MCf,family=pospoisson)
summary(clonenum)
Fclonenum.int=coefficients(clonenum)[1]
Mclonenum.int=coefficients(clonenum)[1]+coefficients(clonenum)[11]
clonenum.slope=c(coefficients(clonenum)[10],coefficients(clonenum)[10]+coefficients(clonenum)[12:19])

# size distribution of clones
clonesize=gls(sizeNext~size,data=MCfs,weights=varExp(form=~size|utrans))
summary(clonesize)
clonesize.int=coefficients(clonesize)[1]
clonesize.slope=coefficients(clonesize)[2]
clonesize.sd2=(summary(clonesize)$sigma)^2
clonesize.varFunc=coef(clonesize$modelStruct$varStruct,"exp")

#recruitment rate
#replace recruit.int with number of recruits per "adult" averaged quad over years then total of quads
recruit.int=0.0905897
#recruit.int=0.01

#size distribution of recruits
recruit.size.mean=mean(MC$sizeNext[MC$fateIPM=="recruit"])
recruit.size.sd=sd(MC$sizeNext[MC$fateIPM=="recruit"])


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
params[11,]<-c(M.int)                  #M growth intercept by site
params[12,]<-c(M.slope)                #M growth slope by site
params[13,]<-c(Pgrow.int)              #Probability of growth (not fission or fusion)intercept by site
params[14,]<-c(Pgrow.slope)          #probability of growth (not fission or fusion) slope, same across sites
params[15,]<-c(Pfiss.int)              #Probability of fission intercept by site
params[16,]<-rep(Pfiss.slope)          #Probability of fission slope, same across sites
params[17,]<-c(Pfuse.int)              #Probability of fusion intercept by site
params[18,]<-c(Pfuse.slope)            #Probability of fusion slope by site
params[19,]<-c(PfuseLG.int)            #Probability of being large fusion intercept given that you survived and fused by site
params[20,]<-rep(PfuseLG.slope)        #probability of being a large fusion slope given that you survived and fused same across sites
params[21,]<-c(PM.int)                 #Probability of M intercept by site
params[22,]<-rep(PM.slope)             #Probability of M slope, same across sites
params[23,]<-rep(clonesize.int)        #clone size intercept, same across sites
params[24,]<-rep(clonesize.slope)      #clone size slope, same across sites
params[25,]<-rep(clonesize.sd2)        #clone size sd, same across sites
params[26,]<-c(clonesize.varFunc)      #clone size variance increases with size by site
params[27,]<-rep(Fclonenum.int)        #number of clones from fissLG intercept, same across sites
params[28,]<-rep(Mclonenum.int)        #number of clones from M intercept, same across sites
params[29,]<-c(clonenum.slope)         #number of clones slope same for fates and across sites
params[30,]<-rep(recruit.int)          #number of recruits (kind of made up), same across sites     
params[31,]<-rep(recruit.size.mean)    #mean size of recruits, same across sites
params[32,]<-rep(recruit.size.sd)      #sd size of recruits, same across sites


## utility functions
for (site in 1:9){
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
  mg.yx=function(xp,x,params) {
    sigmaxp2<-params[5,site]*exp(2*params[6,site]*x)
    sigmaxp<-sqrt(sigmaxp2)
    g<-dnorm(xp,params[11,site]+params[12,site]*x,sigmaxp)
  }
#Probability of growing (not fission or fusing)
  p.grow.x=function(x,params) {
    u=exp(params[13,site]+params[14,site]*x)
    return(u/(1+u))
  }
# probability of fissing     
  p.fiss.x=function(x,params) {
    u=exp(params[15,site]+params[16,site]*x)
    return(u/(1+u))
  }
# probability of fusion
  p.fuse.x=function(x,params) {
    u=exp(params[17,site]+params[18,site]*x)
    return(u/(1+u))
  }
#probability of large fusion  
  p.fuseLG.x=function(x,params) {
    u=exp(params[19,site]+params[20,site]*x)
    return(u/(1+u))
  }
# probability of "M"
  p.M.x=function(x,params) {
    u=exp(params[21,site]+params[22,site]*x)
    return(u/(1+u))
  }
#Fission offspring 
  fc.yx=function(xp,x,params) {
    sigmaxp2<-params[25,site]*exp(2*params[26,site]*x)
    sigmaxp<-sqrt(sigmaxp2)
    p.fiss.x(x,params)*
    exp(params[27,site]+params[29,site]*x)/(1+exp(params[27,site]+params[29,site]*x))*
    dnorm(xp,params[23,site]+params[24,site]*x,sigmaxp)
  }
#"M" offspring 
  mc.yx=function(xp,x,params) {
    sigmaxp2<-params[25,site]*exp(2*params[26,site]+x)
    sigmaxp<-sqrt(sigmaxp2)
    p.M.x(x,params)*
    exp(params[28,site]+params[29,site]*x)/(1+exp(params[28,site]+params[29,site]*x))*
    dnorm(xp,params[23,site]+params[24,site]*x,sigmaxp)
  }
# number of clones per adult for "constant correction" of clonal matrix (no size distribution of offspring)
  c.y=function(x,params) {
    u=p.clone.x(x,params)*
    exp(params$clonenum.int+params$clonenum.slope*x)
    return (u)
  }
#recruitment = sexual offspring 
  f.yx=function(xp,x,params) {
    u=dnorm(xp,params[31,site],params[32,site])*
    exp(params[30,site])
    return (u)
  }
# recruits per adult for "constant correction" of fecundity matrix (no size distribution of offspring)
  f.y=function(xp,params) {
    u=exp(params$recruit.int)
    return (u)
  }
##########################################################discretization/build component kernels##################################
# integration limits - these limits span the range of sizes observed in the data set, and then some.
  min.size<-min(MC$size,na.rm=T)
  max.size<-max(MC$size,na.rm=T)+2
  
# number of cells in the discretized kernel
  n=100 
  
#boundary points (the edges of the cells defining the kernel)
  b=min.size+c(0:n)*(max.size-min.size)/n 
  
# mesh points= midpoints of the cells (could also use cumulative)
  y=0.5*(b[1:n]+b[2:(n+1)])
  
# width of the cells
  h=y[2]-y[1] 
  
# the function outer() evaluates the kernel at all pairwise combinations of the two vectors y and y. 
# P is the growth/suvival kernel.  Pi, Pu, and Pg are growth kernels for 3 growth fates (fission, fusion, and growth)
   
  S=s.x(y,params=params) # survival 
  Gr=p.grow.x(y,params=params)
  I=p.fiss.x(y,params=params)
  U=p.fuse.x(y,params=params)*p.fuseLG.x(y,params=params)
  M=p.M.x(y,params=params)
  
  GG=h*outer(y,y,gg.yx,params=params)
  GI=h*outer(y,y,ig.yx,params=params)
  GU=h*outer(y,y,ug.yx,params=params)
  GM=h*outer(y,y,mg.yx,params=params)
  
  Pg=GG
  for (i in 1:n) Pg[,i]=GG[,i]*S[i]*Gr[i]
  
  Pi=GI
  for(i in 1:n) Pi[,i]=GI[,i]*S[i]*I[i]
  
  Pu=GU
  for(i in 1:n) Pu[,i]=GU[,i]*S[i]*U[i]
  
  Pm=GM
  for(i in 1:n) Pu[,i]=GM[,i]*S[i]*M[i]
  
  
  
#C is the clonal offspring kernel
  fC=h*outer(y,y,fc.yx,params=params)
  mC=h*outer(y,y,mc.yx,params=params)
  
#F is the fecundity/sexually reproduced offspring kernel
  F=h*outer(y,y,f.yx,params=params)
  
  ###############################################################IPM!###############################################################
#Full kernel with 3 growth functions.  One each for fission, fusion and growth
  K=Pg+Pi+Pu+Pm+F+fC+mC
  sn<-levels(MC$utrans)[site]
  rn<-c("Pg","Pi","Pu","Pm","F","fC","mC","K")
  
  filled.contour(y,y,t(K),color=heat.colors,nlevels=20,cex.lab=1.5,main=sn, xlab= "Size at t",ylab="Size at t+1")
  
  retmat<-paste(rn[1:8],sn,sep="")
  assign(retmat[1],Pg)
  assign(retmat[2],Pi)
  assign(retmat[3],Pu)
  assign(retmat[4],Pm)
  assign(retmat[5],F)
  assign(retmat[6],fC)
  assign(retmat[7],mC)
  assign(retmat[8],K)
}

#get lamda,v,w  
lam=as.real(eigen(Kiu)$values[1]) 
w.eigen=as.real(eigen(K)$vectors[,1])
stable.dist=w.eigen/sum(w.eigen) 
v.eigen=as.real(eigen(t(K))$vectors[,1])
repro.val=v.eigen/v.eigen[1] 

#compute elasticity and sensitivity matrices
v.dot.w=sum(stable.dist*repro.val)*h
sens=outer(repro.val,stable.dist)/v.dot.w
elas=matrix(as.vector(sens)*as.vector(K)/lam,nrow=n)

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
P<- t((t(P)/nvals) * S)
nvals <- colSums(Pi, na.rm = TRUE)
Pi<- t((t(Pi)/nvals) * S)
nvals <- colSums(Pu, na.rm = TRUE)
Pu<- t((t(Pu)/nvals) * S)
nvals <- colSums(P, na.rm = TRUE)
Pg<- t((t(Pg)/nvals) * S)
#for F matrix
f<-vector(length=100)
for(i in 1:100) f[i]<- f.y(y,params=params)
correction.here <- f/colSums(F)
F<- t(t(F) * correction.here)
#for Cmatrix
c<-c.y(y,params=params)
correction.here <- c/colSums(C)
C<- t(t(C) * correction.here)
###########################################IPM with Constant Correction#############################################

Kconstant=P+F+C
lam2=as.real(eigen(Kconstant)$values[1])
Kiuconstant=Pg+Pi+Pu+F+C
lam3=as.real(eigen(Kiuconstant)$values[1])

##################some plots for IPM transition probabilities, sensitivities, and elasticities####################
par(mfrow=c(1,1))
filled.contour(y,y,t(K),color=heat.colors,nlevels=20,cex.lab=1.5,main="Kernel", xlab= "Size at t",ylab="Size at t+1")
filled.contour(y,y,t(Kiu),color=heat.colors,nlevels=20,cex.lab=1.5,main="Kernel", xlab= "Size at t",ylab="Size at t+1")
filled.contour(y,y,t(elas),color=heat.colors,nlevels=20,cex.lab=1.5,main="Kernel", xlab= "Size at t",ylab="Size at t+1")
filled.contour(y,y,t(sens),color=heat.colors,nlevels=20,cex.lab=1.5,main="Kernel", xlab= "Size at t",ylab="Size at t+1")
persp(y,y,t(K),phi=30,theta=300,col="lightgrey",border=NA,shade=0.75,main="Pmatrix:survival and growth", xlab= "Size at t",ylab="Size at t+1")
 


plot((MCg$asize)^0.5,(MCg$asizeNext)^0.5,col=c("red","green","blue")[MCg$fateIPM],main="Growth functions by fate",xlab="size",ylab="sizeNext")
legend("bottomright", legend = sprintf(c("growth","fission","fusion")), col = c("blue","red","green"), pch = 1, xjust = 1, bg = "white")
Ilm<-lm(sizeNext~size,data=MCfiss)
Ulm<-lm(sizeNext~size,data=MCfuse)
Glm<-lm(sizeNext~size,data=MCgrowth)
y0<- lines(x,predict(Ilm,newdata=x0,type='response'), col="red")
y0<- lines(x,predict(Ulm,newdata=x0,type='response'), col="green")
y0<- lines(x,predict(Glm,newdata=x0,type='response'), col="blue")


#plot the models over the data - figure 2
#par(mfrow=c(2,3),mar=c(4,4,2,1))
#plot(MCg$size,MCg$sizeNext,main='Growth/Shrinkage/Stasis')  
#lines(x,predict(growth.reg,newdata=x0),col='red',lwd=3)
#plot(MC$size,jitter(MC$surv),main='Survival') # jittered to see easier
#lines(x,predict(surv.reg,newdata=x0,type='response'), col='red',lwd=3)
#plot(MC$size,jitter(MC$fec1),main='Clones') # jittered to see easier
#lines(x,predict(clone,newdata=x0,type='response'), col='red',lwd=3)
#plot(MC$size,MC$fec2,main='Number of clones') # jittered to see easier
#lines(x,predict(clonenum,newdata=x0,type='response'), col='red',lwd=3)
#hist(MC$sizeNext[MC$fateIPM=="fissSM"],main='Size of Clones',freq=FALSE)
#lines(x,dnorm(x,mean(params$clonesize.int+params$clonesize.slope),params$clonesize.sd), col='red',lwd=3)
