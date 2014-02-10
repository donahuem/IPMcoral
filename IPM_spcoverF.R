###Alternate Fmatrix using spcover per quad
###########################################################MeganR IPM Code #############################################################################
library(nlme) # for gls()
library(lme4) # for recruitment
library(VGAM)#For vglmm function with zero truncated poisson
library(aqfig)# for vertical.image.legend
library(popbio) # for eigen.analysis()
library(MASS)# for negative binomial distribution
library(countreg) # for zerotrunc for geometric distribution
library(nnet) # for multinomial distribution
MCr$lnspcover<-log(1+MCr$spcover_t)
###Need this for Hi vs Lo clonenum kernel
MCf$hilo<-ifelse(MCf$numfiss<4,"lo","hi")
MCf$hilo<-as.factor(MCf$hilo)
hinum<-subset(MCf,MCf$numfiss>3)
lonum<-subset(MCf,MCf$numfiss<4)
MCf$Philo<-ifelse(MCf$hilo=="lo",1,0)

MCg<-within(MCg,fate<-relevel(fate,ref="growth")) #makes fate "growth" the zero-level in regression
dat<-MC #data for this species
growdat<-MCg #growth data for this species
fissdat<-MCf #fission only data for this species
fisssizedat<-MCfs #small fission only data for this species
recdat<-MCr # number of recruits per quad per year
recsizedat<-MCrs # recruits only for size distribution
densdat<-MCd #number of colonies per quad per year
denssizedat<-MCds # colonies present at time t (all but recruits and fissSM) for sfd
nx<-10 # number of sites for this species
nf<-3 #number of growth fates for this species
nfp<-6
ncoef<-30 # number of coefficients
n<-300 # number of cuts for meshpoints/discretization into matrix
################################################################Vital Rate Functions###################################################################################

#survival regression
surv.reg=glm(Psurv~size+utrans,data=dat,family=binomial)
summary(surv.reg)
surv.int=c(coefficients(surv.reg)[1], coefficients(surv.reg)[1]+coefficients(surv.reg)[3:(3+nx-1-1)])
surv.slope=coefficients(surv.reg)[2]


# growth regression
growth.reg<-gls(sizeNext~size*utrans+size*fate,weights=varExp(form=~size|utrans),data=growdat)
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
clone=glm(Pfiss~size*utrans,data=dat,family=binomial)
summary(clone)
Pfiss.int=c(coefficients(clone)[1],coefficients(clone)[1]+coefficients(clone)[3:(nx+1)])
Pfiss.slope=c(coefficients(clone)[2],coefficients(clone)[2]+coefficients(clone)[(nx+2):(2*nx)])

# probability of fusion conditional on survival(for fusion growth function)
fuse=glm(Pfuse~size+utrans,data=dat,family=binomial)
summary(fuse)
Pfuse.int=c(coefficients(fuse)[1],coefficients(fuse)[1]+coefficients(fuse)[3:(nx+1)])
Pfuse.slope=coefficients(fuse)[2]

fuseLG=glm(PfuseLG~size+utrans,data=dat,family=binomial)
summary(fuseLG)
PfuseLG.int=c(coefficients(fuseLG)[1],coefficients(fuseLG)[1]+coefficients(fuseLG)[3:(nx+1)])
PfuseLG.slope=coefficients(fuseLG)[2]

##Probability of hi or low fission product numbers conditional on fission
hilo <- glm(Philo~size,data=MCf,family=binomial)
summary(hilo)
Philo.int<-coefficients(hilo)[1]
Philo.slope<-coefficients(hilo)[2]
#number of clones per fission event: zero truncated poisson
#function increases too steeply over time.  Don't want the probability to go to zero or to balloon at higher values of size
clonenumhi=vglm(numfiss~size,data=hinum,family=pospoisson)
summary(clonenumhi)
clonenumhi.int=coefficients(clonenumhi)[1]
clonenumhi.slope=coefficients(clonenumhi)[2]
clonenumlo=vglm(numfiss~size,data=lonum,family=pospoisson)
summary(clonenumlo)
clonenumlo.int=coefficients(clonenumlo)[1]
clonenumlo.slope=coefficients(clonenumlo)[2]


# size distribution of clones
clonesize=gls(sizeNext~size,data=fisssizedat,weights=varExp(form=~size))
summary(clonesize)
clonesize.int=coefficients(clonesize)[1]
clonesize.slope=coefficients(clonesize)[2]
clonesize.sd2=(summary(clonesize)$sigma)^2
clonesize.varFunc=coef(clonesize$modelStruct$varStruct,"exp")

##recruitment rate as in Bruno et al. 2011 page 130
#number of recruits per quad.  Includes all years not just the years including years excluded from growth due to 2 year change etc.
rcrt.nb.SU0 <- glm.nb(recruits~lnspcover+utrans,data=MCr)
summary(rcrt.nb.SU0)
recruit.int=c(coefficients(rcrt.nb.SU0)[1],coefficients(rcrt.nb.SU0)[1]+coefficients(rcrt.nb.SU0)[3:(nx+1)])
recruit.slope=coefficients(rcrt.nb.SU0)[2]
##Recruitment rate using averages of observed recruits and colony density (not using lm just averaged by hand a while ago so data may have changed)
#recruit.int=0.0905897

##Tried different recruitment rates will think about more with sensitivity analyses
#recruit.int=0.01

#size distribution of recruits
rec.sizeU<-lm(sizeNext~utrans,data=recsizedat)
summary(rec.sizeU)
recruit.size.mean=c(coefficients(rec.sizeU)[1],coefficients(rec.sizeU)[1]+coefficients(rec.sizeU)[2:nx])
recruit.size.sd=(summary(rec.sizeU)$sigma)
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
params[12,]<-c(Pfiss.slope)          #Probability of fission slope, same across sites
params[13,]<-c(Pfuse.int)              #Probability of fusion intercept by site
params[14,]<-rep(Pfuse.slope)            #Probability of fusion slope by site
params[15,]<-c(PfuseLG.int)            #Probability of being large fusion intercept given that you survived and fused by site
params[16,]<-rep(PfuseLG.slope)        #probability of being a large fusion slope given that you survived and fused same across sites
params[17,]<-rep(clonesize.int)        #clone size intercept, same across sites
params[18,]<-rep(clonesize.slope)      #clone size slope, same across sites
params[19,]<-rep(clonesize.sd2)        #clone size sd, same across sites
params[20,]<-rep(clonesize.varFunc)    #clone size variance increases with size by site
params[21,]<-rep(clonenumlo.int)     #clonenum intercept lo number of fission events 
params[22,]<-rep(clonenumlo.slope)   #clonenum slope lo number of fission events
params[23,]<-rep(clonenumhi.int)     #clonenum intercept hi number of fission events
params[24,]<-rep(clonenumhi.slope)   #clonenum slope hi number of fission events
params[25,]<-rep(Philo.int)
params[26,]<-rep(Philo.slope)
params[27,]<-c(recruit.size.mean)      #mean size of recruits, same across sites
params[28,]<-rep(recruit.size.sd)      #sd size of recruits, same across sites
params[29,]<-c(recruit.int)
params[30,]<-rep(recruit.slope)


#rownames(params)<-c("surv.int","surv.slope","growth.int","growth.slope","growth.sd2","growth.varFunc","fiss.int","fiss.slope","fuse.int","fuse.slope","Pfiss.int","Pfiss.slope","Pfuse.int","Pfuse.slope","PfuseLG.int","PfuseLG.slope","clonesize.int","clonesize.slope","clonesize.sd2","clonesize.varFunc","clonenum.int","clonenum.slope","recruit.int","recruit.size.mean","recruit.size.sd","dens.int","clonemulti.int.1", "clonemulti.slope.1","clonemulti.int.2","clonemulti.slope.2","clonemulti.int.3","clonemult.slope.3","clonenb.int","clonenb.slope","clonegeo.int","clonegeo.slope")

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
#probability of lo fiss product number  
p.hilo.x=function(x,params) {
  u=exp(params[25,site]+params[26,site]*x)
  return(u/(1+u))
}
#Fission offspring using multinomial response
clo.yx=function(xp,x,params) {
  sigmaxp2<-params[19,site]*exp(2*params[20,site]*x)
  sigmaxp<-sqrt(sigmaxp2)
  nfiss<-exp(params[21,site]+params[22,site]*x)
  ##return: prob of fission * prob of lo number*number of fission products * prob that a fisser of size x prodcues a product of size xp
  p.fiss.x(x,params)*p.hilo.x(x,params)*nfiss*dnorm(xp,params[17,site]+params[18,site]*x,sigmaxp)
}
chi.yx=function(xp,x,params) {
  sigmaxp2<-params[19,site]*exp(2*params[20,site]*x)
  sigmaxp<-sqrt(sigmaxp2)
  nfiss<-exp(params[23,site]+params[24,site]*x)
  ##return: prob of fission * 1-prob of lo number*number of fission products * prob that a fisser of size x prodcues a product of size xp
  p.fiss.x(x,params)*(1-p.hilo.x(x,params))*nfiss*dnorm(xp,params[17,site]+params[18,site]*x,sigmaxp)
}
# number of clones per adult for "constant correction" of clonal matrix (no size distribution of offspring)
c.y=function(x,params) {
  nfiss<-exp(params[21,site]+params[22,site]*x)
  u=p.fiss.x(x,params)*nfiss
  return (u)
}

#recruitment = sexual offspring 
f.yx=function(xp,x,params){
  exp(params[29,site]+params[30,site]*x)*dnorm(xp,params[27,site],params[28,site])
}
# recruits per adult for "constant correction" of fecundity matrix (no size distribution of offspring)
f.y=function(xp,params) {
  r=exp(params[27,site])/exp(params[30,site])
  return (r)
}


##########################################################discretization/build component kernels##################################
# integration limits - these limits span the range of sizes observed in the data set, and then some.
min.size<-min(dat$size,na.rm=T)
max.size<-max(dat$size,na.rm=T)+2
#boundary points (the edges of the cells defining the kernel)
b=min.size+c(0:n)*(max.size-min.size)/n 

# mesh points= midpoints of the cells (could also use cumulative)
y=0.5*(b[1:n]+b[2:(n+1)])

# width of the cells
h=y[2]-y[1]

# the function outer() evaluates the kernel at all pairwise combinations of the two vectors y and y. 
# P is the growth/suvival kernel.  Pi, Pu, Pm,and Pg are growth kernels for 4 growth fates (fission, fusion, and growth)
Kernel<-function(y,n,params){
  
  
  #calculate the proportion of each size class that grows according to I,U,Gr functions
  S=s.x(y,params=params) # survival 
  I=p.fiss.x(y,params=params)-p.fiss.x(y,params=params)*p.fuse.x(y,params=params) #fisser-only = I
  U=(p.fuse.x(y,params=params)*p.fuseLG.x(y,params=params))  #fusion of U AND Ms:  include Ms in fusion growth
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
  
  #with hi vs lo models
  Clo=h*outer(y,y,clo.yx,params=params)
  Chi=h*outer(y,y,chi.yx,params=params)
  C=Clo+Chi
  
  #F is the fecundity/sexually reproduced offspring kernel
  #F is currently using the fecundity constant calculated based on Bruno et al. 2011
  F=h*outer(y,y,f.yx,params=params)
  ###############################################################IPM!###############################################################
  #Full kernel with 3 growth functions.  One each for fission, fusion (including M), and growth
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
  P<- t((t(P)/nvals) * S)*(1-p.fuse.x(y,params)*(1-p.fuseLG.x(y,params)))  # should be:  t((t(P)/nvals)*(S*p.smfuse)
  #want to normalize to the survival and the prob of being a small fuser
  
  #for F matrix
  #f<-vector(length=n)
  #for(i in 1:n) f[i]<- f.y(y,params=params)
  #correction.here <- f/colSums(F)
  #F<- t(t(F) * correction.here)
  
  #for Cmatrix
  #c<-c.y(y,params=params)
  #correction.here <- c/colSums(C)
  #C<- t(t(C) * correction.here)
  
  K=F
  
  return(K)
}  
