###########################################################MeganR IPM Code #############################################################################
library(nlme) # for gls()and gnls()
library(lme4) # for recruitment
library(aqfig)# for vertical.image.legend
library(popbio) # for eigen.analysis()

###Need this for Hi vs Lo clonenum kernel
MCf$hilo<-ifelse(MCf$numfiss<4,"lo","hi")
MCf$hilo<-as.factor(MCf$hilo)
hinum<-subset(MCf,MCf$numfiss>3)
lonum<-subset(MCf,MCf$numfiss<4)
MCf$Philo<-ifelse(MCf$hilo=="lo",1,0)


dat<-MC #data for this species
growdat<-MCg #growth data for this species
fissdat<-MCf #fission only data for this species
fisssizedat<-MCfs #small fission only data for this species
recdat<-MCr # number of recruits per quad per year
recsizedat<-MCrs # recruits only for size distribution
densdat<-MCd #number of colonies per quad per year
denssizedat<-MCds # colonies present at time t (all but recruits and fissSM) for sfd
growthonly<-subset(MCg,MCg$fateIPM=="growth")
fissonly<-subset(MCg,MCg$fateIPM=="fissLG")
fuseonly<-subset(MCg,MCg$fateIPM=="fuseLG"|MCg$fateIPM=="M")

nx<-10 # number of sites for this species
ncoef<-34 # number of coefficients
n<-300 # number of cuts for meshpoints/discretization into matrix

############ integration limits - these limits span the range of sizes observed in the data set, and then some.
min.size<-min(dat$size,na.rm=T)
max.size<-max(dat$size,na.rm=T)+2
#boundary points (the edges of the cells defining the kernel)
b=min.size+c(0:n)*(max.size-min.size)/n 

# mesh points= midpoints of the cells (could also use cumulative)
y=0.5*(b[1:n]+b[2:(n+1)])

# width of the cells
h=y[2]-y[1]
################################################################Vital Rate Functions###################################################################################

#survival regression
surv.reg=glm(Psurv~size+utrans,data=dat,family=binomial)
summary(surv.reg)
surv.int=c(coefficients(surv.reg)[1], coefficients(surv.reg)[1]+coefficients(surv.reg)[3:(3+nx-1-1)])
surv.slope=coefficients(surv.reg)[2]

# growth regression

growth.reg<-gnls(sizeNext~a+b*size+c*exp(-d*size),data=growthonly,params=list(a~utrans,b+c+d~1),start=list(a=c(1,1,1,1,1,1,1,1,1,1),b=1,c=1,d=1),weights=varExp(form=~size))
summary(growth.reg)
growth.int<-c(coefficients(growth.reg)[1],coefficients(growth.reg)[1]+coefficients(growth.reg)[2:10])
growth.slope<-coefficients(growth.reg)[11]
growth.e1<-coefficients(growth.reg)[12]
growth.e2<-coefficients(growth.reg)[13]
growth.sd2=((summary(growth.reg)$sigma))^2
growth.varFunc=coef(growth.reg$modelStruct$varStruct,"exp")


fiss.full<-lm(sizeNext~size,data=fissonly)
fiss.int<-coefficients(fiss.full)[1]
fiss.slope<-coefficients(fiss.full)[2]
fiss.sd<-summary(fiss.full)$sigma

fuse.full.var1<-gls(sizeNext~size,data=fuseonly,weights=varExp(form=~size))
summary(fuse.full.var1)
fuse.int<-coefficients(fuse.full.var1)[1]
fuse.slope<-coefficients(fuse.full.var1)[2]
fuse.varFunc<-coef(fuse.full.var1$modelStruct$varStruct,"exp")
fuse.sd2<-((summary(fuse.full.var1)$sigma))^2


library(VGAM)#For vglmm function with zero truncated poisson

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

#number of clones per fission event: 2 functions for data subset with hi (4-10) and lo (1-3) numfiss using zero truncated poisson

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
rcrt.lmeU.Y <- glmer(recruits~utrans+(1|year),family=poisson,data=recdat)
summary(rcrt.lmeU.Y)
recruit.int= c(fixef(rcrt.lmeU.Y)[1],fixef(rcrt.lmeU.Y)[1]+fixef(rcrt.lmeU.Y)[2:nx])


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

#number of colonies per quadIncludes all years not just the years including years excluded from growth due to 2 year change etc.
dens.lmeU.Y <- glmer(colnum~utrans+(1|year),family=poisson,data=densdat)
summary(dens.lmeU.Y)
dens.int=c(fixef(dens.lmeU.Y)[1],fixef(dens.lmeU.Y)[1]+fixef(dens.lmeU.Y)[2:nx])

# Size distribution of colonies colonies present at time t (all but recruits and fissSM)
col.size<-lm(sizeNext~utrans,data=denssizedat)
summary(col.size)
colsize.int=c(coefficients(col.size)[1],coefficients(col.size)[1]+coefficients(col.size)[2:10])
colsize.sd=(summary(col.size)$sigma)

r<-c(rep(0,nx))
for(i in 1:nx){
  recruitdensity=exp(recruit.int[i])
  colonydensity=exp(dens.int[i])
  sfd= dnorm(y,colsize.int[i],colsize.sd)
  r[i]<-recruitdensity/sum(sfd*y*colonydensity)
}

#make slots for parameters
#params<-data.frame(surv.int=NA,surv.slope=NA,growth.int=NA,growth.slope=NA,growth.sd=NA,growth.varFunc=NA,recruit.int=NA,recruit.size.mean=NA,recruit.size.sd=NA,clone.int=NA,clone.slope=NA,clonenum.int=NA,clonenum.slope=NA,clonesize.int=NA,clonesize.slope=NA,clonesize.sd=NA,fuse.int=NA,fuse.slope=NA,grow.int=NA,grow.slope=NA)

params<-matrix(0, ncoef, nx)  #creates a matrix of zeros with ncoef rows and nx (num of sites) columns

params[1,]<-c(surv.int)                #survival intercept by site
params[2,]<-rep(surv.slope)            #survival slope wrt size, same across sites
params[3,]<-c(growth.int)              #growth intercept by site
params[4,]<-rep(growth.slope)            #growth slope by site
params[5,]<-rep(growth.e1)
params[6,]<-rep(growth.e2)
params[7,]<-rep(growth.sd2)
params[8,]<-rep(growth.varFunc)            #variance in growth, same across sites
params[9,]<-rep(fiss.int)                #fiss growth intercept by site
params[10,]<-rep(fiss.slope)              #fiss growth slope by site
params[11,]<-rep(fiss.sd)
params[12,]<-rep(fuse.int)                #fuse growth intercept by site
params[13,]<-rep(fuse.slope)             #fuse growth slope by site
params[14,]<-rep(fuse.sd2)
params[15,]<-rep(fuse.varFunc)
params[16,]<-c(Pfiss.int)              #Probability of fission intercept by site
params[17,]<-c(Pfiss.slope)          #Probability of fission slope, same across sites
params[18,]<-c(Pfuse.int)              #Probability of fusion intercept by site
params[19,]<-rep(Pfuse.slope)            #Probability of fusion slope by site
params[20,]<-c(PfuseLG.int)            #Probability of being large fusion intercept given that you survived and fused by site
params[21,]<-rep(PfuseLG.slope)        #probability of being a large fusion slope given that you survived and fused same across sites
params[22,]<-rep(clonesize.int)        #clone size intercept, same across sites
params[23,]<-rep(clonesize.slope)      #clone size slope, same across sites
params[24,]<-rep(clonesize.sd2)        #clone size sd, same across sites
params[25,]<-rep(clonesize.varFunc)    #clone size variance increases with size by site
params[26,]<-rep(clonenumlo.int)     #clonenum intercept lo number of fission events 
params[27,]<-rep(clonenumlo.slope)   #clonenum slope lo number of fission events
params[28,]<-rep(clonenumhi.int)     #clonenum intercept hi number of fission events
params[29,]<-rep(clonenumhi.slope)   #clonenum slope hi number of fission events
params[30,]<-rep(Philo.int)
params[31,]<-rep(Philo.slope)
params[32,]<-c(recruit.size.mean)      #mean size of recruits, same across sites
params[33,]<-rep(recruit.size.sd)      #sd size of recruits, same across sites
params[34,]<-c(r)
rownames(params)<-c("surv.int","surv.slope","growth.int","growth.slope","growth.e1","growth.e2","growth.varFunc","growth.sd2","fiss.int","fiss.slope","fiss.sd","fuse.int","fuse.slope","fuse.sd2","fuse.varFunc","Pfiss.int","Pfiss.slope","Pfuse.int","Pfuse.slope","PfuseLG.int","PfuseLG.slope","clonesize.int","clonesize.slope","clonesize.sd2","clonesize.varFunc","clonenumlo.int","clonenumlo.slope","clonenumhi.int","clonenumhi.slope","Philo.int","Philo.slope","recruit.size.mean","recruit.size.sd","r")
## utility functions
#probability of surviving
s.x=function(x,params) {
  u=exp(params[1,site]+params[2,site]*x)
  return(u/(1+u))
}
#growth functions by fateIPM
gg.yx=function(xp,x,params) {
  sigmaxp2<-params[7,site]*exp(2*params[8,site]*x)
  sigmaxp<-sqrt(sigmaxp2)
  g<-dnorm(xp,params[3,site]+params[4,site]*x+params[5,site]*exp((-1)*params[6,site]*x),sigmaxp)
}
ig.yx=function(xp,x,params) {
  g<-dnorm(xp,params[9,site]+params[10,site]*x,params[11,site])
}
ug.yx=function(xp,x,params) {
  sigmaxp2<-params[14,site]*exp(2*params[15,site]*x)
  sigmaxp<-sqrt(sigmaxp2) 
  g<-dnorm(xp,params[12,site]+params[13,site]*x,sigmaxp)
}


# probability of fissing     
p.fiss.x=function(x,params) {
  u=exp(params[16,site]+params[17,site]*x)
  return(u/(1+u))
}
# probability of fusion
p.fuse.x=function(x,params) {
  u=exp(params[18,site]+params[19,site]*x)
  return(u/(1+u))
}
#probability of large fusion  
p.fuseLG.x=function(x,params) {
  u=exp(params[20,site]+params[21,site]*x)
  return(u/(1+u))
}

#probability of lo fiss product number  
p.hilo.x=function(x,params) {
  u=exp(params[30,site]+params[31,site]*x)
  return(u/(1+u))
}
#Fission offspring using separate vglm pospoisson function for hi vs lo number of fission products
clo.yx=function(xp,x,params) {
  sigmaxp2<-params[24,site]*exp(2*params[25,site]*x)
  sigmaxp<-sqrt(sigmaxp2)
  nfiss<-exp(params[26,site]+params[27,site]*x)
  ##return: prob of fission * prob of lo number*number of fission products * prob that a fisser of size x prodcues a product of size xp
  p.fiss.x(x,params)*p.hilo.x(x,params)*nfiss*dnorm(xp,params[22,site]+params[23,site]*x,sigmaxp)
}
chi.yx=function(xp,x,params) {
  sigmaxp2<-params[24,site]*exp(2*params[25,site]*x)
  sigmaxp<-sqrt(sigmaxp2)
  nfiss<-exp(params[28,site]+params[29,site]*x)
  ##return: prob of fission * 1-prob of lo number*number of fission products * prob that a fisser of size x prodcues a product of size xp
  p.fiss.x(x,params)*(1-p.hilo.x(x,params))*nfiss*dnorm(xp,params[22,site]+params[23,site]*x,sigmaxp)
}
# number of clones per adult for "constant correction" of clonal matrix (no size distribution of offspring)
clo.y=function(x,params) {
  nfiss<-exp(params[26,site]+params[27,site]*x)
  u=p.fiss.x(x,params)*p.hilo.x(x,params)*nfiss
  return (u)
}
chi.y=function(x,params) {
  nfiss<-exp(params[28,site]+params[29,site]*x)
  u= p.fiss.x(x,params)*(1-p.hilo.x(x,params))*nfiss
  return (u)
}
#recruitment = sexual offspring 

f.yx=function(xp,x,params){
  params[34,site]*x*dnorm(xp,params[32,site],params[33,site])
}

# recruits per adult for "constant correction" of fecundity matrix (no size distribution of offspring)
f.y=function(xp,params) {
  params[34,site]*xp
  
}


##########################################################discretization/build component kernels##################################


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
  f<- f.y(y,params=params)
  correction.here <- f/colSums(F)
  F<- t(t(F) * correction.here)
  
  #for Cmatrix
  c<-clo.y(y,params=params)+chi.y(y,params)
  correction.here <- c/colSums(C)
  C<- t(t(C) * correction.here)
  
  K=P+C+F
  return(K)
}  
