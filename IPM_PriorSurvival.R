###Alternate survival including fate_prior
###########################################################MeganR IPM Code #############################################################################
library(nlme) # for gls()
library(lme4) # for recruitment
library(VGAM)#For vglmm function with zero truncated poisson
library(aqfig)# for vertical.image.legend
library(popbio) # for eigen.analysis()
library(MASS)# for negative binomial distribution
library(countreg) # for zerotrunc for geometric distribution
library(nnet) # for multinomial distribution

###Need this for Hi vs Lo clonenum kernel
MCf$hilo<-ifelse(MCf$numfiss<4,"lo","hi")
MCf$hilo<-as.factor(MCf$hilo)
hinum<-subset(MCf,MCf$numfiss>3)
lonum<-subset(MCf,MCf$numfiss<4)
MCf$Philo<-ifelse(MCf$hilo=="lo",1,0)

##need this for non-markovian survival
MC$fate_priorSM<-MC$fate_prior
MC$fate_priorSM<-ifelse(MC$allfates=="ID"|MC$allfates=="II"|MC$allfates=="IM"|MC$allfates=="IS"|MC$allfates=="IU","Is",as.character(MC$fate_priorSM))
MC$fate_priorSM<-as.factor(MC$fate_priorSM)
index<-c("I","Is","M","R","S","U")
values<-c(1,0,0,0,0,0)
MC$PI<-values[match(MC$fate_priorSM,index)] # To estimate probability of I fate_priorSM
values<-c(0,1,0,0,0,0)
MC$PIs<-values[match(MC$fate_priorSM,index)] # To estimate probability of I fate_priorSM
values<-c(0,0,1,0,0,0)
MC$PM<-values[match(MC$fate_priorSM,index)] # To estimate probability of I fate_priorSM
values<-c(0,0,0,1,0,0)
MC$PR<-values[match(MC$fate_priorSM,index)] # To estimate probability of I fate_priorSM
values<-c(0,0,0,0,1,0)
MC$PS<-values[match(MC$fate_priorSM,index)] # To estimate probability of I fate_priorSM
values<-c(0,0,0,0,0,1)
MC$PU<-values[match(MC$fate_priorSM,index)] # To estimate probability of I fate_priorSM

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
ncoef<-51 # number of coefficients
n<-300 # number of cuts for meshpoints/discretization into matrix
################################################################Vital Rate Functions###################################################################################

#survival regression
surv.reg=glm(Psurv~utrans+size*fate_priorSM,data=dat,family=binomial)
summary(surv.reg)
survI.int=c(coefficients(surv.reg)[1], coefficients(surv.reg)[1]+coefficients(surv.reg)[2:(nx)])
survI.slope=coefficients(surv.reg)[nx+1]
survIs.int=c(survI.int+coefficients(surv.reg)[nx+2])
survIs.slope=survI.slope+coefficients(surv.reg)[nx+nfp+1]
survM.int=c(survI.int+coefficients(surv.reg)[nx+3])
survM.slope=survI.slope+coefficients(surv.reg)[nx+nfp+2]
survR.int=c(survI.int+coefficients(surv.reg)[nx+4])
survR.slope=survI.slope+coefficients(surv.reg)[nx+nfp+3]
survS.int=c(survI.int+coefficients(surv.reg)[nx+5])
survS.slope=survI.slope+coefficients(surv.reg)[nx+nfp+4]
survU.int=c(survI.int+coefficients(surv.reg)[nx+6])
survU.slope=survI.slope+coefficients(surv.reg)[nx+nfp+5]

PIprior<-glm(PI~size,data=dat,family=binomial)
PIprior.int=coefficients(PIprior)[1]
PIprior.slope=coefficients(PIprior)[2]

PIsprior<-glm(PIs~size+utrans,data=dat,family=binomial)
PIsprior.int=c(coefficients(PIprior)[1],coefficients(PIsprior)[1]+coefficients(PIsprior)[3:(nx+2-1)])
PIsprior.slope=coefficients(PIsprior)[2]

PMprior<-glm(PM~size+utrans,data=dat,family=binomial)
PMprior.int=c(coefficients(PMprior)[1],coefficients(PMprior)[1]+coefficients(PMprior)[3:(nx+2-1)])
PMprior.slope=coefficients(PMprior)[2]

PRprior<-glm(PR~size+utrans,data=dat,family=binomial)
PRprior.int=c(coefficients(PRprior)[1],coefficients(PRprior)[1]+coefficients(PRprior)[3:(nx+2-1)])
PRprior.slope=coefficients(PRprior)[2]

PSprior<-glm(PS~size+utrans,data=dat,family=binomial)
PSprior.int=c(coefficients(PSprior)[1],coefficients(PSprior)[1]+coefficients(PSprior)[3:(nx+2-1)])
PSprior.slope=coefficients(PSprior)[2]

PUprior<-glm(PU~size*utrans,data=dat,family=binomial)
PUprior.int=c(coefficients(PUprior)[1],coefficients(PUprior)[1]+coefficients(PUprior)[3:(nx+2-1)])
PUprior.slope=c(coefficients(PUprior)[2],coefficients(PUprior)[2]+coefficients(PUprior)[(nx+2):(2*nx)])
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
clonenumhi.int=coefficients(clonenum)[1]
clonenumhi.slope=coefficients(clonenum)[2]
clonenumlo=vglm(numfiss~size,data=lonum,family=pospoisson)
summary(clonenumlo)
clonenumlo.int=coefficients(clonenum)[1]
clonenumlo.slope=coefficients(clonenum)[2]


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

params[1,]<-c(survI.int)                #survival intercept by site
params[2,]<-rep(survI.slope)            #survival slope wrt size, same across sites
params[3,]<-c(survIs.int)                #survival intercept by site
params[4,]<-rep(survIs.slope)            #survival slope wrt size, same across sites
params[5,]<-c(survM.int)                #survival intercept by site
params[6,]<-rep(survM.slope)            #survival slope wrt size, same across sites
params[7,]<-c(survR.int)                #survival intercept by site
params[8,]<-rep(survR.slope)            #survival slope wrt size, same across sites
params[9,]<-c(survS.int)                #survival intercept by site
params[10,]<-rep(survS.slope)            #survival slope wrt size, same across sites
params[11,]<-c(survU.int)                #survival intercept by site
params[12,]<-rep(survU.slope)            #survival slope wrt size, same across sites
params[13,]<-c(growth.int)              #growth intercept by site
params[14,]<-c(growth.slope)            #growth slope by site
params[15,]<-rep(growth.sd2)            #variance in growth, same across sites
params[16,]<-c(growth.varFunc)          #additional variance is exponential function of size|utrans
params[17,]<-c(fiss.int)                #fiss growth intercept by site
params[18,]<-c(fiss.slope)              #fiss growth slope by site
params[19,]<-c(fuse.int)                #fuse growth intercept by site
params[20,]<-c(fuse.slope)             #fuse growth slope by site
params[21,]<-c(Pfiss.int)              #Probability of fission intercept by site
params[22,]<-c(Pfiss.slope)          #Probability of fission slope, same across sites
params[23,]<-c(Pfuse.int)              #Probability of fusion intercept by site
params[24,]<-rep(Pfuse.slope)            #Probability of fusion slope by site
params[25,]<-c(PfuseLG.int)            #Probability of being large fusion intercept given that you survived and fused by site
params[26,]<-rep(PfuseLG.slope)        #probability of being a large fusion slope given that you survived and fused same across sites
params[27,]<-rep(clonesize.int)        #clone size intercept, same across sites
params[28,]<-rep(clonesize.slope)      #clone size slope, same across sites
params[29,]<-rep(clonesize.sd2)        #clone size sd, same across sites
params[30,]<-rep(clonesize.varFunc)    #clone size variance increases with size by site
params[31,]<-rep(clonenumlo.int)     #clonenum intercept lo number of fission events 
params[32,]<-rep(clonenumlo.slope)   #clonenum slope lo number of fission events
params[33,]<-rep(clonenumhi.int)     #clonenum intercept hi number of fission events
params[34,]<-rep(clonenumhi.slope)   #clonenum slope hi number of fission events
params[35,]<-rep(Philo.int)
params[36,]<-rep(Philo.slope)
params[37,]<-c(recruit.size.mean)      #mean size of recruits, same across sites
params[38,]<-rep(recruit.size.sd)      #sd size of recruits, same across sites
params[39,]<-c(r)           
params[40,]<-rep(PIprior.int)
params[41,]<-rep(PIprior.slope)
params[42,]<-c(PIsprior.int)
params[43,]<-rep(PIsprior.slope)
params[44,]<-c(PMprior.int)
params[45,]<-rep(PMprior.slope)
params[46,]<-c(PRprior.int)
params[47,]<-rep(PRprior.slope)
params[48,]<-c(PSprior.int)
params[49,]<-rep(PSprior.slope)
params[50,]<-c(PUprior.int)
params[51,]<-c(PUprior.slope)

#rownames(params)<-c("surv.int","surv.slope","growth.int","growth.slope","growth.sd2","growth.varFunc","fiss.int","fiss.slope","fuse.int","fuse.slope","Pfiss.int","Pfiss.slope","Pfuse.int","Pfuse.slope","PfuseLG.int","PfuseLG.slope","clonesize.int","clonesize.slope","clonesize.sd2","clonesize.varFunc","clonenum.int","clonenum.slope","recruit.int","recruit.size.mean","recruit.size.sd","dens.int","clonemulti.int.1", "clonemulti.slope.1","clonemulti.int.2","clonemulti.slope.2","clonemulti.int.3","clonemult.slope.3","clonenb.int","clonenb.slope","clonegeo.int","clonegeo.slope")

## utility functions
#probability of surviving
sI.x=function(x,params) {
  u=exp(params[1,site]+params[2,site]*x)
  return(u/(1+u))
}
sIs.x=function(x,params) {
  u=exp(params[3,site]+params[4,site]*x)
  return(u/(1+u))
}
sM.x=function(x,params) {
  u=exp(params[5,site]+params[6,site]*x)
  return(u/(1+u))
}
sR.x=function(x,params) {
  u=exp(params[7,site]+params[8,site]*x)
  return(u/(1+u))
}
sS.x=function(x,params) {
  u=exp(params[9,site]+params[10,site]*x)
  return(u/(1+u))
}
sU.x=function(x,params) {
  u=exp(params[11,site]+params[12,site]*x)
  return(u/(1+u))
}

p.I.x=function(x,params) {
  u=exp(params[40,site]+params[41,site]*x)
  return(u/(1+u))
}
p.Is.x=function(x,params) {
  u=exp(params[42,site]+params[43,site]*x)
  return(u/(1+u))
}
p.M.x=function(x,params) {
  u=exp(params[44,site]+params[45,site]*x)
  return(u/(1+u))
}
p.R.x=function(x,params) {
  u=exp(params[46,site]+params[47,site]*x)
  return(u/(1+u))
}
p.S.x=function(x,params) {
  u=exp(params[48,site]+params[49,site]*x)
  return(u/(1+u))
}
p.U.x=function(x,params) {
  u=exp(params[50,site]+params[51,site]*x)
  return(u/(1+u))
}

#growth functions by fateIPM
gg.yx=function(xp,x,params) {
  sigmaxp2<-params[15,site]*exp(2*params[16,site]*x)
  sigmaxp<-sqrt(sigmaxp2)
  g<-dnorm(xp,params[13,site]+params[14,site]*x,sigmaxp)
}
ig.yx=function(xp,x,params) {
  sigmaxp2<-params[15,site]*exp(2*params[16,site]*x)
  sigmaxp<-sqrt(sigmaxp2)
  g<-dnorm(xp,params[17,site]+params[18,site]*x,sigmaxp)
}
ug.yx=function(xp,x,params) {
  sigmaxp2<-params[15,site]*exp(2*params[16,site]*x)
  sigmaxp<-sqrt(sigmaxp2)
  g<-dnorm(xp,params[19,site]+params[20,site]*x,sigmaxp)
}


# probability of fissing     
p.fiss.x=function(x,params) {
  u=exp(params[21,site]+params[22,site]*x)
  return(u/(1+u))
}
# probability of fusion
p.fuse.x=function(x,params) {
  u=exp(params[23,site]+params[24,site]*x)
  return(u/(1+u))
}
#probability of large fusion  
p.fuseLG.x=function(x,params) {
  u=exp(params[25,site]+params[26,site]*x)
  return(u/(1+u))
}
#probability of lo fiss product number  
p.hilo.x=function(x,params) {
  u=exp(params[35,site]+params[36,site]*x)
  return(u/(1+u))
}
#Fission offspring using multinomial response
clo.yx=function(xp,x,params) {
  sigmaxp2<-params[29,site]*exp(2*params[30,site]*x)
  sigmaxp<-sqrt(sigmaxp2)
  nfiss<-exp(params[31,site]+params[32,site]*x)
  ##return: prob of fission * prob of lo number*number of fission products * prob that a fisser of size x prodcues a product of size xp
  p.fiss.x(x,params)*p.hilo.x(x,params)*nfiss*dnorm(xp,params[27,site]+params[28,site]*x,sigmaxp)
}
chi.yx=function(xp,x,params) {
  sigmaxp2<-params[29,site]*exp(2*params[30,site]*x)
  sigmaxp<-sqrt(sigmaxp2)
  nfiss<-exp(params[33,site]+params[34,site]*x)
  ##return: prob of fission * 1-prob of lo number*number of fission products * prob that a fisser of size x prodcues a product of size xp
  p.fiss.x(x,params)*(1-p.hilo.x(x,params))*nfiss*dnorm(xp,params[27,site]+params[28,site]*x,sigmaxp)
}
# number of clones per adult for "constant correction" of clonal matrix (no size distribution of offspring)
c.y=function(x,params) {
  nfiss<-exp(params[31,site]+params[32,site]*x)
  u=p.fiss.x(x,params)*nfiss
  return (u)
}

#recruitment = sexual offspring 
f.yx=function(xp,x,params){
  params[39,site]*x*dnorm(xp,params[37,site],params[38,site])
}

# recruits per adult for "constant correction" of fecundity matrix (no size distribution of offspring)
f.y=function(xp,params) {
  r=exp(params[37,site])/exp(params[40,site])
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
  Si=sI.x(y,params=params)*p.I.x(y,params) # survival 
  Sis=sIs.x(y,params=params)*p.Is.x(y,params) # survival 
  Sm=sM.x(y,params=params)*p.M.x(y,params)# survival 
  Sr=sR.x(y,params=params)*p.R.x(y,params) # survival 
  Ss=sS.x(y,params=params)*(1-p.I.x(y,params)-p.Is.x(y,params)-p.M.x(y,params)-p.R.x(y,params)-p.U.x(y,params)) # survival 
  Su=sU.x(y,params=params)*p.U.x(y,params) # survival 
  S=Si+Sis+Sm+Sr+Ss+Su
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
  
  K=P
  
  return(K)
}  
