################MeganR Fixed Effect Vital Rates and Parameters Matrix#############################################################################
library(nlme) # for gls(), lme(), nlmer()
library(lme4) # glmer()
#library(MASS)
library(nnet) #for multinom() for clonenum
nx<-11 # number of sites for this species
ncoef<-35 # number of coefficients

###############Vital Rate Functions###################################################################################

#survival regression
surv.reg.re0<-glmer(Psurv~size+(1|utrans),data=MC,family=binomial)
summary(surv.reg.re0)
surv.int.re0=c(fixed.effects(surv.reg.re0)[1],fixed.effects(surv.reg.re0)[1]+random.effects(surv.reg.re0)[[1]][1:10,1])
surv.slope.re0=fixed.effects(surv.reg.re0)[2]

# growth regression
growth.reg.re.a0<-nlme(sizeNext~a+b*size+c*exp(-d*size),data=GONGS,fixed=a+b+c+d~1,random=c~1|utrans,start=list(fixed=c(1,1,1,1)),control=gnlsControl(nlsTol=0.1))
growth.int.a0<-fixed.effects(growth.reg.re.a0)[1]
growth.slope.a0<-fixed.effects(growth.reg.re.a0)[2]
growth.e1.a0<-c(fixed.effects(growth.reg.re.a0)[3],fixed.effects(growth.reg.re.a0)[3]+random.effects(growth.reg.re.a0)[1:10,1])
growth.e2.a0<-fixed.effects(growth.reg.re.a0)[4]
growth.sd.a0<-summary(growth.reg.re.a0)$sigma

es.reg<-gls(sizeNext~size,data=GOES)
summary(es.reg)
es.int<-coefficients(es.reg)[1]
es.slope<-coefficients(es.reg)[2]
es.sd<-(summary(es.reg)$sigma)

fiss.full<-gls(sizeNext~size,data=fissonly)
fiss.int<-coefficients(fiss.full)[1]
fiss.slope<-coefficients(fiss.full)[2]
fiss.sd<-summary(fiss.full)$sigma

fuse.full<-gls(sizeNext~size,data=fuseonly)
summary(fuse.full)
fuse.int<-coefficients(fuse.full)[1]
fuse.slope<-coefficients(fuse.full)[2]
fuse.sd<-summary(fuse.full)$sigma

#probability of extreme shrinkage conditional on survival
Pes<-glm(PES~1,family=binomial,data=growthonly)
pes.int<-coefficients(Pes)

#probability of fission conditional on survival (for fission growth function and clonal matrix)
fiss.re0=glmer(Pfiss~size+(1|utrans),data=MC,family=binomial)
summary(fiss.re0)
Pfiss.int.re0=c(fixed.effects(fiss.re0)[1],fixed.effects(fiss.re0)[1]+random.effects(fiss.re0)[[1]][1:10,1])
Pfiss.slope.re0=fixed.effects(fiss.re0)[2]

# probability of fusion conditional on survival(for fusion growth function)
fuse.re0<-glmer(Pfuse~size+(1|utrans),data=MC,family=binomial)
summary(fuse.re0)
Pfuse.int.re0=c(fixed.effects(fuse.re0)[1],fixed.effects(fuse.re0)[1]+random.effects(fuse.re0)[[1]][1:10,1])
Pfuse.slope.re0=fixed.effects(fuse.re0)[2]

fuseLG.re0=glmer(PfuseLG~size+(1|utrans),data=MC,family=binomial)
summary(fuseLG.re0)
PfuseLG.int.re0=c(fixed.effects(fuseLG.re0)[1],fixed.effects(fuseLG.re0)[1]+random.effects(fuseLG.re0)[[1]][1:10,1])
PfuseLG.slope.re0=fixed.effects(fuseLG.re0)[2]

#number of clones per fission event: 2 functions for data subset with hi (4-10) and lo (1-3) numfiss using zero truncated poisson
clonenum=multinom(fnumfiss~size,data=MCf)
summary(clonenum)
clonenum1.int=coefficients(clonenum)[1,1]
clonenum1.slope=coefficients(clonenum)[1,2]
clonenum2.int=coefficients(clonenum)[3,1]
clonenum2.slope=coefficients(clonenum)[3,2]
clonenum3.int=coefficients(clonenum)[2,1]
clonenum3.slope=coefficients(clonenum)[2,2]


##############################################################################################################################
#First Attempts at estimating rates of recruitment by density
##Recruitment rate using averages of observed recruits and colony density (not using lm just averaged by hand a while ago so data may have changed)
#recruit.int=0.0905897
##Tried different recruitment rates will think about more with sensitivity analyses
#recruit.int=0.01
#recruit.size.mean=mean(MC$sizeNext[MC$fateIPM=="recruit"])
#recruit.size.sd=sd(MC$sizeNext[MC$fateIPM=="recruit"])
###############################################################################################################################
# Second attempt at estimating recruitment
##based on Bruno et al. 2011 page 130
#Used density, size distribution of 'adult colonies' 
#Basically estimating coral cover.  But we already have this so better to use measured values.  Spoke with MeganD about it.
#number of colonies per quadIncludes all years not just the years including years excluded from growth due to 2 year change etc.
#dens.lmeU.Y <- glmer(colnum~1+(1|utrans),family=poisson,data=MCd)
#summary(dens.lmeU.Y)
#dens.int=c(fixef(dens.lmeU.Y)[1],fixef(dens.lmeU.Y)[1]+ranef(dens.lmeU.Y)[[1]][1:10,1])

# Size distribution of colonies colonies present at time t (all but recruits and fissSM)
#col.size<-lme(size~1,random=~1|utrans,data=MCds)
#summary(col.size)
#colsize.int=c(fixef(col.size)[1],fixef(col.size)[1]+ranef(col.size)[1:10,1])
#colsize.sd=(summary(col.size)$sigma)

#r<-c(rep(0,nx))
#for(i in 1:nx){
# recruitdensity=exp(recruit.int[i])
#colonydensity=exp(dens.int[i])
#sfd= dnorm(y,colsize.int[i],colsize.sd)
#r[i]<-recruitdensity/sum(sfd*y*colonydensity)
#}
#################################################################################################################################
#Third attemp at estimating recruitment rate by colony size using measured coral cover (See notes above)
#number of recruits per quad.  Includes all years not just the years including years excluded from growth due to 2 year change etc.
rcrt.lmeU.Y <- glmer(recruits~1+(1|utrans),family=poisson,data=MCr)
summary(rcrt.lmeU.Y)
recruit.int= c(fixef(rcrt.lmeU.Y)[1],fixef(rcrt.lmeU.Y)[1]+ranef(rcrt.lmeU.Y)[[1]][1:10,1])

#Cover by quad added 4/28/14 by MR.  Using MCr
#want to look at relationship between log transformed area and # of recruits added 5/15/14 by MR
#Transfrom spcover_t
MCr$lnspcov<-log(MCr$spcover_t)
#remove zeros so that we can use log transformation??????????????
#not sure about this have MeganD look it over.
test<-subset(MCr,MCr$spcover_t!=0)
#Used to be spcover_t~1
#Changed to lnspcov~1 1/10/15
#This increases the estimate of r
rcov<-lme(lnspcov~1,random=~1|utrans,data=test)
#rcov<-lme(spcover_t~1,random=~1|utrans,data=test)
summary(rcov)
recruit.cov<-c(fixef(rcov)[1],fixef(rcov)[1]+ranef(rcov)[1:10,1])
#size distribution of recruits
rec.sizeU<-lme(sizeNext~1,random=~1|utrans,data=MCrs)
summary(rec.sizeU)
recruit.size.mean=c(fixef(rec.sizeU)[1],fixef(rec.sizeU)[1]+ranef(rec.sizeU)[1:10,1])
recruit.size.sd=(summary(rec.sizeU)$sigma)

#This is what goes into params
r<-c(rep(0,nx))
for(i in 1:nx){
recruitnumber=exp(recruit.int[i])
cover=recruit.cov[i]
r[i]<-recruitnumber/cover
}
####################################################
#on 2/27/15 MeganD and I had a meeting where we decided to fit size dependent and size independent models.
#If there is no difference then we will pick one and say that the difference did not matter
#Above, we went back and forth about whether we would expect colony fecundity to be a function of surface area or log-transformed surface area
#un-transformed surface area was unrealistic
#We will go with log-transformed cover to calculate closed population r
#for size independent, we still need the number of recruits to be a function of idividual colonies.  
#Size independent will be a function of colony density
#here I will cut and paste code from scenario 2 above.  But we will not be using the sfd as in Bruno
#number of colonies per quadIncludes all years including years excluded from growth due to 2 year change etc.
#dens.lmeU.Y <- glmer(colnum~1+(1|utrans),family=poisson,data=MCd)
#summary(dens.lmeU.Y)
#dens.int=c(fixef(dens.lmeU.Y)[1],fixef(dens.lmeU.Y)[1]+ranef(dens.lmeU.Y)[[1]][1:10,1])

#This is what goes into params
#r<-c(rep(0,nx))
#for(i in 1:nx){
 # recruitnumber=exp(recruit.int[i])
  #density=exp(dens.int[i])
  #r[i]<-recruitnumber/density
#}


##normal distribution for clonal offspring
clonesize.normal<-gls(sizeNext~size,data=MCfs,weights=varFixed(~size))
clone.size.int<-coefficients(clonesize.normal)[1]
clone.size.slope<-coefficients(clonesize.normal)[2]
clone.size.var<-summary(clonesize.normal)$sigma^2

#make slots for parameters
params<-matrix(0, ncoef, nx)  #creates a matrix of zeros with ncoef rows and nx (num of sites) columns

params[1,]<-c(surv.int.re0)               #survival intercept by site
params[2,]<-rep(surv.slope.re0)           #survival slope wrt size, same across sites
params[3,]<-rep(growth.int.a0)             #growth linear intercept by site
params[4,]<-rep(growth.slope.a0)         #growth linear slope by site
params[5,]<-c(growth.e1.a0)            #growth exponential intercept (constant coefficient),same across sites
params[6,]<-rep(growth.e2.a0)            #growth exponential slope, same across sites
params[7,]<-rep(growth.sd.a0)            #growth standard deviation, same across sites
params[8,]<-rep(es.int)               #Extreme negative growth intercept, same across sites
params[9,]<-rep(es.slope)             #Extreme negative growth slope, same across sites
params[10,]<-rep(es.sd)               #Extreme negative growth standard deviation, same across sites
params[11,]<-rep(fiss.int)            #fiss growth intercept by site
params[12,]<-rep(fiss.slope)          #fiss growth slope by site
params[13,]<-rep(fiss.sd)             #fiss growth standard deviation, same across site
params[14,]<-rep(fuse.int)            #fuse growth intercept by site
params[15,]<-rep(fuse.slope)          #fuse growth slope, same across sites
params[16,]<-rep(fuse.sd)             #fuse growth standard deviation, same across sites
params[17,]<-rep(pes.int)             #probability of undergoing extreme negative growth, same across sites
params[18,]<-c(Pfiss.int.re0)             #Probability of fission intercept by site
params[19,]<-rep(Pfiss.slope.re0)           #Probability of fission slopeby site
params[20,]<-c(Pfuse.int.re0)             #Probability of fusion intercept by site
params[21,]<-rep(Pfuse.slope.re0)         #Probability of fusion slope, same across site
params[22,]<-c(PfuseLG.int.re0)           #Probability of being large fusion intercept given that you survived and fused by site
params[23,]<-rep(PfuseLG.slope.re0)       #probability of being a large fusion slope given that you survived and fused same across sites
params[24,]<-rep(clonenum1.int)      #clonenum intercept lo number of fission events, same across sites 
params[25,]<-rep(clonenum1.slope)    #clonenum slope lo number of fission events, same across sites 
params[26,]<-rep(clonenum2.int)      #clonenum intercept lo number of fission events, same across sites 
params[27,]<-rep(clonenum2.slope)    #clonenum slope lo number of fission events, same across sites 
params[28,]<-rep(clonenum3.int)      #clonenum intercept lo number of fission events, same across sites 
params[29,]<-rep(clonenum3.slope)    #clonenum slope lo number of fission events, same across sites 
params[30,]<-c(recruit.size.mean)     #mean size of recruitsby site
params[31,]<-rep(recruit.size.sd)     #sd size of recruits, same across sites
params[32,]<-c(r)                     #estimated number of recruits per colony by site (as in Bruno et al. 2011)
params[33,]<-c(clone.size.int)
params[34,]<-c(clone.size.slope)
params[35,]<-c(clone.size.var)
rownames(params)<-c("surv.int","surv.slope","growth.int","growth.slope","growth.e1","growth.e2","growth.sd","es.int","es.slope","es.sd","fiss.int","fiss.slope","fiss.sd","fuse.int","fuse.slope","fuse.sd","Pes","Pfiss.int","Pfiss.slope","Pfuse.int","Pfuse.slope","PfuseLG.int","PfuseLG.slope","clonenum1.int","clonenum1.slope","clonenum2.int","clonenum2.slope","clonenum3.int","clonenum3.slope","recruit.size.mean","recruit.size.sd","r","clonesize.int","clonesize.slope","clonesize.var")
