###########################################################MeganR Fixed Effect Vital Rates and Parameters Matrix#############################################################################
library(nlme) # for gls()and gnls()
library(lme4) # for recruitment
library(VGAM)#For vglmm function with zero truncated poisson

nx<-10 # number of sites for this species
ncoef<-36 # number of coefficients
n<-300 # number of cuts for meshpoints/discretization into matrix

################################################################Vital Rate Functions###################################################################################

#survival regression
surv.reg=glm(Psurv~size+utrans,data=MC,family=binomial)
summary(surv.reg)
surv.int=c(coefficients(surv.reg)[1], coefficients(surv.reg)[1]+coefficients(surv.reg)[3:(3+nx-1-1)])
surv.slope=coefficients(surv.reg)[2]

# growth regression
growth.reg<-gnls(sizeNext~a+b*size+c*exp(-d*size),data=GONGS,params=list(a~utrans,b+c+d~1),start=list(a=c(1,1,1,1,1,1,1,1,1,1),b=1,c=1,d=1),control=gnlsControl(nlsTol=0.1))
summary(growth.reg)
growth.int<-c(coefficients(growth.reg)[1],coefficients(growth.reg)[1]+coefficients(growth.reg)[2:10])
growth.slope<-coefficients(growth.reg)[11]
growth.e1<-coefficients(growth.reg)[12]
growth.e2<-coefficients(growth.reg)[13]
growth.sd<-summary(growth.reg)$sigma

es.reg<-gls(sizeNext~size,data=GOES)
summary(es.reg)
es.int<-coefficients(es.reg)[1]
es.slope<-coefficients(es.reg)[2]
es.sd<-(summary(es.reg)$sigma)

fiss.full<-lm(sizeNext~size,data=fissonly)
fiss.int<-coefficients(fiss.full)[1]
fiss.slope<-coefficients(fiss.full)[2]
fiss.sd<-summary(fiss.full)$sigma

fuse.full<-lm(sizeNext~size,data=fuseonly)
summary(fuse.full)
fuse.int<-coefficients(fuse.full)[1]
fuse.slope<-coefficients(fuse.full)[2]
fuse.sd<-summary(fuse.full)$sigma

#probability of extreme shrinkage conditional on survival
Pes<-glm(PES~1,family=binomial,data=growthonly)
pes.int<-coefficients(Pes)

#probability of fission conditional on survival (for fission growth function and clonal matrix)
fiss=glm(Pfiss~size*utrans,data=MC,family=binomial)
summary(fiss)
Pfiss.int=c(coefficients(fiss)[1],coefficients(fiss)[1]+coefficients(fiss)[3:(nx+1)])
Pfiss.slope=c(coefficients(fiss)[2],coefficients(fiss)[2]+coefficients(fiss)[(nx+2):(2*nx)])

# probability of fusion conditional on survival(for fusion growth function)
fuse=glm(Pfuse~size+utrans,data=MC,family=binomial)
summary(fuse)
Pfuse.int=c(coefficients(fuse)[1],coefficients(fuse)[1]+coefficients(fuse)[3:(nx+1)])
Pfuse.slope=coefficients(fuse)[2]

fuseLG=glm(PfuseLG~size+utrans,data=MC,family=binomial)
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
clonesize=gls(sizeNext~size,data=MCfs,weights=varExp(form=~size))
summary(clonesize)
clonesize.int=coefficients(clonesize)[1]
clonesize.slope=coefficients(clonesize)[2]
clonesize.sd2=(summary(clonesize)$sigma)^2
clonesize.varFunc=coef(clonesize$modelStruct$varStruct,"exp")

##recruitment rate as in Bruno et al. 2011 page 130
#number of recruits per quad.  Includes all years not just the years including years excluded from growth due to 2 year change etc.
rcrt.lmeU.Y <- glmer(recruits~utrans+(1|year),family=poisson,data=MCr)
summary(rcrt.lmeU.Y)
recruit.int= c(fixef(rcrt.lmeU.Y)[1],fixef(rcrt.lmeU.Y)[1]+fixef(rcrt.lmeU.Y)[2:nx])


##Recruitment rate using averages of observed recruits and colony density (not using lm just averaged by hand a while ago so data may have changed)
#recruit.int=0.0905897

##Tried different recruitment rates will think about more with sensitivity analyses
#recruit.int=0.01

#size distribution of recruits
rec.sizeU<-lm(sizeNext~utrans,data=MCrs)
summary(rec.sizeU)
recruit.size.mean=c(coefficients(rec.sizeU)[1],coefficients(rec.sizeU)[1]+coefficients(rec.sizeU)[2:nx])
recruit.size.sd=(summary(rec.sizeU)$sigma)
#recruit.size.mean=mean(MC$sizeNext[MC$fateIPM=="recruit"])
#recruit.size.sd=sd(MC$sizeNext[MC$fateIPM=="recruit"])

#number of colonies per quadIncludes all years not just the years including years excluded from growth due to 2 year change etc.
dens.lmeU.Y <- glmer(colnum~utrans+(1|year),family=poisson,data=MCd)
summary(dens.lmeU.Y)
dens.int=c(fixef(dens.lmeU.Y)[1],fixef(dens.lmeU.Y)[1]+fixef(dens.lmeU.Y)[2:nx])

# Size distribution of colonies colonies present at time t (all but recruits and fissSM)
col.size<-lm(sizeNext~utrans,data=MCds)
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
params<-matrix(0, ncoef, nx)  #creates a matrix of zeros with ncoef rows and nx (num of sites) columns

params[1,]<-c(surv.int)               #survival intercept by site
params[2,]<-rep(surv.slope)           #survival slope wrt size, same across sites
params[3,]<-c(growth.int)             #growth linear intercept by site
params[4,]<-rep(growth.slope)         #growth linear slope by site
params[5,]<-rep(growth.e1)            #growth exponential intercept (constant coefficient),same across sites
params[6,]<-rep(growth.e2)            #growth exponential slope, same across sites
params[7,]<-rep(growth.sd)            #growth standard deviation, same across sites
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
params[18,]<-c(Pfiss.int)             #Probability of fission intercept by site
params[19,]<-c(Pfiss.slope)           #Probability of fission slopeby site
params[20,]<-c(Pfuse.int)             #Probability of fusion intercept by site
params[21,]<-rep(Pfuse.slope)         #Probability of fusion slope, same across site
params[22,]<-c(PfuseLG.int)           #Probability of being large fusion intercept given that you survived and fused by site
params[23,]<-rep(PfuseLG.slope)       #probability of being a large fusion slope given that you survived and fused same across sites
params[24,]<-rep(clonesize.int)       #clone size intercept, same across sites
params[25,]<-rep(clonesize.slope)     #clone size slope, same across sites
params[26,]<-rep(clonesize.sd2)       #clone size sd, same across sites
params[27,]<-rep(clonesize.varFunc)   #clone size variance increases with size, same across sites
params[28,]<-rep(clonenumlo.int)      #clonenum intercept lo number of fission events, same across sites 
params[29,]<-rep(clonenumlo.slope)    #clonenum slope lo number of fission events, same across sites 
params[30,]<-rep(clonenumhi.int)      #clonenum intercept hi number of fission events, same across sites 
params[31,]<-rep(clonenumhi.slope)    #clonenum slope hi number of fission events, same across sites 
params[32,]<-rep(Philo.int)           #Probability of a fisser producing a lo number of recruits intercept, same across sites
params[33,]<-rep(Philo.slope)         #Probability of a fisser producing a lo number of recruits slope, same across sites
params[34,]<-c(recruit.size.mean)     #mean size of recruitsby site
params[35,]<-rep(recruit.size.sd)     #sd size of recruits, same across sites
params[36,]<-c(r)                     #estimated number of recruits per colony by site (as in Bruno et al. 2011)

rownames(params)<-c("surv.int","surv.slope","growth.int","growth.slope","growth.e1","growth.e2","growth.sd","Pes","es.int","es.slope","es.sd","fiss.int","fiss.slope","fiss.sd","fuse.int","fuse.slope","fuse.sd","Pfiss.int","Pfiss.slope","Pfuse.int","Pfuse.slope","PfuseLG.int","PfuseLG.slope","clonesize.int","clonesize.slope","clonesize.sd2","clonesize.varFunc","clonenumlo.int","clonenumlo.slope","clonenumhi.int","clonenumhi.slope","Philo.int","Philo.slope","recruit.size.mean","recruit.size.sd","r")
