###########################################################MeganR IPM Code #############################################################################
library(nlme) # for gls()
library(lme4) # for recruitment
library(VGAM)#For vglmm function with zero truncated poisson
library(aqfig)# for vertical.image.legend
library(popbio) # for eigen.analysis()
library(MASS)# for negative binomial distribution
library(countreg) # for zerotrunc for geometric distribution
library(nnet) # for multinomial distribution

MCg<-within(MCg,fate<-relevel(fate,ref="growth")) #makes fate "growth" the zero-level in regression
dat<-MC #data for this species
growdat<-MCg #growth data for this species
fissdat<-MCf #fission only data for this species
fisssizedat<-MCfs #small fission only data for this species
recdat<-MCr # number of recruits per quad per year
recsizedat<-MCrs # recruits only for size distribution
densdat<-MCd #number of colonies per quad per year
denssizedat<-MCd # colonies present at time t (all but recruits and fissSM) for sfd
nx<-10 # number of sites for this species
nf<-3 #number of growth fates for this species
ncoef<-50 # number of coefficients
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
#function increases too steeply over time.  Don't want the probability to go to zero or to balloon at higher values of size
clonenum=vglm(numfiss~size,data=fissdat,family=pospoisson)
summary(clonenum)
clonenum.int=c(coefficients(clonenum)[1])#,coefficients(clonenum)[1]+coefficients(clonenum)[2:nx])
clonenum.slope=c(coefficients(clonenum)[2])#[nx+1],coefficients(clonenum)[nx+1]+coefficients(clonenum)[(nx+2):(2*nx)])

#number of clones per fission event: negative binomial distribution
clonenb <- glm.nb(numfiss~size,data=fissdat)
summary (clonenb)
clonenb.int<-coefficients(clonenb)[1]
clonenb.slope<-coefficients(clonenb)[2]

#number of clones per fission event: geometric dtribution
clonegeo=zerotrunc(numfiss~size,data=fissdat,dist="geometric")
summary(clonegeo)
clonegeo.int<-coefficients(clonegeo)[1]
clonegeo.slope<-coefficients(clonegeo)[2]

#Multinomial probability of 1,2,3 or more fission products from a fission event
#Backup in case we can't find an adequate replacement for zero truncated poisson
clonemulti=multinom(fnumfiss~size, data=fissdat)
summary(clonemulti)
clonemulti.int.1=coefficients(clonemulti)[1]
clonemulti.slope.1=coefficients(clonemulti)[4]
clonemulti.int.2=coefficients(clonemulti)[3]
clonemulti.slope.2=coefficients(clonemulti)[6]
clonemulti.int.3=coefficients(clonemulti)[2]
clonemulti.slope.3=coefficients(clonemulti)[5]
##Testing fits for 5-10
clonemulti2=multinom(fnumfiss2~size, data=fissdat)
summary(clonemulti2)
clonemulti2.int.1=coefficients(clonemulti2)[2]
clonemulti2.slope.1=coefficients(clonemulti2)[9]
clonemulti2.int.2=coefficients(clonemulti2)[7]
clonemulti2.slope.2=coefficients(clonemulti2)[14]
clonemulti2.int.3=coefficients(clonemulti2)[6]
clonemulti2.slope.3=coefficients(clonemulti2)[13]
clonemulti2.int.4=coefficients(clonemulti2)[1]
clonemulti2.slope.4=coefficients(clonemulti2)[8]
clonemulti2.int.6=coefficients(clonemulti2)[4]
clonemulti2.slope.6=coefficients(clonemulti2)[11]
clonemulti2.int.7=coefficients(clonemulti2)[3]
clonemulti2.slope.7=coefficients(clonemulti2)[10]
clonemulti2.int.10=coefficients(clonemulti2)[5]
clonemulti2.slope.10=coefficients(clonemulti2)[12]

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

#number of colonies per quadIncludes all years not just the years including years excluded from growth due to 2 year change etc.
dens.lmeU.Y <- glmer(colnum~utrans+(1|year),family=poisson,data=densdat)
summary(dens.lmeU.Y)
dens.int=c(fixef(dens.lmeU.Y)[1],fixef(dens.lmeU.Y)[1]+fixef(dens.lmeU.Y)[2:nx])

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
params[12,]<-rep(Pfiss.slope)          #Probability of fission slope, same across sites
params[13,]<-c(Pfuse.int)              #Probability of fusion intercept by site
params[14,]<-c(Pfuse.slope)            #Probability of fusion slope by site
params[15,]<-c(PfuseLG.int)            #Probability of being large fusion intercept given that you survived and fused by site
params[16,]<-rep(PfuseLG.slope)        #probability of being a large fusion slope given that you survived and fused same across sites
params[17,]<-rep(clonesize.int)        #clone size intercept, same across sites
params[18,]<-rep(clonesize.slope)      #clone size slope, same across sites
params[19,]<-rep(clonesize.sd2)        #clone size sd, same across sites
params[20,]<-rep(clonesize.varFunc)    #clone size variance increases with size by site
params[21,]<-rep(clonenum.int)         #number of clones intercept, same across sites using zero truncated poisson
params[22,]<-rep(clonenum.slope)       #number of clones slope same for fates and across sites using zerotruncated poisson
params[23,]<-c(recruit.int)            #number of recruits (kind of made up), same across sites     
params[24,]<-c(recruit.size.mean)      #mean size of recruits, same across sites
params[25,]<-rep(recruit.size.sd)      #sd size of recruits, same across sites
params[26,]<-c(dens.int)               #colony density by quad   
params[27,]<-rep(clonemulti.int.1)     #clonenum intercept using multinomial distribution: 1clonal product
params[28,]<-rep(clonemulti.slope.1)   #clonenum slope using multinomial distribution: 1 clonal product
params[29,]<-rep(clonemulti.int.2)     #clonenum intercept using multinomial distribution: 2 clonal products
params[30,]<-rep(clonemulti.slope.2)   #clonenum slope using multinomial distribution: 2 clonal products
params[31,]<-rep(clonemulti.int.3)     #clonenum intercept using multinomial distribution: 3clonal products
params[32,]<-rep(clonemulti.slope.3)   #clonenum slope using multinomial distribution: 3 clonal products
params[33,]<-rep(clonenb.int)          #clonenum intercept using multinomial distribution: >3 clonal products
params[34,]<-rep(clonenb.slope)        #clonenum slope using multinomial distribution: >3 clonal products
params[35,]<-rep(clonegeo.int)         #clonenum intercept using geometric distribution
params[36,]<-rep(clonegeo.slope)       #clonenum slope using geometric distribution
params[37,]<-rep(clonemulti2.int.1)     #clonenum intercept using multinomial distribution: 1clonal product
params[38,]<-rep(clonemulti2.slope.1)   #clonenum slope using multinomial distribution: 1 clonal product
params[39,]<-rep(clonemulti2.int.2)     #clonenum intercept using multinomial distribution: 2 clonal products
params[40,]<-rep(clonemulti2.slope.2)   #clonenum slope using multinomial distribution: 2 clonal products
params[41,]<-rep(clonemulti2.int.3)     #clonenum intercept using multinomial distribution: 3clonal products
params[42,]<-rep(clonemulti2.slope.3) 
params[43,]<-rep(clonemulti2.int.4)
params[44,]<-rep(clonemulti2.slope.4)
params[45,]<-rep(clonemulti2.int.6)
params[46,]<-rep(clonemulti2.slope.6)
params[47,]<-rep(clonemulti2.int.7)
params[48,]<-rep(clonemulti2.slope.7)
params[49,]<-rep(clonemulti2.int.10)
params[50,]<-rep(clonemulti2.slope.10)
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


#Fission offspring using zero truncated poisson
  c.yx=function(xp,x,params) {
    sigmaxp2<-params[19,site]*exp(2*params[20,site]*x)
    sigmaxp<-sqrt(sigmaxp2)
    nfiss<-exp(params[21,site]+params[22,site]*x)
    ##return: prob of fission * number of fission products * prob that a fisser of size x prodcues a product of size xp
    p.fiss.x(x,params)*nfiss*dnorm(xp,params[17,site]+params[18,site]*x,sigmaxp)
  }
#Fission offspring using negative binomial distribution
cnb.yx=function(xp,x,params) {
  sigmaxp2<-params[19,site]*exp(2*params[20,site]*x)
  sigmaxp<-sqrt(sigmaxp2)
  nfiss<-exp(params[33,site]+params[34,site]*x)
  ##return: prob of fission * number of fission products * prob that a fisser of size x prodcues a product of size xp
  p.fiss.x(x,params)*nfiss*dnorm(xp,params[17,site]+params[18,site]*x,sigmaxp)
}
#Fission offspring using geometric distribution
cgeo.yx=function(xp,x,params) {
  sigmaxp2<-params[19,site]*exp(2*params[20,site]*x)
  sigmaxp<-sqrt(sigmaxp2)
  nfiss<-exp(params[35,site]+params[36,site]*x)
  ##return: prob of fission * number of fission products * prob that a fisser of size x prodcues a product of size xp
  p.fiss.x(x,params)*nfiss*dnorm(xp,params[17,site]+params[18,site]*x,sigmaxp)
}
#Fission offspring using multinomial response
  p.c1.x=function(x,params) {
    u<-exp(params[27,site]+params[28,site]*x)
    ud<-exp(params[27,site]+params[28,site]*x)+exp(params[29,site]+params[30,site]*x)+exp(params[31,site]+params[32,site]*x)+1    
    return(u/ud)    
  }
  p.c2.x=function(x,params) {
    u<-exp(params[29,site]+params[30,site]*x)
    ud<-exp(params[27,site]+params[28,site]*x)+exp(params[29,site]+params[30,site]*x)+exp(params[31,site]+params[32,site]*x)+1    
    return(u/ud)  
  }
  p.c3.x=function(x,params) {
    u<-exp(params[31,site]+params[32,site]*x)
    ud<-exp(params[27,site]+params[28,site]*x)+exp(params[29,site]+params[30,site]*x)+exp(params[31,site]+params[32,site]*x)+1    
    return(u/ud)  
  }  
####
  p2.c1.x=function(x,params) {
    u<-exp(params[37,site]+params[38,site]*x)
    ud<-exp(params[37,site]+params[38,site]*x)+exp(params[39,site]+params[40,site]*x)+exp(params[41,site]+params[42,site]*x)+exp(params[43,site]+params[44,site]*x)+exp(params[45,site]+params[46,site]*x)+exp(params[47,site]+params[48,site]*x)+exp(params[49,site]+params[50,site]*x)+1    
    return(u/ud)    
  }
  p2.c2.x=function(x,params) {
    u<-exp(params[39,site]+params[40,site]*x)
    ud<-exp(params[37,site]+params[38,site]*x)+exp(params[39,site]+params[40,site]*x)+exp(params[41,site]+params[42,site]*x)+exp(params[43,site]+params[44,site]*x)+exp(params[45,site]+params[46,site]*x)+exp(params[47,site]+params[48,site]*x)+exp(params[49,site]+params[50,site]*x)+1
    return(u/ud)  
  }
  p2.c3.x=function(x,params) {
    u<-exp(params[41,site]+params[42,site]*x)
    ud<-exp(params[37,site]+params[38,site]*x)+exp(params[39,site]+params[40,site]*x)+exp(params[41,site]+params[42,site]*x)+exp(params[43,site]+params[44,site]*x)+exp(params[45,site]+params[46,site]*x)+exp(params[47,site]+params[48,site]*x)+exp(params[49,site]+params[50,site]*x)+1
    return(u/ud)  
  }
  p2.c4.x=function(x,params) {
    u<-exp(params[43,site]+params[44,site]*x)
    ud<-exp(params[37,site]+params[38,site]*x)+exp(params[39,site]+params[40,site]*x)+exp(params[41,site]+params[42,site]*x)+exp(params[43,site]+params[44,site]*x)+exp(params[45,site]+params[46,site]*x)+exp(params[47,site]+params[48,site]*x)+exp(params[49,site]+params[50,site]*x)+1
    return(u/ud)      
  }
  p2.c6.x=function(x,params) {
    u<-exp(params[45,site]+params[46,site]*x)
    ud<-exp(params[37,site]+params[38,site]*x)+exp(params[39,site]+params[40,site]*x)+exp(params[41,site]+params[42,site]*x)+exp(params[43,site]+params[44,site]*x)+exp(params[45,site]+params[46,site]*x)+exp(params[47,site]+params[48,site]*x)+exp(params[49,site]+params[50,site]*x)+1
    return(u/ud)  
  }
  p2.c7.x=function(x,params) {
    u<-exp(params[47,site]+params[48,site]*x)
    ud<-exp(params[37,site]+params[38,site]*x)+exp(params[39,site]+params[40,site]*x)+exp(params[41,site]+params[42,site]*x)+exp(params[43,site]+params[44,site]*x)+exp(params[45,site]+params[46,site]*x)+exp(params[47,site]+params[48,site]*x)+exp(params[49,site]+params[50,site]*x)+1
    return(u/ud)  
  }
  p2.c10.x=function(x,params) {
    u<-exp(params[49,site]+params[50,site]*x)
    ud<-exp(params[37,site]+params[38,site]*x)+exp(params[39,site]+params[40,site]*x)+exp(params[41,site]+params[42,site]*x)+exp(params[43,site]+params[44,site]*x)+exp(params[45,site]+params[46,site]*x)+exp(params[47,site]+params[48,site]*x)+exp(params[49,site]+params[50,site]*x)+1
    return(u/ud)  
  }
####
  c1.yx=function(xp,x,params) {
    sigmaxp2<-params[19,site]*exp(2*params[20,site]*x)
    sigmaxp<-sqrt(sigmaxp2)
    p.c1.x(x,params)*dnorm(xp,params[17,site]+params[18,site]*x,sigmaxp)
  }
  c2.yx=function(xp,x,params) {
    sigmaxp2<-params[19,site]*exp(2*params[20,site]*x)
    sigmaxp<-sqrt(sigmaxp2)
    p.c2.x(x,params)*2*dnorm(xp,params[17,site]+params[18,site]*x,sigmaxp)
  }
  c3.yx=function(xp,x,params) {
    sigmaxp2<-params[19,site]*exp(2*params[20,site]*x)
    sigmaxp<-sqrt(sigmaxp2)
    p.c3.x(x,params)*3*dnorm(xp,params[17,site]+params[18,site]*x,sigmaxp)
  }
  c4.yx=function(xp,x,params) {
    sigmaxp2<-params[19,site]*exp(2*params[20,site]*x)
    sigmaxp<-sqrt(sigmaxp2)
    (1-p.c1.x(x,params)-p.c2.x(x,params)-p.c3.x(x,params))*4*dnorm(xp,params[17,site]+params[18,site]*x,sigmaxp)
  }
######
  c1.2yx=function(xp,x,params) {
    sigmaxp2<-params[19,site]*exp(2*params[20,site]*x)
    sigmaxp<-sqrt(sigmaxp2)
    p.fiss.x(x,params)*p2.c1.x(x,params)*dnorm(xp,params[17,site]+params[18,site]*x,sigmaxp)
  }
  c2.2yx=function(xp,x,params) {
    sigmaxp2<-params[19,site]*exp(2*params[20,site]*x)
    sigmaxp<-sqrt(sigmaxp2)
    p.fiss.x(x,params)*p2.c2.x(x,params)*2*dnorm(xp,params[17,site]+params[18,site]*x,sigmaxp)
  }
  c3.2yx=function(xp,x,params) {
    sigmaxp2<-params[19,site]*exp(2*params[20,site]*x)
    sigmaxp<-sqrt(sigmaxp2)
    p.fiss.x(x,params)*p2.c3.x(x,params)*3*dnorm(xp,params[17,site]+params[18,site]*x,sigmaxp)
  }
  c4.2yx=function(xp,x,params) {
    sigmaxp2<-params[19,site]*exp(2*params[20,site]*x)
    sigmaxp<-sqrt(sigmaxp2)
    p.fiss.x(x,params)*p2.c4.x(x,params)*4*dnorm(xp,params[17,site]+params[18,site]*x,sigmaxp)
  }
  c5.2yx=function(xp,x,params) {
    sigmaxp2<-params[19,site]*exp(2*params[20,site]*x)
    sigmaxp<-sqrt(sigmaxp2)
    p.fiss.x(x,params)*(1-p2.c1.x(x,params)-p2.c2.x(x,params)-p2.c3.x(x,params)-p2.c4.x(x,params)-p2.c6.x(x,params)-p2.c7.x(x,params)-p2.c10.x(x,params))*5*dnorm(xp,params[17,site]+params[18,site]*x,sigmaxp)
  }
  c6.2yx=function(xp,x,params) {
    sigmaxp2<-params[19,site]*exp(2*params[20,site]*x)
    sigmaxp<-sqrt(sigmaxp2)
    p.fiss.x(x,params)*p2.c6.x(x,params)*6*dnorm(xp,params[17,site]+params[18,site]*x,sigmaxp)
  }
  c7.2yx=function(xp,x,params) {
    sigmaxp2<-params[19,site]*exp(2*params[20,site]*x)
    sigmaxp<-sqrt(sigmaxp2)
    p.fiss.x(x,params)*p2.c7.x(x,params)*7*2*dnorm(xp,params[17,site]+params[18,site]*x,sigmaxp)
  }
  c10.2yx=function(xp,x,params) {
    sigmaxp2<-params[19,site]*exp(2*params[20,site]*x)
    sigmaxp<-sqrt(sigmaxp2)
    p.fiss.x(x,params)*p2.c10.x(x,params)*10*3*dnorm(xp,params[17,site]+params[18,site]*x,sigmaxp)
  }
 cs.yx=function(xp,x,params){
   sigmaxp2<-params[19,site]*exp(2*params[20,site]*x)
   sigmaxp<-sqrt(sigmaxp2)
   dnorm(xp,params[17,site]+params[18,site]*x,sigmaxp)
 }
# number of clones per adult for "constant correction" of clonal matrix (no size distribution of offspring)
  c.y=function(x,params) {
    nfiss<-exp(params[21,site]+params[22,site]*x)
    u=p.fiss.x(x,params)*nfiss
    return (u)
  }

#recruitment = sexual offspring 
  f.yx=function(xp,x,params) {
    r=exp(params[23,site])/exp(params[26,site])
    u=dnorm(xp,params[24,site],params[25,site])*r
    return (u)
 }
# recruits per adult for "constant correction" of fecundity matrix (no size distribution of offspring)
  f.y=function(xp,params) {
    r=exp(params[23,site])/exp(params[26,site])
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
  #with zero truncated poisson
  C=h*outer(y,y,c.yx,params=params)
  
  #with negative binomial
  Cnb=h*outer(y,y,cnb.yx,params=params)
  
  #with geometric
  Cgeo=h*outer(y,y,cgeo.yx,params=params)
  
  #with multinomial response
  PI=p.fiss.x(y,params)
  PC1=p.c1.x(y,params)*1
  PC2=p.c2.x(y,params)*2
  PC3=p.c3.x(y,params)*3
  PC4=(1-p.c1.x(y,params)-p.c2.x(y,params)-p.c3.x(y,params))*4
  #C1=h*outer(y,y,c1.yx,params=params)
  PC=PC1+PC2+PC3+PC4
  #C11=C1
  #for (i in 1:n) C11[,i]=C1[,i]*PI[i]
  #C2=h*outer(y,y,c2.yx,params=params)
  #C12=C2
  #for (i in 1:n) C12[,i]=C2[,i]*PI[i]
  #C3=h*outer(y,y,c3.yx,params=params)
  #C13=C3
  #for (i in 1:n) C13[,i]=C3[,i]*PI[i]
  #C4=h*outer(y,y,c4.yx,params=params)
  #C14=C4
  #for (i in 1:n) C14[,i]=C4[,i]*PI[i]
  #Cmulti=C11+C12+C13+C14
  #Cmulti=C1+C2+C3+C4
  Cmulti=h*outer(y,y,cs.yx,params)
  Cmulti1=Cmulti
  for (i in 1:n) Cmulti1[,i]=Cmulti[,i]*PI[i]*PC[i]
  
#F is the fecundity/sexually reproduced offspring kernel
  #R sets minimum reproductive size based on outer 1-2 cm not reproducing used area of a circle to determine smalles reproductive size = pi 2 squared.
  #R=c(rep(0,times=121),rep(1,times=179))
  #F is currently using the fecundity constant calculated based on Bruno et al. 2011
  F=h*outer(y,y,f.yx,params=params)
  #for(i in 1:n) F[,i]=F[,i]*R[i]
  
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
  f<-vector(length=n)
  for(i in 1:n) f[i]<- f.y(y,params=params)
  correction.here <- f/colSums(F)
  F<- t(t(F) * correction.here)
  
  #for Cmatrix
  #c<-c.y(y,params=params)
  #correction.here <- c/colSums(C)
  #C<- t(t(C) * correction.here)
    
  K=P+Cmulti1
  
  return(K)
}  
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
xvals<-barplot(eig)
text(xvals, par("usr")[3] - 0.25, srt = 45, adj = 1,labels = unique(dat$utrans), xpd = TRUE)

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
for (i in c(1,4,9,15)){
  barplot(elasP[,i],horiz=TRUE,xlab="Elasticity",col=cl[1:nx])
}
par(mfrow = c(2, 2))
par(mar = c(7, 2, 2, 2) + 0.1)
for (i in c(1:nx)){
barplot(elasP[i,],main="Elasticity",names.arg=rownames(params),las=2)
}

#legend("bottom",legend=sort(unique(utrans),decreasing=TRUE),horiz=TRUE,pch=1,col=cl[ns:1],cex=1,bty="n")







