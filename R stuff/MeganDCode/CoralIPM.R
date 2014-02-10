##Coral IPM - spatial variation in coral dynamics around Maui
##   by Megan Donahue, June 2012
## Modifield from Syng IPM: harvest of Golden grass scapes in Brazil by Orou G. Gaoue (ogaoue@bio.miami.edu), March 11-12, 2010

rm(list=ls(all=TRUE)) ## Clear all


############################# PART (1) #######################
### I - Data manip and exploration #####
# Run Data Organizing code to pull in the data
source("IPM_DataOrganizing.R")

nx <- 9 #number of sites for this species
sp <- "Montipora capitata"

#Note that the subset of sites with enough data to model is different for each species
  #MC: HonN, HonS, Kahe03, Kahe07,               Molo07, Molo13, Olow03, Olow07,       Puam03
  #PC: HonN,       Kahe03, Kahe07, KanB, KanP10, Molo07, Molo13, Olow03,               Puam03, Puam13
  #PL: HonN, HonS, Kahe03, Kahe07, KanB,         Molo07, Molo13, Olow03,         Papa, Puam03, Puam13


######################### (PART 2) ############################
###  Fiting the statistical models 

### A. Survival (logistic regression for size and harvest)

source("SurvivalAnalysis_IPM.R")

#Best model of survival is surv.U0.Y0:  size(y) + utrans(x) +(1|year); use fit on lnarea_t (rather than centered)
IPM.surv <-surv.A_U.Y0
s.y <- attr(IPM.surv,"fixef")[2]
s.x <- c(attr(IPM.surv,"fixef")[1],attr(IPM.surv,"fixef")[1]+attr(IPM.surv,"fixef")[3:(nx+3-2)])
s.x.vart <- as.numeric(attr(IPM.surv,"ST"))^2

### B. GROWTH
#Best model of growth is grow.U.Q0Y0.vP for MC and PL: y(t+1) ~ g.x*utrans + g.xy*y(t)*utrans + (1|uquad/uyear)
#however, quad-level variation is difficult to manage in the IPM, so refit without quad-level var
#IPM.grow: y(t+1)~g.x*utrans + g.y*y(t) + g.xy*utrans*y(t)+eps
#           g.x ~ N(g.x-hat, g.x.vart), eps~N(0,g.eps.var*(lnarea_t^(2*g.eps.vP)))
source("GrowthAnalysis_Linear_IPM.R")
IPM.grow <- grow.U.Y0.vP
g.x <- c(fixef(IPM.grow)[1], fixef(IPM.grow)[1]+fixef(IPM.grow)[3:(2+nx-1)])
g.xy <- c(fixef(IPM.grow)[2], fixef(IPM.grow)[2]+fixef(IPM.grow)[(3+nx-1):(2*nx)])
g.x.vart <- getVarCov(IPM.grow)[1]
g.eps.var <- IPM.grow$sigma^2
g.eps.vP <- coef(IPM.grow$modelStruct$varStruct,"power")



## C. RECRUITMENT - recruits from larval settlement
#We will use larval recruitment as an external driver 
#Best model for recruitment is rcrt.nbUD
# recruits ~ utrans+lnspcover+lnspcover^2; Negative Binomial
# E(recruits)=NB(utrans + lnspcover + lnspcover^2, theta)

##########for recruitment, some trouble for PL b/c no recruits at oe site -- singular matrix for zinb

source("Recruitment_IPM.R")
# IPM.rcrt <- rcrt.nbU_D2
# r.x <- c(coef(IPM.rcrt)[1], coef(IPM.rcrt)[1]+coef(IPM.rcrt)[2:nx])
# r.d <- coef(IPM.rcrt)[nx+1]
# r.d2 <- coef(IPM.rcrt)[nx+2]
# r.k <- IPM.rcrt$theta
#rcrt.nbU_D2 is the best fit, but for now use easier (non-DD) rcrt.nbU
IPM.rcrt <- rcrt.nbU
#IPM.rcrt<-rcrt.nb0
r.x <- c(coef(IPM.rcrt)[1], coef(IPM.rcrt)[1]+coef(IPM.rcrt)[2:nx])
#r.x <- rep(coef(IPM.rcrt)[1],nx)
r.k <- IPM.rcrt$theta

## C.3. Size distribution of recruits
#mean and variance of recruit size varies by site
# srcrt.glsU.vU:  rcrtsize~utrans where utrans ~ (utrans_hat,sigma^2_utrans)
#  varIdent structure fits each categry (here, each utrans) with a distinct sd
# MC: use varIdent; PC: no varIdent
IPM.rsize <- srcrt.glsU.vU
sr.x <- coef(IPM.rsize)
#sr.x.var <- (IPM.rsize$sigma*IPM.rsize$modelStruct$varStruct)^2   #vector for sites
sr.x.var <- (c(intervals(IPM.rsize)$sigma[2], intervals(IPM.rsize)$sigma[2]*intervals(IPM.rsize)$varStruct[nx:(nx+nx-2)]))^2
if (sp=="Porites compressa") {
  IPM.rsize <- srcrt.glsU
  sr.x <- coef(IPM.rsize)
  sr.x.var <- rep(IPM.rsize$sigma^2,nx)   #make one for each site b/c need it for p.vec later to be same across spp
}



## D. FISSION - Reproduction by fission
source("Fission_IPM.R")

## D.1. Probability of fission
#Fission depends on colony size and transect with some random variation by year
#  P(fiss) ~ Bernoulli(f) where f~N(lnarea_t + utrans + lnarea_t:utrans, sigx_t^2)
# use non-centered lnarea_t
IPM.fiss <- fiss.A_U.Y0
f.x <- c(attr(IPM.fiss,"fixef")[1], attr(IPM.fiss,"fixef")[1]+attr(IPM.fiss,"fixef")[3:(nx+1)])
f.y <- attr(IPM.fiss,"fixef")[2]
f.x.vart <- as.numeric(attr(IPM.fiss,"ST"))^2

## D.2. Number of Fission Products depends on size+location
# Number of Fissees ~ Zero-Trunc-Poisson(utrans + lnarea_prefiss)
IPM.nfiss <- nfiss.ztU_A
nf.x <- c(coef(IPM.nfiss)[1], coef(IPM.nfiss)[1]+coef(IPM.nfiss)[2:nx])
nf.y <- coef(IPM.nfiss)[nx+1]

## D.3. Size distribution of the fissees

IPM.fsize <- sfiss.glsU.vU
sf.x <- c(coef(IPM.fsize)[1], coef(IPM.fsize)[1]+coef(IPM.fsize)[2:nx])
#sf.x.var <- (IPM.fsize$sigma*c(1,IPM.fsize$modelStruct$varStruct))^2            #vector
sf.x.var <- (c(intervals(IPM.fsize)$sigma[2], intervals(IPM.fsize)$sigma[2]*intervals(IPM.fsize)$varStruct[nx:(nx+nx-2)]))^2


##HERE AS OF 25JUNE 4:20AM
########################## PART (3) ###########################
### Matrix of all the regression coefficients
ncoef<-21
p.vec<-matrix(0, ncoef, nx)  #creates a matrix of zeros with ncoef rows and nx (num of sites) columns

p.vec[1,]<-c(s.x)                #survival intercept by site
p.vec[2,]<-rep(s.y,nx)           #survival slope wrt size, same across sites
p.vec[3,]<-rep(s.x.vart,nx)      #temporal variance in surv intercept, same across sites
p.vec[4,]<-c(g.x)                #growth intercept by site
p.vec[5,]<-c(g.xy)               #growth slope by site
p.vec[6,]<-rep(g.x.vart,nx)      #temporal variance in growth intercepts
p.vec[7,]<-rep(g.eps.var,nx)     #additional variance in growth     
p.vec[8,]<-rep(g.eps.vP,nx)    #additional variance is power function of size
p.vec[9,]<-c(r.x)                #recruitment by site
#p.vec[10,]<-rep(r.d,nx)          #effect of density (spcover) on num of recruits
#p.vec[11,]<-rep(r.d2,nx)         #nonlinearity in effect of density (spcover^2) on num of recruits
p.vec[12,]<-rep(r.k,nx)          #clumpling parm for Neg Bin distribution of recruits
p.vec[13,]<-c(sr.x)              #size distribution of recruits varies by site
p.vec[14,]<-c(sr.x.var)          #var of recruit dist varies by site
p.vec[15,]<-c(f.x)               #prob of fission by site
p.vec[16,]<-rep(f.y,nx)          #prob of fission|size (const across sites)
p.vec[17,]<-rep(f.x.vart,nx)     #temporal variation in probability of fission
p.vec[18,]<-c(nf.x)              #number of fission products by site
p.vec[19,]<-rep(nf.y,nx)         #change in num fission products with size
p.vec[20,]<-c(sf.x)              #size dist of fission products - mean by site
p.vec[21,]<-c(sf.x.var)          #size dist of fission products - var by site

dimnames(p.vec)<-list(c("surv.int.bySite", "surv.slope", "surv.int.varT", "growth.int.bySite", "growth.slope.bySite", "growth.int.var", "growth.err.var","growth.err.varPower.wSize","recruit.int.bySite", "recruit.slopeDensity","recruit.slopeDensitySqr", "recruit.clumping","recruitSize.mean.bySite","recruitSize.var.bySite","Pfission.int.bySite","Pfission.slope","Pfission.int.varT","Nfission.int.bySite","Nfission.slope","Sfission.mean.bySite","Sfission.var.bySite"),unique(MC$utrans))

############################# PART (4)#########################
## FIGURES  -- MJD ignoring this section for now bc already have figures

# op<-par(mfrow=c(1,3))
# ## PLOT: Survival
# os<-order(s0) # ordering size 
# os.size<-(s0)[os] # reorganizing size based on that order
# os.surv<-with(dat, (lived07)[os]) # reoganizing surv based on order os
# 
# ## We chop off the last 605 obs of os.size and os.surv to remove the new born that were not included in the survival column.
# os.size<-os.size[!is.na(os.size)]
# os.surv<-os.surv[1:length(os.size)]
# 
# psz<-tapply(os.size,as.numeric(cut(os.size,100)),mean,na.rm=T)
# ps<-tapply((as.numeric(os.surv)-1),as.numeric(cut(os.size,100)),mean,na.rm=T)
# sizef<-seq(min(s0,na.rm=T),max(s0,na.rm=T),length=100)
# 
# plot(psz,ps,pch=19, main="Survival", xlab=expression("Diameter at t (log scale)"), ylab=expression("Probability of surviving"), cex.lab=1.5)
# 
# points(sizef, exp(s.itH + s.slH*sizef)/(1 + exp(s.itH +s.slH*sizef)),type="l", col="red", lty=3)
# 
# points(sizef,exp(s.itUH + s.slUH*sizef)/(1 + exp(s.itUH +s.slUH*sizef)),type="l", col="blue", lty=3)
# 
# 
# ## PLOT for Growth
# ## The issue here is how to fit the regression line when the plot is on log-log scale. I presented here two approaches (use log10 or ggplot2)
# 
# ## 1st method:log10
# ## quartz() ## use x11() for PC!!!
# plot(size.06, size.07, log="xy", main="Growth", col=as.numeric(treat), xlab="Diameter at t (mm, log scale)", ylab="Diameter at t+1 (mm, log scale)", cex.lab=1.5)
# abline(lm(log10(size.07)~log10(size.06)), lty=2, lwd=2, col="blue")
# 
# ## 2d method: ggplot2
# ## install.packages("ggplot2", dependencies=T)
# ## library(ggplot2) 
# ## qplot(size.06, size.07, log="xy") + geom_smooth(method = lm)
# 
# 
# ## PLOT Fecundity
# ## Reorganizing, re-ordering (mainly the fitted values) the data to get a graph that makes sense!!
# os20<-order(size0620)
# size0620n<-size0620[os20]
# fit20<-scape.nb2$fitted.values
# fit20<-fit20[os20]
# 
# ## There is the plot
# plot(size.06,esc.06, log="x", main="Fecundity", xlab="Diameter at t (mm, log scale)", ylab="Number of scapes at t", cex.lab=1.5)
# points(size0620n,fit20, type="l", lwd=3, lty=3, col="red");
# 
# par(op)


############################# PART (5) ##########################
## The IPM FUNCTIONS: s(x), g(y,x), f(y,x), c(y,x), p(y,x), K(y,x)
##For now, don't include temporal variation

### A. SURVIVAL function s(y) - y is size at t, y2 is size at t+1 
sy<-function(y, pvec, site) {
	ybeta<-pvec[1,site] + pvec[2,site]*y;
  s<-exp(ybeta)/(1+exp(ybeta))
	return(s);
}

### B. GROWTH function g(y2,y)
gy2y<-function(y2, y, pvec, site) {
	muy<-pvec[4, site] + pvec[5, site]*y
	sigmay2<-pvec[7, site]*abs(y)^(2*pvec[8, site])
	sigmay<-sqrt(sigmay2)
	g<-dnorm(y2, muy, sigmay)
	return(g)
}

### C.The SURVIVAL-GROWTH function P(y, x)
py2y<-function(y2,y, pvec, site) { 
	p<-sy(y, pvec, site)*gy2y(y2, y, pvec, site)
	return(p) 
	}

### D. FERTILITY (recruit) function f(x,y)

fy2y<-function(y2, y, pvec, site) {
	numrcrts<-exp(pvec[9,site])   #log-link function in NB fit
	scd<-dnorm(y2,pvec[13,site],pvec[14,site])   #size dist of recruits, mean and var by site
	f<-numrcrts*scd                              #note that there is no effect of x
	return(f)
}

### D. REPRO by FISSION (clonal)
cy2y<-function(y2, y, pvec, site) {site
	ybeta<-pvec[15, site] + pvec[16,site]*y;   #prob ratio of fission as function of site (int) and size(slope)
	p_clon<-exp(ybeta)/(1+exp(ybeta))      #prob of fission
	nclone<-pvec[18,site]+pvec[19,site]    #number of clones as a function of site and size (ignore temporal va r for now)
	scd_c<-dnorm(y2, pvec[20,site], pvec[21,site])  #prob of clone of size y2 given the size dist of fissees at each site
	c<-p_clon*nclone*scd_c
	return(c)
}

### The (master) KERNEL: K(y,x)= p(y,x) + f(y,x) + c(y,x)
Ky2y<-function(y2, y, pvec, site) {
	k<-py2y(y2, y, pvec, site)+fy2y(y2, y, pvec, site)+cy2y(y2, y, pvec, site)
	return(k) 
}

######################### PART (6) #############################
### NUMERICAL INTEGRATION 

## The big matrix function

bigmat<-function(bigM, pvec, site){
	## Set matrix size and convergence tolerance 
	min.sz<-0.9*min((range(MC$lnarea_t,na.rm=T)))
	max.sz<-1.1*max((range(MC$lnarea_t,na.rm=T)))

	# Compute meshpoints iteration matrix KD 
	h=(max.sz-min.sz)/(bigM+1); 
	y2=(seq(min.sz,max.sz,length=bigM)+seq(min.sz+h,max.sz+h,length=bigM))/2;  #get midpoints

	## Apply Kyx funct for y and y, pvec (=p.vec), treat (choose!!)
	K=outer(y2,y2,Ky2y, pvec, site);
	KD=h*K;
	return(KD);  		  ## This is your have the matrix
}

############################# PART (7) ##########################
## Eigenvalues and eigenvectors analysis, sensitivity,...

library(popbio)
library(Biodem)

## Checks: Ploting lambda for various big matrix size to see when it level off:

bigm<-c(100, 200, 300, 400, 500, 600, 700)
lam<-matrix(0, length(bigm), nx)


for (i in 1:length(bigm)){
	for (j in 1:nx){
	library(popbio)
	A<-array(0, c(bigm[i], bigm[i], length(bigm), nx))
	A[,,i,j]<-bigmat(bigm[i], pvec=p.vec, site=j)
	lam[i,j]<-eigen.analysis(A[,,i,j])$lambda1
	}
}	


## PLOT FOR THE BIG MATRIX
op<-par(mfrow=c(nx/3,nx/3))

for (j in 1:nx) {
  plot(bigm, lam[,j], main=unique(MC$utrans)[j], type="l", xlab="", ylab="", col=j, lwd=3)
}

par(oma=c(5, 4.5, 0, 0))
title(xlab="Size of the big matrix", ylab=expression(paste("Population growth rate, ", lambda)), outer=T, cex.lab=1.2)

par(op)

## The big matrix is 500 x 500, more than that, lambda does not change much

## Lambda, SSD, RV, Damping ratio, Elasticity matrix
big<-array(0, c(300, 300, nx))
eig<-vector("numeric", nx)

for (i in 1:nx){
	big[,,i]<-bigmat(300, pvec=p.vec, site=i)
  A<-eigen.analysis(big[,,i])
	eig[i]<-A$lambda1
  }

eigsensP<-array(0,c(nx,length(p.vec[,1])))
eigsensN<-array(0,c(nx,length(p.vec[,1])))
#sensN<-array(0, c(300, 300, nx,length(p.vec[,1])))
#sensP<-array(0, c(300, 300, nx,length(p.vec[,1])))

for (i in 1:nx){
  for (j in 1:length(p.vec[,1])) {
    svecP<-p.vec
    svecP[j,]<-1.1*svecP[j,]
    sensP<-bigmat(300, pvec=svecP, site=i)
    eigsensP[i,j]<-eigen.analysis(sensP)$lambda1
    svecN<-p.vec
    svecN[j,]<-0.9*svecN[j,]
    sensN<-bigmat(300, pvec=svecN, site=i)
    eigsensN[i,j]<-eigen.analysis(sensN)$lambda1
}
}

elasP <- (eig-eigsensP)/eig/(0.1)
par(mfrow=c(2,2))
par(mar = c(7, 2, 2, 2) + 0.1)
for (i in c(1,4,9,15)){
  barplot(elasP[,i],horiz=TRUE,xlab="Elasticity",col=cl[1:ns])
}
#legend("bottom",legend=sort(unique(utrans),decreasing=TRUE),horiz=TRUE,pch=1,col=cl[ns:1],cex=1,bty="n")

## To get lambda for each site
eig

par(mar = c(7, 4, 4, 2) + 0.1)
xvals<-barplot(eig)
text(xvals, par("usr")[3] - 0.25, srt = 45, adj = 1,labels = unique(MC$utrans), xpd = TRUE)

#################################

