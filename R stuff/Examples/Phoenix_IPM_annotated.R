##### Model for Phoenix IPM dissertation Ch.3

#### Model subunits from model selection/analysis performed previously
    # See ramet.pop.analysis.R for details


rm(list=ls(all=TRUE)) ## Clear everything from workspace

setwd("C:/Users/lisa/Documents/Dissertation research/Ch. 3 IPM/IPM files/annotated_2_13") #set working directory
getwd()

library(car)
library(lattice)
library(MASS)
library(msm)
library(popbio)
library(gridBase)

################# PART 1. #########################
#### Build matrix of all the regression coefficients
    ##Coefficients from previous model selection
  p.coef<-matrix(0, 57)
  
  #Survival coefs
      p.coef[1]<- 2.527236 #intercept for survival for ramets, F0
      p.coef[2]<- 2.414087 #effect of width for ramets
      p.coef[3]<- -2.832430e-01 #effect of density (number of ramets per genet, log-transformed, centered)
      p.coef[4]<- 3.951971e-01 #intercept for seedlings/singletons, F0
      p.coef[5]<- 6.576427 #effect of width for seedlings/singetons
      p.coef[6]<- 1.156523e+01 #linear effect of grazing (applied only >0.5 prop grazed)
      p.coef[7]<- -2.800432e+01 #quadratic effect of grazing (applied only >0.5 prop grazed)
      p.coef[8]<- 3.930509e-01 #effect of F12
      p.coef[9]<- 6.614860e-01 #effect of NF
  
  #Growth coefs
      p.coef[10]<- 1.613615e-01 #growth intercept for ramet, F0
      p.coef[11]<- 3.115233e-01 #linear effect of ramet-level grazing
      p.coef[12]<- -3.703961e-01 #quadratic effect of ramet-level grazing
      p.coef[13]<- 5.540772e-01 #linear effect of ramet-level harvest
      p.coef[14]<- -2.662937e-01 #quadratic effect of ramet-level harvest
      p.coef[15]<-  1.300230e-03 #effect of canopy openness (cenetered), current models run only at mean canopy openness
      p.coef[16]<-  7.766554e-01 #effect of width, F0
      p.coef[17]<-  -1.569418e-01 #effect of genet-level grazing
      p.coef[18]<- -8.131360e-02 #effect of genet-level harvest
      p.coef[19]<-  3.251101e-03 #grazing-canopy openness interaction (cenetered for canopy), current models run only at mean canopy openness
      p.coef[20]<-  -3.446216e-01 #ramet-level harvest-width interaction
      p.coef[21]<-  9.443623e-03 # intercept for seedling/singleton
      p.coef[22]<- 1.067207 #width effect for seedling/singleton
      p.coef[23]<-  1.805337e-02 #sigma^2 for growth
  
      p.coef[24]<- 9.452016e-02 #intercept for growth,  F12
      p.coef[25]<- 4.277902e-02 #intercept for growth, NF
  
      p.coef[26]<- -1.466495e-01 #effect of width, F12
      p.coef[27]<-  -1.790834e-01 #effect of width, NF
  
  #Flowering/fruiting coef
      p.coef[28]<- -8.023579e+00 #flowering intercept, F12
      p.coef[29]<- -1.043771e+01 #flowering intercept, NF
      p.coef[30]<- -3.313485 #effect of ramet-level harvest
      p.coef[31]<- 9.094276 #effect of width 
      p.coef[32]<- 1.536585 #mean number of flowering stalks 
      p.coef[33]<- 3.448276e-02 #seedlings per stalk, H.F12
      p.coef[34]<- 0.52941176 #seedlings per stalk, NH.F12
      p.coef[35]<- 1 #seedlings per stalk, H.NF
      p.coef[36]<- 4.285714e-01 #seedlings per stalk, NH.NF
      p.coef[37]<- 2.812500e-01 #seedlings per stalk overall 
  
  #Veg repro coefs
      p.coef[38]<- -4.657256 #intercept for probability of vegetative reproduction
      p.coef[39]<- 1.882245e+01 #linear effect of genet-level grazing
      p.coef[40]<-  -1.724860e+02 #quadratic effect of genet-level grazing
      p.coef[41]<- -2.907268e-02 #effect of canopy openness (cenetered), this model only at mean canopy openness
      p.coef[42]<- 1.253968 #mean number of sprouts
      p.coef[43]<- 2.232986e-01 #intercept of sprout size  
      p.coef[44]<- 2.023910e-01 #effect of mother ramet size on new sprout size
      p.coef[45]<- -7.315659e-02 #effect of density (ramets per genet, log transformed,  centered) on new sprout size
      p.coef[46]<- 5.665980e-03 #sigma^2 sprout size
      p.coef[47]<- 1.429351 #effect of mother ramet size on varience, exponential effect    
  
  #Predicting ramet-level grazing and harvest from genet-level intensities (from multinomial model)
      p.coef[48]<- -4.32 #intercept for ramet-level grazed proportion
      p.coef[49]<- 8.88 #effect genet-level grazing on ramet-level grazing
      p.coef[51]<- 0.551 #effect of genet-level harvest on ramet-level grazing
      p.coef[51]<- 1.61 #effect of ramet width on ramet-level proportion grazed
      p.coef[52]<- 0.00676 #interaction between ramet width and genet-level harvest on ramet-level proportion grazed
      p.coef[53]<-  -4.27 #intercept for ramet-level harvest proportion
      p.coef[54]<- -0.664 #effect genet-level grazing on ramet-level harvest
      p.coef[55]<- 4.99 #effect of genet-level harvest on ramet-level harvest
      p.coef[56]<- 2.07 #effect of ramet width on ramet-level proportion harvested
      p.coef[57]<- 2.62 #interaction between ramet width and genet-level harvest on ramet-level proportion harvested
  
  dimnames(p.coef)<-list(c("s.int.ram.F0", "s.width.ram", "s.dens", "s.int.sdlg.F0", 
                           "s.width.sdlg", "s.sl.Ggen", "s.qu.Ggen", "s.int.F12", "s.int.NF",
                           "g.int.F0", "g.sl.Gram", "g.qu.Gram", "g.sl.Hram", "g.qu.Hram",
                           "g.can", "g.width.F0", "g.Ggen","g.Hgen", 
                           "g.Gram_can.cen", "g.Hram_width", 
                           "g.int.sdlg", "g.width.sdlg", "g.sigma2", "g.F12", "g.NF",
                           "g.width.F12", "g.width.NF",
                           "fl.int.F12", "fl.int.NF", 
                           "fl.Hram", "fl.width",
                           "no.stalks.lamba", "sdlgs.stalk.H.F12",
                           "sdlgs.stalk.NH.F12","sdlgs.stalk.H.NF","sdlgs.stalk.NH.NF",
                           "sdlgs.stalk.all",
                           "spr.int", "spr.sl.Ggen", "spr.qu.Ggen", "spr.sl.can",
                           "no.sprts.lambda","ss.int", "ss.mother.size", "ss.ramets.cen",
                           "ss.sig2", "ss.varExp", "g.pre.int", "g.pre.Ggen",
                           "g.pre.Hgen", "g.pre.width","g.pre.Hwidth.int", 
                           "h.pre.int", "h.pre.Ggen",
                           "h.pre.Hgen", "h.pre.width","h.pre.Hwidth.int"
                           ))

################# PART 2. #########################
#### The IPM FUNCTIONS: s(x), g(y,x), f(y,x), c(y,x), p(y,x), K(y,x)

    ### A. SURVIVAL function s(x) for ramets and seedlings
    sx.ramet<-function(x, coef, h.gen, g.gen, fire, density) {
      if(fire=="F0"){ #effect of grazing only when proportion grazed per genet is greater than 0.5
          if(g.gen<0.5) {xbeta<-coef[1]+coef[2]*x+coef[3]*density} else
                         xbeta<-coef[1]+coef[2]*x+coef[3]*density+coef[6]*g.gen+coef[7]*(g.gen^2)}
      if(fire=="F12"){
          if(g.gen<0.5) {xbeta<-coef[1]+coef[8]+coef[2]*x+coef[3]*density} else
                         xbeta<-coef[1]+coef[8]+coef[2]*x+coef[3]*density+coef[6]*g.gen+coef[7]*(g.gen^2)}
      if(fire=="NF"){
          if(g.gen<0.5) {xbeta<-coef[1]+coef[9]+coef[2]*x+coef[3]*density} else
                         xbeta<-coef[1]+coef[9]+coef[2]*x+coef[3]*density+coef[6]*g.gen+coef[7]*(g.gen^2)}
      s<-exp(xbeta)/(1+exp(xbeta))
    	return(s);
    }

    sx.single<-function(x, coef, h.gen, g.gen, fire, density) {
      if(fire=="F0"){  
          if(g.gen<0.5) {xbeta<-coef[4]+coef[5]*x} else
                         xbeta<-coef[4]+coef[5]*x+coef[6]*g.gen+coef[7]*(g.gen^2)}
      if(fire=="F12"){  
          if(g.gen<0.5) {xbeta<-coef[4]+coef[8]+coef[5]*x} else
                         xbeta<-coef[4]+coef[8]+coef[5]*x+coef[6]*g.gen+coef[7]*(g.gen^2)}
      if(fire=="NF"){  
          if(g.gen<0.5) {xbeta<-coef[4]+coef[9]+coef[5]*x} else
                         xbeta<-coef[4]+coef[9]+coef[5]*x+coef[6]*g.gen+coef[7]*(g.gen^2)}
           s<-exp(xbeta)/(1+exp(xbeta))
    	return(s);
    }
    
    ### B. GROWTH function g(y,x) for ramets and seedlings  
    gyx.ramet<-function(y, x, coef, h.gen, g.gen, fire, density) {
        #calculate ramet-level grazing and harvest based on genet-level input specified
        logit.intact=exp(0)
        logit.grazed=exp(coef[48]+coef[49]*g.gen+coef[50]*h.gen+coef[51]*x+coef[52]*h.gen*x)
        logit.harvest=exp(coef[53]+coef[54]*g.gen+coef[55]*h.gen+coef[56]*x+coef[57]*h.gen*x)
        g.ram<-logit.grazed/(logit.intact+logit.grazed+logit.harvest)
        h.ram<-logit.harvest/(logit.intact+logit.grazed+logit.harvest)
        if(h.gen==0) {h.ram<-0}
        if(g.gen==0) {g.ram<-0}

      if(fire=="F0"){  
          mux<-coef[10]+coef[11]*g.ram+coef[12]*(g.ram^2)+coef[13]*h.ram+coef[14]*(h.ram^2)+
            coef[16]*x+coef[17]*g.gen+coef[18]*h.gen+coef[20]*h.ram*x;}
      if(fire=="F12"){  
          mux<-coef[10]+p.coef[24]+coef[11]*g.ram+coef[12]*(g.ram^2)+coef[13]*h.ram+coef[14]*(h.ram^2)+
            (coef[16]+p.coef[26])*x+coef[17]*g.gen+coef[18]*h.gen+coef[20]*h.ram*x;}
      if(fire=="NF"){  
          mux<-coef[10]+p.coef[25]+coef[11]*g.ram+coef[12]*(g.ram^2)+coef[13]*h.ram+coef[14]*(h.ram^2)+
            (coef[16]+p.coef[27])*x+coef[17]*g.gen+coef[18]*h.gen+coef[20]*h.ram*x;}
      
    	sigma<-sqrt(coef[23])
    	g<-dtnorm(y, mux, sigma, lower=min(y), upper=max(y)) 
    	return(g);
    }

    gyx.single<-function(y, x, coef, h.gen, g.gen, fire, density) {
      #calculate ramet-level grazing and harvest based on genet-level input specified
        logit.intact=exp(0)
        logit.grazed=exp(coef[48]+coef[49]*g.gen+coef[50]*h.gen+coef[51]*x+coef[52]*h.gen*x)
        logit.harvest=exp(coef[53]+coef[54]*g.gen+coef[55]*h.gen+coef[56]*x+coef[57]*h.gen*x)
        g.ram<-logit.grazed/(logit.intact+logit.grazed+logit.harvest)
        h.ram<-logit.harvest/(logit.intact+logit.grazed+logit.harvest)
        if(h.gen==0) {h.ram<-0}
        if(g.gen==0) {g.ram<-0}
        
      if(fire=="F0"){
          mux<-coef[21]+coef[11]*g.ram+coef[12]*(g.ram^2)+coef[13]*h.ram+coef[14]*(h.ram^2)+
                coef[22]*x+coef[17]*g.gen+coef[18]*h.gen+coef[20]*h.ram*x;}
      if(fire=="F12"){
          mux<-coef[21]+p.coef[24]+coef[11]*g.ram+coef[12]*(g.ram^2)+coef[13]*h.ram+coef[14]*(h.ram^2)+
                (coef[22]+p.coef[26])*x+coef[17]*g.gen+coef[18]*h.gen+coef[20]*h.ram*x;}
      if(fire=="NF"){
          mux<-coef[21]+p.coef[25]+coef[11]*g.ram+coef[12]*(g.ram^2)+coef[13]*h.ram+coef[14]*(h.ram^2)+
                (coef[22]+p.coef[27])*x+coef[17]*g.gen+coef[18]*h.gen+coef[20]*h.ram*x;}
    	sigma<-sqrt(coef[23])
    	g<-dtnorm(y, mux, sigma, lower=min(y), upper=max(y))
    	return(g);
    }


   ### C.The SURVIVAL-GROWTH function P(y, x) for ramets and seedlings
    pyx.ramet<-function(y,x, coef, h.gen, g.gen, fire, density) { 
    	p<-sx.ramet(x, coef,h.gen,  g.gen, fire, density)*gyx.ramet(y, x, coef, h.gen, g.gen, fire, density)
    	return(p) 
    	}
    
    pyx.single<-function(y,x, coef,h.gen, g.gen, fire, density) { 
      p<-sx.single(x, coef, h.gen,  g.gen, fire, density)*gyx.single(y, x, coef,h.gen, g.gen, fire, density)
    	return(p) 
    	}

    ### D. FERTILITY (SEEDS) function f(y,x) for ramets only
    
    fyx<-function(y, x, coef, h.gen, g.gen, fire, density) {
      
        #calculate ramet-level grazing and harvest based on genet-level input specified
        logit.intact=exp(0)
        logit.grazed=exp(coef[48]+coef[49]*g.gen+coef[50]*h.gen+coef[51]*x+coef[52]*h.gen*x)
        logit.harvest=exp(coef[53]+coef[54]*g.gen+coef[55]*h.gen+coef[56]*x+coef[57]*h.gen*x)
        g.ram<-logit.grazed/(logit.intact+logit.grazed+logit.harvest)
        h.ram<-logit.harvest/(logit.intact+logit.grazed+logit.harvest)
        if(h.gen==0) {h.ram<-0}
        if(g.gen==0) {g.ram<-0}
      
      #Flowering for F0 depends on if it was F12 or NF at start of year, so average of two for simplicity
      if(fire=="F0") {beta<-((coef[28]+coef[29])/2)++coef[31]*x+coef[30]*h.ram;
                     p.flower<-ifelse(x<0.4, 0, (exp(beta)/(1+exp(beta))));} else
      if(fire=="F12") {beta<-coef[28]+coef[31]*x+coef[30]*h.ram;
                     p.flower<-ifelse(x<0.4, 0, (exp(beta)/(1+exp(beta))));} else
      if(fire=="NF") {beta<-coef[29]+coef[31]*x+coef[30]*h.ram;
                      p.flower<-ifelse(x<0.4, 0, (exp(beta)/(1+exp(beta))));
                      }
      
      n.stalks<-coef[32]

      #Can change thresholds for high harvest, originally set to 0.3
      #Divide by 2 with assumption that half of new seedlings are female...
      
      #Choose between modeling mean or harvest/fire-specific rates of seedling production
      sdlgs.stalk <-coef[37]/2 #overall mean
    
      #if(fire=="F12" & h.gen>=0.3) sdlgs.stalk<-coef[33]/2 #H.F12
      #if(fire=="F12" & h.gen<0.3) sdlgs.stalk<-coef[34]/2 #NH.F12
      #if(fire=="NF" & h.gen>=0.3) sdlgs.stalk<-coef[35]/2 #H.NF
      #if(fire=="NF" & h.gen<0.3) sdlgs.stalk<-coef[36]/2 #NH.NF
      #if(fire=="F0") sdlgs.stalk<-0.01
      
      #can divide p.flower by 3 to account for overestimating flower production
      f<-dunif(y, min=min(y), max=ifelse(length(y[y<0.15])>1, max(y[y<0.15]), y[2]))*(p.flower)*n.stalks*sdlgs.stalk #distribute new sdlgs evenly across size classes 0.1 and smaller

      return(f)
    }
    

    ### D. CLONAL FERTILITY c(y,x) for ramets only
    cyx<-function(y, x, coef, h.gen, g.gen, fire, density) {
    	xbeta<-coef[38]+coef[39]*g.gen+coef[40]*(g.gen^2);
    	sprouted<-ifelse(g.gen>0.3, 0, exp(xbeta)/(1+exp(xbeta))) #no sprouting above 0.3 g.gen
    	num.sprouts<-coef[42]
      
      mux<-coef[43]+coef[44]*x+coef[45]*density
      sigmax2<-coef[46]*exp(2*coef[47]*x)
      scd.sprouts<-dtnorm(y, mux, sqrt(sigmax2), lower=0, upper=max(y))
    	c<-sx.ramet(x, coef, h.gen, g.gen, fire, density)*sprouted*num.sprouts*scd.sprouts
    	return(c)
    }

    ### Seedling kernal (only survival-growth)
    single.Kyx<-function(y, x, coef,h.gen, g.gen, fire, density) {
      single.k<-pyx.single(y, x, coef,  h.gen,  g.gen, fire, density)
      return(single.k)
    }

    ### Ramets kernal (survival-growth and clonal reproduction)
    ramet.Kyx<-function(y, x, coef,  h.gen, g.gen, fire, density) {
      ramet.k<-pyx.ramet(y, x, coef, h.gen, g.gen, fire, density)+cyx(y, x, coef, h.gen, g.gen, fire, density)
      return(ramet.k)
      }

    ### Sexual reproduction (fertility) kept separate

################# PART 3. #########################
#### Numerical integration

## The big matrix function

bigmat<-function(bigM, coef,h.gen, g.gen, fire, density){ #big M is matrix size, h is mesh size, y is vector of midpoints
      if(!(fire %in% c("F0", "F12", "NF"))) {stop("Error: invalid fire category (F0, F12 or NF requred)")} else
      if(((h.gen+g.gen)>1)) {stop("Error: invalid harvest-grazing combination (must sum to less than 1)")} else
      if((density>2) | (-2.7>density)) {stop("Error: density outside bounds, remember to log transform & center")} else
  ## Set matrix size and convergence tolerance
  min.sz<-0.95*0.1 ##can play with limits, 95-105 percent of min and max
	max.sz<-1.05*1.5

	# Compute meshpoints iteration matrix KD 
	h=(max.sz-min.sz)/(bigM+1); 
	y.ramets=(seq(min.sz-h,max.sz,length=bigM)+seq(min.sz,max.sz+h,length=bigM))/2;  #midpoints
  y.single<-y.ramets[y.ramets<0.5] #seedlings/singletons defined to be less than 0.5
  
	## Apply Kyx funct for y and y for seedlings and ramets separately
	K.single=outer(y.ramets,y.single, single.Kyx, coef, h.gen, g.gen, fire, density);
	KD.single=h*K.single; #matrix for singleton transitions
  
  K.ramet=outer(y.ramets,y.ramets, ramet.Kyx, coef,  h.gen, g.gen, fire, density);
	KD.ramet=h*K.ramet
  
  ## Calculate sexual reproduction/fertility separately
  F.sdlg=outer(y.single,y.ramets, fyx, coef, h.gen, g.gen, fire, density);
  FD.sdlg=h*F.sdlg
    FD.sdlg[,1:sum(y.ramets<0.4)]<-0  #no seedlings made by sprouts <0.4

  #Bolting -- put together seedlings, ramets and fertility    
  full.matrix<-cbind(rbind(KD.single[1:length(y.single),], #seedlings staying seedlings from KD.single
        matrix(0,length(y.single),length(y.single)),  #seedlings can't become ramets<0.5 (by definition), so assign 0
        KD.single[(length(y.single)+1):length(y.ramets),]), #seedlings becoming ramets from KD.single
    rbind(FD.sdlg, KD.ramet)) #ramets surviving/growing/reproducing by seed and making sprouts from FD.sdlg and KD.ramet
  return(full.matrix);  		  ## This is the full, bolted-together matrix
}

################# PART 4. #########################
#### Eigenvalues and eigenvectors analysis, sensitivity,...

      ## Checks: Calculate and plot lambda for various big matrix sizes to see when it level off:
      
      bigm<-seq(50, 1000, 50)
      lambda<-matrix(0, length(bigm),3)
      
      for (i in 1:length(bigm)){
        for(j in 1:3){
        fire.opts<-c("F0","F12", "NF")
      	lambda[i,j]<-eigen.analysis(bigmat(bigm[i], coef=p.coef, h.gen=0, g.gen=0, fire=fire.opts[j], density=0))$lambda1
      	}}
      
      plot(bigm, lambda[,1], main="Fire", type="l", xlab="Matrix divisions", ylab="Lambda", col="black", lwd=2)
      plot(bigm, lambda[,2], main="Recent fire", type="l", xlab="Matrix divisions", ylab="Lambda", col="black", lwd=2)
      plot(bigm, lambda[,3], main="No fire", type="l", xlab="Matrix divisions", ylab="Lambda", col="black", lwd=2)

    ## Above 500 divisions, estimated lambda does not change much


    ## Effects of density-dependence
        
          densities<-seq(1,80, 2) #max observed is 80
          dens.lambda<-matrix(0, length(densities), 3)
          
          for (i in 1:length(densities)){
            for(j in 1:3){
            fire.opts<-c("F0","F12", "NF")
            dens.calc<-log(densities[i])-log(14.56757) #rescale for model based on log-transformed, centered data
          	dens.lambda[i,j]<-eigen.analysis(bigmat(500, coef=p.coef, h.gen=0, g.gen=0, fire=fire.opts[j], density=dens.calc))$lambda1
          	}
          }
        
        plot(densities, dens.lambda[,3], main="Density dependence", type="l", xlab="Mean ramets/genet", ylab="Lambda", col="black", lwd=1, ylim=c(0.98,1.03))
            abline(h=1, lty=2, col="grey")
            abline(v=14.56, lty=3, col="grey")    
            lines(densities, dens.lambda[,2], main="", lwd=1)
            lines(densities, dens.lambda[,1], main="", lty=2)
            lines(densities, dens.lambda[,3], lwd=2)

    ## Calculate lambda over combinations of harvest, grazing and fire
        h.intensity<-seq(0,0.85, 0.05)
        g.intensity<-seq(0,0.55, 0.05)
        
        ##it takes too much memory to store an array of all matrices for all combinations, so just store lambda and ssd
        hg.lambda<-array(0, c(length(h.intensity), length(g.intensity), 3))
        dimnames(hg.lambda)<-list(h.intensity, g.intensity,c("F0","F12", "NF"))
        hg.ssd<-array(0, c(length(h.intensity), length(g.intensity), 3, 637))
        dimnames(hg.ssd)<-list(h.intensity, g.intensity, c("F0","F12", "NF"), c(seq(0.09352295,0.4976948, length=137),seq(0.09352295,1.576477,length=500)))
        
        for(k in 1:length(h.intensity)){
          for (i in 1:length(g.intensity)){
            for(j in 1:3){
            fire.opts<-c("F0","F12", "NF")
            if((h.intensity[k]+g.intensity[i])<=1) {
                  temp<-eigen.analysis(bigmat(500, coef=p.coef, h.gen=h.intensity[k], g.gen=g.intensity[i], fire=fire.opts[j], density=0))
                  hg.lambda[k,i,j]<-temp$lambda1
                  hg.ssd[k,i,j,]<-temp$stable.stage
                  } else hg.lambda[k,i,j]<-NA         
            }}}
    
        #Plot harvest and grazing combined effects
            #Contour plots
        
          panel.filledcontour <- function(x, y, z, subscripts, at, col.regions =
         gray(seq(28,6,-2)/30),  col = gray(seq(28,6,-2)/30), ...)
        {
          stopifnot(require("gridBase"))
          z <- matrix(z[subscripts],
                      nrow = length(unique(x[subscripts])),
                      ncol = length(unique(y[subscripts])))
          if (!is.double(z)) storage.mode(z) <- "double"
          opar <- par(no.readonly = TRUE)
          on.exit(par(opar))
          if (panel.number() > 1) par(new = TRUE)
          par(fig = gridFIG(), omi = c(0, 0, 0, 0), mai = c(0, 0, 0, 0))
          cpl <- current.panel.limits()
          plot.window(xlim = cpl$xlim, ylim = cpl$ylim,
                      log = "", xaxs = "i", yaxs = "i")
          # paint the color contour regions
          .Internal(filledcontour(as.double(do.breaks(cpl$xlim, nrow(z) - 1)),
                                  as.double(do.breaks(cpl$ylim, ncol(z) - 1)),
                                  z, levels = as.double(at), col = col))
          # add contour lines
          contour(as.double(do.breaks(cpl$xlim, nrow(z) - 1)),
                  as.double(do.breaks(cpl$ylim, ncol(z) - 1)),
                  z, levels = as.double(at), add=T,
                  col = "black", # color of the lines
                  drawlabels=T,  # add labels or not
                  labcex=0.8
                 )
        } 
            x11()
            plot.new()
            levelplot(hg.lambda[,,"NF"], panel=panel.filledcontour, plot.args=list(newpage=FALSE),
                      cuts=13, at=c(0.8,0.85,0.90,0.95, 0.975,0.99,1,1.01,1.02,1.03,1.04,1.05,1.06), 
                          colorkey=list(col=gray(seq(28,6,-2)/30)),
                        labels=list(labels=T, col="black", cex=0.8, label.style="align"),
                        xlab="Harvest intensity", ylab="Grazing intensity")
            x11()
            plot.new()
            levelplot(hg.lambda[,,"F12"], panel=panel.filledcontour, plot.args=list(newpage=FALSE),
                      cuts=13, at=c(0.8,0.85,0.90,0.95, 0.975,0.99,1,1.01,1.02,1.03,1.04,1.05,1.06), 
                          colorkey=list(col=gray(seq(26,4,-2)/30)),
                        labels=list(labels=T, col="black", cex=0.8, label.style="align"),
                        xlab="Harvest intensity", ylab="Grazing intensity")
            x11()
            plot.new()
            levelplot(hg.lambda[,,"F0"], panel=panel.filledcontour, plot.args=list(newpage=FALSE),
                      cuts=13, at=c(0.8,0.85,0.90,0.95, 0.975,0.99,1,1.01,1.02,1.03,1.04,1.05,1.06), 
                          colorkey=list(col=gray(seq(26,4,-2)/30)),
                        labels=list(labels=T, col="black", cex=0.8, label.style="align"),
                        xlab="Harvest intensity", ylab="Grazing intensity")
        
        #Plot effects of harvest
            plot(h.intensity, hg.lambda[,g.intensity==0,"NF"], main="", type="l", xlab="Proportion leaves harvested per genet", ylab="Lambda", col="black", lwd=1, ylim=c(0.98,1.06))
                abline(h=1, lty=2, col="grey")
                lines(h.intensity, hg.lambda[,g.intensity==0,"F12"], main="", lwd=1)
                lines(h.intensity, hg.lambda[,g.intensity==0,"F0"], lty=2)
                lines(h.intensity, hg.lambda[,g.intensity==0,"NF"], lwd=2)
        
        #Plot effects of grazing
            plot(g.intensity, hg.lambda[h.intensity==0,,"NF"], main="", type="l", xlab="Proportion leaves grazed per genet", ylab="Lambda", col="black", lwd=1, ylim=c(0.88,1.06))
                abline(h=1, lty=2, col="grey")
                lines(g.intensity, hg.lambda[h.intensity==0,,"F12"], main="",  lwd=1)
                lines(g.intensity, hg.lambda[h.intensity==0,,"F0"], lty=2)
                lines(g.intensity, hg.lambda[h.intensity==0,,"NF"], lwd=2)
        