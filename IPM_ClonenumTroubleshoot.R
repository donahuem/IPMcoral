###################################NumfissTroubleShooting###############################
###For plots
cl<-c("cyan","mediumturquoise","dodgerblue","green3","purple","olivedrab2","orchid1","blue","darkblue","darkorchid4")
x<-seq(from=-2.4,to=10.5,by=0.01)#Range of values for lines added to plots
library(VGAM) # for zero truncated poisson
library(vcd)# for mosaic plot
dsize<-data.frame(size=x)


##Standard distributions for count data were not adequately describing our observations See IPM_Plots.R
##What if we fit a poisson (not zero truncated) using all of the adult colonies
##the bypasses the probability of fission because only large fissers would make new colonies
##Code MC so that fissLG and M have numfiss and all colonies that are present at t and don't die (growth, fuseLG,fuseSM) are zero's
MC$numfissALL<-MC$numfiss
MC$numfissALL<-ifelse(MC$fateIPM=="growth"|MC$fateIPM=="fuseLG"|MC$fateIPM=="fuseSM",0,MC$numfissALL)
clonenum<-glm(numfiss~size*utrans,data=MC,family=poisson)
drop1(clonenum,test="Chisq")
#drop size:utrans p-value=0.1842
clonenum<-glm(numfiss~size+utrans,data=MC,family=poisson)
drop1(clonenum,test="Chisq")
#drop utrans p-value=0.3153
clonenum<-glm(numfiss~size,data=MC,family=poisson)
drop1(clonenum,test="Chisq")
#keep size p-value=1.083e-07
plot(MC$size,MC$numfiss,pch=18,xlab="size t",ylab="# of fission products per fission event",main="# of fission products for all adult colonies (not just fission)")
lines(x,predict(clonenum,newdata=dsize, "response"))
# Similar to zero truncated poisson without including transect

###Need to look into other ways of describing our numfiss relationship for Cmatrix

###does numfiss vary with prior fate?  Kind of...The higher numbers are in S,U,and M but this will track with larger size
test<-xtabs(~fate_prior+numfiss,MCf)
test<-as.data.frame(test)
mosaic(~fate_prior+numfiss,data=test,rot_labels=c(90,0,90,0))
###does numfiss vary with year?
test<-xtabs(~fyear+numfiss,MCf)
test<-as.data.frame(test)
mosaic(~fyear+numfiss,data=test,rot_labels=c(90,0,90,0))
##not really
###Here is what we have:
#fec1=1 for all fissLG and M 0 for growth, fuseLG,fuseSM,M nothing for dead, recruit or fissSM => same as Pfiss
#fec2=# of fission products for fissLG, fissSM, M
#fec3=sum of fission products for fissLg, fissSM,M
#fec4=# of fission products for fissLG and M = same as numfiss
#numfiss = # of fission products for fissLG and M only => not the same as fec2
#sumfiss = sum of fission products for fissLG, fissSM,M=>same as fec3

###Make some New Columns
MCf$lnsum<-log(MCf$sumfiss)# log transformed sum of fissSM per fission event. This is for fissLg,fissSM, and M (see sumfiss above)
MCf$pSMsum<-MCf$sumfiss/MCf$asize # proportion of colony area at t put into fission offspring
MCf$lnpSMsum<-log(MCf$pSMsum)
MCf$gensize<-MCf$asizeNext+MCf$sumfiss # genet size at t+1 = asize (largest ramet area)+sumfiss(sum of smaller ramet area)
MCf$lngensize<-log(MCf$gensize) # log transformed genet area at t+1
MCf$dsum<-MCf$asize-MCf$gensize #area lost from t to t+1.  cannot take log because negative values where genet larger at t+1 than at t
MCf$psum<-MCf$gensize/MCf$asize #proportional change in untransformed area  
MCf$lnpsum<-log(MCf$psum) # log transformation of proportional change in area
MCf$parealost<-MCf$dsum/MCf$asize # proportion of area lost (untransformed because some colonies grew despite fission)
MCf$drsum<-MCf$asize-MCf$asizeNext # area lost in largest ramet only
MCf$prsum<-MCf$asizeNext/MCf$asize # proportional change in area for largest ramet only
MCf$lnprsum<-log(MCf$prsum) # log transformed proportional change in area for largest ramet only

###Looking at genet level changes in area in relation to number of fission products###
plot(MCf$size,MCf$lngensize,col=cl[MCf$fnumfiss2],pch=18,xlab="size t",ylab="log(genet size at t+1)",main="Genet growth by numfiss")
legend("topleft", legend = sprintf(levels(MCf$fnumfiss2)), col = cl, pch= 18, xjust = 1, bg = "white")
lines(c(-2,10),c(-2,10),lty=2)
plot(MCf$size,MCf$lnpsum,col=cl[MCf$fnumfiss2],pch=18,xlab="size t",ylab="log of proportional growth = log (genet size at t+1/size t)",main="Proportional Genet growth by numfiss")
legend("topleft", legend = sprintf(levels(MCf$fnumfiss2)), col = cl, pch= 18, xjust = 1, bg = "white")
lines(c(-2,10),c(0,0),lty=2)
plot(MCf$size,MCf$parealost,col=cl[MCf$fnumfiss2],pch=18,xlab="size t",ylab="(size t-genet size at t+1)/size t",main="Proportion of tissue lost by numfiss")
legend("bottomright", legend = sprintf(levels(MCf$fnumfiss2)), col = cl, pch= 18, xjust = 1, bg = "white")
lines(c(-2,10),c(0,0),lty=2)
###higher numbers of fission products (the ones that we are not appropriately describing) seem to occur in larger colonies size t. And in colonies that shrink over time

plot(MCf$size,MCf$lnsum,col=cl[MCf$fnumfiss2],pch=18,xlab="size t",ylab="log(sum of small ramet area t+1)",main="genet size at t vs sum of small ramets at t+1")
legend("bottomright", legend = sprintf(levels(MCf$fnumfiss2)), col = cl, pch= 18, xjust = 1, bg = "white")
lines(c(-2,10),c(-2,10),lty=2)
##Sum of small ramet sizes increases with initial genet size
plot(MCf$size,MCf$pSMsum,col=cl[MCf$fnumfiss2],pch=18,xlab="size t",ylab="sum of small ramets t+1/size t",main="Proportion of area t in sum of small ramets at t+1")
legend("topright", legend = sprintf(levels(MCf$fnumfiss2)), col = cl, pch= 18, xjust = 1, bg = "white")
plot(MCf$size,MCf$lnpSMsum,col=cl[MCf$fnumfiss2],pch=18,xlab="size t",ylab="log(sum of small ramets t+1)/size t",main="Proportion of tissue into small fission products")
legend("bottomleft", legend = sprintf(levels(MCf$fnumfiss2)), col = cl, pch= 18, xjust = 1, bg = "white")
# no strong relationship between proportion of initial area contributed to small fission products and numfiss

###Some genets grow and others shrink.  What do functions describing growth only and shrinkage only look like?###
###Are there different processes ie. genet growth = reproduction but genet shrinkage = partial mortality?###
MCf$Net<-ifelse(MCf$lnpsum<0,"shrink","grow")
MCf$Net<-as.factor(MCf$Net)
netgrowth<-subset(MCf,MCf$lnpsum>0)
netshrink<-subset(MCf,MCf$lnpsum<0)
plot(MCf$size,MCf$numfiss,pch=18,col=c("blue","green")[MCf$Net],xlab="size t",ylab="# of fission products per fission event",main="# of fission products by size at t")
legend("topleft", legend =c("Net Growth","Net Shrinkage") ,lwd=1, col = c("blue","green"), xjust = 1, bg = "white")
clonenum<-vglm(numfiss~size,data=netgrowth,family=pospoisson)
lines(x,predict(clonenum,newdata=dsize, "response"), col="blue")
clonenum<-vglm(numfiss~size,data=netshrink,family=pospoisson)
lines(x,predict(clonenum,newdata=dsize, "response"), col="green")

###What would a function for only higher numbers of fission products look like?###
MCf$hilo<-ifelse(MCf$numfiss<4,"lo","hi")
MCf$hilo<-as.factor(MCf$hilo)
hinum<-subset(MCf,MCf$numfiss>3)
lonum<-subset(MCf,MCf$numfiss<4)
plot(MCf$size,MCf$numfiss,pch=18,col=c("purple","magenta")[MCf$hilo],xlab="size t",ylab="# of fission products per fission event",main="# of fission products by size at t: 2 data subsets")
legend("topleft", legend =c("hi","lo") ,lwd=1, col = c("purple","magenta"), xjust = 1, bg = "white")
clonenum<-vglm(numfiss~size,data=hinum,family=pospoisson)
lines(x,predict(clonenum,newdata=dsize, "response"), col="purple")
clonenum<-vglm(numfiss~size,data=lonum,family=pospoisson)
lines(x,predict(clonenum,newdata=dsize, "response"), col="magenta")

###What if we added net growth/shrinkage to the lo model (cannot add to hi as all but one have net shrinkage)?  
logrow<-subset(lonum,lonum$Net=="grow")
loshrink<-subset(lonum,lonum$Net=="shrink")
plot(MCf$size,MCf$numfiss,pch=18,col=c("purple","magenta")[MCf$Net],xlab="size t",ylab="# of fission products per fission event",main="# of fission products by size at t")
legend("topleft", legend =c("logrow","loshrink") ,lwd=1, col = c("purple","magenta"), xjust = 1, bg = "white")
clonenum<-vglm(numfiss~size,data=logrow,family=pospoisson)
lines(x,predict(clonenum,newdata=dsize, "response"), col="purple")
clonenum<-vglm(numfiss~size,data=loshrink,family=pospoisson)
lines(x,predict(clonenum,newdata=dsize, "response"), col="magenta")
###Including net in model does not improve AIC and does not look very different in plot.

###Seems kind of arbitrary to treat higher numbers of fission products as outliers with no data based explanation
###There is obviously a relationship between size and number of fission products###
###Why does this relationship exist?###
###makes sense that if you have more polyps you can break into more parts###
###Does proportional tissue loss or surface area explain some of the variation in number of fission products?###
plot(MCf$lnpsum,MCf$numfiss,pch=18)
netshrink$lndsum<-log(netshrink$dsum)
plot(netshrink$lndsum,netshrink$numfiss,pch=18)
plot(MCf$lnsum,MCf$numfiss,pch=18)
###Because we are looking at net changes we are not able to tell how much of the colony from time t died.   
### New growth is included in SA.
###You could have a colony that had a lot of partial mortality but also a lot of growth.

###What if we used "outliers" based on lower quantiles of the proportional change in surface area
par(mfrow=c(2,3))
for (i in c(0.05,0.1,0.15,0.2,0.25)){
dl<-qnorm(i,mean(MCf$lnpsum),sd(MCf$lnpsum))
MCf$ES<-ifelse(MCf$lnprsum<dl,"ES","NGS")
MCf$ES<-as.factor(MCf$ES)
outliers<-subset(MCf,MCf$ES=="ES")
notout<-subset(MCf,MCf$ES=="NGS")
plot(MCf$size,MCf$numfiss,col=c("red","green")[MCf$ES],pch=18,xlab="size t",ylab="# of fisison products per fission event",main=as.character(i))
clonenum<-vglm(numfiss~size,data=outliers,family=pospoisson)
lines(x,predict(clonenum,newdata=dsize, "response"), col="red")
clonenum<-vglm(numfiss~size,data=notout,family=pospoisson)
lines(x,predict(clonenum,newdata=dsize, "response"), col="green")
}
legend("topleft", legend =c("outliers","not outliers") ,lwd=1, col = c("red","green"), xjust = 1, bg = "white")

##What if we used "outliers"based on upper quantiles of size?
for (i in c(0.95,0.9,0.85,0.8,0.75,0.7)){
dh<-qnorm(i,mean(MCf$size),sd(MCf$size))
MCf$S<-ifelse(MCf$size>dh,"big","small")
MCf$S<-as.factor(MCf$S)
big<-subset(MCf,MCf$S=="big")
small<-subset(MCf,MCf$S=="small")
plot(MCf$size,MCf$numfiss,col=c("red","green")[MCf$S],pch=18,xlab="size t",ylab="# of fisison products per fission event",main=as.character(i))
clonenum<-vglm(numfiss~size,data=big,family=pospoisson)
lines(x,predict(clonenum,newdata=dsize, "response"), col="red")
clonenum<-vglm(numfiss~size,data=small,family=pospoisson)
lines(x,predict(clonenum,newdata=dsize, "response"), col="green")
}
legend("topleft", legend =c("big","small") ,lwd=1, col = c("red","green"), xjust = 1, bg = "white")



###The hi vs lo model seems to fit best visually
###What would we need to make this work in our kernel
###Probability of being a "normal" fission event with <4 offspring
MCf$Philo<-ifelse(MCf$hilo=="lo",1,0)
###cannot use size*utrans interaction low n
hilo.U <- glm(Philo~size+utrans,data=MCf,family=binomial)
drop1(hilo.U,test="Chisq")
#drop utrans p-value=0.1676
hilo <- glm(Philo~size,data=MCf,family=binomial)
drop1(hilo,test="Chisq")
#keep size p-value=0.0002727
###Probably better to fit single model with hilo as an explantory variable
clonenum<-vglm(numfiss~size*hilo,data=MCf,family=pospoisson)
###get error message: 
#Error in if ((temp <- sum(wz[, 1:M, drop = FALSE] < wzepsilon))) warning(paste(temp,  : argument is not interpretable as logical
#Found Ben Bolker on line says outliers may prevent vglm from runninghttp://web.archiveorange.com/archive/v/cu5nv1EBuHyDQyIMbqzT
#Options are 1) fit 2 zt poissons 2) fit one poisson with outliers 3) fit zt poisson minus numfiss=10 outliers
#1) See above
#2)
clonenum<-glm(numfiss~size*hilo*utrans,data=MCf,family=poisson)
drop1(clonenum,test="Chisq")
#drop size:hilo:utrans p-value=0.3844
clonenum<-glm(numfiss~size*hilo+size*utrans+utrans*hilo,data=MCf,family=poisson)
drop1(clonenum,test="Chisq")
#drop size: utrans p-value=0.9010
clonenum<-glm(numfiss~size*hilo+utrans*hilo,data=MCf,family=poisson)
drop1(clonenum,test="Chisq")
#drop size:hilo p-value=0.5450
clonenum<-glm(numfiss~size+utrans*hilo,data=MCf,family=poisson)
drop1(clonenum,test="Chisq")
#drop hilo:utrans p-value=0.5873
clonenum<-glm(numfiss~size+utrans+hilo,data=MCf,family=poisson)
drop1(clonenum,test="Chisq")
#drop utrans p-value=0.89385
clonenum<-glm(numfiss~size+hilo,data=MCf,family=poisson)
drop1(clonenum,test="Chisq")
#keep size p-value=0.03539 and hilo p-value=2e-16
plot(MCf$size,MCf$numfiss,pch=18,col=c("purple","magenta")[MCf$hilo],xlab="size t",ylab="# of fission products per fission event",main="# of fission products by size at t: glm w/Poisson")
legend("topleft", legend =c("hi","lo") ,lwd=1, col = c("purple","magenta"), xjust = 1, bg = "white")
lines(x,exp(coefficients(clonenum)[1]+coefficients(clonenum)[2]*x),col="purple")
lines(x,exp(coefficients(clonenum)[1]+coefficients(clonenum)[3]+coefficients(clonenum)[2]*x),col="magenta")
#3)
outliers<-subset(MCf,MCf$numfiss<10)
clonenum<-vglm(numfiss~size*hilo*utrans,data=outliers,family=pospoisson)
clonenum<-vglm(numfiss~size*hilo,data=outliers,family=pospoisson)
#Get error:
#Error in if ((temp <- sum(wz[, 1:M, drop = FALSE] < wzepsilon))) warning(paste(temp,  : argument is not interpretable as logical
#Even without outliers with any interactions higher than main effect of hilo
#test for effects of utrans in zero truncated poisson without hilo
clonenum.U<-vglm(numfiss~size*utrans,data=MCf,family=pospoisson)
clonenum.U0<-vglm(numfiss~size+utrans,data=MCf,family=pospoisson)
lrtest_vglm(clonenum.U,clonenum.U0)
AIC(clonenum.U)
AIC(clonenum.U0)
# utrans*size interactions significantly improves fit
#We already knew this and chose not to use it because the estimates were unrealistic
##back to hilo with zero truncated
clonenum<-vglm(numfiss~size+hilo,data=outliers,family=pospoisson)
plot(MCf$size,MCf$numfiss,pch=18,col=c("purple","magenta")[MCf$hilo],xlab="size t",ylab="# of fission products per fission event",main="# of fission products by size at t: vglm w/pospoisson")
legend("topleft", legend =c("hi","lo") ,lwd=1, col = c("purple","magenta"), xjust = 1, bg = "white")
lines(x,exp(coefficients(clonenum)[1]+coefficients(clonenum)[2]*x),col="purple")
lines(x,exp(coefficients(clonenum)[1]+coefficients(clonenum)[3]+coefficients(clonenum)[2]*x),col="magenta")
###See IPM_HILOClonal.R for full Cmatrix###



#########################################################################################################################
###What about non-markovian mortality?
###If all of the fissSM die the next year then maybe we don't include them
##think of it as an establishment probability
##Let's look at survival with fate_prior

test<-xtabs(~fate_prior+fateIPM,MC)
test<-as.data.frame(test)
mosaic(~fate_prior+fateIPM,data=test,rot_labels=c(90,0,90,0))
test<-xtabs(~fate_prior+Psurv,MC)
test<-as.data.frame(test)
mosaic(~fate_prior+Psurv,data=test,rot_labels=c(90,0,90,0))
##survival does appear to vary with fate_prior.
##Recruitment has higher mortality, fission is slightly higher, M is 0 and U is lower
##What if we break I up into fissLG and fissSM?
MC$fate_priorSM<-MC$fate_prior
MC$fate_priorSM<-ifelse(MC$allfates=="ID"|MC$allfates=="II"|MC$allfates=="IM"|MC$allfates=="IS"|MC$allfates=="IU","Is",as.character(MC$fate_priorSM))
MC$fate_priorSM<-as.factor(MC$fate_priorSM)
test<-xtabs(~fate_priorSM+fateIPM,MC)
test<-as.data.frame(test)
mosaic(~fate_priorSM+fateIPM,data=test,rot_labels=c(90,0,90,0))
test<-xtabs(~fate_priorSM+Psurv,MC)
test<-as.data.frame(test)
mosaic(~fate_priorSM+Psurv,data=test,rot_labels=c(90,0,90,0))
## fissSM and recruitment have higher mortality, fissLG and growth about the same, M is zero and U is low
##Now lets fit survival models including fate_prior to see if it improves the fit
##cannot run with 3-way interaction low n
surv.UP<-glm(Psurv~size*utrans+size*fate_prior,data=MC,family=binomial)
drop1(surv.UP,test="Chisq")
#drop size*utrans p-value=0.27065
surv.U0P<-glm(Psurv~utrans+size*fate_prior,data=MC,family=binomial)
drop1(surv.U0P,test="Chisq")
#keep utrans p-value=0.0001017, size*fate_prior p-value=0.0437441

##What would survival look like using fate_prior as an explanatory variable?
surv.reg=glm(Psurv~utrans+size*fate_prior,data=MC,family=binomial)
summary(surv.reg)
nx=10
nfp=5
survI.int=c(coefficients(surv.reg)[1], coefficients(surv.reg)[1]+coefficients(surv.reg)[2:(nx)])
survI.slope=coefficients(surv.reg)[nx+1]
survM.int=c(survI.int+coefficients(surv.reg)[nx+2])
survM.slope=survI.slope+coefficients(surv.reg)[nx+nfp+1]
survR.int=c(survI.int+coefficients(surv.reg)[nx+3])
survR.slope=survI.slope+coefficients(surv.reg)[nx+nfp+2]
survS.int=c(survI.int+coefficients(surv.reg)[nx+4])
survS.slope=survI.slope+coefficients(surv.reg)[nx+nfp+3]
survU.int=c(survI.int+coefficients(surv.reg)[nx+5])
survU.slope=survI.slope+coefficients(surv.reg)[nx+nfp+4]
par(mfrow=c(3,4))
for (i in 1:10){
  plot(MC$size,MC$Psurv,col=cl[MC$fate_prior],type="p",xlab="Size (t)",ylab="Survival Probability")
  lines(x,exp(survI.int[i]+survI.slope*x)/(exp(survI.int[i]+survI.slope*x)+1),col="red")
  lines(x,exp(survM.int[i]+survM.slope*x)/(exp(survM.int[i]+survM.slope*x)+1),col="blue")
  lines(x,exp(survR.int[i]+survR.slope*x)/(exp(survR.int[i]+survR.slope*x)+1),col="green")
  lines(x,exp(survS.int[i]+survS.slope*x)/(exp(survS.int[i]+survS.slope*x)+1),col="purple")
  lines(x,exp(survU.int[i]+survU.slope*x)/(exp(survU.int[i]+survU.slope*x)+1),col="magenta")
}
legend("right", legend =sprintf(levels(MCf$fate_prior)) ,lwd=1, col = c("red","blue","green","purple","magenta"), xjust = 1, bg = "white")

surv.UP<-glm(Psurv~size*utrans+size*fate_priorSM,data=MC,family=binomial)
drop1(surv.UP,test="Chisq")
#drop size*utrans p-value=0.26130
surv.U0P<-glm(Psurv~utrans+size*fate_priorSM,data=MC,family=binomial)
drop1(surv.U0P,test="Chisq")
#drop size*fate p-value=0.646648
surv.U0P0<-glm(Psurv~utrans+size+fate_priorSM,data=MC,family=binomial)
drop1(surv.U0P0,test="Chisq")       
#drop fate_priorSM p-value=0.1186
surv.U0<-glm(Psurv~utrans+size,data=MC,family=binomial)
drop1(surv.U0,test="Chisq") 
#Keep size p-value=2.2e-16 and utrans p-value=1.278e-06
       
       
###Try with fissLG and fissSM separate anyway just to see what happens
surv.reg=glm(Psurv~utrans+size*fate_priorSM,data=MC,family=binomial)
summary(surv.reg)
nx=10
nfp=6
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
survS.int=c(survI.int+coefficients(surv.reg)[nx+6])
survS.slope=survI.slope+coefficients(surv.reg)[nx+nfp+5]      
par(mfrow=c(3,4))
for (i in 1:10){
 plot(MC$size,MC$Psurv,col=cl[MC$fate_prior],type="p",xlab="Size (t)",ylab="Survival Probability")
 lines(x,exp(survI.int[i]+survI.slope*x)/(exp(survI.int[i]+survI.slope*x)+1),col="red")
 lines(x,exp(survIs.int[i]+survIs.slope*x)/(exp(survIs.int[i]+survIs.slope*x)+1),col="black")
 lines(x,exp(survM.int[i]+survM.slope*x)/(exp(survM.int[i]+survM.slope*x)+1),col="blue")
 lines(x,exp(survR.int[i]+survR.slope*x)/(exp(survR.int[i]+survR.slope*x)+1),col="green")
 lines(x,exp(survS.int[i]+survS.slope*x)/(exp(survS.int[i]+survS.slope*x)+1),col="purple")
 lines(x,exp(survU.int[i]+survU.slope*x)/(exp(survU.int[i]+survU.slope*x)+1),col="magenta")
 }
legend("right", legend =sprintf(levels(MC$fate_priorSM)) ,lwd=1, col = c("red","black","blue","green","purple","magenta"), xjust = 1, bg = "white")

###Need probability that a colony of size x in year t underwent fate_prior in year t-1
xtabs(~fate_priorSM,MC)
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

##Try glm with size, utrans, fyear (cannot use 3-way interaction or 3 2-way interactions)
I.UY<-glm(PI~size*utrans+size*fyear,data=MC,family=binomial)
drop1(I.UY,test="Chisq")
I.U0Y<-glm(PI~utrans+size*fyear,data=MC,family=binomial)
drop1(I.U0Y,test="Chisq")
I.U0Y0<-glm(PI~utrans+size+fyear,data=MC,family=binomial)
drop1(I.U0Y0,test="Chisq")
I.Y0<-glm(PI~size+fyear,data=MC,family=binomial)
drop1(I.Y0,test="Chisq")


##Looks like year matters so try year as RE
I.U.Y0 <- glmer(PI~size*utrans+(1|fyear),data=MC,family=binomial)
drop1(I.U.Y0,test="Chisq")
#drop size:utrans p-value=0.6226
I.U0.Y0 <- glmer(PI~size+utrans+(1|fyear),data=MC,family=binomial)
drop1(I.U0.Y0,test="Chisq")
#drop utrans p-value=0.234722
I.Y0 <- glmer(PI~size+(1|fyear),data=MC,family=binomial)
drop1(I.Y0,test="Chisq")
#kep size p-value=0.002615

##Compare RE and FE models
I<-glm(PI~size,data=MC,family=binomial)
AIC(I,I.Y0)
###Year as RE improves AIC by delta=3.9822

Is.UY<-glm(PIs~size*utrans+size*fyear+utrans*fyear,data=MC,family=binomial)
drop1(Is.UY,test="Chisq")
Is.UY<-glm(PIs~size*utrans+size*fyear,data=MC,family=binomial)
drop1(Is.UY,test="Chisq")
Is.U0Y<-glm(PIs~utrans+size*fyear,data=MC,family=binomial)
drop1(Is.U0Y,test="Chisq")

##looks like fyear matters try fyear as RE
Is.U.Y0 <- glmer(PIs~size*utrans+(1|fyear),data=MC,family=binomial)
drop1(Is.U.Y0,test="Chisq")
#drop size:utrans p-value=0.2321
Is.U0.Y0 <- glmer(PIs~size+utrans+(1|fyear),data=MC,family=binomial)
drop1(Is.U0.Y0,test="Chisq")
#keep size p-value=2.2e-16 and utrans p-value=2.2e-16
#compare FE and RE models
Is.U<-glm(PIs~size*utrans,data=MC,family=binomial)
drop1(Is.U,test="Chisq")
#drop size:utrans p-value-0.2222
Is.U0<-glm(PIs~size+utrans,data=MC,family=binomial)
drop1(Is.U0,test="Chisq")
#keep size p-value=2.2e-16 and utrans p-value=2.2e-16
AIC(Is.U0.Y0,Is.U0)
#fyear RE improves AIC by delta=22.397

###This was the fullest model that would converge
M.U0Y0<-glm(PM~size+utrans+fyear,data=MC,family=binomial)
drop1(M.U0Y0,test="Chisq")
#drop fyear p-value=0.40457
M.U0<-glm(PM~size+utrans,data=MC,family=binomial)
drop1(M.U0,test="Chisq")
#keep size p-value=1.513e-09 and utrans p-value=0.01866

##fullest model that will converge
R.UY<-glm(PR~size*utrans+size*fyear,data=MC,family=binomial)
drop1(R.UY,test="Chisq")
#drop size:fyear p-value=0.8080
R.UY0<-glm(PR~size*utrans+fyear,data=MC,family=binomial)
drop1(R.UY0,test="Chisq")
#drop size:utrans p-value=0.50506
R.U0Y0<-glm(PR~size+utrans+fyear,data=MC,family=binomial)
drop1(R.U0Y0,test="Chisq")
#keep size p-value=3.865e-12, utrans p-value= 0.0002154, fyear p-value= 0.0199185

##looks like year matters so try as RE
R.U.Y0 <- glmer(PR~size*utrans+(1|fyear),data=MC,family=binomial)
drop1(R.U.Y0,test="Chisq")
#drop size:utrans p-value=0.5254
R.U0.Y0 <- glmer(PR~size+utrans+(1|fyear),data=MC,family=binomial)
drop1(R.U0.Y0,test="Chisq")
#keep size p-value=1.566e-12 and utrans p-value=0.0004346

##Compare RE and FE models
R.U<-glm(PR~size*utrans,data=MC,family=binomial)
drop1(R.U,test="Chisq")
#drop size:utrans p-value=0.5798
R.U0<-glm(PR~size+utrans,data=MC,family=binomial)
drop1(R.U0,test="Chisq")
#keep size p-value=5.097e-13 and utrans p-value=0.0006322
AIC(R.U0,R.U0.Y0)
##fyear RE does not significantly improve AIC with delta=0.2487



S.UY<-glm(PS~size*utrans*fyear,data=MC,family=binomial)
drop1(S.UY,test="Chisq")
#drop size:utrans:fyear p-value= 0.3445
S.UY<-glm(PS~size*utrans+size*fyear+utrans*fyear,data=MC,family=binomial)
drop1(S.UY,test="Chisq")
#drop size:fyear p-vaule=0.9001
S.UY<-glm(PS~size*utrans+utrans*fyear,data=MC,family=binomial)
drop1(S.UY,test="Chisq")
#drop size:utrans p-value=0.4407
S.U0Y<-glm(PS~size+utrans*fyear,data=MC,family=binomial)
drop1(S.U0Y,test="Chisq")
#keep size p-value= 0.0006388 and utrans:fyear p-value=4.204e-06

##looks like fyear matters let's try with fyear RE
S.U.Y0 <- glmer(PS~size*utrans+(1|fyear),data=MC,family=binomial)
drop1(S.U.Y0,test="Chisq")
#drop size:utrans p=0.5244
S.U0.Y0 <- glmer(PS~size+utrans+(1|fyear),data=MC,family=binomial)
drop1(S.U0.Y0,test="Chisq")
#keep size p-value=0.0004909 and utrans p-value=2.032e-05

##compare RE and FE models
S.U<-glm(PS~size*utrans,data=MC, family=binomial)
drop1(S.U,test="Chisq")
#drop size:utrans p-value=0.5351
S.U0<-glm(PS~size+utrans,data=MC, family=binomial)
drop1(S.U0,test="Chisq")
#keep size p-value=0.00015 and utrans p-value=6.925e-06
AIC(S.U0.Y0,S.U0)
#fyear RE improves AIC by delta=5.359




U.UY<-glm(PU~size*utrans+size*fyear,data=MC,family=binomial)
drop1(U.UY,test="Chisq")
#drop size:fyear p-value=0.34336
U.UY0<-glm(PU~size*utrans+fyear,data=MC,family=binomial)
drop1(U.UY0,test="Chisq")
#drop fyear p-value=0.2423
U.U<-glm(PU~size*utrans,data=MC,family=binomial)
drop1(U.U,test="Chisq")
#keep size:utrans p-value=0.02551

## fyear does not appear to matter
U.U<-glm(PU~size*utrans,data=MC,family=binomial)
#keep size:utrans p-value=0.02551



###What if we lumped I and S, M and U the two groups are fairly similar in mosaic plot let recruits and small fission vary
index<-c("I","Is","M","R","S","U")
values<-c("S","Is","U","R","S","U")
MC$fate_priorL<-values[match(MC$fate_priorSM,index)] # To estimate probability of I fate_priorSM
MC$fate_priorL<-as.factor(MC$fate_priorL)
index<-c("S","Is","R","U")
values<-c(0,1,0,0)
MC$PLIs<-values[match(MC$fate_priorL,index)]
values<-c(0,0,1,0)
MC$PLR<-values[match(MC$fate_priorL,index)]
values<-c(0,0,0,1)
MC$PLU<-values[match(MC$fate_priorL,index)]

survL.UF<-glm(Psurv~size*fate_priorL+size*utrans,data=MC,family=binomial)
drop1(survL.UF,test="Chisq")
#drop size:utrans p-value=0.26664
survL.U0F<-glm(Psurv~utrans+size*fate_priorL,data=MC,family=binomial)
drop1(survL.U0F,test="Chisq")
#keep size:fate_priorL p-value=0.0242320 and utrans p-value=0.0001483

summary(survL.U0F)
nx=10
nfp=4
survLIs.int=c(coefficients(survL.U0F)[1], coefficients(survL.U0F)[1]+coefficients(survL.U0F)[2:nx])
survLIs.slope=coefficients(survL.U0F)[nx+1]
survLR.int=c(survLIs.int+coefficients(survL.U0F)[nx+2])
survLR.slope=c(survLIs.slope+coefficients(survL.U0F)[nx+nfp+1])
survLS.int=c(survLIs.int+coefficients(survL.U0F)[nx+3])
survLS.slope=survLIs.slope+coefficients(survL.U0F)[nx+nfp+2]
survLU.int=c(survLIs.int+coefficients(survL.U0F)[nx+4])
survLU.slope=survLIs.slope+coefficients(survL.U0F)[nx+nfp+3]


par(mfrow=c(3,4))
par(mfrow=c(3,4))
for (i in 1:10){
  plot(MC$size,MC$Psurv,col=cl[MC$fate_priorL],type="p",xlab="Size (t)",ylab="Survival Probability")
  lines(x,exp(survLIs.int[i]+survLIs.slope*x)/(exp(survLIs.int[i]+survLIs.slope*x)+1),col="red")
  lines(x,exp(survLR.int[i]+survLR.slope*x)/(exp(survLR.int[i]+survLR.slope*x)+1),col="black")
  lines(x,exp(survLS.int[i]+survLS.slope*x)/(exp(survLS.int[i]+survLS.slope*x)+1),col="blue")
  lines(x,exp(survLU.int[i]+survLU.slope*x)/(exp(survLU.int[i]+survLU.slope*x)+1),col="green")
  
}
legend("right", legend =sprintf(levels(MC$fate_priorL)) ,lwd=1, col = c("red","black","blue","green"), xjust = 1, bg = "white")


PLIs.UY<-glm(PLIs~size*utrans+size*fyear,data=MC,family=binomial)
drop1(PLIs.UY,test="Chisq")
#drop size:utrans p-value=0.108640
PLIs.U0Y<-glm(PLIs~utrans+size*fyear,data=MC,family=binomial)
drop1(PLIs.U0Y,test="Chisq")
#keep utrans p-value<2e-16 and size:fyear p-value=0.03045

#fyear seems to matter so try year as FE
PLIs.U.Y0 <- glmer(PLIs~size*utrans+(1|fyear),data=MC,family=binomial)
drop1(PLIs.U.Y0,test="Chisq")
#drop size:utrans p-value=0.231
PLIs.U0.Y0 <- glmer(PLIs~size+utrans+(1|fyear),data=MC,family=binomial)
drop1(PLIs.U0.Y0,test="Chisq")
#keep size p-value<2.2e-16 and utrans p-value<2.2e-16
PLIs.U0<-glm(PLIs~size+utrans,data=MC,family=binomial)
AIC(PLIs.U0,PLIs.U0.Y0)
##RE improves AIC by delta=22.397


PLR.UY<-glm(PLR~size*utrans+size*fyear,data=MC,family=binomial)
drop1(PLR.UY,test="Chisq")
#drop size:fyear p-value=0.8080
PLR.UY0<-glm(PLR~size*utrans+fyear,data=MC,family=binomial)
drop1(PLR.UY0,test="Chisq")
#drop size:utrans p-value=0.50506
PLR.U0Y0<-glm(PLR~size+utrans+fyear,data=MC,family=binomial)
drop1(PLR.U0Y0,test="Chisq")
#keep size p-value=3.865e-12, utrans p-value=0.0002154, fyear p-value=0.0199185

PLR.U.Y0 <- glmer(PLR~size*utrans+(1|fyear),data=MC,family=binomial)
drop1(PLR.U.Y0,test="Chisq")
#drop size:utrans p-value=0.5254
PLR.U0.Y0 <- glmer(PLR~size+utrans+(1|fyear),data=MC,family=binomial)
drop1(PLR.U0.Y0,test="Chisq")
#keep size p-value=1.566e-12, utrans p-value=0.0004346
PLR.U0<-glm(PLR~size+utrans,data=MC,family=binomial)
AIC(PLR.U0.Y0,PLR.U0)
#RE of year does not significantly improve AIC delta=-0.2487

PLU.UY<-glm(PLU~size*utrans+size*fyear,data=MC,family=binomial)
drop1(PLU.UY,test="Chisq")
#drop size:utrans p-value=0.0.06324
PLU.U0Y<-glm(PLU~utrans+size*fyear,data=MC,family=binomial)
drop1(PLU.U0Y,test="Chisq")
#drop size:fyear p-value=0.0131266 because it will not converge
PLU.U0Y0<-glm(PLU~size+utrans+fyear,data=MC,family=binomial)
drop1(PLU.U0Y0,test="Chisq")
#drop fyear p-value=0.4781988
PLU.U0<-glm(PLU~size+utrans,data=MC,family=binomial)
drop1(PLU.U0,test="Chisq")
#keep size p-value<2.2e-16, utrans p-value=0.0001723

###year does not appear to mattter

plot(MCf$size,MCf$PLIs,col="red")








       
plot(MC$size,MC$PLIs,col="red",pch=18)
points(MC$size,MC$PLR,col="blue",pch=18)
points(MC$size,MC$PLS,col="green",pch=18)
points(MC$size,MC$PLU,col="purple",pch=18)
MCf$matching<-paste(MCf$uyear,MCf$colony,sep=".")
MCfiss<-subset(MC,MC$fateIPM=="fissLG"|MC$fateIPM=="fissSM"|MC$fateIPM=="M")
MCfiss$colony1<-substr(MCfiss$colony,0,2)
MCfiss$matching<-paste(MCfiss$uyear,MCfiss$colony1,sep=".")
bleh<-subset(MCfiss,MCfiss$lnpgensize<0)
plot(bleh$size,bleh$fec4,col=cl[bleh$fateIPM],pch=18,xlab="# of fission products",ylab="log(genet size at t+1/size at t)",main="Proportional growth by fate")
legend("topleft", legend = sprintf(levels(MCfiss$fateIPM)), col = cl, pch= 18, xjust = 1, bg = "white")

clonenum<-vglm(fec4~size,data=bleh,family=pospoisson)
plot(bleh$size,bleh$fec4,col=cl[bleh$utrans],type="p",pch=18,xlab="Size (t)",ylab="# of fission products/'parent colony'",ylim=c(0,10),xlim=c(-2,10))

lines(x,predict(clonenum,newdata=dsize, "response"), col="red")
clonenum<-vglm(numfiss~size,data=MCf,family=pospoisson)
lines(x,predict(clonenum,newdata=dsize, "response"), col="red")
bleh<-subset(MCfiss,MCfiss$fec4<4)
clonenum<-vglm(numfiss~size,data=bleh,family=pospoisson)lines(x,predict(clonenum,newdata=dsize, "response"), col="red")
lines(x,predict(clonenum,newdata=dsize, "response"), col="red")
bleh<-subset(MCfiss,MCfiss$lnpgensize>0)
plot(bleh$size,bleh$fec4,col=cl[bleh$utrans],type="p",pch=18,xlab="Size (t)",ylab="# of fission products/'parent colony'",ylim=c(0,10),xlim=c(-2,10))
clonenum<-vglm(fec4~size,data=bleh,family=pospoisson)
lines(x,predict(clonenum,newdata=dsize, "response"), col="red")



index<-MCf$matching
values<-MCf$gensize
MCfiss$gensize<-values[match(MCfiss$matching,index)]
values<-MCf$lngensize
MCfiss$lngensize<-values[match(MCfiss$matching,index)]
values<-MCf$psum
MCfiss$pgensize<-values[match(MCfiss$matching,index)]
values<-MCf$lnpsum
MCfiss$lnpgensize<-values[match(MCfiss$matching,index)]
MCfiss<-droplevels(MCfiss)
plot(MC$size,MC$sizeNext,col=cl[MC$fateIPM],pch=18,xlab="size t",ylab="size t+1",main="Growth by Fate")
lines(c(0,10),c(0,10),lty=2)
legend("topleft", legend = sprintf(levels(MC$fateIPM)), col = cl, pch= 18, xjust = 1, bg = "white")
plot(MCfiss$size,MCfiss$sizeNext,col=cl[MCfiss$fateIPM],pch=18,xlab="size t",ylab="Size t+1",main="Fission Growth by Fate")
lines(c(0,10),c(0,10),lty=2)
legend("topleft", legend = sprintf(levels(MCfiss$fateIPM)), col = cl, pch= 18, xjust = 1, bg = "white")
plot(MCfiss$lngensize,MCfiss$sizeNext,col=cl[MCfiss$fateIPM],pch=18,xlab="Total genet size t+1",ylab="Ramet size t+1")
lines(c(0,10),c(0,10),lty=2)
legend("topleft", legend = sprintf(levels(MCfiss$fateIPM)), col = cl, pch= 18, xjust = 1, bg = "white")

plot(MCfiss$sizeNext,MCfiss$lnpgensize,col=cl[MCfiss$fateIPM],pch=18,xlab="size t+1",ylab="log(genet size at t+1/size at t)",main="Proportional growth by numfiss")
legend("topleft", legend = sprintf(levels(MCfiss$fateIPM)), col = cl, pch= 18, xjust = 1, bg = "white")
lines(c(-2,10),c(0,0),lty=2)
plot(MCfiss$sizeNext,MCfiss$fec2,col=cl[MCfiss$fateIPM],pch=18,xlab="size t+1",ylab="# of fission products",main="Proportional growth by utrans")
legend("topleft", legend = sprintf(levels(MCfiss$fateIPM)), col = cl, pch= 18, xjust = 1, bg = "white")
plot(MCfiss$fec2,MCfiss$lnpgensize,col=cl[MCfiss$fateIPM],pch=18,xlab="# of fission products",ylab="log(genet size at t+1/size at t)",main="Proportional growth by fate")
legend("topleft", legend = sprintf(levels(MC$fateIPM)), col = cl, pch= 18, xjust = 1, bg = "white")
lines(c(-2,10),c(0,0),lty=2)

plot(MC$lnpsum,MC$numfiss,col=cl[MC$utrans],pch=18)

MC$Iprior<-as.factor(MC$fate_prior)
MC$Iprior<-ifelse(MC$allfates=="ID"|MC$allfates=="IS"|MC$allfates=="II"|MC$allfates=="IU"|MC$allfates=="IM","Is",MC$Iprior)
test<-xtabs(~Iprior+fateIPM,MC)
mosaic(~Iprior+fateIPM,data=test,rot_labels=c(90,0,90,0))
library(vcd)



MC$fate_priorSM<-ifelse(MC$fateIPM=="fissSM","Is",MC$fate_prior)
MC$fate_priorSM<-as.factor(MC$fate_priorSM)
View(MC)
