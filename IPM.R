###########################################################MeganR IPM Code #############################################################################
##Read Data
rm(list=ls(all=TRUE))
setwd("C://Users/Megan/IPMcoral")
data<- read.delim("sortIPMpack.txt")
############################################################################################################################################
##Clean up data
#remove "M" = both fission and fusion in one time step
data<-subset(data,fate_t!="M")

#replace blank cells with NA
data$fate_prior[data$fate_prior==""]<-NA
data$offspringNext[data$offspringNext==""]<-NA
data$fuseLGname[data$fuseLGname==""]<-NA

#Add columns with log and square transformations
data$size<-log(data$asize)
data$sizeNext<-log(data$asizeNext)
data$size2<-(data$size^2)

#needed column called surv instead of survIPM for IPMpack and now everything uses surv
data$surv<-data$survIPM

# needed column with number of fission products for fissLG.  fec2 has number of fission products for both fissLG and fissSM
data$fec4<-ifelse(data$fec1==1,data$fec2,NA)

# used to subset by site
data$utrans<-factor(paste(data$sitename,data$trans,sep="."))

#Made columns to calculate probability of fission, fusion or growth for growth functions
index<-c("growth","dead","fissSM","fissLG","fuseSM","fuseLG","recruit")
values<-c(0,NA,NA,0,NA,1,NA)
data$fusprob<-values[match(data$fateIPM,index)]
values<-c(1,NA,NA,0,NA,0,NA)
data$growprob<-values[match(data$fateIPM,index)]

#subset by species and site with high enough n for analysis
MC<-data[data$species=="MC" & data$sitename!="KanB" & data$sitename!="KanP" & data$sitename!="Papa" & data$utrans!="Puam.13m",]

#subset of data with growth only for growth function
MCg<-subset(MC,fateIPM=="growth"|fateIPM=="fuseLG"|fateIPM=="fissLG")

#subset of data for plotting
#renamed as factors to remove unused levels
MCg$fateIPM<-factor(MCg$fateIPM)
MCg$sitename<-factor(MCg$sitename)
#subset MC by growth fates
MCgrowth<-subset(MC,fateIPM=="growth")
MCfuse<-subset(MC,fateIPM=="fuseLG")
MCfiss<-subset(MC,fateIPM=="fissLG")

#Range of values for lines added to plots
x<-seq(from=-2,to=10,by=0.01)
x0<-data.frame(size=x,size2=x*x)
################################################################Exploratory Plots###################################################################################
#Visualize data
par(mfrow=c(2,2),mar=c(4,4,2,1))
plot(MC$size,jitter(MC$surv),main='Survival') # jittered to see easier
plot(MC$size,MC$sizeNext,main='Growth/Shrinkage/Stasis')  
plot(MC$size,MC$fec1,main='Fission Events') # jittered to see easier
hist(MC$sizeNext[is.na(MC$size)],main='Size of Recruits')

#function for comparison plots like IPMpack (a silly function that only works for me ;)
complot<- function(y,data,family,type,colors="black"){
  lm1<-glm(y~1,data=data,family=family)
  lm2<-glm(y~size,data=data,family=family)
  lm3<-glm(y~size+size2,data=data,family=family)
  
  if (type=="growth"){
    plot(data$size,data$sizeNext,pch=19,cex=1,ylim=c(-2,10),main=paste("Model Comparison",type),xlab="size",ylab="sizeNext",col=colors)
    lines(x,predict(lm1,newdata=x0,type="response"),col="2")
    lines(x,predict(lm2,newdata=x0,type="response"),col="3")
    lines(x,predict(lm3,newdata=x0,type="response"),col="4")
    legend("topleft", legend = sprintf("%s: %s = %.1f",c("1","size","size+size2"), c("AIC"), c(AIC(lm1),AIC(lm2),AIC(lm3))), col = c(2:4), lty = 1, xjust = 1, bg = "white")
  }
  else{
    if (type=="clone prob"|type=="rec prob"|type=="survival"){
      os <- order(data$size)  
      osClon <- (y)[os] 
      osSize<-(data$size)[os] 
      binnedSize <- tapply(osSize, as.numeric(cut(osSize, breaks=20)), mean, na.rm = TRUE); # bin Size data
      binnedClon <- tapply(osClon, as.numeric(cut(osSize, breaks=20)), mean, na.rm = TRUE) #bin Survival probabilities 
      plot(binnedSize, binnedClon, pch = 19, cex=2, ylim=c(0,1.5),xlab = "Size at t", ylab = "Producing clones at t+1", main = paste("Model Comparison",type))
      y0<-predict(lm1,newdata=x0);y0<-exp(y0)/(exp(y0)+1);lines(x,y0,col="2")
      y0<-predict(lm2,newdata=x0);y0<-exp(y0)/(exp(y0)+1);lines(x,y0,col="3")
      y0<-predict(lm3,newdata=x0);y0<-exp(y0)/(exp(y0)+1);lines(x,y0,col="4")
      legend("topleft", legend = sprintf("%s: %s = %.1f",c("1","size","size+size2"), c("AIC"), c(AIC(lm1),AIC(lm2),AIC(lm3))), col = c(2:4), lty = 1, xjust = 1, bg = "white")
    }
    else{
      os <- order(data$size)  
      osClon <- (y)[os] 
      osSize<-(data$size)[os] 
      binnedSize <- tapply(osSize, as.numeric(cut(osSize, breaks=20)), mean, na.rm = TRUE); # bin Size data
      binnedClon <- tapply(osClon, as.numeric(cut(osSize, breaks=20)), mean, na.rm = TRUE) #bin Survival probabilities
      plot(binnedSize, binnedClon, pch = 19, cex=2, xlab = "Size at t", ylab = "number of clones produced", main = paste("Model comparision",type))
      y0<- lines(x,predict(lm1,newdata=x0,type='response'), col="2")
      y0<- lines(x,predict(lm2,newdata=x0,type='response'), col="3")
      y0<-lines(x,predict(lm3,newdata=x0,type='response'), col="4")
      legend("topleft", legend = sprintf("%s: %s = %.1f",c("1","size","size+size2"), c("AIC"), c(AIC(lm1),AIC(lm2),AIC(lm3))), col = c(2:4), lty = 1, xjust = 1, bg = "white")

    }
  }
}

#Comparision plots of vital rate functions: All data
complot(MC$surv,MC,binomial,"survival")
complot(MCg$sizeNext,MCg,gaussian,"growth")
complot(MC$fec1,MC,binomial,"clone prob")
complot(MC$fec4,MC,poisson,"clone num")
complot(MC$fusprob,MC,binomial,"clone prob")

#Comparison plots of growth functions to see if fission, fusion, growth should have thier own parameters
par(mfrow=c(2,2))
complot(MCg$sizeNext,MCg,gaussian,"growth",c("red","green","blue")[MCg$fateIPM])
legend("bottomright", legend = sprintf(c("growth","fission","fusion")), col = c("blue","red","green"), pch = 19, xjust = 1, bg = "white")
index<-list(MCgrowth,MCfiss,MCfuse)
for (i in 1:length(index)){
  complot(index[[i]]$sizeNext,index[[i]],gaussian,"growth")
}
par(mfrow=c(1,1))
#all fission, fusion, growth functions on a single plot
plot(MCg$size,MCg$sizeNext,col=c("red","green","blue")[MCg$fateIPM],main="Growth functions by fate",xlab="size",ylab="sizeNext")
legend("bottomright", legend = sprintf(c("growth","fission","fusion")), col = c("blue","red","green"), pch = 1, xjust = 1, bg = "white")
Ilm<-lm(sizeNext~size,data=MCfiss)
Ulm<-lm(sizeNext~size,data=MCfuse)
Glm<-lm(sizeNext~size,data=MCgrowth)
y0<- lines(x,predict(Ilm,newdata=x0,type='response'), col="red")
y0<- lines(x,predict(Ulm,newdata=x0,type='response'), col="green")
y0<- lines(x,predict(Glm,newdata=x0,type='response'), col="blue")
legend("topleft", legend = c("growth","fission","fusion"), col = c("blue","red","green"), lty = 1, xjust = 1, bg = "white")
################################################################Parameters#########################################################################

#make slots for parameters
params<-data.frame(surv.int=NA,surv.slope=NA,growth.int=NA,growth.slope=NA,growth.sd=NA,recruit.int=NA,recruit.size.mean=NA,recruit.size.sd=NA,clone.int=NA,clone.slope=NA,clonenum.int=NA,clonenum.slope=NA,clonesize.int=NA,clonesize.slope=NA,clonesize.sd=NA,fuse.int=NA,fuse.slope=NA,grow.int=NA,grow.slope=NA)

#survival regression
surv.reg=glm(surv~size,data=MC,family=binomial())
summary(surv.reg)
params$surv.int=coefficients(surv.reg)[1]
params$surv.slope=coefficients(surv.reg)[2]

# growth regression (All fates together)
growth.reg=lm(sizeNext~size,data=MCg)
summary(growth.reg)
params$growth.int=coefficients(growth.reg)[1]
params$growth.slope=coefficients(growth.reg)[2]
params$growth.sd=summary(growth.reg)$sigma

#probability of growth (for fission, fusion, growth with diff growth functions)
grow=glm(growprob~size,data=MC,family=binomial())
summary(grow)
params$grow.int=coefficients(grow)[1]
params$grow.slope=coefficients(grow)[2]

#probability of cloning (for fission growth function and clonal matrix)
clone=glm(fec1~size,data=MC,family=binomial())
summary(clone)
params$clone.int=coefficients(clone)[1]
params$clone.slope=coefficients(clone)[2]

#number of clones per fission event 
clonenum=glm(fec4~size,data=MC,family=poisson())
summary(clonenum)
params$clonenum.int=coefficients(clonenum)[1]
params$clonenum.slope=coefficients(clonenum)[2]

# size distribution of clones
clonesize=lm(sizeNext~size,data=subset(MC,MC$fateIPM=="fissSM"))
summary(clonesize)
params$clonesize.int=coefficients(clonesize)[1]
params$clonesize.slope=coefficients(clonesize)[2]
params$clonesize.sd=summary(clonesize)$sigma

#recruitment rate
#replace recruit.int with number of recruits per "adult" averaged quad over years then total of quads
params$recruit.int=0.0905897

#size distribution of recruits
params$recruit.size.mean=mean(MC$sizeNext[MC$fateIPM=="recruit"])
params$recruit.size.sd=sd(MC$sizeNext[MC$fateIPM=="recruit"])

# probability of fusion (for fusion growth function)
fuse=glm(fusprob~size,data=MC,family=binomial())
params$fuse.int=coefficients(fuse)[1]
params$fuse.slope=coefficients(fuse)[2]

#plot the models over the data - figure 2
par(mfrow=c(2,3),mar=c(4,4,2,1))
plot(MCg$size,MCg$sizeNext,main='Growth/Shrinkage/Stasis')  
lines(x,predict(growth.reg,newdata=x0),col='red',lwd=3)
plot(MC$size,jitter(MC$surv),main='Survival') # jittered to see easier
lines(x,predict(surv.reg,newdata=x0,type='response'), col='red',lwd=3)
plot(MC$size,jitter(MC$fec1),main='Clones') # jittered to see easier
lines(x,predict(clone,newdata=x0,type='response'), col='red',lwd=3)
plot(MC$size,MC$fec2,main='Number of clones') # jittered to see easier
lines(x,predict(clonenum,newdata=x0,type='response'), col='red',lwd=3)
hist(MC$sizeNext[MC$fateIPM=="fissSM"],main='Size of Clones',freq=FALSE)
lines(x,dnorm(x,mean(params$clonesize.int+params$clonesize.slope),params$clonesize.sd), col='red',lwd=3)


## utility functions

#probability of surviving
s.x=function(x,params) {
  u=exp(params$surv.int+params$surv.slope*x)
  return(u/(1+u))
}
#growth function (All growth fates together)
g.yx=function(xp,x,params) {     	
  dnorm(xp,mean=params$growth.int+params$growth.slope*x,sd=params$growth.sd)
}
#Probability of growing (not fission or fusing)
p.grow.x=function(x,params) {
  u=exp(params$grow.int+params$grow.slope*x)
  return(u/(1+u))
}
# probability of fissing     
p.clone.x=function(x,params) {
  u=exp(params$clone.int+params$clone.slope*x)
  return(u/(1+u))
}
#clonal offspring 
c.yx=function(xp,x,params) {
  p.clone.x(x,params)*
    dnorm(xp,mean=params$clonesize.int+params$clonesize.slope*x,sd=params$clonesize.sd)*
    exp(params$clonenum.int+params$clonenum.slope*x)
}
# number of clones per adult for "constant correction" of clonal matrix (no size distribution of offspring)
c.y=function(x,params) {
  u=p.clone.x(x,params)*
    exp(params$clonenum.int+params$clonenum.slope*x)
  return (u)
}

#recruitment = sexual offspring 
f.yx=function(xp,x,params) {
  u=dnorm(xp,mean=params$recruit.size.mean,sd=params$recruit.size.sd)*
    exp(params$recruit.int)
  return (u)
}
# recruits per adult for "constant correction" of fecundity matrix (no size distribution of offspring)
f.y=function(xp,params) {
  u=exp(params$recruit.int)
  return (u)
}

# probability of fusion for fusion growth function      
p.fuse.x=function(x,params) {
  u=exp(params$fuse.int+params$fuse.slope*x)
  return(u/(1+u))
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
G=h*outer(y,y,g.yx,params=params) 
S=s.x(y,params=params) # survival 
I=p.clone.x(y,params=params)
Pi=G
for(i in 1:n) Pi[,i]=G[,i]*S[i]*I[i]
U=p.fuse.x(y,params=params)
Pu=G
for(i in 1:n) Pu[,i]=G[,i]*S[i]*U[i]
Gr=p.grow.x(y,params=params)
Pg=G
for (i in 1:n) Pg[,i]=G[,i]*S[i]*Gr[i]
P=G 
for(i in 1:n) P[,i]=G[,i]*S[i] 

#C is the clonal offspring kernel
C=h*outer(y,y,c.yx,params=params)
#F is the fecundity/sexually reproduced offspring kernel
F=h*outer(y,y,f.yx,params=params)

###############################################################IPM!###############################################################
#Full kernel with a single growth function (fission, fusion, growth all together)
K=P+F+C

#Full kernel with 3 growth functions.  One each for fission, fusion and growth
Kiu=Pg+Pi+Pu+F+C

#get lamda,v,w  
lam=as.real(eigen(K)$values[1]) 
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
 