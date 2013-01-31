######Basic IPM with Cory Code Olowalu, MC, all years, cheat Fecundity, colonal Matrix######
####The data read and organization is the same in BasicIPMpack code######
rm(list=ls(all=TRUE))
setwd("C://Users/Megan/IPMcoral")
data<- read.delim("sortIPMpack.txt")
######Clean up data: remove M, add NA, make recruit fecs,add stage and stageNext######
data<-subset(data,fate_t!="M")
data$fate_prior[data$fate_prior==""]<-NA
data$offspringNext[data$offspringNext==""]<-NA
data$fuseLGname[data$fuseLGname==""]<-NA
data$size<-log(data$asize)
data$sizeNext<-log(data$asizeNext)
data$surv<-data$survIPM
index<-c("growth","dead","fissSM","fissLG","fuseSM","fuseLG","recruit")
values<-c(1,NA,NA,1,NA,1,NA)
data$rec<-values[match(data$fateIPM,index)]
values2<-c("continuous","continuous","continuous","continuous","continuous","dead","seedbank")
data$stage<-as.factor(values2[match(data$fateIPM,index)])
values3<-c("continuous","dead","continuous","continuous","dead","continuous","continuous")
data$stageNext<-as.factor(values3[match(data$fateIPM,index)])
data$size2<-(data$size^2)
data$fec4<-ifelse(data$fec1==1,data$fec2,NA)
######Use Olowalu MC because best n######
OMC<-subset(data,sitename=="Olow"&species=="MC")
######growth,fuseLg, fissLg only for growth makeGrowthObj only removes values where size or sizeNext are NA so keeps fuseSM and fissSM. We don't want those for growth######
OMCg<-subset(OMC,fateIPM=="growth"|fateIPM=="fuseLG"|fateIPM=="fissLG")
#####################################################################################################
# make some plots -  figure 1
par(mfrow=c(2,2),mar=c(4,4,2,1))
plot(OMC$size,jitter(OMC$surv),main='Survival') # jittered to see easier
plot(OMC$size,OMC$sizeNext,main='Growth/Shrinkage/Stasis')  
plot(OMC$size,OMC$fec1,main='Fission Events') # jittered to see easier
hist(OMC$sizeNext[is.na(OMC$size)],main='Size of Recruits')

x<-seq(from=-2,to=10,length=1001)
x0<-data.frame(size=x,size2=x*x)
xx=seq(-2,10,by=.01)
####################survival model selection modified from Rob's complex life cycle exercise###############
par(mfrow=c(1,2),mar=c(4,4,2,1))
so1<-glm(surv~1,data=OMC,family=binomial())
so2<-glm(surv~size,data=OMC,family=binomial())
so3<-glm(surv~size+(size2),data=OMC,family=binomial())
fs <- order(OMC$size)
fs.fec <- (OMC$surv)[fs]
fs.size <- (OMC$size)[fs]
pfz <- tapply(fs.size, as.numeric(cut(fs.size, 21)), mean, na.rm = TRUE) # 6=ncuts
ps <- tapply(fs.fec, as.numeric(cut(fs.size, 21)), mean, na.rm = TRUE)
plot(as.numeric(pfz), as.numeric(ps), pch = 19,cex=2,ylim=c(0,1),xlab="size",ylab="proportion surviving",main="Survival")
y0<-predict(so1,newdata=x0);y0<-exp(y0)/(exp(y0)+1);lines(x,y0,col="red")
y0<-predict(so2,newdata=x0);y0<-exp(y0)/(exp(y0)+1);lines(x,y0,col="green")
y0<-predict(so3,newdata=x0);y0<-exp(y0)/(exp(y0)+1);lines(x,y0,col="blue")
c(AIC(so1),AIC(so2),AIC(so3))
legend("topleft", legend = sprintf("%s: %s = %.1f",c("1","size","size+size2"), c("AIC"), c(AIC(so1),AIC(so2),AIC(so3))), col = c(2:4), lty = 1, xjust = 1, bg = "white")
########make slots for parameters######
params<-data.frame(surv.int=NA,surv.slope=NA,growth.int=NA,growth.slope=NA,growth.sd=NA,recruit.int=NA,recruit.size.mean=NA,recruit.size.sd=NA,clone.int=NA,clone.slope=NA,clonenum.int=NA,clonenum.slope=NA,clonesize.int=NA,clonesize.slope=NA,clonesize.sd=NA)
# 1. survival regression from Cory's code
surv.reg=glm(surv~size,data=OMC,family=binomial())
summary(surv.reg)
params$surv.int=coefficients(surv.reg)[1]
params$surv.slope=coefficients(surv.reg)[2]

####################growth model selection modified from Rob's complex life cycle exercise###############
go1<-lm(sizeNext~1,data=OMCg)
go2<-lm(sizeNext~size,data=OMCg)
go3<-lm(sizeNext~size+(size2),data=OMCg)

plot(OMCg$size, OMCg$sizeNext, pch = 19,cex=1,ylim=c(-2,10),xlab="size",ylab="sizeNext",main="Growth")
y0<-predict(go1,newdata=x0);lines(x,y0,col="red")
y0<-predict(go2,newdata=x0);lines(x,y0,col="green")
y0<-predict(go3,newdata=x0);lines(x,y0,col="blue")
c(AIC(go1),AIC(go2),AIC(go3))
legend("topleft", legend = sprintf("%s: %s = %.1f",c("1","size","size+size2"), c("AIC"), c(AIC(go1),AIC(go2),AIC(go3))), col = c(2:4), lty = 1, xjust = 1, bg = "white")

# 2. growth regression from Cory's code
growth.reg=lm(sizeNext~size,data=OMCg)
summary(growth.reg)
params$growth.int=coefficients(growth.reg)[1]
params$growth.slope=coefficients(growth.reg)[2]
params$growth.sd=summary(growth.reg)$sigma

######Clonal Matrix: modified fecundity from Cory's code######
# whether or not a colony fissions
co1<-glm(fec1~1, data=OMC,family = binomial())
co2<-glm(fec1~size, data=OMC,family = binomial())
co3<-glm(fec1~size+size2, data=OMC,family = binomial())

# plot these models for comparison
par(mfrow=c(1,2),mar=c(4,4,2,1))
os <- order(OMC$size)  
osClon <- (OMC$fec1)[os] 
osSize<-(OMC$size)[os] 
binnedSize <- tapply(osSize, as.numeric(cut(osSize, breaks=20)), mean, na.rm = TRUE); # bin Size data
binnedClon <- tapply(osClon, as.numeric(cut(osSize, breaks=20)), mean, na.rm = TRUE) #bin Survival probabilities
plot(binnedSize, binnedClon, pch = 19, cex=2, xlab = "Size at t", ylab = "Producing clones at t+1", main = "Probability of Fission")
y0<-predict(co1,newdata=x0);y0<-exp(y0)/(exp(y0)+1);lines(x,y0,col="red")
y0<-predict(co2,newdata=x0);y0<-exp(y0)/(exp(y0)+1);lines(x,y0,col="green")
y0<-predict(co3,newdata=x0);y0<-exp(y0)/(exp(y0)+1);lines(x,y0,col="blue")
legend("topleft", legend = sprintf("%s: %s = %.1f",c("1","size","size+size2"), c("AIC"), c(AIC(co1),AIC(co2),AIC(co3))), col = c(2:4), lty = 1, xjust = 1, bg = "white")

clone=glm(fec1~size,data=OMC,family=binomial())
summary(clone)
params$clone.int=coefficients(clone)[1]
params$clone.slope=coefficients(clone)[2]

##How many colonies a fission produces
co1<-glm(fec4~1, data=OMC[OMC$fec1==1,],family = poisson())
co2<-glm(fec4~size, data=OMC[OMC$fec1==1],family = poisson())
co3<-glm(fec4~size+size2, data=OMC[OMC$fec1==1,],family = poisson())


# plot these models for comparison
os <- order(OMC$size)  
osClon <- (OMC$fec2)[os] 
osSize<-(OMC$size)[os] 
binnedSize <- tapply(osSize, as.numeric(cut(osSize, breaks=20)), mean, na.rm = TRUE); # bin Size data
binnedClon <- tapply(osClon, as.numeric(cut(osSize, breaks=20)), mean, na.rm = TRUE) #bin Survival probabilities
plot(binnedSize, binnedClon, pch = 19, cex=2, xlab = "Size at t", ylab = "number of clones produced", main = "How many clones produced")
y0<- lines(xx,predict(co1,data.frame(size=xx),type='response'), col='red')
y0<- lines(xx,predict(co2,data.frame(size=xx),type='response'), col='green')
y0<-lines(xx,predict(co3,data.frame(size=xx,size2=xx),data=OMC,type='response'), col='blue')
legend("topleft", legend = sprintf("%s: %s = %.1f",c("1","size","size+size2"), c("AIC"), c(AIC(co1),AIC(co2),AIC(co3))), col = c(2:4), lty = 1, xjust = 1, bg = "white")
###From Cory's exercise patge 7: d. update the regression for seed number to include only the individuals that flowered. this is easily done by using the code for the seed regression in section B.3. and changing the data argument to data=d[d$fec.flower==1,] 
clonenum=glm(fec4~size,data=OMC,family=poisson())
summary(clonenum)
params$clonenum.int=coefficients(clonenum)[1]
params$clonenum.slope=coefficients(clonenum)[2]

# 4. size distribution of clones
# in the dataframe, clones are those individuals who have "fissSM" for fateIPM
clonesize=lm(sizeNext~size,data=subset(OMC,OMC$fateIPM=="fissSM"))
summary(clonesize)
params$clonesize.int=coefficients(clonesize)[1]
params$clonesize.slope=coefficients(clonesize)[2]
params$clonesize.sd=summary(clonesize)$sigma


# 3. recruit rate
#replace recruit.int with number of recruits per "adult" averaged quad over years then total of quads
params$recruit.int=0.0905897


# 4. size distribution of recruits
# in the dataframe, recruits are those individuals who have a value for sizeNext but not for size
params$recruit.size.mean=mean(OMC$sizeNext[OMC$fateIPM=="recruit"])
params$recruit.size.sd=sd(OMC$sizeNext[OMC$fateIPM=="recruit"])

# 6. plot the models over the data - figure 2
par(mfrow=c(2,3),mar=c(4,4,2,1))
plot(OMCg$size,OMCg$sizeNext,main='Growth/Shrinkage/Stasis')  
lines(xx,predict(growth.reg,data.frame(size=xx)),col='red',lwd=3)
plot(OMC$size,jitter(OMC$surv),main='Survival') # jittered to see easier
lines(xx,predict(surv.reg,data.frame(size=xx),type='response'), col='red',lwd=3)
plot(OMC$size,jitter(OMC$fec1),main='Clones') # jittered to see easier
lines(xx,predict(clone,data.frame(size=xx),type='response'), col='red',lwd=3)
plot(OMC$size,OMC$fec2,main='Number of clones') # jittered to see easier
lines(xx,predict(clonenum,data.frame(size=xx),type='response'), col='red',lwd=3)
hist(OMC$sizeNext[OMC$fateIPM=="fissSM"],main='Size of Clones',freq=FALSE)
lines(xx,dnorm(xx,mean(params$clonesize.int+params$clonesize.slope),params$clonesize.sd), col='red',lwd=3)
## vital rate functions

# 1. probability of surviving
s.x=function(x,params) {
  u=exp(params$surv.int+params$surv.slope*x)
  return(u/(1+u))
}
# 2. growth function
g.yx=function(xp,x,params) {   		
  dnorm(xp,mean=params$growth.int+params$growth.slope*x,sd=params$growth.sd)
}
gs<-function(xp,x,g.yx,s.x){
  g.yx(xp,x,params)*s.x(x,params)
}
# 3. Clonal function      
p.clone.x=function(x,params) {
  u=exp(params$clone.int+params$clone.slope*x)
  return(u/(1+u))
}
c.yx=function(xp,x,params) {
  p.clone.x(x,params)*
    dnorm(xp,mean=params$clonesize.int+params$clonesize.slope*x,sd=params$clonesize.sd)*
    exp(params$clonenum.int+params$clonenum.slope*x)
}
c.y=function(x,params) {
  u=p.clone.x(x,params)*
    exp(params$clonenum.int+params$clonenum.slope*x)
  return (u)
}
# 3. recruitment function 
f.yx=function(xp,x,params) {
      u=dnorm(xp,mean=params$recruit.size.mean,sd=params$recruit.size.sd)*
      exp(params$recruit.int)
      return (u)
}
f.y=function(xp,params) {
  u=exp(params$recruit.int)
  return (u)
}
# 1. boundary points b, mesh points y and step size h
# integration limits - these limits span the range of sizes observed in the data set, and then some.
min.size<-min(OMC$size,na.rm=T)
max.size<-max(OMC$size,na.rm=T)+2
#From cory's code.  I used the one from Rob's code instead#min.size=.9*min(c(d$size,d$sizeNext),na.rm=T)
#From cory's code.  I used the one from Rob's code instead#max.size=1.1*max(c(d$size,d$sizeNext),na.rm=T)
# number of cells in the discretized kernel
n=100 
# #boundary points (the edges of the cells defining the kernel)
b=min.size+c(0:n)*(max.size-min.size)/n 
# mesh points (midpoints of the cells)
y=0.5*(b[1:n]+b[2:(n+1)])
# width of the cells
h=y[2]-y[1] 

# 2. make component kernels
# the function outer() evaluates the kernal at all pairwise combinations of the two vectors y and y. for the numerical integration, we're estimating the area of a rectangle under the curve. the hieghts of the rectangles are given by the outer function and the width of the rectangles is h. 
G=h*outer(y,y,g.yx,params=params) # growth kernel
S=s.x(y,params=params) # survival 
P=G # placeholder; we're about to redefine P on the next line
for(i in 1:n) P[,i]=G[,i]*S[i]  # growth/survival kernel

P2=outer(y,y,gs,g.yx=g.yx,s.x=s.x)*h
C=h*outer(y,y,c.yx,params=params) # reproduction kernel
F=h*outer(y,y,f.yx,params=params) # reproduction kernel

K=P+F+C #full kernel

# 1. get lamda,v,w  
lam=as.real(eigen(K)$values[1]) 
w.eigen=as.real(eigen(K)$vectors[,1])
stable.dist=w.eigen/sum(w.eigen) 
v.eigen=as.real(eigen(t(K))$vectors[,1])
repro.val=v.eigen/v.eigen[1] 

# 2. compute elasticity and sensitivity matrices
v.dot.w=sum(stable.dist*repro.val)*h
sens=outer(repro.val,stable.dist)/v.dot.w
elas=matrix(as.vector(sens)*as.vector(K)/lam,nrow=n)

# 3. plot results
par(mfrow=c(1,1))
filled.contour(y,y,t(K),color=heat.colors,nlevels=20,cex.lab=1.5,main="Kernel", xlab= "Size at t",ylab="Size at t+1")
filled.contour(y,y,t(elas),color=heat.colors,nlevels=20,cex.lab=1.5,main="Kernel", xlab= "Size at t",ylab="Size at t+1")
filled.contour(y,y,t(sens),color=heat.colors,nlevels=20,cex.lab=1.5,main="Kernel", xlab= "Size at t",ylab="Size at t+1")

######Cory's plots######
par(mfrow=c(2,3))
image(y,y,t(K), xlab="Size (t)",ylab="Size (t+1)",col=heat.colors(100), main="Kernel")
contour(y,y,t(K), add = TRUE, drawlabels = TRUE)
plot(y,stable.dist,xlab="Size",type="l",main="Stable size distribution")
plot(y,repro.val,xlab="Size",type="l",main="Reproductive values") 
image(y,y,t(elas),xlab="Size (t)",ylab="Size (t+1)",main="Elasticity")
image(y,y,t(sens),xlab="Size (t)",ylab="Size (t+1)", main="Sensitivity")
persp(y,y,t(K),phi=30,theta=300,col="lightgrey",border=NA,shade=0.75,main="Pmatrix:survival and growth", xlab= "Size at t",ylab="Size at t+1")
######constantcorrection for P#####
nvals <- colSums(P, na.rm = TRUE)
loc0 <- which(nvals == 0, arr.ind = TRUE)
if (length(loc0) > 0) {
  print("warnings - columns that sum to 0 or that have NAs - assuming survival is along the diagonal; plot your Pmatrix to check it")
  P[, loc0] <- 0
  P[cbind(loc0, loc0)] <- S
}
nvals <- colSums(P, na.rm = TRUE)
P<- t((t(P)/nvals) * S)
######################################################################################
###########################constantcorrection for F####################################
f<-vector(length=100)
for(i in 1:100) f[i]<- f.y(y,params=params)
correction.here <- f/colSums(F)
F<- t(t(F) * correction.here)
#######################################################################################
########################ConstantcorrectionforC#########################################
c<-c.y(y,params=params)
correction.here <- c/colSums(C)
C<- t(t(C) * correction.here)
########################################################################################

Kconstant=P+F+C
lam=as.real(eigen(Kconstant)$values[1]) 