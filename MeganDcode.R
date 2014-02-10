library(nlme)
growthonly<- read.csv("growthonly.csv")
min.size<-(-1.90897)
max.size<-9.152484
n=300
b=min.size+c(0:n)*(max.size-min.size)/n 
y=0.5*(b[1:n]+b[2:(n+1)])

##full no varFunc
growth.gnls.abc<-gnls(sizeNext~a+b*size+c*exp(-d*size),data=growthonly,params=list(a+b+c~utrans,d~1),start=list(rep(1,31)))
summary(growth.gnls.abc)
growth.int.abc<-c(coefficients(growth.gnls.abc)[1],coefficients(growth.gnls.abc)[1]+coefficients(growth.gnls.abc)[2:10])
growth.slope.abc<-c(coefficients(growth.gnls.abc)[11],coefficients(growth.gnls.abc)[11]+coefficients(growth.gnls.abc)[12:20])
growth.e1.abc<-c(coefficients(growth.gnls.abc)[21],coefficients(growth.gnls.abc)[21]+coefficients(growth.gnls.abc)[22:30])
growth.e2.abc<-coefficients(growth.gnls.abc)[31]

par(mfrow=c(3,4))
for (i in 1:10){
  x<-subset(growthonly,growthonly$utrans==levels(growthonly$utrans)[i])
  plot(x$size,x$sizeNext,xlim=c(-2,10),ylim=c(-2,10))
  lines(y,growth.int.abc[i]+growth.slope.abc[i]*y+growth.e1.abc[i]*exp(-growth.e2.abc*y))
}


##full with varExp
growth.gnls.ve.abc<-gnls(sizeNext~a+b*size+c*exp(-d*size),data=growthonly,params=list(a+b+c~utrans,d~1),start=list(rep(1,31)),weights=varExp(form=~size))
summary(growth.gnls.ve.abc)
growth.int.ve.abc<-c(coefficients(growth.gnls.ve.abc)[1],coefficients(growth.gnls.ve.abc)[1]+coefficients(growth.gnls.ve.abc)[2:10])
growth.slope.ve.abc<-c(coefficients(growth.gnls.ve.abc)[11],coefficients(growth.gnls.ve.abc)[11]+coefficients(growth.gnls.ve.abc)[12:20])
growth.e1.ve.abc<-c(coefficients(growth.gnls.ve.abc)[21],coefficients(growth.gnls.ve.abc)[21]+coefficients(growth.gnls.ve.abc)[22:30])
growth.e2.ve.abc<-coefficients(growth.gnls.ve.abc)[31]

par(mfrow=c(3,4))
for (i in 1:10){
  x<-subset(growthonly,growthonly$utrans==levels(growthonly$utrans)[i])
  plot(x$size,x$sizeNext,xlim=c(-2,10),ylim=c(-2,10))
  lines(y,growth.int.ve.abc[i]+growth.slope.ve.abc[i]*y+growth.e1.ve.abc[i]*exp(-growth.e2.ve.abc*y))
}
anova(growth.gnls.abc,growth.gnls.ve.abc)
#varExp significantly improves fit deltaAIC=15.809, p-value<0.0001





##drop utrans:c

growth.gnls.ve.ab<-gnls(sizeNext~a+b*size+c*exp(-d*size),data=growthonly,params=list(a+b~utrans,c+d~1),start=list(rep(1,22)),weights=varExp(form=~size))
summary(growth.gnls.ve.ab)
growth.int.ve.ab<-c(coefficients(growth.gnls.ve.ab)[1],coefficients(growth.gnls.ve.ab)[1]+coefficients(growth.gnls.ve.ab)[2:10])
growth.slope.ve.ab<-c(coefficients(growth.gnls.ve.ab)[11],coefficients(growth.gnls.ve.ab)[11]+coefficients(growth.gnls.ve.ab)[12:20])
growth.e1.ve.ab<-coefficients(growth.gnls.ve.ab)[21]
growth.e2.ve.ab<-coefficients(growth.gnls.ve.ab)[22]

par(mfrow=c(3,4))
for (i in 1:10){
  x<-subset(growthonly,growthonly$utrans==levels(growthonly$utrans)[i])
  plot(x$size,x$sizeNext,xlim=c(-2,10),ylim=c(-2,10))
  lines(y,growth.int.ve.ab[i]+growth.slope.ve.ab[i]*y+growth.e1.ve.ab*exp(-growth.e2.ve.ab*y))
}

anova(growth.gnls.ve.abc,growth.gnls.ve.ab)
##keep c:utrans deltaAIC=5.055 p-value=0.0061


##drop a:utrans
growth.gnls.ve.bc<-gnls(sizeNext~a+b*size+c*exp(-d*size),data=growthonly,params=list(b+c~utrans,a+d~1),start=list(rep(1,22)),weights=varExp(form=~size))
summary(growth.gnls.ve.bc)
growth.int.ve.bc<-coefficients(growth.gnls.ve.bc)[21]
growth.slope.ve.bc<-c(coefficients(growth.gnls.ve.bc)[1],coefficients(growth.gnls.ve.bc)[1]+coefficients(growth.gnls.ve.bc)[2:10])
growth.e1.ve.bc<-c(coefficients(growth.gnls.ve.bc)[11],coefficients(growth.gnls.ve.bc)[11]+coefficients(growth.gnls.ve.bc)[12:20])
growth.e2.ve.bc<-coefficients(growth.gnls.ve.bc)[22]

par(mfrow=c(3,4))
for (i in 1:10){
  x<-subset(growthonly,growthonly$utrans==levels(growthonly$utrans)[i])
  plot(x$size,x$sizeNext,xlim=c(-2,10),ylim=c(-2,10))
  lines(y,growth.int.ve.bc+growth.slope.ve.bc[i]*y+growth.e1.ve.bc[i]*exp(-growth.e2.ve.bc*y))
}

anova(growth.gnls.ve.abc,growth.gnls.ve.bc)
##keep a:utrans deltaAIC=12.481 p-value=4e-04

##a:utrans only
growth.gnls.ve.a<-gnls(sizeNext~a+b*size+c*exp(-d*size),data=growthonly,params=list(a~utrans,b+c+d~1),start=list(rep(1,13)),weights=varExp(form=~size))
summary(growth.gnls.ve.a)
growth.int.ve.a<-c(coefficients(growth.gnls.ve.a)[1],coefficients(growth.gnls.ve.a)[1]+coefficients(growth.gnls.ve.a)[2:10])
growth.slope.ve.a<-coefficients(growth.gnls.ve.a)[11]
growth.e1.ve.a<-coefficients(growth.gnls.ve.a)[12]
growth.e2.ve.a<-coefficients(growth.gnls.ve.a)[13]

par(mfrow=c(3,4))
for (i in 1:10){
  x<-subset(growthonly,growthonly$utrans==levels(growthonly$utrans)[i])
  plot(x$size,x$sizeNext,xlim=c(-2,10),ylim=c(-2,10))
  lines(y,growth.int.ve.a[i]+growth.slope.ve.a*y+growth.e1.ve.a*exp(-growth.e2.ve.a*y))
}
anova(growth.gnls.ve.ab,growth.gnls.ve.a)
##keep b deltaAIC=20.291 p-value<0.0001

##b:utrans only
growth.gnls.ve.b<-gnls(sizeNext~a+b*size+c*exp(-d*size),data=growthonly,params=list(b~utrans,a+c+d~1),start=list(rep(1,13)),weights=varExp(form=~size))
summary(growth.gnls.ve.b)
growth.int.ve.b<-coefficients(growth.gnls.ve.b)[11]
growth.slope.ve.b<-c(coefficients(growth.gnls.ve.b)[1],coefficients(growth.gnls.ve.b)[1]+coefficients(growth.gnls.ve.b)[2:10])
growth.e1.ve.b<-coefficients(growth.gnls.ve.b)[12]
growth.e2.ve.b<-coefficients(growth.gnls.ve.b)[13]

par(mfrow=c(3,4))
for (i in 1:10){
  x<-subset(growthonly,growthonly$utrans==levels(growthonly$utrans)[i])
  plot(x$size,x$sizeNext,xlim=c(-2,10),ylim=c(-2,10))
  lines(y,growth.int.ve.b+growth.slope.ve.b[i]*y+growth.e1.ve.b*exp(-growth.e2.ve.b*y))
}
anova(growth.gnls.ve.ab,growth.gnls.ve.b)
#keep a deltaAIC=47.827 p-value<0.0001
anova(growth.gnls.ve.bc,growth.gnls.ve.b)
#keep c deltaAIC=40.689 p-value<0.0001


##c:utrans only
growth.gnls.ve.c<-gnls(sizeNext~a+b*size+c*exp(-d*size),data=growthonly,params=list(c~utrans,a+b+d~1),start=list(rep(1,13)),weights=varExp(form=~size))
summary(growth.gnls.ve.c)
growth.int.ve.c<-coefficients(growth.gnls.ve.c)[11]
growth.slope.ve.c<-coefficients(growth.gnls.ve.c)[12]
growth.e1.ve.c<-c(coefficients(growth.gnls.ve.c)[1],coefficients(growth.gnls.ve.c)[1]+coefficients(growth.gnls.ve.c)[2:10])
growth.e2.ve.c<-coefficients(growth.gnls.ve.c)[13]

par(mfrow=c(3,4))
for (i in 1:10){
  x<-subset(growthonly,growthonly$utrans==levels(growthonly$utrans)[i])
  plot(x$size,x$sizeNext,xlim=c(-2,10),ylim=c(-2,10))
  lines(y,growth.int.ve.c+growth.slope.ve.c*y+growth.e1.ve.c[i]*exp(-growth.e2.ve.c*y))
}
anova(growth.gnls.ve.bc,growth.gnls.ve.c)
#drop b:utrans deltaAIC=2.506 p-value=0.0782

AIC(growth.gnls.ve.abc,growth.gnls.ve.ab,growth.gnls.ve.bc,growth.gnls.ve.a,growth.gnls.ve.b,growth.gnls.ve.c)



##Look for best varFunc = varExp vs formula to varFixed
                         
                               
                               
par(mfrow=c(1,1))
plot(growthonly$size,residuals(growth.gnls.ve.abc),xlim=c(-2,10),ylim=c(-4,4),main="",xlab="Size t",ylab="Size t+1")
sd2<-(summary(growth.gnls.ve.abc)$sigma)^2
varF<-coef(growth.gnls.ve.abc$modelStruct$varStruct,"exp")
lines(y,sqrt(sd2*exp(2*varF*y)))

growth.exp.gnls<-gnls(sizeNext~a+b*size+c*exp(-d*size),data=growthonly,params=a+b+c+d~1,start=list(a=1,b=1,c=1,d=1),weights=varExp(form=~size))
growth.int<-coefficients(growth.exp.gnls)[1]
growth.slope<-coefficients(growth.exp.gnls)[2]
growth.e1<-coefficients(growth.exp.gnls)[3]
growth.e2<-coefficients(growth.exp.gnls)[4]


par(mfrow=c(3,4))
for (i in 1:10){
  x<-subset(growthonly,growthonly$utrans==levels(growthonly$utrans)[i])
  #splinex<-smooth.spline(x$size,x$sizeNext)
  #plot(splinex,type="l",lwd=2,xlim=c(-2,10),ylim=c(-2,10),xlab="Size t",ylab="Size t+1",main=levels(growthonly$utrans)[i])
  plot(x$size,x$sizeNext,pch=18,xlim=c(-2,10),ylim=c(-2,10),xlab="Size t",ylab="Size t+1",main=levels(growthonly$utrans)[i])
  #lines(y,growth.int.ve.abc[i]+growth.slope.ve.abc[i]*y+growth.e1.ve.abc[i]*exp(-growth.e2.ve.abc*y),col="red")
  #lines(y,growth.int.ve.ab[i]+growth.slope.ve.ab[i]*y+growth.e1.ve.ab*exp(-growth.e2.ve.ab*y),col="blue")
  #lines(y,growth.int.ve.bc+growth.slope.ve.bc[i]*y+growth.e1.ve.bc[i]*exp(-growth.e2.ve.bc*y),col="green")
  lines(y,growth.int.ve.a[i]+growth.slope.ve.a*y+growth.e1.ve.a*exp(-growth.e2.ve.a*y))
  #lines(y,growth.int.ve.b+growth.slope.ve.b[i]*y+growth.e1.ve.b*exp(-growth.e2.ve.b*y),col="magenta")
  #lines(y,growth.int.ve.c+growth.slope.ve.c*y+growth.e1.ve.c[i]*exp(-growth.e2.ve.c*y),col="orange")
  #lines(y,growth.int+growth.slope*y+growth.e1*exp(-growth.e2*y),col="red")
  
}
                               