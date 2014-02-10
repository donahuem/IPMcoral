x<-seq(from=-2.4,to=10.5,by=0.01)#Range of values for lines added to plots

plot(MCg$size,MCg$sizeNext,pch=18,col=c("black","red","green")[MCg$fate])
test<-lm(sizeNext~size*size2*fate,data=MCg)
drop1(test,test="Chisq")
#drop size:size2:fate p-value=0.6476
test<-lm(sizeNext~size*size2+size*fate+size2*fate,data=MCg)
drop1(test,test="Chisq")
#drop size2:fate p-value=0.9414
test<-lm(sizeNext~size*size2+size*fate,data=MCg)
drop1(test,test="Chisq")
#drop size:size2 p-value=0.0944908
test<-lm(sizeNext~size*fate+size2,data=MCg)
drop1(test,test="Chisq")
#keep size 2 p-value=3.122e-06 and size:fate p-value=0.0002565
plot(MCg$size,MCg$sizeNext,pch=18,col=c("black","red","green")[MCg$fate])
lines(y,(coefficients(test)[1]+coefficients(test)[2]*y)+coefficients(test)[5]*y^2,lwd=2)
lines(y,coefficients(test)[1]+coefficients(test)[3]+(coefficients(test)[2]+coefficients(test)[6])*y+coefficients(test)[5]*y^2,lwd=2,col="red")
lines(y,coefficients(test)[1]+coefficients(test)[4]+(coefficients(test)[2]+coefficients(test)[7])*y+coefficients(test)[5]*y^2,lwd=2,col="green")
##Compare linear model 
fategrow<-lm(sizeNext~size+fate,data=growdat)
growth<-subset(MCg,MCg$fate=="growth")
plot(growth$size,growth$sizeNext)
test1<-lm(sizeNext~size+size2,data=growth)
bleh<-data.frame(size=y,size2=y^2)
lines(y,predict(test1,newdata=bleh),col="purple")
lines(y,(coefficients(test)[1]+coefficients(test)[2]*y)+coefficients(test)[5]*y^2)
lines(x,coefficients(fategrow)[1]+coefficients(fategrow)[2]*x,col="green",lwd=1)
legend("bottomright",legend=c("single","byfate","linear"),col=c("black","purple","green"),lwd=2)

fission<-subset(MCg,MCg$fate=="fissLG")
plot(fission$size,fission$sizeNext)
test2<-lm(sizeNext~size+size2,data=fission)
lines(y,predict(test2,newdata=bleh),col="blue")
lines(y,coefficients(test)[1]+coefficients(test)[3]+(coefficients(test)[2]+coefficients(test)[6])*y+coefficients(test)[5]*y^2,col="red")
lines(x,coefficients(fategrow)[1]+coefficients(fategrow)[3]+coefficients(fategrow)[2]*x,col="green",lwd=1)
legend("bottomright",legend=c("single","byfate","linear"),col=c("red","blue","green"),lwd=2)

fusion<-subset(MCg,MCg$fate=="fuseLG")
plot(fusion$size,fusion$sizeNext)
test3<-lm(sizeNext~size+size2,data=fusion)
lines(y,predict(test3,newdata=bleh),col="magenta")
lines(y,coefficients(test)[1]+coefficients(test)[4]+(coefficients(test)[2]+coefficients(test)[7])*y+coefficients(test)[5]*y^2,col="green")
lines(x,coefficients(fategrow)[1]+coefficients(fategrow)[4]+coefficients(fategrow)[2]*x,col="blue",lwd=1)
legend("bottomright",legend=c("single","byfate","linear"),col=c("green","magenta","blue"),lwd=2)



###Include Transect
test<-lm(sizeNext~size*utrans*fate*size2,data=MCg)
drop1(test,test="Chisq")
#drop size:utrans:fate:size2 p-value=0.3683
test<-lm(sizeNext~size*utrans+size*fate+size*size2+size2*fate+size2*utrans+fate*utrans,data=MCg)
drop1(test,test="Chisq")
#drop fate:size2 p-value=0.85252
test<-lm(sizeNext~size*utrans+size*fate+size*size2+size2*utrans+size*utrans+fate*utrans,data=MCg)
drop1(test,test="Chisq")
#drop utrans:fate p-value=0.781489
test<-lm(sizeNext~size*utrans+size*fate+size*size2+size2*utrans,data=MCg)
drop1(test,test="Chisq")
#keep everything else


##Try to plot it
##All together
summary(test)
growth.int<-c(coefficients(test)[1],coefficients(test)[1]+coefficients(test)[3:(nx+1)])
growth.slope1<-c(coefficients(test)[2],coefficients(test)[2]+coefficients(test)[15:23])
growth.slope2<-c(coefficients(test)[14],coefficients(test)[14]+coefficients(test)[27:35])
growth.slope3<-coefficients(test)[26]
fiss.int<-c(growth.int+coefficients(test)[12])
fiss.slope1<-c(growth.slope1+coefficients(test)[24])
fuse.int<-c(growth.int+coefficients(test)[13])
fuse.slope1<-c(growth.slope1+coefficients(test)[25])

par(mfrow=c(3,4))
for (i in 1:10){
  plot(MCg$size[MCg$utrans==levels(MCg$utrans)[i]],MCg$sizeNext[MCg$utrans==levels(MCg$utrans)[i]],pch=18,col=c("black","red","green")[MCg$fate],main="",ylab="Size (t+1)",xlab="Size (t)",xlim=c(-2,10),ylim=c(-2,10))
  lines(y,growth.int[i]+growth.slope1[i]*y+growth.slope2[i]*y^2+growth.slope3*y*y^2,col="black")
  lines(y,fiss.int[i]+fiss.slope1[i]*y+growth.slope2[i]*y^2+growth.slope3*y*y^2,col="red")
  lines(y,fuse.int[i]+fuse.slope1[i]*y+growth.slope2[i]*y^2+growth.slope3*y*y^2,col="green")
}
  

### try without size:size2 interaction
test<-lm(sizeNext~size*utrans+size*fate+size2*utrans,data=MCg)
drop1(test,test="Chisq")
#keep everything else


##Try to plot it
##All together
summary(test)
growth.int<-c(coefficients(test)[1],coefficients(test)[1]+coefficients(test)[3:(nx+1)])
growth.slope1<-c(coefficients(test)[2],coefficients(test)[2]+coefficients(test)[15:23])
growth.slope2<-c(coefficients(test)[14],coefficients(test)[14]+coefficients(test)[26:34])
fiss.int<-c(growth.int+coefficients(test)[12])
fiss.slope1<-c(growth.slope1+coefficients(test)[24])
fuse.int<-c(growth.int+coefficients(test)[13])
fuse.slope1<-c(growth.slope1+coefficients(test)[25])

par(mfrow=c(3,4))
for (i in 1:10){
  plot(MCg$size[MCg$utrans==levels(MCg$utrans)[i]],MCg$sizeNext[MCg$utrans==levels(MCg$utrans)[i]],pch=18,col=c("black","red","green")[MCg$fate],main="",ylab="Size (t+1)",xlab="Size (t)",xlim=c(-2,10),ylim=c(-2,10))
  lines(y,growth.int[i]+growth.slope1[i]*y+growth.slope2[i]*y^2,col="black")
  lines(y,fiss.int[i]+fiss.slope1[i]*y+growth.slope2[i]*y^2,col="red")
  lines(y,fuse.int[i]+fuse.slope1[i]*y+growth.slope2[i]*y^2,col="green")
}


##Try without utrans:size2
test<-lm(sizeNext~size*utrans+size*fate+size2,data=MCg)
summary(test)

growth.int<-c(coefficients(test)[1],coefficients(test)[1]+coefficients(test)[3:(nx+1)])
growth.slope1<-c(coefficients(test)[2],coefficients(test)[2]+coefficients(test)[15:23])
growth.slope2<-coefficients(test)[14]
fiss.int<-c(growth.int+coefficients(test)[12])
fiss.slope1<-c(growth.slope1+coefficients(test)[24])
fuse.int<-c(growth.int+coefficients(test)[13])
fuse.slope1<-c(growth.slope1+coefficients(test)[25])

par(mfrow=c(3,4))
for (i in 1:10){
  plot(MCg$size[MCg$utrans==levels(MCg$utrans)[i]],MCg$sizeNext[MCg$utrans==levels(MCg$utrans)[i]],pch=18,col=c("black","red","green")[MCg$fate],main="",ylab="Size (t+1)",xlab="Size (t)",xlim=c(-2,10),ylim=c(-2,10))
  lines(y,growth.int[i]+growth.slope1[i]*y+growth.slope2*y^2,col="black")
  lines(y,fiss.int[i]+fiss.slope1[i]*y+growth.slope2*y^2,col="red")
  lines(y,fuse.int[i]+fuse.slope1[i]*y+growth.slope2*y^2,col="green")
}




###Should be still include varFunc?
test<-gls(sizeNext~size*utrans+size*fate+size2,data=MCg)
test.exp<-gls(sizeNext~size*utrans+size*fate+size2,weights=varExp(form=~size),data=MCg)
anova(test, test.exp)
#varExp significantly improves fit p-value=0.0001, deltaAIC=30.092

##Compare residuals to model with varExp
grow.UF.ves<-gls(sizeNext~size*utrans+size*fate,data=growdat,weights=varExp(form=~size))
E<-resid(grow.UF.ves,type="normalized")
coplot(E~size|utrans,data=growdat,ylab="normalized residuals")
E<-resid(test.exp,type="normalized")
coplot(E~size|utrans,data=growdat,ylab="normalized residuals")
##residuals look better with size2 included




###Fit 3 models for fates
growthonly<-subset(MCg,MCg$fateIPM=="growth")
fissonly<-subset(MCg,MCg$fateIPM=="fissLG")
fuseonly<-subset(MCg,MCg$fateIPM=="fuseLG"|MCg$fateIPM=="M")
growth<-gls(sizeNext~size*utrans+size2*utrans,data=growthonly)
growthvar<-gls(sizeNext~size*utrans+size2*utrans,data=growthonly,weights=varExp(form=~size))
anova(growth,growthvar)
#keep varExp p-value=1e-04
growthvar<-gls(sizeNext~size*utrans+size2*utrans,data=growthonly,weights=varExp(form=~size),method="ML")
growthvar1<-gls(sizeNext~size*utrans+size2,data=growthonly,weights=varExp(form=~size),method="ML")
anova(growthvar,growthvar1)
#keep size2:utrans p-value=4e-04

#model with size2:utrans would not converge
fiss<-gls(sizeNext~size*utrans+size2,data=fissonly)
fissvar<-gls(sizeNext~size*utrans+size2,data=fissonly,weights=varExp(form=~size))
anova(fiss, fissvar)
#varExp does not improve fit p-value=0.1036
fissvar<-gls(sizeNext~size*utrans+size2,data=fissonly,weights=varExp(form=~size),method="ML")
fissvar1<-gls(sizeNext~size*utrans,data=fissonly,weights=varExp(form=~size),method="ML")
anova(fissvar,fissvar1)
#drop size p-value=0.4356
fissvar2<-gls(sizeNext~size+utrans,data=fissonly,weights=varExp(form=~size),method="ML")
anova(fissvar1,fissvar2)
#drop size:utrans p-value=0.0522
fissvar3<-gls(sizeNext~size,data=fissonly,weights=varExp(form=~size),method="ML")
anova(fissvar2,fissvar3)
#drop utrans p-value=0.7344
fissvar4<-gls(sizeNext~1,data=fissonly,weights=varExp(form=~size),method="ML")
anova(fissvar3,fissvar4)
#keep size p-value<0.0001

fuse<-gls(sizeNext~size+utrans+size2,data=fuseonly)
fusevar<-gls(sizeNext~size+utrans+size2,data=fuseonly,weights=varExp(form=~size))
anova(fuse,fusevar)
#varExp improves fit p-value=0.0108
fusevar<-gls(sizeNext~size+utrans+size2,data=fuseonly,weights=varExp(form=~size),method="ML")
fusevar1<-gls(sizeNext~size+utrans,data=fuseonly,weights=varExp(form=~size),method="ML")
anova(fusevar,fusevar1)
#drop size2 p-value=0.0514
fusevar2<-gls(sizeNext~size,data=fuseonly,weights=varExp(form=~size),method="ML")
anova(fusevar1,fusevar2)
#drop utrans p-value=0.0586
fusevar3<-gls(sizeNext~1,data=fuseonly,weights=varExp(form=~size),method="ML")
anova(fusevar2,fusevar3)
#keep size p-value<0.0001

##Refit everything using REML
growthvar<-gls(sizeNext~size*utrans+size2,data=growthonly,weights=varExp(form=~size))
fissvar3<-gls(sizeNext~size,data=fissonly,weights=varExp(form=~size))
fusevar2<-gls(sizeNext~size,data=fuseonly,weights=varExp(form=~size))
growth.int<-c(coefficients(growthvar)[1],coefficients(growthvar)[1]+coefficients(growthvar)[3:11])
growth.slope1<-c(coefficients(growthvar)[2],coefficients(growthvar)[2]+coefficients(growthvar)[13:21])
growth.slope2<-coefficients(growthvar)[12]
fiss.int<-coefficients(fissvar3)[1]
fiss.slope<-coefficients(fissvar3)[2]
fuse.int<-coefficients(fusevar2)[1]
fuse.slope<-coefficients(fusevar2)[2]
par(mfrow=c(3,4))
for (i in 1:10){
  plot(MCg$size[MCg$utrans==levels(MCg$utrans)[i]],MCg$sizeNext[MCg$utrans==levels(MCg$utrans)[i]],col=c("black","red","green")[MCg$fate],xlim=c(-2,10),ylim=c(-2,10),xlab="Size t",ylab="Size t+1")
  lines(y,growth.int[i]+growth.slope1[i]*y+growth.slope2*y^2,col="black")
  lines(y,fiss.int+fiss.slope*y,col="red")
  lines(y,fuse.int+fiss.slope*y,col="green")
}



#nonlinear

growthonly<-subset(MCg,MCg$fateIPM=="growth")

#fis nonlinear least squares regression
growth.exp<-nls(sizeNext~a+b*size+c*exp(-d*size),data=growthonly,start=list(a=1,b=1,c=1,d=1),nls.control(maxiter=10000,minFactor=1/(2^12)))
plot(growthonly$size,growthonly$sizeNext)
dsize<-data.frame(size=y)
lines(y,predict(growth.exp,newdata=dsize),lwd=2)

par(mfrow=c(3,4))
for (i in 1:10){
  x<-subset(growthonly,growthonly$utrans==levels(growthonly$utrans)[i])
  plot(x$size,x$sizeNext,xlim=c(-2,10),ylim=c(-2,10))
  dsize<-data.frame(size=y)
  lines(y,predict(growth.exp,newdata=dsize),lwd=2)
}


##make dummy variables to fit main effects of utrans and utrans:size interaction in nls()
index<-levels(growthonly$utrans)
values<-c(1,0,0,0,0,0,0,0,0,0)
growthonly$D1<-values[match(growthonly$utrans,index)]
values<-c(0,1,0,0,0,0,0,0,0,0)
growthonly$D2<-values[match(growthonly$utrans,index)]
values<-c(0,0,1,0,0,0,0,0,0,0)
growthonly$D3<-values[match(growthonly$utrans,index)]
values<-c(0,0,0,1,0,0,0,0,0,0)
growthonly$D4<-values[match(growthonly$utrans,index)]
values<-c(0,0,0,0,1,0,0,0,0,0)
growthonly$D5<-values[match(growthonly$utrans,index)]
values<-c(0,0,0,0,0,1,0,0,0,0)
growthonly$D6<-values[match(growthonly$utrans,index)]
values<-c(0,0,0,0,0,0,1,0,0,0)
growthonly$D7<-values[match(growthonly$utrans,index)]
values<-c(0,0,0,0,0,0,0,1,0,0)
growthonly$D8<-values[match(growthonly$utrans,index)]
values<-c(0,0,0,0,0,0,0,0,1,0)
growthonly$D9<-values[match(growthonly$utrans,index)]
values<-c(0,0,0,0,0,0,0,0,0,1)
growthonly$D10<-values[match(growthonly$utrans,index)]

growth.exp.U0<-nls(sizeNext~a+b*size+c*D2+d*D3+e*D4+f*D5+g*D6+h*D7+i*D8+j*D9+k*D10+l*exp(-m*size),data=growthonly,start=list(a=1,b=1,c=1,d=1,e=1,f=1,g=1,h=1,i=1,j=1,k=1,l=1,m=1),nls.control(maxiter=10000,minFactor=1/(2^20)))
summary(growth.exp.U0)

growth.int<-c(coefficients(growth.exp.U0)[1],coefficients(growth.exp.U0)[1]+coefficients(growth.exp.U0)[3:11])
growth.slope<-coefficients(growth.exp.U0)[2]
growth.e1<-coefficients(growth.exp.U0)[12]
growth.e2<-coefficients(growth.exp.U0)[13]

par(mfrow=c(3,4))
for (i in 1:10){
  x<-subset(growthonly,growthonly$utrans==levels(growthonly$utrans)[i])
  plot(x$size,x$sizeNext,xlim=c(-2,10),ylim=c(-2,10))
  lines(y,growth.int[i]+growth.slope*y+growth.e1*exp(-growth.e2*y))
}

growth.exp.U<-nls(sizeNext~a+b*size+c*D2+d*D3+e*D4+f*D5+g*D6+h*D7+i*D8+j*D9+k*D10+l*exp(-m*size)+n*D2*size+o*D3*size+p*D4*size+q*D5*size+r*D6*size+s*D7*size+t*D8*size+u*D9*size+v*D10*size,data=growthonly,start=list(a=1,b=1,c=1,d=1,e=1,f=1,g=1,h=1,i=1,j=1,k=1,l=1,m=1,n=1,o=1,p=1,q=1,r=1,s=1,t=1,u=1,v=1),nls.control(maxiter=10000,minFactor=1/(2^20)))
summary(growth.exp.U)

growth.int1<-c(coefficients(growth.exp.U)[1],coefficients(growth.exp.U)[1]+coefficients(growth.exp.U)[3:11])
growth.slope1<-c(coefficients(growth.exp.U)[2],coefficients(growth.exp.U)[2]+coefficients(growth.exp.U)[14:22])
growth.e11<-coefficients(growth.exp.U)[12]
growth.e21<-coefficients(growth.exp.U)[13]

par(mfrow=c(3,4))
for (i in 1:10){
  x<-subset(growthonly,growthonly$utrans==levels(growthonly$utrans)[i])
  plot(x$size,x$sizeNext,xlim=c(-2,10),ylim=c(-2,10))
  lines(y,growth.int1[i]+growth.slope1[i]*y+growth.e11*exp(-growth.e21*y))
}

###Visualize all together
par(mfrow=c(3,4))
for (i in 1:10){
  x<-subset(growthonly,growthonly$utrans==levels(growthonly$utrans)[i])
  plot(x$size,x$sizeNext,xlim=c(-2,10),ylim=c(-2,10))
  lines(y,predict(growth.exp,newdata=dsize),lwd=2,col="black")
  lines(y,growth.int[i]+growth.slope*y+growth.e1*exp(-growth.e2*y),col="green")
  lines(y,growth.int1[i]+growth.slope1[i]*y+growth.e11*exp(-growth.e21*y),col="red")
  lines(y,growth.int2[i]+growth.slope2[i]*y+growth.e12*exp(-growth.e22*y),col="purple")
}



##Add variance component cannot use varFunc in nls() can use weights=
summary(growth.exp.U)$sigma
0.774043^2
growthonly$resid<-residuals(growth.exp.U)
plot(growthonly$size,growthonly$resid)
vw<-nls(resid~0.5991426*exp(2*a*size),data=growthonly,start=list(a=1))
plot(growthonly$size,growthonly$resid)
lines(y,predict(vw,dsize))
summary(vw)
ve<-varExp(coefficients(vw),form=~size)
ve<-Initialize(ve,growthonly)
ve<-varWeights(ve)

growth.exp.Uvar<-nls(sizeNext~a+b*size+c*D2+d*D3+e*D4+f*D5+g*D6+h*D7+i*D8+j*D9+k*D10+l*exp(-m*size)+n*D2*size+o*D3*size+p*D4*size+q*D5*size+r*D6*size+s*D7*size+t*D8*size+u*D9*size+v*D10*size,data=growthonly,start=list(a=1,b=1,c=1,d=1,e=1,f=1,g=1,h=1,i=1,j=1,k=1,l=1,m=1,n=1,o=1,p=1,q=1,r=1,s=1,t=1,u=1,v=1),nls.control(maxiter=10000,minFactor=1/(2^20)),weights=ve)
growth.int2<-c(coefficients(growth.exp.Uvar)[1],coefficients(growth.exp.Uvar)[1]+coefficients(growth.exp.Uvar)[3:11])
growth.slope2<-c(coefficients(growth.exp.Uvar)[2],coefficients(growth.exp.Uvar)[2]+coefficients(growth.exp.Uvar)[14:22])
growth.e12<-coefficients(growth.exp.Uvar)[12]
growth.e22<-coefficients(growth.exp.Uvar)[13]

par(mfrow=c(3,4))
for (i in 1:10){
x<-subset(growthonly,growthonly$utrans==levels(growthonly$utrans)[i])
plot(x$size,x$sizeNext,xlim=c(-2,10),ylim=c(-2,10))
lines(y,growth.int2[i]+growth.slope2[i]*y+growth.e12*exp(-growth.e22*y))
}
library(lmtest)
lrtest(growth.exp.Uvar,growth.exp.U)


##I can do all of this in gnls just had to figure out syntax for params
growth.exp.gnls<-gnls(sizeNext~a+b*size+c*exp(-d*size),data=growthonly,params=a+b+c+d~1,start=list(a=1,b=1,c=1,d=1))
growth.exp.U.gnls<-gnls(sizeNext~a+b*size+c*exp(-d*size),data=growthonly,params=list(a+b~utrans,c+d~1),start=list(a=c(1,1,1,1,1,1,1,1,1,1),b=c(1,1,1,1,1,1,1,1,1,1),c=1,d=1))
growth.exp.UvarExp.gnls<-gnls(sizeNext~a+b*size+c*exp(-d*size),data=growthonly,params=list(a+b~utrans,c+d~1),start=list(a=c(1,1,1,1,1,1,1,1,1,1),b=c(1,1,1,1,1,1,1,1,1,1),c=1,d=1),weights=varExp(form=~size))

growth.int<-c(coefficients(growth.exp.UvarExp.gnls)[1],coefficients(growth.exp.UvarExp.gnls)[1]+coefficients(growth.exp.UvarExp.gnls)[2:10])
growth.slope<-c(coefficients(growth.exp.UvarExp.gnls)[11],coefficients(growth.exp.UvarExp.gnls)[11]+coefficients(growth.exp.UvarExp.gnls)[12:20])
growth.e1<-coefficients(growth.exp.UvarExp.gnls)[21]
growth.e2<-coefficients(growth.exp.UvarExp.gnls)[22]

par(mfrow=c(3,4))
for (i in 1:10){
  x<-subset(growthonly,growthonly$utrans==levels(growthonly$utrans)[i])
  plot(x$size,x$sizeNext,xlim=c(-2,10),ylim=c(-2,10))
  lines(y,growth.int[i]+growth.slope[i]*y+growth.e1*exp(-growth.e2*y))
}

growth.exp.UvarExp.gnls1<-gnls(sizeNext~a+b*size+c*exp(-d*size),data=growthonly,params=list(b+c~utrans,a+d~1),start=list(a=1,b=c(1,1,1,1,1,1,1,1,1,1),c=c(1,1,1,1,1,1,1,1,1,1),d=1),weights=varExp(form=~size))

growth.int1<-coefficients(growth.exp.UvarExp.gnls1)[21]
growth.slope1<-c(coefficients(growth.exp.UvarExp.gnls1)[1],coefficients(growth.exp.UvarExp.gnls1)[1]+coefficients(growth.exp.UvarExp.gnls1)[2:10])
growth.e11<-c(coefficients(growth.exp.UvarExp.gnls1)[11],coefficients(growth.exp.UvarExp.gnls1)[11]+coefficients(growth.exp.UvarExp.gnls1)[12:20])
growth.e21<-coefficients(growth.exp.UvarExp.gnls1)[22]

par(mfrow=c(3,4))
for (i in 1:10){
  x<-subset(growthonly,growthonly$utrans==levels(growthonly$utrans)[i])
  plot(x$size,x$sizeNext,xlim=c(-2,10),ylim=c(-2,10))
  lines(y,growth.int1+growth.slope1[i]*y+growth.e11[i]*exp(-growth.e21*y))
}


growth.exp.UvarExp.gnls2<-gnls(sizeNext~a+b*size+c*exp(-d*size),data=growthonly,params=list(a+b+c~utrans,d~1),start=list(a=c(1,1,1,1,1,1,1,1,1,1),b=c(1,1,1,1,1,1,1,1,1,1),c=c(1,1,1,1,1,1,1,1,1,1),d=1),weights=varExp(form=~size))
summary(growth.exp.UvarExp.gnls2)
growth.int2<-c(coefficients(growth.exp.UvarExp.gnls2)[1],coefficients(growth.exp.UvarExp.gnls2)[1]+coefficients(growth.exp.UvarExp.gnls2)[2:10])
growth.slope2<-c(coefficients(growth.exp.UvarExp.gnls2)[11],coefficients(growth.exp.UvarExp.gnls2)[11]+coefficients(growth.exp.UvarExp.gnls2)[12:20])
growth.e12<-c(coefficients(growth.exp.UvarExp.gnls2)[21],coefficients(growth.exp.UvarExp.gnls2)[21]+coefficients(growth.exp.UvarExp.gnls2)[22:30])
growth.e22<-coefficients(growth.exp.UvarExp.gnls2)[31]

par(mfrow=c(3,4))
for (i in 1:10){
  x<-subset(growthonly,growthonly$utrans==levels(growthonly$utrans)[i])
  plot(x$size,x$sizeNext,xlim=c(-2,10),ylim=c(-2,10))
  lines(y,growth.int2[i]+growth.slope2[i]*y+growth.e12[i]*exp(-growth.e22*y))
}

growth.exp.UvarExp.gnls3<-gnls(sizeNext~a+b*size+c*exp(-d*size),data=growthonly,params=list(a~utrans,b+c+d~1),start=list(a=c(1,1,1,1,1,1,1,1,1,1),b=1,c=1,d=1),weights=varExp(form=~size))
summary(growth.exp.UvarExp.gnls3)
growth.int3<-c(coefficients(growth.exp.UvarExp.gnls3)[1],coefficients(growth.exp.UvarExp.gnls3)[1]+coefficients(growth.exp.UvarExp.gnls3)[3:10])
growth.slope3<-coefficients(growth.exp.UvarExp.gnls3)[11]
growth.e13<-coefficients(growth.exp.UvarExp.gnls3)[12]
growth.e23<-coefficients(growth.exp.UvarExp.gnls3)[13]

par(mfrow=c(3,4))
for (i in 1:10){
  x<-subset(growthonly,growthonly$utrans==levels(growthonly$utrans)[i])
  plot(x$size,x$sizeNext,xlim=c(-2,10),ylim=c(-2,10))
  lines(y,growth.int3[i]+growth.slope3*y+growth.e13*exp(-growth.e23*y))
}

growth.exp.UvarExp.gnls4<-gnls(sizeNext~a+b*size+c*exp(-d*size),data=growthonly,params=list(b~utrans,a+c+d~1),start=list(a=1,b=c(1,1,1,1,1,1,1,1,1,1),c=1,d=1),weights=varExp(form=~size))
summary(growth.exp.UvarExp.gnls4)
growth.int4<-coefficients(growth.exp.UvarExp.gnls4)[11]
growth.slope4<-c(coefficients(growth.exp.UvarExp.gnls4)[1],coefficients(growth.exp.UvarExp.gnls4)[1]+coefficients(growth.exp.UvarExp.gnls4)[2:10])
growth.e14<-coefficients(growth.exp.UvarExp.gnls4)[12]
growth.e24<-coefficients(growth.exp.UvarExp.gnls4)[13]

par(mfrow=c(3,4))
for (i in 1:10){
  x<-subset(growthonly,growthonly$utrans==levels(growthonly$utrans)[i])
  plot(x$size,x$sizeNext,xlim=c(-2,10),ylim=c(-2,10))
  lines(y,growth.int4+growth.slope4[i]*y+growth.e14*exp(-growth.e24*y))
}

growth.exp.UvarExp.gnls5<-gnls(sizeNext~a+b*size+c*exp(-d*size),data=growthonly,params=list(c~utrans,a+b+d~1),start=list(a=1,b=1,c=c(1,1,1,1,1,1,1,1,1,1),d=1),weights=varExp(form=~size))
summary(growth.exp.UvarExp.gnls5)
growth.int5<-coefficients(growth.exp.UvarExp.gnls5)[11]
growth.slope5<-coefficients(growth.exp.UvarExp.gnls5)[12]
growth.e15<-c(coefficients(growth.exp.UvarExp.gnls5)[1],coefficients(growth.exp.UvarExp.gnls5)[1]+coefficients(growth.exp.UvarExp.gnls5)[2:10])
growth.e25<-coefficients(growth.exp.UvarExp.gnls5)[13]

par(mfrow=c(3,4))
for (i in 1:10){
  x<-subset(growthonly,growthonly$utrans==levels(growthonly$utrans)[i])
  plot(x$size,x$sizeNext,xlim=c(-2,10),ylim=c(-2,10))
  lines(y,growth.int5+growth.slope5*y+growth.e15[i]*exp(-growth.e25*y))
}


growth.exp.UvarExp.gnls6<-gnls(sizeNext~a+b*size+c*exp(-d*size),data=growthonly,params=list(b+d~utrans,a+c~1),start=list(a=1,b=c(1,1,1,1,1,1,1,1,1,1),c=1,d=c(1,1,1,1,1,1,1,1,1,1)),weights=varExp(form=~size))
growth.int6<-c(coefficients(growth.exp.UvarExp.gnls6)[1],coefficients(growth.exp.UvarExp.gnls6)[1]+coefficients(growth.exp.UvarExp.gnls6)[2:10])
growth.slope6<-coefficients(growth.exp.UvarExp.gnls6)[22]
growth.e16<-coefficients(growth.exp.UvarExp.gnls6)[21]
growth.e26<-c(coefficients(growth.exp.UvarExp.gnls6)[11],coefficients(growth.exp.UvarExp.gnls6)[11]+coefficients(growth.exp.UvarExp.gnls6)[12:20])

par(mfrow=c(3,4))
for (i in 1:10){
  x<-subset(growthonly,growthonly$utrans==levels(growthonly$utrans)[i])
  plot(x$size,x$sizeNext,xlim=c(-2,10),ylim=c(-2,10))
  lines(y,growth.int6[i]+growth.slope6*y+growth.e16*exp(-growth.e26[i]*y))
}
anova(growth.exp.UvarExp.gnls,growth.exp.UvarExp.gnls2)

##Look for best varFunc = varExp vs formula to varFixed
growth.exp.UvarExp.gnls1<-gnls(sizeNext~a+b*size+c*exp(-d*size),data=growthonly,params=list(b+c~utrans,a+d~1),start=list(a=c(1,1,1,1,1,1,1,1,1,1),b=c(1,1,1,1,1,1,1,1,1,1),c=1,d=1),weights=varExp(form=~size))
growth.exp.U.gnls1<-gnls(sizeNext~a+b*size+c*exp(-d*size),data=growthonly,params=list(b+c~utrans,a+d~1),start=list(a=c(1,1,1,1,1,1,1,1,1,1),b=c(1,1,1,1,1,1,1,1,1,1),c=1,d=1))
growth.exp.UvarExp.gnls1<-gnls(sizeNext~a+b*size+c*exp(-d*size),data=growthonly,params=list(b+c~utrans,a+d~1),start=list(a=c(1,1,1,1,1,1,1,1,1,1),b=c(1,1,1,1,1,1,1,1,1,1),c=1,d=1),weights=~
                                 
x<-(summary(growth.exp.UvarExp.gnls1)$sigma)^2*exp(2*coef(growth.exp.UvarExp.gnls1$modelStruct$varStruct,"exp")*y)