
#Spline for MCg, and growth, fission and fusion alone
bleh<-smooth.spline(MCg$size,MCg$sizeNext)
plot(bleh,xlab="Size t",ylab="Size t+1")
points(MCg$size,MCg$sizeNext,col=c(1:3)[MCg$fate])

growthspline<-smooth.spline(growthonly$size,growthonly$sizeNext)
plot(growthspline,xlab="Size t",ylab="Size t+1")
points(growthonly$size,growthonly$sizeNext)

bleh<-smooth.spline(fissonly$size,fissonly$sizeNext)
plot(bleh,xlab="Size t",ylab="Size t+1")
points(fissonly$size,fissonly$sizeNext)

bleh<-smooth.spline(fuseonly$size,fuseonly$sizeNext)
plot(bleh,xlab="Size t",ylab="Size t+1")
points(fuseonly$size,fuseonly$sizeNext)

##What would a step function look like
##Use max size at t for extreme growth value for cut off

y1<-subset(y,y<2.758039)
y2<-subset(y,y>2.758039)

dsize1<-data.frame(size=y1)
dsize2<-data.frame(size=y2)

steplo<-subset(growthonly,growthonly$size<max(growthonly$size[growthonly$ES=="EG"]))
test<-lm(sizeNext~size,data=steplo)
stephi<-subset(growthonly,growthonly$size>max(growthonly$size[growthonly$ES=="EG"]))
test2<-lm(sizeNext~size,data=stephi)

plot(growthspline,xlab="Size t",ylab="Size t+1")
points(growthonly$size,growthonly$sizeNext)
lines(y1,predict(test,newdata=dsize1),xlim=c(-2,2.758039),col="green",lwd=2)
lines(y2,predict(test2,newdata=dsize2),col="blue",lwd=2)

###
y3<-subset(y,y<1.5)
y4<-subset(y,y>1.5)
dsize3<-data.frame(size=y3)
dsize4<-data.frame(size=y4)

steplo1<-subset(growthonly,growthonly$size<1.5)
test3<-lm(sizeNext~size,data=steplo1)
stephi1<-subset(growthonly,growthonly$size>1.5)
test4<-lm(sizeNext~size,data=stephi1)

plot(growthspline,xlab="Size t",ylab="Size t+1")
points(growthonly$size,growthonly$sizeNext)
lines(y3,predict(test3,newdata=dsize3),col="purple",lwd=2)
lines(y4,predict(test4,newdata=dsize4),col="magenta",lwd=2)

###
dsize5<-data.frame(size=y1,size2=y1^2)
test5<-lm(sizeNext~size+size2,data=steplo)
plot(growthspline,xlab="Size t",ylab="Size t+1")
points(growthonly$size,growthonly$sizeNext)

lines(y1,predict(test5,newdata=dsize5),lwd=2,col="green")
lines(y2,predict(test2,newdata=dsize2),col="blue",lwd=2)
