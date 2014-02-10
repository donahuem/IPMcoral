
#Are elaticities correlated with lambda?
eig<-as.data.frame(eig)
eig$utrans<-levels(MC$utrans)

#What are the highest elasticities by site
#rank elasticities by site
bleh<-matrix(0,nrows=10,ncols=34)
for(i in 1:10){
  bleh[i,]<-as.numeric(factor(abs(elasP)[i,]))}

#Sum columns in ranking to see which elasticities are most important across sites
s<-colSums(bleh)
barplot(s)

#rank columns sums
sr<-as.numeric(factor(s))



for (i in c(2:3,5:7,9:11,16:18,20,21,28:31,33,34)){
  x<-cor.test(eig$eig,elasP[,i])
  print(x)
}
#plot parameters with significant correlation with lambda
par(mfrow=c(2,2))
for (i in c(2:3,5:7,9:11,16:18,20,21,28:31,33,34)){
  plot(eig$eig,elasP[,i],ylab=rownames(params)[i],col=cl,pch=18,xlab="lambda")
}
#of these the following do not look reasonably linear:
# they all have the same outlier below I list one's where the outlier seems to have most influence
#r looks linear but I think that the outlier is driving it
#recruit.size.sd outlier
#Philo.slope doesn't look too linear
#Philo.int outlier
#Pfiss.int outlier
#fiss.int not too linear

for (i in c(1:34)){
  x<-cor.test(slopes$slopes,elasP[,i])
  print(x)
}

par(mfrow=c(2,2))
for (i in c(2:3,18,28:31)){
  plot(slopes$slopes,elasP[,i],ylab=rownames(params)[i],col=cl,pch=18)
}

par(mfrow=c(2,2))
for (i in c(2,18,28:29)){
  plot(slopes$slopes,elasP[,i],ylab=rownames(params)[i],col=cl,pch=18)
}