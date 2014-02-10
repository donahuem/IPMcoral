fusever<-subset(MC,MC$fateIPM=="fuseLG"|MC$fateIPM=="fuseSM"|MC$fateIPM=="M")
Sfuse<-subset(MC,MC$fateIPM=="fuseSM")
Sfuse$test<-paste(Sfuse$uyear,Sfuse$fuseLGname)
fuseonly$test<-paste(fuseonly$uyear, fuseonly$colony)
index<-fuseonly$test
values<-fuseonly$asize

Sfuse$LGsize<-values[match(Sfuse$test,index)]
Sfuse$test2<-Sfuse$asize-Sfuse$LGsize
