#GET AND ORGANIZE THE DATA
rm(list=ls(all=TRUE)) ## Clear all
setwd("C://Users/Megan/IPMcoral")# for MeganR
getwd() #check working directory is set properly
data<- read.delim("sortIPMpack_28Mar2013.txt")

#create site names
sitenametable<-c("HonN",
                 "HonS",
                 "Molo",
                 "Puam",
                 "Kahe",
                 "KanB",
                 "KanP",
                 "Maal",
                 "Olow",
                 "Papa")
sitetable<-c("Honolua",
             "Honolua",
             "Molokini",
             "Puamana",
             "Kahekili",
             "KanahenaBay",
             "KanahenaPt",
             "Maalaea",
             "Olowalu",
             "Papaula")
data$site<-sitetable[charmatch(data$sitename,sitenametable)]
data$site<-factor(data$site)
data$quad<-factor(paste("Q",data$quad,sep=""))

#create transect names
data$trans <- as.character(data$depthname)
data <- within(data,trans[sitename=="HonS"]<-"S")
data <- within(data,trans[sitename=="HonN"]<-"N")
data$trans <-factor(data$trans)

#create unique names for transect, quad, and colony
data$utrans<-factor(paste(data$site,data$trans,sep="."))
data$uquad <-factor(paste(data$utrans,data$quad,sep="."))
data$uyear <-factor(paste(data$utrans,data$quad,data$year,sep="."))
data$tyear <-factor(paste(data$utrans,data$year,sep="."))
data$ucolony <-factor(paste(data$utrans,data$quad,data$species,data$colony,sep="."))

#make year a factor
data$fyear<-factor(data$year)
data$fyearNext<-factor(data$yearNext)                          

#Add some NAs, change some column names, add probability of growth and fusion
#replace blank cells with NA
data$fate_prior[data$fate_prior==""]<-NA
data$offspringNext[data$offspringNext==""]<-NA
data$fuseLGname[data$fuseLGname==""]<-NA  
data$fec4<-ifelse(data$fec1==1,data$fec2,NA)
# change some column names
data$numfiss<-data$fec4
data$sumfiss<-data$fec3

#Add probability of fission, all & large only fusion, survival
index<-c("growth","dead","fissSM","fissLG","fuseSM","fuseLG","recruit","M")
values<-c(1,NA,NA,0,0,0,NA,0)
data$Pgrowth<-values[match(data$fateIPM,index)] # Keep this.  We need it for plots. Estimate the probability of growth as the difference between 1 and the probability of fission and fusion. 
values<-c(0,NA,NA,0,1,1,NA,1)
data$Pfuse<-values[match(data$fateIPM,index)] # To estimate probability of All fusion. Includes M
values<-c(NA,NA,NA,NA,0,1,NA,1)
data$PfuseLG<-values[match(data$fateIPM,index)] # To estimate probability of Large fusion given fusion. Includes M 
values<-c(0,NA,NA,1,0,0,NA,1)
data$Pfiss<-values[match(data$fateIPM,index)] # To estiamte probability of fission. Includes M (need it for clonal reproduction)
values<-c(1,0,NA,1,1,1,NA,1)
data$Psurv<-values[match(data$fateIPM,index)] # To estimate probability of survival
values<-c("growth","dead","fissSM","fissLG","fuseSM","fuseLG","recruit","fuseLG")
data$fate<-values[match(data$fateIPM,index)]# Make fate factor to include M in fuse LG for growth function
data$fate<-as.factor(data$fate)


#create some size & growth measures
data$size<-log(data$asize)
data$sizeNext<-log(data$asizeNext)
#data$lnpsize <-log(data$psize)
#data$lnpsizeNext <-log(data$psizeNext)
#data$logsize<-log10(data$asize)
#data$logasizeNext<-log10(data$asizeNext)
data$size2<-(data$size^2)
#data$sqrtsize<-(data$asize)^(1/2)
#data$sqrtsizeNext<-(data$asizeNext)^(1/2)

#data$PAratio<-data$psize/data$asize
data$grArea<-log(data$asizeNext/data$asize)
#data$grPeri<-log(data$psizeNext/data$psize)
#data$grPAR<-log((data$asizeNext/data$psizeNext)/(data$asize/data$psize))
#data$grDel <- data$asizeNext-data$asize
#data$grPct<-100*data$grDel/data$asize
#data$grRatio<-data$asizeNext/data$asize


#check format of data
str(data) #check structure of data and variables
#exclude 2-year observations  MIGHT WANT TO RECONSIDER THIS CHOICE LATER
data <- subset(data,data$twoyear==0)
data<-droplevels(data)
#check which sites should be excluded for which species b/c too few colonies
  table(data$utrans,data$species)
#                 MC  MF  MP  PB  PC  PD  PL  PM  PR  PV
#Honolua.N       391  39  13   0 197   0 568   9   0  11
#Honolua.S       101 239  52   0   0   0 132   9   0  51
#Kahekili.03m    165   0  65   0 246   0 365  12   0   0
#Kahekili.07m    245   0  11   0 513   0 354  24   0   0
#KanahenaBay.03m  25   0 118   0 153   0  76  26   0 220
#KanahenaPt.03m    0   0  29   0   0  22   5  27   0  39
#KanahenaPt.10m    7   0  36   0 244  25   0   3 159 106
#Molokini.07m    144   0 192   0 367   0 323  15   0  28
#Molokini.13m    194   0 139   8 137   0  45  21   0  13
#Olowalu.03m     400   0  42   0 355 241 246  11   0   0
#Olowalu.07m     260   0 373   0   0   0  29  37   0   0
#Papaula.03m       4 386  25   0   0   0 130   6   0   0
#Puamana.03m     275   0  33   0  95   6 276 204   0   0
#Puamana.13m      19   0   2   0 381   0  89  12   0   0
 
#exclude transects with fewer than 20 growth observations
MCi2<-data[data$species=="MC" & data$site!="KanahenaBay"  & data$site!="KanahenaPt" & data$utrans!="Papaula.03m" & data$utrans!="Puamana.13m",]
  MCi2$utrans<-factor(MCi2$utrans)
  MCi2$site<-factor(MCi2$site)
  MCi2$uquad<-factor(MCi2$uquad)
#create a variable:  spcoverD = spcover_t-average transect cover over years and quads
attach(MCi2)
coverTQY <- tapply(spcover,list(utrans,uquad,year),mean,na.rm=T)
coverT <- rowMeans(coverTQY,na.rm=T)
MCi2$spcoverD<-(spcover-coverT[as.character(utrans)])/coverT[as.character(utrans)]
detach(MCi2)

MC<-MCi2

#Need this centered variable later
#MC$csize <- MC$size-mean(MC$size,na.rm=TRUE)

MCg<-subset(MC,fate=="growth"|fate=="fuseLG"|fate=="fissLG")
amc={};
for (i in 1:length(MCg)) {
  if (!(sum(is.na(MCg[i]))>0)) {
    amc<-c(amc,i);
    #str(MCa[i])
  }
}
MCg<-droplevels(MCg[,amc])

MCg<-droplevels(MCg)

MCf<-subset(MC,fateIPM=="fissLG"|fateIPM=="M")
xtabs(~numfiss,MCf)
index<-c(1,2,3,4,5,6,7,10)
values<-c("one","two","three",">three",">three",">three",">three",">three")
MCf$fnumfiss<-values[match(MCf$numfiss,index)]
values<-c("one","two","three","four","five","six","seven","ten")
MCf$fnumfiss2<-values[match(MCf$numfiss,index)]
MCf$fnumfiss<-as.factor(MCf$fnumfiss)
MCf$fnumfiss2<-as.factor(MCf$fnumfiss2)
xtabs(~fnumfiss,MCf)
MCf<-droplevels(MCf)

MCfs<-subset(MC,fateIPM=="fissSM")
MCfs<-droplevels(MCfs)
#For recruit size, exclude all but recruits
MCrs <- subset(MC,MC$fate_t=="R")


#################################################
#Data for Recruitment
library(gdata)

#Use alldata.txt (just a list of colonies and sizes with photoquad info) to get the recruit info
rawdata<-data<- read.delim("alldata_28Mar2013.txt")
rawdata$site<-sitetable[charmatch(rawdata$sitename,sitenametable)]
rawdata$site<-factor(rawdata$site)
rawdata$quad<-factor(paste("Q",rawdata$quad,sep=""))

#create transect names just like for sortdata->data
rawdata$trans <- as.character(rawdata$depthname)
rawdata <- within(rawdata,trans[sitename=="HonS"]<-"S")
rawdata <- within(rawdata,trans[sitename=="HonN"]<-"N")
rawdata$trans <-factor(rawdata$trans)

#create unique names for transect, quad, and colony just like for sortdata->data
rawdata$utrans<-factor(paste(rawdata$site,rawdata$trans,sep="."))
rawdata$uquad <-factor(paste(rawdata$utrans,rawdata$quad,sep="."))
rawdata$uyear <-factor(paste(rawdata$utrans,rawdata$quad,rawdata$year_t,sep="."))
rawdata$tyear <-factor(paste(rawdata$utrans,rawdata$year_t,sep="."))

#check for missing quads
#checking with MR, but for now
#No Data for Papaula 03m Q2 2003, Puamana 13m Q1 (exclude 2004 and 2005; earlier years missing), Puamana.13m.Q5.2001
#Puamana 13m Q4 does not exist for all years,
table(rawdata$uquad,rawdata$year_t)
RFtab<-xtabs(fate_t=="R"~uyear+species,rawdata,drop.unused.levels=FALSE)     
RF<-as.data.frame(RFtab)
#RF <- RF[(RF$uyear != "Papaula.03m.Q2.2003" & RF$uyear != "Puamana.13m.Q5.2001"),]
#RF <- RF[(RF$uyear != "Puamana.13m.Q1.2004"& RF$uyear != "Puamana.13m.Q1.2005"),]
#RF <- RF[(RF$uyear != "Papaula.03m.Q2.2003" & RF$uyear != "Puamana.13m.Q5.2001"),]

a<-unlist(strsplit(as.character(RF$uyear),"\\."))
newcol<-matrix(a,ncol=4,byrow=TRUE)
RF<-cbind(RF,newcol)
RF$utrans<-paste(RF[,4],RF[,5])
RF$quad<-RF[,6]
RF$year<-RF[,7]
RF<-subset(RF,select=c(1,2,3,8,9,10))
RF$year<-as.factor(RF$year)
RF$species<-as.factor(RF$species)
RF$utrans<-as.factor(RF$utrans)
RF$quad<-as.factor(RF$quad)
#create density variable
spcoverlookup <- as.data.frame.table(tapply(rawdata$area_t,list(rawdata$uyear,rawdata$species),sum,na.rm=TRUE))
spcoverlookup$Freq[is.na(spcoverlookup$Freq)] <- 0  #these NAs indicated that none of that species was present, so it's a zero
allcoverlookup <- as.data.frame.table(tapply(rawdata$area_t,list(rawdata$uyear),sum,na.rm=TRUE))
RF$allcover_t <- allcoverlookup$Freq[allcoverlookup$Var1==as.character(RF$uyear)]
RF$spcover_t <- spcoverlookup$Freq
RF <- rename.vars(RF,from="Freq",to="recruits")

table(RF$utrans,RF$year)
#OK, looking at this table, we can see that we have the correct number of observations in each quad/year:  
#Missing Quads:
#KanP 03m Q2 in 2004, 2005 (blank log and layer files)
#Papa 03m Q2 in 2003, 2004 (missing output file)
#Puam 13m Q1 in all years (blank log and layer files)
#Puam 13m Q4 in 1999-2003 (blank log and layer files)
#Puam 13m Q5 in 2001      (missing photo file)
#For these observations, zeros are real zeros.

#This approach includes the two year transitions, which I think is ok b/c obs are so low.
#   this really only effects Honolua (for all three spp) and Papaula.03m for PL.  

#Remember to exclude the first year of observation  b/c those are false zeros -- all "N" fates
#   1st years:  1999(HonN,HonS,Olow,Papa,Puam), 2000(Molo), 2001 (Kahe,Kana)

# exclude the first photos from each site in appropriate year -- all NEW; 

RF <-subset(RF,!(year==1999))
RF <-subset(RF,!(year==2000 & (utrans=="Molokini 07m" | utrans=="Molokini 13m")))
RF <-subset(RF,!(year==2001 & (utrans=="KanahenaBay 03m" | utrans=="KanahenaPt 03m" | utrans=="KanahenaPt 10m" | utrans=="Kahekili 03m" | utrans=="Kahekili 07m")))
table(RF$utrans,RF$year)

RFMC <- droplevels(subset(RF,RF$species=="MC"))
RFPC <- droplevels(subset(RF,RF$species=="PC"))
RFPL <- droplevels(subset(RF,RF$species=="PL"))

#For the IPM, subset to species, and only look at sites where there are enough observations to measure growth, surv, and fission, too
MCr <- droplevels(subset(RF,RF$species=="MC" & (RF$utrans!="KanahenaBay 03m" & RF$utrans!="KanahenaPt 03m"& RF$utrans!="KanahenaPt 10m"& RF$utrans!="Papaula 03m" & RF$utrans!="Puamana 13m")))
#MCr <- droplevels(subset(RF,RF$species=="MC" & (RF$utrans!="KanahenaBay 03m" & RF$utrans!="Honolua.S" & RF$utrans!="KanahenaPt 03m"& RF$utrans!="KanahenaPt 10m"& RF$utrans!="Papaula 03m" & RF$utrans!="Puamana 13m")))
PLr <- droplevels(subset(RF,RF$species=="PL" & (RF$utrans!="KanahenaPt 10m" & RF$utrans!="KanahenaPt 03m" & RF$utrans!="Olowalu 07m")))
PCr <- droplevels(subset(RF,RF$species=="PC" & (RF$utrans!="Honolua S" & RF$utrans!="KanahenaPt 03m"  & RF$utrans!="Olowalu 07m" & RF$utrans!="Papaula 03m")))
##########################################################################
#Data for Density
MCds<-subset(MC,fateIPM!="recruit"&fateIPM!="fissSM")
MCds<-droplevels(MCds)
MCd<-xtabs(~uyear,MCds)
MCd<-as.data.frame(MCd)
MCd$uyear<-as.character(MCd$uyear)
temp<-matrix(nrow=190,ncol=4)
for (i in 1:190){
  temp[i,]<-unlist(strsplit(MCd[i,1],"[.]"))

}
for(i in 1:190){
  temp2<-paste(temp[,1],temp[,2])
}
MCd<-cbind(MCd,temp2,temp[,3],temp[,4])
colnames(MCd)<-c("uyear","colnum","utrans","quad","year")




###################################################################
MCf$hilo<-ifelse(MCf$numfiss<4,"lo","hi")
MCf$hilo<-as.factor(MCf$hilo)
hinum<-subset(MCf,MCf$numfiss>3)
lonum<-subset(MCf,MCf$numfiss<4)
MCf$Philo<-ifelse(MCf$hilo=="lo",1,0)

MCg$sizediff<-(MCg$sizeNext-MCg$size)
growthonly<-subset(MCg,MCg$fateIPM=="growth")

dl<-qnorm(0.05,mean(growthonly$sizediff),sd(growthonly$sizediff))
growthonly$ES<-ifelse(growthonly$sizediff<dl,"ES","NGS")
growthonly$ES<-as.factor(growthonly$ES)
GONGS<-subset(growthonly,growthonly$ES=="NGS")
GOES<-subset(growthonly,growthonly$ES=="ES")
index<-c("ES","NGS")
values<-c(1,0)
growthonly$PES<-values[match(growthonly$ES,index)]
fissonly<-subset(MCg,MCg$fateIPM=="fissLG")
fuseonly<-subset(MCg,MCg$fateIPM=="fuseLG"|MCg$fateIPM=="M")
