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

#add max wave height to each transect
wavedata<-read.table("WaveValuesMaui_Import.txt", sep="\t", header=TRUE)#read in wave data
data$avgWH <- wavedata$avgWH[match(data$utrans,wavedata$transect)]
data$maxWH <- wavedata$maxWH[match(data$utrans,wavedata$transect)]

#make year a factor
data$fyear<-factor(data$year)
data$fyearNext<-factor(data$yearNext)                          

#Add some NAs, change some column names, add probability of growth and fusion
data$fec4<-ifelse(data$fec1==1,data$fec2,NA)
index<-c("growth","dead","fissSM","fissLG","fuseSM","fuseLG","recruit","M")
values<-c(0,NA,NA,0,1,1,NA,0)
data$Pfuse<-values[match(data$fateIPM,index)]
values<-c(NA,NA,NA,NA,0,1,NA,NA)
data$PfuseLG<-values[match(data$fateIPM,index)]
values<-c(1,NA,NA,0,0,0,NA,0)
data$Pgrowth<-values[match(data$fateIPM,index)]
values<-c(0,NA,NA,0,0,0,NA,1)
data$PM<-values[match(data$fateIPM,index)]
values<-c(0,NA,NA,1,0,0,NA,0)
data$Pfiss<-values[match(data$fateIPM,index)]
values<-c(1,0,NA,1,1,1,NA,1)
data$Psurv<-values[match(data$fateIPM,index)]
data$numfiss<-data$fec4
data$sumfiss<-data$fec3
values<-c("growth","dead","fissSM","fissLG","fuseSM","fuseLG","recruit","growth")
data$fate<-values[match(data$fateIPM,index)]
data$fate<-as.factor(data$fate)


#create some size & growth measures
data$size<-log(data$asize)
data$sizeNext<-log(data$asizeNext)
data$lnpsize <-log(data$psize)
data$lnpsizeNext <-log(data$psizeNext)
data$logsize<-log10(data$asize)
data$logasizeNext<-log10(data$asizeNext)
data$size2<-(data$size^2)
data$sqrtsize<-(data$asize)^(1/2)
data$sqrtsizeNext<-(data$asizeNext)^(1/2)

data$PAratio<-data$psize/data$asize
data$grArea<-log(data$asizeNext/data$asize)
data$grPeri<-log(data$psizeNext/data$psize)
data$grPAR<-log((data$asizeNext/data$psizeNext)/(data$asize/data$psize))
data$grDel <- data$asizeNext-data$asize
data$grPct<-100*data$grDel/data$asize
data$grRatio<-data$asizeNext/data$asize

#replace blank cells with NA
data$fate_prior[data$fate_prior==""]<-NA
data$offspringNext[data$offspringNext==""]<-NA
data$fuseLGname[data$fuseLGname==""]<-NA  
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
MC$csize <- MC$size-mean(MC$size,na.rm=TRUE)

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
MCf<-droplevels(MCf)

MCfs<-subset(MC,fateIPM=="fissSM")
MCfs<-droplevels(MCfs)
#For recruit size, exclude all but recruits
MCr <- subset(MC,MC$fate_t=="R")
