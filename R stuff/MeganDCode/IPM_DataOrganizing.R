#GET AND ORGANIZE THE DATA
#21June2012 - Cleaned up and annotated 
#18June2012 - created MCr/PLr/PCr from alldata.txt to coutn recruits into each observed quad
#7June2012 - fixed the MCa/PLa/PCa to include fate_t1=M; cleaned up the $fiss, $fuse definitions.
#31May2012 - large fissers and large fusers get included in growth (MCa), small fusers get included with death
#31May2012 - fixed site exclusion names so MC, PL, and PC all exclude sites with no (or very few) observations with the correct site names

#rm(list=ls(all=TRUE)) ## Clear all

setwd("~/MJD/Projects/CAMEO/CAMEO-R") #set working directory to location of this script & data file

getwd() #check working directory is set properly

#read in data
data<-read.table("C:/Documents and Settings/megan/My Documents/MJD/Projects/CAMEO/Parsing w Perl/sortdata_18June2012.txt", sep="\t", header=TRUE)#read in data

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
data$uyear <-factor(paste(data$utrans,data$quad,data$year_t,sep="."))
data$tyear <-factor(paste(data$utrans,data$year_t,sep="."))
data$ucolony <-factor(paste(data$utrans,data$quad,data$species,data$colony,sep="."))

#add max wave height to each transect
wavedata<-read.table("WaveValuesMaui_Import.txt", sep="\t", header=TRUE)#read in wave data
data$avgWH <- wavedata$avgWH[match(data$utrans,wavedata$transect)]
data$maxWH <- wavedata$maxWH[match(data$utrans,wavedata$transect)]

#make year a factor
data$fyear_t<-factor(data$year_t)
data$fyear_t1<-factor(data$year_t1)                          

#create some additional variables using logical conditions on fate
  #recall that fate_t1 exists only for colonies that existed at fate_t (therefore, R and N are never in column fate_t1)
data$surv = NA   #SURVIVAL -- consider fusing death for smaller fusers
  data <- within(data,surv[fate_t1=="D"]<-0) #if you're dead, you're dead
  data$surv[(data$fate_t1=="U" & data$fusee_t1!="")]<-0 #if you fused and were the smaller one (i.e., you have a fusee) then dead
  data$surv[(data$fate_t1=="U" & data$fusee_t1=="")]<-1 #if you fused and were the larger one (i.e., no fusee) then survived
  data <- within(data,surv[fate_t1=="S" | fate_t1=="M" | fate_t1=="I"]<-1) #if you grew, fissed, or fissed and fused, then you survived'
  #fate_t1 is never R
data$fiss = NA   #FISSION
  data <- within(data,fiss[fate_t1=="I"|fate_t1=="M"]<-1)
  #if fate_t1=I or M then fissed.  Otherwise, no.  Want cases (I, M, D, S, U) in the denominator
  #leave as NA those cases where fate_t1==""; this includes fate_t=R(in last year of photo)
  data <- within(data,fiss[fate_t1=="D"|fate_t1=="S"|fate_t1=="U"]<-0)
data$fuse = NA   #FUSION
  data <- within(data,fuse[fate_t1=="U"|fate_t1=="M"]<-1)
  data <- within(data,fuse[fate_t1=="D"|fate_t1=="S"|fate_t1=="I"]<-0)
  #leave as NA those cases where fate_t1==""; this includes fate_t=R(in last year of photo)


#create some size & growth measures
data$lnarea_t<-log(data$area_t)
data$lnarea_t1<-log(data$area_t1)
data$lnperi_t <-log(data$peri_t)
data$lnperi_t1 <-log(data$peri_t1)
data$logarea_t<-log10(data$area_t)
data$logarea_t1<-log10(data$area_t1)
data$lnareafiss <- log(data$areafiss)
data$lnarea_prefiss <- log(data$area_prefiss)


data$PAratio<-data$peri_t/data$area_t
data$grArea<-log(data$area_t1/data$area_t)
data$grPeri<-log(data$peri_t1/data$peri_t)
data$grPAR<-log((data$area_t1/data$peri_t1)/(data$area_t/data$peri_t))
data$grDel <- data$area_t1-data$area_t
data$grPct<-100*data$grDel/data$area_t
data$grRatio<-data$area_t1/data$area_t
  
#check format of data
  str(data) #check structure of data and variables
#check which sites should be excluded for which species b/c too few colonies
  table(data$utrans,data$species)

#                 MC  MF  MP  PB  PC  PD  PL  PM  PR  PV
#Honolua.N       324  31  11   0 150   0 485   8   0  10
#Honolua.S        80 186  46   0   0   0 119   8   0  39
#Kahekili.03m    127   0  49   0 181   0 275  11   0   0
#Kahekili.07m    190   0   8   0 408   0 275  20   0   0
#KanahenaBay.03m  14   0  93   0 113   0  57  22   0 169
#KanahenaPt.03m    0   0  19   0   0  16   4  25   0  29
#KanahenaPt.10m    6   0  34   0 190  19   0   2 106  89
#Molokini.07m    109   0 162   0 299   0 289  14   0  25
#Molokini.13m    155   0 121   8 105   0  41  21   0  11
#Olowalu.03m     306   0  33   0 255 176 190  10   0   0
#Olowalu.07m     195   0 277   0   0   0  25  34   0   0
#Papaula.03m       3 264  14   0   0   0 102   4   0   0
#Puamana.03m     217   0  23   0  75   3 236 166   0   0
#Puamana.13m      12   0   1   0 269   0  71   8   0   0

#create separate dataframes for three most abundant species

#exclude transects with fewer than 20 growth observations
MCi2<-data[data$species=="MC" & data$site!="KanahenaBay"  & data$site!="KanahenaPt" & data$site!="Papaula" & data$utrans!="Puamana.13m",]
#MC<-data[data$species=="MC" & !data$skipyear & data$site!="KanahenaBay" & data$utrans != "Honolua.S" & data$site!="KanahenaPt" & data$site!="Papaula" & data$utrans!="Puamana.13m",]
  MCi2$utrans<-factor(MCi2$utrans)
  MCi2$site<-factor(MCi2$site)
  MCi2$uquad<-factor(MCi2$uquad)
PLi2<-data[data$species=="PL" & data$utrans!="KanahenaPt.10m" & data$utrans!="KanahenaPt.03m" & data$utrans!="Olowalu.07m",]
  PLi2$utrans<-factor(PLi2$utrans)
  PLi2$site<-factor(PLi2$site)
  PLi2$uquad<-factor(PLi2$uquad)
PCi2<-data[data$species=="PC" & data$utrans!="Honolua.S" & data$utrans!="KanahenaPt.03m"  & data$utrans!="Olowalu.07m" & data$utrans!="Papaula.03m",]
  PCi2$utrans<-factor(PCi2$utrans)
  PCi2$site<-factor(PCi2$site)
  PCi2$uquad<-factor(PCi2$uquad)

#create a variable:  spcoverD = spcover_t-average transect cover over years and quads
attach(MCi2)
coverTQY <- tapply(spcover_t,list(utrans,uquad,year_t),mean,na.rm=T)
coverT <- rowMeans(coverTQY,na.rm=T)
MCi2$spcoverD<-(spcover_t-coverT[as.character(utrans)])/coverT[as.character(utrans)]
detach(MCi2)
attach(PCi2)
coverTQY <- tapply(spcover_t,list(utrans,uquad,year_t),mean,na.rm=T)
coverT <- rowMeans(coverTQY,na.rm=T)
PCi2$spcoverD<-(spcover_t-coverT[as.character(utrans)])/coverT[as.character(utrans)]
detach(PCi2)
attach(PLi2)
coverTQY <- tapply(spcover_t,list(utrans,uquad,year_t),mean,na.rm=T)
coverT <- rowMeans(coverTQY,na.rm=T)
PLi2$spcoverD<-(spcover_t-coverT[as.character(utrans)])/coverT[as.character(utrans)]
detach(PLi2)

#exclude 2-year observations  MIGHT WANT TO RECONSIDER THIS CHOICE LATER
MC <- subset(MCi2,MCi2$skipyear==0)
PL <- subset(PLi2,PLi2$skipyear==0)
PC <- subset(PCi2,PCi2$skipyear==0)

#For recruit size, exclude all but recruits
MCrs <- subset(MCi2,MCi2$fate_t=="R")
PCrs <- subset(PCi2,PCi2$fate_t=="R")
PLrs <- subset(PLi2,PLi2$fate_t=="R")

#Need this centered variable later
MC$clnarea_t <- MC$lnarea_t-mean(MC$lnarea_t)
PC$clnarea_t <- PC$lnarea_t-mean(PC$lnarea_t)
PL$clnarea_t <- PL$lnarea_t-mean(PL$lnarea_t)

#data for GROWTH model, 
  #at first, I only included S in growth
    #MCs<-MC[MC$fate_t1=="S",]
    #PLs<-PL[PL$fate_t1=="S",]
    #PCs<-PC[PC$fate_t1=="S",]
#But we decided the cleanest way to break up the growth/fission/fusion was to consider growth inclusive of 
  #fate_t1=S is grow from t to t1
  #fate_t1=I means colony fissed (and survived) from t to t1
  #fate_t1=U means colony fused from t to t1, and need to distinguish between fuser (=dead) or fusee (=growth)
#NOTE(6-7-12): also need to include fate_t1=M if no fusee
MCa<-MC[(MC$fate_t1=="S" | MC$fate_t1 =="I" | (MC$fate_t1=="U" & MC$fusee_t1 == "") | (MC$fate_t1=="M" & MC$fusee_t1 == "")),]
PCa<-PC[(PC$fate_t1=="S" | PC$fate_t1 =="I" | (PC$fate_t1=="U" & PC$fusee_t1 == "") | (PC$fate_t1=="M" & PC$fusee_t1 == "")),]
PLa<-PL[(PL$fate_t1=="S" | PL$fate_t1 =="I" | (PL$fate_t1=="U" & PL$fusee_t1 == "") | (PL$fate_t1=="M" & PL$fusee_t1 == "")),]


#For growth analysis, excludes columns with empty cells which seems to confuse NLME
amc={};
for (i in 1:length(MCa)) {
  if (!(sum(is.na(MCa[i]))>0)) {
    amc<-c(amc,i);
    #str(MCa[i])
  }
}
apl={};
for (i in 1:length(PLa)) {
  if (!(sum(is.na(PLa[i]))>0)) {
    apl<-c(apl,i);
    #str(PLa[i])
  }
}
apc={};
for (i in 1:length(PCa)) {
  if (!(sum(is.na(PCa[i]))>0)) {
    apc<-c(apc,i);
    #str(PCa[i])
  }
}

MCa<-droplevels(MCa[,amc])
PCa<-droplevels(PCa[,apc])
PLa<-droplevels(PLa[,apl])

#################################################
#Data for Recruitment
library(gdata)

#Use alldata.txt (just a list of colonies and sizes with photoquad info) to get the recruit info
rawdata<-read.table("C:/Documents and Settings/megan/My Documents/MJD/Projects/CAMEO/Parsing w Perl/alldata_18June2012.txt", sep="\t", header=TRUE)#read in data
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

